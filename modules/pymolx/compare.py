'''
Structure comparison (command "compare_structures", Structure Comparison
dialog) of closely related proteins (allelic variants, homologues from
other species, models): all structures are aligned and superimposed on a
reference, one multiple alignment object is made, and each residue gets a
structural similarity score, shown as colors from blue (similar) through
green and yellow to red (different); residues without a partner are grey.

Scores, per aligned residue (C-alpha atoms):

- lDDT (default; Mariani et al. 2013, as in CASP and AlphaFold's pLDDT):
  the fraction of the residue's distances to the reference residues
  within 15 A that are kept within 0.5, 1, 2 and 4 A (averaged), 0-100.
  Needs no superposition, so a domain moved as a whole stays blue; only
  local changes (loops, mutation sites) show.
- distance: C-alpha distance after the superposition, A.

Compared to the reference: each structure is colored by its similarity to
the reference, the reference by the mean of the others. All against all:
every structure is colored by the mean over all pairs at that alignment
position (how conserved the position is across the set).
'''

import math

import numpy as np

from pymol import cmd

METHODS = ('align', 'super', 'cealign')
METRICS = ('lddt', 'distance')
MODES = ('reference', 'all')
LDDT_RADIUS = 15.0
LDDT_THRESHOLDS = (0.5, 1.0, 2.0, 4.0)
# default color ranges: score at red and at blue (close relatives keep
# lDDT above 80 nearly everywhere: 70 gives contrast)
RANGES = {'lddt': (70.0, 100.0), 'distance': (3.0, 0.0)}
UNALIGNED_COLOR = 'grey50'
RAMP = ['red', 'yellow', 'green', 'blue']  # low to high similarity
RAMP_STEPS = 100


def _objects(objects, _self):
    import fnmatch
    names = [n for n in (_self.get_object_list('all') or [])
             if _self.get_type(n) == 'object:molecule']
    if not str(objects).strip():
        return names
    found = []
    for pattern in str(objects).replace(',', ' ').split():
        found += [n for n in names
                  if fnmatch.fnmatchcase(n, pattern) and n not in found]
    return found


def _ca(obj, selection, state, _self):
    '''
    :return: {atom index: (chain, resi, resn, xyz)} of C-alpha atoms
    '''
    atoms = {}
    _self.iterate_state(
        state, '%%%s and (%s) and name CA and alt \'\'+A and present' % (
            obj, selection),
        'atoms[index] = (chain, resi, resn, (x, y, z))',
        space={'atoms': atoms})
    return atoms


def lddt(ref_xyz, model_xyz, ref_all):
    '''
    Per-residue C-alpha lDDT of the model, against the reference.

    :param ref_xyz: (n, 3) the reference's aligned C-alphas
    :param model_xyz: (n, 3) the model's partners, rows matching
    :param ref_all: (m, 3) all the reference's C-alphas (the aligned ones
        first, in the same order): neighbors the model lacks count as
        not kept
    :return: (n,) scores 0-100 (nan: no neighbors)
    '''
    n = len(ref_xyz)
    if not n:
        return np.zeros(0)
    d_ref = np.linalg.norm(ref_xyz[:, None] - ref_all[None], axis=-1)
    d_model = np.linalg.norm(model_xyz[:, None] - model_xyz[None], axis=-1)
    near = d_ref < LDDT_RADIUS
    near[np.arange(n), np.arange(n)] = False
    counts = near.sum(axis=1)
    diff = np.abs(d_ref[:, :n] - d_model)
    kept = 0.0
    for t in LDDT_THRESHOLDS:
        kept = kept + ((diff < t) & near[:, :n]).sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        return 100.0 * kept / (len(LDDT_THRESHOLDS) * counts)


class Result:
    '''
    The alignment and the scores

    - objects, reference, alignment (object name)
    - columns: [{object: atom index}] in alignment order
    - residues: {object: {atom index: (chain, resi, resn, xyz)}}
    - scores: {object: {atom index: score}}
    - stats: {object: {'aligned', 'identity', 'rmsd', 'score'}}
    '''

    def __init__(self):
        self.objects = []
        self.reference = ''
        self.alignment = ''
        self.metric = 'lddt'
        self.mode = 'reference'
        self.columns = []
        self.residues = {}
        self.scores = {}
        self.stats = {}

    def profile(self, obj):
        '''
        :return: [(column number, score)] along the alignment
        '''
        return [(k, self.scores[obj][col[obj]])
                for k, col in enumerate(self.columns)
                if obj in col and col[obj] in self.scores[obj]]

    def table(self):
        label = 'lDDT' if self.metric == 'lddt' else 'dist (A)'
        lines = [' %-24s %8s %9s %9s %9s' % (
            'Structure', 'Aligned', 'Identity', 'RMSD (A)', 'Mean ' + label)]
        for obj in self.objects:
            s = self.stats[obj]
            lines.append(' %-24s %8s %8s %9s %9s' % (
                obj + (' (reference)' if obj == self.reference else ''),
                s['aligned'],
                '%.0f%%' % s['identity'] if s['identity'] is not None
                else '-',
                '%.2f' % s['rmsd'] if s['rmsd'] is not None else '-',
                '%.1f' % s['score'] if s['score'] is not None else '-'))
        return '\n'.join(lines)

    def save_csv(self, path):
        import csv
        with open(path, 'w', newline='') as handle:
            out = csv.writer(handle)
            out.writerow(['object', 'chain', 'resi', 'resn', 'column',
                          self.metric])
            for obj in self.objects:
                column_of = {col[obj]: k for k, col in enumerate(self.columns)
                             if obj in col}
                for index, (chain, resi, resn, _) in sorted(
                        self.residues[obj].items()):
                    score = self.scores[obj].get(index)
                    out.writerow([obj, chain, resi, resn,
                                  column_of.get(index, ''),
                                  '' if score is None else '%.2f' % score])


def compute(objects='', reference='', selection='polymer.protein',
            method='align', metric='lddt', mode='reference', state=1,
            alignment='compare_aln', _self=cmd):
    '''
    Align and superimpose all on the reference, score every residue.

    :return: Result
    '''
    method, metric, mode = (str(v).strip().lower()
                            for v in (method, metric, mode))
    for value, allowed, what in ((method, METHODS, 'method'),
                                 (metric, METRICS, 'metric'),
                                 (mode, MODES, 'mode')):
        if value not in allowed:
            raise ValueError('%s must be one of: %s' % (what,
                                                        ', '.join(allowed)))
    names = _objects(objects, _self)
    if len(names) < 2:
        raise ValueError('need at least two structures (objects)')
    reference = str(reference).strip() or names[0]
    if reference not in names:
        raise ValueError('reference "%s" is not among the structures' %
                         reference)
    names = [reference] + [n for n in names if n != reference]
    state = int(state)
    selection = str(selection).strip() or 'polymer.protein'

    result = Result()
    result.objects, result.reference = names, reference
    result.metric, result.mode = metric, mode
    for obj in names:
        result.residues[obj] = _ca(obj, selection, state, _self)
        if len(result.residues[obj]) < 3:
            raise ValueError('fewer than 3 C-alpha atoms in "%s and (%s)"' %
                             (obj, selection))

    # one multiple alignment with all residue pairs (no outlier rejection:
    # that would drop the moved regions), then all superimposed on the
    # reference by their common core (with outlier rejection)
    result.alignment = _self.get_legal_name(alignment)
    _self.delete(result.alignment)
    atoms = '(%s) and (%s) and name CA and alt \'\'+A' % (
        ' or '.join('%' + n for n in names), selection)
    every = {'cycles': 0} if method in ('align', 'super') else {}
    _self.extra_fit(atoms, reference, method=method, object=result.alignment,
                    transform=0, quiet=1, **every)
    _self.extra_fit(atoms, reference, method=method, quiet=1)
    _self.disable(result.alignment)  # (its lines are clutter here)
    for column in _self.get_raw_alignment(result.alignment) or []:
        entry = {obj: index for obj, index in column
                 if index in result.residues.get(obj, {})}
        if len(entry) > 1:
            result.columns.append(entry)
    # coordinates after the superposition
    for obj in names:
        result.residues[obj] = _ca(obj, selection, state, _self)

    def xyz(obj, index):
        return result.residues[obj][index][3]

    pair_scores = {}

    def scores_of(model, ref):
        '''
        {model atom index: score} for the model against ref
        '''
        if (model, ref) in pair_scores:
            return pair_scores[(model, ref)]
        cols = [c for c in result.columns if model in c and ref in c]
        if metric == 'distance':
            values = {c[model]: math.dist(xyz(model, c[model]),
                                          xyz(ref, c[ref])) for c in cols}
        else:
            aligned = [c[ref] for c in cols]
            others = [i for i in result.residues[ref] if i not in aligned]
            ref_all = np.array([xyz(ref, i) for i in aligned + others])
            values = dict(zip(
                (c[model] for c in cols),
                lddt(ref_all[:len(cols)],
                     np.array([xyz(model, c[model]) for c in cols]),
                     ref_all)))
        values = {i: float(v) for i, v in values.items() if not np.isnan(v)}
        pair_scores[(model, ref)] = values
        return values

    def symmetric(a, b):
        '''
        {column number: score} of the pair (a, b), both ways for lDDT
        '''
        ab, ba = scores_of(a, b), scores_of(b, a)
        out = {}
        for k, c in enumerate(result.columns):
            if a in c and b in c:
                values = [v for v in (ab.get(c[a]), ba.get(c[b]))
                          if v is not None]
                if values:
                    out[k] = sum(values) / len(values)
        return out

    result.scores = {obj: {} for obj in names}
    if mode == 'reference':
        sums = {}
        for model in names[1:]:
            values = scores_of(model, reference)
            result.scores[model] = values
            for k, c in enumerate(result.columns):
                if model in c and reference in c and c[model] in values:
                    sums.setdefault(c[reference], []).append(
                        values[c[model]])
        result.scores[reference] = {i: sum(v) / len(v)
                                    for i, v in sums.items()}
    else:
        per_column = {}
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                for k, value in symmetric(a, b).items():
                    per_column.setdefault(k, []).append(value)
        for k, values in per_column.items():
            mean = sum(values) / len(values)
            for obj, index in result.columns[k].items():
                result.scores[obj][index] = mean

    # per structure: against the reference
    for obj in names:
        cols = [c for c in result.columns if obj in c and reference in c]
        values = list(result.scores[obj].values())
        stats = {'aligned': len(cols) if obj != reference else
                 len(result.residues[obj]),
                 'score': sum(values) / len(values) if values else None,
                 'rmsd': None, 'identity': None}
        if obj != reference and cols:
            d = [math.dist(xyz(obj, c[obj]), xyz(reference, c[reference]))
                 for c in cols]
            stats['rmsd'] = math.sqrt(sum(x * x for x in d) / len(d))
            same = sum(result.residues[obj][c[obj]][2] ==
                       result.residues[reference][c[reference]][2]
                       for c in cols)
            stats['identity'] = 100.0 * same / len(cols)
        result.stats[obj] = stats
    return result


def ramp_color(fraction):
    '''
    (r, g, b) from red (0, different) through yellow and green to blue
    (1, similar)
    '''
    fraction = min(1.0, max(0.0, fraction))
    stops = [(1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 0.8, 0.0),
             (0.0, 0.0, 1.0)]
    position = fraction * (len(stops) - 1)
    k = min(int(position), len(stops) - 2)
    t = position - k
    return tuple(a + (b - a) * t for a, b in zip(stops[k], stops[k + 1]))


def similarity_fraction(score, metric, low=None, high=None):
    '''
    0 (red) to 1 (blue): where a score is in its range
    '''
    red, blue = RANGES[metric]
    red = red if low is None else low
    blue = blue if high is None else high
    if blue == red:
        return 1.0
    return min(1.0, max(0.0, (score - red) / (blue - red)))


def apply(result, low=None, high=None, legend=True, _self=cmd):
    '''
    Color the structures by their scores; a color bar (legend)

    :param low: score shown red {default: RANGES}
    :param high: score shown blue
    '''
    names = []
    for k in range(RAMP_STEPS + 1):
        name = 'cmp_%03d' % k
        _self.set_color(name, list(ramp_color(k / RAMP_STEPS)))
        names.append(name)
    for obj in result.objects:
        _self.color(UNALIGNED_COLOR, '%%%s and polymer' % obj)
        groups = {}
        for index, score in result.scores[obj].items():
            k = round(RAMP_STEPS * similarity_fraction(
                score, result.metric, low, high))
            groups.setdefault(k, []).append(index)
        for k, indices in groups.items():
            for start in range(0, len(indices), 500):
                _self.color(names[k], 'byres (%%%s and index %s)' % (
                    obj, '+'.join(map(str, indices[start:start + 500]))))
    if legend:
        red, blue = RANGES[result.metric]
        red = red if low is None else low
        blue = blue if high is None else high
        legend_name = _self.get_legal_name('compare_legend')
        _self.delete(legend_name)
        lo, hi = sorted((red, blue))
        colors = RAMP if red < blue else list(reversed(RAMP))
        _self.ramp_new(legend_name, 'none', [lo, hi], colors, quiet=1)
        result.legend = legend_name


def select_dissimilar(result, threshold, name='dissimilar', _self=cmd):
    '''
    Residues scoring worse than the threshold (lDDT below, distance above)

    :return: number of residues
    '''
    parts = []
    count = 0
    for obj in result.objects:
        bad = [i for i, s in result.scores[obj].items()
               if (s < threshold if result.metric == 'lddt'
                   else s > threshold)]
        count += len(bad)
        for start in range(0, len(bad), 500):
            parts.append('(%%%s and index %s)' % (
                obj, '+'.join(map(str, bad[start:start + 500]))))
    _self.select(name, 'byres (%s)' % (' or '.join(parts) or 'none'),
                 enable=1)
    return count


def compare_structures(objects='', reference='', selection='polymer.protein',
                       method='align', metric='lddt', mode='reference',
                       low=None, high=None, file='', state=1, quiet=0, *,
                       _self=cmd):
    '''
DESCRIPTION

    "compare_structures" aligns closely related proteins on a reference
    and colors them by local structural similarity, from blue (similar)
    through green and yellow to red (different); grey: no partner.

USAGE

    compare_structures [ objects [, reference [, selection [, method
        [, metric [, mode ]]]]]]

ARGUMENTS

    objects = str: names or patterns {default: all molecular objects}

    reference = str {default: the first}

    selection = str: what to compare, e.g. "chain A" {default:
    polymer.protein}

    method = align (sequence-guided), super or cealign {default: align}

    metric = lddt (local distance difference test, 0-100, independent
    of the superposition) or distance (C-alpha, after superposition, A)
    {default: lddt}

    mode = reference (each to the reference) or all (mean over all
    pairs at each alignment position) {default: reference}

    low, high = float: scores shown red and blue {default: lDDT 70 and
    100; distance 3 and 0 A}

    file = str: save the scores per residue as CSV

EXAMPLE

    compare_structures wt mut1 mut2, reference=wt, selection=chain A

SEE ALSO

    align, super, cealign, extra_fit
    '''
    result = compute(objects, reference, selection, method, metric, mode,
                     state, _self=_self)
    low = None if low in (None, '') else float(low)
    high = None if high in (None, '') else float(high)
    apply(result, low, high, _self=_self)
    if file:
        result.save_csv(_self.exp_path(file))
    if not int(quiet):
        print(' compare_structures: %d structures on %s (%s), %d aligned '
              'positions, %s' % (len(result.objects), result.reference,
                                 method, len(result.columns),
                                 'lDDT' if result.metric == 'lddt' else
                                 'C-alpha distance'))
        print(result.table())
    return result


def extend(_self):
    _self.extend('compare_structures', compare_structures)
    _self.compare_structures = compare_structures
