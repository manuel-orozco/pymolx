'''
Sequence comparison (command "compare_sequences", Sequence Comparison
dialog) of related proteins: the chains of the structures are aligned by
sequence (Biopython's global aligner, BLOSUM62, gaps -10/-0.5, free end
gaps) to the reference's chains with the same chain ID (or one chain to
one chain), the pairwise alignments merged on the reference (a star
alignment: right for close relatives), and each residue colored like the
Structure Comparison: blue (identical), green (similar: BLOSUM62 > 0),
yellow (neutral: 0), red (different: < 0), grey (no partner).

Residues PyMOL knows without coordinates (e.g. PDB REMARK 465) are part
of the sequences (lower case in the alignment), but can't be colored.

Compared to the reference: each chain's residues by their pair with the
reference; the reference's by the mean over the others (a deletion inside
another chain counts as different, missing ends don't). All against all:
every residue by the mean over all pairs at its alignment position.
'''

from pymol import cmd

from .cluster import THREE_TO_ONE

MODES = ('reference', 'all')
# score per residue pair, 0-100: identical, similar, neutral, different
IDENTICAL, SIMILAR, NEUTRAL, DIFFERENT = 100.0, 200.0 / 3, 100.0 / 3, 0.0
CATEGORIES = [(IDENTICAL, 'identical'), (SIMILAR, 'similar (BLOSUM62 > 0)'),
              (NEUTRAL, 'neutral (BLOSUM62 0)'),
              (DIFFERENT, 'different (BLOSUM62 < 0)')]
# modified amino acids as their parents (phosphorylated, methylated,
# oxidized, ...)
MODIFIED = {'SEP': 'S', 'TPO': 'T', 'PTR': 'Y', 'HYP': 'P', 'MLY': 'K',
            'M3L': 'K', 'MLZ': 'K', 'ALY': 'K', 'KCX': 'K', 'LLP': 'K',
            'CSO': 'C', 'CSD': 'C', 'CME': 'C', 'CSX': 'C', 'OCS': 'C',
            'PCA': 'Q', 'NLE': 'L', 'ABA': 'A', 'FME': 'M', 'SAC': 'S',
            'CGU': 'E', 'HIC': 'H', 'DAH': 'F', 'TYS': 'Y', 'AGM': 'R'}
ONE_LETTER = dict(THREE_TO_ONE, SEC='U', PYL='O', ASX='B', GLX='Z',
                  **MODIFIED)
GAP_OPEN, GAP_EXTEND = -10.0, -0.5

_matrix = None


def blosum62():
    global _matrix
    if _matrix is None:
        from Bio.Align import substitution_matrices
        _matrix = substitution_matrices.load('BLOSUM62')
    return _matrix


def _scorable(aa):
    return {'U': 'C', 'O': 'K'}.get(aa, aa if aa in 'ARNDCQEGHILKMFPSTWYVBZX'
                                    else 'X')


def pair_score(a, b):
    '''
    0-100 for two residues (one-letter codes)
    '''
    if a == b and a != 'X':
        return IDENTICAL
    value = blosum62()[_scorable(a)][_scorable(b)]
    return SIMILAR if value > 0 else NEUTRAL if value == 0 else DIFFERENT


class Entry:
    '''
    One chain of one structure: its residues in order
    '''

    def __init__(self, obj, chain):
        self.obj = obj
        self.chain = chain
        self.residues = []   # (resi, resn, one letter, CA index or None)

    @property
    def name(self):
        return '%s/%s' % (self.obj, self.chain or '-')

    @property
    def sequence(self):
        return ''.join(r[2] for r in self.residues)


def _entries(obj, selection, _self):
    '''
    :return: [Entry] per chain, residues in order (with residues without
        coordinates)
    '''
    rows = []
    _self.iterate('%%%s and (%s) and polymer.protein and name CA' % (
        obj, selection), 'rows.append((chain, resi, resn, index))',
        space={'rows': rows})
    present = set()
    _self.iterate('%%%s and (%s) and polymer.protein and name CA and present'
                  % (obj, selection), 'present.add(index)',
                  space={'present': present})
    entries = {}
    seen = set()
    for chain, resi, resn, index in rows:
        if (chain, resi) in seen:
            continue  # (alternate locations)
        seen.add((chain, resi))
        entry = entries.setdefault(chain, Entry(obj, chain))
        entry.residues.append((resi, resn, ONE_LETTER.get(resn, 'X'),
                               index if index in present else None))
    return list(entries.values())


def _align_pair(reference, other):
    '''
    :return: {reference position: other position} for aligned residues
    '''
    from Bio import Align
    aligner = Align.PairwiseAligner(mode='global',
                                    substitution_matrix=blosum62(),
                                    open_gap_score=GAP_OPEN,
                                    extend_gap_score=GAP_EXTEND)
    aligner.open_end_gap_score = 0.0
    aligner.extend_end_gap_score = 0.0
    a = ''.join(_scorable(c) for c in reference)
    b = ''.join(_scorable(c) for c in other)
    alignment = aligner.align(a, b)[0]
    mapping = {}
    for (r0, r1), (o0, o1) in zip(*alignment.aligned):
        for k in range(r1 - r0):
            mapping[r0 + k] = o0 + k
    return mapping


class Group:
    '''
    The reference chain and the matching chains: a star alignment

    columns: [{entry name: position}]
    '''

    def __init__(self, reference, others):
        self.reference = reference
        self.entries = [reference] + others
        n = len(reference.residues)
        maps = {e.name: _align_pair(reference.sequence, e.sequence)
                for e in others}
        # each other chain: residues before reference position k (inserts)
        inserts = {}
        for e in others:
            mapping = maps[e.name]
            mapped = sorted(mapping.items())
            inserts[e.name] = blocks = [[] for _ in range(n + 1)]
            done = -1
            for ref_pos, pos in mapped:
                blocks[ref_pos].extend(range(done + 1, pos))
                done = pos
            blocks[n].extend(range(done + 1, len(e.residues)))
        self.columns = []
        for k in range(n + 1):
            width = max((len(inserts[e.name][k]) for e in others), default=0)
            for j in range(width):
                column = {}
                for e in others:
                    block = inserts[e.name][k]
                    # (right-aligned before the reference position, so a
                    # leading extension ends where the reference starts)
                    offset = j - (width - len(block))
                    if offset >= 0:
                        column[e.name] = block[offset]
                self.columns.append(column)
            if k < n:
                column = {reference.name: k}
                for e in others:
                    if k in maps[e.name]:
                        column[e.name] = maps[e.name][k]
                self.columns.append(column)
        # each chain's first and last aligned column (end gaps outside)
        self.span = {}
        for e in self.entries:
            ks = [k for k, c in enumerate(self.columns) if e.name in c]
            self.span[e.name] = (min(ks), max(ks)) if ks else (0, -1)

    def letter(self, entry, column):
        pos = self.columns[column].get(entry.name)
        return None if pos is None else entry.residues[pos][2]

    def internal_gap(self, entry, column):
        first, last = self.span[entry.name]
        return entry.name not in self.columns[column] and \
            first < column < last


class Result:
    '''
    groups, scores {entry name: {position: score}}, stats per entry,
    identity matrix
    '''

    def __init__(self):
        self.groups = []
        self.reference = ''
        self.mode = 'reference'
        self.scores = {}
        self.stats = {}
        self.alignment = ''

    @property
    def entries(self):
        return [e for g in self.groups for e in g.entries]

    def identity_matrix(self):
        '''
        :return: (names, {(name, name): percent identity or None}) over
            the aligned positions of each pair (same chain group)
        '''
        names, values = [], {}
        for group in self.groups:
            for a in group.entries:
                names.append(a.name)
                for b in group.entries:
                    same = aligned = 0
                    for k, c in enumerate(group.columns):
                        if a.name in c and b.name in c:
                            aligned += 1
                            same += group.letter(a, k) == group.letter(b, k)
                    values[(a.name, b.name)] = 100.0 * same / aligned \
                        if aligned else None
        return names, values

    def table(self):
        lines = [' %-28s %6s %8s %9s %10s %5s' % (
            'Chain', 'Length', 'Aligned', 'Identity', 'Similarity', 'Gaps')]
        for e in self.entries:
            s = self.stats[e.name]
            lines.append(' %-28s %6d %8s %8s %10s %5s' % (
                e.name + (' (reference)' if s['reference'] else ''),
                s['length'], s['aligned'],
                '-' if s['identity'] is None else '%.1f%%' % s['identity'],
                '-' if s['similarity'] is None else
                '%.1f%%' % s['similarity'],
                '-' if s['gaps'] is None else s['gaps']))
        return '\n'.join(lines)

    def fasta(self):
        '''
        The alignments, one block per chain group (residues without
        coordinates in lower case)
        '''
        lines = []
        for group in self.groups:
            for e in group.entries:
                text = []
                for k in range(len(group.columns)):
                    pos = group.columns[k].get(e.name)
                    if pos is None:
                        text.append('-')
                    else:
                        resi, resn, letter, index = e.residues[pos]
                        text.append(letter if index is not None
                                    else letter.lower())
                lines.append('>' + e.name)
                text = ''.join(text)
                lines += [text[i:i + 60] for i in range(0, len(text), 60)]
        return '\n'.join(lines) + '\n'

    def save_csv(self, path):
        import csv
        with open(path, 'w', newline='') as handle:
            out = csv.writer(handle)
            out.writerow(['chain', 'resi', 'resn', 'column', 'score'])
            for group in self.groups:
                for e in group.entries:
                    column_of = {c[e.name]: k for k, c in
                                 enumerate(group.columns) if e.name in c}
                    for pos, (resi, resn, _, _) in enumerate(e.residues):
                        score = self.scores[e.name].get(pos)
                        out.writerow([e.name, resi, resn,
                                      column_of.get(pos, ''),
                                      '' if score is None else
                                      '%.1f' % score])


def _objects(objects, _self):
    from .compare import _objects as objects_of
    return objects_of(objects, _self)


def compute(objects='', reference='', selection='polymer.protein',
            mode='reference', _self=cmd):
    '''
    :return: Result
    '''
    mode = str(mode).strip().lower()
    if mode not in MODES:
        raise ValueError('mode must be one of: ' + ', '.join(MODES))
    names = _objects(objects, _self)
    if len(names) < 2:
        raise ValueError('need at least two structures (objects)')
    reference = str(reference).strip() or names[0]
    if reference not in names:
        raise ValueError('reference "%s" is not among the structures' %
                         reference)
    selection = str(selection).strip() or 'polymer.protein'
    chains = {n: _entries(n, selection, _self) for n in names}
    if not chains[reference]:
        raise ValueError('no protein chains in "%s and (%s)"' % (
            reference, selection))

    result = Result()
    result.reference, result.mode = reference, mode
    single = len(chains[reference]) == 1
    for ref_entry in chains[reference]:
        others = []
        for n in names:
            if n == reference:
                continue
            same = [e for e in chains[n] if e.chain == ref_entry.chain]
            if not same and single and len(chains[n]) == 1:
                same = chains[n]
            others += same[:1]
        if others:
            result.groups.append(Group(ref_entry, others))
    if not result.groups:
        raise ValueError('no chains to compare with the reference\'s '
                         '(matched by chain ID)')

    for group in result.groups:
        ref = group.reference
        for e in group.entries:
            result.scores[e.name] = {}
        for k, column in enumerate(group.columns):
            if mode == 'reference':
                a = group.letter(ref, k)
                values = []
                for e in group.entries[1:]:
                    b = group.letter(e, k)
                    if a is not None and b is not None:
                        score = pair_score(a, b)
                        result.scores[e.name][column[e.name]] = score
                        values.append(score)
                    elif a is not None and group.internal_gap(e, k):
                        values.append(DIFFERENT)
                if a is not None and values:
                    result.scores[ref.name][column[ref.name]] = \
                        sum(values) / len(values)
            else:
                values = []
                entries = group.entries
                for i, x in enumerate(entries):
                    for y in entries[i + 1:]:
                        a, b = group.letter(x, k), group.letter(y, k)
                        if a is not None and b is not None:
                            values.append(pair_score(a, b))
                        elif (a is not None and group.internal_gap(y, k)) or \
                                (b is not None and group.internal_gap(x, k)):
                            values.append(DIFFERENT)
                if values:
                    mean = sum(values) / len(values)
                    for name, pos in column.items():
                        result.scores[name][pos] = mean

        for e in group.entries:
            stats = {'reference': e is ref, 'length': len(e.residues),
                     'aligned': len(e.residues) if e is ref else 0,
                     'identity': None, 'similarity': None, 'gaps': None}
            if e is not ref:
                aligned = same = similar = gaps = 0
                for k, column in enumerate(group.columns):
                    a, b = group.letter(ref, k), group.letter(e, k)
                    if a is not None and b is not None:
                        aligned += 1
                        score = pair_score(a, b)
                        same += score == IDENTICAL
                        similar += score >= SIMILAR
                    elif (a is None) != (b is None):
                        first, last = group.span[e.name]
                        rf, rl = group.span[ref.name]
                        if max(first, rf) < k < min(last, rl):
                            gaps += 1
                stats.update(aligned=aligned, gaps=gaps,
                             identity=100.0 * same / aligned
                             if aligned else None,
                             similarity=100.0 * similar / aligned
                             if aligned else None)
            result.stats[e.name] = stats
    return result


def apply(result, superpose=False, alignment='seq_aln', _self=cmd):
    '''
    Color the residues (with coordinates) by their scores, make an
    alignment object, superimpose on the reference if asked
    '''
    from .compare import RAMP_STEPS, UNALIGNED_COLOR, ramp_color
    names = []
    for k in range(RAMP_STEPS + 1):
        name = 'cmp_%03d' % k
        _self.set_color(name, list(ramp_color(k / RAMP_STEPS)))
        names.append(name)
    for group in result.groups:
        for e in group.entries:
            present = [r[3] for r in e.residues if r[3] is not None]
            for start in range(0, len(present), 500):
                _self.color(UNALIGNED_COLOR, 'byres (%%%s and index %s)' % (
                    e.obj, '+'.join(map(str, present[start:start + 500]))))
            by_color = {}
            for pos, score in result.scores[e.name].items():
                index = e.residues[pos][3]
                if index is not None:
                    by_color.setdefault(
                        round(RAMP_STEPS * score / 100.0), []).append(index)
            for k, indices in by_color.items():
                for start in range(0, len(indices), 500):
                    _self.color(names[k], 'byres (%%%s and index %s)' % (
                        e.obj, '+'.join(map(str,
                                            indices[start:start + 500]))))
    if superpose:
        _superpose(result, _self)
    # an alignment object (residues with coordinates), made after the
    # superposition (its lines would join the old places) and hidden
    raw = []
    for group in result.groups:
        for column in group.columns:
            atoms = []
            for e in group.entries:
                pos = column.get(e.name)
                if pos is not None and e.residues[pos][3] is not None:
                    atoms.append((e.obj, e.residues[pos][3]))
            if len(atoms) > 1:
                raw.append(atoms)
    result.alignment = _self.get_legal_name(alignment)
    _self.delete(result.alignment)
    if raw:
        _self.set_raw_alignment(result.alignment, raw)
        _self.disable(result.alignment)


def _superpose(result, _self):
    '''
    Each structure on the reference, by its aligned C-alphas
    '''
    import numpy as np
    from .cluster import kabsch
    pairs = {}
    for group in result.groups:
        ref = group.reference
        for column in group.columns:
            a = column.get(ref.name)
            if a is None or ref.residues[a][3] is None:
                continue
            for e in group.entries[1:]:
                b = column.get(e.name)
                if b is not None and e.residues[b][3] is not None:
                    pairs.setdefault(e.obj, []).append(
                        (ref.residues[a][3], e.residues[b][3]))
    ref_obj = result.groups[0].reference.obj
    ref_xyz = {}
    _self.iterate_state(1, '%' + ref_obj, 'xyz[index] = (x, y, z)',
                        space={'xyz': ref_xyz})
    for obj, matched in pairs.items():
        if len(matched) < 3:
            continue
        xyz = {}
        _self.iterate_state(1, '%' + obj, 'xyz[index] = (x, y, z)',
                            space={'xyz': xyz})
        rotation, translation = kabsch(
            np.array([xyz[b] for _, b in matched]),
            np.array([ref_xyz[a] for a, _ in matched]))
        coords = _self.get_coords('%' + obj, 1)
        _self.load_coords(coords @ rotation.T + translation, '%' + obj,
                          state=1)


def compare_sequences(objects='', reference='', selection='polymer.protein',
                      mode='reference', superpose=0, file='', quiet=0, *,
                      _self=cmd):
    '''
DESCRIPTION

    "compare_sequences" aligns the chains of related proteins by sequence
    (BLOSUM62) on a reference and colors the residues: blue identical,
    green similar, yellow neutral, red different, grey no partner.

USAGE

    compare_sequences [ objects [, reference [, selection [, mode
        [, superpose [, file ]]]]]]

ARGUMENTS

    objects = str: names or patterns {default: all molecular objects}

    reference = str {default: the first}

    selection = str: what to compare, e.g. "chain H" {default:
    polymer.protein}; chains are matched by chain ID

    mode = reference (each to the reference) or all (mean over all
    pairs at each position) {default: reference}

    superpose = 0/1: superimpose on the reference by the aligned residues
    {default: 0}

    file = str: save the alignment as FASTA

EXAMPLE

    compare_sequences wt mut*, reference=wt

SEE ALSO

    compare_structures, align
    '''
    result = compute(objects, reference, selection, mode, _self=_self)
    apply(result, int(superpose), _self=_self)
    if file:
        with open(_self.exp_path(file), 'w') as handle:
            handle.write(result.fasta())
    if not int(quiet):
        print(' compare_sequences: %d chains on %s, %s' % (
            len(result.entries), result.reference,
            'each to the reference' if result.mode == 'reference' else
            'all against all'))
        print(result.table())
    return result


def extend(_self):
    _self.extend('compare_sequences', compare_sequences)
    _self.compare_sequences = compare_sequences
