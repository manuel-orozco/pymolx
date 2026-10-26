'''
Clustering of structures by RMSD (command "cluster_rmsd", RMSD
Clustering dialog), made for sets of predicted complexes (e.g. antibody-
antigen models from different seeds): each pair of models is superimposed
on one part (the antigen) and the RMSD is measured on another part (the
antibody) without refitting, like the ligand RMSD of docking poses: small
when the antibody binds the same epitope in the same orientation.

Atoms are matched by chain, residue and name; only atoms all models have
are used. Clusters: hierarchical (average, complete or single linkage)
cut at an RMSD, or GROMOS (Daura et al. 1999: the model with the most
neighbors within the cutoff and its neighbors form a cluster, repeated).
Each cluster's representative is the member with the lowest mean RMSD to
the other members; its best model the one with the highest confidence
score (ipTM by default), read from the predictor's JSON next to the model
file (AlphaFold 3 "<name>_summary_confidences.json", Boltz
"confidence_<name>.json", ColabFold "<name with scores>.json"), or found
in a given folder by the object's name, or from a CSV file with a row per
model (its scores take precedence; any of its columns can pick the best
model).
'''

import csv

import json
import os

import fnmatch
import re

from pymol import cmd

METHODS = ('average', 'complete', 'single', 'gromos')
SCORES = {'iptm': 'ipTM', 'ptm': 'pTM', 'ranking': 'ranking score',
          'plddt': 'pLDDT'}
# CSV column names meaning the scores above (lower case)
SCORE_COLUMNS = {
    'iptm': 'iptm', 'ptm': 'ptm', 'plddt': 'plddt', 'mean_plddt': 'plddt',
    'ranking_score': 'ranking', 'ranking': 'ranking',
    'confidence_score': 'ranking',
}
# CSV columns naming the model
NAME_COLUMNS = ('model', 'name', 'object', 'file', 'filename', 'path',
                'structure', 'id')
# scores where lower is better (by their names)
LOWER_IS_BETTER = ('rank', 'energy', 'rmsd', 'dg', 'kcal', 'clash')
ATOMS = {'ca': 'name CA', 'backbone': 'name N+CA+C+O',
         'heavy': 'not hydro'}
# conserved motifs of antibody variable domains (heavy, light, VHH, scFv;
# human and mouse): the framework 1 Cys 9-16 residues before the framework
# 2 Trp (W-V/I/Y/F-x-Q/K/R, e.g. WVRQ, WIRK, WYQQ), later the framework 4
# W-G-x-G (heavy) or F-G-x-G (light)
VDOMAIN_MOTIFS = (r'C.{9,16}W[VIFLYM].[QKRE]', r'[WF]G.G')


THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V', 'MSE': 'M'}


def _sequences(obj, _self):
    '''
    :return: {chain: one-letter sequence}
    '''
    seq = {}
    _self.iterate('%%%s and polymer.protein and name CA' % obj,
                  '_ = seq.setdefault(chain, []).append(resn)',
                  space={'seq': seq})
    return {chain: ''.join(THREE_TO_ONE.get(r, 'X') for r in resns)
            for chain, resns in seq.items()}


def antibody_chains(obj, _self=cmd):
    '''
    Chains with an antibody variable domain (heavy, light, VHH, scFv), by
    their conserved framework motifs

    :return: sorted chain IDs
    '''
    chains = []
    for chain, seq in _sequences(obj, _self).items():
        first = re.search(VDOMAIN_MOTIFS[0], seq)
        if first and re.search(VDOMAIN_MOTIFS[1], seq[first.end():]):
            chains.append(chain)
    return sorted(chains)


def default_selections(obj, _self=cmd):
    '''
    :return: (fit selection, RMSD selection): the antigen and the
        antibody, or the whole structure twice if no antibody is found
    '''
    antibody = antibody_chains(obj, _self)
    chains = sorted(_sequences(obj, _self))
    antigen = [c for c in chains if c not in antibody]
    if antibody and antigen:
        return ('chain %s' % '+'.join(antigen),
                'chain %s' % '+'.join(antibody))
    return 'polymer', 'polymer'


def _stem(path):
    name = os.path.basename(path)
    for ext in ('.gz', '.bz2'):
        if name.endswith(ext):
            name = name[:-len(ext)]
    return os.path.splitext(name)[0]


def _score_names(stem):
    '''
    Score file names a predictor writes for a model file "stem"
    '''
    names = []
    if stem.endswith('_model'):                       # AlphaFold 3
        names.append(stem[:-len('_model')] + '_summary_confidences.json')
    if stem == 'model':                               # AlphaFold 3 sample
        names.append('summary_confidences.json')
    names.append('confidence_%s.json' % stem)         # Boltz
    for kind in ('_unrelaxed_', '_relaxed_'):         # ColabFold
        if kind in stem:
            names.append(stem.replace(kind, '_scores_') + '.json')
    names += [stem + '_summary_confidences.json', stem + '_scores.json']
    return names


def read_scores(path):
    '''
    :return: {'iptm', 'ptm', 'ranking': float} found in a score JSON
    '''
    with open(path) as handle:
        data = json.load(handle)
    scores = {}
    for key, names in (('iptm', ('iptm', 'protein_iptm')),
                       ('ptm', ('ptm', 'protein_ptm')),
                       ('ranking', ('ranking_score', 'confidence_score'))):
        for name in names:
            value = data.get(name)
            if isinstance(value, list) and len(value) == 1:
                value = value[0]
            if isinstance(value, (int, float)):
                scores[key] = float(value)
                break
    return scores


def find_scores(obj, folder=''):
    '''
    The confidence scores of a predicted model: from the JSON next to the
    file it was loaded from, or found in "folder" by the object's name

    :return: ({'iptm', 'ptm', 'ranking': float}, path) or ({}, None)
    '''
    from . import files
    places = []
    path = files.source(obj)
    if path:
        places.append((os.path.dirname(path), _score_names(_stem(path))))
    if folder:
        names = _score_names(obj)
        folder = os.path.expanduser(folder)
        for root, _, found in os.walk(folder):
            places.append((root, [n for n in names if n in found]))
    for directory, names in places:
        for name in names:
            candidate = os.path.join(directory, name)
            if os.path.isfile(candidate):
                try:
                    return read_scores(candidate), candidate
                except (OSError, ValueError):
                    continue
    return {}, None


def read_scores_csv(path):
    '''
    Scores per model from a CSV file: a header row; the model in a column
    named model, name, object, file, ... (else the first column), or in
    "seed" and "sample" columns (AlphaFold 3's <job>_ranking_scores.csv:
    models <job>_seed-<seed>_sample-<sample>_model); numeric columns as
    scores (iptm, ptm, ranking_score, ... under their usual names, others
    by their lower-case header)

    :return: {model name, file or "seed-1_sample-0": {score: float}}
    '''
    with open(path, newline='') as handle:
        text = handle.read()
    try:
        dialect = csv.Sniffer().sniff(text.splitlines()[0], ',;\t')
    except (csv.Error, IndexError):
        dialect = csv.excel
    rows = list(csv.reader(text.splitlines(), dialect))
    if len(rows) < 2:
        raise ValueError('%s: a header row and a row per model' % path)
    header = [h.strip().lower() for h in rows[0]]
    named = [header.index(n) for n in NAME_COLUMNS if n in header]
    seed_sample = not named and {'seed', 'sample'} <= set(header)
    key = named[0] if named else 0
    keys = {header.index('seed'), header.index('sample')} if seed_sample \
        else {key}
    scores = {}
    for row in rows[1:]:
        if len(row) <= max(keys):
            continue
        if seed_sample:
            name = 'seed-%s_sample-%s' % (row[header.index('seed')].strip(),
                                          row[header.index('sample')].strip())
        else:
            name = row[key].strip()
        if not name:
            continue
        values = {}
        for column, (title, cell) in enumerate(zip(header, row)):
            if column in keys:
                continue
            try:
                value = float(cell)
            except ValueError:
                continue
            values[SCORE_COLUMNS.get(title, title)] = value
        scores[name] = values
    return scores


def read_scores_csvs(pattern, _self=cmd):
    '''
    Scores from one or more CSV files (a path, or a pattern like
    "results/*_ranking_scores.csv"); each AlphaFold 3 file's rows only
    for the models of its job (<job>_ranking_scores.csv)

    :return: {(job prefix or '', row name): {score: float}}
    '''
    import glob
    paths = sorted(glob.glob(_self.exp_path(pattern))) or \
        [_self.exp_path(pattern)]
    table = {}
    for path in paths:
        job = os.path.basename(path)
        job = job[:-len('_ranking_scores.csv')] \
            if job.endswith('_ranking_scores.csv') else ''
        for name, values in read_scores_csv(path).items():
            table[(job, name)] = values
    return table


def match_scores(models, table):
    '''
    Rows of a scores table for the models: by object name, by file name
    without extension, or by a part of the name only one model has, as a
    whole word (e.g. "seed-10_sample-2" in "job_seed-10_sample-2_model",
    not in "job_seed-10_sample-21_model"), among the models of the row's
    job if it has one

    :param table: {(job, row name): scores} (read_scores_csvs) or
        {row name: scores}
    :return: ({model: scores}, [row names not matched])
    '''
    found = {}
    unmatched = []
    for entry, values in table.items():
        job, name = entry if isinstance(entry, tuple) else ('', entry)
        candidates = [m for m in models if m.startswith(job)] if job \
            else models
        stem = _stem(name)
        if name in candidates:
            found[name] = values
        elif stem in candidates:
            found[stem] = values
        else:
            word = re.compile(r'(?<![A-Za-z0-9])%s(?![A-Za-z0-9])' %
                              re.escape(stem)) if stem else None
            hits = [m for m in candidates if word and word.search(m)]
            if len(hits) == 1:
                found[hits[0]] = values
            else:
                unmatched.append(name)
    return found, unmatched


def higher_is_better(score):
    name = str(score).lower()
    if name in SCORES or 'score' in name and 'rank' in name:
        return True
    return not any(word in name for word in LOWER_IS_BETTER)


def _objects(objects, _self):
    names = _self.get_object_list('all') or []
    molecules = [n for n in names
                 if _self.get_type(n) == 'object:molecule']
    if not str(objects).strip():
        return molecules
    found = []
    for pattern in str(objects).split():
        found += [n for n in molecules
                  if fnmatch.fnmatchcase(n, pattern) and n not in found]
    return found


def _atoms(obj, selection, state, _self):
    '''
    :return: {(chain, resi, name): (x, y, z)}
    '''
    atoms = {}
    _self.iterate_state(state, '%%%s and (%s) and present' % (obj, selection),
                        'atoms[(chain, resi, name)] = (x, y, z)',
                        space={'atoms': atoms})
    return atoms


def kabsch(mobile, target):
    '''
    :return: (rotation, translation) moving "mobile" onto "target"
    '''
    import numpy as np
    pc, tc = mobile.mean(axis=0), target.mean(axis=0)
    h = (mobile - pc).T @ (target - tc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rotation, tc - rotation @ pc


def rmsd_matrix(fit, measure):
    '''
    :param fit: [array (n, 3)] per model, atoms to superimpose on
    :param measure: [array (m, 3)] per model, atoms to measure
    :return: symmetric matrix (models x models), Angstrom
    '''
    import numpy as np
    n = len(fit)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            rotation, translation = kabsch(fit[j], fit[i])
            moved = measure[j] @ rotation.T + translation
            d = np.sqrt(((moved - measure[i]) ** 2).sum(axis=1).mean())
            matrix[i, j] = matrix[j, i] = d
    return matrix


def gromos(matrix, cutoff):
    '''
    :return: cluster number (0-based, by size) per model
    '''
    import numpy as np
    n = len(matrix)
    left = set(range(n))
    labels = [-1] * n
    cluster = 0
    while left:
        members = sorted(left)
        counts = [sum(matrix[i, j] <= cutoff for j in members)
                  for i in members]
        center = members[int(np.argmax(counts))]
        for j in members:
            if matrix[center, j] <= cutoff:
                labels[j] = cluster
                left.discard(j)
        cluster += 1
    return labels


def hierarchical(matrix, cutoff, method='average'):
    '''
    :return: cluster number (0-based, by size) per model
    '''
    import numpy as np
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    if len(matrix) < 2:
        return [0] * len(matrix)
    tree = linkage(squareform(matrix, checks=False), method=method)
    labels = fcluster(tree, t=cutoff, criterion='distance')
    # renumber by size (largest first), ties by first member
    order = sorted(set(labels), key=lambda c: (-np.sum(labels == c),
                                               list(labels).index(c)))
    rank = {c: k for k, c in enumerate(order)}
    return [rank[c] for c in labels]


class Result:
    '''
    models, matrix, labels and per-cluster summaries
    '''

    def __init__(self, models, matrix, labels, scores, fit_atoms,
                 rmsd_atoms, confidences=None, best_by='iptm', higher=None):
        self.models = models
        # per model {'iptm', 'ptm', 'ranking', 'plddt'} (those known)
        self.confidences = confidences or [{} for _ in models]
        for conf, plddt in zip(self.confidences, scores):
            if plddt is not None:
                conf.setdefault('plddt', plddt)
        self.best_by = best_by
        # None: by the score's name (e.g. lower is better for a rank)
        self.higher = higher
        self.matrix = matrix
        self.labels = labels
        self.scores = scores            # mean pLDDT/B per model (or None)
        self.fit_atoms = fit_atoms      # number of atoms
        self.rmsd_atoms = rmsd_atoms
        self.clusters = []
        for c in sorted(set(labels)):
            members = [i for i, label in enumerate(labels) if label == c]
            if len(members) > 1:
                mean = [sum(matrix[i, j] for j in members if j != i) /
                        (len(members) - 1) for i in members]
                representative = members[mean.index(min(mean))]
                spread = sum(matrix[i, j] for i in members
                             for j in members if i < j) / (
                    len(members) * (len(members) - 1) / 2)
            else:
                representative, spread = members[0], 0.0
            known = [scores[i] for i in members if scores[i] is not None]
            best, best_score, key = self._best(members)
            self.clusters.append({
                'cluster': c + 1,
                'size': len(members),
                'representative': models[representative],
                'members': [models[i] for i in members],
                'mean_rmsd': spread,
                'score': sum(known) / len(known) if known else None,
                'best': models[best] if best is not None else None,
                'best_score': best_score,
                'best_by': key,
            })

    def _best(self, members):
        '''
        The member with the highest score: best_by, else the first score
        all members have (ipTM, ranking, pTM, pLDDT)

        :return: (index, score, score name) or (None, None, None)
        '''
        for key in [self.best_by] + [k for k in ('iptm', 'ranking', 'ptm',
                                                 'plddt')
                                     if k != self.best_by]:
            values = [self.confidences[i].get(key) for i in members]
            if all(v is not None for v in values):
                higher = higher_is_better(key) if (
                    self.higher is None or key != self.best_by) \
                    else self.higher
                pick = max if higher else min
                k = pick(range(len(members)), key=lambda n: values[n])
                return members[k], values[k], key
        return None, None, None

    def table(self):
        width = max([len('Best')] + [len(format_best(c))
                                     for c in self.clusters])
        lines = [' %-7s %4s  %-9s %-6s  %-*s  %s' % (
            'Cluster', 'Size', 'Mean RMSD', 'pLDDT', width, 'Best',
            'Representative')]
        for c in self.clusters:
            lines.append(' %-7d %4d  %7.2f A  %-6s  %-*s  %s' % (
                c['cluster'], c['size'], c['mean_rmsd'],
                '%.1f' % c['score'] if c['score'] is not None else '-',
                width, format_best(c), c['representative']))
        for c in self.clusters:
            if c['best']:
                lines.append(' cluster %d best: %s' % (c['cluster'],
                                                      c['best']))
        return '\n'.join(lines)

    def save_csv(self, path):
        import csv
        with open(path, 'w', newline='') as handle:
            out = csv.writer(handle)
            extra = [] if self.best_by in SCORES else [self.best_by]
            out.writerow(['model', 'cluster', 'representative', 'best',
                          'pLDDT', 'ipTM', 'pTM', 'ranking'] + extra +
                         self.models)
            reps = {m: c['representative'] for c in self.clusters
                    for m in c['members']}
            best = {c['best'] for c in self.clusters}

            def value(conf, key, digits):
                v = conf.get(key)
                return '' if v is None else '%.*f' % (digits, v)

            for i, model in enumerate(self.models):
                conf = self.confidences[i]
                out.writerow([model, self.labels[i] + 1, reps[model] == model,
                              model in best, value(conf, 'plddt', 2),
                              value(conf, 'iptm', 3), value(conf, 'ptm', 3),
                              value(conf, 'ranking', 3)] +
                             [value(conf, k, 3) for k in extra] +
                             ['%.3f' % v for v in self.matrix[i]])


def short_names(names):
    '''
    The parts of the names that differ (common prefix and suffix cut at
    "_"), e.g. "seed-10_sample-2" for "job_seed-10_sample-2_model"

    :return: {name: short name}
    '''
    import os.path
    if len(names) < 2:
        return {n: n for n in names}
    prefix = os.path.commonprefix(names)
    suffix = os.path.commonprefix([n[::-1] for n in names])[::-1]
    cut = prefix.rfind('_') + 1
    end = suffix.find('_') if '_' in suffix else len(suffix)
    shorts = {}
    for n in names:
        short = n[cut:len(n) - (len(suffix) - end)]
        shorts[n] = short or n
    if len(set(shorts.values())) < len(names):
        return {n: n for n in names}
    return shorts


def format_best(c):
    '''
    "model_3 (ipTM 0.86)" for a cluster, or "-"
    '''
    if not c.get('best'):
        return '-'
    return '%s (%s %s)' % (c['best'], score_label(c['best_by']),
                           format_value(c['best_score']))


def score_label(key):
    return SCORES.get(key, key)


def format_value(value):
    # (AlphaFold 3 ranking scores of a job often differ in the 3rd digit)
    return '%.3f' % value if abs(value) < 10 else '%.1f' % value


def compute(objects='', fit='', rmsd='', atoms='ca', cutoff=3.0,
            method='average', state=1, best_by='iptm', scores='',
            scores_csv='', higher=None, _self=cmd):
    '''
    The clustering, without changing the session.

    :return: Result
    '''
    import numpy as np
    best_by = str(best_by).strip().lower()
    table = {}
    if str(scores_csv).strip():
        table = read_scores_csvs(str(scores_csv).strip(), _self=_self)
    columns = {k for values in table.values() for k in values}
    if best_by not in SCORES and best_by not in columns:
        raise ValueError('best_by must be one of: ' + ', '.join(
            list(SCORES) + sorted(columns - set(SCORES))))
    method = str(method).strip().lower()
    if method not in METHODS:
        raise ValueError('method must be one of: ' + ', '.join(METHODS))
    atoms = str(atoms).strip().lower()
    if atoms not in ATOMS:
        raise ValueError('atoms must be one of: ' + ', '.join(ATOMS))
    models = _objects(objects, _self)
    if len(models) < 2:
        raise ValueError('need at least two structures (objects)')
    if not str(fit).strip() or not str(rmsd).strip():
        auto_fit, auto_rmsd = default_selections(models[0], _self)
        fit = str(fit).strip() or auto_fit
        rmsd = str(rmsd).strip() or auto_rmsd
    state = int(state)

    fit_sel = '(%s) and %s' % (fit, ATOMS[atoms])
    rmsd_sel = '(%s) and %s' % (rmsd, ATOMS[atoms])
    fit_sets = [_atoms(m, fit_sel, state, _self) for m in models]
    rmsd_sets = [_atoms(m, rmsd_sel, state, _self) for m in models]
    fit_keys = sorted(set.intersection(*(set(a) for a in fit_sets)))
    rmsd_keys = sorted(set.intersection(*(set(a) for a in rmsd_sets)))
    if len(fit_keys) < 3:
        raise ValueError('fewer than 3 atoms in "%s" common to all '
                         'structures' % fit_sel)
    if not rmsd_keys:
        raise ValueError('no atoms in "%s" common to all structures' %
                         rmsd_sel)
    fit_xyz = [np.array([a[k] for k in fit_keys]) for a in fit_sets]
    rmsd_xyz = [np.array([a[k] for k in rmsd_keys]) for a in rmsd_sets]
    matrix = rmsd_matrix(fit_xyz, rmsd_xyz)
    cutoff = float(cutoff)
    labels = gromos(matrix, cutoff) if method == 'gromos' else \
        hierarchical(matrix, cutoff, method)

    # model confidence: mean B-factor of the measured part, if it looks
    # like pLDDT (predicted models store it there)
    plddt = []
    for m in models:
        values = []
        _self.iterate_state(state, '%%%s and (%s) and name CA' % (m, rmsd),
                            'values.append(b)', space={'values': values})
        mean = sum(values) / len(values) if values else None
        plddt.append(mean if mean is not None and 0 < mean <= 100 else None)
    confidences = [find_scores(m, scores)[0] for m in models]
    unmatched = []
    if table:
        matched, unmatched = match_scores(models, table)
        for m, conf in zip(models, confidences):
            conf.update(matched.get(m, {}))  # (the CSV first)
    result = Result(models, matrix, labels, plddt, len(fit_keys),
                    len(rmsd_keys), confidences, best_by, higher)
    result.unmatched_scores = unmatched
    result.fit, result.rmsd, result.method, result.cutoff = \
        fit, rmsd, method, cutoff
    return result


CLUSTER_COLORS = ['tv_blue', 'tv_orange', 'forest', 'raspberry', 'teal',
                  'olive', 'violetpurple', 'deepsalmon', 'marine', 'limon']


def apply(result, superpose=True, group=True, color=True, prefix='cluster',
          state=1, _self=cmd):
    '''
    Show the clustering: all models superimposed on the fit part of the
    largest cluster's representative, one group per cluster (the
    representative first), the measured part colored by cluster
    '''
    import numpy as np
    reference = result.clusters[0]['representative']
    fit_sel = '(%s) and %s' % (result.fit, 'name CA')
    if superpose:
        target = _atoms(reference, fit_sel, state, _self)
        for model in result.models:
            if model == reference:
                continue
            mobile = _atoms(model, fit_sel, state, _self)
            keys = sorted(set(target) & set(mobile))
            if len(keys) < 3:
                continue
            rotation, translation = kabsch(
                np.array([mobile[k] for k in keys]),
                np.array([target[k] for k in keys]))
            xyz = _self.get_coords('%' + model, state)
            _self.load_coords(xyz @ rotation.T + translation, '%' + model,
                              state=state)
    if group:
        # groups of an earlier clustering: emptied, then removed (deleting
        # a group would delete the models in it)
        pattern = re.compile(r'%s_\d+$' % re.escape(prefix))
        for name in _self.get_names('group_objects') or []:
            members = _self.get_object_list('(%s)' % name) or []
            if pattern.match(name) and set(members) <= set(result.models):
                _self.ungroup(name)
                _self.delete(name)
    for c in result.clusters:
        if group:
            name = _self.get_legal_name('%s_%d' % (prefix, c['cluster']))
            first = [m for m in (c['best'], c['representative']) if m]
            members = list(dict.fromkeys(first)) + [
                m for m in c['members'] if m not in first]
            _self.group(name, ' '.join(members))
            _self.order(' '.join(members), location='current')
        if color:
            colour = CLUSTER_COLORS[(c['cluster'] - 1) % len(CLUSTER_COLORS)]
            for m in c['members']:
                _self.color(colour, '%%%s and (%s) and elem C' % (
                    m, result.rmsd))
    if group:
        _self.order(' '.join('%s_%d' % (prefix, c['cluster'])
                             for c in result.clusters), location='current')


def cluster_rmsd(objects='', fit='', rmsd='', atoms='ca', cutoff=3.0,
                 method='average', superpose=1, group=1, color=1,
                 file='', best_by='iptm', scores='', scores_csv='',
                 higher='auto', state=1, quiet=0, *, _self=cmd):
    '''
DESCRIPTION

    "cluster_rmsd" clusters structures (objects) by RMSD, e.g. predicted
    antibody-antigen complexes: each pair is superimposed on "fit" (the
    antigen) and the RMSD is measured on "rmsd" (the antibody), without
    refitting. Atoms are matched by chain, residue and name.

USAGE

    cluster_rmsd [ objects [, fit [, rmsd [, atoms [, cutoff
        [, method ]]]]]]

ARGUMENTS

    objects = str: object names or patterns (e.g. "model_*") {default:
    all molecular objects}

    fit, rmsd = str: selections to superimpose on and to measure
    {default: antigen chains and antibody chains, found by their
    sequences; else all polymer atoms}

    atoms = ca, backbone or heavy {default: ca}

    cutoff = float: RMSD, Angstrom {default: 3.0}

    method = average, complete, single (hierarchical) or gromos
    {default: average}

    superpose = 0/1: superimpose all models on the fit part of the
    largest cluster's representative {default: 1}

    group = 0/1: one group per cluster, cluster_1 the largest {default: 1}

    color = 0/1: color the measured part's carbons by cluster
    {default: 1}

    file = str: also save the clusters and the RMSD matrix as CSV

    best_by = iptm, ptm, ranking or plddt: the score that picks each
    cluster's best model (from the predictor's JSON next to the model
    file; pLDDT from the B-factors) {default: iptm}

    scores = str: folder to search for score files by object name, for
    models whose files pymolx doesn't know {default: none}

    scores_csv = str: CSV file (or a pattern for several) with a row per
    model: a model/name/file column (else the first), or seed and sample
    columns (AlphaFold 3's <job>_ranking_scores.csv); score columns such
    as iptm, ptm, ranking_score, or any other (then usable as best_by)
    {default: none}

    higher = auto, 1 or 0: whether a higher best_by score is better
    {default: auto: lower for names like rank, energy, rmsd}

EXAMPLE

    cluster_rmsd *_model, cutoff=2.5

SEE ALSO

    align, super, rms_cur
    '''
    higher = str(higher).strip().lower()
    higher = None if higher in ('auto', '', '-1') else \
        higher not in ('0', 'no', 'false', 'off')
    result = compute(objects, fit, rmsd, atoms, cutoff, method, state,
                     best_by, scores, scores_csv, higher, _self=_self)
    apply(result, int(superpose), int(group), int(color), state=int(state),
          _self=_self)
    if file:
        result.save_csv(_self.exp_path(file))
    if not int(quiet):
        print(' cluster_rmsd: %d structures, %d clusters (%s, cutoff %.1f A);'
              ' fit on %s (%d atoms), RMSD of %s (%d atoms)' % (
                  len(result.models), len(result.clusters), result.method,
                  result.cutoff, result.fit, result.fit_atoms, result.rmsd,
                  result.rmsd_atoms))
        print(result.table())
        if result.unmatched_scores:
            print(' cluster_rmsd: scores without a matching model: %s' %
                  ', '.join(result.unmatched_scores))
    return result


def extend(_self):
    _self.extend('cluster_rmsd', cluster_rmsd)
    _self.cluster_rmsd = cluster_rmsd
