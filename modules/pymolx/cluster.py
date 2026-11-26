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

import pymol
from pymol import cmd

METHODS = ('average', 'complete', 'single', 'gromos')
# what the distance between two models means
MEASURES = {
    'ligand': 'ligand RMSD: fit the antigen, measure the antibody',
    'interface': 'interface RMSD: fit and measure the contacting '
                 'residues of both (CAPRI\'s I-RMSD, blind to how the '
                 'rest of the antibody sits)',
    'epitope': 'epitope overlap: the antigen residues each model '
               'contacts, as 1 - their Jaccard overlap (0 the same '
               'patch, 1 no residue in common)',
}
DEFAULT_MEASURE = 'ligand'
CONTACT = 5.0          # A between heavy atoms: a residue contact
# a distance in the epitope measure is a fraction, not an Angstrom
EPITOPE_CUTOFF = 0.5
SCORES = {'ab_ag_iptm': 'antibody-antigen ipTM', 'iptm': 'ipTM',
          'ptm': 'pTM', 'ranking': 'ranking score', 'plddt': 'pLDDT',
          'ab_ag_pae': 'antibody-antigen PAE'}
# AlphaFold 3's own advice for an antibody-antigen interface is the ipTM
# of that chain pair (docs/output.md on "chain_pair_iptm": "Can be used
# for ranking a specific interface between two chains, when you know that
# they interact, e.g. for antibody-antigen interactions"), not the whole
# complex's ipTM, which the always-confident heavy-light interface lifts
AB_AG_SCORE = 'ab_ag_iptm'
# confident / failed, from the same document's reading of ipTM
IPTM_CONFIDENT = 0.8
IPTM_FAILED = 0.6
# scores of SCORES where a lower value is the better one
LOWER_SCORES = ('ab_ag_pae',)
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
LOWER_IS_BETTER = ('rank', 'energy', 'rmsd', 'dg', 'kcal', 'clash', 'pae')
ATOMS = {'ca': 'name CA', 'backbone': 'name N+CA+C+O',
         'heavy': 'not hydro'}
# the same by atom name (None: all heavy atoms)
ATOM_NAMES = {'ca': ('CA',), 'backbone': ('N', 'CA', 'C', 'O'),
              'heavy': None}
# numbers (coordinates) worked on at a time, when many models are
# superimposed at once: 16 MB each
BATCH = 2000000
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


def chain_pair_value(data, key, chains, antibody, best):
    """
    One value out of a [chains x chains] array of a confidence JSON, for
    the antibody-antigen chain pairs

    :param chains: the structure's chain IDs, in the order the array uses
    :param antibody: the antibody's chain IDs
    :param best: max (ipTM) or min (PAE)
    :return: float, or None
    """
    matrix = data.get(key)
    if not isinstance(matrix, list) or not chains:
        return None
    if len(matrix) != len(chains):
        return None             # (not this structure's chains)
    found = []
    for i, chain_i in enumerate(chains):
        for j, chain_j in enumerate(chains):
            if i == j or (chain_i in antibody) == (chain_j in antibody):
                continue        # (both the antibody's, or both not)
            row = matrix[i]
            if isinstance(row, list) and j < len(row) and \
                    isinstance(row[j], (int, float)):
                found.append(float(row[j]))
    return best(found) if found else None


def read_scores(path, chains=(), antibody=()):
    '''
    :param chains: the structure's chain IDs in order, and which of them
        are the antibody's: for the antibody-antigen chain pair scores
    :return: {'iptm', 'ptm', 'ranking', 'ab_ag_iptm', 'ab_ag_pae': float}
        found in a score JSON
    '''
    with open(path) as handle:
        data = json.load(handle)
    scores = {}
    if chains and antibody and not set(chains) <= set(antibody):
        pair = chain_pair_value(data, 'chain_pair_iptm', chains, antibody,
                                max)
        if pair is not None:
            scores['ab_ag_iptm'] = pair
        pae = chain_pair_value(data, 'chain_pair_pae_min', chains,
                               antibody, min)
        if pae is not None:
            scores['ab_ag_pae'] = pae
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


def folder_listing(folder):
    '''
    :return: [(directory, {file names})] of a folder and those inside it
    '''
    if not folder:
        return []
    return [(root, set(found))
            for root, _, found in os.walk(os.path.expanduser(folder))]


def find_scores(obj, folder='', chains=(), antibody=(), listing=None):
    '''
    The confidence scores of a predicted model: from the JSON next to the
    file it was loaded from, or found in "folder" by the object's name

    :param chains, antibody: the structure's chain IDs in order and the
        antibody's, for the antibody-antigen chain pair scores
    :param listing: folder_listing(folder), when many models are looked
        up (default: the folder is listed here)

    :return: ({'iptm', 'ptm', 'ranking', 'ab_ag_iptm', ...}, path) or
        ({}, None)
    '''
    from . import files
    places = []
    path = files.source(obj)
    if path:
        places.append((os.path.dirname(path), _score_names(_stem(path))))
    if folder:
        names = _score_names(obj)
        if listing is None:
            listing = folder_listing(folder)
        for root, found in listing:
            places.append((root, [n for n in names if n in found]))
    for directory, names in places:
        for name in names:
            candidate = os.path.join(directory, name)
            if os.path.isfile(candidate):
                try:
                    return (read_scores(candidate, chains, antibody),
                            candidate)
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
    if name in LOWER_SCORES:
        return False
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


class Side:
    '''
    One part of one model (e.g. its antigen), as read_sides reads it:

    xyz: array (atoms, 3) of its heavy atoms
    residue: for each of them, the index of its residue in "residues"
    residues: [(chain, resi)]
    keys, coords: [(chain, resi, name)] and array (n, 3) of the atoms to
    compare (those with the wanted names; all of them without names)
    '''

    __slots__ = ('xyz', 'residue', 'residues', 'keys', 'coords')


def _side(chains, resis, atom_names, flat, names):
    import numpy as np
    n = len(chains)
    side = Side()
    side.xyz = np.array(flat, dtype=float).reshape(n, 3)
    side.residues = residues = []
    index = []
    last = None
    for key in zip(chains, resis):
        if key != last:
            residues.append(key)
            last = key
        index.append(len(residues) - 1)
    side.residue = np.array(index, dtype=np.intp)
    # (an atom given twice, as with alternate locations: the last one)
    rows = {(chains[i], resis[i], atom_names[i]): i for i in range(n)
            if names is None or atom_names[i] in names}
    side.keys = list(rows)
    side.coords = side.xyz if len(rows) == n else \
        side.xyz[list(rows.values())]
    return side


def read_sides(models, selection, state=1, names=None, _self=cmd):
    '''
    The heavy atoms of a selection in each model, all read in one pass.
    (Asking PyMOL model by model means a pass over every atom of the
    session each time, and a "within" on top for the contacts: the time
    grows with the square of the number of models, minutes for a few
    hundred. This way it grows with their number.)

    :param names: names of the atoms to compare, e.g. ("CA",) (default:
        all heavy atoms)
    :return: {model: Side}, for every model (empty where nothing of the
        selection is)
    '''
    wanted = set(models)
    names = None if names is None else frozenset(names)
    sides = {}
    chains, resis, atom_names, flat = [], [], [], []
    current = [None, False]     # the model being read, and if it is wanted

    def flush():
        if current[1] and chains:
            sides[current[0]] = _side(chains, resis, atom_names, flat,
                                      names)
        del chains[:], resis[:], atom_names[:], flat[:]

    def add(model, chain, resi, name, x, y, z):
        # (PyMOL goes through the atoms object by object)
        if model != current[0]:
            flush()
            current[:] = model, model in wanted
        if current[1]:
            chains.append(chain)
            resis.append(resi)
            atom_names.append(name)
            flat.append(x)
            flat.append(y)
            flat.append(z)

    _self.iterate_state(state, '(%s) and not hydro and present' % selection,
                        'add(model, chain, resi, name, x, y, z)',
                        space={'add': add})
    flush()
    return {m: sides[m] if m in sides else _side([], [], [], [], names)
            for m in models}


def close_residues(side, other, distance=CONTACT):
    '''
    The residues of "side" with a heavy atom within "distance" of one of
    "other" (two Side of one model)

    :return: {(chain, resi)}
    '''
    import numpy as np
    from scipy.spatial import cKDTree
    if not len(side.xyz) or not len(other.xyz):
        return set()
    nearest, _ = cKDTree(other.xyz).query(side.xyz,
                                          distance_upper_bound=distance)
    return {side.residues[k]
            for k in np.unique(side.residue[np.isfinite(nearest)])}


def joined(first, second, residues=None):
    '''
    The atoms to compare of two Side of one model as one set, e.g. the
    antigen's and the antibody's

    :param residues: {(chain, resi)} to keep (default: all)
    :return: (keys, array (n, 3))
    '''
    import numpy as np
    rows = {}
    xyz = np.concatenate([first.coords, second.coords])
    for i, key in enumerate(first.keys + second.keys):
        if residues is None or key[:2] in residues:
            rows[key] = i
    return list(rows), xyz[list(rows.values())]


def matched_atoms(atoms):
    '''
    The atoms all models have, matched by chain, residue and name

    :param atoms: [(keys, array (n, 3))] per model
    :return: (keys, array (models, atoms, 3))
    '''
    import numpy as np
    first = atoms[0][0]
    if all(keys == first for keys, _ in atoms[1:]):
        # (the usual case, models of one complex: nothing to look up)
        return list(first), np.array([xyz for _, xyz in atoms], dtype=float)
    keys = sorted(set(first).intersection(*(keys for keys, _ in atoms[1:])))
    out = np.zeros((len(atoms), len(keys), 3))
    for n, (model_keys, xyz) in enumerate(atoms):
        row = {key: i for i, key in enumerate(model_keys)}
        if keys:
            out[n] = xyz[[row[key] for key in keys]]
    return keys, out


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


def rotations(mobile, target):
    '''
    kabsch for many structures at once

    :param mobile: array (models, atoms, 3), each centered on its origin
    :param target: array (atoms, 3) or (models, atoms, 3), centered too
    :return: array (models, 3, 3), the rotations of "mobile" onto "target"
    '''
    import numpy as np
    u, _, vt = np.linalg.svd(mobile.transpose(0, 2, 1) @ target)
    v = vt.transpose(0, 2, 1).copy()
    ut = u.transpose(0, 2, 1)
    v[:, :, 2] *= np.sign(np.linalg.det(v @ ut))[:, None]
    return v @ ut


def batches(first, end, atoms):
    '''
    (start, stop) over first..end, so many models at a time that their
    atoms' coordinates stay below BATCH numbers
    '''
    size = max(1, BATCH // max(1, 3 * atoms))
    return [(start, min(end, start + size))
            for start in range(first, end, size)]


def rmsd_matrix(fit, measure):
    '''
    :param fit: array (models, n, 3), or a list of arrays (n, 3): the
        atoms to superimpose on
    :param measure: the same for the atoms to measure (m per model)
    :return: symmetric matrix (models x models), Angstrom
    '''
    import numpy as np
    fit = np.asarray(fit, dtype=float)
    measure = np.asarray(measure, dtype=float)
    n = len(fit)
    matrix = np.zeros((n, n))
    # each model with its fit atoms' centroid at the origin
    centers = fit.mean(axis=1)[:, None, :]
    fit = fit - centers
    measure = measure - centers
    for i in range(n - 1):
        # one model against all after it, in one go
        for start, stop in batches(i + 1, n, measure.shape[1]):
            rotation = rotations(fit[start:stop], fit[i])
            moved = measure[start:stop] @ rotation.transpose(0, 2, 1)
            moved -= measure[i]
            d = np.sqrt(np.einsum('jak,jak->j', moved, moved) /
                        measure.shape[1])
            matrix[i, start:stop] = matrix[start:stop, i] = d
    return matrix


def contact_residues(obj, side, other, state, cutoff=CONTACT, _self=cmd):
    '''
    The residues of "side" with a heavy atom within cutoff of "other"

    :return: {(chain, resi)}
    '''
    found = set()
    _self.iterate('(%%%s and (%s) and not hydro) within %g of '
                  '(%%%s and (%s) and not hydro)'
                  % (obj, side, cutoff, obj, other),
                  'found.add((chain, resi))', space={'found': found})
    return found


def _residues_selection(keys):
    """
    "(chain A and resi 5+6) or (chain B and resi 9)" for {(chain, resi)}
    """
    by_chain = {}
    for chain, resi in sorted(keys):
        by_chain.setdefault(chain, []).append(str(resi))
    return ' or '.join('(chain "%s" and resi %s)' % (chain, '+'.join(res))
                       for chain, res in sorted(by_chain.items())) or 'none'


def jaccard_matrix(sets):
    '''
    1 - the Jaccard overlap of each pair of sets (0: the same residues)

    :return: symmetric matrix
    '''
    import numpy as np
    n = len(sets)
    column = {item: k for k, item in enumerate(set().union(*sets))}
    member = np.zeros((n, len(column)))
    for row, items in enumerate(sets):
        member[row, [column[item] for item in items]] = 1.0
    shared = member @ member.T
    sizes = member.sum(axis=1)
    union = sizes[:, None] + sizes[None, :] - shared
    matrix = np.where(union > 0, 1.0 - shared / np.maximum(union, 1.0), 1.0)
    np.fill_diagonal(matrix, 0.0)
    return matrix


def gromos(matrix, cutoff):
    '''
    :return: cluster number (0-based, by size) per model
    '''
    import numpy as np
    within = np.asarray(matrix) <= cutoff
    n = len(within)
    labels = np.full(n, -1)
    left = np.ones(n, dtype=bool)
    # neighbors among the models still without a cluster
    counts = within.sum(axis=1)
    cluster = 0
    while left.any():
        # (the first of those with the most neighbors)
        center = int(np.argmax(np.where(left, counts, -1)))
        members = left & within[center]
        members[center] = True
        labels[members] = cluster
        left &= ~members
        counts -= within[:, members].sum(axis=1)
        cluster += 1
    return labels.tolist()


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


def set_aside(labels, min_size):
    '''
    Clusters of fewer than min_size models put together as one last
    "cluster": the models that found no company

    :return: (labels, the label of that pool or None)
    '''
    labels = [int(c) for c in labels]
    small = {c for c in set(labels) if labels.count(c) < min_size}
    if not small:
        return labels, None
    kept = sorted(set(labels) - small)      # (still largest first)
    rank = {c: k for k, c in enumerate(kept)}
    return [rank.get(c, len(kept)) for c in labels], len(kept)


def cluster_label(c):
    '''
    "3" for a cluster, "-" for the models left without one
    '''
    return '-' if c.get('unclustered') else str(c['cluster'])


def group_name(prefix, c):
    return '%s_unclustered' % prefix if c.get('unclustered') else \
        '%s_%d' % (prefix, c['cluster'])


class Result:
    '''
    models, matrix, labels and per-cluster summaries

    min_size: clusters with fewer models are no clusters; their models
    come last as one entry marked "unclustered", without a representative
    or a best model
    '''

    def __init__(self, models, matrix, labels, scores, fit_atoms,
                 rmsd_atoms, confidences=None, best_by='iptm', higher=None,
                 min_size=1):
        self.min_size = max(int(min_size), 1)
        labels, pool = set_aside(labels, self.min_size)
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
        self.measure = DEFAULT_MEASURE
        self.clusters = []
        for c in sorted(set(labels)):
            members = [i for i, label in enumerate(labels) if label == c]
            if c == pool:
                known = [scores[i] for i in members if scores[i] is not None]
                self.clusters.append({
                    'unclustered': True, 'mean_ab_ag_iptm': None,
                    'cluster': c + 1, 'size': len(members),
                    'representative': None,
                    'members': [models[i] for i in members],
                    'mean_rmsd': None,
                    'score': sum(known) / len(known) if known else None,
                    'best': None, 'best_score': None, 'best_by': None})
                continue
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
            confident = [self.confidences[i].get(AB_AG_SCORE)
                         for i in members]
            confident = [v for v in confident if isinstance(v, (int, float))]
            self.clusters.append({
                'mean_ab_ag_iptm': (sum(confident) / len(confident)
                                    if confident else None),
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

    def ranked(self):
        '''
        The clusters, most trustworthy first: a populated cluster of
        confident models is the signal, a confident singleton rarely is.
        Ranked by the number of models and the mean antibody-antigen
        ipTM of their interfaces, when that is known.

        :return: [cluster], a copy of the list
        '''
        known = [c for c in self.clusters
                 if c.get('mean_ab_ag_iptm') is not None]
        if not known:
            return list(self.clusters)      # (by size, as they are)
        most = max(c['size'] for c in self.clusters)
        return sorted(self.clusters, key=lambda c: (
            bool(c.get('unclustered')),
            -(c['size'] / most * (c.get('mean_ab_ag_iptm') or 0.0))))

    def table(self):
        pair = any(c.get('mean_ab_ag_iptm') is not None
                   for c in self.clusters)
        unit = '' if self.measure == 'epitope' else ' A'
        spread = 'Shared' if self.measure == 'epitope' else 'Mean RMSD'
        width = max([len('Best')] + [len(format_best(c))
                                     for c in self.clusters])
        lines = [' %-7s %4s  %-9s %-6s %s %-*s  %s' % (
            'Cluster', 'Size', spread, 'pLDDT',
            'Ab-Ag ipTM ' if pair else '', width, 'Best',
            'Representative')]
        for c in self.clusters:
            value = c.get('mean_ab_ag_iptm')
            if c.get('unclustered'):
                within = '%7s%-2s' % ('-', '')
            else:
                within = '%7.2f%-2s' % (
                    1 - c['mean_rmsd'] if self.measure == 'epitope'
                    else c['mean_rmsd'], unit)
            lines.append(' %-7s %4d  %s %-6s %s%-*s  %s' % (
                cluster_label(c), c['size'], within,
                '%.1f' % c['score'] if c['score'] is not None else '-',
                ('%-11s' % ('%.2f' % value if value is not None else '-'))
                if pair else '',
                width, format_best(c), c['representative'] or '-'))
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
            # 0: a model left without a cluster (min_size)
            number = {m: 0 if c.get('unclustered') else c['cluster']
                      for c in self.clusters for m in c['members']}

            def value(conf, key, digits):
                v = conf.get(key)
                return '' if v is None else '%.*f' % (digits, v)

            for i, model in enumerate(self.models):
                conf = self.confidences[i]
                out.writerow([model, number[model], reps[model] == model,
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


def scores_table(best_by, scores_csv, _self=cmd):
    '''
    The scores CSV, read, and the score to pick the best model by,
    checked against what there is to pick from

    :return: (best_by, {(job, row name): scores})
    '''
    best_by = str(best_by).strip().lower()
    table = {}
    if str(scores_csv).strip():
        table = read_scores_csvs(str(scores_csv).strip(), _self=_self)
    columns = {k for values in table.values() for k in values}
    if best_by not in SCORES and best_by not in columns:
        raise ValueError('best_by must be one of: ' + ', '.join(
            list(SCORES) + sorted(columns - set(SCORES))))
    return best_by, table


def checked(method, atoms):
    '''
    :return: (method, atoms) in lower case, or ValueError
    '''
    method = str(method).strip().lower()
    if method not in METHODS:
        raise ValueError('method must be one of: ' + ', '.join(METHODS))
    atoms = str(atoms).strip().lower()
    if atoms not in ATOMS:
        raise ValueError('atoms must be one of: ' + ', '.join(ATOMS))
    return method, atoms


def model_scores(models, part, state, scores, table, _self=cmd):
    '''
    What is known about each model's quality: the mean B-factor of
    "part" if it looks like pLDDT (predicted models store it there), the
    predictor's confidence scores, and the rows of a scores table (which
    come first)

    :return: ([pLDDT or None], [{score: float}], [rows without a model])
    '''
    # (all models in one pass: see read_sides)
    found = {m: [] for m in models}

    def add(model, b):
        values = found.get(model)
        if values is not None:
            values.append(b)

    _self.iterate_state(state, '(%s) and name CA' % part, 'add(model, b)',
                        space={'add': add})
    plddt = []
    for m in models:
        values = found[m]
        mean = sum(values) / len(values) if values else None
        plddt.append(mean if mean is not None and 0 < mean <= 100 else None)
    order = sorted(_self.get_chains(models[0]) or [])
    antibody_ids = antibody_chains(models[0], _self)
    listing = folder_listing(scores)
    confidences = [find_scores(m, scores, order, antibody_ids, listing)[0]
                   for m in models]
    unmatched = []
    if table:
        matched, unmatched = match_scores(models, table)
        for m, conf in zip(models, confidences):
            conf.update(matched.get(m, {}))  # (the CSV first)
    return plddt, confidences, unmatched


def compute(objects='', fit='', rmsd='', atoms='ca', cutoff=3.0,
            method='average', state=1, best_by='iptm', scores='',
            scores_csv='', higher=None, measure=DEFAULT_MEASURE,
            _self=cmd):
    '''
    The clustering, without changing the session.

    :return: Result
    '''
    import numpy as np
    best_by, table = scores_table(best_by, scores_csv, _self)
    method, atoms = checked(method, atoms)
    measure = str(measure).strip().lower()
    if measure not in MEASURES:
        raise ValueError('measure must be one of: ' + ', '.join(MEASURES))
    models = _objects(objects, _self)
    if len(models) < 2:
        raise ValueError('need at least two structures (objects)')
    if not str(fit).strip() or not str(rmsd).strip():
        auto_fit, auto_rmsd = default_selections(models[0], _self)
        fit = str(fit).strip() or auto_fit
        rmsd = str(rmsd).strip() or auto_rmsd
    state = int(state)

    cutoff = float(cutoff)
    # both parts of every model, each read in one pass
    fits = read_sides(models, fit, state, ATOM_NAMES[atoms], _self)
    parts = read_sides(models, rmsd, state, ATOM_NAMES[atoms], _self)
    epitopes = [close_residues(fits[m], parts[m]) for m in models]
    if measure == 'epitope':
        # where each model binds, not how the antibody is turned
        empty = [m for m, e in zip(models, epitopes) if not e]
        if empty:
            raise ValueError('nothing of "%s" touches "%s" in %s: no '
                             'epitope to compare' % (rmsd, fit,
                                                     ', '.join(empty[:3])))
        matrix = jaccard_matrix(epitopes)
        fit_keys = rmsd_keys = sorted(set.union(*epitopes))
    else:
        if measure == 'interface':
            # the contacting residues of both sides, in any model: fitted
            # on and measured, as CAPRI's interface RMSD
            paratopes = [close_residues(parts[m], fits[m]) for m in models]
            residues = set.union(*epitopes) | set.union(*paratopes)
            if not residues:
                raise ValueError('"%s" and "%s" do not touch in any '
                                 'structure: no interface to compare'
                                 % (fit, rmsd))
            where = _residues_selection(residues)
            fit_sel = rmsd_sel = '(%s) and %s' % (where, ATOMS[atoms])
            fit_keys, fit_xyz = rmsd_keys, rmsd_xyz = matched_atoms(
                [joined(fits[m], parts[m], residues) for m in models])
        else:
            fit_sel = '(%s) and %s' % (fit, ATOMS[atoms])
            rmsd_sel = '(%s) and %s' % (rmsd, ATOMS[atoms])
            fit_keys, fit_xyz = matched_atoms(
                [(fits[m].keys, fits[m].coords) for m in models])
            rmsd_keys, rmsd_xyz = matched_atoms(
                [(parts[m].keys, parts[m].coords) for m in models])
        if len(fit_keys) < 3:
            raise ValueError('fewer than 3 atoms in "%s" common to all '
                             'structures' % fit_sel)
        if not rmsd_keys:
            raise ValueError('no atoms in "%s" common to all structures' %
                             rmsd_sel)
        matrix = rmsd_matrix(fit_xyz, rmsd_xyz)
    labels = gromos(matrix, cutoff) if method == 'gromos' else \
        hierarchical(matrix, cutoff, method)
    plddt, confidences, unmatched = model_scores(models, rmsd, state,
                                                 scores, table, _self)
    result = Result(models, matrix, labels, plddt, len(fit_keys),
                    len(rmsd_keys), confidences, best_by, higher)
    result.unmatched_scores = unmatched
    result.fit, result.rmsd, result.method, result.cutoff = \
        fit, rmsd, method, cutoff
    result.measure = measure
    result.epitopes = epitopes
    return result


CLUSTER_COLORS = ['tv_blue', 'tv_orange', 'forest', 'raspberry', 'teal',
                  'olive', 'violetpurple', 'deepsalmon', 'marine', 'limon']
UNCLUSTERED_COLOR = 'grey50'    # models left without a cluster


def apply(result, superpose=True, group=True, color=True, prefix='cluster',
          state=1, _self=cmd):
    '''
    Show the clustering: all models superimposed on the fit part of the
    largest cluster's representative, one group per cluster (the
    representative first), the measured part colored by cluster
    '''
    # (no representative: every model was left without a cluster)
    reference = result.clusters[0]['representative'] or result.models[0]
    if superpose:
        # (all models read in one pass, and each moved by its name, not
        # as a selection: see read_sides)
        sides = read_sides(result.models, '(%s) and name CA' % result.fit,
                           state, _self=_self)
        target = sides[reference]
        for model in result.models:
            if model == reference:
                continue
            keys, (mobile, onto) = matched_atoms([
                (sides[model].keys, sides[model].coords),
                (target.keys, target.coords)])
            if len(keys) < 3:
                continue
            rotation, translation = kabsch(mobile, onto)
            try:
                xyz = _self.get_coordset(model, state)
            except pymol.CmdException:
                xyz = None      # (no such state)
            if xyz is None:
                continue
            _self.load_coordset(xyz @ rotation.T + translation, model, state)
    if group:
        # groups of an earlier clustering: emptied, then removed (deleting
        # a group would delete the models in it)
        pattern = re.compile(r'%s_(\d+|unclustered)$' % re.escape(prefix))
        for name in _self.get_names('group_objects') or []:
            members = _self.get_object_list('(%s)' % name) or []
            if pattern.match(name) and set(members) <= set(result.models):
                _self.ungroup(name)
                _self.delete(name)
    colors = {}     # model: color index
    for c in result.clusters:
        if group:
            name = _self.get_legal_name(group_name(prefix, c))
            first = [m for m in (c['best'], c['representative']) if m]
            members = list(dict.fromkeys(first)) + [
                m for m in c['members'] if m not in first]
            _self.group(name, ' '.join(members))
            _self.order(' '.join(members), location='current')
        if color:
            colour = UNCLUSTERED_COLOR if c.get('unclustered') else \
                CLUSTER_COLORS[(c['cluster'] - 1) % len(CLUSTER_COLORS)]
            for m in c['members']:
                colors[m] = _self.get_color_index(colour)
    if colors:
        # (every model's carbons in one pass, not a "color" per model)
        _self.alter('(%s) and elem C' % result.rmsd,
                    'color = colors.get(model, color)',
                    space={'colors': colors})
        _self.recolor()
    if group:
        _self.order(' '.join(group_name(prefix, c)
                             for c in result.clusters), location='current')


def cluster_rmsd(objects='', fit='', rmsd='', atoms='ca', cutoff=3.0,
                 method='average', superpose=1, group=1, color=1,
                 file='', best_by='iptm', scores='', scores_csv='',
                 higher='auto', state=1, measure=DEFAULT_MEASURE, quiet=0,
                 *, _self=cmd):
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

    measure = what the distance between two models is {default: ligand}:
    ligand (fit the antigen, measure the antibody: CAPRI's ligand RMSD),
    interface (fit and measure the residues the two sides contact:
    CAPRI's interface RMSD, which the antibody's distal domains don't
    sway) or epitope (the antigen residues each model contacts, compared
    as 1 - their overlap: where it binds, not how it is turned)

    cutoff = float: RMSD in Angstrom, or for measure=epitope a fraction
    of the contacted residues not shared, 0 to 1 {default: 3.0, and 0.5
    for epitope}

    method = average, complete, single (hierarchical) or gromos
    {default: average}

    superpose = 0/1: superimpose all models on the fit part of the
    largest cluster's representative {default: 1}

    group = 0/1: one group per cluster, cluster_1 the largest {default: 1}

    color = 0/1: color the measured part's carbons by cluster
    {default: 1}

    file = str: also save the clusters and the RMSD matrix as CSV

    best_by = ab_ag_iptm, iptm, ptm, ranking, plddt or ab_ag_pae: the
    score that picks each cluster's best model (from the predictor's
    JSON next to the model file; pLDDT from the B-factors). ab_ag_iptm
    is AlphaFold 3's ipTM for the antibody-antigen chain pair, which is
    what its own documentation suggests for such an interface, rather
    than the whole complex's {default: iptm}

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
    measure = str(measure).strip().lower()
    if measure == 'epitope' and float(cutoff) == 3.0:
        cutoff = EPITOPE_CUTOFF     # (a fraction, not an Angstrom)
    result = compute(objects, fit, rmsd, atoms, cutoff, method, state,
                     best_by, scores, scores_csv, higher, measure,
                     _self=_self)
    apply(result, int(superpose), int(group), int(color), state=int(state),
          _self=_self)
    if file:
        result.save_csv(_self.exp_path(file))
    if not int(quiet):
        if result.measure == 'epitope':
            print(' cluster_rmsd: %d structures, %d clusters (%s, cutoff '
                  '%.2f of the epitope not shared); the residues of %s '
                  'that %s touches' % (
                      len(result.models), len(result.clusters),
                      result.method, result.cutoff, result.fit,
                      result.rmsd))
        else:
            print(' cluster_rmsd: %d structures, %d clusters (%s, cutoff '
                  '%.1f A); %s, fit on %d atoms, RMSD of %d' % (
                      len(result.models), len(result.clusters),
                      result.method, result.cutoff,
                      MEASURES[result.measure].split(':')[0],
                      result.fit_atoms, result.rmsd_atoms))
        print(result.table())
        if result.unmatched_scores:
            print(' cluster_rmsd: scores without a matching model: %s' %
                  ', '.join(result.unmatched_scores))
    return result


def extend(_self):
    _self.extend('cluster_rmsd', cluster_rmsd)
    _self.cluster_rmsd = cluster_rmsd
