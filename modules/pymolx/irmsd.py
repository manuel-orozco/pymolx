'''
Clustering of complexes by interface RMSD (command "cluster_irmsd",
Interface RMSD Clustering dialog), made for the many models an epitope
prediction gives (pymolx.epitope: seeds x samples of one antibody-antigen
complex): which of them agree on how the two meet?

The interface RMSD (i-RMSD) is CAPRI's (Lensink et al.): the residues of
either side with an atom within 10 A of the other side are the interface;
two structures are superimposed on the backbone atoms of those residues
and the RMSD is taken over the same atoms. It looks only at the place
where the two meet, so a constant domain that swings, or a floppy
terminus, does not blur it. CAPRI reads it as high (up to 1 A), medium
(2 A) or acceptable (4 A) agreement.

CAPRI has a reference structure to take the interface from. Among models
there is none, so there are two ways to choose the residues ("scope"):

- pair: for each two models, the residues that are interface in either
  of them. Two models on the same epitope are compared over that epitope
  alone; two on different epitopes over both, where they cannot agree.
- all: the residues that are interface in any model, the same for every
  pair (what pymolx.cluster's measure=interface does, and HADDOCK's RMSD
  clustering). With models spread over the antigen this approaches the
  RMSD of the whole complex.

Models whose two sides do not touch have no interface and are left out.
Clustering, the clusters' representatives and best models, the scores and
the way the result is shown are pymolx.cluster's; clusters smaller than
"min_size" are no clusters, and their models are kept together as
"unclustered".
'''

from pymol import cmd

from . import cluster

# A between heavy atoms of the two sides: an interface residue (CAPRI)
INTERFACE = 10.0
# how far apart two models of one binding mode may be: CAPRI's
# "acceptable" (high 1 A, medium 2 A)
CUTOFF = 4.0
CAPRI = (('high', 1.0), ('medium', 2.0), ('acceptable', 4.0))
SCOPES = {
    'pair': 'each pair of models: the residues that are interface in '
            'either of the two',
    'all': 'all models: the residues that are interface in any of them',
}
DEFAULT_SCOPE = 'pair'
DEFAULT_ATOMS = 'backbone'      # CAPRI's choice
MIN_SIZE = 1
# A: pairs closer than this are measured atom by atom (irmsd_matrix)
EXACT_BELOW = 0.01
# fewer atoms than this cannot be superimposed
MIN_ATOMS = 3
# the atoms compared, by name (the keys of cluster.ATOMS); None: all heavy
ATOM_NAMES = cluster.ATOM_NAMES


def capri_class(value):
    '''
    "high", "medium", "acceptable" or "" for an i-RMSD
    '''
    for name, limit in CAPRI:
        if value <= limit:
            return name
    return ''


def interface_residues(obj, side, other, state=1, distance=INTERFACE,
                       _self=cmd):
    '''
    The residues of both sides with a heavy atom within "distance" of the
    other side

    :return: {(chain, resi)}
    '''
    return (cluster.contact_residues(obj, side, other, state, distance,
                                     _self=_self) |
            cluster.contact_residues(obj, other, side, state, distance,
                                     _self=_self))


def irmsd(first, second):
    '''
    The RMSD of two sets of matched atoms after their best superposition

    :param first, second: arrays (n, 3)
    '''
    import numpy as np
    rotation, translation = cluster.kabsch(second, first)
    moved = second @ rotation.T + translation
    return float(np.sqrt(((moved - first) ** 2).sum(axis=1).mean()))


def _pair_rows(xyz, masks, i, others):
    '''
    The i-RMSD of model i to each of "others" (indices), each pair over
    the atoms that are interface in either: superimposed and measured
    atom by atom

    :return: array (others)
    '''
    import numpy as np
    out = np.zeros(len(others))
    for start, stop in cluster.batches(0, len(others), masks.shape[1]):
        rows = others[start:stop]
        # (atoms outside a pair's interface weigh nothing)
        use = (masks[rows] | masks[i])[:, :, None]
        counts = use.sum(axis=(1, 2))
        mobile = np.where(use, xyz[rows], 0.0)
        target = np.where(use, xyz[i], 0.0)
        mobile -= np.where(
            use, (mobile.sum(axis=1) / counts[:, None])[:, None, :], 0.0)
        target -= np.where(
            use, (target.sum(axis=1) / counts[:, None])[:, None, :], 0.0)
        rotation = cluster.rotations(mobile, target)
        moved = mobile @ rotation.transpose(0, 2, 1)
        moved -= target
        out[start:stop] = np.sqrt(
            np.einsum('jak,jak->j', moved, moved) / counts)
    return out


def irmsd_matrix(xyz, masks, scope=DEFAULT_SCOPE):
    '''
    :param xyz: array (models, atoms, 3), the same atoms in each model
    :param masks: bool array (models, atoms): the atoms of each model's
        interface residues
    :return: (symmetric matrix (models x models) in Angstrom, the fewest
        and the most atoms a pair was compared over)
    '''
    import numpy as np
    xyz = np.asarray(xyz, dtype=float)
    masks = np.asarray(masks, dtype=bool)
    n = len(xyz)
    if n < 2:
        return np.zeros((n, n)), (0, 0)
    # only atoms that are interface in some model can count
    shared = masks.any(axis=0)
    xyz, masks = xyz[:, shared], masks[:, shared]
    if scope == 'all':
        # the same atoms for every pair
        count = xyz.shape[1]
        return cluster.rmsd_matrix(xyz, xyz), (count, count)

    # All pairs at once. A pair (i, j) is compared over the atoms in
    # mask i or mask j, and a sum over those is the sum over i's, plus
    # the sum over j's, minus the sum over the atoms in both: each of
    # them a product of two (models x atoms) matrices. From the sums
    # come each pair's centroids, the 3 x 3 matrix that Kabsch's method
    # takes apart, and the RMSD itself (the summed squares less twice
    # the singular values), without a pair of structures ever being
    # laid on each other.
    x = xyz - xyz.reshape(-1, 3).mean(axis=0)   # (smaller numbers)
    m = masks.astype(float)
    axes = [x[:, :, k] for k in range(3)]
    masked = [m * axis for axis in axes]

    def over_pairs(values, in_mask):
        # [i, j]: the sum of model j's values over the atoms of pair
        # (i, j); in_mask is m * values
        return m @ values.T + in_mask.sum(axis=1)[None, :] - m @ in_mask.T

    own = m.sum(axis=1)
    counts = own[:, None] + own[None, :] - m @ m.T
    # [i, j, k]: the sum of coordinate k of model j, then of model i
    sum_j = np.stack([over_pairs(axes[k], masked[k]) for k in range(3)],
                     axis=2)
    sum_i = sum_j.transpose(1, 0, 2)
    squares = (x * x).sum(axis=2)
    square_j = over_pairs(squares, m * squares)
    spread = square_j + square_j.T - (
        (sum_j * sum_j).sum(axis=2) + (sum_i * sum_i).sum(axis=2)) / counts
    rows, columns = np.triu_indices(n, 1)
    # [pair, k, l]: coordinate k of the later model (j, the one moved)
    # times coordinate l of the earlier one (i), summed over the atoms
    h = np.empty((len(rows), 3, 3))
    for k in range(3):
        for l in range(3):
            cross = masked[l] @ axes[k].T + axes[l] @ masked[k].T - \
                masked[l] @ masked[k].T
            h[:, k, l] = cross[rows, columns]
    h -= sum_j[rows, columns][:, :, None] * sum_i[rows, columns][:, None, :] \
        / counts[rows, columns][:, None, None]
    u, s, vt = np.linalg.svd(h)
    s[:, 2] *= np.sign(np.linalg.det(u) * np.linalg.det(vt))
    squared = (spread[rows, columns] - 2.0 * s.sum(axis=1)) / \
        counts[rows, columns]
    values = np.sqrt(np.maximum(squared, 0.0))
    matrix = np.zeros((n, n))
    matrix[rows, columns] = matrix[columns, rows] = values
    # A difference of large sums says little about a pair of (nearly)
    # the same structure: those are superimposed and measured.
    for i in np.unique(rows[values < EXACT_BELOW]):
        others = np.flatnonzero(matrix[i, i + 1:] < EXACT_BELOW) + i + 1
        matrix[i, others] = matrix[others, i] = _pair_rows(
            xyz, masks, i, others)
    compared = counts[rows, columns]
    return matrix, (int(round(compared.min())), int(round(compared.max())))


def compute(objects='', antigen='', antibody='', atoms=DEFAULT_ATOMS,
            cutoff=CUTOFF, method='average', distance=INTERFACE,
            scope=DEFAULT_SCOPE, min_size=MIN_SIZE, state=1, best_by='iptm',
            scores='', scores_csv='', higher=None, _self=cmd):
    '''
    The clustering, without changing the session.

    :return: cluster.Result, with "interfaces" ({(chain, resi)} per
        model), "skipped" (models without an interface, left out) and
        "atoms_compared" (fewest, most per pair)
    '''
    import numpy as np
    best_by, table = cluster.scores_table(best_by, scores_csv, _self)
    method, atoms = cluster.checked(method, atoms)
    scope = str(scope).strip().lower()
    if scope not in SCOPES:
        raise ValueError('scope must be one of: ' + ', '.join(SCOPES))
    cutoff, distance, min_size = float(cutoff), float(distance), int(min_size)
    if cutoff <= 0:
        raise ValueError('cutoff: an RMSD above 0 A')
    if distance <= 0:
        raise ValueError('distance: above 0 A')
    if min_size < 1:
        raise ValueError('min_size: 1 or more structures')
    models = cluster._objects(objects, _self)
    if len(models) < 2:
        raise ValueError('need at least two structures (objects)')
    if not str(antigen).strip() or not str(antibody).strip():
        auto_antigen, auto_antibody = cluster.default_selections(models[0],
                                                                 _self)
        antigen = str(antigen).strip() or auto_antigen
        antibody = str(antibody).strip() or auto_antibody
    if antigen == antibody:
        raise ValueError('an interface has two sides: give the antigen '
                         'and the antibody (no antibody chains were found '
                         'in %s)' % models[0])
    state = int(state)

    # (both sides of every model, each read in one pass)
    antigens = cluster.read_sides(models, antigen, state, ATOM_NAMES[atoms],
                                  _self)
    antibodies = cluster.read_sides(models, antibody, state,
                                    ATOM_NAMES[atoms], _self)
    keys, xyz = cluster.matched_atoms(
        [cluster.joined(antigens[m], antibodies[m]) for m in models])
    if len(keys) < MIN_ATOMS:
        raise ValueError('fewer than %d %s atoms of "%s" and "%s" common '
                         'to all structures' % (MIN_ATOMS, atoms, antigen,
                                                antibody))
    epitopes = [cluster.close_residues(antigens[m], antibodies[m], distance)
                for m in models]
    interfaces = [epitope | cluster.close_residues(antibodies[m],
                                                   antigens[m], distance)
                  for m, epitope in zip(models, epitopes)]
    residues = [key[:2] for key in keys]
    masks = np.array([[r in interface for r in residues]
                      for interface in interfaces], dtype=bool)
    # a model whose sides do not touch has nothing to compare
    touching = masks.sum(axis=1) >= MIN_ATOMS
    skipped = [m for m, ok in zip(models, touching) if not ok]
    if touching.sum() < 2:
        raise ValueError('"%s" and "%s" come within %g A of each other in '
                         '%d structure(s): no interfaces to compare'
                         % (antigen, antibody, distance, touching.sum()))
    models = [m for m, ok in zip(models, touching) if ok]
    interfaces = [i for i, ok in zip(interfaces, touching) if ok]
    epitopes = [e for e, ok in zip(epitopes, touching) if ok]
    masks = masks[touching]
    xyz = xyz[touching]

    matrix, compared = irmsd_matrix(xyz, masks, scope)
    labels = cluster.gromos(matrix, cutoff) if method == 'gromos' else \
        cluster.hierarchical(matrix, cutoff, method)
    plddt, confidences, unmatched = cluster.model_scores(
        models, antibody, state, scores, table, _self)
    result = cluster.Result(models, matrix, labels, plddt, compared[0],
                            compared[1], confidences, best_by, higher,
                            min_size=min_size)
    result.unmatched_scores = unmatched
    # (cluster.apply superimposes on "fit" and colors "rmsd")
    result.fit, result.rmsd, result.method, result.cutoff = \
        antigen, antibody, method, cutoff
    result.measure = 'interface'
    result.scope, result.distance, result.atoms = scope, distance, atoms
    result.interfaces = interfaces
    result.epitopes = epitopes      # (the antigen's share of them)
    result.skipped = skipped
    result.atoms_compared = compared
    return result


def summary(result):
    '''
    One line on a result: what was clustered, and over what
    '''
    import numpy as np
    clusters = [c for c in result.clusters if not c.get('unclustered')]
    alone = sum(c['size'] for c in result.clusters if c.get('unclustered'))
    low, high = result.atoms_compared
    line = '%d structures, %d cluster%s' % (
        len(result.models), len(clusters), '' if len(clusters) == 1 else 's')
    if alone:
        line += ', %d unclustered (in clusters of fewer than %d)' % (
            alone, result.min_size)
    line += '; i-RMSD over %s %s atoms of residues within %g A' % (
        '%d' % low if low == high else '%d to %d' % (low, high),
        result.atoms, result.distance)
    # (what there is to cut: a cutoff outside this range changes nothing)
    apart = result.matrix[np.triu_indices(len(result.models), 1)]
    line += ', pairs %.1f to %.1f A apart' % (apart.min(), apart.max())
    if result.skipped:
        line += '; %d left out, their sides do not touch: %s' % (
            len(result.skipped), ', '.join(result.skipped[:5]) +
            (', ...' if len(result.skipped) > 5 else ''))
    return line


def cluster_irmsd(objects='', antigen='', antibody='', cutoff=CUTOFF,
                  min_size=MIN_SIZE, distance=INTERFACE, scope=DEFAULT_SCOPE,
                  atoms=DEFAULT_ATOMS, method='average', superpose=1,
                  group=1, color=1, file='', best_by='iptm', scores='',
                  scores_csv='', higher='auto', state=1, quiet=0,
                  *, _self=cmd):
    '''
DESCRIPTION

    "cluster_irmsd" clusters complexes (objects) by interface RMSD, e.g.
    the many antibody-antigen models of an epitope prediction: the
    residues of either side close to the other side are superimposed and
    their RMSD is taken, as CAPRI does, so only the place where the two
    meet counts. Atoms are matched by chain, residue and name.

USAGE

    cluster_irmsd [ objects [, antigen [, antibody [, cutoff
        [, min_size [, distance [, scope [, atoms [, method ]]]]]]]]]

ARGUMENTS

    objects = str: object names or patterns (e.g. "*_model") {default:
    all molecular objects}

    antigen, antibody = str: selections of the two sides {default: the
    antibody chains, found by their sequences, and the other chains}

    cutoff = float: the largest interface RMSD, in Angstrom, at which
    models still belong together. CAPRI calls an agreement within 1 A
    high, 2 A medium and 4 A acceptable {default: 4.0}

    min_size = int: clusters of fewer models are no clusters; their
    models are kept together as "unclustered" {default: 1}

    distance = float: a residue is interface if one of its heavy atoms
    is within this many Angstrom of the other side {default: 10.0,
    CAPRI's; 5.0 for the residues in direct contact}

    scope = pair or all: where the interface residues come from. pair:
    for each two models, those of either of them. all: those of any
    model, the same for every pair {default: pair}

    atoms = backbone, ca or heavy {default: backbone, CAPRI's}

    method = average, complete, single (hierarchical) or gromos
    {default: average}

    superpose = 0/1: superimpose all models on the antigen of the
    largest cluster's representative {default: 1}

    group = 0/1: one group per cluster, cluster_1 the largest, and
    cluster_unclustered {default: 1}

    color = 0/1: color the antibody's carbons by cluster {default: 1}

    file = str: also save the clusters and the i-RMSD matrix as CSV

    best_by, scores, scores_csv, higher: the score that picks each
    cluster's best model and where to find it, as for "cluster_rmsd"

EXAMPLE

    cluster_irmsd *_model, cutoff=2.5, min_size=3

NOTES

    Models whose antibody and antigen do not come within "distance" of
    each other have no interface: they are left out, and named.

SEE ALSO

    cluster_rmsd
    '''
    higher = str(higher).strip().lower()
    higher = None if higher in ('auto', '', '-1') else \
        higher not in ('0', 'no', 'false', 'off')
    result = compute(objects, antigen, antibody, atoms, cutoff, method,
                     distance, scope, min_size, state, best_by, scores,
                     scores_csv, higher, _self=_self)
    cluster.apply(result, int(superpose), int(group), int(color),
                  state=int(state), _self=_self)
    if file:
        result.save_csv(_self.exp_path(file))
    if not int(quiet):
        print(' cluster_irmsd: %s (%s, cutoff %.1f A)' % (
            summary(result), result.method, result.cutoff))
        print(result.table())
        if result.unmatched_scores:
            print(' cluster_irmsd: scores without a matching model: %s' %
                  ', '.join(result.unmatched_scores))
    return result


def extend(_self):
    _self.extend('cluster_irmsd', cluster_irmsd)
    _self.cluster_irmsd = cluster_irmsd
