'''
Hydrogen-bond optimization for pymolx.prep, like Reduce (MolProbity) and
Maestro's "Optimize H-bond assignments". On OpenMM topologies, no PyMOL
calls.

- flips (heavy atoms, before hydrogens are added): Asn and Gln amides and
  His rings turned 180 degrees where that makes better hydrogen bonds
  (X-ray data can't tell their N from O or C), His tautomers (HID or HIE;
  HIP when charged) chosen at the same time. Interacting groups are
  optimized together: every combination for small clusters, one group
  at a time (until nothing changes) for big ones.
- rotatable hydrogens (after hydrogens are added): Ser, Thr and Tyr OH,
  Cys SH and Lys NH3+ turned, waters reoriented, for hydrogen bonds and
  away from clashes; one group at a time, twice over.

Scores are simple counts. Heavy atoms: a donor and an acceptor 2.4 to 3.3
A apart +1 (to 0 at 3.6 A); two donors or two acceptors closer than 3.2 A
-1 (their hydrogens or lone pairs would clash), any pair closer than
2.4 A -1. A flip must win by FLIP_MARGIN. Hydrogens: H to an acceptor 1.5
to 2.5 A with the donor-H-acceptor angle 120 degrees or more +1; H closer
than 1.6 A to another H, or than 2.2 A to a heavy atom that isn't an
acceptor, -1.
'''

import itertools
import math

import numpy as np

DONOR, ACCEPTOR = 1, 2
BOTH = DONOR | ACCEPTOR
WATER = {'HOH', 'WAT', 'H2O', 'DOD', 'TIP3', 'SOL'}
AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
}
SIDE_CHAIN = {
    'SER': {'OG': BOTH}, 'THR': {'OG1': BOTH}, 'TYR': {'OH': BOTH},
    'ASN': {'OD1': ACCEPTOR, 'ND2': DONOR},
    'GLN': {'OE1': ACCEPTOR, 'NE2': DONOR},
    'ASP': {'OD1': ACCEPTOR, 'OD2': ACCEPTOR},
    'GLU': {'OE1': ACCEPTOR, 'OE2': ACCEPTOR},
    'LYS': {'NZ': DONOR}, 'ARG': {'NE': DONOR, 'NH1': DONOR, 'NH2': DONOR},
    'TRP': {'NE1': DONOR}, 'HIS': {'ND1': BOTH, 'NE2': BOTH},
    'MET': {'SD': ACCEPTOR}, 'CYS': {'SG': BOTH},
}
VARIANT_ROLES = {
    'ASH': {'OD1': BOTH, 'OD2': BOTH}, 'GLH': {'OE1': BOTH, 'OE2': BOTH},
    'LYN': {'NZ': BOTH}, 'HID': {'ND1': DONOR, 'NE2': ACCEPTOR},
    'HIE': {'ND1': ACCEPTOR, 'NE2': DONOR}, 'HIP': {'ND1': DONOR,
                                                   'NE2': DONOR},
    'CYX': {'SG': 0},
}
# flip groups: atom pairs whose positions swap
FLIPS = {'ASN': [('OD1', 'ND2')], 'GLN': [('OE1', 'NE2')],
         'HIS': [('ND1', 'CD2'), ('CE1', 'NE2')]}
FLIP_MARGIN = 0.5
EXHAUSTIVE_LIMIT = 4096
REACH = 3.6
# rotatable hydrogens: residue, parent, axis atom, hydrogens, symmetry
ROTORS = [('SER', 'OG', 'CB', ('HG',), 1), ('THR', 'OG1', 'CB', ('HG1',), 1),
          ('TYR', 'OH', 'CZ', ('HH',), 1), ('CYS', 'SG', 'CB', ('HG',), 1),
          ('LYS', 'NZ', 'CE', ('HZ1', 'HZ2', 'HZ3'), 3)]
WATER_ORIENTATIONS = 300


def _xyz(positions):
    from openmm import unit
    return np.array(positions.value_in_unit(unit.angstrom), dtype=float)


def _positions(xyz):
    from openmm import Vec3, unit
    return unit.Quantity([Vec3(*map(float, v / 10.0)) for v in xyz],
                         unit.nanometer)


def roles(topology, variants=None):
    '''
    :param variants: {residue: protonation state (ASH, HID, ...)}
    :return: {atom index: DONOR, ACCEPTOR or BOTH} of heavy atoms
    '''
    variants = variants or {}
    result = {}
    for residue in topology.residues():
        atoms = {a.name: a for a in residue.atoms()}
        name = residue.name
        if name in WATER:
            for atom in atoms.values():
                if atom.element is not None and atom.element.symbol == 'O':
                    result[atom.index] = BOTH
            continue
        if name in AMINO_ACIDS or name in ('ACE', 'NME'):
            if 'N' in atoms and name != 'PRO':
                result[atoms['N'].index] = DONOR
            for o in ('O', 'OXT'):
                if o in atoms:
                    result[atoms[o].index] = ACCEPTOR
            table = dict(SIDE_CHAIN.get(name, {}))
            table.update(VARIANT_ROLES.get(variants.get(residue), {}))
            for atom_name, role in table.items():
                if atom_name in atoms and role:
                    result[atoms[atom_name].index] = role
            continue
        for atom in atoms.values():  # ligands etc.: unknown, either way
            if atom.element is not None and atom.element.symbol in 'NO':
                result[atom.index] = BOTH
    return result


def pair_score(role1, role2, d):
    if d > REACH:
        return 0.0
    if d < 2.4:
        return -1.0
    if (role1 & DONOR and role2 & ACCEPTOR) or \
            (role1 & ACCEPTOR and role2 & DONOR):
        return 1.0 if d <= 3.3 else (REACH - d) / (REACH - 3.3)
    return -1.0 if d < 3.2 else 0.0


class _Group:
    '''
    A flippable Asn, Gln or His; states: (flipped, His tautomer)
    '''

    def __init__(self, residue, atoms, his_charged, flippable=True,
                 tautomer=None):
        self.residue = residue
        self.atoms = atoms  # name: index
        if residue.name == 'HIS':
            tautomers = [tautomer] if tautomer else \
                ['HIP'] if his_charged else ['HIE', 'HID']
        else:
            tautomers = [None]
        self.states = [(flip, t) for flip in ((False, True) if flippable
                                              else (False,))
                       for t in tautomers]
        self.mobile = {atoms[name] for pair in FLIPS[residue.name]
                       for name in pair}

    def sites(self, xyz, state):
        '''
        :return: [(position, role)] of the polar atoms in a state
        '''
        flip, tautomer = state
        a = self.atoms
        if self.residue.name == 'HIS':
            nd1 = xyz[a['CD2']] if flip else xyz[a['ND1']]
            ne2 = xyz[a['CE1']] if flip else xyz[a['NE2']]
            role = VARIANT_ROLES[tautomer]
            return [(nd1, role['ND1']), (ne2, role['NE2'])]
        acceptor, donor = FLIPS[self.residue.name][0]
        if flip:
            return [(xyz[a[acceptor]], DONOR), (xyz[a[donor]], ACCEPTOR)]
        return [(xyz[a[acceptor]], ACCEPTOR), (xyz[a[donor]], DONOR)]


def optimize_flips(topology, positions, his_charged=None, variants=None,
                   unflippable=(), tautomers=None):
    '''
    Choose Asn/Gln/His flips and His tautomers.

    :param his_charged: {residue: True if the His is charged (HIP)}
    :param variants: {residue: protonation state} of other residues
    :param unflippable: residues not to flip (e.g. built ones: their
        orientation means nothing)
    :param tautomers: {His residue: tautomer} decided already (e.g. by a
        metal it binds)
    :return: (positions, {His residue: tautomer}, [(residue, 'flip')])
    '''
    from scipy.spatial import cKDTree

    his_charged = his_charged or {}
    fixed_tautomers = tautomers or {}
    xyz = _xyz(positions)
    groups = []
    for residue in topology.residues():
        if residue.name not in FLIPS:
            continue
        atoms = {a.name: a.index for a in residue.atoms()}
        needed = {n for pair in FLIPS[residue.name] for n in pair}
        if needed <= set(atoms):
            groups.append(_Group(residue, atoms,
                                 his_charged.get(residue, False),
                                 residue not in unflippable,
                                 fixed_tautomers.get(residue)))
    if not groups:
        return positions, {}, []

    role = roles(topology, variants)
    mobile = set().union(*(g.mobile for g in groups))
    env = [i for i in role if i not in mobile and role[i]]
    env_xyz = xyz[env] if env else np.zeros((0, 3))
    tree = cKDTree(env_xyz) if env else None
    residue_of = {a.index: a.residue for a in topology.atoms()}

    # each group against the fixed environment, per state
    env_score = []
    for g in groups:
        scores = []
        for state in g.states:
            total = -FLIP_MARGIN if state[0] else 0.0
            for position, r in g.sites(xyz, state):
                if tree is None or not r:
                    continue
                for k in tree.query_ball_point(position, REACH):
                    atom = env[k]
                    if residue_of[atom] is g.residue:
                        continue
                    total += pair_score(r, role[atom], float(
                        np.linalg.norm(env_xyz[k] - position)))
            scores.append(total)
        env_score.append(scores)

    # groups against each other
    centers = np.array([xyz[list(g.mobile)].mean(axis=0) for g in groups])
    near = cKDTree(centers).query_pairs(2 * REACH + 3.0)
    pair = {}
    for i, j in near:
        table = np.zeros((len(groups[i].states), len(groups[j].states)))
        for si, sa in enumerate(groups[i].states):
            sites_a = groups[i].sites(xyz, sa)
            for sj, sb in enumerate(groups[j].states):
                table[si, sj] = sum(
                    pair_score(ra, rb, float(np.linalg.norm(pa - pb)))
                    for pa, ra in sites_a
                    for pb, rb in groups[j].sites(xyz, sb) if ra and rb)
        if table.any():
            pair[(i, j)] = table

    # clusters of interacting groups
    parent = list(range(len(groups)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in pair:
        parent[root(i)] = root(j)
    clusters = {}
    for i in range(len(groups)):
        clusters.setdefault(root(i), []).append(i)

    choice = [0] * len(groups)
    for members in clusters.values():
        def total(assign):
            score = sum(env_score[i][assign[i]] for i in members)
            for (i, j), table in pair.items():
                if i in assign and j in assign:
                    score += table[assign[i], assign[j]]
            return score

        sizes = [len(groups[i].states) for i in members]
        if math.prod(sizes) <= EXHAUSTIVE_LIMIT:
            best, best_score = None, None
            for combo in itertools.product(*(range(n) for n in sizes)):
                assign = dict(zip(members, combo))
                score = total(assign)
                if best_score is None or score > best_score + 1e-9:
                    best, best_score = assign, score
        else:
            best = {i: 0 for i in members}
            for _ in range(20):
                changed = False
                for i in members:
                    current = total(best)
                    for s in range(len(groups[i].states)):
                        trial = dict(best, **{i: s})
                        if total(trial) > current + 1e-9:
                            best, current, changed = trial, total(trial), True
                if not changed:
                    break
        for i, s in best.items():
            choice[i] = s

    flips = []
    tautomers = {}  # (the result)
    for g, s in zip(groups, choice):
        flip, tautomer = g.states[s]
        if tautomer:
            tautomers[g.residue] = tautomer
        if flip:
            flips.append(g.residue)
            for a, b in FLIPS[g.residue.name]:
                ia, ib = g.atoms[a], g.atoms[b]
                xyz[[ia, ib]] = xyz[[ib, ia]]
    return _positions(xyz), tautomers, flips


def _rotation(axis, angle):
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    c, s = math.cos(angle), math.sin(angle)
    return np.array([
        [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
        [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
        [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)]])


def _random_rotations(n, seed=1):
    '''
    n rotation matrices, roughly uniform (random unit quaternions)
    '''
    rng = np.random.default_rng(seed)
    q = rng.normal(size=(n, 4))
    q /= np.linalg.norm(q, axis=1)[:, None]
    w, x, y, z = q.T
    return np.stack([
        np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w),
                  2 * (x * z + y * w)], -1),
        np.stack([2 * (x * y + z * w), 1 - 2 * (x * x + z * z),
                  2 * (y * z - x * w)], -1),
        np.stack([2 * (x * z - y * w), 2 * (y * z + x * w),
                  1 - 2 * (x * x + y * y)], -1)], 1)


def orient_rotatable_hydrogens(topology, positions, variants=None):
    '''
    Turn OH, SH and NH3+ hydrogens and reorient waters for hydrogen bonds.

    :return: (positions, number of groups moved, number of groups)
    '''
    from scipy.spatial import cKDTree

    xyz = _xyz(positions)
    role = roles(topology, variants)
    hydrogens = np.array([a.element is not None and a.element.symbol == 'H'
                          for a in topology.atoms()])
    parent_of = {}
    for bond in topology.bonds():
        a, b = bond.atom1, bond.atom2
        if hydrogens[a.index]:
            parent_of[a.index] = b.index
        if hydrogens[b.index]:
            parent_of[b.index] = a.index

    # rotors: (hydrogen indices, parent index, candidate positions fn)
    rotors = []
    for residue in topology.residues():
        atoms = {a.name: a.index for a in residue.atoms()}
        for name, parent, axis, hs, symmetry in ROTORS:
            if residue.name == name and parent in atoms and axis in atoms \
                    and all(h in atoms for h in hs):
                rotors.append(('torsion', [atoms[h] for h in hs],
                               atoms[parent], atoms[axis], symmetry))
        if residue.name in WATER:
            o = [i for n, i in atoms.items()
                 if not hydrogens[i] and n.startswith('O')]
            hs = [i for i in atoms.values() if hydrogens[i]]
            if len(o) == 1 and len(hs) == 2:
                rotors.append(('water', hs, o[0], None, None))
    if not rotors:
        return positions, 0, 0

    angles = np.radians(np.arange(0, 360, 10))
    water_rotations = np.concatenate(
        [np.eye(3)[None], _random_rotations(WATER_ORIENTATIONS)])

    def candidates(rotor):
        kind, hs, p, axis, symmetry = rotor
        h0 = xyz[hs]
        if kind == 'water':
            return xyz[p] + np.einsum('rij,hj->rhi', water_rotations,
                                      h0 - xyz[p])
        steps = angles[angles < 2 * math.pi / symmetry]
        rel = h0 - xyz[p]
        return np.stack([xyz[p] + rel @ _rotation(xyz[p] - xyz[axis],
                                                 t).T for t in steps])

    def scores(rotor, cands, tree):
        kind, hs, p, axis, symmetry = rotor
        own = set(hs) | {p}
        # (the tree is from the start of the pass: hydrogens move a little)
        near = [k for k in tree.query_ball_point(xyz[p], 5.0)
                if k not in own]
        if not near:
            return np.zeros(len(cands))
        nxyz = xyz[near]
        n_h = hydrogens[near]
        n_acc = np.array([bool(role.get(k, 0) & ACCEPTOR) for k in near])
        d = np.linalg.norm(cands[:, :, None, :] - nxyz[None, None], axis=-1)
        # donor-H...acceptor angle
        to_h = cands - xyz[p]
        to_a = nxyz[None, None] - cands[:, :, None, :]
        cos = np.einsum('rhk,rhnk->rhn', -to_h, to_a) / (
            np.linalg.norm(to_h, axis=-1)[:, :, None] *
            np.linalg.norm(to_a, axis=-1) + 1e-9)
        hbond = (n_acc[None, None] & (d >= 1.5) & (d <= 2.5) &
                 (cos <= math.cos(math.radians(120))))
        clash_h = n_h[None, None] & (d < 1.6)
        clash_heavy = ~n_h[None, None] & ~n_acc[None, None] & (d < 2.2)
        return (hbond.sum(axis=(1, 2)) - clash_h.sum(axis=(1, 2)) -
                clash_heavy.sum(axis=(1, 2))).astype(float)

    start = xyz.copy()
    for _ in range(2):
        tree = cKDTree(xyz)
        for rotor in rotors:
            cands = candidates(rotor)
            s = scores(rotor, cands, tree)
            best = int(np.argmax(s))
            if s[best] > s[0]:  # candidate 0: as it is
                xyz[rotor[1]] = cands[best]
    moved = sum(1 for rotor in rotors
                if np.abs(xyz[rotor[1]] - start[rotor[1]]).max() > 0.1)
    return _positions(xyz), moved, len(rotors)
