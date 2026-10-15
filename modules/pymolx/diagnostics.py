'''
Structure diagnostics, like the Diagnostics tab of Maestro's Protein
Preparation Wizard (Structure Preparation dialog, command
"structure_diagnostics"). Each issue has a selection, so the dialog can
show it.

Checks: missing residues and chain breaks, missing side-chain atoms,
steric clashes (overlap of van der Waals radii of 0.4 A or more between
atoms more than 3 bonds apart, like MolProbity's "serious" clashes;
hydrogen bonds and metal contacts are not clashes), Ramachandran and
rotamer outliers (MolProbity's Top8000 reference data and cutoffs, in
data/pymolx/molprobity, CC BY 4.0, Richardson Lab), C-beta deviations,
cis peptides, D-amino acids, stretched or compressed bonds, alternate
locations, and whether hydrogens are present. "Scores" sums them up like
MolProbity: clashscore (serious clashes per 1000 atoms), Ramachandran
favored and outliers, rotamer outliers, C-beta deviations.
'''

import itertools

import math
from collections import namedtuple

from pymol import cmd

Issue = namedtuple('Issue', 'category severity residue detail selection')

# category: description, in report order
CATEGORIES = {
    'Scores': 'MolProbity-style summary',
    'Missing residues': 'residues known from the file but not modeled, '
                        'and chain breaks',
    'Missing atoms': 'heavy atoms missing from standard residues',
    'Steric clashes': 'van der Waals overlaps of 0.4 A or more',
    'Ramachandran outliers': 'phi/psi rare in the Top8000 reference data',
    'Rotamer outliers': 'side-chain conformations rare in the Top8000 data',
    'C-beta deviations': 'C-beta 0.25 A or more from its ideal position',
    'Cis peptides': 'peptide bonds with omega near 0',
    'Chirality': 'D-amino acids in a chain of L-amino acids',
    'Bond lengths': 'bonds much longer or shorter than covalent radii',
    'Alternate locations': 'residues with alternate conformations',
    'Hydrogens': 'hydrogens present or not',
}
SEVERITIES = ('error', 'warning', 'info')
# reported without a count (one line per item)
INFO_CATEGORIES = ('Scores', 'Hydrogens')

# MolProbity's Top8000 grids (data/pymolx/molprobity) and cutoffs
TOP8000 = '$PYMOL_DATA/pymolx/molprobity/top8000.npz'
FAVORED = 0.02
RAMA_ALLOWED = {'general': 0.0005, 'cispro': 0.002}  # others: 0.001
ROTAMER_ALLOWED = 0.003
CBETA_DEVIATION = 0.25
CHI = {  # side-chain dihedrals, atom names
    'ARG': ['N CA CB CG', 'CA CB CG CD', 'CB CG CD NE', 'CG CD NE CZ'],
    'ASN': ['N CA CB CG', 'CA CB CG OD1'],
    'ASP': ['N CA CB CG', 'CA CB CG OD1'],
    'CYS': ['N CA CB SG'],
    'GLN': ['N CA CB CG', 'CA CB CG CD', 'CB CG CD OE1'],
    'GLU': ['N CA CB CG', 'CA CB CG CD', 'CB CG CD OE1'],
    'HIS': ['N CA CB CG', 'CA CB CG ND1'],
    'ILE': ['N CA CB CG1', 'CA CB CG1 CD1'],
    'LEU': ['N CA CB CG', 'CA CB CG CD1'],
    'LYS': ['N CA CB CG', 'CA CB CG CD', 'CB CG CD CE', 'CG CD CE NZ'],
    'MET': ['N CA CB CG', 'CA CB CG SD', 'CB CG SD CE'],
    'PHE': ['N CA CB CG', 'CA CB CG CD1'],
    'PRO': ['N CA CB CG'],
    'SER': ['N CA CB OG'],
    'THR': ['N CA CB OG1'],
    'TRP': ['N CA CB CG', 'CA CB CG CD1'],
    'TYR': ['N CA CB CG', 'CA CB CG CD1'],
    'VAL': ['N CA CB CG1'],
}
_top8000 = None

HEAVY_ATOMS = {
    'ALA': 'CB', 'ARG': 'CB CG CD NE CZ NH1 NH2', 'ASN': 'CB CG OD1 ND2',
    'ASP': 'CB CG OD1 OD2', 'CYS': 'CB SG', 'GLN': 'CB CG CD OE1 NE2',
    'GLU': 'CB CG CD OE1 OE2', 'GLY': '', 'HIS': 'CB CG ND1 CD2 CE1 NE2',
    'ILE': 'CB CG1 CG2 CD1', 'LEU': 'CB CG CD1 CD2', 'LYS': 'CB CG CD CE NZ',
    'MET': 'CB CG SD CE', 'PHE': 'CB CG CD1 CD2 CE1 CE2 CZ',
    'PRO': 'CB CG CD', 'SER': 'CB OG', 'THR': 'CB OG1 CG2',
    'TRP': 'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
    'TYR': 'CB CG CD1 CD2 CE1 CE2 CZ OH', 'VAL': 'CB CG1 CG2',
}
HEAVY_ATOMS = {resn: {'N', 'CA', 'C', 'O'} | set(atoms.split())
               for resn, atoms in HEAVY_ATOMS.items()}

# van der Waals radii (Bondi; O and carbonyl C as MolProbity's probe, so
# O(i)...C(i+1) of ordinary peptides doesn't count) and covalent radii,
# Angstrom
VDW = {'H': 1.10, 'C': 1.70, 'N': 1.55, 'O': 1.40, 'S': 1.80, 'P': 1.80,
       'F': 1.47, 'CL': 1.75, 'BR': 1.85, 'I': 1.98, 'SE': 1.90}
CARBONYL_C = 1.65
COVALENT = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05,
            'P': 1.07, 'F': 0.57, 'CL': 1.02, 'BR': 1.20, 'I': 1.39,
            'SE': 1.20}
POLAR = {'N', 'O'}
WATER = {'HOH', 'WAT', 'H2O', 'DOD', 'TIP3', 'SOL'}
HBOND = {'N', 'O', 'S'}  # donors (of polar H) and acceptors
POLAR_H = 1.0            # van der Waals radius of polar H (MolProbity)
CLASH_OVERLAP = 0.4
SEVERE_OVERLAP = 0.9
CHAIN_BREAK = 2.0   # C-N distance, Angstrom


def _label(atom):
    chain = atom['chain'] or '-'
    return '%s/%s %s' % (chain, atom['resn'], atom['resi'])


def _residue_key(atom):
    return (atom['model'], atom['segi'], atom['chain'], atom['resi'])


def _byres(atoms):
    groups = {}
    for atom in atoms:
        groups.setdefault(atom['model'], []).append(str(atom['index']))
    return ' or '.join('byres (%%%s and index %s)' % (model, '+'.join(idx))
                       for model, idx in groups.items())


def _dihedral(p0, p1, p2, p3):
    b0 = [a - b for a, b in zip(p0, p1)]
    b1 = [a - b for a, b in zip(p2, p1)]
    b2 = [a - b for a, b in zip(p3, p2)]
    norm = math.sqrt(sum(c * c for c in b1))
    b1 = [c / norm for c in b1]

    def dot(u, v):
        return sum(a * b for a, b in zip(u, v))

    def cross(u, v):
        return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0]]

    v = [a - dot(b0, b1) * b for a, b in zip(b0, b1)]
    w = [a - dot(b2, b1) * b for a, b in zip(b2, b1)]
    x = dot(v, w)
    y = dot(cross(b1, v), w)
    return math.degrees(math.atan2(y, x))


def _read(selection, state, _self):
    '''
    :return: (atoms, bonds): atom dicts with coordinates and index pairs
    '''
    model = _self.get_model('(%s) and present' % selection, state=state)
    atoms = []
    for i, a in enumerate(model.atom):
        atoms.append({
            'i': i, 'index': a.index, 'model': a.model, 'segi': a.segi,
            'chain': a.chain, 'resn': a.resn, 'resi': a.resi,
            'resv': a.resi_number, 'name': a.name,
            'elem': (a.symbol or '').upper(), 'alt': a.alt,
            'hetatm': a.hetatm, 'q': a.q, 'coord': tuple(a.coord)})
    bonds = [tuple(b.index) for b in model.bond]
    return atoms, bonds


def _residues(atoms):
    '''
    :return: {residue key: {atom name: atom}}, of residues with alternate
        locations the first one (like most validation reports)
    '''
    first = {}
    for atom in atoms:
        if atom['alt']:
            key = _residue_key(atom)
            first[key] = min(first.get(key, atom['alt']), atom['alt'])
    residues = {}
    for atom in atoms:
        key = _residue_key(atom)
        if atom['alt'] and atom['alt'] != first[key]:
            continue
        residues.setdefault(key, {})[atom['name']] = atom
    return residues


def _protein_chains(residues):
    '''
    :return: {(model, segi, chain): [residue atoms dicts in order]}
    '''
    chains = {}
    for key, named in residues.items():
        if 'CA' in named and 'N' in named and 'C' in named:
            chains.setdefault(key[:3], []).append(named)
    return chains


def check_missing_residues(selection, atoms, residues, _self):
    issues = []
    # known but unmodeled (no coordinates)
    missing = []
    _self.iterate('(%s) and polymer and not present' % selection,
                  'missing.append((model, chain, resi, resv, resn))',
                  space={'missing': missing})
    runs = []
    for model, chain, resi, resv, resn in missing:
        if runs and runs[-1][0] == (model, chain) and \
                resv - runs[-1][2] <= 1:
            runs[-1][2] = resv
            runs[-1][3].append(resi)
        else:
            runs.append([(model, chain), resv, resv, [resi]])
    for (model, chain), first, last, resis in runs:
        span = '%s' % first if first == last else '%s-%s' % (first, last)
        sel = '%%%s and chain %s and resi %s' % (
            model, chain or '""', '+'.join(resis))
        issues.append(Issue(
            'Missing residues', 'warning', '%s/%s' % (chain or '-', span),
            '%d residue%s without coordinates' % (
                len(resis), 's' if len(resis) > 1 else ''),
            'byres (%s) or (%s)' % (sel, 'neighbor (%s)' % sel)))

    # chain breaks between modeled residues (not where residues are known
    # to be missing)
    unmodeled = {}
    for model, chain, resi, resv, resn in missing:
        unmodeled.setdefault((model, chain), set()).add(resv)
    for chain_residues in _protein_chains(residues).values():
        for a, b in zip(chain_residues, chain_residues[1:]):
            d = math.dist(a['C']['coord'], b['N']['coord'])
            known = unmodeled.get((a['CA']['model'], a['CA']['chain']), ())
            if any(a['CA']['resv'] < v < b['CA']['resv'] for v in known):
                continue
            if d > CHAIN_BREAK:
                issues.append(Issue(
                    'Missing residues', 'warning', '%s - %s' % (
                        _label(a['CA']), _label(b['CA'])),
                    'chain break, C-N %.1f A' % d,
                    _byres([a['CA'], b['CA']])))
    return issues


def check_missing_atoms(residues):
    issues = []
    for named in residues.values():
        atom = next(iter(named.values()))
        expected = HEAVY_ATOMS.get(atom['resn'])
        if expected is None or atom['hetatm']:
            continue
        absent = expected - set(named)
        if absent:
            issues.append(Issue(
                'Missing atoms', 'error', _label(atom),
                'missing %s' % ' '.join(sorted(absent)), _byres([atom])))
    return issues


def check_clashes(atoms, bonds):
    import numpy as np
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        return [Issue('Steric clashes', 'info', 'all',
                      'not checked: needs scipy (pip install scipy)', '')], None

    usable = [a for a in atoms if a['elem'] in VDW]
    if len(usable) < 2:
        return [], 0
    neighbors = {a['i']: set() for a in atoms}
    for i, j in bonds:
        neighbors[i].add(j)
        neighbors[j].add(i)

    def within_bonds(i, n):
        near = set(neighbors[i])
        for _ in range(n - 1):
            near |= {k for j in near for k in neighbors[j]}
        return near

    xyz = np.array([a['coord'] for a in usable])
    tree = cKDTree(xyz)
    pairs = tree.query_pairs(2 * max(VDW.values()) - CLASH_OVERLAP,
                             output_type='ndarray')
    hydrogen_parent = {}
    for i, j in bonds:
        if atoms[i]['elem'] == 'H':
            hydrogen_parent[i] = atoms[j]['elem']
        if atoms[j]['elem'] == 'H':
            hydrogen_parent[j] = atoms[i]['elem']

    def radius(atom):
        if atom['elem'] == 'H' and \
                hydrogen_parent.get(atom['i']) in HBOND:
            return POLAR_H
        if atom['name'] == 'C' and atom['elem'] == 'C':
            return CARBONYL_C
        return VDW[atom['elem']]

    excluded = {}
    worst = {}
    clashing_pairs = 0  # for the clashscore, MolProbity's way
    for p, q in pairs:
        a, b = usable[p], usable[q]
        if a['alt'] and b['alt'] and a['alt'] != b['alt']:
            continue  # alternate locations never coexist
        d = float(np.linalg.norm(xyz[p] - xyz[q]))
        overlap = radius(a) + radius(b) - d
        if overlap < CLASH_OVERLAP:
            continue
        # hydrogen bonds: polar H to N/O/S, and N/O to N/O
        ea, eb = a['elem'], b['elem']
        if (ea == 'H' and hydrogen_parent.get(a['i']) in HBOND and
                eb in HBOND) or (eb == 'H' and ea in HBOND and
                                 hydrogen_parent.get(b['i']) in HBOND):
            continue
        if ea in POLAR and eb in POLAR and d >= 2.5:
            continue
        if a['i'] not in excluded:
            excluded[a['i']] = (within_bonds(a['i'], 3),
                                within_bonds(a['i'], 4))
        if b['i'] in excluded[a['i']][0]:
            continue
        # the clashscore (like probe -4H): not water, occupancy 0.33 or
        # more, atoms with a hydrogen 4 bonds apart don't count
        if not (a['resn'] in WATER or b['resn'] in WATER or
                a['q'] < 0.33 or b['q'] < 0.33 or
                ('H' in (ea, eb) and b['i'] in excluded[a['i']][1])):
            clashing_pairs += 1
        key = tuple(sorted((_residue_key(a), _residue_key(b))))
        if key not in worst or overlap > worst[key][0]:
            worst[key] = (overlap, d, a, b)

    issues = []
    for overlap, d, a, b in sorted(worst.values(), key=lambda w: -w[0]):
        issues.append(Issue(
            'Steric clashes',
            'error' if overlap >= SEVERE_OVERLAP else 'warning',
            '%s - %s' % (_label(a), _label(b)),
            '%s-%s %.2f A, overlap %.2f A' % (a['name'], b['name'], d,
                                              overlap),
            _byres([a, b])))
    return issues, clashing_pairs


def top8000(_self=cmd):
    '''
    :return: {table name: array}, empty if the data file is missing
    '''
    global _top8000
    if _top8000 is None:
        import numpy as np
        try:
            with np.load(_self.exp_path(TOP8000)) as data:
                _top8000 = {key: data[key] for key in data.files}
        except OSError:
            _top8000 = {}
    return _top8000


def lookup(tables, name, angles):
    '''
    Value of a Top8000 table at the angles (degrees): linear
    interpolation between bin centers, all axes wrapping (like
    MolProbity's NDimTable)
    '''
    grid = tables[name]
    dims = tables[name + '_dims']
    index = []
    fraction = []
    for angle, (lo, hi, n) in zip(angles, dims):
        n = int(n)
        position = ((angle - lo) % (hi - lo)) / ((hi - lo) / n) - 0.5
        first = math.floor(position)
        index.append((first % n, (first + 1) % n))
        fraction.append(position - first)
    value = 0.0
    for corner in itertools.product((0, 1), repeat=len(index)):
        weight = 1.0
        for k, c in enumerate(corner):
            weight *= fraction[k] if c else 1 - fraction[k]
        if weight:
            value += weight * float(grid[tuple(
                index[k][c] for k, c in enumerate(corner))])
    return value


def _phi_psi_triples(residues):
    '''
    :return: [(previous, current, next)] residue atoms of bonded amino
        acids: the current one with N, CA and C, the previous one with C,
        the next one with N (as much as phi and psi need)
    '''
    chains = {}
    for key, named in residues.items():
        if {'N', 'CA', 'C'} & set(named):
            chains.setdefault(key[:3], []).append(named)
    triples = []
    for chain_residues in chains.values():
        for prev, cur, nxt in zip(chain_residues, chain_residues[1:],
                                  chain_residues[2:]):
            if not ({'N', 'CA', 'C'} <= set(cur) and 'C' in prev and
                    'N' in nxt):
                continue
            if math.dist(prev['C']['coord'], cur['N']['coord']) > \
                    CHAIN_BREAK or math.dist(
                        cur['C']['coord'], nxt['N']['coord']) > CHAIN_BREAK:
                continue
            if cur['CA']['resn'] in HEAVY_ATOMS:
                triples.append((prev, cur, nxt))
    return triples


def rama_class(prev, cur, nxt):
    '''
    MolProbity's Ramachandran residue class
    '''
    resn = cur['CA']['resn']
    if resn == 'GLY':
        return 'glycine'
    if resn == 'PRO':
        if 'CA' not in prev:
            return 'transpro'  # (like MolProbity without omega)
        omega = _dihedral(prev['CA']['coord'], prev['C']['coord'],
                          cur['N']['coord'], cur['CA']['coord'])
        return 'cispro' if -90 < omega < 90 else 'transpro'
    if nxt['N']['resn'] == 'PRO':
        return 'prepro'
    if resn in ('ILE', 'VAL'):
        return 'ileval'
    return 'general'


def check_ramachandran(residues, tables):
    '''
    :return: (issues, {'favored': n, 'allowed': n, 'outlier': n})
    '''
    issues = []
    counts = {'favored': 0, 'allowed': 0, 'outlier': 0}
    for prev, cur, nxt in _phi_psi_triples(residues):
        phi = _dihedral(prev['C']['coord'], cur['N']['coord'],
                        cur['CA']['coord'], cur['C']['coord'])
        psi = _dihedral(cur['N']['coord'], cur['CA']['coord'],
                        cur['C']['coord'], nxt['N']['coord'])
        kind = rama_class(prev, cur, nxt)
        value = lookup(tables, 'rama_' + kind, (phi, psi))
        if value >= FAVORED:
            counts['favored'] += 1
        elif value >= RAMA_ALLOWED.get(kind, 0.001):
            counts['allowed'] += 1
        else:
            counts['outlier'] += 1
            issues.append(Issue(
                'Ramachandran outliers', 'warning', _label(cur['CA']),
                'phi %.0f, psi %.0f (%s, %.2f%%)' % (
                    phi, psi, kind, 100 * value), _byres([cur['CA']])))
    return issues, counts


def check_rotamers(residues, tables):
    '''
    :return: (issues, number of residues evaluated)
    '''
    issues = []
    evaluated = 0
    for named in residues.values():
        atom = next(iter(named.values()))
        chis = CHI.get(atom['resn'])
        if not chis or atom['hetatm']:
            continue
        needed = {name for chi in chis for name in chi.split()}
        if not needed <= set(named):
            continue
        angles = [_dihedral(*(named[name]['coord'] for name in chi.split()))
                  for chi in chis]
        table = 'rota_' + {'PHE': 'phetyr', 'TYR': 'phetyr'}.get(
            atom['resn'], atom['resn'].lower())
        value = lookup(tables, table, angles)
        evaluated += 1
        if value < ROTAMER_ALLOWED:
            issues.append(Issue(
                'Rotamer outliers', 'warning', _label(atom),
                '%s (%.2f%%)' % (', '.join(
                    'chi%d %.0f' % (k, a % 360)
                    for k, a in enumerate(angles, 1)), 100 * value),
                _byres([atom])))
    return issues, evaluated


def ideal_cbeta(n, ca, c):
    '''
    C-beta position from the backbone, ideal L-amino acid geometry
    '''
    b = [p - q for p, q in zip(ca, n)]
    d = [p - q for p, q in zip(c, ca)]
    a = [b[1] * d[2] - b[2] * d[1], b[2] * d[0] - b[0] * d[2],
         b[0] * d[1] - b[1] * d[0]]
    return [-0.58273431 * x + 0.56802827 * y - 0.54067466 * z + w
            for x, y, z, w in zip(a, b, d, ca)]


def check_cbeta(residues):
    issues = []
    for named in residues.values():
        if not {'N', 'CA', 'C', 'CB'} <= set(named):
            continue
        atom = named['CA']
        if atom['resn'] not in HEAVY_ATOMS or atom['resn'] == 'GLY':
            continue
        ideal = ideal_cbeta(named['N']['coord'], atom['coord'],
                            named['C']['coord'])
        deviation = math.dist(ideal, named['CB']['coord'])
        if deviation >= CBETA_DEVIATION:
            issues.append(Issue(
                'C-beta deviations', 'warning', _label(atom),
                'C-beta %.2f A from its ideal position' % deviation,
                _byres([atom])))
    return issues


def check_backbone(residues):
    '''
    Coarse phi/psi check, without the reference data: clearly disallowed
    regions only
    '''
    issues = []
    for chain_residues in _protein_chains(residues).values():
        for prev, cur, nxt in zip(chain_residues, chain_residues[1:],
                                  chain_residues[2:]):
            if math.dist(prev['C']['coord'], cur['N']['coord']) > \
                    CHAIN_BREAK or math.dist(
                        cur['C']['coord'], nxt['N']['coord']) > CHAIN_BREAK:
                continue
            phi = _dihedral(prev['C']['coord'], cur['N']['coord'],
                            cur['CA']['coord'], cur['C']['coord'])
            psi = _dihedral(cur['N']['coord'], cur['CA']['coord'],
                            cur['C']['coord'], nxt['N']['coord'])
            resn = cur['CA']['resn']
            if resn == 'GLY':
                continue
            if resn == 'PRO':
                outlier = not -120 <= phi <= -30
            else:
                left = 25 <= phi <= 110 and -70 <= psi <= 120
                outlier = phi > 0 and not left
            if outlier:
                issues.append(Issue(
                    'Ramachandran outliers', 'warning', _label(cur['CA']),
                    'phi %.0f, psi %.0f' % (phi, psi), _byres([cur['CA']])))
    return issues


def check_cis_peptides(residues):
    issues = []
    for chain_residues in _protein_chains(residues).values():
        for a, b in zip(chain_residues, chain_residues[1:]):
            if math.dist(a['C']['coord'], b['N']['coord']) > CHAIN_BREAK:
                continue
            omega = _dihedral(a['CA']['coord'], a['C']['coord'],
                              b['N']['coord'], b['CA']['coord'])
            if abs(omega) < 30:
                proline = b['CA']['resn'] == 'PRO'
                issues.append(Issue(
                    'Cis peptides', 'info' if proline else 'warning',
                    '%s - %s' % (_label(a['CA']), _label(b['CA'])),
                    '%s, omega %.0f' % ('cis-proline' if proline else
                                        'cis peptide', omega),
                    _byres([a['CA'], b['CA']])))
    return issues


def _chirality(named):
    n, ca, c, cb = (named[k]['coord'] for k in ('N', 'CA', 'C', 'CB'))
    u = [a - b for a, b in zip(n, ca)]
    v = [a - b for a, b in zip(c, ca)]
    w = [a - b for a, b in zip(cb, ca)]
    return (u[0] * (v[1] * w[2] - v[2] * w[1]) -
            u[1] * (v[0] * w[2] - v[2] * w[0]) +
            u[2] * (v[0] * w[1] - v[1] * w[0]))


def check_chirality(residues):
    issues = []
    for named in residues.values():
        if not {'N', 'CA', 'C', 'CB'} <= set(named):
            continue
        atom = named['CA']
        if atom['resn'] not in HEAVY_ATOMS or atom['resn'] == 'GLY':
            continue
        # L-amino acids: (N-CA) . ((C-CA) x (CB-CA)) > 0
        if _chirality(named) < -0.5:
            issues.append(Issue('Chirality', 'warning', _label(atom),
                                'D-amino acid (CA)', _byres([atom])))
    return issues


def check_bonds(atoms, bonds):
    issues = []
    for i, j in bonds:
        a, b = atoms[i], atoms[j]
        if a['elem'] not in COVALENT or b['elem'] not in COVALENT:
            continue
        if a['alt'] and b['alt'] and a['alt'] != b['alt']:
            continue
        d = math.dist(a['coord'], b['coord'])
        ideal = COVALENT[a['elem']] + COVALENT[b['elem']]
        if d > ideal + 0.45 or d < ideal - 0.4:
            issues.append(Issue(
                'Bond lengths', 'error' if d > ideal + 1.0 else 'warning',
                _label(a) if _residue_key(a) == _residue_key(b) else
                '%s - %s' % (_label(a), _label(b)),
                '%s-%s %.2f A (expected about %.2f)' % (
                    a['name'], b['name'], d, ideal), _byres([a, b])))
    return issues


def check_alternates(atoms):
    alts = {}
    for atom in atoms:
        if atom['alt']:
            alts.setdefault(_residue_key(atom), (atom, set()))[1].add(
                atom['alt'])
    return [Issue('Alternate locations', 'info', _label(atom),
                  'alternates %s' % ', '.join(sorted(letters)),
                  _byres([atom]))
            for atom, letters in alts.values()]


def check_hydrogens(atoms, selection):
    protein = [a for a in atoms if a['resn'] in HEAVY_ATOMS]
    hydrogens = sum(a['elem'] == 'H' for a in protein)
    if not protein:
        return []
    if hydrogens:
        return [Issue('Hydrogens', 'info', 'all', '%d hydrogens on protein '
                      'residues' % hydrogens, '(%s) and hydro' % selection)]
    return [Issue('Hydrogens', 'warning', 'all', 'no hydrogens: protonate '
                  'before minimizing or scoring', '(%s)' % selection)]


def run(selection='all', state=-1, *, _self=cmd):
    '''
    :return: list of Issue
    '''
    state = int(state)
    if state < 1:
        state = _self.get_state()
    atoms, bonds = _read(selection, state, _self)
    if not atoms:
        raise ValueError('no atoms with coordinates in "%s"' % selection)
    residues = _residues(atoms)
    issues = []
    issues += check_missing_residues(selection, atoms, residues, _self)
    issues += check_missing_atoms(residues)
    clashes, clashing_pairs = check_clashes(atoms, bonds)
    issues += clashes
    tables = top8000(_self)
    rama = rotamers = None
    if tables:
        found, rama = check_ramachandran(residues, tables)
        issues += found
        found, evaluated = check_rotamers(residues, tables)
        issues += found
        rotamers = (len(found), evaluated)
    else:
        issues += check_backbone(residues)
    cbeta = check_cbeta(residues)
    issues += cbeta
    issues += check_cis_peptides(residues)
    issues += check_chirality(residues)
    issues += check_bonds(atoms, bonds)
    issues += check_alternates(atoms)
    issues += check_hydrogens(atoms, selection)
    hydrogens = any(a['elem'] == 'H' for a in atoms)
    counted = sum(a['resn'] not in WATER for a in atoms)
    return _scores(counted, clashing_pairs, hydrogens, rama, rotamers,
                   len(cbeta), selection) + issues


def _scores(n_atoms, clashing_pairs, hydrogens, rama, rotamers, cbeta,
            selection):
    '''
    MolProbity-style summary lines (category "Scores")
    '''
    scores = []

    def add(detail):
        scores.append(Issue('Scores', 'info', 'all', detail, ''))

    if clashing_pairs is not None:
        add('clashscore %.1f (serious clashes per 1000 atoms, waters '
            'excluded%s)' % (1000.0 * clashing_pairs / max(1, n_atoms),
            '' if hydrogens else '; no hydrogens, not comparable'))
    if rama:
        total = sum(rama.values())
        if total:
            add('Ramachandran favored %.1f%%, outliers %.1f%% (%d '
                'residues)' % (100.0 * rama['favored'] / total,
                               100.0 * rama['outlier'] / total, total))
    if rotamers and rotamers[1]:
        add('rotamer outliers %.1f%% (%d of %d)' % (
            100.0 * rotamers[0] / rotamers[1], rotamers[0], rotamers[1]))
    add('C-beta deviations: %d' % cbeta)
    return scores


def summary(issues):
    '''
    :return: {category: count}, in CATEGORIES order (all categories)
    '''
    counts = {category: 0 for category in CATEGORIES}
    for issue in issues:
        if issue.category not in INFO_CATEGORIES:
            counts[issue.category] += 1
    return counts


def format_report(issues, limit=10):
    lines = ['Diagnostics']
    counts = summary(issues)
    for category in CATEGORIES:
        found = [i for i in issues if i.category == category]
        if category in INFO_CATEGORIES:
            for issue in found:
                lines.append('  %s: %s' % (category, issue.detail))
            continue
        lines.append('  %s: %s' % (category, counts[category] or 'none'))
        for issue in found[:limit]:
            lines.append('    %-7s %s: %s' % (issue.severity, issue.residue,
                                              issue.detail))
        if len(found) > limit:
            lines.append('    ... %d more' % (len(found) - limit))
    return '\n'.join(lines)


def structure_diagnostics(selection='all', state=-1, quiet=0, *, _self=cmd):
    '''
DESCRIPTION

    "structure_diagnostics" reports problems of a structure: missing
    residues and atoms, chain breaks, steric clashes, Ramachandran and
    rotamer outliers (MolProbity's Top8000 data), C-beta deviations, cis
    peptides, D-amino acids, unusual bond lengths, alternate locations
    and missing hydrogens, with MolProbity-style scores (clashscore, ...).

USAGE

    structure_diagnostics [ selection [, state ]]

SEE ALSO

    prepare_structure
    '''
    issues = run(selection, state, _self=_self)
    if not int(quiet):
        print(format_report(issues))
    return issues


def extend(_self):
    _self.extend('structure_diagnostics', structure_diagnostics)
    _self.structure_diagnostics = structure_diagnostics
