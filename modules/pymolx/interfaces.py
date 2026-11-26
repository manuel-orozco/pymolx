'''
PISA-like interface analysis (pymolx).

For each pair of polymer chains in contact within one object:

- interface area: half of the solvent-accessible area buried on complex
  formation (PISA's definition)
- dG_int: solvation free energy gain on forming the interface, from
  atomic solvation parameters fitted to PISA (see ASP below)
- hydrogen bonds, salt bridges and disulfide bonds across the interface,
  with distance criteria tuned to match PISA's counts (no hydrogens
  needed)
- interface residues with buried area, dG contribution and bond flags

No P-value, CSS or assembly prediction (these need PISA itself).

Also the epitope (the residues of chosen target chains buried by the
other chains) and the proximity profile: how close each residue of the
targets comes to the other chains.
'''

import itertools
from collections import defaultdict

import numpy

from pymol import cmd

# donor-acceptor distance; distance-only, so shorter than PISA's own
# criterion: 3.4 A gives counts with near zero bias against PISA
# (tools/pymolx/calibrate_interfaces.py)
HBOND_CUTOFF = 3.4
SALT_BRIDGE_CUTOFF = 4.0   # as PISA
DISULFIDE_CUTOFF = 2.5
CONTACT_CUTOFF = 5.0       # chains closer than this may form an interface
SASA_DOT_DENSITY = 3

# atomic solvation parameters (cal/mol/A^2), Eisenberg & McLachlan 1986
ASP_EISENBERG = {'C': 16.0, 'N': -6.0, 'O': -6.0, 'O-': -24.0, 'N+': -50.0,
                 'S': 21.0}

# atomic solvation parameters (cal/mol/A^2) fitted to PDBe PISA's
# int_solv_en: 179 interfaces of 65 PDB entries (2026-10-06,
# tools/pymolx/calibrate_interfaces.py). Fitted on half of the entries,
# the other half had a mean absolute error of 0.66 kcal/mol (R^2 0.988),
# against 3.9 kcal/mol (R^2 0.59) for ASP_EISENBERG.
ASP_PISA_FIT = {'C': 16.3, 'N': -19.1, 'O': -4.4, 'O-': -10.6, 'N+': -26.0,
                'S': 45.0}

ASP = ASP_PISA_FIT

METHOD_NOTE = (
    'Computed locally by pymolx with PISA-like definitions. Interface area: '
    'half the buried solvent-accessible area. \u0394iG: solvation free '
    'energy gain on forming the interface, with atomic solvation parameters '
    'fitted to PISA (typically within 1 kcal/mol of PISA), without hydrogen '
    'bond and salt bridge contributions. No P-value, CSS or assembly '
    'prediction.')

CARBOXYLATE = {('ASP', 'OD1'), ('ASP', 'OD2'), ('GLU', 'OE1'),
               ('GLU', 'OE2')}
CATIONIC = {('LYS', 'NZ'), ('ARG', 'NE'), ('ARG', 'NH1'), ('ARG', 'NH2')}

# heavy atom donors/acceptors of standard amino acids (backbone N/O added
# below); other residues: any N or O counts as both
DONORS = {('ARG', 'NE'), ('ARG', 'NH1'), ('ARG', 'NH2'), ('ASN', 'ND2'),
          ('GLN', 'NE2'), ('HIS', 'ND1'), ('HIS', 'NE2'), ('LYS', 'NZ'),
          ('SER', 'OG'), ('THR', 'OG1'), ('TYR', 'OH'), ('TRP', 'NE1')}
ACCEPTORS = {('ASP', 'OD1'), ('ASP', 'OD2'), ('GLU', 'OE1'), ('GLU', 'OE2'),
             ('ASN', 'OD1'), ('GLN', 'OE1'), ('HIS', 'ND1'), ('HIS', 'NE2'),
             ('SER', 'OG'), ('THR', 'OG1'), ('TYR', 'OH')}
AMINO_ACIDS = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
               'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
               'TYR', 'VAL', 'MSE', 'SEC'}


class Atom:
    __slots__ = ('index', 'chain', 'resn', 'resi', 'name', 'elem', 'alt',
                 'xyz', 'asa', 'bsa')

    def __init__(self, index, chain, resn, resi, name, elem, alt):
        self.index = index
        self.chain = chain
        self.resn = resn
        self.resi = resi
        self.name = name
        self.elem = elem
        self.alt = alt
        self.xyz = None
        self.asa = 0.0  # in the isolated chain
        self.bsa = 0.0  # buried in this interface

    @property
    def key(self):
        return (self.chain, self.resi, self.resn, self.name, self.alt)

    @property
    def residue(self):
        return (self.chain, self.resi, self.resn)

    def asp_type(self):
        if self.elem == 'C':
            return 'C'
        if self.elem == 'S':
            return 'S'
        if (self.resn, self.name) in CARBOXYLATE or self.name == 'OXT':
            return 'O-'
        if (self.resn, self.name) in CATIONIC:
            return 'N+'
        if self.elem in ('N', 'O'):
            return self.elem
        return None

    def is_donor(self):
        if self.elem not in ('N', 'O'):
            return False
        if self.resn not in AMINO_ACIDS:
            return True
        if self.name == 'N':
            return self.resn != 'PRO'
        return (self.resn, self.name) in DONORS

    def is_acceptor(self):
        if self.elem not in ('N', 'O'):
            return False
        if self.resn not in AMINO_ACIDS:
            return True
        if self.name in ('O', 'OXT'):
            return True
        return (self.resn, self.name) in ACCEPTORS

    def is_anion(self):
        return (self.resn, self.name) in CARBOXYLATE or self.name == 'OXT'

    def is_cation(self):
        return (self.resn, self.name) in CATIONIC

    def label(self):
        return '%s/%s %s/%s' % (self.chain, self.resn, self.resi, self.name)


def chain_selection(obj, chain):
    return '%%%s and chain "%s" and polymer and not hydro' % (obj, chain)


def _load_atoms(sele, state, _self):
    atoms = []
    _self.iterate_state(
        state, sele,
        'atoms.append(Atom(index, chain, resn, resi, name, elem, alt))',
        space={'atoms': atoms, 'Atom': Atom})
    xyz = _self.get_coords(sele, state)
    for atom, coord in zip(atoms, xyz if xyz is not None else ()):
        atom.xyz = coord
    return atoms


def _isolated_areas(sele, state, _self):
    '''
    Per-atom solvent-accessible area of the selection on its own.

    :return: {atom key: area}
    '''
    tmp = _self.get_unused_name('_pymolx_iface')
    _self.create(tmp, sele, state, 1, zoom=0)
    try:
        _self.get_area(tmp, load_b=1)
        areas = {}
        _self.iterate(tmp, 'areas[(chain, resi, resn, name, alt)] = b',
                      space={'areas': areas})
        return areas
    finally:
        _self.delete(tmp)


def _pairs_within(atoms1, atoms2, cutoff):
    '''
    :return: list of (atom1, atom2, distance) closer than cutoff
    '''
    if not atoms1 or not atoms2:
        return []
    xyz1 = numpy.array([a.xyz for a in atoms1])
    xyz2 = numpy.array([a.xyz for a in atoms2])
    result = []
    for i in range(0, len(xyz1), 512):
        block = xyz1[i:i + 512]
        dist = numpy.sqrt(((block[:, None, :] - xyz2[None, :, :])**2).sum(-1))
        for a, b in zip(*numpy.nonzero(dist <= cutoff)):
            result.append((atoms1[i + a], atoms2[b], float(dist[a, b])))
    return result


class Interface:
    '''
    Result for one chain pair (see module docstring).
    '''

    def __init__(self, obj, chain1, chain2):
        self.obj = obj
        self.chains = (chain1, chain2)
        self.area = 0.0
        self.dg = 0.0
        self.hbonds = []        # (donor/acceptor atom, atom, distance)
        self.salt_bridges = []  # (atom, atom, distance)
        self.disulfides = []    # (atom, atom, distance)
        self.residues = []      # dicts, see analyze_interface
        self.bsa_by_type = defaultdict(float)  # ASP type -> buried area

    def counts(self):
        return (len(self.hbonds), len(self.salt_bridges),
                len(self.disulfides))

    def residue_count(self, side):
        chain = self.chains[side]
        return sum(1 for r in self.residues if r['chain'] == chain)


def analyze_interface(obj, chain1, chain2, atoms1, atoms2, state, _self):
    '''
    :param atoms1, atoms2: Atom lists with isolated-chain areas
    '''
    iface = Interface(obj, chain1, chain2)
    complex_areas = _isolated_areas(
        '(%s) or (%s)' % (chain_selection(obj, chain1),
                          chain_selection(obj, chain2)), state, _self)

    bsa_total = 0.0
    dg = 0.0
    residues = defaultdict(lambda: {'asa': 0.0, 'bsa': 0.0, 'dg': 0.0})
    for atom in itertools.chain(atoms1, atoms2):
        res = residues[atom.residue]
        res['asa'] += atom.asa
        atom.bsa = max(0.0, atom.asa - complex_areas.get(atom.key, atom.asa))
        if atom.bsa <= 0.0:
            continue
        bsa_total += atom.bsa
        iface.bsa_by_type[atom.asp_type()] += atom.bsa
        sigma = ASP.get(atom.asp_type(), 0.0)
        contribution = -sigma * atom.bsa / 1000.0  # kcal/mol
        dg += contribution
        res['bsa'] += atom.bsa
        res['dg'] += contribution

    iface.area = bsa_total / 2.0
    iface.dg = dg

    # bonds between atoms near the interface
    near1 = [a for a in atoms1 if a.elem in ('N', 'O', 'S')]
    near2 = [a for a in atoms2 if a.elem in ('N', 'O', 'S')]
    for a, b, d in _pairs_within(near1, near2, HBOND_CUTOFF):
        if ((a.is_donor() and b.is_acceptor()) or
                (a.is_acceptor() and b.is_donor())):
            iface.hbonds.append((a, b, d))
    for a, b, d in _pairs_within(near1, near2, SALT_BRIDGE_CUTOFF):
        if (a.is_anion() and b.is_cation()) or (a.is_cation() and b.is_anion()):
            iface.salt_bridges.append((a, b, d))
    sg1 = [a for a in near1 if a.resn == 'CYS' and a.name == 'SG']
    sg2 = [a for a in near2 if a.resn == 'CYS' and a.name == 'SG']
    iface.disulfides = _pairs_within(sg1, sg2, DISULFIDE_CUTOFF)

    flags = defaultdict(set)
    for kind, bonds in (('H', iface.hbonds), ('S', iface.salt_bridges),
                        ('D', iface.disulfides)):
        for a, b, _ in bonds:
            flags[a.residue].add(kind)
            flags[b.residue].add(kind)

    for (chain, resi, resn), res in residues.items():
        if res['bsa'] > 0.0:
            iface.residues.append({
                'chain': chain, 'resi': resi, 'resn': resn,
                'asa': res['asa'], 'bsa': res['bsa'], 'dg': res['dg'],
                'flags': ''.join(k for k in 'HSD' if k in flags[(chain, resi, resn)]),
            })
    return iface


def analyze(obj, state=-1, *, _self=cmd):
    '''
    PISA-like analysis of all chain-chain interfaces of an object.

    :return: list of Interface, largest interface first
    '''
    if state < 1:
        state = _self.get_state()

    chains = [c for c in _self.get_chains('%%%s and polymer' % obj)]
    if len(chains) < 2:
        return []

    settings = {name: _self.get(name) for name in ('dot_solvent',
                                                   'dot_density')}
    view = _self.get_view()
    with _self.UndoPauseCM() if hasattr(_self, 'UndoPauseCM') else _nullcm():
        try:
            _self.set('dot_solvent', 1, quiet=1)
            _self.set('dot_density', SASA_DOT_DENSITY, quiet=1)

            atoms = {}
            for chain in chains:
                sele = chain_selection(obj, chain)
                atoms[chain] = _load_atoms(sele, state, _self)
                areas = _isolated_areas(sele, state, _self)
                for atom in atoms[chain]:
                    atom.asa = areas.get(atom.key, 0.0)

            interfaces = []
            for chain1, chain2 in itertools.combinations(chains, 2):
                contact = _self.count_atoms(
                    '(%s) within %s of (%s)' % (
                        chain_selection(obj, chain1), CONTACT_CUTOFF,
                        chain_selection(obj, chain2)), state=state)
                if not contact:
                    continue
                iface = analyze_interface(obj, chain1, chain2, atoms[chain1],
                                          atoms[chain2], state, _self)
                if iface.area > 0.0:
                    interfaces.append(iface)
        finally:
            for name, value in settings.items():
                _self.set(name, value, quiet=1)
            _self.set_view(view)

    interfaces.sort(key=lambda i: -i.area)
    return interfaces


class _nullcm:
    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def format_table(obj, interfaces):
    lines = [
        ' Interfaces of "%s" (PISA-like, computed by pymolx)' % obj,
        '  %2s  %-9s %9s %13s %4s %4s %4s  %s' % (
            '#', 'Chains', 'Area(A^2)', 'dG(kcal/mol)', 'HB', 'SB', 'DS',
            'Residues'),
    ]
    for n, iface in enumerate(interfaces, 1):
        hb, sb, ds = iface.counts()
        lines.append('  %2d  %-9s %9.1f %13.1f %4d %4d %4d  %d+%d' % (
            n, '%s-%s' % iface.chains, iface.area, iface.dg, hb, sb, ds,
            iface.residue_count(0), iface.residue_count(1)))
    if not interfaces:
        lines.append('  (no interfaces between polymer chains)')
    lines.append(' dG: solvation energy gain, atomic solvation parameters '
                 'fitted to PISA (estimate)')
    return '\n'.join(lines)


def _resi(resi):
    # negative residue numbers need escaping in selections
    return '\\' + resi if resi.startswith('-') else resi


def show_interface(iface, *, prefix=None, _self=cmd):
    '''
    Visualize an interface: residues as sticks, bonds as dashes, all in a
    group named <prefix> (default iface_<obj>_<chain1><chain2>).

    :return: group name
    '''
    obj = iface.obj
    c1, c2 = iface.chains
    group = prefix or _self.get_legal_name('iface_%s_%s%s' % (obj, c1, c2))

    residues = ' or '.join(
        '(chain "%s" and resi %s)' % (r['chain'], _resi(r['resi']))
        for r in iface.residues) or 'none'
    _self.select(group + '_res', '%%%s and (%s)' % (obj, residues), enable=0)
    _self.show('sticks', '%s_res and not name N+C+O' % group)

    for suffix, bonds, color in (('hbonds', iface.hbonds, 'yellow'),
                                 ('saltbr', iface.salt_bridges, 'magenta'),
                                 ('disulf', iface.disulfides, 'orange')):
        name = '%s_%s' % (group, suffix)
        _self.delete(name)
        for a, b, _ in bonds:
            _self.distance(name, '%%%s and index %d' % (obj, a.index),
                           '%%%s and index %d' % (obj, b.index), label=0)
        if bonds:
            _self.color(color, name)

    members = ' '.join(n for n in _self.get_names('all')
                       if n.startswith(group + '_'))
    _self.group(group, members)
    _self.zoom(group + '_res', 4, animate=1)
    return group


def interface_analysis(selection='', state=-1, quiet=0, *, _self=cmd):
    '''
DESCRIPTION

    PISA-like analysis of the interfaces between the polymer chains of
    an object: interface area, solvation energy gain estimate (dG),
    hydrogen bonds, salt bridges and disulfides.

    dG uses atomic solvation parameters fitted to PISA's dG_int, so it
    is typically within 1 kcal/mol of PISA. For P-values, CSS and
    assembly prediction, use PISA itself.

USAGE

    interface_analysis object [, state ]

ARGUMENTS

    object = str: object name {default: first object}

    state = int: state {default: -1 (current)}

SEE ALSO

    distance, get_area
    '''
    obj = selection.strip() or (_self.get_object_list() or [''])[0]
    interfaces = analyze(obj, int(state), _self=_self)
    if not int(quiet):
        print(format_table(obj, interfaces))
    return interfaces


# epitope: the target's residues buried by the other chains


def default_targets(obj, _self=cmd):
    '''
    Target chains: the chains that aren't antibody chains (found by their
    framework motifs, pymolx.cluster), else none

    :return: sorted chain IDs
    '''
    from .cluster import antibody_chains
    chains = sorted(_self.get_chains('%%%s and polymer' % obj))
    antibody = antibody_chains(obj, _self)
    if not antibody or len(antibody) == len(chains):
        return []
    return [c for c in chains if c not in antibody]


def _resi_sort_key(resi):
    digits = ''.join(ch for ch in resi if ch.isdigit() or ch == '-')
    try:
        return (int(digits), resi)
    except ValueError:
        return (0, resi)


def contact_residues(interface_list, targets, min_bsa=0.0, side='target'):
    '''
    Residues buried between the target chains and the others: those of
    the targets (side "target": the epitope) or of the others ("partner":
    the paratope). Interfaces between two target chains, or between two
    others, don't count.

    :param min_bsa: buried area a residue needs, A^2 (more than this)
    :return: {chain: [(resi, resn, buried area)]}, residues in order
    '''
    targets = set(targets)
    found = {}
    for iface in interface_list:
        a, b = iface.chains
        if (a in targets) == (b in targets):
            continue
        for r in iface.residues:
            if (r['chain'] in targets) != (side == 'target'):
                continue
            key = (r['chain'], r['resi'])
            entry = found.setdefault(key, [r['resn'], 0.0])
            entry[1] += r['bsa']
    result = {}
    for (chain, resi), (resn, bsa) in found.items():
        if bsa > min_bsa:
            result.setdefault(chain, []).append((resi, resn, bsa))
    for residues in result.values():
        residues.sort(key=lambda r: _resi_sort_key(r[0]))
    return dict(sorted(result.items()))


NAME_STYLES = ('one', 'three', 'none')


def residue_label(resi, resn, names='one'):
    '''
    "D52" (one), "ASP52" (three) or "52" (none)
    '''
    if names == 'one':
        from .cluster import THREE_TO_ONE
        return THREE_TO_ONE.get(resn, 'X') + resi
    if names == 'three':
        return resn + resi
    return resi


def format_residues(residues, names='one'):
    '''
    "D52, Y53, W100A" (names "one"), "ASP52, ..." ("three") or "52, 53,
    100A" ("none") from contact_residues' list for one chain
    '''
    return ', '.join(residue_label(resi, resn, names)
                     for resi, resn, _ in residues)


def format_epitope(by_chain, names='one'):
    '''
    One chain: its residues; several: one line per chain ("A: D52, Y53")
    '''
    if len(by_chain) == 1:
        return format_residues(next(iter(by_chain.values())), names)
    return '\n'.join('%s: %s' % (chain or '-',
                                 format_residues(residues, names))
                     for chain, residues in by_chain.items())


# proximity profile: how close each residue of the targets comes to the
# other chains

# atoms closer than this (clashes of a predicted model) score no higher
PROXIMITY_MIN_DISTANCE = 1.0
CONTACT_DISTANCE = 4.0     # the default "in contact" line of the profile
# the score at which the profile's colors end (2.5 A): fixed, so that
# two profiles, or two colored structures, can be compared
PROXIMITY_SCORE_TOP = 0.4
# the viewer's colors for a profile: this many steps along the ramp, named
PROXIMITY_COLOR_STEPS = 100
PROXIMITY_COLOR_NAME = 'proximity_%03d'


class ChainAtoms:
    '''
    Heavy atoms of an object's polymer chains with their coordinates,
    read once (residue_atoms) so that proximity() needs no PyMOL calls
    '''

    def __init__(self):
        self.residues = []  # (chain, resi, resn), in the object's order
        self.owner = []     # per atom: index of its residue
        self.xyz = []       # per atom: coordinates

    def chains(self):
        return sorted({chain for chain, _, _ in self.residues})


def residue_atoms(obj, state=-1, *, _self=cmd):
    '''
    :return: ChainAtoms of the object's polymer chains
    '''
    if state < 1:
        state = _self.get_state()
    rows = []
    _self.iterate_state(
        state, '%%%s and polymer and not hydro' % obj,
        'rows.append((chain, resi, resn, x, y, z))', space={'rows': rows})
    atoms = ChainAtoms()
    index = {}
    for chain, resi, resn, x, y, z in rows:
        key = (chain, resi, resn)
        if key not in index:
            index[key] = len(atoms.residues)
            atoms.residues.append(key)
        atoms.owner.append(index[key])
        atoms.xyz.append((x, y, z))
    return atoms


def _nearest(xyz1, xyz2, tree=True):
    '''
    For each point of xyz1, the closest point of xyz2

    :param tree: use scipy's k-d tree if it is installed
    :return: (distances, indices into xyz2)
    '''
    cKDTree = None
    if tree:
        try:
            from scipy.spatial import cKDTree
        except ImportError:
            pass
    if cKDTree is not None:
        return cKDTree(xyz2).query(xyz1)
    distance = numpy.empty(len(xyz1))
    index = numpy.empty(len(xyz1), dtype=int)
    square2 = (xyz2**2).sum(1)
    for i in range(0, len(xyz1), 256):
        block = xyz1[i:i + 256]
        squared = ((block**2).sum(1)[:, None] + square2[None, :] -
                   2.0 * block @ xyz2.T)
        best = squared.argmin(1)
        index[i:i + 256] = best
        distance[i:i + 256] = numpy.sqrt(numpy.maximum(
            squared[numpy.arange(len(block)), best], 0.0))
    return distance, index


def proximity(atoms, targets):
    '''
    How close each residue of the target chains comes to the other
    chains: the distance between the closest pair of heavy atoms, one
    of the residue and one of any chain that isn't a target.

    :param atoms: ChainAtoms (residue_atoms)
    :return: {chain: [(resi, resn, distance, partner)]}: every residue of
        the target chains in the object's order, partner being the
        (chain, resi, resn) of the residue it is closest to; {} without
        targets or without other chains
    '''
    targets = set(targets)
    is_target = numpy.array([atoms.residues[k][0] in targets
                             for k in atoms.owner], dtype=bool)
    if not is_target.any() or is_target.all():
        return {}
    owner = numpy.array(atoms.owner)
    xyz = numpy.array(atoms.xyz, dtype=float)
    distance, index = _nearest(xyz[is_target], xyz[~is_target])
    partner_owner = owner[~is_target]
    closest = {}  # residue index -> (distance, partner residue index)
    for k, d, j in zip(owner[is_target].tolist(), distance.tolist(),
                       index.tolist()):
        if k not in closest or d < closest[k][0]:
            closest[k] = (d, int(partner_owner[j]))
    result = {}
    for k in sorted(closest):
        chain, resi, resn = atoms.residues[k]
        d, partner = closest[k]
        result.setdefault(chain, []).append(
            (resi, resn, d, atoms.residues[partner]))
    return dict(sorted(result.items()))


def proximity_score(distance):
    '''
    The profile's score, inversely proportional to the distance: 1/A
    (0.25 at 4 A, 0.1 at 10 A)
    '''
    return 1.0 / max(distance, PROXIMITY_MIN_DISTANCE)


def proximity_fraction(distance):
    '''
    Where a distance lies on the profile's color ramp: 0 (far) to 1
    (PROXIMITY_SCORE_TOP and above: 2.5 A and closer)
    '''
    return min(1.0, proximity_score(distance) / PROXIMITY_SCORE_TOP)


def ramp_rgb(ramp, fraction):
    '''
    The color at a fraction (0 to 1) of a ramp

    :param ramp: its colors from 0 to 1, each "#rrggbb"
    :return: (r, g, b), each 0 to 1
    '''
    stops = [tuple(int(color.lstrip('#')[k:k + 2], 16) / 255.0
                   for k in (0, 2, 4)) for color in ramp]
    position = min(1.0, max(0.0, fraction)) * (len(stops) - 1)
    k = min(int(position), len(stops) - 2)
    return tuple(a + (b - a) * (position - k)
                 for a, b in zip(stops[k], stops[k + 1]))


def color_proximity(obj, by_chain, ramp, *, _self=cmd):
    '''
    Color the residues of a profile in the viewer by how close they come
    to the other chains, as the profile's graph does: all their atoms,
    and nothing else of the object.

    :param by_chain: proximity()'s {chain: [(resi, resn, distance,
        partner)]}
    :param ramp: the colors from far to close, each "#rrggbb" (the
        graph's are theme.COLORS' proximity_0 to proximity_4)
    :return: the number of residues colored
    '''
    indices = []
    for step in range(PROXIMITY_COLOR_STEPS + 1):
        name = PROXIMITY_COLOR_NAME % step
        _self.set_color(name, list(ramp_rgb(
            ramp, step / PROXIMITY_COLOR_STEPS)))
        indices.append(_self.get_color_index(name))
    colors = {}
    for chain, rows in by_chain.items():
        for resi, _, distance, _ in rows:
            colors[(chain, resi)] = indices[round(
                PROXIMITY_COLOR_STEPS * proximity_fraction(distance))]
    if colors:
        # (one pass over the object, not a "color" per residue)
        _self.alter('%%%s and polymer' % obj,
                    'color = colors.get((chain, resi), color)',
                    space={'colors': colors})
        _self.recolor('%' + obj)
    return len(colors)


def format_proximity(by_chain):
    '''
    The profile as a tab-separated table, one line per residue
    '''
    lines = ['\t'.join(('chain', 'residue', 'number', 'distance (A)',
                        'score (1/A)', 'closest chain', 'closest residue',
                        'closest number'))]
    for chain, residues in by_chain.items():
        for resi, resn, distance, partner in residues:
            lines.append('\t'.join((
                chain, resn, resi, '%.2f' % distance,
                '%.3f' % proximity_score(distance), partner[0], partner[2],
                partner[1])))
    return '\n'.join(lines)


def epitope(selection='', target='', min_bsa=0.0, names='one', state=-1,
            quiet=0, *, _self=cmd):
    '''
DESCRIPTION

    "epitope" lists the residues of the target chain(s) buried by the
    other chains (e.g. an antigen's epitope on an antibody), from the
    PISA-like interface analysis, separated by commas: "D52, Y53".

USAGE

    epitope object [, target [, min_bsa [, names [, state ]]]]

ARGUMENTS

    object = str: object name {default: first object}

    target = str: target chain(s), e.g. "A" or "A+B" {default: the chains
    that aren't antibody chains}

    min_bsa = float: buried area a residue needs, A^2 {default: 0}

    names = one (D52), three (ASP52) or none (52) {default: one}

EXAMPLE

    epitope model_1, A, names=three

SEE ALSO

    interface_analysis
    '''
    names = str(names).strip().lower()
    if names not in NAME_STYLES:
        raise ValueError('names must be one, three or none')
    obj = selection.strip() or (_self.get_object_list() or [''])[0]
    targets = [c for c in str(target).replace('+', ' ').split() if c] or \
        default_targets(obj, _self)
    if not targets:
        raise ValueError('give the target chain(s): no antibody chains to '
                         'tell the target from')
    chains = set(_self.get_chains('%%%s and polymer' % obj))
    missing = [c for c in targets if c not in chains]
    if missing:
        raise ValueError('no chain %s in %s' % (', '.join(missing), obj))
    interface_list = analyze(obj, int(state), _self=_self)
    by_chain = contact_residues(interface_list, targets, float(min_bsa))
    if not int(quiet):
        if not by_chain:
            print(' epitope: no residues of %s buried by other chains' %
                  '+'.join(targets))
        for chain, residues in by_chain.items():
            print(' epitope, chain %s (%d residues, %.0f A^2 buried):' % (
                chain, len(residues), sum(r[2] for r in residues)))
            print(' ' + format_residues(residues, names))
    return format_epitope(by_chain, names)
