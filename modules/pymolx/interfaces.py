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
