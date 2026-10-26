'''
Binding pocket (command "binding_pocket", Binding Pocket dialog): the
residues with an atom within a distance of a ligand (small molecule,
ion, or any selection), as a selection and as a list ("D52, Y53" or
"ASP52, TYR53" or "52, 53"; one line per chain when there are several).

Heavy atoms only, unless asked (results then don't depend on whether
hydrogens were added). Waters are left out unless asked. Each residue
gets its closest distance and atom pair, and "polar" when one of its N/O
atoms is within 3.5 A of a ligand N/O (a likely hydrogen bond or salt
bridge).
'''

import numpy

from pymol import cmd

from .interfaces import NAME_STYLES, format_epitope, _resi, _resi_sort_key

POLAR_CUTOFF = 3.5


def ligands(obj, _self=cmd):
    '''
    Small molecules and ions of an object (not polymer, not water)

    :return: [{'chain', 'resi', 'resn', 'atoms', 'selection'}]
    '''
    found = {}
    _self.iterate('%%%s and not polymer and not solvent' % obj,
                  '_ = found.setdefault((chain, resi, resn), []).append(1)',
                  space={'found': found})
    result = []
    for (chain, resi, resn), atoms in found.items():
        result.append({
            'chain': chain, 'resi': resi, 'resn': resn,
            'atoms': len(atoms),
            'selection': '%%%s and chain "%s" and resi %s and resn %s' % (
                obj, chain, _resi(resi), resn)})
    result.sort(key=lambda r: (r['atoms'] == 1, r['chain'],
                               _resi_sort_key(r['resi'])))
    return result


def _atoms(selection, state, _self):
    atoms = []
    _self.iterate_state(
        state, selection,
        'atoms.append((model, chain, resi, resn, name, elem, (x, y, z)))',
        space={'atoms': atoms})
    return atoms


class Pocket:
    '''
    residues: [{'model', 'chain', 'resi', 'resn', 'distance', 'atom',
    'ligand_atom', 'polar'}] in chain/number order
    '''

    def __init__(self, ligand, cutoff):
        self.ligand = ligand
        self.cutoff = cutoff
        self.residues = []
        self.selection = ''

    def by_chain(self):
        '''
        {chain: [(resi, resn, distance)]} for interfaces.format_epitope
        '''
        chains = {}
        for r in self.residues:
            chains.setdefault(r['chain'], []).append(
                (r['resi'], r['resn'], r['distance']))
        return chains

    def text(self, names='one'):
        return format_epitope(self.by_chain(), names) if self.residues \
            else ''

    def save_csv(self, path):
        import csv
        with open(path, 'w', newline='') as handle:
            out = csv.writer(handle)
            out.writerow(['object', 'chain', 'resi', 'resn', 'distance',
                          'atom', 'ligand_atom', 'polar'])
            for r in self.residues:
                out.writerow([r['model'], r['chain'], r['resi'], r['resn'],
                              '%.2f' % r['distance'], r['atom'],
                              r['ligand_atom'], int(r['polar'])])


def find(ligand, cutoff=4.0, protein='', hydrogens=False, waters=False,
         name='pocket', state=-1, _self=cmd):
    '''
    The pocket around a ligand, and a selection of its residues.

    :param ligand: selection of the ligand atoms
    :param protein: selection of the residues to consider {default: the
        polymer of the ligand's object(s)}
    :return: Pocket
    '''
    ligand = str(ligand).strip()
    cutoff = float(cutoff)
    if cutoff <= 0:
        raise ValueError('the distance must be positive')
    state = int(state)
    if state < 1:
        state = _self.get_state()
    heavy = '' if hydrogens else ' and not hydro'
    ligand_sel = '(%s)%s' % (ligand, heavy)
    if not _self.count_atoms(ligand_sel):
        raise ValueError('no atoms in the ligand selection "%s"' % ligand)
    objects = _self.get_object_list(ligand_sel) or []
    protein = str(protein).strip() or '(%s) and polymer' % ' or '.join(
        '%' + o for o in objects)
    if waters:
        protein = '(%s) or ((%s) and solvent)' % (protein, ' or '.join(
            '%' + o for o in objects))
    candidates = '((%s)%s) and not (%s)' % (protein, heavy, ligand)

    near = _atoms('(%s) within %g of (%s)' % (candidates, cutoff, ligand_sel),
                  state, _self)
    keys = {(a[0], a[1], a[2]) for a in near}
    lig = _atoms(ligand_sel, state, _self)
    lig_xyz = numpy.array([b[6] for b in lig], float)
    lig_polar = numpy.array([b[5] in ('N', 'O') for b in lig])
    result = Pocket(ligand, cutoff)
    if keys:
        residue_atoms = _atoms('(byres ((%s) within %g of (%s)))%s' % (
            candidates, cutoff, ligand_sel, heavy), state, _self)
        by_residue = {}
        for a in residue_atoms:
            if (a[0], a[1], a[2]) in keys:
                by_residue.setdefault((a[0], a[1], a[2]), []).append(a)
        for (model, chain, resi), atoms in by_residue.items():
            xyz = numpy.array([a[6] for a in atoms], float)
            dist = numpy.linalg.norm(xyz[:, None] - lig_xyz[None], axis=2)
            i, j = numpy.unravel_index(numpy.argmin(dist), dist.shape)
            d, a, b = float(dist[i, j]), atoms[i], lig[j]
            polar = numpy.array([x[5] in ('N', 'O') for x in atoms])
            polar = bool((dist[polar][:, lig_polar] <= POLAR_CUTOFF).any())
            result.residues.append({
                'model': model, 'chain': chain, 'resi': resi,
                'resn': a[3], 'distance': d, 'atom': a[4],
                'ligand_atom': '%s%s %s' % (b[3], b[2], b[4]),
                'polar': polar})
    result.residues.sort(key=lambda r: (r['model'], r['chain'],
                                        _resi_sort_key(r['resi'])))
    if name:
        name = _self.get_legal_name(str(name).strip())
        parts = ['(%%%s and chain "%s" and resi %s)' % (
            r['model'], r['chain'], _resi(r['resi'])) for r in result.residues]
        _self.select(name, ' or '.join(parts) or 'none', enable=0)
        result.selection = name
    return result


def show(result, _self=cmd):
    '''
    The ligand and the pocket's residues as sticks, zoomed
    '''
    _self.show('sticks', '(%s) and not hydro' % result.ligand)
    if result.selection:
        _self.show('sticks', '(%s) and not name N+C+O and not hydro' %
                   result.selection)
        _self.zoom('(%s) or (%s)' % (result.ligand, result.selection), 3,
                   animate=0.5)


def binding_pocket(ligand='organic', cutoff=4.0, names='one', name='pocket',
                   protein='', hydrogens=0, waters=0, file='', state=-1,
                   quiet=0, *, _self=cmd):
    '''
DESCRIPTION

    "binding_pocket" lists the residues with an atom within "cutoff" of
    the ligand (heavy atoms) and selects them.

USAGE

    binding_pocket [ ligand [, cutoff [, names [, name ]]]]

ARGUMENTS

    ligand = str: the ligand's atoms {default: organic}

    cutoff = float: maximum distance, A {default: 4.0}

    names = one (D52), three (ASP52) or none (52) {default: one}

    name = str: the selection made {default: pocket}

    protein = str: the residues to consider {default: the polymer of the
    ligand's object}

    hydrogens = 0/1: hydrogens count too {default: 0}

    waters = 0/1: waters can be pocket residues {default: 0}

    file = str: also save the residues as CSV

EXAMPLE

    binding_pocket resn ATP, 5
    binding_pocket 1abc and resn HEM, names=three

SEE ALSO

    epitope, interface_analysis
    '''
    names = str(names).strip().lower()
    if names not in NAME_STYLES:
        raise ValueError('names must be one, three or none')
    result = find(ligand, cutoff, protein, int(hydrogens), int(waters),
                  name, state, _self=_self)
    if file:
        result.save_csv(_self.exp_path(file))
    if not int(quiet):
        polar = sum(r['polar'] for r in result.residues)
        print(' binding_pocket: %d residues within %g A of %s (%d with '
              'polar contacts); selection "%s"' % (
                  len(result.residues), result.cutoff, result.ligand, polar,
                  result.selection))
        if result.residues:
            print(' ' + result.text(names).replace('\n', '\n '))
    return result.text(names)


def extend(_self):
    _self.extend('binding_pocket', binding_pocket)
    _self.binding_pocket = binding_pocket
