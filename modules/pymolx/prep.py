'''
Protein preparation workflow, like Maestro's Protein Preparation Wizard
(Structure Preparation dialog, command "prepare_structure"):

1. structure: a selection (one state), read from PyMOL
2. preprocess (PDBFixer): delete waters / other heterogens, replace
   nonstandard residues, fill in missing residues (from the structure's
   own record of unmodeled residues, or a sequence file / PDB ID), fill
   in missing side-chain atoms, cap chain termini with ACE / NME
3. protonation at a pH: pKa values from PROPKA (or the standard ones);
   His tautomers from hydrogen bonds (OpenMM Modeller)
4. restrained minimization (OpenMM): force field, implicit / explicit /
   no solvent, restraints, convergence tolerance, iterations
5. diagnostics of the result (pymolx.diagnostics), by the caller

PrepJob.run does steps 2-4 without PyMOL calls (it may run in a worker
thread); prepare reads PyMOL and PrepJob.apply loads the result as a new
object. Waters and ligands are exported as chains of their own (so
sequences, caps and terminal hydrogens see only the polymer) and get
their chains back on loading.
'''

import io
import logging
import math
import threading
import time

from pymol import cmd

from . import mm

AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
}
PKA_METHODS = ('propka', 'standard')
SOLVENTS = mm.IMPLICIT_SOLVENTS + ('explicit', 'vacuum')
RESTRAINTS = ('heavy', 'backbone', 'ca', 'none')
BACKBONE = {'N', 'CA', 'C', 'O'}
EXPLICIT_PADDING = 1.0  # nm of water around the structure
IONIC_STRENGTH = 0.15   # molar (Na+ Cl-), plus neutralizing ions
BUILD_K = 1000.0        # kJ/mol/rad^2: built residues' chirality, trans
MAX_ROUNDS = 5          # minimizer restarts while some force is large
MAX_FORCE_FACTOR = 50   # ... larger than this times the tolerance
L_IMPROPER = 2.138      # rad, dihedral N-C-CA-CB of L-amino acids (122.5)
L_IMPROPER_HA = -2.062  # rad, dihedral N-C-CA-HA (-118.1)

# protonation variants (OpenMM Modeller names): (protonated, deprotonated)
TITRATABLE = {'ASP': ('ASH', 'ASP'), 'GLU': ('GLH', 'GLU'),
              'HIS': ('HIP', None), 'LYS': ('LYS', 'LYN')}
CHARGE = {'ASP': -1, 'GLU': -1, 'LYS': 1, 'ARG': 1, 'HIP': 1, 'CYM': -1,
          'ASH': 0, 'GLH': 0, 'LYN': 0, 'HID': 0, 'HIE': 0,
          'NA': 1, 'K': 1, 'CL': -1, 'MG': 2, 'CA': 2, 'ZN': 2, 'MN': 2,
          'FE2': 2, 'FE': 3, 'CO': 2, 'NI': 2, 'CU': 2, 'CD': 2}
# metal ions (element): coordination distance cutoff, A
METALS = {'ZN': 2.6, 'MG': 2.6, 'CA': 2.8, 'MN': 2.6, 'FE': 2.6, 'CO': 2.6,
          'NI': 2.6, 'CU': 2.6, 'CD': 2.8, 'NA': 2.9, 'K': 3.2}
# protonation states of metal-coordinating residues (by coordinating atom)
METAL_STATES = {('HIS', 'ND1'): 'HIE', ('HIS', 'NE2'): 'HID',
                ('CYS', 'SG'): 'CYX', ('ASP', 'OD1'): 'ASP',
                ('ASP', 'OD2'): 'ASP', ('GLU', 'OE1'): 'GLU',
                ('GLU', 'OE2'): 'GLU'}
ONE_LETTER = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
              'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
              'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
              'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}


class Cancelled(Exception):
    pass


def propka_available():
    try:
        import propka.run  # noqa: F401
        return True
    except ImportError:
        return False


def _chain_label(chain):
    return chain if chain else '(blank)'


# reading PyMOL


def _choose_alternates(rows, altloc):
    '''
    One alternate location per residue: the one with the highest mean
    occupancy (ties: the first letter), or the letter "altloc" where the
    residue has it.

    :return: ({(chain, resi): letter}, [(chain, resi, resn, letter,
        {letter: occupancy})])
    '''
    occupancy = {}
    resn = {}
    for (record, elem, name, alt, rname, chain, resv, resi, x, y, z, q,
         b) in rows:
        if alt:
            occupancy.setdefault((chain, resi), {}).setdefault(
                alt, []).append(q)
            resn[(chain, resi)] = rname
    chosen = {}
    report = []
    for key, letters in occupancy.items():
        means = {k: sum(v) / len(v) for k, v in letters.items()}
        if altloc != 'occupancy' and altloc in means:
            letter = altloc
        else:
            letter = sorted(means, key=lambda k: (-means[k], k))[0]
        chosen[key] = letter
        report.append((key[0], key[1], resn[key], letter, means))
    return chosen, report


def _export(selection, state, _self, hydrogens=False, altloc='occupancy'):
    '''
    mmCIF text of the selection (with or without its hydrogens), one
    alternate location per residue (_choose_alternates). Each chain's
    polymer residues and its other residues (waters, ligands, ions)
    become chains of their own, with generated IDs.

    :return: (text, {chain ID: PyMOL chain}, {PyMOL chain: polymer chain
        ID}), b-factors {(chain, resi, name): b}, alternates report
    '''
    atoms = '(%s)%s' % (selection, '' if hydrogens else ' and not hydro')
    polymer = set()
    _self.iterate_state(state, '%s and polymer' % atoms,
                        'polymer.add((chain, resi))',
                        space={'polymer': polymer})
    rows = []
    _self.iterate_state(
        state, atoms,
        'rows.append((type, elem, name, alt, resn, chain, resv, resi, '
        'x, y, z, q, b))', space={'rows': rows})
    if not rows:
        raise ValueError('no atoms with coordinates in "%s"' % selection)
    chosen, alternates = _choose_alternates(rows, altloc)
    rows = [row[:3] + ('',) + row[4:] for row in rows
            if not row[3] or row[3] == chosen[(row[5], row[7])]]

    # chains in order of appearance; polymer first, then the rest
    chain_order = []
    for row in rows:
        if row[5] not in chain_order:
            chain_order.append(row[5])
    rows.sort(key=lambda row: (chain_order.index(row[5]),
                               (row[5], row[7]) not in polymer))

    ids = {}
    chain_map = {}
    polymer_ids = {}
    seen = set()
    duplicates = []
    lines = ['data_pymolx', 'loop_'] + ['_atom_site.' + column for column in (
        'group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id',
        'label_comp_id', 'label_asym_id', 'label_seq_id',
        'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy',
        'B_iso_or_equiv', 'auth_seq_id', 'auth_asym_id',
        'pdbx_PDB_model_num')]
    b_factors = {}
    for n, (record, elem, name, alt, resn, chain, resv, resi, x, y, z, q,
            b) in enumerate(rows, 1):
        is_polymer = (chain, resi) in polymer
        group = (chain, is_polymer)
        if group not in ids:
            ids[group] = 'c%d' % (len(ids) + 1)
            chain_map[ids[group]] = chain
            if is_polymer:
                polymer_ids[chain] = ids[group]
        key = (chain, resi, name, alt)
        if key in seen:
            duplicates.append(key)
        seen.add(key)
        b_factors[(chain, resi, name)] = b
        inscode = resi[len(str(resv)):] if resi.startswith(str(resv)) else ''
        lines.append(' '.join(mm._cif_value(v) for v in (
            record, n, elem, name, alt, resn, ids[group],
            resv, inscode or '?', '%.3f' % x, '%.3f' % y, '%.3f' % z,
            '%.2f' % q, '%.2f' % b, resv, ids[group], 1)))
    if duplicates:
        raise ValueError(
            'atoms are not unique by chain, residue and name (e.g. %s); '
            'rename chains or residues first' % (
                ', '.join('%s/%s/%s' % k[:3] for k in duplicates[:3])))
    return ('\n'.join(lines) + '\n', chain_map, polymer_ids, b_factors,
            alternates)


def _unmodeled(selection, _self):
    '''
    Each chain's polymer residues in order, with those PyMOL knows of but
    has no coordinates for (PDB REMARK 465, mmCIF unobserved residues)
    marked as missing.

    :return: {chain: [(resi, resn, present)]}
    '''
    order = {}
    # (an assignment: iterate prints the value of an expression)
    _self.iterate('(%s) and polymer' % selection,
                  '_ = order.setdefault(chain, {}).setdefault(resi, resn)',
                  space={'order': order})
    present = set()
    _self.iterate('(%s) and polymer and present' % selection,
                  'present.add((chain, resi))', space={'present': present})
    return {chain: [(resi, resn, (chain, resi) in present)
                    for resi, resn in residues.items()]
            for chain, residues in order.items()}


# the job


class PrepJob:
    '''
    One preparation: prepare_job (reads PyMOL), run() (no PyMOL calls),
    apply() (loads the result).
    '''

    def __init__(self):
        self.selection = 'all'
        self.state = 1
        self.objects = []
        self.name = ''
        # preprocess
        self.preprocess = True
        self.altloc = 'occupancy'  # or a letter
        self.delete_hydrogens = True
        self.cap_termini = False
        self.fill_loops = True
        self.fill_termini = False
        self.source = ''
        self.fill_side_chains = True
        self.replace_nonstandard = True
        self.delete_waters = False
        self.water_distance = 0.0   # A: delete only waters farther from
                                    # het groups (0: all waters)
        self.water_min_hbonds = 0   # delete waters with fewer H-bonds to
                                    # non-waters (0: off)
        self.delete_heterogens = False
        self.mutations = []          # [(chain, resi, new name)], resolved
                                     # by prepare_job: (chain, resi, old,
                                     # new)
        self.disulfides = True       # keep (create) disulfide bonds
        self.metals = True           # metal sites: states, bonds, restraints
        # protonation
        self.protonate = True
        self.ph = 7.0
        self.pka_method = 'propka'
        self.ligand_hydrogens = True   # from the PDB's CCD (downloaded)
        self.optimize_hbonds = True    # flips, tautomers, rotatable H
        # minimization
        self.minimize = True
        self.forcefield = 'amber14'
        self.solvent = 'gbn2'
        self.restrain = 'heavy'
        self.restraint_k = 1000.0   # kJ/mol/nm^2
        self.tolerance = 10.0       # kJ/mol/nm
        self.max_iterations = 0     # 0: until converged
        self.platform = 'auto'
        self.diagnostics = True     # (for the caller)

        self.input_text = ''
        self.chain_map = {}
        self.polymer_ids = {}
        self.residues = {}
        self.b_factors = {}
        self.input_hydrogens = 0
        self.alternates = []
        self.het_bonds = []
        self.metal_bonds = []
        self._metal_pairs = []
        self._metal_states = {}
        self.topology = None
        self.positions = None
        self.report = {}
        self.log = []
        self.status = ''
        self.progress = None        # callable(iteration, energy kJ/mol)
        self._cancel = threading.Event()

    def cancel(self):
        self._cancel.set()

    def _check(self):
        if self._cancel.is_set():
            raise Cancelled()

    def _say(self, text):
        self.log.append(text)

    def _step(self, text):
        self.status = text
        self._check()

    def check_options(self):
        if self.forcefield not in mm.FORCEFIELD_INFO:
            raise ValueError('forcefield must be one of: ' +
                             ', '.join(mm.FORCEFIELD_INFO))
        info = mm.FORCEFIELD_INFO[self.forcefield]
        if self.solvent not in SOLVENTS:
            raise ValueError('solvent must be one of: ' + ', '.join(SOLVENTS))
        if self.restrain not in RESTRAINTS:
            raise ValueError('restrain must be one of: ' +
                             ', '.join(RESTRAINTS))
        if self.altloc != 'occupancy' and (len(str(self.altloc)) != 1 or
                                           not str(self.altloc).isalnum()):
            raise ValueError('altloc must be "occupancy" or one letter')
        if self.pka_method not in PKA_METHODS:
            raise ValueError('pka_method must be propka or standard')
        if self.minimize and self.solvent in mm.IMPLICIT_SOLVENTS and \
                not info['implicit']:
            raise ValueError('%s has no implicit solvent here: use explicit '
                             'or vacuum' % info['label'])
        if self.minimize and self.preprocess and self.cap_termini and \
                not info['caps']:
            raise ValueError('%s has no ACE/NME caps here: use an AMBER force '
                             'field or leave the termini uncapped'
                             % info['label'])
        if self.protonate and self.pka_method == 'propka' and \
                not propka_available():
            raise ValueError('PROPKA is not installed (pip install propka); '
                             'use standard pKa values')

    def run(self):
        '''
        Preprocess, protonate and minimize (no PyMOL calls).

        :raises Cancelled: if cancel() was called
        '''
        _, app, _, PDBFixer = mm._require()
        self.check_options()
        self._step('Reading the structure...')
        if self.alternates:
            self._say_alternates()
        fixer = self._fixer = PDBFixer(
            pdbxfile=io.StringIO(self.input_text))
        fixer.missingResidues = {}
        self._unknown_hydrogens = 0
        if self.input_hydrogens and not (self.preprocess and
                                         self.delete_hydrogens):
            self._keep_known_hydrogens(fixer, app)
        # atoms of the input: only these are restrained (built ones move)
        self.input_atoms = {self._atom_key(a)
                            for a in fixer.topology.atoms()}

        if self.preprocess:
            self._preprocess(fixer, app)
        self.topology, self.positions = fixer.topology, fixer.positions
        self._sites(app)

        if self.protonate:
            self._protonate(app)
        elif self.minimize:
            self._say('Hydrogens (standard pKa values, pH %.1f), needed for '
                      'the minimization' % self.ph)
            self._add_hydrogens(app, [None] * self.topology.getNumResidues(),
                                self.ph)

        if self.minimize:
            self._minimize()
        self.status = 'Done'
        return self.report

    def _keep_known_hydrogens(self, fixer, app):
        '''
        Existing hydrogens of polymer residues are kept if OpenMM knows
        their names (e.g. from an earlier preparation, or standard PDB
        names); others (e.g. PyMOL's "1HB") are deleted and added again.
        '''
        app.Modeller._loadStandardHydrogenDefinitions()
        known = {name: {h.name for h in data.hydrogens}
                 for name, data in app.Modeller._residueHydrogens.items()}
        # (terminal hydrogens of amino acids)
        terminal = {'H1', 'H2', 'H3', 'HXT'}
        polymer = set(self.polymer_ids.values())
        unknown = [a for a in fixer.topology.atoms()
                   if a.element is not None and a.element.symbol == 'H' and
                   a.residue.chain.id in polymer and
                   a.name not in known.get(a.residue.name, set()) | terminal]
        if unknown:
            modeller = app.Modeller(fixer.topology, fixer.positions)
            modeller.delete(unknown)
            fixer.topology, fixer.positions = (modeller.topology,
                                               modeller.positions)
        self._unknown_hydrogens = len(unknown)

    def _label(self, residue):
        return '%s %s/%s' % (residue.name, _chain_label(self.chain_map.get(
            residue.chain.id, residue.chain.id)), mm._residue_id(residue))

    def _sites(self, app):
        '''
        Disulfide bonds (reported, or broken) and metal sites: what each
        metal ion coordinates (by distance), for the coordinating
        residues' protonation states, the minimization's restraints and
        zero-order bonds in the result
        '''
        import numpy as np
        from openmm import unit
        from scipy.spatial import cKDTree

        say = self._say
        lines = []
        bonds = [b for b in self.topology.bonds()
                 if b.atom1.name == 'SG' and b.atom2.name == 'SG' and
                 b.atom1.residue.name == b.atom2.residue.name == 'CYS']
        if bonds and not self.disulfides:
            modeller = app.Modeller(self.topology, self.positions)
            modeller.delete(bonds)
            self.topology, self.positions = (modeller.topology,
                                             modeller.positions)
            lines.append('  disulfide bonds broken (free cysteines): %d' %
                         len(bonds))
            bonds = []
        elif bonds:
            pairs = sorted((sorted((b.atom1.residue, b.atom2.residue),
                                   key=lambda r: r.index) for b in bonds),
                           key=lambda pair: pair[0].index)
            lines.append('  disulfide bonds: %d (%s)' % (len(bonds), ', '.join(
                '%s-%s' % (self._label(a), self._label(b))
                for a, b in pairs)))
        self.report['disulfides'] = len(bonds)

        xyz = np.array(self.positions.value_in_unit(unit.angstrom))
        atoms = list(self.topology.atoms())
        # close cysteines without a bond
        bonded = {a for b in bonds for a in (b.atom1, b.atom2)}
        sg = [a for a in atoms if a.name == 'SG' and
              a.residue.name == 'CYS' and a not in bonded]
        close = []
        for i, a in enumerate(sg):
            for b in sg[i + 1:]:
                d = float(np.linalg.norm(xyz[a.index] - xyz[b.index]))
                if d < 3.5:
                    close.append('%s-%s %.1f A' % (
                        self._label(a.residue), self._label(b.residue), d))
        if close:
            lines.append('  close cysteines without a disulfide bond (SG-SG):'
                         ' %s' % ', '.join(close))

        self.metal_bonds = []
        self._metal_pairs = []
        self._metal_states = {}
        polymer = set(self.polymer_ids.values())
        metals = [a for a in atoms if a.element is not None and
                  a.element.symbol.upper() in METALS and
                  a.residue.chain.id not in polymer]
        if self.metals and metals:
            ligands = [a for a in atoms if a.element is not None and
                       a.element.symbol in ('N', 'O', 'S')]
            tree = cKDTree(xyz[[a.index for a in ligands]])
            sites = {}
            for metal in metals:
                cutoff = METALS[metal.element.symbol.upper()]
                contacts = []
                for k in tree.query_ball_point(xyz[metal.index], cutoff):
                    other = ligands[k]
                    if other.residue is metal.residue:
                        continue
                    d = float(np.linalg.norm(xyz[other.index] -
                                             xyz[metal.index]))
                    contacts.append((other, d))
                    self._metal_pairs.append((self._atom_key(metal),
                                              self._atom_key(other)))
                    self.metal_bonds.append((self._author_key(metal),
                                             self._author_key(other), 0))
                    state = METAL_STATES.get((other.residue.name,
                                              other.name))
                    if state:
                        self._metal_states[other.residue] = state
                if contacts:
                    label = self._label(metal.residue)
                    sites[label] = len(contacts)
                    lines.append('  metal %s: %s' % (label, ', '.join(
                        '%s %s %.2f A' % ('water' if other.residue.name
                                          in mm.WATER else
                                          self._label(other.residue),
                                          other.name, d)
                        for other, d in sorted(contacts,
                                               key=lambda c: c[1]))))
            self.report['metal_sites'] = sites
        if lines:
            say('Disulfides and metal sites')
            for line in lines:
                say(line)

    def _author_key(self, atom):
        '''
        (chain, resi, name) as in PyMOL
        '''
        return (self.chain_map.get(atom.residue.chain.id,
                                   atom.residue.chain.id),
                mm._residue_id(atom.residue), atom.name)

    def _metal_dashes(self, _self):
        dashes = _self.get_legal_name(self.name + '_metal')
        _self.delete(dashes)
        for a, b, _ in self.metal_bonds:
            _self.distance(dashes, _select_atom(self.name, a),
                           _select_atom(self.name, b), quiet=1)
        _self.hide('labels', dashes)
        self.report['metal_dashes'] = dashes

    def _say_alternates(self):
        kept = {}
        for _, _, _, letter, _ in self.alternates:
            kept[letter] = kept.get(letter, 0) + 1
        self._say('Alternate locations: %d residues, kept %s (%s)' % (
            len(self.alternates), ', '.join(
                '%s in %d' % item for item in sorted(kept.items())),
            'highest occupancy' if self.altloc == 'occupancy' else
            '%s where present, else highest occupancy' % self.altloc))
        self.report['alternates'] = {
            '%s/%s %s' % (_chain_label(chain), resn, resi): letter
            for chain, resi, resn, letter, _ in self.alternates}

    @staticmethod
    def _atom_key(atom):
        return (atom.residue.chain.id, mm._residue_id(atom.residue),
                atom.residue.name, atom.name)

    # preprocess

    def _polymer_chains(self, fixer):
        polymer = set(self.polymer_ids.values())
        return [(index, chain) for index, chain in
                enumerate(fixer.topology.chains()) if chain.id in polymer]

    def _preprocess(self, fixer, app):
        self._step('Preprocessing...')
        say = self._say
        say('Preprocess')
        polymer = set(self.polymer_ids.values())
        if self.input_hydrogens and self.delete_hydrogens:
            say('  existing hydrogens deleted: %d' % self.input_hydrogens)
        elif self.input_hydrogens:
            say('  existing hydrogens kept: %d%s' % (
                self.input_hydrogens - self._unknown_hydrogens,
                ' (%d with names OpenMM doesn\'t know replaced)' %
                self._unknown_hydrogens if self._unknown_hydrogens else ''))

        if self.delete_waters or self.delete_heterogens or \
                self.water_min_hbonds:
            waters = self._doomed_waters(fixer, polymer)
            doomed = waters + [
                r for r in fixer.topology.residues()
                if r.chain.id not in polymer and r.name not in mm.WATER and
                self.delete_heterogens]
            if doomed:
                modeller = app.Modeller(fixer.topology, fixer.positions)
                modeller.delete(doomed)
                fixer.topology, fixer.positions = (modeller.topology,
                                                   modeller.positions)
            waters = len(waters)
            if self.delete_waters or self.water_min_hbonds:
                why = []
                if self.delete_waters and self.water_distance:
                    why.append('beyond %g A from het groups' %
                               self.water_distance)
                if self.water_min_hbonds:
                    why.append('with fewer than %d hydrogen bonds to '
                               'non-waters' % self.water_min_hbonds)
                say('  waters deleted: %d%s' % (
                    waters, ' (%s)' % '; '.join(why) if why else ''))
            if self.delete_heterogens:
                say('  other heterogens deleted: %d%s' % (
                    len(doomed) - waters, ' (%s)' % ', '.join(sorted(
                        {r.name for r in doomed if r.name not in mm.WATER}))
                    if len(doomed) > waters else ''))

        mutated = set()
        if self.mutations:
            for chain, resi, old, new in self.mutations:
                fixer.applyMutations(['%s-%s-%s' % (old, int(resi), new)],
                                     self.polymer_ids[chain])
            mutated = {(self.polymer_ids[chain], str(int(resi)))
                       for chain, resi, _, _ in self.mutations}
            say('  mutations: %s' % ', '.join(
                '%s/%s %s -> %s' % (_chain_label(chain), old, resi, new)
                for chain, resi, old, new in self.mutations))
            self.report['mutations'] = len(self.mutations)

        if self.replace_nonstandard:
            fixer.findNonstandardResidues()
            replaced = [r for r, _ in fixer.nonstandardResidues
                        if r.chain.id in polymer]
            fixer.nonstandardResidues = [
                (r, new) for r, new in fixer.nonstandardResidues
                if r.chain.id in polymer]
            say('  nonstandard residues replaced: %s' % (', '.join(sorted(
                {'%s->%s' % (r.name, new)
                 for r, new in fixer.nonstandardResidues})) or 'none'))
            if replaced:
                fixer.replaceNonstandardResidues()
        self._check()

        # missing residues: from the structure (or a sequence source)
        missing = {}
        if self.fill_loops:
            missing = self._missing_residues(fixer)
            chains = dict(self._polymer_chains(fixer))
            sizes = {index: len(list(chain.residues()))
                     for index, chain in chains.items()}
            loops = {key: value for key, value in missing.items()
                     if 0 < key[1] < sizes.get(key[0], 0)}
            ends = {key: value for key, value in missing.items()
                    if key not in loops}
            say('  missing residues filled in: %d in %d gaps' % (
                sum(map(len, loops.values())), len(loops)))
            if self.fill_termini:
                say('  missing residues at chain ends filled in: %d' %
                    sum(map(len, ends.values())))
            else:
                if ends:
                    say('  missing residues at chain ends left out: %d '
                        '(floppy; "including chain ends" to build them)' %
                        sum(map(len, ends.values())))
                missing = loops

        # ACE/NME caps at the ends of protein chains
        capped = []
        if self.cap_termini:
            capped, missing = self._caps(fixer, app, missing)
            say('  termini capped (ACE/NME): %s' % (
                ', '.join(capped) or 'none'))
        fixer.missingResidues = missing

        self._step('Building missing atoms...')
        fixer.findMissingAtoms()
        side_chain_atoms = sum(len(v) for v in fixer.missingAtoms.values())
        if not self.fill_side_chains:
            # (mutated residues get their new side chains anyway)
            fixer.missingAtoms = {
                r: atoms for r, atoms in fixer.missingAtoms.items()
                if (r.chain.id, r.id) in mutated}
        if self.cap_termini:
            ends = {residue for _, chain in self._polymer_chains(fixer)
                    for residue in list(chain.residues())[-1:]}
            fixer.missingTerminals = {
                r: atoms for r, atoms in fixer.missingTerminals.items()
                if r not in ends or r.name not in AMINO_ACIDS}
        if self.fill_side_chains:
            say('  missing side-chain atoms filled in: %d in %d residues' % (
                side_chain_atoms, len(fixer.missingAtoms)))
        elif side_chain_atoms:
            say('  missing side-chain atoms left out: %d' % side_chain_atoms)
        self._build(fixer)
        self.report['missing_residues'] = sum(map(len, missing.values()))
        self.report['missing_atoms'] = side_chain_atoms
        self.report['capped'] = capped

    def _built(self, topology):
        '''
        :return: set of heavy atoms which aren't in the input
        '''
        return {a for a in topology.atoms()
                if a.element is not None and a.element.symbol != 'H' and
                self._atom_key(a) not in self.input_atoms}

    def _build_problems(self, topology, positions):
        '''
        D-amino acids and cis peptides involving built atoms (PDBFixer
        builds missing residues roughly)
        '''
        from openmm import unit
        from .diagnostics import _chirality, _dihedral

        built = self._built(topology)
        xyz = positions.value_in_unit(unit.angstrom)
        problems = []
        for residue, named in _backbone(topology):
            if residue.name in AMINO_ACIDS and 'CB' in named and \
                    built & {named[k] for k in ('N', 'CA', 'C', 'CB')} and \
                    _chirality({k: {'coord': tuple(xyz[a.index])}
                                for k, a in named.items()}) < -0.5:
                problems.append('D-%s %s' % (residue.name, residue.id))
        for a, b in _peptide_bonds(topology):
            atoms = (a['CA'], a['C'], b['N'], b['CA'])
            if built & set(atoms) and abs(_dihedral(
                    *(tuple(xyz[atom.index]) for atom in atoms))) < 90:
                problems.append('cis %s-%s' % (a['C'].residue.id,
                                               b['N'].residue.id))
        from .diagnostics import COVALENT
        for bond in topology.bonds():
            a, b = bond.atom1, bond.atom2
            if not (built & {a, b}) or a.element is None or \
                    b.element is None:
                continue
            ea, eb = a.element.symbol.upper(), b.element.symbol.upper()
            if 'H' in (ea, eb) or ea not in COVALENT or eb not in COVALENT:
                continue
            d = math.dist(xyz[a.index], xyz[b.index])
            ideal = COVALENT[ea] + COVALENT[eb]
            if d > ideal + 0.45 or d < ideal - 0.4:
                problems.append('bond %s %s %s-%s %.2f A' % (
                    a.residue.name, a.residue.id, a.name, b.name, d))
        return problems

    def _build(self, fixer):
        '''
        Add missing atoms and residues (one seed: reproducible). PDBFixer
        builds some residues as D-amino acids: their side chains are
        mirrored here; cis peptides of built residues are turned trans in
        the minimization.
        '''
        fixer.addMissingAtoms(seed=1)
        fixer.positions = self._fix_chirality(fixer.topology,
                                              fixer.positions)
        problems = self._build_problems(fixer.topology, fixer.positions)
        if problems and not self.minimize:
            self._say('  built residues with problems (minimize to correct '
                      'them): %s' % ', '.join(problems))

    def _fix_chirality(self, topology, positions):
        '''
        Built D-amino acids made L: the side chain is reflected through
        the N-CA-C plane (then CG2 of Ile/Thr back through its CB plane,
        whose chirality the reflection inverted too). Heavy atoms only:
        hydrogens come later.
        '''
        import numpy as np
        from openmm import unit
        from .diagnostics import _chirality

        built = self._built(topology)
        xyz = np.array(positions.value_in_unit(unit.nanometer))

        def reflect(indices, a, b, c):
            normal = np.cross(xyz[b] - xyz[a], xyz[c] - xyz[a])
            normal /= np.linalg.norm(normal)
            for i in indices:
                xyz[i] -= 2 * np.dot(xyz[i] - xyz[a], normal) * normal

        fixed = 0
        for residue, named in _backbone(topology):
            if residue.name not in AMINO_ACIDS or 'CB' not in named or \
                    not built & {named[k] for k in ('N', 'CA', 'C', 'CB')}:
                continue
            coords = {k: {'coord': tuple(xyz[a.index])}
                      for k, a in named.items()}
            if _chirality(coords) >= 0:
                continue
            side = [a.index for a in residue.atoms()
                    if a.name not in ('N', 'CA', 'C', 'O', 'OXT') and
                    a.element is not None and a.element.symbol != 'H']
            reflect(side, named['N'].index, named['CA'].index,
                    named['C'].index)
            atoms = {a.name: a.index for a in residue.atoms()}
            other = {'ILE': 'CG1', 'THR': 'OG1'}.get(residue.name)
            if other and {'CB', 'CG2', other} <= set(atoms):
                reflect([atoms['CG2']], atoms['CA'], atoms['CB'],
                        atoms[other])
            fixed += 1
        if fixed:
            self._say('  built D-amino acids made L: %d' % fixed)
        return unit.Quantity([openmm_vec(v) for v in xyz], unit.nanometer)

    def _doomed_waters(self, fixer, polymer):
        '''
        Waters to delete: all, or those beyond water_distance from het
        groups (non-water, non-polymer residues), and those with fewer
        than water_min_hbonds hydrogen bonds to non-waters (N or O 2.5 to
        3.3 A from the water oxygen)
        '''
        import numpy as np
        from openmm import unit

        xyz = np.array(fixer.positions.value_in_unit(unit.angstrom))
        waters = [r for r in fixer.topology.residues()
                  if r.chain.id not in polymer and r.name in mm.WATER]
        oxygen = {r: [a.index for a in r.atoms()
                      if a.element is not None and a.element.symbol == 'O']
                  for r in waters}
        waters = [r for r in waters if oxygen[r]]
        doomed = set()
        if self.delete_waters:
            het = [a.index for a in fixer.topology.atoms()
                   if a.residue.chain.id not in polymer and
                   a.residue.name not in mm.WATER]
            if self.water_distance > 0 and het:
                het_xyz = xyz[het]
                for r in waters:
                    d = np.linalg.norm(het_xyz - xyz[oxygen[r][0]], axis=1)
                    if d.min() > self.water_distance:
                        doomed.add(r)
            else:
                doomed.update(waters)
        if self.water_min_hbonds:
            polar = [a.index for a in fixer.topology.atoms()
                     if a.residue.name not in mm.WATER and
                     a.element is not None and a.element.symbol in 'NO']
            polar_xyz = xyz[polar] if polar else np.zeros((0, 3))
            for r in waters:
                if r in doomed:
                    continue
                d = np.linalg.norm(polar_xyz - xyz[oxygen[r][0]], axis=1)
                if np.count_nonzero((d >= 2.5) & (d <= 3.3)) < \
                        self.water_min_hbonds:
                    doomed.add(r)
        return [r for r in waters if r in doomed]

    def _missing_residues(self, fixer):
        '''
        :return: PDBFixer missingResidues {(chain index, residue index):
            [names]}
        '''
        if self.source:
            fixer.sequences = mm._sequences_from_source(self.source)
            for sequence in fixer.sequences:
                sequence.chainId = self.polymer_ids.get(sequence.chainId,
                                                        '-')
            fixer.findMissingResidues()
            polymer = {index for index, _ in self._polymer_chains(fixer)}
            return {key: value for key, value in fixer.missingResidues.items()
                    if key[0] in polymer}

        missing = {}
        for index, chain in self._polymer_chains(fixer):
            residues = self.residues.get(self.chain_map[chain.id], [])
            if all(present for _, _, present in residues):
                continue
            position = {mm._residue_id(r): i
                        for i, r in enumerate(chain.residues())}
            pending = []
            for resi, resn, present in residues:
                if not present:
                    pending.append(resn)
                elif resi in position:
                    if pending:
                        missing[(index, position[resi])] = pending
                    pending = []
            if pending:
                missing[(index, len(position))] = pending
        return missing

    def _caps(self, fixer, app, missing):
        missing = {key: list(value) for key, value in missing.items()}
        capped = []
        oxt = []
        for index, chain in self._polymer_chains(fixer):
            residues = list(chain.residues())
            label = _chain_label(self.chain_map[chain.id])
            ends = []
            first, last = residues[0], residues[-1]
            if first.name in AMINO_ACIDS:
                missing[(index, 0)] = ['ACE'] + missing.get((index, 0), [])
                ends.append('N')
            if last.name in AMINO_ACIDS:
                key = (index, len(residues))
                missing[key] = missing.get(key, []) + ['NME']
                oxt += [a for a in last.atoms() if a.name == 'OXT']
                ends.append('C')
            if ends:
                capped.append('%s (%s)' % (label, '+'.join(ends)))
        if oxt:
            # the C-terminal carboxylate becomes a peptide bond to NME
            modeller = app.Modeller(fixer.topology, fixer.positions)
            modeller.delete(oxt)
            fixer.topology, fixer.positions = (modeller.topology,
                                               modeller.positions)
        return capped, missing

    # protonation

    def _pkas(self):
        '''
        PROPKA pKa values of the polymer chains' titratable residues.

        :return: {residue: pKa}
        '''
        import propka.run
        from openmm import unit

        polymer = set(self.polymer_ids.values())
        letters = {}
        by_key = {}
        lines = []
        xyz = self.positions.value_in_unit(unit.angstrom)
        serial = 0
        for chain in self.topology.chains():
            if chain.id not in polymer:
                continue
            letter = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'[len(letters) % 36]
            letters[chain.id] = letter
            for residue in chain.residues():
                if residue.name not in AMINO_ACIDS:
                    continue
                number = int(residue.id)
                icode = residue.insertionCode.strip()[:1]
                by_key[(letter, number, icode or ' ')] = residue
                for atom in residue.atoms():
                    if atom.element is None or atom.element.symbol == 'H':
                        continue
                    serial += 1
                    name = atom.name if len(atom.name) == 4 else \
                        ' %-3s' % atom.name
                    x, y, z = xyz[atom.index]
                    lines.append(
                        'ATOM  %5d %-4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f'
                        '  1.00  0.00          %2s' % (
                            serial % 100000, name, residue.name, letter,
                            number % 10000, icode or ' ', x, y, z,
                            atom.element.symbol))
        if not lines:
            return {}
        lines.append('END')
        # (PROPKA logs warnings about every unusual group)
        logger = logging.getLogger('propka')
        level = logger.level
        logger.setLevel(logging.CRITICAL)
        try:
            molecule = propka.run.single(
                'pymolx.pdb', optargs=['--quiet'],
                stream=io.StringIO('\n'.join(lines) + '\n'),
                write_pka=False)
        finally:
            logger.setLevel(level)
        conformation = molecule.conformations[molecule.conformation_names[0]]
        pkas = {}
        for group in conformation.groups:
            if group.residue_type not in TITRATABLE:
                continue
            atom = group.atom
            key = (atom.chain_id, atom.res_num,
                   (getattr(atom, 'icode', '') or ' ').strip() or ' ')
            residue = by_key.get(key)
            if residue is not None and residue.name == group.residue_type:
                pkas[residue] = group.pka_value
        return pkas

    def _ligand_variants(self, app):
        '''
        Hydrogens of residues OpenMM has no hydrogen definitions for
        (ligands, kept nonstandard residues), from the PDB's Chemical
        Component Dictionary (downloaded by PDBFixer; its private helpers)

        :return: {residue index: [(hydrogen name, parent name)]}
        '''
        app.Modeller._loadStandardHydrogenDefinitions()
        known = app.Modeller._residueHydrogens
        residues = [r for r in self.topology.residues()
                    if r.name not in known and r.name not in mm.WATER]
        if not residues:
            return {}
        self._step('Getting ligand definitions (PDB Chemical Component '
                   'Dictionary)...')
        fixer = self._fixer
        fixer.topology, fixer.positions = self.topology, self.positions
        try:
            definitions = fixer._downloadNonstandardDefinitions()
            variants = {r.index: fixer._describeVariant(r, definitions)
                        for r in residues}
        except Exception as e:
            self._say('  ligand hydrogens not added: %s' % e)
            return {}
        added = {}
        missing = set()
        for r in residues:
            variant = variants[r.index]
            if variant is None:
                missing.add(r.name)
            else:
                added.setdefault(r.name, []).append(len(variant))
        if added:
            self._say('  ligand hydrogens (CCD): %s' % ', '.join(
                '%s %d%s' % (name, n[0], ' x%d' % len(n) if len(n) > 1
                             else '') for name, n in sorted(added.items())))
        if missing:
            self._say('  no hydrogens for %s (no matching CCD definition, '
                      'or no internet)' % ', '.join(sorted(missing)))
        self.report['ligand_hydrogens'] = {
            name: sum(n) for name, n in added.items()}
        return {i: v for i, v in variants.items() if v is not None}

    def _flip(self, app, residues, variants, pkas, fixed=None):
        '''
        Asn/Gln/His flips and His tautomers (pymolx.hbonds) before the
        hydrogens are added; His variants are set in "variants"

        :return: ["Asn A/45 flipped", ...]
        '''
        from . import hbonds

        charged = {}
        for residue in residues:
            if residue.name == 'HIS':
                pka = pkas.get(residue, 6.5)  # (standard: Modeller's rule)
                charged[residue] = pka > self.ph if residue in pkas else \
                    self.ph <= 6.5
        states = {r: v for r, v in zip(residues, variants) if v}
        built = {a.residue for a in self._built(self.topology)}
        fixed = {r: s for r, s in (fixed or {}).items() if r.name == 'HIS'}
        self.positions, tautomers, flipped = hbonds.optimize_flips(
            self.topology, self.positions, charged, states,
            built | set(fixed), fixed)
        index = {r: i for i, r in enumerate(residues)}
        for residue, tautomer in tautomers.items():
            variants[index[residue]] = tautomer
        # kept hydrogens of flipped groups would be in the wrong places
        names = {'ASN': ('HD21', 'HD22'), 'GLN': ('HE21', 'HE22'),
                 'HIS': ('HD1', 'HE2', 'HD2', 'HE1')}
        stale = [a for r in flipped for a in r.atoms()
                 if a.name in names[r.name]]
        if stale:
            modeller = app.Modeller(self.topology, self.positions)
            modeller.delete(stale)
            self.topology, self.positions = (modeller.topology,
                                             modeller.positions)
        return ['%s %s/%s flipped' % (
            r.name.capitalize(), _chain_label(self.chain_map.get(
                r.chain.id, r.chain.id)), mm._residue_id(r))
            for r in flipped]

    def _add_hydrogens(self, app, variants, ph):
        if self.ligand_hydrogens:
            variants = list(variants)
            for index, variant in self._ligand_variants(app).items():
                if variants[index] is None:
                    variants[index] = variant
        modeller = app.Modeller(self.topology, self.positions)
        actual = modeller.addHydrogens(pH=ph, variants=variants)
        self.topology, self.positions = modeller.topology, modeller.positions
        return actual

    def _protonate(self, app):
        say = self._say
        residues = list(self.topology.residues())
        variants = [None] * len(residues)
        pkas = {}
        skipped = []
        if self.pka_method == 'propka':
            self._step('Computing pKa values (PROPKA)...')
            pkas = self._pkas()
            # side chains built here are placed arbitrarily: their PROPKA
            # pKa means little, standard values are used
            built = {a.residue for a in self._built(self.topology)
                     if a.name not in ('N', 'CA', 'C', 'O', 'OXT')}
            skipped = sorted((r for r in pkas if r in built),
                             key=lambda r: r.index)
            for residue in skipped:
                del pkas[residue]
            for i, residue in enumerate(residues):
                if residue in pkas:
                    protonated, deprotonated = TITRATABLE[residue.name]
                    variants[i] = (protonated if pkas[residue] > self.ph
                                   else deprotonated)
        metal_states = {}
        for i, residue in enumerate(residues):
            state = self._metal_states.get(residue)
            if state:
                variants[i] = state
                metal_states[residue] = state
        say('Protonation at pH %.1f (%s)' % (
            self.ph, 'PROPKA 3 pKa values' if pkas else
            'standard pKa values'))
        if metal_states:
            say('  metal-coordinating residues: %s' % ', '.join(
                '%s %s' % (self._label(r), 'CYM' if s == 'CYX' else s)
                for r, s in metal_states.items()))
        if self.pka_method == 'propka' and skipped:
            say('  standard pKa values for built side chains: %s' % ', '.join(
                '%s %s/%s' % (r.name, _chain_label(self.chain_map.get(
                    r.chain.id, r.chain.id)), mm._residue_id(r))
                for r in skipped))
        flips = []
        if self.optimize_hbonds:
            self._step('Optimizing hydrogen bonds (flips)...')
            flips = self._flip(app, residues, variants, pkas, metal_states)
        self._step('Adding hydrogens...')
        # (Modeller picks the His tautomer from hydrogen bonds only above
        # pH 6.5; PROPKA has decided which His are charged already)
        model_ph = max(self.ph, 7.0) if pkas else self.ph
        self._add_hydrogens(app, variants, model_ph)
        # (from the hydrogens: Modeller reports no variant for residues it
        # decided itself)
        disulfide = {a.residue for b in self.topology.bonds()
                     for a in (b.atom1, b.atom2)
                     if b.atom1.name == b.atom2.name == 'SG'}
        actual = [_protonation_state(r, r in disulfide)
                  for r in self.topology.residues()]
        if self.optimize_hbonds:
            self._step('Optimizing hydrogen bonds (rotatable hydrogens)...')
            from . import hbonds
            states = dict(zip(self.topology.residues(), actual))
            self.positions, moved, rotors = \
                hbonds.orient_rotatable_hydrogens(self.topology,
                                                  self.positions, states)
            say('  hydrogen bonds optimized: %s; rotatable hydrogens '
                'turned: %d of %d (OH, SH, NH3+, waters)' % (
                    ', '.join(flips) or 'no flips', moved, rotors))
            self.report['flips'] = flips
            self.report['rotated'] = moved
        names = [r.name for r in residues]
        states = {}
        for name, variant in zip(names, actual):
            if name in ('ASP', 'GLU', 'HIS', 'LYS', 'CYS'):
                states.setdefault(name, {}).setdefault(variant or name, 0)
                states[name][variant or name] += 1
        for name in sorted(states):
            say('  %s: %s' % (name, ', '.join(
                '%d %s' % (n, v) for v, n in sorted(states[name].items()))))
        unusual = []
        for residue, variant in zip(residues, actual):
            if residue not in pkas:
                continue
            # states the residue's usual (model) pKa wouldn't give
            usual = {'ASP': 'ASP', 'GLU': 'GLU', 'LYS': 'LYS'}.get(
                residue.name)
            if residue.name == 'HIS':
                usual = 'HIP' if self.ph < 6.5 else None
                odd = (variant == 'HIP') != (usual == 'HIP')
            else:
                odd = variant != usual
            if odd:
                unusual.append('%s %s%s (pKa %.1f): %s' % (
                    residue.name, _chain_label(self.chain_map.get(
                        residue.chain.id, residue.chain.id)) + '/',
                    mm._residue_id(residue), pkas[residue], variant))
        if pkas:
            say('  unusual protonation states: %s' % (
                '; '.join(unusual) if unusual else 'none'))
        charge = 0
        polymer = set(self.polymer_ids.values())
        for residue, variant in zip(residues, actual):
            if residue.chain.id in polymer or residue.name in CHARGE:
                charge += CHARGE.get(variant or residue.name, 0)
        for chain in self.topology.chains():
            if chain.id not in polymer:
                continue
            chain_residues = list(chain.residues())
            if chain_residues[0].name in AMINO_ACIDS:
                charge += 1
            if chain_residues[-1].name in AMINO_ACIDS:
                charge -= 1
        say('  net charge (protein and ions): %+d' % charge)
        self.report['pka'] = {(self.chain_map.get(r.chain.id, r.chain.id),
                               mm._residue_id(r), r.name): v
                              for r, v in pkas.items()}
        self.report['variants'] = states
        self.report['unusual'] = unusual
        self.report['charge_estimate'] = charge

    # minimization

    def _context(self, openmm, system, positions):
        names = mm.PLATFORMS if self.platform == 'auto' else (self.platform,)
        errors = []
        for name in names:
            try:
                platform = openmm.Platform.getPlatformByName(name)
                integrator = openmm.VerletIntegrator(0.001)
                # mixed precision: in single precision the minimizer's line
                # search could stop at once
                properties = {'Precision': 'mixed'} if name in (
                    'CUDA', 'OpenCL') else {}
                context = openmm.Context(system, integrator, platform,
                                         properties)
                context.setPositions(positions)
                return context, name
            except Exception as e:
                errors.append('%s: %s' % (name, e))
        raise RuntimeError('no OpenMM platform worked: ' + '; '.join(errors))

    def _descend(self, openmm, system, platform, positions, steps=500):
        '''
        Steepest descent (at most 0.02 A per step and coordinate; massless
        particles stay fixed), in a context of its own.

        :return: positions, nm
        '''
        integrator = openmm.CustomIntegrator(0.0)
        integrator.addGlobalVariable('h', 1e-6)  # nm per kJ/mol/nm
        integrator.addComputePerDof(
            'x', 'x + select(m, max(-0.002, min(0.002, h*f)), 0)')
        context = openmm.Context(system, integrator,
                                 openmm.Platform.getPlatformByName(platform))
        context.setPositions(positions)
        for _ in range(steps // 50):
            integrator.step(50)
            self._check()
        return context.getState(getPositions=True).getPositions(asNumpy=True)

    def _minimize(self):
        openmm, app, unit, _ = mm._require()
        say = self._say
        info = mm.FORCEFIELD_INFO[self.forcefield]
        self._step('Setting up the minimization...')

        water_file, water_model = info['water']
        files = list(info['files']) + [water_file]
        if self.solvent in mm.IMPLICIT_SOLVENTS:
            files.append('implicit/%s.xml' % self.solvent)
        forcefield = app.ForceField(*files)

        # residues without parameters (ligands) stay, but outside the system
        topology, positions = self.topology, self.positions
        # (only non-polymer residues: getUnmatchedResidues ignores CHARMM's
        # terminal patches; createSystem reports polymer residues it can't
        # handle)
        polymer = set(self.polymer_ids.values())
        excluded = [r for r in forcefield.getUnmatchedResidues(topology)
                    if r.chain.id not in polymer]
        excluded_atoms = {a.index for r in excluded for a in r.atoms()}
        kept = [a.index for a in topology.atoms()
                if a.index not in excluded_atoms]
        if excluded:
            modeller = app.Modeller(topology, positions)
            modeller.delete(excluded)
            topology, positions = modeller.topology, modeller.positions
        n_structure = topology.getNumAtoms()

        if self.solvent == 'explicit':
            self._step('Adding water...')
            modeller = app.Modeller(topology, positions)
            modeller.addSolvent(forcefield, model=water_model,
                                padding=EXPLICIT_PADDING * unit.nanometer,
                                ionicStrength=IONIC_STRENGTH * unit.molar,
                                neutralize=True)
            topology, positions = modeller.topology, modeller.positions
            system = forcefield.createSystem(
                topology, nonbondedMethod=app.PME,
                nonbondedCutoff=1.0 * unit.nanometer, rigidWater=False)
        else:
            many = topology.getNumAtoms() > 5000
            system = forcefield.createSystem(
                topology,
                nonbondedMethod=app.CutoffNonPeriodic if many else
                app.NoCutoff, rigidWater=False,
                **({'nonbondedCutoff': 2.0 * unit.nanometer} if many else {}))
        self._check()

        charge = 0.0
        for force in system.getForces():
            force.setForceGroup(0)
            if isinstance(force, openmm.NonbondedForce):
                charge = sum(force.getParticleParameters(i)[0].value_in_unit(
                    unit.elementary_charge)
                    for i in range(n_structure))

        restrained = 0
        if self.restrain != 'none':
            periodic = self.solvent == 'explicit'
            restraint = openmm.CustomExternalForce(
                '0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2' if periodic
                else '0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
            restraint.addGlobalParameter('k', self.restraint_k)
            for name in ('x0', 'y0', 'z0'):
                restraint.addPerParticleParameter(name)
            for atom in list(topology.atoms())[:n_structure]:
                if atom.element is None or atom.element.symbol == 'H':
                    continue
                if self.restrain == 'backbone' and atom.name not in BACKBONE:
                    continue
                if self.restrain == 'ca' and atom.name != 'CA':
                    continue
                if self._atom_key(atom) not in self.input_atoms:
                    continue  # built: free to move
                restraint.addParticle(atom.index, positions[atom.index]
                                      .value_in_unit(unit.nanometer))
                restrained += 1
            restraint.setForceGroup(1)
            system.addForce(restraint)

        # metal sites: the metal and the atoms it binds stay where they are
        # (a nonbonded ion model distorts them: its van der Waals sphere
        # pushes partners away, without it charges collapse onto it)
        site_atoms = set()
        if self._metal_pairs:
            index = {self._atom_key(a): a.index for a in topology.atoms()}
            for a, b in self._metal_pairs:
                if a in index and b in index:
                    site_atoms |= {index[a], index[b]}
            for i in site_atoms:
                system.setParticleMass(i, 0.0)
            self.report['metal_fixed_atoms'] = len(site_atoms)

        # first stage: input heavy atoms held fixed (massless), so that
        # hydrogens, built atoms and solvent clashing at the start can't
        # push them around (e.g. invert a chiral center)
        fixed = {}
        for atom in list(topology.atoms())[:n_structure]:
            if atom.element is not None and atom.element.symbol != 'H' and \
                    self._atom_key(atom) in self.input_atoms and \
                    atom.index not in site_atoms:
                fixed[atom.index] = system.getParticleMass(atom.index)
                system.setParticleMass(atom.index, 0.0)
        # built residues and the ones joining them are guided to L
        # chirality and trans peptide bonds, in both stages (PDBFixer builds
        # residues roughly, and a strained loop may pull its neighbors)
        guide = None
        built = self._built(topology)
        if built:
            # residues with built atoms and their neighbors in the chain
            built_residues = {a.residue for a in built}
            near = set(built_residues)
            for a, b in _peptide_bonds(topology):
                if a['C'].residue in built_residues or \
                        b['N'].residue in built_residues:
                    near |= {a['C'].residue, b['N'].residue}
            # (harmonic in the angle difference, wrapped so that the cut is
            # opposite the target: an exactly cis peptide still gets a push,
            # unlike with 1 - cos)
            guide = openmm.CustomTorsionForce(
                '0.5*k_build*w^2; '
                'w = d - 6.28318531*floor((d + 3.14159265)/6.28318531); '
                'd = theta - theta0')
            guide.addGlobalParameter('k_build', BUILD_K)
            guide.addPerTorsionParameter('theta0')
            # (hydrogens too: they were added to the built geometry)
            for residue, named in _backbone(topology):
                if 'CB' in named and residue.name in AMINO_ACIDS and \
                        residue in near:
                    guide.addTorsion(*(named[k].index for k in
                                       ('N', 'C', 'CA', 'CB')),
                                     [L_IMPROPER])
                    if 'HA' in named:
                        guide.addTorsion(*(named[k].index for k in
                                           ('N', 'C', 'CA', 'HA')),
                                         [L_IMPROPER_HA])
            for a, b in _peptide_bonds(topology):
                atoms = (a['CA'], a['C'], b['N'], b['CA'])
                if a['C'].residue in near or b['N'].residue in near:
                    guide.addTorsion(*(atom.index for atom in atoms),
                                     [3.14159265])
                    if 'O' in a and 'H' in b:
                        guide.addTorsion(a['O'].index, a['C'].index,
                                         b['N'].index, b['H'].index,
                                         [3.14159265])
            if guide.getNumTorsions():
                guide.setForceGroup(2)
                system.addForce(guide)
            else:
                guide = None
        started = time.time()
        context, platform = self._context(openmm, system, positions)
        device = ''
        if platform in ('CUDA', 'OpenCL'):
            try:
                device = context.getPlatform().getPropertyValue(
                    context, 'DeviceName')
            except Exception:
                pass

        def energy():
            return context.getState(getEnergy=True, groups={0}) \
                .getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

        job = self

        class Reporter(openmm.MinimizationReporter):
            def report(self, iteration, x, grad, args):
                if job.progress is not None:
                    job.progress(iteration, args['system energy']
                                 if 'system energy' in args else 0.0)
                return job._cancel.is_set()

        def positions_now():
            return context.getState(getPositions=True).getPositions(
                asNumpy=True).value_in_unit(unit.nanometer)

        def lbfgs():
            before = positions_now()
            openmm.LocalEnergyMinimizer.minimize(
                context, self.tolerance * unit.kilojoule_per_mole /
                unit.nanometer, int(self.max_iterations), Reporter())
            self._check()
            return abs(positions_now() - before).max() > 1e-4

        def max_force():
            import numpy as np
            forces = context.getState(getForces=True).getForces(
                asNumpy=True).value_in_unit(unit.kilojoule_per_mole /
                                            unit.nanometer)
            mobile = [system.getParticleMass(i).value_in_unit(unit.dalton)
                      > 0 for i in range(system.getNumParticles())]
            return float(np.linalg.norm(forces[mobile], axis=1).max())

        def minimize():
            # OpenMM's minimizer stops when the root-mean-square force is
            # below the tolerance; in a large system that can happen with a
            # few atoms still strained, and sometimes it stops at once
            # (seen on CUDA). So: again while some atom's force is large,
            # after steepest descent steps if it didn't move.
            for _ in range(MAX_ROUNDS):
                if not lbfgs():
                    context.setPositions(self._descend(
                        openmm, system, platform, positions_now()))
                    self._check()
                    lbfgs()
                if max_force() < MAX_FORCE_FACTOR * self.tolerance or \
                        self.max_iterations:
                    break

        e0 = energy()
        if fixed:
            self._step('Minimizing hydrogens and built atoms...')
            minimize()
            relaxed = context.getState(getPositions=True).getPositions()
            for index, mass in fixed.items():
                system.setParticleMass(index, mass)
            # a new context: after context.reinitialize() the CUDA platform's
            # minimizer moved nothing (OpenMM 8.6)
            del context
            context, platform = self._context(openmm, system, relaxed)
        e_relaxed = energy()
        self._step('Minimizing...')
        minimize()
        e1 = energy()

        final = context.getState(getPositions=True).getPositions(asNumpy=True)
        final = final.value_in_unit(unit.nanometer)[:n_structure]
        start = self.positions.value_in_unit(unit.nanometer)
        new = list(start)
        squares = []
        atoms = list(self.topology.atoms())
        for index, xyz in zip(kept, final):
            new[index] = openmm.Vec3(*map(float, xyz))
            if atoms[index].element is not None and \
                    atoms[index].element.symbol != 'H' and \
                    self._atom_key(atoms[index]) in self.input_atoms:
                d = new[index] - start[index]
                squares.append(d[0]**2 + d[1]**2 + d[2]**2)
        self.positions = unit.Quantity(new, unit.nanometer)
        rmsd = 10 * (sum(squares) / len(squares)) ** 0.5 if squares else 0.0

        solvent = {'explicit': 'explicit water (%s), removed afterwards' %
                   water_model.upper(), 'vacuum': 'vacuum'}.get(
            self.solvent, 'implicit %s' % self.solvent.upper())
        say('Minimization (%s, %s, %s%s)' % (
            info['label'], solvent, platform,
            ': %s' % device if device else ''))
        say('  time: %.1f s' % (time.time() - started))
        say('  restraints: %s, %d atoms, k = %g kJ/mol/nm^2' % (
            self.restrain, restrained, self.restraint_k)
            if restrained else '  restraints: none')
        if fixed:
            say('  energy: %.1f, hydrogens and built atoms relaxed %.1f, '
                'final %.1f kJ/mol (%.1f kcal/mol)' % (
                    e0, e_relaxed, e1, e1 / mm.KJ_PER_KCAL))
        else:
            say('  energy: %.1f -> %.1f kJ/mol (%.1f -> %.1f kcal/mol)' % (
                e0, e1, e0 / mm.KJ_PER_KCAL, e1 / mm.KJ_PER_KCAL))
        say('  RMSD of input heavy atoms: %.2f A' % rmsd)
        if site_atoms:
            say('  metal sites held fixed: %d atoms (metals and the atoms '
                'they bind)' % len(site_atoms))
        if built:
            problems = self._build_problems(self.topology, self.positions)
            say('  built atoms: %d; %s' % (
                len([a for a in built if a.element is not None and
                     a.element.symbol != 'H']),
                'chirality and peptide bonds good' if not problems else
                'still with problems: ' + ', '.join(problems)))
        say('  net charge (force field): %+.2f' % charge)
        if excluded:
            say('  left out, no force field parameters (kept unchanged): '
                '%s' % ', '.join(sorted({r.name for r in excluded})))
        self.report.update(energy_before=e0, energy_after=e1, rmsd=rmsd,
                           platform=platform, device=device,
                           restrained_atoms=restrained,
                           charge=charge,
                           excluded_residues=sorted({r.name
                                                     for r in excluded}))

    # PyMOL part

    def apply(self, _self=cmd):
        '''
        Load the result as a new object; the input objects are hidden.
        '''
        if self.topology is None:
            raise RuntimeError('run() first')
        _self.delete(self.name)
        mm._load_topology(self.topology, self.positions, self.name,
                          self.b_factors, _self, self.chain_map)
        mm._restore_bonds(self.name, self.het_bonds, _self)
        # metal coordination: zero-order bonds (like Maestro), and dashes to
        # see them (PyMOL hardly draws zero-order bonds)
        mm._restore_bonds(self.name, self.metal_bonds, _self)
        if self.metal_bonds:
            self._metal_dashes(_self)
        _self.disable(' '.join(o for o in self.objects if o != self.name))
        self.report['object'] = self.name
        self.report['atoms'] = _self.count_atoms('%' + self.name)
        self._say('New object "%s": %d atoms' % (self.name,
                                                 self.report['atoms']))
        return self.report


def _select_atom(name, key):
    chain, resi, atom = key
    return '%%%s and chain "%s" and resi \\%s and name "%s"' % (
        name, chain, resi, atom)


def prepare_job(selection='all', name='', state=-1, _self=cmd, **options):
    '''
    :param options: PrepJob attributes
    :return: PrepJob with its input read from PyMOL
    '''
    mm._require()
    job = PrepJob()
    for key, value in options.items():
        if not hasattr(job, key) or key.startswith('_'):
            raise TypeError('unknown option: %s' % key)
        setattr(job, key, value)
    job.check_options()
    job.state = mm._state(state, _self)
    job.objects = mm._objects(selection, _self)
    job.selection = selection
    job.name = name or _self.get_legal_name(job.objects[0] + '_prep')
    (job.input_text, job.chain_map, job.polymer_ids, job.b_factors,
     job.alternates) = _export(
         selection, job.state, _self,
         hydrogens=not (job.preprocess and job.delete_hydrogens),
         altloc=job.altloc)
    job.input_hydrogens = _self.count_atoms(
        '(%s) and hydro and present' % selection, state=job.state)
    job.het_bonds = mm._het_bonds(selection, job.state, _self)
    if job.mutations and job.preprocess:
        resolved = []
        chains = sorted(job.polymer_ids)
        for chain, resi, new in job.mutations:
            if not chain:
                if len(chains) != 1:
                    raise ValueError('mutation %s:%s: give the chain (%s)' % (
                        resi, new, ', '.join(map(_chain_label, chains))))
                chain = chains[0]
            names = []
            _self.iterate_state(
                job.state, '(%s) and polymer and chain "%s" and resi \\%s '
                'and name CA' % (selection, chain, resi),
                'names.append(resn)', space={'names': names})
            if not names:
                raise ValueError('mutation %s:%s:%s: no such residue in "%s"'
                                 % (chain, resi, new, selection))
            resolved.append((chain, resi, names[0], new))
        job.mutations = resolved
    if job.preprocess and job.fill_loops and not job.source:
        job.residues = _unmodeled(selection, _self)
    return job


def _protonation_state(residue, disulfide=False):
    '''
    :return: variant name (ASH, HID, LYN, ...) from the residue's hydrogens
    '''
    names = {a.name for a in residue.atoms()}
    name = residue.name
    if name == 'ASP':
        return 'ASH' if names & {'HD1', 'HD2'} else 'ASP'
    if name == 'GLU':
        return 'GLH' if names & {'HE1', 'HE2'} else 'GLU'
    if name == 'HIS':
        d, e = 'HD1' in names, 'HE2' in names
        return 'HIP' if d and e else 'HID' if d else 'HIE' if e else 'HIN'
    if name == 'LYS':
        return 'LYS' if 'HZ3' in names else 'LYN'
    if name == 'CYS':
        return 'CYS' if 'HG' in names else 'CYX' if disulfide else 'CYM'
    return name


def openmm_vec(v):
    from openmm import Vec3
    return Vec3(float(v[0]), float(v[1]), float(v[2]))


def _backbone(topology):
    '''
    :return: [(residue, {'N', 'CA', 'C' (and 'CB', 'O', 'H', 'HA'):
        atom})] of residues with a backbone
    '''
    result = []
    for residue in topology.residues():
        named = {a.name: a for a in residue.atoms()
                 if a.name in ('N', 'CA', 'C', 'CB', 'O', 'H', 'HA')}
        if {'N', 'CA', 'C'} <= set(named):
            result.append((residue, named))
    return result


def _peptide_bonds(topology):
    '''
    :return: [(atoms of residue i, atoms of residue i+1)] bonded C-N
    '''
    named = {residue: atoms for residue, atoms in _backbone(topology)}
    pairs = []
    for bond in topology.bonds():
        a, b = bond.atom1, bond.atom2
        if a.name == 'N' and b.name == 'C':
            a, b = b, a
        if a.name == 'C' and b.name == 'N' and a.residue in named and \
                b.residue in named and a.residue is not b.residue:
            pairs.append((named[a.residue], named[b.residue]))
    return pairs


def parse_mutations(text):
    '''
    "A:45:GLY, A:102:W" (chain:residue:new; the chain can be left out for
    a structure with one chain; one- or three-letter codes)

    :return: [(chain or '', resi, new name)]
    '''
    import re
    found = []
    for token in re.split(r'[,;\s]+', str(text).strip()):
        if not token:
            continue
        parts = re.split(r'[:/]', token)
        if len(parts) == 2 and not re.fullmatch(r'-?\d+', parts[1]):
            parts = [''] + parts
        if len(parts) != 3:
            raise ValueError('mutation "%s": write chain:residue:new, e.g. '
                             'A:45:GLY' % token)
        chain, resi, new = parts
        new = ONE_LETTER.get(new.upper(), new.upper())
        if new not in AMINO_ACIDS:
            raise ValueError('mutation "%s": %s is not a standard amino acid'
                             % (token, new))
        if not re.fullmatch(r'-?\d+', resi):
            raise ValueError('mutation "%s": residue numbers only (no '
                             'insertion codes)' % token)
        found.append((chain, resi, new))
    return found


def _bool(value):
    if isinstance(value, str):
        return value.strip().lower() not in ('0', 'off', 'no', 'false', '')
    return bool(value)


def prepare_structure(selection='all', name='', preprocess=1,
                      altloc='occupancy', delete_hydrogens=1, cap_termini=0,
                      fill_loops=1, fill_termini=0, source='',
                      fill_side_chains=1, replace_nonstandard=1,
                      delete_waters=0, water_distance=0.0,
                      water_min_hbonds=0, delete_heterogens=0, mutations='',
                      disulfides=1, metals=1, protonate=1,
                      ph=7.0, pka_method='propka', ligand_hydrogens=1,
                      optimize_hbonds=1, minimize=1,
                      forcefield='amber14', solvent='gbn2', restrain='heavy',
                      restraint_k=1000.0, tolerance=10.0, max_iterations=0,
                      platform='auto', diagnostics=1, state=-1, quiet=0, *,
                      _self=cmd):
    '''
DESCRIPTION

    "prepare_structure" prepares a protein structure in a new object, like
    Maestro's Protein Preparation Wizard: preprocess (PDBFixer), protonate
    (PROPKA pKa values) and restrained minimization (OpenMM), then
    diagnostics.

USAGE

    prepare_structure [ selection [, name [, ... ]]]

ARGUMENTS

    selection = str: atoms to prepare {default: all}

    name = str: new object {default: <object>_prep}

    preprocess = 0/1: the preprocessing step {default: 1}

    altloc = occupancy or a letter: which alternate location to keep
    {default: occupancy (the highest)}

    delete_hydrogens = 0/1: delete existing hydrogens first (else they
    are kept and only missing ones added) {default: 1}

    cap_termini = 0/1: ACE/NME caps on protein chain ends {default: 0}

    fill_loops = 0/1: build missing residues (known from the file, or
    from "source") {default: 1}

    fill_termini = 0/1: also missing residues at chain ends {default: 0}

    source = str: sequence file or PDB ID for missing residues
    {default: the structure's own record}

    fill_side_chains = 0/1: build missing side-chain atoms {default: 1}

    replace_nonstandard = 0/1: e.g. MSE -> MET {default: 1}

    delete_waters, delete_heterogens = 0/1: delete waters, delete other
    non-polymer residues (ligands, ions) {default: 0}

    water_distance = float: with delete_waters, only waters farther than
    this from het groups (A) {default: 0 = all waters}

    water_min_hbonds = int: delete waters with fewer hydrogen bonds to
    non-waters {default: 0 = off}

    mutations = str: e.g. "A:45:GLY A:102:W" (chain:residue:new; one- or
    three-letter codes) {default: none}

    disulfides = 0/1: keep the disulfide bonds OpenMM recognizes (else
    free cysteines) {default: 1}

    metals = 0/1: metal sites: coordinating residues' protonation states
    (e.g. thiolate Cys, His bound by NE2: HID), coordination held during
    the minimization, zero-order bonds in the result {default: 1}

    protonate = 0/1: the protonation step {default: 1}

    ph = float {default: 7.0}

    pka_method = propka or standard {default: propka}

    ligand_hydrogens = 0/1: hydrogens for ligands from the PDB's Chemical
    Component Dictionary (needs internet) {default: 1}

    optimize_hbonds = 0/1: flip Asn/Gln/His, choose His tautomers, orient
    OH, SH, NH3+ and water hydrogens for hydrogen bonds {default: 1}

    minimize = 0/1: the minimization step {default: 1}

    forcefield = amber14, amber19, amber99sbildn, charmm36_2024 or
    charmm36 {default: amber14}

    solvent = gbn2, obc2, hct (implicit, AMBER only), explicit (water
    box, removed afterwards) or vacuum {default: gbn2}

    restrain = heavy, backbone, ca or none {default: heavy}

    restraint_k = float: kJ/mol/nm^2 {default: 1000}

    tolerance = float: convergence, kJ/mol/nm {default: 10}

    max_iterations = int: 0 = until converged {default: 0}

    platform = auto, CUDA, OpenCL, CPU or Reference {default: auto}

    diagnostics = 0/1: report problems of the result {default: 1}

EXAMPLE

    fetch 1oky
    prepare_structure 1oky, cap_termini=1, ph=6.5

SEE ALSO

    structure_diagnostics, fix_structure, minimize
    '''
    options = dict(
        preprocess=_bool(preprocess),
        altloc=str(altloc).strip(),
        delete_hydrogens=_bool(delete_hydrogens),
        cap_termini=_bool(cap_termini),
        fill_loops=_bool(fill_loops), fill_termini=_bool(fill_termini),
        source=str(source).strip(), fill_side_chains=_bool(fill_side_chains),
        replace_nonstandard=_bool(replace_nonstandard),
        delete_waters=_bool(delete_waters),
        water_distance=float(water_distance),
        water_min_hbonds=int(water_min_hbonds),
        delete_heterogens=_bool(delete_heterogens),
        mutations=parse_mutations(mutations),
        disulfides=_bool(disulfides), metals=_bool(metals),
        protonate=_bool(protonate), ph=float(ph),
        pka_method=str(pka_method).strip(),
        ligand_hydrogens=_bool(ligand_hydrogens),
        optimize_hbonds=_bool(optimize_hbonds), minimize=_bool(minimize),
        forcefield=str(forcefield).strip(), solvent=str(solvent).strip(),
        restrain=str(restrain).strip(), restraint_k=float(restraint_k),
        tolerance=float(tolerance), max_iterations=int(max_iterations),
        platform=str(platform).strip(), diagnostics=_bool(diagnostics))
    job = prepare_job(selection, name, state, _self=_self, **options)
    job.run()
    report = job.apply(_self)
    if job.diagnostics:
        from . import diagnostics
        report['diagnostics'] = diagnostics.run(
            '%' + job.name, _self=_self)
        job.log.append(diagnostics.format_report(report['diagnostics']))
    if not int(quiet):
        print('\n'.join(' ' + line for line in job.log))
    return report


def extend(_self):
    _self.extend('prepare_structure', prepare_structure)
    _self.prepare_structure = prepare_structure
