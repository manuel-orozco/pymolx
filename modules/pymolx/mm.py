'''
Structure preparation (PDBFixer) and energy minimization (OpenMM).

Both are optional dependencies (pip install openmm pdbfixer); they are
imported on first use.

- fix_structure: PDBFixer fixes into a new object (missing atoms,
  nonstandard residues, heterogens, hydrogens; missing residues when the
  original file or PDB ID is given, as PyMOL doesn't keep SEQRES)
- minimize: OpenMM energy minimization of a selection, updating the
  coordinates in place. Atoms outside the selection (in the same objects)
  are held fixed. Hydrogens are added internally; only existing atoms
  are written back.

Minimization is split into prepare (reads PyMOL), MinimizationJob.run
(OpenMM only, may run in a worker thread) and MinimizationJob.apply
(writes PyMOL), so a GUI can keep its event loop running.
'''

import io
import os
import threading

from pymol import cmd

KJ_PER_KCAL = 4.184

# residue names OpenMM's Amber force fields have templates for
STANDARD_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT',
}
WATER = {'HOH', 'WAT', 'H2O', 'TIP3', 'SOL'}

FORCEFIELDS = {
    'amber14': ['amber14-all.xml'],
    'amber99sbildn': ['amber99sbildn.xml'],
    'charmm36': ['charmm36.xml'],
}
SOLVENT_MODELS = {
    # implicit solvent files per force field family
    'gbn2': {'amber14': 'implicit/gbn2.xml', 'amber99sbildn': 'implicit/gbn2.xml'},
    'obc2': {'amber14': 'implicit/obc2.xml', 'amber99sbildn': 'implicit/obc2.xml'},
    'vacuum': {},
}
RESTRAINTS = ('none', 'backbone', 'heavy')
BACKBONE = {'N', 'CA', 'C', 'O'}
PLATFORMS = ('CUDA', 'OpenCL', 'CPU', 'Reference')


class MissingDependency(Exception):
    pass


def _require():
    '''
    :return: (openmm, openmm.app, openmm.unit, PDBFixer)
    '''
    try:
        import openmm
        from openmm import app, unit
        from pdbfixer import PDBFixer
    except ImportError as e:
        raise MissingDependency(
            'Structure preparation needs OpenMM and PDBFixer: '
            'pip install openmm pdbfixer (missing: %s)' % e.name) from e
    return openmm, app, unit, PDBFixer


def available():
    try:
        _require()
        return True
    except MissingDependency:
        return False


def _chain_id(chain):
    # PDBFixer writes blank chain IDs as "."
    return '' if chain.id == '.' else chain.id


def _residue_id(residue):
    return residue.id + residue.insertionCode.strip()


def _key(chain, resi, name):
    return (chain, str(resi), name)


def _openmm_key(atom):
    residue = atom.residue
    return _key(_chain_id(residue.chain), _residue_id(residue), atom.name)


def _objects(selection, _self):
    from pymol import CmdException
    try:
        objects = _self.get_object_list('(%s)' % selection)
    except CmdException:
        objects = []
    if not objects:
        raise ValueError('no atoms in selection "%s"' % selection)
    return objects


def _state(state, _self):
    state = int(state)
    return _self.get_state() if state < 1 else state


def _cif_value(value):
    value = str(value)
    if value == '':
        return '.'
    if "'" in value or ' ' in value:
        return '"%s"' % value
    return value


def _export_cif(selection, state, _self):
    '''
    mmCIF text for OpenMM/PDBFixer. Written here rather than with
    cmd.get_str: label_asym_id is set to the chain (author chain ID),
    otherwise OpenMM may take the label IDs (one per molecule) as chains,
    and auth_seq_id carries the residue number.

    :raises ValueError: if (chain, residue, atom name) isn't unique, as
        coordinates are matched back by these
    '''
    rows = []
    _self.iterate_state(
        state, '(%s)' % selection,
        'rows.append((type, elem, name, alt, resn, chain, resv, resi, '
        'x, y, z, q, b))', space={'rows': rows})

    seen = set()
    duplicates = []
    lines = ['data_pymolx', 'loop_'] + ['_atom_site.' + column for column in (
        'group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id',
        'label_comp_id', 'label_asym_id', 'label_seq_id',
        'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy',
        'B_iso_or_equiv', 'auth_seq_id', 'auth_asym_id',
        'pdbx_PDB_model_num')]
    for n, (record, elem, name, alt, resn, chain, resv, resi, x, y, z, q,
            b) in enumerate(rows, 1):
        key = (chain, resi, name, alt)
        if key in seen:
            duplicates.append(key)
        seen.add(key)
        inscode = resi[len(str(resv)):] if resi.startswith(str(resv)) else ''
        lines.append(' '.join(_cif_value(v) for v in (
            record, n, elem, name, alt, resn, chain,
            resv, inscode or '?', '%.3f' % x, '%.3f' % y, '%.3f' % z,
            '%.2f' % q, '%.2f' % b, resv, chain, 1)))
    if duplicates:
        raise ValueError(
            'atoms are not unique by chain, residue and name (e.g. %s); '
            'rename chains or residues first' % (
                ', '.join('%s/%s/%s' % k[:3] for k in duplicates[:3])))
    return '\n'.join(lines) + '\n'


def _load_topology(topology, positions, name, b_factors, _self):
    '''
    Create a PyMOL object from an OpenMM topology with its exact bonds
    (loading a file would make PyMOL guess bonds from distances, which
    goes wrong for newly built atoms that start out close together).
    '''
    from chempy import models, Atom, Bond
    from openmm import unit

    xyz = positions.value_in_unit(unit.angstrom)
    model = models.Indexed()
    for atom in topology.atoms():
        residue = atom.residue
        a = Atom()
        a.name = atom.name
        a.symbol = atom.element.symbol if atom.element is not None else 'X'
        a.resn = residue.name
        a.resi = _residue_id(residue)
        digits = residue.id.lstrip('-')
        a.resi_number = int(residue.id) if digits.isdigit() else 0
        a.chain = a.segi = _chain_id(residue.chain)
        a.hetatm = int(residue.name not in STANDARD_RESIDUES)
        a.coord = [float(c) for c in xyz[atom.index]]
        a.q = 1.0
        a.b = b_factors.get((a.chain, a.resi, a.name), 0.0)
        model.add_atom(a)
    for bond in topology.bonds():
        b = Bond()
        b.index = [bond.atom1.index, bond.atom2.index]
        order = getattr(bond, 'order', None)
        b.order = order if order in (1, 2, 3) else 1
        model.add_bond(b)
    _self.load_model(model, name, zoom=0)


def _sequences_from_source(source):
    '''
    Polymer sequences per author chain from a structure file or PDB ID,
    as PDBFixer Sequence objects.
    '''
    from pdbfixer.pdbfixer import Sequence
    import gemmi

    path = source
    if not os.path.exists(source):
        import tempfile
        import urllib.request
        url = 'https://files.rcsb.org/download/%s.cif' % source.strip().lower()
        with urllib.request.urlopen(url, timeout=60) as response:
            data = response.read()
        handle, path = tempfile.mkstemp(suffix='.cif')
        with os.fdopen(handle, 'wb') as out:
            out.write(data)
    try:
        structure = gemmi.read_structure(path)
    finally:
        if path != source:
            os.unlink(path)
    structure.setup_entities()
    sequences = []
    for chain in structure[0]:
        polymer = chain.get_polymer()
        entity = structure.get_entity_of(polymer) if len(polymer) else None
        if entity is not None and entity.full_sequence:
            sequences.append(Sequence(chain.name, [
                gemmi.Entity.first_mon(item) for item in entity.full_sequence]))
    return sequences


# PDBFixer


def fix_structure(selection='all', name='', add_missing_residues=0,
                  add_missing_atoms=1, replace_nonstandard=1,
                  heterogens='keep_water', add_hydrogens=1, ph=7.0,
                  source='', state=-1, quiet=0, *, _self=cmd):
    '''
DESCRIPTION

    Fix a structure with PDBFixer into a new object: add missing heavy
    atoms, replace nonstandard residues, remove heterogens, add
    hydrogens. Missing residues need the original file or PDB ID
    ("source"), because PyMOL doesn't keep the sequence (SEQRES).

USAGE

    fix_structure selection [, name [, add_missing_residues
        [, add_missing_atoms [, replace_nonstandard [, heterogens
        [, add_hydrogens [, ph [, source ]]]]]]]]

ARGUMENTS

    selection = str: atoms to fix {default: all}

    name = str: new object {default: <object>_fixed}

    heterogens = none, keep_water or all: which heterogens to remove
    {default: keep_water}

    ph = float: pH for protonation {default: 7.0}

    source = str: original structure file or PDB ID, for the sequence
    of missing residues {default: none}

SEE ALSO

    minimize
    '''
    _, app, _, PDBFixer = _require()
    state = _state(state, _self)
    objects = _objects(selection, _self)
    name = name or _self.get_legal_name(objects[0] + '_fixed')
    quiet = int(quiet)

    text = _export_cif(selection, state, _self)
    fixer = PDBFixer(pdbxfile=io.StringIO(text))
    report = {'object': name, 'atoms_before': _self.count_atoms(
        '(%s)' % selection, state=state)}

    # missing residues: sequence from the original file or PDB ID
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    if int(add_missing_residues):
        if not source:
            raise ValueError('adding missing residues needs the original '
                             'structure file or PDB ID ("source")')
        fixer.sequences = _sequences_from_source(source)
        fixer.findMissingResidues()
    report['missing_residues'] = sum(
        len(v) for v in fixer.missingResidues.values())

    report['nonstandard'] = []
    if int(replace_nonstandard):
        fixer.findNonstandardResidues()
        report['nonstandard'] = sorted(
            {'%s->%s' % (r.name, new) for r, new in fixer.nonstandardResidues})
        fixer.replaceNonstandardResidues()

    report['heterogens_removed'] = 0
    if heterogens in ('keep_water', 'all'):
        before = fixer.topology.getNumResidues()
        fixer.removeHeterogens(keepWater=(heterogens == 'keep_water'))
        report['heterogens_removed'] = before - fixer.topology.getNumResidues()
    elif heterogens != 'none':
        raise ValueError('heterogens must be none, keep_water or all')

    fixer.findMissingAtoms()
    report['missing_atoms'] = sum(
        len(v) for v in fixer.missingAtoms.values()) + sum(
        len(v) for v in fixer.missingTerminals.values())
    if not int(add_missing_atoms):
        fixer.missingAtoms = {}
        fixer.missingTerminals = {}
    fixer.addMissingAtoms()

    n_before_h = fixer.topology.getNumAtoms()
    if int(add_hydrogens):
        fixer.addMissingHydrogens(float(ph))
    report['hydrogens_added'] = fixer.topology.getNumAtoms() - n_before_h

    # B-factors (e.g. AlphaFold pLDDT) of atoms which were there
    b_factors = {}
    _self.iterate_state(state, '(%s)' % selection,
                        'b_factors[(chain, resi, name)] = b',
                        space={'b_factors': b_factors})
    _self.delete(name)
    _load_topology(fixer.topology, fixer.positions, name, b_factors, _self)
    _self.disable(' '.join(objects))

    report['atoms_after'] = _self.count_atoms('%' + name)
    if not quiet:
        print(format_fix_report(report))
    return report


def format_fix_report(report):
    lines = [' fix_structure: new object "%s"' % report['object'],
             '  missing residues added: %d' % report['missing_residues'],
             '  missing heavy atoms added: %d' % report['missing_atoms'],
             '  nonstandard residues replaced: %s' % (
                 ', '.join(report['nonstandard']) or 'none'),
             '  heterogen residues removed: %d' % report['heterogens_removed'],
             '  hydrogens added: %d' % report['hydrogens_added'],
             '  atoms: %d -> %d' % (report['atoms_before'],
                                     report['atoms_after'])]
    return '\n'.join(lines)


# OpenMM minimization


class MinimizationCancelled(Exception):
    pass


class MinimizationJob:
    '''
    One minimization: created by prepare_minimization (reads PyMOL), then
    run() (OpenMM only), then apply() (writes PyMOL).
    '''

    def __init__(self):
        self.objects = []
        self.selection = ''
        self.state = 1
        self.forcefield = 'amber14'
        self.solvent = 'gbn2'
        self.restrain = 'backbone'
        self.restraint_k = 1000.0   # kJ/mol/nm^2
        self.max_iterations = 0      # 0: until converged
        self.tolerance = 10.0        # kJ/mol/nm
        self.platform = 'auto'
        self.ph = 7.0
        self.input_text = ''
        self.mobile_keys = set()
        self.existing_keys = set()   # all heavy atoms of the objects
        self.excluded = []           # residue names left out of the system
        self.result = None           # {key: (x, y, z)} in Angstrom
        self.report = {}
        self.progress = None         # callable(iteration, energy kJ/mol)
        self._cancel = threading.Event()

    def cancel(self):
        self._cancel.set()

    # OpenMM part

    def _build(self, openmm, app, unit, PDBFixer):
        fixer = PDBFixer(pdbxfile=io.StringIO(self.input_text))

        # leave out residues the force field has no templates for
        excluded = [r for r in fixer.topology.residues()
                    if r.name not in STANDARD_RESIDUES]
        self.excluded = sorted({r.name for r in excluded})
        if excluded:
            modeller = app.Modeller(fixer.topology, fixer.positions)
            modeller.delete(excluded)
            fixer.topology, fixer.positions = (modeller.topology,
                                               modeller.positions)

        fixer.findMissingResidues()
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        self.report['atoms_added_internally'] = sum(
            len(v) for v in fixer.missingAtoms.values())
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(self.ph)

        files = list(FORCEFIELDS[self.forcefield])
        implicit = SOLVENT_MODELS[self.solvent].get(self.forcefield)
        if implicit:
            files.append(implicit)
        elif self.solvent != 'vacuum':
            raise ValueError('%s implicit solvent is not available with %s' % (
                self.solvent, self.forcefield))
        forcefield = app.ForceField(*files)

        n_atoms = fixer.topology.getNumAtoms()
        method = app.NoCutoff if n_atoms <= 5000 else app.CutoffNonPeriodic
        kwargs = {} if method is app.NoCutoff else {
            'nonbondedCutoff': 1.2 * unit.nanometer}
        system = forcefield.createSystem(fixer.topology,
                                         nonbondedMethod=method,
                                         constraints=None, **kwargs)
        return fixer.topology, fixer.positions, system

    def _context(self, openmm, system, positions):
        names = PLATFORMS if self.platform == 'auto' else (self.platform,)
        errors = []
        for name in names:
            try:
                platform = openmm.Platform.getPlatformByName(name)
                integrator = openmm.VerletIntegrator(0.001)
                context = openmm.Context(system, integrator, platform)
                context.setPositions(positions)
                return context, name
            except Exception as e:
                errors.append('%s: %s' % (name, e))
        raise RuntimeError('no OpenMM platform worked: ' + '; '.join(errors))

    def run(self):
        '''
        Minimize (OpenMM only, no PyMOL calls). Fills result and report.

        :raises MinimizationCancelled: if cancel() was called
        '''
        openmm, app, unit, PDBFixer = _require()
        topology, positions, system = self._build(openmm, app, unit, PDBFixer)

        # fixed atoms: not in the selection (hydrogens follow their parent)
        heavy_parent = {}
        for bond in topology.bonds():
            a, b = bond.atom1, bond.atom2
            if a.element is not None and a.element.symbol == 'H':
                heavy_parent[a.index] = b
            elif b.element is not None and b.element.symbol == 'H':
                heavy_parent[b.index] = a
        # atoms added internally (e.g. a missing OXT) move with their residue
        mobile_residues = {key[:2] for key in self.mobile_keys}
        mobile = []
        for atom in topology.atoms():
            key = _openmm_key(heavy_parent.get(atom.index, atom))
            is_mobile = key in self.mobile_keys or (
                key not in self.existing_keys and key[:2] in mobile_residues)
            mobile.append(is_mobile)
            if not is_mobile:
                system.setParticleMass(atom.index, 0.0)
        self.report['mobile_atoms'] = sum(mobile)
        self.report['fixed_atoms'] = len(mobile) - sum(mobile)
        if not any(mobile):
            raise ValueError('nothing to minimize: no mobile atoms')

        for force in system.getForces():
            force.setForceGroup(0)
        restrained = 0
        if self.restrain != 'none':
            restraint = openmm.CustomExternalForce(
                '0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
            restraint.addGlobalParameter('k', self.restraint_k)
            for name in ('x0', 'y0', 'z0'):
                restraint.addPerParticleParameter(name)
            for atom, is_mobile in zip(topology.atoms(), mobile):
                if not is_mobile or atom.element.symbol == 'H':
                    continue
                if self.restrain == 'backbone' and atom.name not in BACKBONE:
                    continue
                restraint.addParticle(
                    atom.index, positions[atom.index].value_in_unit(
                        unit.nanometer))
                restrained += 1
            restraint.setForceGroup(1)
            system.addForce(restraint)
        self.report['restrained_atoms'] = restrained

        context, platform = self._context(openmm, system, positions)
        self.report['platform'] = platform

        def energy():
            return context.getState(getEnergy=True, groups={0}) \
                .getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

        self.report['energy_before'] = energy()

        job = self

        class Reporter(openmm.MinimizationReporter):
            def report(self, iteration, x, grad, args):
                if job.progress is not None:
                    # args is a C++ map: no .get()
                    energy = (args['system energy']
                              if 'system energy' in args else 0.0)
                    job.progress(iteration, energy)
                return job._cancel.is_set()

        openmm.LocalEnergyMinimizer.minimize(
            context, self.tolerance * unit.kilojoule_per_mole / unit.nanometer,
            int(self.max_iterations), Reporter())
        if self._cancel.is_set():
            raise MinimizationCancelled()

        self.report['energy_after'] = energy()
        final = context.getState(getPositions=True).getPositions(asNumpy=True)
        final = final.value_in_unit(unit.angstrom)
        self.result = {_openmm_key(atom): tuple(map(float, final[atom.index]))
                       for atom, is_mobile in zip(topology.atoms(), mobile)
                       if is_mobile}
        self.report['excluded_residues'] = self.excluded
        return self.report

    # PyMOL part

    def apply(self, _self=cmd):
        '''
        Write minimized coordinates of existing atoms back (one undo step).
        '''
        if self.result is None:
            raise RuntimeError('run() first')
        new = self.result
        before = {}
        _self.iterate_state(self.state, '(%s) and not hydro' % self.selection,
                            'before[(chain, resi, name)] = (x, y, z)',
                            space={'before': before})
        _self.alter_state(
            self.state, '(%s)' % self.selection,
            '(x, y, z) = new.get((chain, resi, name), (x, y, z))',
            space={'new': new})
        # existing hydrogens which weren't matched follow their atoms
        if _self.count_atoms('(%s) and hydro' % self.selection):
            _self.h_fix('(%s) and hydro' % self.selection)

        moved = [(before[k], new[k]) for k in before if k in new]
        if moved:
            msd = sum((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2
                      for a, b in moved) / len(moved)
            self.report['rmsd'] = msd ** 0.5
        self.report['updated_atoms'] = len(moved)
        return self.report


def prepare_minimization(selection='all', *, forcefield='amber14',
                         solvent='gbn2', restrain='backbone',
                         restraint_k=1000.0, max_iterations=0,
                         tolerance=10.0, platform='auto', ph=7.0, state=-1,
                         _self=cmd):
    '''
    :return: MinimizationJob with the input read from PyMOL
    '''
    _require()
    if forcefield not in FORCEFIELDS:
        raise ValueError('forcefield must be one of: ' + ', '.join(FORCEFIELDS))
    if solvent not in SOLVENT_MODELS:
        raise ValueError('solvent must be one of: ' + ', '.join(SOLVENT_MODELS))
    if restrain not in RESTRAINTS:
        raise ValueError('restrain must be one of: ' + ', '.join(RESTRAINTS))

    job = MinimizationJob()
    job.state = _state(state, _self)
    job.objects = _objects(selection, _self)
    job.selection = selection
    job.forcefield = forcefield
    job.solvent = solvent
    job.restrain = restrain
    job.restraint_k = float(restraint_k)
    job.max_iterations = int(max_iterations)
    job.tolerance = float(tolerance)
    job.platform = platform
    job.ph = float(ph)

    system = '(%s) and not hydro' % ' or '.join(
        '%' + obj for obj in job.objects)
    job.input_text = _export_cif(system, job.state, _self)
    _self.iterate_state(job.state, system, 'keys.add((chain, resi, name))',
                        space={'keys': job.existing_keys})
    _self.iterate_state(job.state, '(%s) and not hydro' % selection,
                        'keys.add((chain, resi, name))',
                        space={'keys': job.mobile_keys})
    return job


def minimize(selection='all', forcefield='amber14', solvent='gbn2',
             restrain='backbone', restraint_k=1000.0, max_iterations=0,
             tolerance=10.0, platform='auto', ph=7.0, state=-1, quiet=0,
             *, _self=cmd):
    '''
DESCRIPTION

    Energy-minimize a selection with OpenMM and update its coordinates.
    Other atoms of the same objects are held fixed. Hydrogens are added
    internally; residues without force field parameters (ligands,
    ions, waters) are left out.

USAGE

    minimize [ selection [, forcefield [, solvent [, restrain
        [, restraint_k [, max_iterations [, tolerance [, platform ]]]]]]]]

ARGUMENTS

    selection = str: atoms to minimize {default: all}

    forcefield = amber14, amber99sbildn or charmm36 {default: amber14}

    solvent = gbn2, obc2 (implicit) or vacuum {default: gbn2}

    restrain = none, backbone or heavy: harmonic restraints to the
    starting positions {default: backbone}

    restraint_k = float: restraint force constant, kJ/mol/nm^2
    {default: 1000}

    max_iterations = int: 0 = until converged {default: 0}

    tolerance = float: convergence, kJ/mol/nm {default: 10}

    platform = auto, CUDA, OpenCL, CPU or Reference {default: auto}

SEE ALSO

    fix_structure, sculpt_activate
    '''
    job = prepare_minimization(
        selection, forcefield=forcefield, solvent=solvent, restrain=restrain,
        restraint_k=restraint_k, max_iterations=max_iterations,
        tolerance=tolerance, platform=platform, ph=ph, state=state,
        _self=_self)
    job.run()
    report = job.apply(_self)
    if not int(quiet):
        print(format_minimize_report(report))
    return report


def format_minimize_report(report):
    e0, e1 = report['energy_before'], report['energy_after']
    lines = [
        ' minimize: %d atoms moved, %d held fixed (%s)' % (
            report['mobile_atoms'], report['fixed_atoms'],
            report['platform']),
        '  energy: %.1f -> %.1f kJ/mol (%.1f -> %.1f kcal/mol)' % (
            e0, e1, e0 / KJ_PER_KCAL, e1 / KJ_PER_KCAL),
        '  RMSD of updated atoms: %.3f A' % report.get('rmsd', 0.0),
    ]
    if report.get('excluded_residues'):
        lines.append('  left out (no force field parameters): %s' %
                     ', '.join(report['excluded_residues']))
    return '\n'.join(lines)


def extend(_self):
    for func in (fix_structure, minimize):
        _self.extend(func.__name__, func)
        setattr(_self, func.__name__, func)
