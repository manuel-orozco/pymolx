'''
Structure preparation and minimization (pymolx.mm: fix_structure,
minimize) and the Structure Preparation dialog. Skipped without OpenMM
and PDBFixer. Uses the CPU platform, so no GPU is needed.
'''

import os

import pytest

from pymol import cmd
from pymolx import mm, undo

pytestmark = pytest.mark.skipif(not mm.available(),
                                reason='needs openmm and pdbfixer')

DATA = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
SEQUENCE = 'ACDEFGHIKLMNPQRSTVWY'


@pytest.fixture(autouse=True)
def session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


def peptide(name='pep'):
    cmd.fab(SEQUENCE, name, ss=1)
    cmd.remove('%s and hydro' % name)


def coords(selection):
    return cmd.get_coords(selection).copy()


# minimize

def test_minimize_lowers_energy():
    peptide()
    before = coords('pep')
    report = cmd.minimize('pep', platform='CPU', quiet=1)
    assert report['energy_after'] < report['energy_before']
    assert report['mobile_atoms'] > 0 and report['fixed_atoms'] == 0
    assert report['updated_atoms'] == cmd.count_atoms('pep')
    after = coords('pep')
    assert abs(after - before).max() > 0.01
    assert 0 < report['rmsd'] < 2.0  # backbone restrained
    assert cmd.count_atoms('pep and hydro') == 0  # no atoms added


def test_minimize_energy_drops_without_restraints():
    peptide()
    report = cmd.minimize('pep', restrain='none', solvent='vacuum',
                          platform='CPU', quiet=1)
    assert report['energy_after'] < report['energy_before']
    assert report['restrained_atoms'] == 0


def test_unselected_atoms_stay_fixed():
    peptide()
    fixed_before = coords('pep and resi 1-10')
    report = cmd.minimize('pep and resi 11-20', platform='CPU', quiet=1)
    assert report['fixed_atoms'] > 0
    assert abs(coords('pep and resi 1-10') - fixed_before).max() == 0.0
    assert report['energy_after'] < report['energy_before']


def test_minimize_is_one_undo_step():
    peptide()
    before = coords('pep')
    cmd.undo_enable()
    try:
        cmd.minimize('pep', platform='CPU', max_iterations=50, quiet=1)
        assert len(undo.stack.undo_stack) == 1
        cmd.undo()
        assert abs(coords('pep') - before).max() < 1e-3
    finally:
        cmd.undo_disable()


def test_ligands_and_waters_are_left_out():
    cmd.load(os.path.join(DATA, '1rx1.pdb'), 'p')
    ligand_before = coords('p and resn NAP')
    report = cmd.minimize('p and resi 1-10', platform='CPU',
                          max_iterations=20, quiet=1)
    assert {'NAP', 'HOH'} <= set(report['excluded_residues'])
    assert abs(coords('p and resn NAP') - ligand_before).max() == 0.0


def test_cancel():
    peptide()
    job = mm.prepare_minimization('pep', platform='CPU')
    job.cancel()
    with pytest.raises(mm.MinimizationCancelled):
        job.run()


def test_progress_reports_energies():
    peptide()
    job = mm.prepare_minimization('pep', platform='CPU', max_iterations=30)
    seen = []
    job.progress = lambda iteration, energy: seen.append((iteration, energy))
    job.run()
    assert seen and all(isinstance(e, float) for _, e in seen)


def test_bad_options():
    peptide()
    for kwargs in ({'forcefield': 'nope'}, {'solvent': 'nope'},
                   {'restrain': 'nope'}):
        with pytest.raises(ValueError):
            mm.prepare_minimization('pep', **kwargs)
    with pytest.raises(ValueError):
        mm.prepare_minimization('nothing_here')


# fix_structure

def test_fix_structure_adds_atoms_and_hydrogens():
    peptide()
    cmd.alter('pep', 'b = 77.0')
    cmd.remove('pep and resi 5 and not name N+CA+C+O+CB')  # truncated PHE
    report = cmd.fix_structure('pep', quiet=1)
    assert report['object'] == 'pep_fixed'
    assert report['missing_atoms'] >= 6
    assert report['hydrogens_added'] > 0
    assert cmd.count_atoms('pep_fixed and resi 5 and not hydro') == 11
    assert cmd.get_chains('pep_fixed') == ['']
    # existing atoms keep their B-factors (e.g. AlphaFold pLDDT)
    assert cmd.count_atoms('pep_fixed and b=77') == cmd.count_atoms('pep')
    assert cmd.get_names(enabled_only=1) == ['pep_fixed']


def test_fix_structure_heterogens():
    cmd.load(os.path.join(DATA, '1rx1.pdb'), 'p')
    cmd.fix_structure('p', 'keep', heterogens='keep_water',
                      add_hydrogens=0, quiet=1)
    cmd.fix_structure('p', 'none', heterogens='all', add_hydrogens=0,
                      quiet=1)
    assert cmd.count_atoms('keep and solvent') == cmd.count_atoms(
        'p and solvent')
    assert cmd.count_atoms('keep and resn NAP') == 0
    assert cmd.count_atoms('none and solvent') == 0


def test_fix_missing_residues_needs_source():
    peptide()
    with pytest.raises(ValueError):
        cmd.fix_structure('pep', add_missing_residues=1, quiet=1)


def test_fix_then_minimize_keeps_names():
    peptide()
    cmd.fix_structure('pep', quiet=1)
    report = cmd.minimize('pep_fixed', platform='CPU', max_iterations=50,
                          quiet=1)
    # hydrogens from PDBFixer have OpenMM names: all atoms written back
    assert report['updated_atoms'] == cmd.count_atoms(
        'pep_fixed and not hydro')
    assert report['energy_after'] < report['energy_before']


# dialog

def test_dialog():
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import pymol.plugins  # noqa: F401 (adds the startup plugin path)
    from pmg_tk.startup.structure_prep_gui.main import StructurePrepDialog

    peptide()
    dialog = StructurePrepDialog(None, _self=cmd)
    try:
        dialog.refresh()
        dialog.selection.setEditText('pep')
        dialog.fix()
        assert 'pep_fixed' in cmd.get_names()
        assert dialog.selection.currentText() == 'pep_fixed'

        before = coords('pep_fixed')
        dialog.platform.setCurrentText('CPU')
        dialog.max_iterations.setValue(50)
        dialog.minimize()
        assert dialog.thread is not None
        assert not dialog.minimize_button.isEnabled()
        dialog.thread.wait()
        dialog._finished()  # normally delivered by the event loop
        assert dialog.progress.text() == 'Done'
        assert abs(coords('pep_fixed') - before).max() > 0.01
        assert 'energy:' in dialog.log.toPlainText()
        assert dialog.minimize_button.isEnabled()

        # CHARMM36 has no implicit solvent here: vacuum is forced
        dialog.forcefield.setCurrentText('charmm36')
        assert dialog.solvent.currentData() == 'vacuum'
    finally:
        dialog.deleteLater()


# regressions: chains renamed by OpenMM's chain-ID heuristic, and bonds
# guessed from distances for newly built atoms (stretched "starbursts")

def two_chains():
    '''
    Chains H and L whose mmCIF label IDs (segi) outnumber the chains, the
    case where OpenMM would otherwise use label IDs as chains.
    '''
    cmd.fab('EVQLVESGGG', 'h', ss=1)
    cmd.fab('DIQMTQSPSS', 'l', ss=1)
    cmd.alter('h', 'chain="H"; segi="B"')
    cmd.alter('l', 'chain="L"; segi="C"; resv += 100')
    cmd.create('ab', 'h or l')
    cmd.delete('h l')
    cmd.alter('ab and chain L', 'segi="D"', space={})
    cmd.alter('ab and chain L and resi 106-110', 'segi="E"')
    cmd.remove('ab and hydro')
    cmd.translate([0, 15, 0], 'ab and chain L')


def bond_problems(obj):
    import math
    from collections import Counter
    model = cmd.get_model(obj)
    stretched = sum(
        1 for b in model.bond
        if math.dist(model.atom[b.index[0]].coord,
                     model.atom[b.index[1]].coord) > 2.2)
    degree = Counter(i for b in model.bond for i in b.index)
    crowded = sum(1 for n in degree.values() if n > 4)
    return stretched, crowded


def test_chains_kept_with_many_label_ids():
    two_chains()
    before = coords('ab and chain L')
    report = cmd.minimize('ab', platform='CPU', max_iterations=100, quiet=1)
    assert report['updated_atoms'] == cmd.count_atoms('ab')
    assert abs(coords('ab and chain L') - before).max() > 0.01

    cmd.fix_structure('ab', quiet=1)
    assert cmd.get_chains('ab_fixed') == ['H', 'L']


def test_fixed_structure_has_exact_bonds():
    two_chains()
    cmd.remove('ab and chain H and resi 3 and not name N+CA+C+O+CB')
    cmd.fix_structure('ab', quiet=1)
    assert bond_problems('ab_fixed') == (0, 0)
    # one bond per pair; peptide bonds across residues are there
    assert cmd.count_atoms('ab_fixed and chain H and resi 1 and name C '
                           'and bound_to (resi 2 and name N)') == 1
    cmd.minimize('ab_fixed', platform='CPU', max_iterations=200, quiet=1)
    assert bond_problems('ab_fixed') == (0, 0)


def test_duplicate_atoms_are_refused():
    peptide()
    cmd.copy('pep2', 'pep')
    cmd.create('both', 'pep or pep2')  # same chain, residues and names
    with pytest.raises(ValueError, match='not unique'):
        mm.prepare_minimization('both', platform='CPU')
