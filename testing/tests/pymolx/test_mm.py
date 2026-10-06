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
