'''
Binding pocket (pymolx.pocket: binding_pocket) and the Binding Pocket
dialog.
'''

import csv
import os

import numpy as np
import pytest

from pymol import cmd
from pymolx import pocket, undo

DATA = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
# PDK1 with staurosporine (STU 1368): hinge S160 O and A162 N bind it
STU = '1oky and resn STU'


@pytest.fixture(autouse=True)
def session():
    cmd.reinitialize()
    cmd.load(os.path.join(DATA, '1oky.pdb.gz'), '1oky')
    yield
    cmd.reinitialize()


def brute_force(ligand, cutoff, protein='1oky and polymer'):
    '''
    {(chain, resi): closest distance} over heavy atoms
    '''
    lig = cmd.get_coords('(%s) and not hydro' % ligand)
    atoms = []
    cmd.iterate_state(1, '(%s) and not hydro' % protein,
                      'atoms.append((chain, resi, x, y, z))',
                      space={'atoms': atoms})
    found = {}
    for chain, resi, *xyz in atoms:
        d = float(np.linalg.norm(lig - np.array(xyz), axis=1).min())
        if d <= cutoff:
            key = (chain, resi)
            found[key] = min(d, found.get(key, d))
    return found


def test_ligands():
    found = pocket.ligands('1oky')
    names = [x['resn'] for x in found]
    assert set(names) == {'STU', 'GOL', 'SO4'}  # no water, no SEP
    stu = found[names.index('STU')]
    assert stu['atoms'] == 35 and stu['resi'] == '1368'
    assert cmd.count_atoms(stu['selection']) == 35


def test_pocket_residues_and_distances():
    result = pocket.find(STU, 4.0)
    expected = brute_force(STU, 4.0)
    assert {(r['chain'], r['resi']) for r in result.residues} == \
        set(expected)
    assert len(result.residues) == 16
    for r in result.residues:
        assert r['distance'] == pytest.approx(
            expected[(r['chain'], r['resi'])], abs=1e-3)
        assert r['distance'] <= 4.0
    by_resi = {r['resi']: r for r in result.residues}
    for resi in ('160', '162', '166'):
        assert by_resi[resi]['polar']
    assert by_resi['160']['atom'] == 'O' and by_resi['162']['atom'] == 'N'
    assert not by_resi['96']['polar']  # Val CG2 to carbon
    assert by_resi['160']['ligand_atom'] == 'STU1368 N1'
    # in sequence order
    numbers = [int(r['resi']) for r in result.residues]
    assert numbers == sorted(numbers)


def test_selection_and_formats():
    result = pocket.find(STU, 4.0, name='site')
    assert result.selection == 'site'
    residues = set()
    cmd.iterate('site', 'residues.add(resi)', space={'residues': residues})
    assert residues == {r['resi'] for r in result.residues}
    # whole residues
    assert cmd.count_atoms('site') == cmd.count_atoms('byres site')
    assert result.text('one').startswith('L88, G89, V96, A109, K111,')
    assert result.text('three').startswith('LEU88, GLY89, VAL96,')
    assert result.text('none').startswith('88, 89, 96, 109,')
    assert result.text('one').endswith('T222, D223')


def test_cutoff():
    small = pocket.find(STU, 3.0, name='')
    large = pocket.find(STU, 6.0, name='')
    key = lambda result: {(r['chain'], r['resi']) for r in result.residues}
    assert key(small) < key(pocket.find(STU, 4.0, name='')) < key(large)
    assert key(large) == set(brute_force(STU, 6.0))
    with pytest.raises(ValueError, match='positive'):
        pocket.find(STU, 0)


def test_waters_and_hydrogens():
    plain = pocket.find(STU, 4.0, name='')
    wet = pocket.find(STU, 4.0, waters=True, name='')
    waters = [r for r in wet.residues if r['resn'] == 'HOH']
    assert waters and len(wet.residues) == len(plain.residues) + len(waters)
    # hydrogens added: the default (heavy atoms) doesn't change
    cmd.h_add('1oky')
    assert pocket.find(STU, 4.0, name='').text() == plain.text()
    with_h = pocket.find(STU, 4.0, hydrogens=True, name='')
    assert len(with_h.residues) > len(plain.residues)


def test_several_chains():
    cmd.alter('1oky and polymer and resi 150-400', 'chain="B"')
    cmd.alter(STU, 'chain="B"')
    text = pocket.find(STU, 4.0, name='').text('one')
    lines = text.split('\n')
    assert lines[0] == 'A: L88, G89, V96, A109, K111'
    assert lines[1].startswith('B: L159, S160, Y161,')


def test_ion_and_other_selections():
    cmd.load(os.path.join(DATA, '1rx1.pdb'), 'rx')
    result = pocket.find('rx and resn CA', 3.0, name='')
    assert result.residues
    assert all(r['model'] == 'rx' for r in result.residues)
    # a protein selection limits the pocket
    part = pocket.find(STU, 4.0, protein='1oky and resi 1-150', name='')
    assert [r['resi'] for r in part.residues] == ['88', '89', '96', '109',
                                                  '111']
    with pytest.raises(ValueError, match='no atoms'):
        pocket.find('1oky and resn XYZ')


def test_command(tmp_path, capsys):
    path = tmp_path / 'pocket.csv'
    text = cmd.binding_pocket(STU, 4.0, 'three', file=str(path))
    assert text.startswith('LEU88, GLY89')
    out = capsys.readouterr().out
    assert '16 residues within 4 A' in out and 'selection "pocket"' in out
    assert cmd.count_atoms('pocket')
    with open(path) as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 16 and rows[0]['resn'] == 'LEU'
    assert {row['resi'] for row in rows if row['polar'] == '1'} >= {
        '160', '162', '166'}
    with pytest.raises(ValueError, match='names'):
        cmd.binding_pocket(STU, names='four')


def test_one_undo_step():
    cmd.undo_enable()
    try:
        steps = len(undo.stack.undo_stack)
        cmd.binding_pocket(STU, quiet=1)
        assert len(undo.stack.undo_stack) == steps + 1
    finally:
        cmd.undo_disable()


def test_dialog():
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import pymol.plugins  # noqa: F401 (adds the startup plugin path)
    from pmg_tk.startup.pocket_gui.main import PocketDialog

    dialog = PocketDialog(None, _self=cmd)
    try:
        dialog.refresh()
        assert dialog.object.currentText() == '1oky'
        labels = [dialog.ligands.item(k).text()
                  for k in range(dialog.ligands.count())]
        assert labels[0].startswith('GOL 1360')
        assert [item.text() for item in dialog.ligands.selectedItems()] == [
            'STU 1368  (chain A, 35 atoms)']  # the largest
        dialog.run()
        assert dialog.text.toPlainText().startswith('L88, G89')
        assert dialog.table.rowCount() == 16
        dialog.names.setCurrentIndex(1)
        assert dialog.text.toPlainText().startswith('LEU88, GLY89')
        assert cmd.count_atoms('pocket')
        dialog.cutoff.setValue(3.0)
        dialog.custom.setText('resn STU')
        dialog.run()
        assert dialog.table.rowCount() < 16
        assert dialog.show_button.isEnabled()
        dialog.show_pocket()
        assert not dialog.grab().isNull()
    finally:
        dialog.deleteLater()
