'''
Residue labels (pymolx.labels) and the toolbar's Labels window
(pymolx.gui.label_dialog).
'''

import os

import pytest

from pymol import cmd
from pymolx import labels, undo

DATA = os.path.join(os.path.dirname(__file__), '..', '..', 'data')


@pytest.fixture(autouse=True)
def clean_session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


def shown(selection='all'):
    '''{(model, resi, name): label} of the atoms showing one'''
    found = {}
    cmd.iterate('(%s) and rep labels' % selection,
                '_ = found.update({(model, resi, name): label}) '
                'if label else None', space={'found': found})
    return found


def test_residue_text():
    assert labels.residue_text('LYS', '593', 'K') == 'K593'
    assert labels.residue_text('LYS', '593', 'K', 'three') == 'LYS593'
    assert labels.residue_text('LYS', '593', 'K', 'none') == '593'
    assert labels.residue_text('GLY', '28A', 'G') == 'G28A'
    # no one-letter code: the name
    assert labels.residue_text('NAG', '601', 'X') == 'NAG601'
    assert labels.residue_text('ZN', '301', '') == 'ZN301'
    # ... set apart from the number where it ends in a digit
    assert labels.residue_text('SO4', '502', 'X') == 'SO4 502'
    assert labels.residue_text('SO4', '502', 'X', 'three') == 'SO4 502'
    # PyMOL's letter for water is not an amino acid's
    assert labels.residue_text('HOH', '812', 'O') == 'HOH812'
    assert labels.residue_text('DA', '5', 'X') == 'A5'
    assert labels.residue_text('U', '7', 'X') == 'U7'
    assert labels.residue_text('DA', '5', 'X', 'three') == 'DA5'


def test_one_label_per_residue():
    cmd.fab('AKG', 'pep')
    assert cmd.label_residues('pep') == 3
    assert shown() == {('pep', '1', 'CA'): 'A1', ('pep', '2', 'CA'): 'K2',
                       ('pep', '3', 'CA'): 'G3'}
    # another format and another atom: the labels move, none is left
    assert cmd.label_residues('pep', 'three', 'CB') == 3
    assert shown() == {('pep', '1', 'CB'): 'ALA1', ('pep', '2', 'CB'): 'LYS2',
                       # glycine has no beta carbon
                       ('pep', '3', 'CA'): 'GLY3'}
    assert cmd.label_residues('pep', 'none', 'NZ') == 3
    assert shown() == {('pep', '1', 'CA'): '1', ('pep', '2', 'NZ'): '2',
                       ('pep', '3', 'CA'): '3'}
    with pytest.raises(ValueError):
        cmd.label_residues('pep', names='two')


def test_whole_residues_of_the_selection():
    cmd.fab('AKG', 'pep')
    cmd.fab('WW', 'other')
    # one atom picks its residue, and no other
    assert cmd.label_residues('pep and resi 2 and name NZ') == 1
    assert shown() == {('pep', '2', 'CA'): 'K2'}
    # a label the user wrote on another residue stays
    cmd.label('other and resi 1 and name CZ2', '"ring"')
    cmd.label('pep and resi 2 and name NZ', '"amine"')
    assert cmd.label_residues('pep and resi 2+3') == 2
    assert shown() == {('pep', '2', 'CA'): 'K2', ('pep', '3', 'CA'): 'G3',
                       ('other', '1', 'CZ2'): 'ring'}
    assert cmd.label_residues('none') == 0

    assert labels.clear('pep and resi 2') == 1
    assert shown() == {('pep', '3', 'CA'): 'G3', ('other', '1', 'CZ2'): 'ring'}
    assert labels.clear('pep and resi 1+2') == 0


def test_font_size():
    cmd.fab('AK', 'pep')
    before = cmd.get_setting_float('label_size')
    cmd.label_residues('pep')
    assert cmd.get_setting_float('label_size') == before
    cmd.label_residues('pep', size=20)
    assert cmd.get_setting_float('label_size') == 20


def test_residues_without_an_alpha_carbon():
    cmd.fnab('AT', 'dna', mode='DNA', form='B', dbl_helix=0)
    cmd.pseudoatom('ion', name='ZN', resn='ZN', resi='301', chain='A',
                   segi='', elem='Zn', pos=(0, 0, 0))
    cmd.pseudoatom('wat', name='O', resn='HOH', resi='812', chain='A',
                   segi='', elem='O', pos=(3, 0, 0))
    cmd.label_residues('all', atom='CB')
    found = shown()
    # a nucleotide: its letter, on the sugar's C1'
    assert {key: text for key, text in found.items() if key[0] == 'dna'} \
        == {('dna', '1', "C1'"): 'A1', ('dna', '2', "C1'"): 'T2'}
    assert found[('ion', '301', 'ZN')] == 'ZN301'
    assert found[('wat', '812', 'O')] == 'HOH812'


def test_a_real_structure():
    # residues absent from the model are atoms without coordinates
    # (pymolx.missing), which could not show a label
    cmd.load(os.path.join(DATA, '1oky.pdb.gz'))
    absent = cmd.count_atoms('1oky and not present')
    assert absent == 27
    residues = set()
    cmd.iterate('1oky and present', 'residues.add((chain, resi))',
                space={'residues': residues})
    assert cmd.label_residues('1oky') == len(residues)
    labeled = []
    cmd.iterate('1oky and not present', 'labeled.append(label)',
                space={'labeled': labeled})
    assert labeled == [''] * absent
    found = shown('1oky')
    assert len(found) == len(residues)
    assert found[('1oky', '72', 'CA')] == 'P72'
    assert found[('1oky', '1368', 'C1')] == 'STU1368'     # the ligand
    assert found[('1oky', '1363', 'O1')] == 'SO4 1363'


def test_one_undo_step():
    cmd.fab('AKG', 'pep')
    cmd.undo_enable()
    try:
        steps = len(undo.stack.undo_stack)
        cmd.label_residues('pep', size=20)
        assert len(undo.stack.undo_stack) == steps + 1
        assert len(shown()) == 3
        cmd.undo()
        assert shown() == {}
        assert cmd.get_setting_float('label_size') != 20
    finally:
        cmd.undo_disable()


# the window

try:
    from pymol.Qt import QtCore, QtWidgets
    import shiboken6
except ImportError:
    QtWidgets = None

needs_qt = pytest.mark.skipif(QtWidgets is None, reason='no Qt')


@pytest.fixture
def dialog():
    from pymolx.gui import label_dialog
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    d = label_dialog.LabelDialog(_self=cmd)
    yield d
    shiboken6.delete(d)


def choices(d):
    return [(d.target.itemText(i), d.target.itemData(i))
            for i in range(d.target.count())]


@needs_qt
def test_dialog_targets(dialog):
    d = dialog
    # nothing loaded: everything, which is nothing
    assert choices(d) == [('All residues', 'all')]
    d.apply()
    assert d.status.text() == 'No residues to label in All residues.'

    cmd.fab('AKG', 'pep')
    cmd.select('lys', 'pep and resi 2', enable=0)
    cmd.pseudoatom('_hidden')
    d.refresh()
    assert choices(d) == [('(lys)', 'lys'), ('pep', 'pep'),
                          ('All residues', 'all')]
    assert d.selection() == 'all'
    # what was clicked in the viewer comes first, once there is some
    cmd.select('sele', 'pep and resi 3')
    d.refresh()
    assert d.selection() == 'sele'
    # ... unless the user asked for something else, while that exists
    row = d.target.findData('lys')
    d.target.setCurrentIndex(row)
    d.target.activated.emit(row)
    cmd.select('later', 'pep and resi 1')
    d.refresh()
    assert d.selection() == 'lys'
    cmd.delete('lys')
    d.refresh()
    assert d.selection() == 'sele'
    # the list is read again whenever it opens or the window is back
    cmd.delete('sele')
    d.target.about_to_open.emit()
    assert d.selection() == 'all'
    cmd.select('sele', 'pep and resi 3')
    QtWidgets.QApplication.sendEvent(
        d, QtCore.QEvent(QtCore.QEvent.Type.WindowActivate))
    assert d.selection() == 'sele'


@needs_qt
def test_dialog_labels(dialog):
    d = dialog
    cmd.fab('AKG', 'pep')
    cmd.select('sele', 'pep and resi 1+2')
    d.refresh()
    assert [d.names.itemText(i) for i in range(d.names.count())] == [
        '1-letter code (K593)', '3-letter code (LYS593)',
        'Number only (593)']
    assert [d.atom.itemText(i) for i in range(d.atom.count())] == \
        list(labels.ATOMS)
    assert d.atom.isEditable()

    d.label_button.click()
    assert shown() == {('pep', '1', 'CA'): 'A1', ('pep', '2', 'CA'): 'K2'}
    assert d.status.text() == 'Labeled 2 residues of (sele).'

    d.names.setCurrentIndex(d.names.findData('three'))
    d.atom.setCurrentText('NZ')         # any atom name can be typed
    d.label_button.click()
    assert shown() == {('pep', '1', 'CA'): 'ALA1', ('pep', '2', 'NZ'): 'LYS2'}

    # the font size is the viewer's, there and back
    assert d.size.value() == round(cmd.get_setting_float('label_size'))
    d.size.setValue(22)
    assert cmd.get_setting_float('label_size') == 22
    cmd.set('label_size', 9)
    d.refresh()
    assert d.size.value() == 9
    cmd.set('label_size', -0.8)         # in Angstrom: shown as the default
    d.refresh()
    assert d.size.value() == 14
    assert cmd.get_setting_float('label_size') == pytest.approx(-0.8)

    cmd.select('sele', 'pep and resi 2')
    d.clear_button.click()
    assert shown() == {('pep', '1', 'CA'): 'ALA1'}
    assert d.status.text() == 'Took the labels off 1 residue of (sele).'
    d.clear_button.click()
    assert d.status.text() == 'No labels in (sele).'
