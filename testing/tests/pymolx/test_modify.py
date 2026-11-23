'''
Modify Sequences (pymolx.modify: renumber_residues, set_chain_id,
set_segment_id, swap_chains, merge_chains, combine_chains) and its
dialog.
'''

import os

import pytest

from pymol import cmd
from pymolx import modify, undo

DATA = os.path.join(os.path.dirname(__file__), '..', '..', 'data')


@pytest.fixture(autouse=True)
def session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


def peptide(name, sequence, chain='A', first=1, segi=''):
    cmd.fab(sequence, name, chain=chain, resi=first, segi=segi, hydro=0)
    return name


def residues(selection):
    '''
    [(segi, chain, resi, one-letter code)] in the order of the atoms
    '''
    found = []
    cmd.iterate('(%s) and name CA' % selection,
                'found.append((segi, chain, resi, oneletter))',
                space={'found': found})
    return found


def numbers(selection):
    return [r[2] for r in residues(selection)]


def count_residues(selection):
    '''
    Residues with any atom in the selection (1oky's last one has no CA)
    '''
    found = set()
    cmd.iterate(selection, 'found.add((segi, chain, resi))',
                space={'found': found})
    return len(found)


def names():
    '''
    Objects and selections, without the "_pk..." ones cmd.fab leaves
    '''
    return sorted(n for n in cmd.get_names('all') if not n.startswith('_pk'))


def letters(selection):
    return ''.join(r[3] for r in residues(selection))


def antibody():
    '''
    One object with a heavy chain numbered 1-6 (28A style insertion at
    2A) and a light chain numbered 1-5
    '''
    peptide('h', 'EVQLVE', 'H')
    cmd.alter('h and resi 3', 'resi = "2A"')
    cmd.sort('h')
    peptide('l', 'DIQMT', 'L')
    cmd.create('fab_1', 'h or l')
    cmd.delete('h')
    cmd.delete('l')
    assert numbers('fab_1 and chain H') == ['1', '2', '2A', '4', '5', '6']
    return 'fab_1'


def test_split_resi():
    assert modify.split_resi('28A') == (28, 'A')
    assert modify.split_resi('-3') == (-3, '')
    assert modify.split_resi('-3B') == (-3, 'B')


def test_resi_selection():
    assert modify.resi_selection('') == ''
    assert modify.resi_selection('  ') == ''
    assert modify.resi_selection('100-150') == 'resi 100-150'
    assert modify.resi_selection('100 - 150, 200A') == 'resi 100-150+200A'
    assert modify.resi_selection('5 to 9') == 'resi 5-9'
    assert modify.resi_selection('5:9+12') == 'resi 5-9+12'
    assert modify.resi_selection('-5-3') == 'resi \\-5-3'
    assert modify.resi_selection('-5--3') == 'resi \\-5-\\-3'
    with pytest.raises(ValueError, match='neither a number nor a range'):
        modify.resi_selection('chain A')
    # PyMOL reads them as meant
    peptide('m', 'ACDEFGHIK', first=-5)
    assert numbers('m and ' + modify.resi_selection('-5--3, 1 to 2')) == \
        ['-5', '-4', '-3', '1', '2']


def test_chain_table():
    antibody()
    cmd.load(os.path.join(DATA, '1oky.pdb.gz'), '1oky')
    table = modify.chain_table()
    assert [(r['object'], r['chain']) for r in table] == [
        ('fab_1', 'H'), ('fab_1', 'L'), ('1oky', 'A')]
    heavy, light, kinase = table
    assert heavy['residues'] == 6 and heavy['sequence'] == 'EVQLVE'
    assert (heavy['first'], heavy['last'], heavy['gaps']) == ('1', '6', 1)
    assert heavy['others'] == 0 and heavy['segi'] == ''
    assert (light['first'], light['last'], light['gaps']) == ('1', '5', 0)
    assert cmd.count_atoms(light['selection'] + ' and name CA') == 5
    # the kinase: residues absent from the model count (REMARK 465),
    # ligands and waters are "others"
    assert kinase['residues'] == count_residues('1oky and polymer')
    assert kinase['residues'] > count_residues('1oky and polymer and present')
    assert kinase['others'] == count_residues('1oky and not polymer') > 100
    assert kinase['sequence'].count('?') == 0
    assert len(kinase['sequence']) == kinase['residues']


def test_offset_keeps_gaps_and_insertion_codes():
    antibody()
    xyz = cmd.get_coords('fab_1 and chain H').copy()
    change = modify.renumber('fab_1 and chain H', offset=10)
    assert numbers('fab_1 and chain H') == ['11', '12', '12A', '14', '15',
                                            '16']
    assert numbers('fab_1 and chain L') == ['1', '2', '3', '4', '5']
    assert letters('fab_1 and chain H') == 'EVQLVE'
    assert cmd.get_coords('fab_1 and chain H') == pytest.approx(xyz)
    assert change.count == 6 and change.chains == [('fab_1', '', 'H')]
    assert 'chain H of fab_1, 1 to 6 is now 11 to 16' in change.text
    # back, across its own old numbers
    modify.renumber('fab_1 and chain H', offset=-10)
    assert numbers('fab_1 and chain H') == ['1', '2', '2A', '4', '5', '6']
    # negative numbers
    modify.renumber('fab_1 and chain L', offset=-3)
    assert numbers('fab_1 and chain L') == ['-2', '-1', '0', '1', '2']
    assert letters('fab_1 and chain L') == 'DIQMT'


def test_start_moves_each_chain_to_a_number():
    antibody()
    change = modify.renumber('fab_1', start=101)
    assert numbers('fab_1 and chain H') == ['101', '102', '102A', '104',
                                            '105', '106']
    assert numbers('fab_1 and chain L') == ['101', '102', '103', '104',
                                            '105']
    assert change.count == 11
    assert change.chains == [('fab_1', '', 'H'), ('fab_1', '', 'L')]
    # nothing to do: said so, and no residues counted
    again = modify.renumber('fab_1', start=101)
    assert again.count == 0 and 'already' in again.text


def test_one_by_one_closes_gaps():
    antibody()
    modify.renumber('fab_1 and chain H', start=3, gaps=False)
    # new numbers overlap the old ones (4 was 5's neighbour's number):
    # done in one pass, no two residues become one
    assert numbers('fab_1 and chain H') == ['3', '4', '5', '6', '7', '8']
    assert letters('fab_1 and chain H') == 'EVQLVE'
    assert cmd.count_atoms('fab_1 and chain H and resi 5 and name CA') == 1


def test_part_of_a_chain():
    antibody()
    modify.renumber('fab_1 and chain L and resi 4-5', offset=10)
    assert numbers('fab_1 and chain L') == ['1', '2', '3', '14', '15']
    # atoms of the selection stand for their residues
    modify.renumber('fab_1 and chain L and resi 14 and name CA', start=20)
    assert numbers('fab_1 and chain L') == ['1', '2', '3', '15', '20']
    assert cmd.count_atoms('fab_1 and chain L and resi 20') == \
        cmd.count_atoms('byres (fab_1 and chain L and resi 20 and name CA)')
    assert cmd.count_atoms('fab_1 and chain L and resi 14') == 0


def test_collisions_change_nothing():
    antibody()
    before = residues('fab_1')
    with pytest.raises(ValueError) as error:
        modify.renumber('fab_1 and chain L and resi 4-5', offset=-2)
    assert str(error.value) == (
        'fab_1: residue 4 of chain L would become residue 2 of chain L, '
        'which is there already')
    assert residues('fab_1') == before
    with pytest.raises(ValueError, match='would become residue 1 of '
                       'chain H, which is there already'):
        modify.set_ids('fab_1 and chain L', chain='H')
    assert residues('fab_1') == before
    # two residues to the same place: both are number 4 of chain L, in
    # different segments
    cmd.alter('fab_1 and chain L and resi 5', 'segi, resi = "S", "4"')
    cmd.sort()
    before = residues('fab_1')
    with pytest.raises(ValueError) as error:
        modify.set_ids('fab_1 and chain L and resi 4', segi='T')
    assert str(error.value) == (
        'fab_1: residue 4 of chain L and residue 4 of chain L (segment S) '
        'would both become residue 4 of chain L (segment T)')
    assert residues('fab_1') == before
    with pytest.raises(ValueError, match='no atoms'):
        modify.renumber('fab_1 and chain X', offset=1)
    with pytest.raises(ValueError, match='either the start or the offset'):
        modify.renumber('fab_1')
    with pytest.raises(ValueError, match='either the start or the offset'):
        modify.renumber('fab_1', start=1, offset=1)


def test_chain_and_segment_ids():
    antibody()
    atoms = cmd.count_atoms('fab_1 and chain H')
    change = modify.set_ids('fab_1 and chain H', chain='A')
    assert cmd.get_chains('fab_1') == ['A', 'L']
    assert cmd.count_atoms('fab_1 and chain A') == atoms
    assert change.count == 6 and change.chains == [('fab_1', '', 'A')]
    assert change.text == '6 residues now have chain A'
    # a segment ID, and both at once
    modify.set_ids('fab_1 and chain L', segi='LC')
    assert {r[0] for r in residues('fab_1 and chain L')} == {'LC'}
    change = modify.set_ids('fab_1 and chain L', chain='B', segi='')
    assert {r[:2] for r in residues('fab_1')} == {('', 'A'), ('', 'B')}
    assert change.text == '5 residues now have chain B and no segment ID'
    # part of a chain becomes a chain of its own
    modify.set_ids('fab_1 and chain A and resi 4-6', chain='C')
    assert letters('fab_1 and chain A') == 'EVQ'
    assert letters('fab_1 and chain C') == 'LVE'
    # and joins another, where its numbers are free
    modify.renumber('fab_1 and chain C', start=11)
    modify.set_ids('fab_1 and chain C', chain='B')
    assert numbers('fab_1 and chain B') == ['1', '2', '3', '4', '5', '11',
                                            '12', '13']
    # no chain ID at all
    modify.set_ids('fab_1 and chain A', chain='')
    assert cmd.get_chains('fab_1') == ['', 'B']
    for wrong in ('A B', 'A"', "A'", '(A)', 'A,B'):
        with pytest.raises(ValueError, match='cannot have'):
            modify.set_ids('fab_1', chain=wrong)
    with pytest.raises(ValueError, match='chain ID or a segment ID'):
        modify.set_ids('fab_1')


def test_swap():
    antibody()
    change = modify.swap('fab_1', 'H', 'L')
    assert letters('fab_1 and chain H') == 'DIQMT'
    assert letters('fab_1 and chain L') == 'EVQLVE'
    assert numbers('fab_1 and chain L') == ['1', '2', '2A', '4', '5', '6']
    assert change.count == 11
    assert 'chains H and L swapped' in change.text
    with pytest.raises(ValueError, match='the same'):
        modify.swap('fab_1', 'H', 'H')
    with pytest.raises(ValueError, match='no chain X or Y'):
        modify.swap('fab_1', 'X', 'Y')


def test_merge_continues_the_numbers():
    antibody()
    change = modify.merge('fab_1')
    assert cmd.get_chains('fab_1') == ['H']
    assert numbers('fab_1') == ['1', '2', '2A', '4', '5', '6', '7', '8',
                                '9', '10', '11']
    assert letters('fab_1') == 'EVQLVEDIQMT'
    assert change.chains == [('fab_1', '', 'H')]
    assert change.text == 'chains H + L of fab_1 are now chain H, 1 to 11'


def test_merge_settings():
    antibody()
    modify.merge('fab_1', chain='A', gap=10, order=['L', 'H'])
    assert cmd.get_chains('fab_1') == ['A']
    # the light chain first, the heavy chain 10 numbers on, as it was
    assert numbers('fab_1') == ['1', '2', '3', '4', '5', '16', '17', '17A',
                                '19', '20', '21']
    assert letters('fab_1') == 'DIQMTEVQLVE'

    cmd.reinitialize()
    antibody()
    modify.merge('fab_1', numbering='consecutive', start=5)
    assert numbers('fab_1') == [str(n) for n in range(5, 16)]
    assert letters('fab_1') == 'EVQLVEDIQMT'

    cmd.reinitialize()
    antibody()
    before = residues('fab_1')
    with pytest.raises(ValueError, match='there already'):
        modify.merge('fab_1', numbering='keep')
    assert residues('fab_1') == before
    modify.renumber('fab_1 and chain L', start=201)
    modify.merge('fab_1', numbering='keep')
    assert numbers('fab_1')[-5:] == ['201', '202', '203', '204', '205']
    assert cmd.get_chains('fab_1') == ['H']


def test_merge_gives_one_segment_and_takes_two_chains():
    antibody()
    cmd.alter('fab_1 and chain L', 'segi = "LC"')
    cmd.sort()
    modify.merge('fab_1')
    assert {r[:2] for r in residues('fab_1')} == {('', 'H')}
    assert letters('fab_1') == 'EVQLVEDIQMT'
    with pytest.raises(ValueError, match='two or more chains'):
        modify.merge('fab_1')
    with pytest.raises(ValueError, match='numbering must be'):
        modify.merge('fab_1', numbering='random')
    with pytest.raises(ValueError, match='gap'):
        modify.merge('fab_1', gap=-1)
    # each object by itself; one with a single chain is left alone
    peptide('one', 'ACD', 'H')
    peptide('two_a', 'GGG', 'H')
    peptide('two_b', 'KKK', 'L')
    cmd.create('two', 'two_a or two_b')
    change = modify.merge('one or two', chain='Z')
    assert cmd.get_chains('one') == ['H']
    assert cmd.get_chains('two') == ['Z']
    assert numbers('two') == ['1', '2', '3', '4', '5', '6']
    assert change.chains == [('two', '', 'Z')]


def test_merge_numbers_ligands_and_waters_after_the_polymers():
    # as in a PDB entry: each chain's waters follow its last residue
    antibody()
    for chain, number in (('H', 7), ('H', 8), ('L', 6)):
        cmd.pseudoatom('fab_1', name='O', resn='HOH', resi=number,
                       chain=chain, segi='', elem='O', pos=(number, 0, 0))
    cmd.alter('fab_1 and resn HOH', 'type = "HETATM"')
    cmd.sort('fab_1')
    assert count_residues('fab_1 and polymer') == 11
    change = modify.merge('fab_1')
    # the light chain goes on where the heavy chain's residues end,
    # not after its waters
    assert numbers('fab_1 and polymer') == ['1', '2', '2A', '4', '5', '6',
                                            '7', '8', '9', '10', '11']
    assert letters('fab_1') == 'EVQLVEDIQMT'
    waters = []
    cmd.iterate_state(1, 'fab_1 and resn HOH',
                      'waters.append((chain, resi, x))',
                      space={'waters': waters})
    assert waters == [('H', '12', 7.0), ('H', '13', 8.0), ('H', '14', 6.0)]
    assert change.text == ('chains H + L of fab_1 are now chain H, 1 to 11, '
                           'its ligands and waters 12 to 14')
    assert change.count == 8


def test_combine_into_one_object():
    peptide('heavy', 'EVQLVE', 'A')
    peptide('light', 'DIQMT', 'A')
    peptide('antigen', 'GSGSG', 'B')
    atoms = cmd.count_atoms('heavy or light or antigen')
    change = modify.combine('complex', 'heavy or light or antigen')
    assert change.name == 'complex'
    assert cmd.count_atoms('complex') == atoms
    # the second "A" is not merged into the first: it gets a free ID
    assert cmd.get_chains('complex') == ['A', 'B', 'C']
    assert letters('complex and chain A') == 'EVQLVE'
    assert letters('complex and chain C') == 'DIQMT'
    assert letters('complex and chain B') == 'GSGSG'
    assert change.renamed == [('light', 'A', 'C')]
    assert 'chain A of light is chain C there' in change.text
    assert change.count == 16
    assert change.chains == [('complex', '', 'A'), ('complex', '', 'B'),
                             ('complex', '', 'C')]
    # copied: the originals are as they were, and no helpers are left
    assert cmd.get_chains('light') == ['A']
    assert names() == ['antigen', 'complex', 'heavy', 'light']
    with pytest.raises(ValueError, match='"complex" already'):
        modify.combine('complex', 'heavy')
    with pytest.raises(ValueError, match='no atoms'):
        modify.combine('other', 'chain X')
    with pytest.raises(ValueError, match='name'):
        modify.combine(' ', 'heavy')
    assert names() == ['antigen', 'complex', 'heavy', 'light']


def test_combine_moves_when_asked():
    peptide('heavy', 'EVQLVE', 'H')
    peptide('light', 'DIQMT', 'L')
    xyz = cmd.get_coords('light').copy()
    # "all" must not take in the new object itself
    modify.combine('fv', 'heavy or (light and resi 1-3)', delete=True)
    assert cmd.get_names('objects') == ['light', 'fv']
    assert letters('light') == 'MT'
    assert letters('fv and chain L') == 'DIQ'
    assert letters('fv and chain H') == 'EVQLVE'
    assert cmd.get_coords('fv and chain L') == pytest.approx(
        xyz[:cmd.count_atoms('fv and chain L')])
    modify.combine('everything', 'all')
    assert cmd.count_atoms('everything') == cmd.count_atoms('light or fv')


def test_commands():
    antibody()
    for func in modify.COMMANDS:
        name = func.__name__
        assert name in cmd.keyword and name in undo.UNDOABLE
        assert 'DESCRIPTION' in func.__doc__ and 'USAGE' in func.__doc__
    assert cmd.renumber_residues('fab_1 and chain H', 11, quiet=1) == 6
    assert numbers('fab_1 and chain H')[0] == '11'
    assert cmd.renumber_residues('fab_1 and chain H', offset='-10',
                                 quiet=1) == 6
    # one by one: only 2A (now 3) gets another number
    assert cmd.renumber_residues('fab_1 and chain H', 1, gaps=0,
                                 quiet=1) == 1
    assert numbers('fab_1 and chain H') == ['1', '2', '3', '4', '5', '6']
    with pytest.raises(ValueError, match='whole number'):
        cmd.renumber_residues('fab_1', 'one')
    assert cmd.swap_chains('H', 'L', quiet=1) == 11
    assert cmd.set_chain_id('chain L', 'X', quiet=1) == 6
    assert cmd.set_segment_id('chain X', 'SEG1', quiet=1) == 6
    assert cmd.merge_chains('fab_1', 'M', order='X H', gap='4',
                            quiet=1) == 11
    assert numbers('fab_1') == [str(n) for n in range(1, 7)] + [
        str(n) for n in range(11, 16)]
    assert {r[:2] for r in residues('fab_1')} == {('SEG1', 'M')}
    assert cmd.combine_chains('copy_1', 'fab_1', quiet=1) == 11
    assert cmd.count_atoms('copy_1') == cmd.count_atoms('fab_1')


def test_commands_from_the_command_line(capsys):
    antibody()
    cmd.do('renumber_residues fab_1 and chain L, offset=100')
    cmd.do('merge_chains fab_1, A, order=L H')
    cmd.do('set_chain_id fab_1, LONG')
    cmd.sync()
    assert cmd.get_chains('fab_1') == ['LONG']
    assert numbers('fab_1') == ['101', '102', '103', '104', '105', '106',
                                '107', '107A', '109', '110', '111']
    out = capsys.readouterr().out
    assert 'renumber_residues: 5 residues renumbered: chain L of fab_1, ' \
        '1 to 5 is now 101 to 105' in out
    assert 'a PDB file holds 1 character of a chain ID (LONG)' in out


def test_one_undo_step_each():
    antibody()
    before = residues('fab_1')
    cmd.undo_enable()
    try:
        steps = len(undo.stack.undo_stack)
        cmd.merge_chains('fab_1', quiet=1)
        assert len(undo.stack.undo_stack) == steps + 1
        cmd.combine_chains('copy_1', 'fab_1', delete=1, quiet=1)
        assert len(undo.stack.undo_stack) == steps + 2
        assert cmd.get_names('objects') == ['copy_1']
        cmd.undo()
        assert cmd.get_names('objects') == ['fab_1']
        cmd.undo()
        assert residues('fab_1') == before
        # a change that fails is no step
        with pytest.raises(ValueError):
            cmd.set_chain_id('fab_1 and chain L', 'H', quiet=1)
        assert len(undo.stack.undo_stack) == steps
    finally:
        cmd.undo_disable()


def test_a_real_structure():
    cmd.load(os.path.join(DATA, '1oky.pdb.gz'), '1oky')
    sequence = letters('1oky and polymer')
    bonds = len(cmd.get_model('1oky').bond)
    xyz = cmd.get_coords('1oky').copy()
    ligands = numbers('1oky and not polymer')
    count = count_residues('1oky and polymer')
    # the protein from 1, residues absent from the model included;
    # ligands and waters keep their numbers
    change = modify.renumber('1oky and polymer', start=1, gaps=False)
    # (the last residue has no alpha carbon)
    assert numbers('1oky and polymer') == [str(n) for n in range(1, count)]
    assert cmd.count_atoms('1oky and polymer and resi %d' % count) == 1
    assert change.count == count
    assert letters('1oky and polymer') == sequence
    assert numbers('1oky and not polymer') == ligands
    assert len(cmd.get_model('1oky').bond) == bonds
    assert cmd.get_coords('1oky') == pytest.approx(xyz)
    # the waters as a chain of their own
    waters = cmd.count_atoms('1oky and solvent')
    modify.set_ids('1oky and solvent', chain='W')
    assert cmd.count_atoms('1oky and chain W') == waters
    assert cmd.count_atoms('1oky and chain A and solvent') == 0


def test_dialog():
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import shiboken6
    import pymol.plugins  # noqa: F401 (adds the startup plugin path)
    from pmg_tk.startup.modify_gui.main import ModifyDialog

    antibody()
    peptide('antigen', 'GSGSG', 'H')
    dialog = ModifyDialog(None, _self=cmd)

    def cells(row):
        return [dialog.table.item(row, col).text()
                for col in range(dialog.table.columnCount())]

    def select(*rows):
        dialog.table.clearSelection()
        mode = dialog.table.selectionModel().SelectionFlag
        for row in rows:
            dialog.table.selectionModel().select(
                dialog.table.model().index(row, 0), mode.Select | mode.Rows)

    try:
        assert dialog.windowTitle() == 'Modify Sequences'
        dialog.refresh()
        assert dialog.table.rowCount() == 3
        assert cells(0) == ['fab_1', 'H', '', '6', '1 to 6 (1 gap)', '',
                            'EVQLVE']
        assert cells(2)[:2] == ['antigen', 'H']
        # nothing selected: nothing changes, and it says why
        assert dialog.renumber() is None
        assert dialog.status.text() == ('Nothing changed: select one or '
                                        'more chains in the table first')

        # renumber: the three ways
        select(0)
        dialog.number.setValue(101)
        dialog.renumber()
        assert numbers('fab_1 and chain H') == ['101', '102', '102A', '104',
                                                '105', '106']
        assert numbers('antigen') == ['1', '2', '3', '4', '5']
        assert cells(0)[4] == '101 to 106 (1 gap)'
        assert dialog.status.text().startswith('6 residues renumbered')
        assert [r['chain'] for r in dialog.selected_rows()] == ['H']
        dialog.mode.setCurrentIndex(dialog.mode.findData('offset'))
        dialog.number.setValue(-100)
        dialog.renumber()
        assert numbers('fab_1 and chain H')[:3] == ['1', '2', '2A']
        dialog.mode.setCurrentIndex(dialog.mode.findData('count'))
        dialog.number.setValue(1)
        dialog.renumber()
        assert numbers('fab_1 and chain H') == ['1', '2', '3', '4', '5', '6']
        assert cells(0)[4] == '1 to 6'

        # only some residues
        dialog.residues.setText('5-6')
        dialog.mode.setCurrentIndex(dialog.mode.findData('start'))
        dialog.number.setValue(10)
        dialog.renumber()
        assert numbers('fab_1 and chain H') == ['1', '2', '3', '4', '10',
                                                '11']
        dialog.residues.setText('chain L')
        assert dialog.renumber() is None
        assert 'neither a number nor a range' in dialog.status.text()
        dialog.residues.setText('')

        # chain and segment IDs
        select(2)
        dialog.chain_id.setText('G')
        dialog.segi_id.setText('AG')
        dialog.segi_id.textEdited.emit('AG')
        assert dialog.segi_check.isChecked()
        dialog.change_ids()
        assert {r[:2] for r in residues('antigen')} == {('AG', 'G')}
        assert cells(2)[:3] == ['antigen', 'G', 'AG']
        assert [r['chain'] for r in dialog.selected_rows()] == ['G']
        dialog.chain_check.setChecked(False)
        dialog.segi_check.setChecked(False)
        assert dialog.change_ids() is None
        assert 'check Chain ID' in dialog.status.text()
        # a clash is reported, not made
        select(1)
        dialog.chain_check.setChecked(True)
        dialog.chain_id.setText('H')
        assert dialog.change_ids() is None
        assert dialog.status.text().startswith(
            'Nothing changed: fab_1: residue 1 of chain L would become')
        assert cmd.get_chains('fab_1') == ['H', 'L']

        # swap two
        select(0)
        assert not dialog.swap_button.isEnabled()
        select(0, 1)
        assert dialog.swap_button.isEnabled()
        assert dialog.order.text() == 'H L'
        dialog.swap()
        assert letters('fab_1 and chain H') == 'DIQMT'
        dialog.swap()
        assert letters('fab_1 and chain H') == 'EVQLVE'

        # merge into one chain, in the order typed
        assert [r['chain'] for r in dialog.selected_rows()] == ['H', 'L']
        dialog.order.setText('L, H')
        dialog.merge_chain.setText('F')
        dialog.gap.setValue(4)
        assert dialog.gap.isEnabled()
        dialog.merge()
        assert cmd.get_chains('fab_1') == ['F']
        assert numbers('fab_1') == ['1', '2', '3', '4', '5', '10', '11',
                                    '12', '13', '19', '20']
        assert dialog.table.rowCount() == 2
        assert cells(0)[:2] == ['fab_1', 'F']
        dialog.numbering.setCurrentIndex(
            dialog.numbering.findData('consecutive'))
        assert not dialog.gap.isEnabled()
        select(0)
        assert dialog.merge() is None
        assert 'two or more chains' in dialog.status.text()

        # into one object
        select(0, 1)
        dialog.new_name.setText('complex')
        dialog.combine()
        assert cmd.get_names('objects') == ['fab_1', 'antigen', 'complex']
        assert dialog.table.rowCount() == 4
        assert cmd.get_chains('complex') == ['F', 'G']
        assert {r['object'] for r in dialog.selected_rows()} == {'complex'}
        assert dialog.combine() is None
        assert 'already' in dialog.status.text()

        # the filter hides objects, and what is hidden is not changed
        dialog.filter.setText('anti')
        assert [dialog.table.isRowHidden(k) for k in range(4)] == [
            True, False, True, True]
        dialog._select_all(True)
        assert [r['object'] for r in dialog.selected_rows()] == ['antigen']
        dialog.filter.setText('')

        # with no undo to take back, it says so
        dialog.undo()
        assert dialog.status.text() == 'Nothing to undo'
        dialog.resize(1040, 640)
        assert not dialog.grab().isNull()
    finally:
        shiboken6.delete(dialog)


def test_dialog_undo():
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import shiboken6
    import pymol.plugins  # noqa: F401 (adds the startup plugin path)
    from pmg_tk.startup.modify_gui.main import ModifyDialog

    antibody()
    dialog = ModifyDialog(None, _self=cmd)
    cmd.undo_enable()
    try:
        dialog.refresh()
        dialog._select_all(True)
        steps = len(undo.stack.undo_stack)
        dialog.merge()
        assert cmd.get_chains('fab_1') == ['H']
        assert len(undo.stack.undo_stack) == steps + 1
        # a change refused is no step
        dialog.new_name.setText('fab_1')
        assert dialog.combine() is None
        assert len(undo.stack.undo_stack) == steps + 1
        dialog.undo()
        assert dialog.status.text() == 'Undone'
        assert cmd.get_chains('fab_1') == ['H', 'L']
        assert dialog.table.rowCount() == 2
    finally:
        cmd.undo_disable()
        shiboken6.delete(dialog)
