'''
Residues missing from the model (INCENTIVE_PARITY.md L-07): PDB REMARK
465 kept by "load", mmCIF residues without coordinates, and how the
sequence viewer data shows them.
'''

import os

import pytest

from pymol import cmd
from pymolx import missing, sequence, undo

DATA = os.path.join(os.path.dirname(__file__), '..', '..', 'data')


@pytest.fixture(autouse=True)
def clean_session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


def chain_text(chain, model=None):
    '''modeled residues upper case, missing ones lower case'''
    row = [r for r in sequence.get_rows()
           if r.chain == chain and (model is None or r.model == model)][0]
    return ''.join(r.text if r.present else r.text.lower()
                   for r in row.residues if r.polymer)


def test_parse_remark465():
    lines = [
        'REMARK 465 MISSING RESIDUES',
        'REMARK 465   M RES C SSSEQI',
        'REMARK 465     MET A    51',
        'REMARK 465     ASP A    52A',
        'REMARK 465     GLY B    -3',
        'REMARK 465     MET A    51',  # repeated (NMR models)
        'REMARK 470 MISSING ATOM',
        'REMARK 465     SER A    99',  # after the table: ignored
    ]
    assert missing.parse_remark465(lines) == [
        ('A', '51', 51, 'MET'), ('A', '52A', 52, 'ASP'),
        ('B', '-3', -3, 'GLY')]
    assert missing.parse_remark465(['REMARK 465 NONE']) == []


def missing_atoms(selection):
    return cmd.count_atoms('(%s) and not present' % selection)


def test_pdb_load_adds_remark465_residues():
    cmd.load(os.path.join(DATA, '1oky.pdb.gz'))
    # one CA atom without coordinates per missing residue
    assert missing_atoms('1oky') == 27
    assert cmd.count_atoms('1oky and not present and not name CA') == 0

    # in residue order: the missing N-terminus (51-71), then the model
    resis = []
    cmd.iterate('1oky and chain A and name CA',
                'resis.append(int(resv))', space={'resis': resis})
    assert resis[:22] == list(range(51, 73))
    assert resis == sorted(resis)

    text = chain_text('A')
    assert text.startswith('mdgtaaeprpgagslqhaqppPQPRKK')
    assert sum(c.islower() for c in text) == 27

    # atoms: they follow renames and sessions
    cmd.set_name('1oky', 'kinase')
    session = cmd.get_session()
    cmd.reinitialize()
    cmd.set_session(session)
    assert missing_atoms('kinase') == 27


def test_pdb_load_into_named_object_and_command_language():
    cmd.do('load %s, mykinase' % os.path.join(DATA, '1oky.pdb.gz'))
    cmd.sync()
    assert missing_atoms('mykinase') == 27


def test_other_formats_untouched():
    cmd.load(os.path.join(DATA, '1hbb_entity_poly_seq.cif'), 'hb')
    assert missing_atoms('hb') == 227  # PyMOL's own, not more
    cmd.load(os.path.join(DATA, '1rx1.pdb'), 'dhfr')  # no REMARK 465
    assert missing_atoms('dhfr') == 0


def test_retain_order_is_respected():
    # like the mmCIF reader: no added atoms when atom order is kept
    cmd.set('retain_order')
    cmd.load(os.path.join(DATA, '1oky.pdb.gz'))
    assert missing_atoms('1oky') == 0


def test_mmcif_residues_without_coordinates():
    # PyMOL adds missing residues as CA atoms without coordinates
    cmd.load(os.path.join(DATA, '1hbb_entity_poly_seq.cif'), 'hb')
    text = chain_text('A')
    assert text.startswith('vlSPADKTNVKAAWGKV')
    assert 'Ktyfphfdlshgsaqvkghgkkvadaltnavahvddmpnalsalsdlhahklrvdpvnfkl' \
        in text
    assert len(text) == 141  # the whole sequence


def test_missing_residues_can_be_selected():
    cmd.load(os.path.join(DATA, '1oky.pdb.gz'))
    row = [r for r in sequence.get_rows() if r.chain == 'A'][0]
    absent = [r for r in row.residues if not r.present]
    sequence.select_residues(absent[:3])
    assert sequence.get_selected() == {r.key for r in absent[:3]}

    # a range across the gap: missing and modeled residues
    cmd.delete('sele')
    first = row.residues.index(absent[0])
    span = row.residues[first:first + 25]  # 21 missing + 4 modeled
    sequence.select_residues(span)
    assert sequence.get_selected() == {r.key for r in span}


def test_load_is_still_one_undo_step():
    cmd.undo_enable()
    try:
        before = len(undo.stack.undo_stack)
        cmd.load(os.path.join(DATA, '1oky.pdb.gz'))
        assert len(undo.stack.undo_stack) == before + 1
        cmd.undo()
        assert cmd.get_names() == []
    finally:
        cmd.undo_disable()


def test_unreadable_extra_data_never_breaks_load(monkeypatch, capsys):
    def broken(*args):
        raise OSError('disk on fire')

    monkeypatch.setattr(missing, '_pdb_lines', broken)
    cmd.load(os.path.join(DATA, '1oky.pdb.gz'))
    assert cmd.count_atoms('1oky') > 0
    assert 'could not read missing residues' in capsys.readouterr().out


def test_messages_note_missing_residues():
    cmd.load(os.path.join(DATA, '1oky.pdb.gz'))
    row = [r for r in sequence.get_rows() if r.chain == 'A'][0]
    i = [k for k, r in enumerate(row.residues) if not r.present][-1]
    grey, mixed = row.residues[i - 2:i + 1], row.residues[i - 2:i + 4]
    assert sequence.describe_residues(grey[:1]) == \
        'You clicked %s (no coordinates)' % grey[0].macro
    assert grey[0].macro.startswith('/1oky//A/')
    assert sequence.describe_residues(grey).endswith('(no coordinates)')
    assert sequence.describe_residues(mixed).endswith(
        '(3 without coordinates)')
    assert 'coordinates' not in sequence.describe_residues(
        row.residues[i + 1:i + 3])
    assert sequence.describe_residues(row.residues, whole=True) == \
        'You clicked /1oky//A/ (27 without coordinates)'


def test_wizards_only_get_residues_with_coordinates():
    from pymol.wizard import Wizard

    picked = []

    class Picker(Wizard):
        def do_select(self, name):
            picked.append(cmd.count_atoms(name + ' and present'))

    cmd.load(os.path.join(DATA, '1oky.pdb.gz'))
    row = [r for r in sequence.get_rows() if r.chain == 'A'][0]
    i = [k for k, r in enumerate(row.residues) if not r.present][-1]
    cmd.set_wizard(Picker(_self=cmd))
    try:
        sequence.select_residues(row.residues[i - 2:i + 1])  # grey only
        assert picked == []
        assert cmd.count_atoms('sele') == 3  # still selected
        sequence.select_residues(row.residues[i - 2:i + 4])  # mixed
        assert len(picked) == 1 and picked[0] > 0
    finally:
        cmd.set_wizard()


def test_mutagenesis_unharmed_by_missing_residues():
    cmd.load(os.path.join(DATA, '1oky.pdb.gz'))
    row = [r for r in sequence.get_rows() if r.chain == 'A'][0]
    grey = [r for r in row.residues if not r.present][0]
    modeled = [r for r in row.residues if r.present and r.polymer][3]
    cmd.wizard('mutagenesis')
    try:
        sequence.select_residues([grey])  # used to raise in the wizard
        assert cmd.count_atoms('sele') == 1
        # a modeled residue: the wizard takes (and deletes) the selection;
        # the click is still reported
        sequence.select_residues([modeled])
        assert cmd.count_states('mutation') > 1
        import io, contextlib
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            sequence.report_selection(sequence.describe_residues([modeled]),
                                      [modeled], True)
        assert 'You clicked ' + modeled.macro in out.getvalue()
    finally:
        cmd.set_wizard()
