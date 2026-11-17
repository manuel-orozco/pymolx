'''
CDR loops (pymolx.cdr, data/startup/cdr_gui, parity item L-31): ANARCI
numbering under each scheme, where the loops fall, the command on a
structure, and the dialog.
'''

import os

import pytest

from pymol import cmd
import pymol.plugins  # noqa: F401 (puts $PYMOL_DATA/startup on the path)
from pymolx import cdr

# trastuzumab, whose CDRs are published under each scheme
VH = ('EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG'
      'RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS')
VL = ('DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR'
      'SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK')
NOT_AN_ANTIBODY = ('MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAF'
                   'AQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVA')

needs_anarci = pytest.mark.skipif(bool(cdr.available()),
                                  reason=cdr.available() or '')


@pytest.fixture(autouse=True)
def session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


def test_schemes_listed():
    assert [v for _, v in cdr.SCHEMES] == ['imgt', 'kabat', 'chothia',
                                           'martin', 'aho']
    assert set(cdr.REGIONS) == {v for _, v in cdr.SCHEMES}
    for scheme, kinds in cdr.REGIONS.items():
        assert set(kinds) == {'H', 'L'}
        for ranges in kinds.values():
            assert len(ranges) == 3            # CDR1, CDR2, CDR3
            assert all(a <= b for a, b in ranges)
            # the loops don't overlap and come in order
            assert ranges[0][1] < ranges[1][0] < ranges[1][1] < ranges[2][0]
    assert cdr.scheme_label('kabat') == 'Kabat'


@needs_anarci
def test_published_cdrs():
    # the CDRs of trastuzumab as each scheme defines them
    expected = {
        'kabat': {'CDR1': 'DTYIH', 'CDR2': 'RIYPTNGYTRYADSVKG',
                  'CDR3': 'WGGDGFYAMDY'},
        'chothia': {'CDR1': 'GFNIKDT', 'CDR2': 'YPTNGY',
                    'CDR3': 'WGGDGFYAMDY'},
        'imgt': {'CDR1': 'GFNIKDTY', 'CDR2': 'IYPTNGYT',
                 'CDR3': 'SRWGGDGFYAMDY'},
    }
    for scheme, loops in expected.items():
        domain, = cdr.domains(VH, scheme)
        assert domain.chain_type == 'H'
        assert domain.kind == 'H'
        assert domain.is_antibody
        assert domain.cdrs() == loops, scheme
    # the light chain, and that it is recognised as kappa
    domain, = cdr.domains(VL, 'kabat')
    assert domain.chain_type == 'K' and domain.kind == 'L'
    assert domain.cdrs() == {'CDR1': 'RASQDVNTAVA', 'CDR2': 'SASFLYS',
                             'CDR3': 'QQHYTTPPT'}


@needs_anarci
def test_regions_cover_the_domain():
    domain, = cdr.domains(VH, 'imgt')
    regions = domain.regions()
    assert [name for name, _ in regions] == [
        'FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']
    # every numbered residue lands in exactly one region, in order
    joined = ''.join(''.join(letter for _, letter in residues)
                     for _, residues in regions)
    assert joined == domain.sequence()
    assert domain.sequence() in VH
    # the numbering runs forward
    numbers = [position for _, residues in regions
               for position, _ in residues]
    assert numbers == sorted(numbers)


@needs_anarci
def test_not_an_antibody():
    assert cdr.domains(NOT_AN_ANTIBODY, 'imgt') == []
    with pytest.raises(ValueError, match='scheme'):
        cdr.domains(VH, 'nonsense')
    with pytest.raises(ValueError, match='no sequence'):
        cdr.domains('', 'imgt')


@needs_anarci
def test_two_domains_in_an_scfv():
    # an scFv: ANARCI finds the two variable domains in one chain
    found = cdr.domains(VH + 'GGGGS' * 3 + VL, 'imgt')
    assert len(found) == 2
    assert [d.chain_type for d in found] == ['H', 'K']
    assert found[0].end < found[1].start      # in order, not overlapping


@needs_anarci
def test_find_cdrs_on_a_structure(tmp_path):
    cmd.fab(VH, 'ab', chain='H', ss=1)
    found = cmd.find_cdrs('ab', scheme='kabat', quiet=1)
    assert len(found) == 1
    chain_id, domain, residues = found[0]
    assert chain_id == 'H'
    assert len(residues) == len(domain.positions)
    # a selection per loop, and one of all of them
    for name in cdr.CDR_NAMES:
        selection = 'ab_H_%s' % name
        assert selection in cmd.get_names('selections')
        letters = []
        cmd.iterate('%s and name CA' % selection,
                    'letters.append(resn)', space={'letters': letters})
        assert len(letters) == len(domain.cdrs()[name])
    assert cmd.count_atoms('ab_cdrs') > 0
    # the loops are coloured apart from the framework
    colours = set()
    cmd.iterate('ab_H_CDR3', 'colours.add(color)', space={'colours': colours})
    assert colours == {cmd.get_color_index(cdr.CDR_COLORS[2])}
    # the report names the scheme and the loops
    text = cdr.report(found, 'kabat')
    assert 'Kabat' in text and 'CDR3' in text
    assert domain.cdrs()['CDR3'] in text


@needs_anarci
def test_dialog():
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import shiboken6
    from pmg_tk.startup.cdr_gui.main import CDRDialog, _read_sequences

    assert _read_sequences('>a\nEVQL\nVESG\n>b\nDIQM') == [
        ('a', 'EVQLVESG'), ('b', 'DIQM')]
    assert _read_sequences('EVQLVESG') == [('sequence', 'EVQLVESG')]
    assert _read_sequences('   ') == []

    cmd.fab(VH, 'ab', chain='H', ss=1)
    dialog = CDRDialog(None, _self=cmd)
    try:
        dialog.refresh()
        dialog.object.setCurrentText('ab')
        dialog.scheme.setCurrentIndex(dialog.scheme.findData('chothia'))
        dialog.run()
        assert dialog.table.rowCount() == 3          # the three loops
        rows = {dialog.table.item(r, 2).text():
                dialog.table.item(r, 3).text()
                for r in range(3)}
        assert rows['CDR1'] == 'GFNIKDT'
        assert 'Chothia' in dialog.status.text()
        # the lengths column
        assert dialog.table.item(0, 4).text() == str(len(rows['CDR1']))

        # a pasted sequence instead
        dialog.source.button(1).setChecked(True)
        dialog.sequence.setPlainText('>mine\n' + VL)
        dialog.scheme.setCurrentIndex(dialog.scheme.findData('imgt'))
        dialog.run()
        assert dialog.table.rowCount() == 3
        assert dialog.table.item(0, 0).text() == 'mine'
        assert 'IMGT' in dialog.status.text()

        # something that isn't an antibody is said so, not crashed on
        dialog.sequence.setPlainText(NOT_AN_ANTIBODY)
        dialog.run()
        assert 'No antibody' in dialog.status.text()
    finally:
        dialog.close()
        shiboken6.delete(dialog)
        QtWidgets.QApplication.processEvents()
