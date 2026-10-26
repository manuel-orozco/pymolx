'''
PISA-like interface analysis (pymolx.interfaces) and its plugin dialog
(data/startup/interface_gui). Uses the 1tii demo structure (AB5 toxin).
'''

import csv
import os

import pytest

from pymol import cmd
from pymolx import interfaces, undo

DEMO = cmd.exp_path('$PYMOL_DATA/demo/1tii.pdb')


@pytest.fixture(scope='module')
def result():
    cmd.reinitialize()
    cmd.load(DEMO, 'm')
    return {tuple(sorted(i.chains)): i for i in interfaces.analyze('m')}


@pytest.fixture
def session():
    cmd.reinitialize()
    cmd.load(DEMO, 'm')
    yield
    cmd.reinitialize()


def test_interfaces_found(result):
    # B5 ring (D-E, E-F, F-G, G-H, D-H) and A1-A2 (A-C) are the large ones
    large = {k for k, i in result.items() if i.area > 900}
    assert large == {('D', 'E'), ('E', 'F'), ('F', 'G'), ('G', 'H'),
                     ('D', 'H'), ('A', 'C')}
    assert len(result) == 14  # as PISA finds in the asymmetric unit


def test_ring_interfaces_are_alike(result):
    ring = [result[k] for k in [('D', 'E'), ('E', 'F'), ('F', 'G'),
                                ('G', 'H'), ('D', 'H')]]
    areas = [i.area for i in ring]
    assert max(areas) - min(areas) < 0.1 * min(areas)
    assert all(i.dg < 0 for i in ring)  # hydrophobic burial: favorable


def test_bonds(result):
    # A1 and A2 subunits are linked by a disulfide
    assert result[('A', 'C')].counts()[2] == 1
    a, b, d = result[('A', 'C')].disulfides[0]
    assert (a.resn, a.name, b.resn, b.name) == ('CYS', 'SG', 'CYS', 'SG')
    for iface in result.values():
        for a, b, d in iface.hbonds:
            assert d <= interfaces.HBOND_CUTOFF
            assert a.chain != b.chain
        for a, b, d in iface.salt_bridges:
            assert d <= interfaces.SALT_BRIDGE_CUTOFF
            assert {a.is_anion(), b.is_anion()} == {True, False}


def test_area_and_dg_are_consistent(result):
    for iface in result.values():
        bsa = sum(r['bsa'] for r in iface.residues)
        assert bsa == pytest.approx(2 * iface.area, rel=1e-6)
        dg = -sum(interfaces.ASP.get(t, 0.0) * area
                  for t, area in iface.bsa_by_type.items()) / 1000.0
        assert iface.dg == pytest.approx(dg, abs=1e-6)
        assert sum(r['dg'] for r in iface.residues) == pytest.approx(
            iface.dg, abs=1e-6)
        for r in iface.residues:
            assert 0 < r['bsa'] <= r['asa'] + 1e-6
            assert set(r['flags']) <= set('HSD')


def test_atom_typing():
    def atom(resn, name, elem):
        return interfaces.Atom(1, 'A', resn, '1', name, elem, '')
    assert atom('ASP', 'OD1', 'O').asp_type() == 'O-'
    assert atom('LYS', 'NZ', 'N').asp_type() == 'N+'
    assert atom('LEU', 'CD1', 'C').asp_type() == 'C'
    assert atom('MET', 'SD', 'S').asp_type() == 'S'
    assert atom('SER', 'OG', 'O').asp_type() == 'O'
    assert atom('PRO', 'N', 'N').is_donor() is False
    assert atom('GLY', 'N', 'N').is_donor()
    assert atom('GLY', 'O', 'O').is_acceptor()
    assert atom('PHE', 'CZ', 'C').is_donor() is False
    assert atom('ARG', 'NH2', 'N').is_cation()
    assert atom('GLU', 'OE2', 'O').is_anion()


def test_no_side_effects(session):
    settings = {n: cmd.get(n) for n in ('dot_solvent', 'dot_density')}
    view = cmd.get_view()
    cmd.undo_enable()
    try:
        interfaces.analyze('m')
        assert not undo.stack.can_undo()
    finally:
        cmd.undo_disable()
    assert cmd.get_names('all') == ['m']
    assert {n: cmd.get(n) for n in settings} == settings
    assert cmd.get_view() == pytest.approx(view)


def test_single_chain(session):
    cmd.create('a_only', 'm and chain A')
    assert interfaces.analyze('a_only') == []


def test_command(session, capsys):
    result = cmd.interface_analysis('m')
    assert len(result) == 14
    out = capsys.readouterr().out
    assert 'Interfaces of "m"' in out and 'G-H' in out


def test_show_interface(session):
    iface = [i for i in interfaces.analyze('m') if i.chains == ('A', 'C')][0]
    group = interfaces.show_interface(iface)
    assert group == 'iface_m_AC'
    names = cmd.get_names('all')
    assert 'iface_m_AC_res' in names and 'iface_m_AC_disulf' in names
    assert cmd.count_atoms('iface_m_AC_res and rep sticks') > 0
    residues = {(r['chain'], r['resi']) for r in iface.residues}
    selected = set()
    cmd.iterate('iface_m_AC_res', 'selected.add((chain, resi))',
                space={'selected': selected})
    assert selected == residues


def test_write_csv(result, tmp_path):
    pytest.importorskip('pymol.Qt')
    import pymol.plugins  # noqa: F401 (adds the startup plugin path)
    from pmg_tk.startup.interface_gui.main import write_csv
    path = tmp_path / 'ifaces.csv'
    write_csv(str(path), 'm', list(result.values()))
    rows = list(csv.DictReader(open(path)))
    assert len(rows) == sum(len(i.residues) for i in result.values())
    assert {'chains', 'area', 'dG', 'resn', 'bsa', 'bonds'} <= set(rows[0])


def test_dialog(session):
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import pymol.plugins  # noqa: F401
    from pmg_tk.startup.interface_gui.main import InterfaceDialog

    dialog = InterfaceDialog(None, _self=cmd)
    try:
        dialog.refresh_objects()
        assert dialog.object_combo.currentText() == 'm'
        dialog.analyze()
        assert dialog.interface_table.rowCount() == 14
        assert dialog.selected_interface() is dialog.interfaces[0]
        assert dialog.residue_table.rowCount() == len(
            dialog.interfaces[0].residues)
        assert dialog.bond_table.rowCount() == sum(
            dialog.interfaces[0].counts())
        dialog.show_selected()
        assert any(n.startswith('iface_m_') for n in cmd.get_names('all'))
    finally:
        dialog.deleteLater()


# epitope

def buried(result, target, partners):
    '''
    {resi} of the target chain buried in its interfaces with the partners
    '''
    found = set()
    for chains, iface in result.items():
        if target in chains and set(chains) - {target} <= set(partners) \
                and len(set(chains)) == 2:
            found |= {r['resi'] for r in iface.residues
                      if r['chain'] == target and r['bsa'] > 0}
    return found


def test_epitope_residues(result):
    chains = sorted({c for pair in result for c in pair})
    target = chains[0]
    others = chains[1:]
    by_chain = interfaces.contact_residues(list(result.values()), [target])
    assert list(by_chain) == [target]
    resis = [resi for resi, _, _ in by_chain[target]]
    assert set(resis) == buried(result, target, others)
    assert resis == sorted(resis, key=lambda r: int(r))  # in order
    # more buried area asked: fewer residues
    fewer = interfaces.contact_residues(list(result.values()), [target], 20)
    assert 0 < len(fewer[target]) < len(resis)
    assert all(bsa > 20 for _, _, bsa in fewer[target])
    # the partners' side: only residues buried by the target
    partners = interfaces.contact_residues(list(result.values()), [target],
                                           side='partner')
    assert target not in partners
    for chain, residues in partners.items():
        assert {r for r, _, _ in residues} == buried(result, chain, [target])


def test_interfaces_between_targets_dont_count(result):
    chains = sorted({c for pair in result for c in pair})
    # two touching chains both as targets: their mutual interface is left
    # out, only what the others bury counts
    pair = next(p for p in result if len(set(p)) == 2)
    by_chain = interfaces.contact_residues(list(result.values()), list(pair))
    others = [c for c in chains if c not in pair]
    for chain in pair:
        assert {r for r, _, _ in by_chain.get(chain, [])} == \
            buried(result, chain, others)


def test_epitope_format():
    residues = [('-3', 'GLY', 5.0), ('52', 'ASP', 30.0),
                ('100A', 'TYR', 12.0), ('7', 'HOH', 1.0)]
    assert interfaces.format_residues(residues) == 'G-3, D52, Y100A, X7'
    assert interfaces.format_residues(residues, 'three') == \
        'GLY-3, ASP52, TYR100A, HOH7'
    assert interfaces.format_residues(residues, 'none') == '-3, 52, 100A, 7'
    assert interfaces.format_epitope({'A': residues[:3]}, 'none') == \
        '-3, 52, 100A'
    assert interfaces.format_epitope({'A': residues[:1],
                                      'B': residues[1:3]}) == \
        'A: G-3\nB: D52, Y100A'


def test_default_targets_are_the_non_antibody_chains():
    # trastuzumab's variable domains (heavy, light), and an antigen
    VH = ('EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADS'
          'VKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS')
    VL = ('DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFS'
          'GSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK')
    ANTIGEN = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVK'
    cmd.reinitialize()
    try:
        cmd.fab(VH, 'h', chain='H')
        cmd.fab(VL, 'l', chain='L')
        cmd.fab(ANTIGEN, 'a', chain='A')
        cmd.create('cpx', 'h or l or a')
        assert interfaces.default_targets('cpx') == ['A']
        assert interfaces.default_targets('a') == []
    finally:
        cmd.reinitialize()


def test_epitope_command(session, capsys):
    chains = sorted(cmd.get_chains('m and polymer'))
    text = cmd.epitope('m', chains[0])
    out = capsys.readouterr().out
    assert 'epitope, chain %s' % chains[0] in out
    assert text == out.splitlines()[-1].strip()
    assert cmd.epitope('m', '%s+%s' % tuple(chains[:2]), quiet=1).startswith(
        '%s: ' % chains[0])
    one = cmd.epitope('m', chains[0], quiet=1).split(', ')
    three = cmd.epitope('m', chains[0], names='three', quiet=1).split(', ')
    plain = cmd.epitope('m', chains[0], names='none', quiet=1).split(', ')
    assert [n[1:] for n in one] == [n[3:] for n in three] == plain
    assert one[0][0].isalpha() and three[0][:3].isalpha()
    with pytest.raises(ValueError, match='names'):
        cmd.epitope('m', chains[0], names='four')
    with pytest.raises(ValueError, match='no chain Z'):
        cmd.epitope('m', 'Z')
    with pytest.raises(ValueError, match='give the target'):
        cmd.epitope('m')  # no antibody chains here


def test_dialog_epitope_tab(session):
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import pymol.plugins  # noqa: F401
    from pmg_tk.startup.interface_gui.main import InterfaceDialog
    Qt = pytest.importorskip('pymol.Qt').QtCore.Qt

    dialog = InterfaceDialog(None, _self=cmd)
    try:
        dialog.refresh_objects()
        dialog.analyze()
        assert [dialog.tabs.tabText(i) for i in range(dialog.tabs.count())] \
            == ['Residues', 'Bonds', 'Epitope']
        chains = sorted(cmd.get_chains('m and polymer'))
        items = [dialog.target_list.item(i)
                 for i in range(dialog.target_list.count())]
        assert [i.text() for i in items] == chains
        assert dialog.targets() == []          # no antibody: none chosen
        assert dialog.epitope_text.toPlainText() == ''
        items[0].setCheckState(Qt.CheckState.Checked)
        expected = cmd.epitope('m', chains[0], quiet=1)
        assert dialog.epitope_text.toPlainText() == expected
        assert dialog.name_style.currentData() == 'one'
        dialog.name_style.setCurrentIndex(dialog.name_style.findData('three'))
        assert dialog.epitope_text.toPlainText() == cmd.epitope(
            'm', chains[0], names='three', quiet=1)
        dialog.name_style.setCurrentIndex(dialog.name_style.findData('one'))
        assert dialog.paratope_text.toPlainText()
        assert dialog.epitope_summary.text().startswith(
            '%d epitope residues' % len(expected.split(', ')))
        dialog.min_bsa.setValue(20)
        assert len(dialog.epitope_text.toPlainText().split(', ')) < \
            len(expected.split(', '))
        dialog.select_epitope()
        assert cmd.count_atoms('epitope') > 0
        assert cmd.count_atoms('paratope') > 0
        assert cmd.count_atoms('epitope and not chain %s' % chains[0]) == 0
    finally:
        dialog.deleteLater()
