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


# proximity profile

def test_proximity(session):
    import numpy
    atoms = interfaces.residue_atoms('m')
    chains = sorted(cmd.get_chains('m and polymer'))
    assert atoms.chains() == chains
    target = chains[0]
    profile = interfaces.proximity(atoms, [target])
    assert list(profile) == [target]
    rows = profile[target]
    # every residue of the chain, in the object's order
    resis = []
    cmd.iterate('m and polymer and chain %s and name CA' % target,
                'resis.append(resi)', space={'resis': resis})
    assert [r[0] for r in rows] == resis
    # the distance is the closest pair of heavy atoms
    others = cmd.get_coords('m and polymer and not hydro and not chain %s' %
                            target)
    for resi, resn, distance, partner in rows[::15]:
        own = cmd.get_coords('m and polymer and not hydro and chain %s and '
                             'resi %s' % (target, resi))
        expected = numpy.sqrt(((own[:, None] - others[None])**2).sum(-1))
        assert distance == pytest.approx(expected.min(), abs=1e-3)
        assert partner[0] != target
        assert cmd.count_atoms('m and chain %s and resi %s and resn %s' %
                               partner) > 0
    # the residues buried in an interface are the close ones
    result = interfaces.analyze('m')
    epitope = {r[0] for r in interfaces.contact_residues(
        result, [target])[target]}
    close = {r[0] for r in rows if r[2] <= 4.0}
    assert close and close <= epitope
    assert all(r[2] < 7.0 for r in rows if r[0] in epitope)
    # two targets: measured to the rest only
    both = interfaces.proximity(atoms, chains[:2])
    assert list(both) == chains[:2]
    assert all(r[3][0] not in chains[:2] for rows in both.values()
               for r in rows)
    assert all(a[2] >= b[2] - 1e-9 for a, b in zip(both[target], rows))
    # nothing to measure
    assert interfaces.proximity(atoms, []) == {}
    assert interfaces.proximity(atoms, chains) == {}
    assert interfaces.proximity(atoms, ['?']) == {}


def test_nearest_without_scipy():
    import numpy
    random = numpy.random.default_rng(1)
    a = random.uniform(-40, 40, (700, 3))
    b = random.uniform(-40, 40, (900, 3))
    tree = interfaces._nearest(a, b)
    plain = interfaces._nearest(a, b, tree=False)
    assert numpy.allclose(tree[0], plain[0], atol=1e-6)
    assert (numpy.asarray(tree[1]) == plain[1]).all()


def test_proximity_score_and_table():
    assert interfaces.proximity_score(4.0) == pytest.approx(0.25)
    assert interfaces.proximity_score(10.0) == pytest.approx(0.1)
    assert interfaces.proximity_score(2.0) > interfaces.proximity_score(3.0)
    # atoms on top of each other (clashes) score no higher than at 1 A
    assert interfaces.proximity_score(0.0) == interfaces.proximity_score(1.0)
    table = interfaces.format_proximity({
        'A': [('52', 'ASP', 3.125, ('H', '100A', 'TYR')),
              ('53', 'GLY', 12.5, ('L', '32', 'SER'))]}).split('\n')
    assert table[0].split('\t')[:5] == ['chain', 'residue', 'number',
                                        'distance (A)', 'score (1/A)']
    assert table[1].split('\t') == ['A', 'ASP', '52', '3.12', '0.320', 'H',
                                    'TYR', '100A']
    assert table[2].split('\t') == ['A', 'GLY', '53', '12.50', '0.080', 'L',
                                    'SER', '32']


def test_color_proximity(session):
    from pymolx.gui import theme
    ramp = [theme.COLORS['proximity_%d' % k] for k in range(5)]
    # the ramp: its ends, a stop and halfway between two
    assert interfaces.ramp_rgb(ramp, 0) == pytest.approx(
        [int(ramp[0][k:k + 2], 16) / 255 for k in (1, 3, 5)])
    assert interfaces.ramp_rgb(ramp, 7) == interfaces.ramp_rgb(ramp, 1)
    assert interfaces.ramp_rgb(ramp, 0.5) == pytest.approx(
        [int(ramp[2][k:k + 2], 16) / 255 for k in (1, 3, 5)])
    between = interfaces.ramp_rgb(['#000000', '#ff8000'], 0.5)
    assert between == pytest.approx((0.5, 64 / 255, 0.0))
    # the colors end at 2.5 A, as the graph's do
    assert interfaces.proximity_fraction(2.5) == pytest.approx(1.0)
    assert interfaces.proximity_fraction(1.0) == 1.0
    assert interfaces.proximity_fraction(10.0) == pytest.approx(0.25)

    chains = sorted(cmd.get_chains('m and polymer'))
    target = chains[0]
    profile = interfaces.proximity(interfaces.residue_atoms('m'), [target])
    rows = profile[target]

    def colors(selection):
        found = set()
        cmd.iterate(selection, 'found.add(color)', space={'found': found})
        return found

    others = 'm and not chain %s' % target
    before = {}
    cmd.iterate(others, 'before[(chain, resi, name)] = color',
                space={'before': before})
    assert interfaces.color_proximity('m', {}, ramp) == 0
    assert interfaces.color_proximity('m', profile, ramp) == len(rows)
    # every residue of the target in the color of its distance, all its
    # atoms; the closest in the ramp's last color, the farthest low on it
    for resi, resn, distance, partner in rows[::9] + [
            min(rows, key=lambda r: r[2]), max(rows, key=lambda r: r[2])]:
        step = round(100 * interfaces.proximity_fraction(distance))
        assert colors('m and chain %s and resi %s' % (target, resi)) == {
            cmd.get_color_index('proximity_%03d' % step)}
    closest = min(r[2] for r in rows)
    assert closest < 2.6 and max(r[2] for r in rows) > 15
    assert cmd.get_color_tuple('proximity_100') == pytest.approx(
        interfaces.ramp_rgb(ramp, 1))
    assert cmd.get_color_tuple('proximity_000') == pytest.approx(
        interfaces.ramp_rgb(ramp, 0))
    assert len(colors('m and chain %s and polymer' % target)) > 10
    # the other chains as they were
    after = {}
    cmd.iterate(others, 'after[(chain, resi, name)] = color',
                space={'after': after})
    assert after == before


def _chart_module():
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import pymol.plugins  # noqa: F401
    from pmg_tk.startup.interface_gui import chart
    return chart


def _brightness(color):
    return (0.2126 * color.redF() + 0.7152 * color.greenF() +
            0.0722 * color.blueF())


def test_chart_axis():
    chart = _chart_module()
    assert chart.residue_number('100A') == (100, 'A')
    assert chart.residue_number('-3') == (-3, '')
    assert chart.residue_number('X') == (None, 'X')
    # one slot per number; an insertion code takes the next one
    assert chart.slots(['1', '2', '3']) == ([0, 1, 2], 3)
    assert chart.slots(['-2', '-1', '1']) == ([0, 1, 3], 4)
    assert chart.slots(['99', '100', '100A', '100B', '101']) == \
        ([0, 1, 2, 3, 4], 5)
    # missing residues leave a gap, but a jump in the numbers a short one
    assert chart.slots(['5', '6', '10']) == ([0, 1, 5], 6)
    assert chart.slots(['5', '1005']) == ([0, chart.MAX_GAP],
                                          chart.MAX_GAP + 1)
    # numbered out of order: next slot
    assert chart.slots(['50', '51', '1', '2']) == ([0, 1, 2, 3], 4)
    assert chart.slots([]) == ([], 0)
    assert chart.number_step(60) == 1
    assert chart.number_step(2.2) == 50
    assert chart.number_step(0.5) == 200
    # the colors run from the theme's first to its last
    from pymolx.gui import theme
    assert chart.ramp_color(0).name() == theme.COLORS[chart.RAMP[0]]
    assert chart.ramp_color(1).name() == theme.COLORS[chart.RAMP[-1]]
    assert chart.ramp_color(7).name() == theme.COLORS[chart.RAMP[-1]]
    assert chart.ramp_color(0.5).name() == theme.COLORS[chart.RAMP[2]]
    # brighter all the way: closer residues stand out on the dark chart
    brightness = [_brightness(chart.ramp_color(k / 20)) for k in range(21)]
    assert brightness == sorted(brightness)
    assert len(set(brightness)) == 21


def test_chart(session):
    import shiboken6
    chart = _chart_module()
    rows = [('10', 'ASP', 3.0, ('H', '100A', 'TYR')),
            ('11', 'GLY', 12.0, ('H', '33', 'SER')),
            ('12', 'LYS', 2.0, ('L', '50', 'ASP')),
            ('20', 'ALA', 40.0, ('L', '50', 'ASP'))]
    widget = chart.ProximityChart()
    try:
        widget.resize(600, 220)
        assert widget.residue_at(100) is None        # nothing shown yet
        assert widget.image().width() == 1200
        widget.set_profile(rows, 'A', 'caption')
        assert widget.slot_count == 11
        assert widget.top_score() == pytest.approx(0.5)  # 2 A: above 0.4
        plot = widget.plot_rect()
        bars = [widget.bar_rect(k) for k in range(4)]
        # closer: taller; bars stand on the axis, inside the plot
        assert bars[2].height() > bars[0].height() > bars[1].height() > \
            bars[3].height() > 0
        assert all(b.bottom() == pytest.approx(plot.bottom()) for b in bars)
        assert all(plot.left() <= b.left() and b.right() <= plot.right() + 1
                   and b.top() >= plot.top() - 1e-6 for b in bars)
        assert bars[0].height() == pytest.approx(
            plot.height() * (1 / 3.0) / 0.5)
        assert bars[0].right() <= bars[1].left() < bars[2].left() < \
            bars[3].left()
        # the residue under the pointer
        for k, bar in enumerate(bars):
            assert widget.residue_at(bar.center().x()) == k
        assert widget.residue_at((bars[2].right() + bars[3].left()) / 2) \
            is None                                    # residues 13 to 19
        assert widget.residue_at(2) is None
        assert widget.describe(0) == ('ASP 10, chain A\n3.0 Å from TYR 100A, '
                                      'chain H\nScore 0.33')
        assert widget.within_cutoff() == 2
        widget.set_cutoff(2.5)
        assert widget.within_cutoff() == 1
        # 47 pixels per residue: every second number is labelled
        assert [t for _, t in widget.number_ticks(plot)] == ['10', '12', '20']
        # a click names the residue
        clicked = []
        widget.residue_clicked.connect(lambda *a: clicked.append(a))
        from pymol.Qt import QtCore, QtGui, QtWidgets
        point = QtCore.QPointF(bars[2].center().x(), 100)
        left = QtCore.Qt.MouseButton.LeftButton
        QtWidgets.QApplication.sendEvent(widget, QtGui.QMouseEvent(
            QtCore.QEvent.Type.MouseButtonPress, point, point, left, left,
            QtCore.Qt.KeyboardModifier.NoModifier))
        assert clicked == [('A', '12')]
        # hovering marks it, leaving the chart ends that
        QtWidgets.QApplication.sendEvent(widget, QtGui.QMouseEvent(
            QtCore.QEvent.Type.MouseMove, point, point,
            QtCore.Qt.MouseButton.NoButton, QtCore.Qt.MouseButton.NoButton,
            QtCore.Qt.KeyboardModifier.NoModifier))
        assert widget.hovered == 2
        QtWidgets.QApplication.sendEvent(
            widget, QtCore.QEvent(QtCore.QEvent.Type.Leave))
        assert widget.hovered is None
        # the picture: bars in the colors of the ramp, close to far
        picture = widget.image(1)
        assert (picture.width(), picture.height()) == (600, 220)

        def color_of(k):
            bar = bars[k]
            return picture.pixelColor(int(bar.center().x()),
                                      int(bar.bottom() - 2))
        assert color_of(2).name() == chart.ramp_color(1).name()
        assert color_of(0).name() == chart.ramp_color(
            (1 / 3.0) / chart.SCORE_TOP).name()
        assert _brightness(color_of(2)) > _brightness(color_of(0)) > \
            _brightness(color_of(1))
        # the short chain fills the plot: nothing to scroll, bars as wide
        # as they get
        assert widget.scroll_range() == 0 and widget.shown() == [0, 1, 2, 3]
        assert not widget.can_zoom_in() and not widget.can_zoom_out()
        widget.set_profile([])
        widget.set_message('nothing')
        assert widget.residue_at(300) is None
        assert not widget.image(1).isNull()
    finally:
        shiboken6.delete(widget)


def test_chart_scrolls_a_long_chain(session):
    import shiboken6
    chart = _chart_module()
    from pymol.Qt import QtCore, QtGui, QtWidgets
    peak_color = chart.ramp_color((1 / 3.0) / chart.SCORE_TOP).name()
    many = [(str(k), 'ALA', 30.0 if k != 700 else 3.0, ('H', '1', 'GLY'))
            for k in range(1, 2001)]
    widget = chart.ProximityChart()
    row = chart.ChartScroll(widget)
    try:
        widget.resize(600, 220)
        widget.set_profile(many, 'A')
        plot = widget.plot_rect()
        width = int(round(plot.width()))
        assert width < 600
        # a pixel per residue, from the start: the rest is scrolled to
        assert widget.slot_width() == chart.SLOT_WIDTH == 1.0
        assert widget.offset == 0
        assert widget.scroll_range() == 2000 - width
        shown = widget.shown()
        assert shown[0] == 0 and 0 <= len(shown) - width <= 1
        assert widget.residue_at(plot.left() + 10.5) == 10
        assert 699 not in shown
        assert all(int(t) <= width + 1 for _, t in widget.number_ticks(plot))
        changes = []
        widget.view_changed.connect(lambda: changes.append(1))
        widget.set_offset(600)
        assert widget.offset == 600 and len(changes) == 1
        peak = widget.bar_rect(699)
        assert plot.left() < peak.left() < plot.right()
        assert widget.residue_at(peak.center().x()) == 699
        assert widget.image(1).pixelColor(
            int(peak.left()), int(peak.top() + 3)).name() == peak_color
        ticks = widget.number_ticks(plot)
        assert ticks and all(600 <= int(t) <= 601 + width for _, t in ticks)
        assert all(plot.left() <= x <= plot.right() for x, _ in ticks)
        # nothing is drawn outside the plot
        left_of_plot = widget.image(1).pixelColor(int(plot.left()) - 3,
                                                  int(plot.bottom()) - 2)
        from pymolx.gui import theme
        assert left_of_plot.name() == theme.COLORS['base']
        # no further than the ends
        widget.set_offset(10**6)
        assert widget.offset == widget.scroll_range()
        assert widget.shown()[-1] == 1999
        widget.set_offset(-5)
        assert widget.offset == 0
        # wider bars: the residue in the middle stays there
        widget.set_offset(600)
        middle = widget.residue_at(plot.center().x())
        widget.zoom_in()
        assert widget.slot_width() == pytest.approx(1.5)
        assert abs(widget.residue_at(plot.center().x()) - middle) <= 1
        assert widget.scroll_range() == 3000 - width
        widget.zoom_out()
        assert widget.slot_width() == pytest.approx(1.0)
        assert abs(widget.residue_at(plot.center().x()) - middle) <= 1
        # the whole chain: residues share pixels, the closest shows
        assert widget.can_zoom_out()
        widget.show_whole()
        assert widget.scroll_range() == 0 and widget.offset == 0
        assert len(widget.shown()) == 2000
        assert not widget.can_zoom_out() and widget.can_zoom_in()
        peak = widget.bar_rect(699)
        assert widget.residue_at(peak.center().x()) == 699
        assert widget.image(1).pixelColor(
            int(peak.left()), int(peak.top() + 3)).name() == peak_color
        # the widest
        for _ in range(30):
            widget.zoom_in()
        assert widget.slot_width() == chart.MAX_SLOT_WIDTH
        assert not widget.can_zoom_in()
        # the same chain again (say, another line): the view stays
        offset = widget.offset
        assert offset > 0
        widget.set_profile(many, 'A', 'caption')
        assert widget.slot_width() == chart.MAX_SLOT_WIDTH
        assert widget.offset == offset
        # another chain: from its start
        widget.set_profile(many[:1500], 'B')
        assert widget.slot_width() == 1.0 and widget.offset == 0
        assert widget.scroll_range() == 1500 - width

        # the mouse wheel scrolls, with Ctrl it changes the bars' width
        def wheel(amount, modifiers=QtCore.Qt.KeyboardModifier.NoModifier):
            point = QtCore.QPointF(plot.center().x(), 100)
            QtWidgets.QApplication.sendEvent(widget, QtGui.QWheelEvent(
                point, point, QtCore.QPoint(0, 0), QtCore.QPoint(0, amount),
                QtCore.Qt.MouseButton.NoButton, modifiers,
                QtCore.Qt.ScrollPhase.NoScrollPhase, False))
        wheel(-120)
        assert widget.offset == 120
        wheel(120)
        assert widget.offset == 0
        wheel(120, QtCore.Qt.KeyboardModifier.ControlModifier)
        assert widget.slot_width() == pytest.approx(1.5)
        wheel(-120, QtCore.Qt.KeyboardModifier.ControlModifier)
        assert widget.slot_width() == pytest.approx(1.0)

        # the slider below the chart follows it, and moves it
        assert row.slider.maximum() == widget.scroll_range() == 1500 - width
        assert row.slider.pageStep() == int(plot.width())
        row.slider.setValue(300)
        assert widget.offset == 300
        widget.set_offset(50)
        assert row.slider.value() == 50
        assert row.wider.isEnabled() and row.narrower.isEnabled()
        row.wider.click()
        assert widget.slot_width() == pytest.approx(1.5)
        assert row.slider.maximum() == 2250 - width
        row.narrower.click()
        assert widget.slot_width() == pytest.approx(1.0)
        # the whole chain: the slider's handle reaches from end to end and
        # moves nothing
        row.whole.click()
        assert widget.scroll_range() == 0
        assert row.slider.maximum() == 1
        assert row.slider.pageStep() == int(widget.plot_rect().width())
        row.slider.setValue(1)
        assert widget.offset == 0
        assert not row.whole.isEnabled() and not row.narrower.isEnabled()
        assert row.wider.isEnabled()
        # resizing keeps the view inside the chain
        row.wider.click()
        widget.set_offset(10**6)
        assert widget.offset > 0
        old = widget.size()
        widget.resize(900, 220)
        # (a widget that isn't shown gets no resize event of its own)
        QtWidgets.QApplication.sendEvent(
            widget, QtGui.QResizeEvent(widget.size(), old))
        assert widget.offset <= widget.scroll_range()
        assert row.slider.maximum() == max(1, widget.scroll_range())
        assert row.slider.value() == widget.offset
        # nothing shown: nothing to move or to widen
        widget.set_profile([])
        assert row.slider.maximum() == 1 and widget.scroll_range() == 0
        assert not row.wider.isEnabled() and not row.whole.isEnabled()
    finally:
        shiboken6.delete(row)
        shiboken6.delete(widget)


def test_dialog_profile(session):
    import shiboken6
    chart = _chart_module()
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    Qt = pytest.importorskip('pymol.Qt').QtCore.Qt
    from pmg_tk.startup.interface_gui.main import InterfaceDialog

    dialog = InterfaceDialog(None, _self=cmd)
    try:
        assert isinstance(dialog.chart, chart.ProximityChart)
        assert dialog.chart_scroll.chart is dialog.chart
        assert dialog.chart.message == 'Analyze an object to see its profile'
        assert not dialog.copy_profile_button.isEnabled()
        assert not dialog.save_profile_button.isEnabled()
        assert not dialog.color_profile_button.isEnabled()
        dialog.refresh_objects()
        dialog.analyze()
        chains = sorted(cmd.get_chains('m and polymer'))
        assert dialog.chart.rows == []                 # no target chosen
        assert dialog.chart.message == 'Choose the target chain(s)'
        items = [dialog.target_list.item(i)
                 for i in range(dialog.target_list.count())]
        items[0].setCheckState(Qt.CheckState.Checked)
        expected = interfaces.proximity(interfaces.residue_atoms('m'),
                                        [chains[0]])
        assert dialog.proximity == expected
        assert [dialog.profile_chain.itemText(i)
                for i in range(dialog.profile_chain.count())] == [chains[0]]
        assert dialog.chart.rows == expected[chains[0]]
        assert dialog.chart.chain == chains[0]
        assert dialog.chart.caption.startswith(
            'Proximity of m chain %s to chains %s' % (
                chains[0], '+'.join(chains[1:])))
        assert dialog.copy_profile_button.isEnabled()
        assert dialog.save_profile_button.isEnabled()
        within = sum(1 for r in expected[chains[0]] if r[2] <= 4.0)
        assert within > 0
        assert dialog.profile_note.text().startswith(
            '%d of %d residues within 4.0 Å' % (
                within, len(expected[chains[0]])))
        # the slider under the graph: this chain fits, wider bars don't
        # (a dialog that isn't shown has no size yet: give the graph one)
        row = dialog.chart_scroll
        dialog.chart.resize(800, 240)
        row.follow()
        assert dialog.chart.scroll_range() == 0 and row.wider.isEnabled()
        assert not row.whole.isEnabled()
        while not dialog.chart.scroll_range():
            row.wider.click()
        assert row.whole.isEnabled()
        assert row.slider.maximum() == dialog.chart.scroll_range()
        row.slider.setValue(row.slider.maximum())
        assert dialog.chart.shown()[-1] == len(expected[chains[0]]) - 1
        assert dialog.chart.shown()[0] > 0
        # the line: the view stays where it is
        view = (dialog.chart.slot_width(), dialog.chart.offset)
        dialog.contact_distance.setValue(8.0)
        assert (dialog.chart.slot_width(), dialog.chart.offset) == view
        assert dialog.chart.cutoff == 8.0
        more = sum(1 for r in expected[chains[0]] if r[2] <= 8.0)
        assert more > within
        assert dialog.profile_note.text().startswith(
            '%d of %d residues within 8.0 Å' % (
                more, len(expected[chains[0]])))
        # a click on a bar zooms in on its residue
        view = cmd.get_view()
        dialog.chart.residue_clicked.emit(chains[0],
                                          expected[chains[0]][5][0])
        assert cmd.get_view() != view
        # the table for a spreadsheet
        dialog.copy_profile()
        copied = QtWidgets.QApplication.clipboard().text()
        assert copied == interfaces.format_proximity(expected)
        assert len(copied.split('\n')) == len(expected[chains[0]]) + 1
        assert dialog.status.text().startswith('Copied the proximity of')
        # the structure in the graph's colors
        assert dialog.color_profile_button.isEnabled()
        dialog.color_profile_button.click()
        assert dialog.status.text() == (
            'Colored %d residues of chain %s by their proximity to the '
            'other chains' % (len(expected[chains[0]]), chains[0]))
        near = min(expected[chains[0]], key=lambda r: r[2])
        found = set()
        cmd.iterate('m and chain %s and resi %s' % (chains[0], near[0]),
                    'found.add(color)', space={'found': found})
        assert found == {cmd.get_color_index('proximity_%03d' % round(
            100 * interfaces.proximity_fraction(near[2])))}
        assert cmd.get_color_tuple('proximity_100') == pytest.approx(
            chart.ramp_color(1).getRgbF()[:3], abs=1e-6)
        # a second target: its own profile, chosen in the chain list
        items[1].setCheckState(Qt.CheckState.Checked)
        assert [dialog.profile_chain.itemText(i)
                for i in range(dialog.profile_chain.count())] == chains[:2]
        assert dialog.profile_chain.currentText() == chains[0]  # kept
        dialog.profile_chain.setCurrentIndex(1)
        assert dialog.chart.chain == chains[1]
        assert dialog.chart.rows == dialog.proximity[chains[1]]
        assert dialog.chart.caption.startswith(
            'Proximity of m chain %s to chains %s' % (
                chains[1], '+'.join(chains[2:])))
        # no targets left
        for item in items[:2]:
            item.setCheckState(Qt.CheckState.Unchecked)
        assert dialog.chart.rows == [] and dialog.proximity == {}
        assert dialog.profile_chain.count() == 0
        assert not dialog.copy_profile_button.isEnabled()
        assert dialog.profile_note.text() == ''
        # every chain a target: nothing to measure to
        for item in items:
            item.setCheckState(Qt.CheckState.Checked)
        assert dialog.chart.message == 'No other chains to measure to'
    finally:
        shiboken6.delete(dialog)
