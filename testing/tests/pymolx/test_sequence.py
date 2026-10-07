'''
pymolx sequence viewer (INCENTIVE_PARITY.md L-07): rows, colors and
selections from pymolx.sequence, and the Qt viewer (headless, Qt
offscreen platform).
'''

import os
from collections import defaultdict

import pytest

from pymol import cmd
from pymolx import sequence

DATA = os.path.join(os.path.dirname(__file__), '..', '..', 'data')


@pytest.fixture(autouse=True)
def clean_session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


def row_text(row):
    '''cells with their gaps, e.g. "ACDG NAP CA"'''
    text = ''
    for res in row.residues:
        text += ' ' * (res.col - len(text)) + res.text
    return text


def by_label(rows):
    return {row.label: row for row in rows}


# rows

def test_rows_per_chain():
    cmd.fab('ACDEFGHIKLMNPQ', 'p', chain='A')
    cmd.fab('WY', 'p2', chain='B', resi=101)
    cmd.create('p', 'p or p2')
    cmd.delete('p2')
    rows = sequence.get_rows()
    assert [row.label for row in rows] == ['p/A', 'p/B']
    assert row_text(rows[0]) == 'ACDEFGHIKLMNPQ'
    assert row_text(rows[1]) == 'WY'
    assert rows[0].ncols == 14
    assert all(row.has_polymer for row in rows)


def test_residue_numbers_and_ticks():
    cmd.fab('A' * 25, 'p')
    row, = sequence.get_rows()
    # first residue, then every 10th; ticks every 5th
    assert row.numbers == [(0, '1'), (9, '10'), (19, '20')]
    assert row.ticks == [4, 9, 14, 19, 24]


def test_numbers_dont_overlap():
    cmd.fab('AAAA', 'p', resi=998)  # 998..1001
    row, = sequence.get_rows()
    # "998" covers columns 0-2, so "1000" at column 2 is left out
    assert row.numbers == [(0, '998')]


def test_ligands_ions_and_water():
    cmd.load(os.path.join(DATA, '1rx1.pdb'), 'dhfr')
    row, = sequence.get_rows()
    assert row.label == 'dhfr/A'  # segment "A" is not repeated
    text = row_text(row)
    assert text.startswith('MISLIAALAVDRVIG')
    # non-polymer residues: names with a blank column before; no water
    assert text.endswith('ERR NAP CA')
    assert 'HOH' not in text
    nap = [r for r in row.residues if r.resn == 'NAP'][0]
    assert not nap.polymer and nap.width == 3
    assert row.residue_at(nap.col + 2) is nap
    assert row.residue_at(nap.col - 1) is None


def test_ligands_after_polymer():
    cmd.load(os.path.join(DATA, '1rx1.pdb'), 'dhfr')
    cmd.alter('resn NAP', 'resi="0"')  # before the protein in atom order
    cmd.sort()
    row, = sequence.get_rows()
    assert row_text(row).endswith('ERR NAP CA')


def two_chains():
    cmd.fab('ACD', 'p', chain='A')
    cmd.fab('EF', 'p2', chain='B', resi=10)
    cmd.create('p', 'p or p2')
    cmd.delete('p2')


def test_rows_by_object():
    two_chains()
    cmd.fab('GH', 'other')
    rows = sequence.get_rows(by_object=True)
    assert [row.label for row in rows] == ['p', 'other']
    row = rows[0]
    assert row.chain is None
    # like the OpenGL viewer: "/A/" markers above blank cells
    assert [(col, text) for col, text, _ in row.markers] == [
        (0, '/A/'), (8, '/B/')]
    assert row_text(row) == '    ACD     EF'
    # each chain numbered from its first residue
    assert row.numbers == [(4, '1'), (12, '10')]
    assert [r.chain for r in row.marker_at(9)] == ['B', 'B']
    assert row.marker_at(3) is None  # the blank after the marker
    assert row.residue_at(12).resn == 'GLU'


def test_rows_by_object_with_ligand_chain():
    cmd.load(os.path.join(DATA, '1rx1.pdb'), 'dhfr')
    cmd.alter('resn NAP', 'chain="B"')
    row, = sequence.get_rows(by_object=True)
    assert [text for _, text, _ in row.markers] == ['/A/', '/B/']
    assert row_text(row).endswith('ERR CA     NAP')


def test_nucleic_acids():
    cmd.fnab('ACGT', mode='DNA', name='dna')
    rows = by_label(sequence.get_rows())
    # terminal residues included, all one-letter codes
    assert row_text(rows['dna/A']) == 'ACGT'
    assert row_text(rows['dna/B']) == 'ACGT'


def test_residue_names_format():
    cmd.fab('ACD', 'p')
    cmd.set('seq_view_format', 1)
    row, = sequence.get_rows()
    assert row_text(row) == 'ALA CYS ASP'
    assert row.numbers[0] == (0, '1')


def test_objects_shown():
    cmd.fab('ACD', 'p1')
    cmd.fab('EFG', 'p2')
    cmd.fab('HIK', '_hidden')
    cmd.fab('LMN', 'p3')
    cmd.disable('p2')
    cmd.set('seq_view', 'off', 'p3')
    assert [row.model for row in sequence.get_rows()] == ['p1']
    cmd.pseudoatom('ps')
    cmd.enable('p2')
    rows = sequence.get_rows()
    assert [row.model for row in rows] == ['p1', 'p2', 'ps']
    assert rows[2].label == 'ps/PSDO/P'  # segment and chain
    assert not rows[2].has_polymer


def test_colors_like_opengl_viewer():
    cmd.fab('ACD', 'p')
    cmd.color('green')
    cmd.color('red', 'resi 2')
    cmd.color('blue', 'resi 3 and not name CA')  # CA is the guide atom
    colors = [res.color for res in sequence.get_rows()[0].residues]
    assert colors == [cmd.get_color_index(c) for c in ('green', 'red', 'green')]


def test_nearest_index():
    cmd.load(os.path.join(DATA, '1rx1.pdb'), 'dhfr')
    row, = sequence.get_rows()
    nap = [i for i, r in enumerate(row.residues) if r.resn == 'NAP'][0]
    assert row.nearest_index(-5) == 0
    assert row.nearest_index(row.residues[nap].col - 1) in (nap - 1, nap)
    assert row.nearest_index(10 ** 6) == len(row.residues) - 1


# change counters

def test_change_counts():
    cmd.fab('ACD', 'p')
    changed, dirty = sequence.get_change_counts()
    cmd.color('red', 'resi 2')
    assert sequence.get_change_counts()[0] > changed
    changed, dirty = sequence.get_change_counts()
    cmd.select('s1', 'resi 1')
    assert sequence.get_change_counts()[1] > dirty


# selections

def keys(*resis, model='p', chain=''):
    return {(model, '', chain, str(r)) for r in resis}


def residues(row, *resis):
    return [r for r in row.residues if r.resi in map(str, resis)]


def test_active_selection():
    assert sequence.active_selection() is None
    assert sequence.active_selection(create=True) == 'sele'
    cmd.fab('ACD', 'p')
    cmd.select('s1', 'resi 1')
    cmd.select('s2', 'resi 2')
    assert sequence.active_selection() == 's2'
    cmd.disable('s2')
    assert sequence.active_selection() in ('s1', None)
    cmd.set('auto_number_selections')
    cmd.delete('s*')
    assert sequence.active_selection(create=True) == 'sel01'


def test_select_residues():
    cmd.fab('ACDEF', 'p')
    row, = sequence.get_rows()
    name = sequence.select_residues(residues(row, 2))
    assert name == 'sele'
    assert sequence.get_selected() == keys(2)
    sequence.select_residues(residues(row, 4, 5))
    assert sequence.get_selected() == keys(2, 4, 5)
    sequence.select_residues(residues(row, 4), add=False)
    assert sequence.get_selected() == keys(2, 5)
    # temporary selections are gone
    assert cmd.get_names('public_selections') == ['sele']
    assert not [n for n in cmd.get_names('selections')
                if n.startswith('_seqview')]


def test_select_follows_mouse_selection_mode():
    cmd.fab('ACD', 'p', chain='A')
    cmd.fab('EF', 'p2', chain='B')
    cmd.create('p', 'p or p2')
    cmd.delete('p2')
    row = sequence.get_rows()[0]

    cmd.set('mouse_selection_mode', 2)  # chains
    sequence.select_residues(residues(row, 2))
    assert sequence.get_selected() == keys(1, 2, 3, chain='A')

    cmd.set('mouse_selection_mode', 0)  # atoms: still the residue's atoms
    cmd.delete('sele')
    sequence.select_residues(residues(row, 2))
    assert cmd.count_atoms('sele') == cmd.count_atoms('p and chain A and resi 2')


def test_drag_from_base():
    cmd.fab('ACDEFGH', 'p')
    row, = sequence.get_rows()
    cmd.select('sele', 'resi 7')
    base = sequence.save_base()
    sequence.select_residues(row.residues[0:4], base=base)
    assert sequence.get_selected() == keys(1, 2, 3, 4, 7)
    # dragging back: the range shrinks
    sequence.select_residues(row.residues[0:2], base=base)
    assert sequence.get_selected() == keys(1, 2, 7)
    sequence.clear_base()
    assert sequence.BASE_SELE not in cmd.get_names('all')


def test_deselect_all():
    cmd.fab('ACD', 'p')
    cmd.select('sele', 'resi 1')
    assert sequence.deselect_all() == 'sele'
    assert cmd.count_atoms('sele') == 0


def test_wizard_gets_selection():
    from pymol.wizard import Wizard

    class Picker(Wizard):
        picked = None

        def do_select(self, name):
            Picker.picked = (name, cmd.count_atoms(name))

    cmd.fab('ACD', 'p')
    cmd.set_wizard(Picker(_self=cmd))
    try:
        row, = sequence.get_rows()
        sequence.select_residues(residues(row, 3))
        assert Picker.picked == ('sele', cmd.count_atoms('p and resi 3'))
    finally:
        cmd.set_wizard()


# Qt viewer

try:
    from pymol.Qt import QtCore, QtGui, QtWidgets
except ImportError:
    QtWidgets = None

needs_qt = pytest.mark.skipif(QtWidgets is None, reason='no Qt')


@pytest.fixture(scope='module')
def app():
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    return (QtWidgets.QApplication.instance() or
            QtWidgets.QApplication(['pymolx-test']))


@pytest.fixture
def window(app):
    from pymolx.gui import content_panel, sequence_viewer

    win = QtWidgets.QMainWindow()
    win.cmd = cmd
    win.setting_callbacks = defaultdict(list)
    display = QtWidgets.QMenu('Display', win)
    win.menudict = {'Wizard': QtWidgets.QMenu('Wizard', win),
                    'Display': display}
    # like pmg_qt's SettingAction
    action = display.addAction('Sequence')
    action.setCheckable(True)
    action.triggered.connect(
        lambda: cmd.set('seq_view', int(action.isChecked())))
    win.setting_callbacks[cmd.setting._get_index('seq_view')].append(
        lambda v: action.setChecked(bool(v)))
    win.sequence_action = action

    win.browser = QtWidgets.QPlainTextEdit()
    win.lineedit = QtWidgets.QLineEdit()
    frame = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout(frame)
    layout.addWidget(win.browser)
    layout.addWidget(win.lineedit)
    win.ext_window = QtWidgets.QDockWidget(win)
    win.ext_window.setWidget(frame)
    win.addDockWidget(QtCore.Qt.DockWidgetArea.BottomDockWidgetArea,
                      win.ext_window)
    win.setCentralWidget(QtWidgets.QWidget())

    mode = QtWidgets.QMenu('Sequence Mode', display)
    display.addMenu(mode)
    mode.addAction('Residue Codes')
    win.sequence_mode_menu = mode

    win.content_panel = content_panel.setup(win)
    win.content_panel.timer.stop()
    win.prefs = sequence_viewer.MemoryPrefs()
    win.sequence_viewer = sequence_viewer.setup(win, prefs=win.prefs)
    win.sequence_viewer.timer.stop()
    win.resize(900, 600)
    win.show()
    yield win
    cmd.set('internal_gui', 1)
    cmd.set('internal_feedback', 1)
    win.close()
    win.deleteLater()


def fire_setting_callbacks(win):
    for setting in cmd.get_setting_updates() or ():
        for callback in win.setting_callbacks.get(setting, ()):
            callback(cmd.get_setting_tuple(setting)[1][0])


def process():
    # twice: layouts settle, then zero timers (dock height correction) run
    QtWidgets.QApplication.processEvents()
    QtWidgets.QApplication.processEvents()


def mouse(widget, kind, pos, button=None, modifiers=None):
    Qt = QtCore.Qt
    button = button or Qt.MouseButton.LeftButton
    modifiers = modifiers or Qt.KeyboardModifier.NoModifier
    types = {'press': QtCore.QEvent.Type.MouseButtonPress,
             'move': QtCore.QEvent.Type.MouseMove,
             'release': QtCore.QEvent.Type.MouseButtonRelease,
             'double': QtCore.QEvent.Type.MouseButtonDblClick}
    buttons = Qt.MouseButton.NoButton if kind == 'release' else button
    pos = QtCore.QPointF(pos)
    event = QtGui.QMouseEvent(types[kind], pos, widget.mapToGlobal(pos),
                              button, buttons, modifiers)
    QtWidgets.QApplication.sendEvent(widget, event)


def cell_center(viewer, row, res):
    return viewer.residue_rect(row, res).center()


@needs_qt
def test_seq_toggle_shows_qt_viewer(window):
    viewer = window.sequence_viewer
    dock = window.sequence_dock
    button = window.content_panel.toggles.seq_button
    cmd.fab('ACDEF', 'p')

    assert dock.isHidden()
    button.click()
    assert viewer.shown and not dock.isHidden()
    assert cmd.get_setting_int('seq_view') == 0  # no OpenGL bar
    assert window.sequence_action.isChecked()
    button.click()
    assert not viewer.shown and dock.isHidden()
    assert not window.sequence_action.isChecked()


@needs_qt
def test_seq_view_setting_and_menu(window):
    viewer = window.sequence_viewer
    button = window.content_panel.toggles.seq_button
    cmd.fab('ACDEF', 'p')

    cmd.set('seq_view', 1)
    fire_setting_callbacks(window)
    assert viewer.shown and button.isChecked()
    assert window.sequence_action.isChecked()
    assert cmd.get_setting_int('seq_view') == 0
    fire_setting_callbacks(window)  # the reset to 0 doesn't hide it
    assert viewer.shown

    cmd.set('seq_view', 0)
    fire_setting_callbacks(window)
    assert not viewer.shown and not button.isChecked()

    window.sequence_action.trigger()  # Display > Sequence
    assert viewer.shown and button.isChecked()
    assert cmd.get_setting_int('seq_view') == 0


@needs_qt
def test_hidden_without_polymers(window):
    viewer = window.sequence_viewer
    dock = window.sequence_dock
    viewer.set_shown(True)
    assert dock.isHidden()
    cmd.pseudoatom('ps')
    viewer.refresh()
    assert dock.isHidden()
    cmd.load(os.path.join(DATA, '1rx1.pdb'), 'p')
    viewer.refresh()
    assert not dock.isHidden()
    assert [row.label for row in viewer.rows] == ['ps/PSDO/P', 'p/A']
    cmd.delete('p')
    viewer.refresh()
    assert dock.isHidden()


@needs_qt
def test_refresh_follows_changes(window):
    viewer = window.sequence_viewer
    cmd.fab('ACD', 'p')
    viewer.set_shown(True)
    rows = viewer.rows
    counts = viewer._counts
    viewer.refresh()
    assert viewer.rows is rows  # nothing changed, no rebuild
    assert viewer._counts == counts  # its own queries don't count

    cmd.select('sele', 'resi 2')  # e.g. a click in the 3D viewer
    viewer.refresh()
    assert viewer.rows is rows
    assert viewer.selected == keys(2)

    cmd.color('red', 'resi 1')
    viewer.refresh()
    assert viewer.rows is not rows
    assert viewer.rows[0].residues[0].color == cmd.get_color_index('red')

    cmd.disable('p')
    viewer.refresh()
    assert viewer.rows == []


@needs_qt
def test_click_drag_shift_click(window):
    viewer = window.sequence_viewer
    vp = viewer.viewport()
    cmd.fab('ACDEFGHIKL', 'p')
    viewer.set_shown(True)
    process()
    row = viewer.rows[0]
    res = row.residues

    # click selects, click again deselects; not an undo step
    mouse(vp, 'press', cell_center(viewer, row, res[1]))
    mouse(vp, 'release', cell_center(viewer, row, res[1]))
    assert sequence.get_selected() == keys(2)
    assert viewer.selected == keys(2)
    mouse(vp, 'press', cell_center(viewer, row, res[1]))
    mouse(vp, 'release', cell_center(viewer, row, res[1]))
    assert sequence.get_selected() == set()

    # drag 3..6, back to 5
    mouse(vp, 'press', cell_center(viewer, row, res[2]))
    mouse(vp, 'move', cell_center(viewer, row, res[5]))
    assert sequence.get_selected() == keys(3, 4, 5, 6)
    mouse(vp, 'move', cell_center(viewer, row, res[4]))
    mouse(vp, 'release', cell_center(viewer, row, res[4]))
    assert sequence.get_selected() == keys(3, 4, 5)
    assert sequence.BASE_SELE not in cmd.get_names('all')

    # shift-click extends from the start of the last range
    mouse(vp, 'press', cell_center(viewer, row, res[7]),
          modifiers=QtCore.Qt.KeyboardModifier.ShiftModifier)
    mouse(vp, 'release', cell_center(viewer, row, res[7]))
    assert sequence.get_selected() == keys(3, 4, 5, 6, 7, 8)

    # double click on empty space clears
    empty = QtCore.QPoint(viewer.width() - 5,
                          cell_center(viewer, row, res[0]).y())
    mouse(vp, 'double', empty)
    assert sequence.get_selected() == set()


@needs_qt
def test_selection_clicks_are_not_undo_steps(window):
    from pymolx import undo
    viewer = window.sequence_viewer
    vp = viewer.viewport()
    cmd.fab('ACD', 'p')
    viewer.set_shown(True)
    process()
    row = viewer.rows[0]
    undo.stack.enable()
    try:
        before = len(undo.stack.undo_stack)
        mouse(vp, 'press', cell_center(viewer, row, row.residues[0]))
        mouse(vp, 'release', cell_center(viewer, row, row.residues[0]))
        assert sequence.get_selected() == keys(1)
        assert len(undo.stack.undo_stack) == before
    finally:
        undo.stack.disable()


@needs_qt
def test_label_click_selects_chain(window):
    viewer = window.sequence_viewer
    vp = viewer.viewport()
    two_chains()
    viewer.set_shown(True)
    process()
    row_b = viewer.rows[1]
    y = cell_center(viewer, row_b, row_b.residues[0]).y()
    mouse(vp, 'press', QtCore.QPoint(5, y))
    assert sequence.get_selected() == keys(10, 11, chain='B')
    mouse(vp, 'press', QtCore.QPoint(5, y))
    assert sequence.get_selected() == set()


@needs_qt
def test_residue_menu(window):
    viewer = window.sequence_viewer
    cmd.fab('ACD', 'p')
    viewer.set_shown(True)
    row = viewer.rows[0]

    menu = viewer._residue_menu(row.residues[1], QtCore.QPoint(0, 0))
    labels = [a.text() for a in menu.actions()]
    assert 'zoom' in labels and 'color' in labels
    assert cmd.count_atoms(sequence_viewer_menu_sele()) == \
        cmd.count_atoms('resi 2')
    menu.close()

    cmd.select('sele', 'resi 2')
    viewer.refresh()
    menu = viewer._residue_menu(row.residues[1], QtCore.QPoint(0, 0))
    labels = [a.text() for a in menu.actions()]
    assert 'actions' in labels  # the selection's menu
    menu.close()


def sequence_viewer_menu_sele():
    from pymolx.gui import sequence_viewer
    return sequence_viewer.MENU_SELE


@needs_qt
def test_paint_colors_and_highlight(window):
    viewer = window.sequence_viewer
    cmd.fab('ACDEF', 'p')
    cmd.color('white')
    cmd.color('red', 'resi 3')
    viewer.set_shown(True)
    process()
    row = viewer.rows[0]

    def cell_pixels(res):
        image = viewer.viewport().grab().toImage()
        rect = viewer.residue_rect(row, res)
        return [image.pixelColor(x, y)
                for x in range(rect.left(), rect.right() + 1)
                for y in range(rect.top(), rect.bottom() + 1)]

    def reddish(pixels):
        # glyphs are antialiased: count clearly red pixels
        return sum(1 for c in pixels if c.redF() > 0.6 and
                   c.greenF() < 0.3 and c.blueF() < 0.3)

    assert reddish(cell_pixels(row.residues[2])) > 0      # red letter
    assert reddish(cell_pixels(row.residues[1])) == 0
    letter = reddish(cell_pixels(row.residues[2]))
    cmd.select('sele', 'resi 3')
    viewer.refresh()
    pixels = cell_pixels(row.residues[2])
    # inverted: red box with a dark letter
    assert reddish(pixels) > max(letter, len(pixels) // 2)
    assert reddish(pixels) < len(pixels)


@needs_qt
def test_horizontal_scroll(window):
    viewer = window.sequence_viewer
    cmd.fab('ACDEFGHIKLMNPQRSTVWY' * 10, 'p')
    viewer.set_shown(True)
    process()
    hbar = viewer.horizontalScrollBar()
    assert hbar.maximum() > 0
    row = viewer.rows[0]
    last = row.residues[-1]
    viewer.scroll_to(row, last)
    rect = viewer.residue_rect(row, last)
    assert rect.left() >= viewer._label_width()
    assert rect.right() <= viewer.viewport().width()
    assert viewer.hit(rect.center())[1] is last


@needs_qt
def test_row_layout_by_object(window):
    from pymolx.gui.sequence_viewer import PREF_BY_OBJECT
    viewer = window.sequence_viewer
    vp = viewer.viewport()
    two_chains()
    viewer.set_shown(True)
    process()
    assert [row.label for row in viewer.rows] == ['p/A', 'p/B']

    # Display > Sequence Mode has the layouts
    by_object = window.sequence_layout_actions[True]
    assert by_object in window.sequence_mode_menu.actions()
    by_object.trigger()
    assert viewer.by_object and window.prefs[PREF_BY_OBJECT] is True
    assert [row.label for row in viewer.rows] == ['p']

    # clicking the "/B/" marker selects chain B
    row = viewer.rows[0]
    col, text, _ = row.markers[1]
    x = viewer._col_x(col) + viewer.char_w
    y = viewer._row_top(0) + viewer.line_h // 2
    mouse(vp, 'press', QtCore.QPoint(x, y))
    assert sequence.get_selected() == keys(10, 11, chain='B')

    # the label column's menu switches back, and the Display menu follows
    menu = viewer.layout_menu()
    by_chain = [a for a in menu.actions() if a.text() == 'One Row per Chain']
    by_chain[0].trigger()
    assert not viewer.by_object
    assert window.sequence_layout_actions[False].isChecked()
    assert [row.label for row in viewer.rows] == ['p/A', 'p/B']


@needs_qt
def test_label_width_drag(window):
    viewer = window.sequence_viewer
    vp = viewer.viewport()
    cmd.fab('ACDEF', 'p')
    viewer.set_shown(True)
    process()
    auto = viewer._label_width()
    y = viewer._row_top(0) + viewer.line_h

    mouse(vp, 'press', QtCore.QPoint(auto, y))
    mouse(vp, 'move', QtCore.QPoint(auto + 60, y))
    assert viewer._label_width() == auto + 60
    mouse(vp, 'release', QtCore.QPoint(auto + 60, y))
    assert viewer.label_width == auto + 60
    assert not [k for k in window.prefs if 'width' in k]  # this run only
    assert sequence.get_selected() == set()  # no residue clicked
    # cells move with the label column
    assert viewer.residue_rect(viewer.rows[0],
                               viewer.rows[0].residues[0]).left() > auto + 60

    # can't squeeze the cells away or the labels to nothing
    viewer.set_label_width(10 ** 6)
    assert viewer._label_width() < vp.width() - 5 * viewer.char_w
    viewer.set_label_width(1)
    assert viewer._label_width() >= 3 * viewer.char_w

    # double-click on the edge fits the labels again
    mouse(vp, 'double', QtCore.QPoint(viewer._label_width(), y))
    assert viewer.label_width is None
    assert viewer._label_width() == auto


@needs_qt
def test_layout_prefs_are_loaded(app):
    from pymolx.gui import sequence_viewer as sv
    # a width saved by an earlier version is not used: labels fit
    prefs = sv.MemoryPrefs({sv.PREF_BY_OBJECT: True,
                            'pymolx_seq_label_width': 150})
    viewer = sv.SequenceViewer(cmd, prefs=prefs)
    viewer.timer.stop()
    viewer.resize(800, 100)
    two_chains()
    viewer.set_shown(True)
    assert [row.label for row in viewer.rows] == ['p']
    assert viewer.label_width is None
    viewer.deleteLater()


def wheel(widget, dy, modifiers=None):
    Qt = QtCore.Qt
    pos = QtCore.QPointF(widget.width() / 2, widget.height() / 2)
    event = QtGui.QWheelEvent(
        pos, widget.mapToGlobal(pos), QtCore.QPoint(0, 0),
        QtCore.QPoint(0, dy), Qt.MouseButton.NoButton,
        modifiers or Qt.KeyboardModifier.NoModifier,
        Qt.ScrollPhase.NoScrollPhase, False)
    QtWidgets.QApplication.sendEvent(widget, event)


@needs_qt
def test_font_size_ctrl_wheel(window):
    from pymolx.gui import sequence_viewer as sv
    from pymolx.gui import theme
    ctrl = QtCore.Qt.KeyboardModifier.ControlModifier
    viewer = window.sequence_viewer
    vp = viewer.viewport()
    cmd.fab('ACDEFGHIKLMNPQRSTVWY' * 10, 'p')
    viewer.set_shown(True)
    process()
    size = viewer.font_size
    assert size == theme.FONT_SIZE
    row_h, char_w, height = viewer.row_h, viewer.char_w, viewer.height()

    wheel(vp, 240, ctrl)  # two notches up
    process()
    assert viewer.font_size == size + 2
    assert not [k for k in window.prefs if 'font' in k]  # this run only
    assert viewer.row_h > row_h and viewer.char_w > char_w
    assert viewer.height() > height  # the panel grows with the rows

    wheel(vp, -120, ctrl)
    assert viewer.font_size == size + 1
    for _ in range(4):  # touchpad: quarter notches add up
        wheel(vp, -30, ctrl)
    assert viewer.font_size == size

    wheel(vp, -120 * 50, ctrl)
    assert viewer.font_size == sv.FONT_SIZES[0]
    wheel(vp, 120 * 50, ctrl)
    assert viewer.font_size == sv.FONT_SIZES[1]

    wheel(vp, 120)  # no Ctrl: scrolls, font unchanged
    assert viewer.font_size == sv.FONT_SIZES[1]


@needs_qt
def test_font_size_keeps_position_and_label_width(window):
    viewer = window.sequence_viewer
    cmd.fab('ACDEFGHIKLMNPQRSTVWY' * 10, 'p')
    viewer.set_shown(True)
    process()
    row = viewer.rows[0]
    viewer.set_font_size(8)  # doubled below, within FONT_SIZES
    viewer.set_label_width(120)
    hbar = viewer.horizontalScrollBar()
    hbar.setValue(100 * viewer.char_w)
    left_col = hbar.value() // viewer.char_w

    char_w = viewer.char_w
    viewer.set_font_size(viewer.font_size * 2)
    assert abs(hbar.value() // viewer.char_w - left_col) <= 1
    # the dragged width follows the (whole pixel) character width
    assert viewer.label_width == pytest.approx(
        120 * viewer.char_w / char_w, abs=1)
    assert viewer.char_w > char_w
    assert viewer.residue_rect(row, row.residues[0]).height() == \
        viewer.line_h


@needs_qt
def test_font_size_menu(window):
    from pymolx.gui import theme
    viewer = window.sequence_viewer
    cmd.fab('ACD', 'p')
    viewer.set_shown(True)

    def actions():
        return {a.text(): a for a in viewer.layout_menu().actions()}

    default = 'Default Font Size (%d pt)' % theme.FONT_SIZE
    assert not actions()[default].isEnabled()
    actions()['Larger Font'].trigger()
    actions()['Larger Font'].trigger()
    assert viewer.font_size == theme.FONT_SIZE + 2
    actions()['Smaller Font'].trigger()
    assert viewer.font_size == theme.FONT_SIZE + 1
    actions()[default].trigger()
    assert viewer.font_size == theme.FONT_SIZE


@needs_qt
def test_font_size_matches_gui(app):
    # saved sizes from earlier versions are ignored: GUI's size, and the
    # monospace letters as tall as the interface font's
    from pymolx.gui import sequence_viewer as sv
    prefs = sv.MemoryPrefs({'pymolx_seq_font_size': 9})
    viewer = sv.SequenceViewer(cmd, prefs=prefs)
    viewer.timer.stop()
    assert viewer.font_size == sv.theme.FONT_SIZE
    ui = QtGui.QFont(QtWidgets.QApplication.font())
    ui.setPointSizeF(sv.theme.FONT_SIZE)
    cap = QtGui.QFontMetricsF(viewer.font()).capHeight()
    assert cap == pytest.approx(QtGui.QFontMetricsF(ui).capHeight(), abs=1)
    viewer.deleteLater()


def chains(n, name='p'):
    letters = 'ABCDEFGH'[:n]
    for i, chain in enumerate(letters):
        cmd.fab('ACDEFGHIKL', '_part%d' % i, chain=chain)
    cmd.create(name, ' '.join('_part%d' % i for i in range(n)))
    cmd.delete('_part*')


@needs_qt
def test_label_column_fits_long_names(window):
    viewer = window.sequence_viewer
    chains(2, 'migg1_fc_g0f_proteinA')
    viewer.set_shown(True)
    process()
    fm = QtGui.QFontMetrics(viewer.font())
    label_w = viewer._label_width()
    for row in viewer.rows:
        assert fm.horizontalAdvance(row.label) <= label_w - 2 * 4

    # absurdly long names: the cells keep most of the width
    cmd.set_name('migg1_fc_g0f_proteinA', 'x' * 150)
    viewer.refresh()
    assert viewer._label_width() <= viewer.viewport().width() * 0.5 + 1


@needs_qt
def test_default_height_three_rows(window):
    viewer = window.sequence_viewer
    chains(2)
    viewer.set_shown(True)
    process()
    assert viewer.height() == viewer.rows_height(2)

    chains(5)
    viewer.refresh()
    process()
    assert len(viewer.rows) == 5
    assert viewer.height() == viewer.rows_height(3)
    assert viewer.verticalScrollBar().maximum() > 0
    lo, hi = viewer.height_limits()
    assert viewer.minimumHeight() == lo == viewer.rows_height(1)
    assert hi == viewer.rows_height(5)
    # at most a dock separator's width of slack
    separator = window.style().pixelMetric(
        QtWidgets.QStyle.PixelMetric.PM_DockWidgetSeparatorExtent, None,
        window)
    assert viewer.maximumHeight() - hi == separator


def separator_drag(win, dock, dy):
    '''drag the edge between the dock and the central widget'''
    Qt = QtCore.Qt
    x = dock.geometry().center().x()
    y = (dock.geometry().bottom() + win.centralWidget().geometry().top()) // 2
    for kind, pos, buttons in [
            (QtCore.QEvent.Type.MouseButtonPress, y, Qt.MouseButton.LeftButton),
            (QtCore.QEvent.Type.MouseMove, y + dy // 2, Qt.MouseButton.LeftButton),
            (QtCore.QEvent.Type.MouseMove, y + dy, Qt.MouseButton.LeftButton),
            (QtCore.QEvent.Type.MouseButtonRelease, y + dy, Qt.MouseButton.NoButton)]:
        p = QtCore.QPointF(x, pos)
        event = QtGui.QMouseEvent(kind, p, win.mapToGlobal(p),
                                  Qt.MouseButton.LeftButton, buttons,
                                  Qt.KeyboardModifier.NoModifier)
        QtWidgets.QApplication.sendEvent(win, event)
        process()


@needs_qt
def test_user_resizes_height(window):
    viewer = window.sequence_viewer
    chains(5)
    viewer.set_shown(True)
    process()
    three = viewer.rows_height(3)
    assert viewer.height() == three

    separator_drag(window, window.sequence_dock, viewer.row_h)
    assert viewer.height() > three
    assert viewer.user_height == viewer.height()
    user = viewer.height()

    # kept while the rows change; never more than all rows
    cmd.delete('p')
    chains(4)
    viewer.refresh()
    process()
    assert viewer.height() == min(user, viewer.rows_height(4))

    # back to the default
    menu = {a.text(): a for a in viewer.layout_menu().actions()}
    menu['Default Height (3 Rows)'].trigger()
    process()
    assert viewer.user_height is None
    assert viewer.height() == viewer.rows_height(3)


@needs_qt
def test_prefs_read_after_window_setup(window):
    # PyMOL reads ~/.pymolpluginsrc.py after the window is built: saved
    # values that show up later still apply, without being saved again
    from pymolx.gui import sequence_viewer as sv
    viewer = window.sequence_viewer
    assert not viewer.by_object
    window.prefs.update({sv.PREF_BY_OBJECT: True})
    process()  # setup's deferred load_prefs has run already
    viewer.load_prefs()
    assert viewer.by_object
    assert window.sequence_layout_actions[True].isChecked()

    class Watch(sv.MemoryPrefs):
        def set(self, key, value):
            raise AssertionError('saved %s' % key)

    viewer.prefs = Watch(window.prefs)
    viewer.load_prefs()  # loading again doesn't save anything
