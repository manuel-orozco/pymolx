'''
Hover highlighting (like Maestro): the unit under the mouse pointer drawn
with lighter colors, same style. Picking itself needs an on-screen OpenGL
view (checked in the GUI); here: what gets highlighted, and how it looks
in a ray traced image (same line and stick geometry).
'''

import os

import pytest

from pymol import cmd
from pymolx import hover, sequence, undo

DATA = os.path.join(os.path.dirname(__file__), '..', '..', 'data')


@pytest.fixture(autouse=True)
def clean_session():
    cmd.reinitialize()
    cmd.load(os.path.join(DATA, '1rx1.pdb'), 'dhfr')
    cmd.show_as('lines', 'polymer')
    yield
    hover.clear()
    cmd.reinitialize()


def ca(resi):
    found = []
    cmd.iterate('resi %s and name CA' % resi, 'found.append(index)',
                space={'found': found})
    return found[0]


def atoms(selection):
    found = set()
    cmd.iterate(selection, 'found.add(index)', space={'found': found})
    return found


def first_atom(selection):
    return min(atoms(selection))


def test_unit_follows_mouse_selection_mode():
    i = ca(50)
    assert atoms(hover.unit_selection('dhfr', i, mode=0)) == {i}
    assert atoms(hover.unit_selection('dhfr', i, mode=1)) == atoms('resi 50')
    assert atoms(hover.unit_selection('dhfr', i, mode=2)) == \
        atoms('chain A')


def test_unit_key():
    # residue mode: any atom of the residue, same unit
    a, b = ca(50), first_atom('resi 50 and not name CA')
    assert hover.unit_key('dhfr', a, 1) == hover.unit_key('dhfr', b, 1)
    assert hover.unit_key('dhfr', a, 1) != hover.unit_key('dhfr', ca(51), 1)
    # atom mode: per atom
    assert hover.unit_key('dhfr', a, 0) != hover.unit_key('dhfr', b, 0)
    # chain mode: the whole chain
    assert hover.unit_key('dhfr', a, 2) == hover.unit_key('dhfr', ca(90), 2)


def brightness(tmp_path, name):
    from PIL import Image, ImageStat
    path = str(tmp_path / name)
    cmd.png(path, width=200, height=150, ray=1)
    image = Image.open(path).convert('RGBA')
    return ImageStat.Stat(image.convert('L'), image.getchannel('A')).mean[0]


@pytest.mark.parametrize('rep', ['lines', 'sticks'])
def test_lighter_colors_same_style(tmp_path, rep):
    cmd.hide('everything')
    cmd.show(rep, 'resi 48-52')
    cmd.orient('resi 48-52')
    # sub-pixel lines (e.g. 0.75 px) don't ray trace pixel-exactly
    # from run to run
    cmd.set('line_width', 1.49)
    colors = []
    cmd.iterate('resi 50', 'colors.append(color)', space={'colors': colors})

    plain = brightness(tmp_path, 'plain.png')
    assert hover.show('resi 50')
    lighter = brightness(tmp_path, 'hover.png')
    hover.clear()
    assert lighter > plain + 3
    assert brightness(tmp_path, 'clear.png') == pytest.approx(plain, abs=0.5)
    # the atoms' colors themselves never change
    after = []
    cmd.iterate('resi 50', 'after.append(color)', space={'after': after})
    assert after == colors


def test_highlight_stays_out_of_the_way():
    names = cmd.get_names('all')
    n_all = cmd.count_atoms('all')
    cmd.undo_enable()
    try:
        steps = len(undo.stack.undo_stack)
        counts = sequence.get_change_counts()
        assert hover.show('resi 50')
        assert hover.show('resi 51')  # replaced
        hover.clear()
        # nothing for the sequence viewer to update
        assert sequence.get_change_counts() == counts
        assert len(undo.stack.undo_stack) == steps  # no undo steps
    finally:
        cmd.undo_disable()
    # no objects or atoms added, ever
    assert hover.show('resi 50')
    assert cmd.get_names('all') == names
    assert cmd.count_atoms('all') == n_all


def test_deleted_object_while_hovered():
    hover.show('resi 50')
    cmd.delete('dhfr')
    hover.clear()  # no crash
    cmd.load(os.path.join(DATA, '1rx1.pdb'), 'dhfr')
    assert hover.show('resi 50')


def test_nothing_to_highlight():
    assert not hover.show('none')


def test_uncounted_nests():
    counts = sequence.get_change_counts()
    with sequence.uncounted():
        with sequence.uncounted():
            cmd.color('red', 'resi 1')
    assert sequence.get_change_counts() == counts
    cmd.color('blue', 'resi 1')  # outside: counted
    assert sequence.get_change_counts()[0] > counts[0]


# the Qt part, with a stand-in for picking

try:
    from pymol.Qt import QtCore, QtGui, QtWidgets
except ImportError:
    QtWidgets = None


@pytest.mark.skipif(QtWidgets is None, reason='no Qt')
def test_highlighter_events(monkeypatch):
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(['t'])
    from pymolx.gui import hover as ghover
    from pymolx.gui.sequence_viewer import MemoryPrefs

    picks = {'at': ('dhfr', ca(50))}

    class Pick:
        @staticmethod
        def pick_at(_cob, x, y):
            return picks['at']

    monkeypatch.setattr(ghover, '_cmd', Pick)
    shown = []
    monkeypatch.setattr(ghover.hover, 'show',
                        lambda sel, _self=cmd: shown.append(sel) or True)
    monkeypatch.setattr(ghover.hover, 'clear',
                        lambda _self=cmd: shown.append(None))

    class GLStandIn(QtWidgets.QWidget):
        fb_scale = 1.0

        def __init__(self):
            super().__init__()
            self.before_draw = []
            self.updates = 0

        def update(self):  # (no paintGL: tests run it with paint())
            self.updates += 1

        def paint(self):
            for func in self.before_draw:
                func()

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    win = QtWidgets.QMainWindow()
    win.cmd = cmd
    win.pymolwidget = GLStandIn()
    display = QtWidgets.QMenu('Display', win)
    display.addAction('Orthoscopic View')
    win.menudict = {'Display': display}
    prefs = MemoryPrefs()
    hl = ghover.setup(win, prefs=prefs)
    assert [a.text() for a in display.actions()] == [
        'Hover Highlight', 'Orthoscopic View']

    def event(kind, buttons=QtCore.Qt.MouseButton.NoButton):
        Type = QtCore.QEvent.Type
        if kind == 'leave':
            ev = QtCore.QEvent(Type.Leave)
        else:
            t = {'move': Type.MouseMove, 'press': Type.MouseButtonPress,
                 'release': Type.MouseButtonRelease}[kind]
            p = QtCore.QPointF(10, 10)
            ev = QtGui.QMouseEvent(t, p, p, QtCore.Qt.MouseButton.LeftButton,
                                   buttons, QtCore.Qt.KeyboardModifier.NoModifier)
        QtWidgets.QApplication.sendEvent(win.pymolwidget, ev)

    cmd.set('mouse_selection_mode', 1)
    event('move')
    hl.update()  # (the timer, without waiting)
    # picking draws into the frame buffer: only in paintGL, never here
    assert shown == [] and win.pymolwidget.updates == 1
    win.pymolwidget.paint()
    assert shown == ['byresi(%%dhfr and index %d)' % ca(50)]
    win.pymolwidget.paint()  # nothing pending
    assert len(shown) == 1

    # another atom of the same residue: nothing new
    picks['at'] = ('dhfr', first_atom('resi 50 and not name CA'))
    event('move')
    hl.update()
    win.pymolwidget.paint()
    assert len(shown) == 1

    event('press', QtCore.Qt.MouseButton.LeftButton)  # rotating: out of the way
    assert shown[-1] is None
    event('release')
    hl.update()
    win.pymolwidget.paint()
    assert shown[-1].startswith('byresi')

    picks['at'] = None  # empty space
    event('move')
    hl.update()
    win.pymolwidget.paint()
    assert shown[-1] is None

    picks['at'] = ('dhfr', ca(50))
    event('move')
    hl.update()
    win.pymolwidget.paint()
    event('leave')
    assert shown[-1] is None

    # the pointer left before the paint: no pick
    count = len(shown)
    event('move')
    hl.update()
    event('leave')
    picks['at'] = ('dhfr', ca(60))
    win.pymolwidget.paint()
    assert all(s is None for s in shown[count:])

    # Display > Hover Highlight: off (kept in the preferences)
    hl.action.trigger()
    count = len(shown)
    event('move')
    hl.update()
    win.pymolwidget.paint()
    assert all(s is None for s in shown[count:])
    assert prefs[ghover.PREF_HOVER] is False
    win.deleteLater()
