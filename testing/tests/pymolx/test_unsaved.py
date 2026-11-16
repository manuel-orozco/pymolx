'''
pymolx.unsaved: has the session changed since it was saved or opened
(asked about before quitting, pymolx.gui.unsaved)
'''
import os

import pytest
from pymol import cmd

from pymolx import unsaved


@pytest.fixture(autouse=True)
def fresh():
    cmd.reinitialize()
    unsaved._baseline = None
    yield
    cmd.reinitialize()
    unsaved._baseline = None


def test_nothing_to_save():
    assert not unsaved.has_unsaved()


def test_never_saved():
    cmd.fab('ACDEF', 'pep')
    assert unsaved.has_unsaved()


def test_fingerprint_is_stable():
    cmd.fab('ACDEF', 'pep')
    # (the first cmd.get_session sets pse_binary_dump: left out)
    assert unsaved.fingerprint() == unsaved.fingerprint()


def test_saved_then_changed(tmp_path):
    cmd.fab('ACDEF', 'pep')
    cmd.save(str(tmp_path / 'one.pse'))
    assert not unsaved.has_unsaved()
    # looking around and GUI state aren't changes
    cmd.turn('y', 40)
    cmd.set('seq_view', 1)
    cmd.set('mouse_selection_mode', 2)
    assert not unsaved.has_unsaved()
    # what the session holds is
    cmd.color('red', 'pep and resi 2')
    assert unsaved.has_unsaved()
    cmd.save(str(tmp_path / 'one.pse'))
    assert not unsaved.has_unsaved()
    cmd.translate([1, 0, 0], 'pep and resi 3', camera=0)     # (a mouse edit)
    assert unsaved.has_unsaved()
    cmd.save(str(tmp_path / 'one.pse'))
    cmd.set('cartoon_transparency', 0.4)
    assert unsaved.has_unsaved()
    cmd.save(str(tmp_path / 'one.pse'))
    cmd.delete('all')
    assert not unsaved.has_unsaved()            # (nothing to save)


def test_other_formats_and_partial(tmp_path):
    cmd.fab('ACDEF', 'pep')
    cmd.save(str(tmp_path / 'pep.pdb'))         # a structure, not the session
    assert unsaved.has_unsaved()
    cmd.save(str(tmp_path / 'part.pse'), 'pep', partial=1)
    assert unsaved.has_unsaved()
    cmd.save(str(tmp_path / 'whole.pze'))
    assert not unsaved.has_unsaved()


def test_opened_session(tmp_path):
    cmd.fab('ACDEF', 'pep')
    cmd.save(str(tmp_path / 'one.pse'))
    cmd.reinitialize()
    unsaved._baseline = None
    cmd.load(str(tmp_path / 'one.pse'))
    assert not unsaved.has_unsaved()
    # merged into another session: that one isn't saved
    cmd.fab('GG', 'other')
    cmd.save(str(tmp_path / 'two.pse'))
    cmd.load(str(tmp_path / 'one.pse'), partial=1)
    assert unsaved.has_unsaved()
    # from the command language too
    cmd.do('save %s' % (tmp_path / 'three.pse'))
    cmd.sync()
    assert os.path.exists(str(tmp_path / 'three.pse'))
    assert not unsaved.has_unsaved()


def test_question(tmp_path):
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import shiboken6
    from pymolx.gui import unsaved as gui
    B = QtWidgets.QMessageBox.StandardButton

    class Window:
        # what gui.answer uses of the main window: Save, as File > Save
        # Session does (to session_file, or asking where: here, nowhere)
        def __init__(self, path):
            self.path = path

        def session_save(self):
            if self.path:
                cmd.save(self.path)

    cmd.fab('ACDEF', 'pep')
    box = gui.question(None, cmd)
    try:
        assert 'has not been saved' in box.informativeText()
        assert box.button(B.Discard).text() == "Don't Save"
        assert box.defaultButton() is box.button(B.Save)
        assert box.escapeButton() is box.button(B.Cancel)
    finally:
        shiboken6.delete(box)
    assert gui.answer(Window(''), cmd, B.Discard)       # quit, unsaved
    assert not gui.answer(Window(''), cmd, B.Cancel)    # stay
    assert not gui.answer(Window(''), cmd, B.Save)      # not saved: stay
    path = str(tmp_path / 'kept.pse')
    assert gui.answer(Window(path), cmd, B.Save)        # saved: quit
    assert os.path.exists(path) and not unsaved.has_unsaved()
    # nothing unsaved: no question
    assert gui.confirm(None, cmd)
    cmd.color('red', 'pep')
    box = gui.question(None, cmd)
    try:
        assert 'since kept.pse was saved' in box.informativeText()
    finally:
        shiboken6.delete(box)


def test_quit_asking(monkeypatch):
    # menus quit through the GUI's handler (which asks first) if there is
    # one; without a GUI it is cmd.quit
    called = []
    monkeypatch.setattr(unsaved, '_quit_handler', lambda: called.append(1))
    unsaved.quit_asking(cmd)
    assert called == [1]
    from pymol import menu
    entry = [e for e in menu.main_menu(cmd, (0.0, 0.0, 0.0), (0, 0))
             if 'quit' in str(e[1])][0]
    assert 'quit_asking' in entry[2]
