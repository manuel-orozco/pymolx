'''
"Scripts" in the menu bar: File > Run Script..., starting in the scripts
folder of the source checkout pymolx was built from.
'''

import os

import pytest

from pymolx import paths

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..',
                                    '..'))


def test_scripts_dir_is_in_the_source_checkout():
    # setup.py records where pymolx was built from
    assert paths.source_dir() == REPO
    assert paths.scripts_dir() == os.path.join(REPO, 'scripts')
    assert os.path.isfile(os.path.join(paths.scripts_dir(),
                                       'save_all_pdb.py'))


def test_scripts_dir_without_checkout(monkeypatch, tmp_path):
    monkeypatch.setattr(paths, 'source_dir', lambda: None)
    monkeypatch.setattr(paths, 'USER_SCRIPTS_DIR', str(tmp_path / 'scripts'))
    assert paths.scripts_dir() == str(tmp_path / 'scripts')
    assert not os.path.exists(paths.scripts_dir())
    paths.scripts_dir(create=True)
    assert os.path.isdir(str(tmp_path / 'scripts'))


try:
    from pymol.Qt import QtCore, QtGui, QtWidgets
except ImportError:
    QtWidgets = None


@pytest.mark.skipif(QtWidgets is None, reason='no Qt')
def test_menu_bar_item_runs_a_script(monkeypatch, tmp_path):
    # same as File > Run Script... (window.file_run), starting in scripts/
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(['t'])
    from pymolx.gui import scripts

    folder = tmp_path / 'scripts'
    other = tmp_path / 'elsewhere'

    def scripts_dir(create=False):
        if create:
            folder.mkdir(exist_ok=True)
        return str(folder)

    monkeypatch.setattr(scripts.paths, 'scripts_dir', scripts_dir)

    class Window(QtWidgets.QMainWindow):
        initialdir = '/home/user/structures'
        chosen = None
        calls = []

        def file_run(self):
            # like pmg_qt: the dialog opens in initialdir; a chosen script
            # moves initialdir to its folder
            self.calls.append(self.initialdir)
            if self.chosen:
                self.initialdir = os.path.dirname(self.chosen)

    win = Window()
    menubar = win.menuBar()
    win.menudict = {}
    for title in ('File', 'Plugin', 'Help'):
        menu = QtWidgets.QMenu(title, menubar)
        menubar.addMenu(menu)
        win.menudict[title] = menu
    action = scripts.setup(win)
    titles = [a.text() for a in menubar.actions()]
    assert titles == ['File', 'Plugin', 'Scripts', 'Help']  # before Help
    assert action.menu() is None  # a plain item

    # cancelled: the dialogs' folder is left as it was
    action.trigger()
    assert win.calls == [str(folder)] and folder.is_dir()
    assert win.initialdir == '/home/user/structures'

    # a script chosen elsewhere: as File > Run Script... leaves it
    win.chosen = str(other / 'x.pml')
    action.trigger()
    assert win.calls[-1] == str(folder)
    assert win.initialdir == str(other)
    win.deleteLater()


def test_real_window_has_file_run():
    # the menu item calls what File > Run Script... calls
    from pmg_qt.pymol_qt_gui import PyMOLQtGUI
    assert callable(PyMOLQtGUI.file_run)
