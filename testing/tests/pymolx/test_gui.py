'''
pymolx Qt GUI: theme, pymol.menu conversion and the main window toolbar
(INCENTIVE_PARITY.md: L-01, L-03). Runs headless on Qt's offscreen platform.
'''

import os
import re
from collections import defaultdict

import pytest

from pymol import cmd

try:
    from pymol.Qt import QtGui, QtWidgets
except ImportError:
    pytest.skip('no Qt', allow_module_level=True)

from pymolx.gui import menus, theme, toolbar


@pytest.fixture(scope='module')
def app():
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    return (QtWidgets.QApplication.instance() or
            QtWidgets.QApplication(['pymolx-test']))


class StubWindow(QtWidgets.QMainWindow):
    '''
    The parts of PyMOLQtGUI the toolbar uses.
    '''

    def __init__(self):
        super().__init__()
        self.cmd = cmd
        self.setting_callbacks = defaultdict(list)
        self.calls = []

    def fire_setting_callbacks(self):
        # like PyMOLQtGUI.update_feedback
        for setting in cmd.get_setting_updates() or ():
            for callback in self.setting_callbacks.get(setting, ()):
                callback(cmd.get_setting_tuple(setting)[1][0])

    def open_builder_panel(self):
        self.calls.append('builder')

    def scene_panel_menu_dialog(self):
        self.calls.append('scenes')

    def render_dialog(self, widget=None):
        self.calls.append('render')

    def get_view(self):
        self.calls.append('get_view')

    def open_props_dialog(self):
        self.calls.append('properties')


@pytest.fixture
def window(app):
    cmd.reinitialize()
    cmd.get_setting_updates()  # discard pending updates
    win = StubWindow()
    win.toolbar = toolbar.setup(win)
    yield win
    cmd.rock(0)
    cmd.reinitialize()
    win.deleteLater()


def toolbar_buttons(win):
    '''tool buttons in toolbar order'''
    tb = win.toolbar
    return [w for w in map(tb.widgetForAction, tb.actions())
            if isinstance(w, QtWidgets.QToolButton)]


def buttons(win):
    return {b.text(): b for b in toolbar_buttons(win) if b.text()}


def actions(menu):
    return {a.text(): a for a in menu.actions() if a.text()}


# theme

def test_stylesheet_tokens_substituted():
    style = theme.stylesheet(cmd)
    without_comments = re.sub(r'/\*.*?\*/', '', style, flags=re.S)
    assert '$' not in without_comments
    assert theme.COLORS['accent'] in style
    icons = cmd.exp_path(theme.ICON_DIR).replace('\\', '/')
    assert 'url(%s/chevron-down.svg)' % icons in style
    assert os.path.isfile(os.path.join(cmd.exp_path(theme.ICON_DIR),
                                       'chevron-down.svg'))


def test_palette(app):
    pal = theme.palette()
    assert pal.color(QtGui.QPalette.ColorRole.Window).name() == \
        theme.COLORS['window']
    assert pal.color(QtGui.QPalette.ColorRole.Highlight).name() == \
        theme.COLORS['selection']


def test_feedback_html():
    html = theme.feedback_html('PyMOL>zoom\n a < b')
    lines = html.split('<br>')
    assert lines[0].startswith('<span style="color:%s">' %
                               theme.COLORS['accent_text'])
    assert 'PyMOL&gt;zoom' in lines[0]
    assert lines[1] == '&nbsp;a&nbsp;&lt;&nbsp;b'


# pymol.menu conversion

def test_strip_label():
    assert menus.strip_label('\\900red ') == 'red'
    assert menus.strip_label('\\272Note:\\559 text') == 'Note: text'


def test_label_segments():
    assert menus.label_segments('plain') == [(None, 'plain')]
    assert menus.label_segments('\\900re\\090d ') == [
        ('#ff0000', 're'), ('#00ff00', 'd')]
    assert menus.label_segments('by \\559x') == [
        (None, 'by '), ('#8e8eff', 'x')]


def test_fill_menu(app):
    data = [
        [2, 'Title:', ''],
        [1, 'scale', 'set sphere_scale, 0.3'],
        [0, '', ''],
        [1, 'sub', [[2, 'Sub:', ''], [1, '\\090green', 'pass']]],
        [2, 'a note', ''],
    ]
    menu = menus.fill_menu(menus.PyMenu(), data, cmd)
    entries = [(a.text(), a.isSeparator(), a.isEnabled())
               for a in menu.actions()]
    assert entries == [
        ('Title:', False, False),
        ('scale', False, True),
        ('', True, True),
        ('sub', False, True),
        ('a note', False, False),
    ]
    assert menu.actions()[0].property('pymol_title')
    submenu = menu.actions()[3].menu()
    assert [a.text() for a in submenu.actions()] == ['Sub:', 'green']
    assert submenu.actions()[1].data() == '\\090green'

    cmd.set('sphere_scale', 1.0)
    menu.actions()[1].trigger()
    assert cmd.get_setting_float('sphere_scale') == pytest.approx(0.3)
    cmd.set('sphere_scale', 1.0)


def test_colored_labels_are_drawn(app):
    window = QtWidgets.QMainWindow()
    window.setStyleSheet(theme.stylesheet(cmd))
    menu = menus.fill_menu(menus.PyMenu('', window), [
        [2, 'Color:', ''],
        [1, '\\900reds', 'pass'],
    ], cmd)
    menu.ensurePolished()
    menu.adjustSize()
    image = menu.grab().toImage()

    def pixels(rect):
        return {image.pixelColor(x, y).getRgb()[:3]
                for x in range(rect.left(), rect.right())
                for y in range(rect.top(), rect.bottom())}

    reds = pixels(menu.actionGeometry(menu.actions()[1]))
    assert any(r > 150 and g < 80 and b < 80 for r, g, b in reds)
    title = pixels(menu.actionGeometry(menu.actions()[0]))
    header = QtGui.QColor(theme.COLORS['selection']).getRgb()[:3]
    assert header in title
    window.deleteLater()


# toolbar

def test_toolbar_buttons(window):
    assert [b.text() for b in toolbar_buttons(window)] == [
        'Residues', '', '', 'Zoom', 'Orient', 'Rock', 'Presets...',
        'Builder...', 'Scenes', 'Draw/Ray', '…']


def test_selection_mode(window):
    btn = window.toolbar.selection_mode_button
    assert btn.text() == 'Residues'

    actions(btn.menu())['Chains'].trigger()
    assert cmd.get_setting_int('mouse_selection_mode') == 2
    window.fire_setting_callbacks()
    assert btn.text() == 'Chains'

    # changes from elsewhere update the button
    cmd.set('mouse_selection_mode', 4)
    window.fire_setting_callbacks()
    assert btn.text() == 'Objects'
    assert actions(btn.menu())['Objects'].isChecked()


def test_rock(window):
    btn = window.toolbar.rock_button
    assert not btn.isChecked()

    btn.click()
    assert cmd.get('rock') == 'on'
    btn.click()
    assert cmd.get('rock') == 'off'

    cmd.rock(1)
    window.fire_setting_callbacks()
    assert btn.isChecked()


def test_presets(window):
    cmd.fragment('ala')
    cmd.hide('everything')
    menu = buttons(window)['Presets...'].menu()
    menu.aboutToShow.emit()
    actions(menu)['ball and stick'].trigger()
    assert cmd.count_atoms('rep spheres') > 0


def test_zoom_menu(window):
    menu = buttons(window)['Zoom'].menu()
    menu.aboutToShow.emit()
    assert not actions(menu)['Selection'].isEnabled()
    cmd.fragment('ala')
    cmd.select('sele', 'name CA')
    menu.aboutToShow.emit()
    assert actions(menu)['Selection'].isEnabled()


def test_window_actions(window):
    btns = buttons(window)
    btns['Builder...'].click()
    btns['Scenes'].click()
    actions(btns['…'].menu())['Get View'].trigger()
    actions(btns['…'].menu())['Properties...'].trigger()
    assert window.calls == ['builder', 'scenes', 'get_view', 'properties']

    movie = actions(btns['…'].menu())['Movie'].menu()
    assert list(actions(movie)) == [
        'Rewind', 'Backward', 'Stop', 'Play', 'Forward', 'End', 'Clear Movie']
