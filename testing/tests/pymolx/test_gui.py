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
    # select tool, selection level | undo redo | zoom all, zoom selection
    # | presets ... (the drawn buttons paint their own label, so no text)
    bar = window.toolbar
    assert [b.text() for b in toolbar_buttons(window)] == [
        '', '', '', '', '', '', '',
        'Builder...', 'Scenes', 'Draw/Ray', '\u2022\u2022\u2022']
    assert toolbar_buttons(window)[:7] == [
        bar.select_button, bar.selection_mode_button,
        bar.undo_button, bar.redo_button,
        bar.zoom_all_button, bar.zoom_selection_button, bar.presets_button]
    assert bar.selection_mode_button.letter == 'R'
    # + and - sit between the arrow and the level, stacked in a column
    column = bar.add_button_plus.parent()
    assert column is bar.add_button_minus.parent()
    assert (bar.add_button_plus.text(), bar.add_button_minus.text()) == \
        ('+', '\u2212')
    # every button in the left half is an icon (or the level's letter)
    # over a caption, so each one fills the bar's height
    assert [b.caption for b in toolbar_buttons(window)[:7]] == [
        'SELECT', 'RESIDUES', 'UNDO', 'REDO', 'ZOOM ALL', 'ZOOM SEL',
        'PRESETS']
    for btn in toolbar_buttons(window)[:7]:
        assert btn.sizeHint().height() == bar.presets_button.sizeHint().height()
        assert not btn.grab().isNull()      # they paint


def test_selection_mode(window):
    # the level's letter with its name underneath, as Maestro shows it
    btn = window.toolbar.selection_mode_button
    assert (btn.letter, btn.caption) == ('R', 'RESIDUES')

    actions(btn.menu())['Chains  (C)'].trigger()
    assert cmd.get_setting_int('mouse_selection_mode') == 2
    window.fire_setting_callbacks()
    assert (btn.letter, btn.caption) == ('C', 'CHAINS')

    # changes from elsewhere update the button
    cmd.set('mouse_selection_mode', 4)
    window.fire_setting_callbacks()
    assert (btn.letter, btn.caption) == ('O', 'OBJECTS')
    assert actions(btn.menu())['Objects  (O)'].isChecked()

    # every level PyMOL has, each with its own letter
    from pymolx.gui import toolbar
    labels = [a.text() for a in btn.menu().actions()]
    assert len(labels) == len(toolbar.SELECTION_MODES) == 7
    letters = [letter for _, _, letter in toolbar.SELECTION_MODES]
    assert letters == ['A', 'R', 'C', 'S', 'O', 'M', 'CA']
    assert len(set(letters)) == len(letters)
    for label, value, letter in toolbar.SELECTION_MODES:
        cmd.set('mouse_selection_mode', value)
        window.fire_setting_callbacks()
        assert (btn.letter, btn.caption) == (letter, label.upper())
    # it is wide enough for the longest of them
    cmd.set('mouse_selection_mode', 1)
    window.fire_setting_callbacks()
    assert btn.sizeHint().width() >= 40
    assert btn.sizeHint().height() > btn.fontMetrics().height()
    assert not btn.grab().isNull()      # it paints


def test_select_tool(window):
    from pymolx import picking
    bar = window.toolbar
    # the arrow is on and nothing else: a click replaces the selection
    assert bar.select_button.isChecked()
    assert picking.click_action() == 'replace'

    bar.add_button_plus.setChecked(True)
    assert picking.click_action() == 'add'
    # + and - are alternatives
    bar.add_button_minus.setChecked(True)
    assert not bar.add_button_plus.isChecked()
    assert picking.click_action() == 'subtract'
    bar.add_button_minus.setChecked(False)
    assert picking.click_action() == 'replace'

    # with the arrow off, clicks leave the selection alone and + and -
    # can't be pressed
    bar.add_button_plus.setChecked(True)
    bar.select_button.setChecked(False)
    assert picking.click_action() == 'off'
    assert not bar.add_button_plus.isEnabled()
    assert not bar.add_button_minus.isEnabled()
    bar.select_button.setChecked(True)
    assert bar.add_button_plus.isEnabled()
    assert picking.click_action() == 'add'
    bar.add_button_plus.setChecked(False)



def test_select_tool_survives_a_mouse_mode_change(window):
    # cmd.mouse rebinds every button, "+/-" included, so the toolbar
    # binds its own choice again when button_mode_name changes
    from pymolx import picking
    bar = window.toolbar
    bar.add_button_plus.setChecked(True)
    assert picking.click_action() == 'add'

    cmd.mouse('three_button_editing')
    window.fire_setting_callbacks()
    assert picking.click_action() == 'add'
    assert bar.add_button_plus.isChecked()
    bar.add_button_plus.setChecked(False)
    cmd.mouse('three_button_viewing')

def test_rock_in_the_more_menu(window):
    menu = buttons(window)['\u2022\u2022\u2022'].menu()
    action = actions(menu)['Rock']
    menu.aboutToShow.emit()
    assert not action.isChecked()

    action.trigger()
    assert cmd.get('rock') == 'on'
    action.trigger()
    assert cmd.get('rock') == 'off'

    # rocking started elsewhere shows up when the menu opens
    cmd.rock(1)
    menu.aboutToShow.emit()
    assert action.isChecked()


def test_presets(window):
    # an icon with its name underneath, like the selection mode's letter
    btn = window.toolbar.presets_button
    assert btn.caption == 'PRESETS'
    assert btn.glyph is not None and not btn.glyph.isNull()
    assert not btn.letter and not btn.text()
    # the two drawn buttons are the same height, so the toolbar lines up
    assert btn.sizeHint().height() == \
        window.toolbar.selection_mode_button.sizeHint().height()
    assert not btn.grab().isNull()

    cmd.fragment('ala')
    cmd.hide('everything')
    menu = btn.menu()
    menu.aboutToShow.emit()
    # the presets PyMOL offers are all there
    labels = [a.text() for a in menu.actions() if a.text()]
    for one in ('simple', 'technical', 'ligands', 'pretty',
                'publication', 'ball and stick', 'default'):
        assert one in labels
    actions(menu)['ball and stick'].trigger()
    assert cmd.count_atoms('rep spheres') > 0


def test_zoom_buttons(window):
    bar = window.toolbar
    cmd.fragment('ala')
    cmd.zoom('none')
    before = cmd.get_view()
    bar.zoom_all_button.click()
    assert cmd.get_view() != before

    # zooming on the selection needs one
    bar.zoom_timer.timeout.emit()
    assert not bar.zoom_selection_button.isEnabled()
    cmd.select('sele', 'name CA')
    bar.zoom_timer.timeout.emit()
    assert bar.zoom_selection_button.isEnabled()
    before = cmd.get_view()
    bar.zoom_selection_button.click()
    assert cmd.get_view() != before
    cmd.delete('sele')
    bar.zoom_timer.timeout.emit()
    assert not bar.zoom_selection_button.isEnabled()


def test_window_actions(window):
    btns = buttons(window)
    btns['Builder...'].click()
    btns['Scenes'].click()
    actions(btns['•••'].menu())['Get View'].trigger()
    actions(btns['•••'].menu())['Properties...'].trigger()
    assert window.calls == ['builder', 'scenes', 'get_view', 'properties']

    movie = actions(btns['•••'].menu())['Movie'].menu()
    assert list(actions(movie)) == [
        'Rewind', 'Backward', 'Stop', 'Play', 'Forward', 'End', 'Clear Movie']


def test_console_font(app):
    font = theme.console_font(['Arial', 'Consolas', 'DejaVu Sans Mono'])
    assert font.family() == 'DejaVu Sans Mono'
    font = theme.console_font(['Arial', 'Consolas'])  # e.g. Windows
    assert font.family() == 'Consolas'
    # letters as tall as the interface font's at the same size (within a
    # pixel: metrics come in whole pixels)
    for size in (9, theme.FONT_SIZE, 16):
        font = theme.console_font(['DejaVu Sans Mono'], size=size)
        ui = QtGui.QFont(QtWidgets.QApplication.font())
        ui.setPointSizeF(size)
        assert QtGui.QFontMetricsF(font).capHeight() == pytest.approx(
            QtGui.QFontMetricsF(ui).capHeight(), abs=1)
    # installed fonts: always a font
    assert theme.console_font().family()


def test_one_text_size_no_bold(app):
    # all GUI text the same size (FONT_SIZE); bold only for titles
    import re
    from pymol import cmd
    qss = theme.stylesheet(cmd)
    sizes = set(re.findall(r'font-size:\s*([^;]+);', qss))
    # symbols, not text: the entry dot, "...", the selecting tool's + and -
    assert sizes == {'%dpt' % theme.FONT_SIZE, '8px', '15px', '14px'}
    bold = [m.group(1).strip().split('\n')[-1] for m in re.finditer(
        r'([^{}]*)\{[^}]*font-weight:\s*bold', qss)]
    # titles: the wizard's, Structure Preparation's workflow steps
    assert bold == ['QLabel#wizard_title', 'QGroupBox#prep_step::title']


def test_gui_startup_has_no_local_pymolx_name():
    # an "import pymolx..." inside execapp() made "pymolx" a local name
    # there and the GUI failed to start (UnboundLocalError)
    from pmg_qt import pymol_qt_gui
    assert 'pymolx' not in pymol_qt_gui.execapp.__code__.co_varnames
