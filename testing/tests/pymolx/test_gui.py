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
    from pymol.Qt import QtCore, QtGui, QtWidgets
except ImportError:
    pytest.skip('no Qt', allow_module_level=True)

import shiboken6

from pymolx.gui import menus, plugin_menu, theme, toolbar


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
    # | presets, labels | plugins | mouse mode | wizard | sequence,
    # timeline, output ... (the drawn buttons paint their own label, so
    # no text)
    bar = window.toolbar
    assert [b.text() for b in toolbar_buttons(window)] == [
        '', '', '', '', '', '', '', '', '', '', '', '', '', '',
        'Builder...', 'Scenes', 'Draw/Ray', '\u2022\u2022\u2022']
    assert toolbar_buttons(window)[:14] == [
        bar.select_button, bar.selection_mode_button,
        bar.undo_button, bar.redo_button,
        bar.zoom_all_button, bar.zoom_selection_button, bar.presets_button,
        bar.labels_button,
        bar.plugins_button, bar.mouse_button, bar.wizard_button,
        bar.seq_button, bar.timeline_button, bar.command_button]
    assert bar.selection_mode_button.letter == 'R'
    # + and - sit between the arrow and the level, stacked in a column
    column = bar.add_button_plus.parent()
    assert column is bar.add_button_minus.parent()
    assert (bar.add_button_plus.text(), bar.add_button_minus.text()) == \
        ('+', '\u2212')
    # every button in the left half is an icon (or the level's letter)
    # over a caption, so each one fills the bar's height
    assert [b.caption for b in toolbar_buttons(window)[:14]] == [
        'SELECT', 'RESIDUES', 'UNDO', 'REDO', 'ZOOM ALL', 'ZOOM SEL',
        'PRESETS', 'LABELS', 'PLUGINS', '3-BUTTON VIEWING', 'WIZARD',
        'SEQUENCE', 'TIMELINE', 'OUTPUT']
    for btn in toolbar_buttons(window)[:14]:
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


def test_plugin_groups():
    # the labels the plugins give themselves, as the Plugin menu has them
    startup = os.path.join(os.path.dirname(__file__), '..', '..', '..',
                           'data', 'startup')
    labels = []
    for name in sorted(os.listdir(startup)):
        path = os.path.join(startup, name, '__init__.py')
        if os.path.exists(path):
            with open(path) as handle:
                labels += re.findall(r"addmenuitem(?:qt)?\('([^']+)'",
                                     handle.read())
    assert len(labels) >= 16
    groups, last = plugin_menu.grouped(
        ['Plugin Manager'] + labels + ['Somebody Else\'s Tool'])
    assert groups == [
        ('Sequence Tools', ['Sequence Comparison', 'Modify Sequences']),
        ('Structure Tools', [
            'Structure Comparison', 'Structure Prediction',
            'Structure Preparation', 'Interface Analysis (PISA-like)',
            'Alignment/Superposition']),
        ('Clustering', ['RMSD Clustering', 'Interface RMSD Clustering']),
        ('Antibody Tools', ['Antibody Builder', 'CDR Loops',
                            'Epitope Prediction', 'Paratope Prediction']),
        ('Ligand Tools', ['Binding Pocket']),
        ('Physicochemical Tools', ['Aggregation Prediction',
                                   'APBS Electrostatics']),
        ('Miscellaneous', ['Lighting Settings', 'Somebody Else\'s Tool']),
    ]
    assert last == ['Plugin Manager']
    # every name in a group is a plugin there is (no misspelt label)
    for _heading, names in plugin_menu.GROUPS:
        assert set(names) <= set(labels)
    # a group nothing is installed of is left out
    assert plugin_menu.grouped(['Binding Pocket']) == (
        [('Ligand Tools', ['Binding Pocket'])], [])
    assert plugin_menu.grouped([]) == ([], [])


def test_plugins_button(window):
    btn = window.toolbar.plugins_button
    assert btn.caption == 'PLUGINS'
    assert btn.glyph is not None and not btn.glyph.isNull()
    assert btn.toolTip() == 'Open a plugin'
    assert not btn.grab().isNull()
    menu = btn.menu()
    # a window without a Plugin menu: said, not an empty list
    menu.aboutToShow.emit()
    assert [a.text() for a in menu.actions()] == ['No plugins']
    assert not menu.actions()[0].isEnabled()

    # the Plugin menu as PyMOLQtGUI.initializePlugins sets it up
    opened = []
    plugin = QtWidgets.QMenu('Plugin', window)
    plugin.addAction('Plugin Manager', lambda: opened.append('manager'))
    plugin.addSeparator()
    legacy = QtWidgets.QMenu('Legacy Plugins', plugin)
    plugin.addMenu(legacy)
    plugin.addSeparator()
    for label in ('Lighting Settings', 'CDR Loops', 'Binding Pocket',
                  'RMSD Clustering', 'Antibody Builder'):
        plugin.addAction(label, lambda label=label: opened.append(label))
    legacy.addAction('Old Tool', lambda: opened.append('Old Tool'))
    window.menudict = {'PluginQt': plugin, 'Plugin': legacy}

    def listed():
        return [None if a.isSeparator() else a.text()
                for a in menu.actions()]

    menu.aboutToShow.emit()
    assert listed() == [
        'Clustering', 'RMSD Clustering', None,
        'Antibody Tools', 'Antibody Builder', 'CDR Loops', None,
        'Ligand Tools', 'Binding Pocket', None,
        'Miscellaneous', 'Old Tool', 'Lighting Settings', None,
        'Plugin Manager']
    # the headings cannot be clicked; the entries open their plugin
    assert [a.text() for a in menu.actions() if not a.isEnabled()
            and not a.isSeparator()] == ['Clustering', 'Antibody Tools',
                                         'Ligand Tools', 'Miscellaneous']
    actions(menu)['CDR Loops'].trigger()
    actions(menu)['Old Tool'].trigger()
    actions(menu)['Plugin Manager'].trigger()
    assert opened == ['CDR Loops', 'Old Tool', 'manager']

    # a plugin installed since is there the next time, and showing the
    # list again does not take the entries out of the Plugin menu
    plugin.addAction('Paratope Prediction', lambda: None)
    menu.aboutToShow.emit()
    assert listed()[3:8] == ['Antibody Tools', 'Antibody Builder',
                             'CDR Loops', 'Paratope Prediction', None]
    assert len(listed()) == 16
    assert [a.text() for a in plugin.actions() if a.text()] == [
        'Plugin Manager', 'Legacy Plugins', 'Lighting Settings',
        'CDR Loops', 'Binding Pocket', 'RMSD Clustering',
        'Antibody Builder', 'Paratope Prediction']
    assert not menu.grab().isNull()

    # before the plugins are loaded the Plugin menu only offers to
    window.menudict = {'Plugin': QtWidgets.QMenu('Plugin', window)}
    window.menudict['Plugin'].addAction('Initialize Plugin System')
    menu.aboutToShow.emit()
    assert listed() == ['Initialize Plugin System']


def test_labels_button(window):
    # opens the Labels window, the same one each time
    from pymolx.gui import label_dialog
    btn = window.toolbar.labels_button
    assert btn.caption == 'LABELS' and btn.menu() is None
    try:
        btn.click()
        dialog = label_dialog._dialog
        assert isinstance(dialog, label_dialog.LabelDialog)
        assert dialog.parent() is window and dialog.isVisible()
        dialog.hide()
        btn.click()
        assert label_dialog._dialog is dialog and dialog.isVisible()
    finally:
        shiboken6.delete(label_dialog._dialog)
        label_dialog._dialog = None


def test_mouse_mode_button(window):
    # the mode's name under the mouse, as the selection level's is under
    # its letter; the menu is the mode list with the controls table
    from pymolx.gui import mouse_modes
    btn = window.toolbar.mouse_button
    assert btn.caption == '3-BUTTON VIEWING'
    assert isinstance(btn.menu(), mouse_modes.MouseModeMenu)
    width = btn.sizeHint().width()
    try:
        for mode, name in mouse_modes.mode_names().items():
            cmd.mouse(mode)
            window.fire_setting_callbacks()
            assert btn.caption == name.upper()
            # wide enough for any mode: the buttons after it stay put
            assert btn.sizeHint().width() == width
            assert not btn.grab().isNull()
    finally:
        cmd.mouse('three_button_viewing')


def test_toggle_buttons(app):
    # the switches that were a strip under the content panel
    cmd.reinitialize()
    cmd.set('movie_panel', 0)       # as content_panel.setup leaves it
    cmd.get_setting_updates()
    win = StubWindow()
    win.menudict = {'Wizard': QtWidgets.QMenu('Wizard', win)}
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
    bar = win.toolbar = toolbar.setup(win)
    try:
        # the Wizard menu is the menu bar's own
        assert bar.wizard_button.menu() is win.menudict['Wizard']

        # switches are pressed in while on, and show their "-on" icon
        for btn in (bar.seq_button, bar.timeline_button,
                    bar.command_button):
            assert btn.isCheckable() and not btn.isChecked()
        off = bar.timeline_button.grab().toImage()
        bar.timeline_button.click()
        assert cmd.get_setting_int('movie_panel') == 1
        assert bar.timeline_button.grab().toImage() != off
        # ... and follow the setting
        cmd.set('movie_panel', 0)
        win.fire_setting_callbacks()
        assert not bar.timeline_button.isChecked()

        # the output pane; the command line stays
        bar.set_output_shown(False)
        assert win.browser.isHidden() and not bar.command_button.isChecked()
        bar.command_button.click()
        assert not win.browser.isHidden()
        assert not win.lineedit.isHidden()
        bar.command_button.click()
        assert win.browser.isHidden() and not win.lineedit.isHidden()
        bar.set_output_shown(True)
        assert not win.browser.isHidden() and bar.command_button.isChecked()
    finally:
        cmd.reinitialize()
        shiboken6.delete(win)


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
