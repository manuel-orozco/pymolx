'''
pymolx content panel (INCENTIVE_PARITY.md L-14, L-15): the panel list
from the core, menu dispatch per object type, and the Qt panel itself
(headless, Qt offscreen platform).
'''

import os
from collections import defaultdict

import pytest

from pymol import cmd
from pymolx.panel import get_panel_list, PanelItem


@pytest.fixture(autouse=True)
def clean_session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


def rows():
    return [(it.name, it.type, it.enabled, it.nest_level)
            for it in get_panel_list()]


# panel list

def test_panel_list():
    cmd.fragment('ala', 'm1')
    cmd.select('s1', 'name CA')
    cmd.fragment('gly', 'm2')
    cmd.disable('m2')
    cmd.fragment('ser', '_hidden')
    # loading hides selections (auto_hide_selections); "all" can be
    # disabled by earlier sessions
    enabled = cmd.get_names('all', enabled_only=1)
    assert rows() == [
        ('all', 'all', rows()[0][2], 0),
        ('m1', 'object:molecule', True, 0),
        ('s1', 'selection', 's1' in enabled, 0),
        ('m2', 'object:molecule', False, 0),
    ]


def test_panel_list_groups():
    cmd.fragment('ala', 'm1')
    cmd.fragment('gly', 'm2')
    cmd.group('g1', 'm1 m2')

    cmd.group('g1', action='close')
    items = get_panel_list()
    assert [(it.name, it.is_group, it.is_open) for it in items] == [
        ('all', False, False), ('g1', True, False)]

    cmd.group('g1', action='open')
    assert rows()[1:] == [
        ('g1', 'object:group', True, 0),
        ('m1', 'object:molecule', True, 1),
        ('m2', 'object:molecule', True, 1),
    ]


# menus per object type

def make_objects():
    cmd.fragment('gly', 'mol')
    cmd.set('gaussian_b_floor', 30)
    cmd.map_new('map')
    cmd.isomesh('mesh', 'map')
    cmd.isosurface('surf', 'map')
    cmd.distance('dist', 'mol and name N', 'mol and name C')
    cmd.pseudoatom('pseudo')
    cmd.load_cgo([1.0, 4.0, 0, 0, 0, 1, 0, 0, 2.0], 'cgo')  # BEGIN LINES END
    cmd.select('sele', 'mol')
    cmd.group('grp', 'pseudo')
    cmd.ramp_new('ramp', 'map')
    cmd.slice_new('slice', 'map')
    cmd.volume('vol', 'map')


def test_menu_data_for_all_types():
    pytest.importorskip('pymol.Qt')
    from pymolx.gui import content_panel

    make_objects()
    items = get_panel_list()
    types = {it.type for it in items}
    assert types >= {
        'all', 'selection', 'object:molecule', 'object:map', 'object:mesh',
        'object:surface', 'object:measurement', 'object:cgo', 'object:group',
        'object:ramp', 'object:slice', 'object:volume'}

    for it in items:
        assert it.type in content_panel.MENUS, it.type
        for column, fname in zip('ASHLC', content_panel.MENUS[it.type]):
            data = content_panel.menu_data(it, column, cmd)
            if fname is None:
                assert data is None
            else:
                assert isinstance(data, list) and data, (it, column)


# Qt panel

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
    from pymolx.gui import content_panel

    win = QtWidgets.QMainWindow()
    win.cmd = cmd
    win.setting_callbacks = defaultdict(list)
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
    win.setCentralWidget(QtWidgets.QWidget())

    win.content_panel = content_panel.setup(win)
    win.content_panel.timer.stop()
    yield win
    cmd.set('internal_gui', 1)
    cmd.set('internal_feedback', 1)
    win.deleteLater()


def fire_setting_callbacks(win):
    for setting in cmd.get_setting_updates() or ():
        for callback in win.setting_callbacks.get(setting, ()):
            callback(cmd.get_setting_tuple(setting)[1][0])


def panel_rows(win):
    panel = win.content_panel
    panel.refresh()
    return {row.item.name: row for row in panel.row_widgets()}


@needs_qt
def test_setup_hides_opengl_panel(window):
    assert cmd.get_setting_int('internal_gui') == 0
    assert cmd.get_setting_int('internal_feedback') == 0


@needs_qt
def test_rows_follow_session(window):
    assert list(panel_rows(window)) == ['all']
    cmd.fragment('ala', 'm1')
    cmd.select('s1', 'm1')
    rows = panel_rows(window)
    assert list(rows) == ['all', 'm1', 's1']
    assert rows['s1'].label.text() == '(s1)'
    assert rows['all'].label.text() == 'All'


def click(widget, button=None):
    '''press on widget; unhandled presses propagate to the parent'''
    Qt = QtCore.Qt
    button = button or Qt.MouseButton.LeftButton
    pos = QtCore.QPointF(2, 2)
    event = QtGui.QMouseEvent(QtCore.QEvent.Type.MouseButtonPress, pos,
                              widget.mapToGlobal(pos), button, button,
                              Qt.KeyboardModifier.NoModifier)
    QtWidgets.QApplication.sendEvent(widget, event)


@needs_qt
def test_click_toggles_enabled(window):
    cmd.fragment('ala', 'm1')
    click(panel_rows(window)['m1'].label)
    assert cmd.get_names(enabled_only=1) == []
    row = panel_rows(window)['m1']
    assert row.property('enabled_item') is False
    click(row.label)
    assert cmd.get_names(enabled_only=1) == ['m1']


@needs_qt
def test_button_menus(window):
    cmd.fragment('ala', 'm1')
    row = panel_rows(window)['m1']
    menu = row.buttons['S'].menu()
    menu.aboutToShow.emit()
    labels = [a.text() for a in menu.actions()]
    assert 'sticks' in labels and 'cartoon' in labels

    cmd.hide('everything')
    [a for a in menu.actions() if a.text() == 'sticks'][0].trigger()
    cmd.sync()
    assert cmd.count_atoms('rep sticks') == cmd.count_atoms('m1')


@needs_qt
def test_group_arrow(window):
    cmd.fragment('ala', 'm1')
    cmd.group('g1', 'm1')
    cmd.group('g1', action='open')
    rows = panel_rows(window)
    assert 'm1' in rows
    arrow = rows['g1'].findChild(QtWidgets.QToolButton, 'group_arrow')
    arrow.click()
    assert 'm1' not in panel_rows(window)


@needs_qt
def test_toggles(window):
    toggles = window.content_panel.toggles

    toggles.seq_button.click()
    assert cmd.get_setting_int('seq_view') == 1
    cmd.set('seq_view', 0)
    fire_setting_callbacks(window)
    assert not toggles.seq_button.isChecked()

    assert toggles.command_button.isChecked()
    toggles.command_button.click()
    assert window.browser.isHidden()
    assert not window.lineedit.isHidden()
    toggles.command_button.click()
    assert not window.browser.isHidden()

    cmd.mouse('three_button_editing')
    fire_setting_callbacks(window)
    assert toggles.mouse_button.text() == cmd.get('button_mode_name')
    cmd.mouse('three_button_viewing')


@needs_qt
def test_forward_key_to_command_line(window):
    from pymolx.gui.content_panel import forward_key_to_command_line
    Qt = QtCore.Qt

    def key(text, modifiers=Qt.KeyboardModifier.NoModifier):
        return QtGui.QKeyEvent(QtCore.QEvent.Type.KeyPress, Qt.Key.Key_A,
                               modifiers, text)

    window.lineedit.clear()
    assert forward_key_to_command_line(window, key('z'))
    assert window.lineedit.text() == 'z'
    assert not forward_key_to_command_line(
        window, key('c', Qt.KeyboardModifier.ControlModifier))
    assert not forward_key_to_command_line(window, key('\x1b'))  # Esc
    assert window.lineedit.text() == 'z'


@needs_qt
def test_long_names_keep_buttons_visible(window):
    name = 'fold_lilrb3_d1d4_his_4cmut_model_1_with_a_very_long_name'
    cmd.fragment('ala', name)
    window.resize(900, 600)
    window.show()
    row = panel_rows(window)[name]
    QtWidgets.QApplication.processEvents()
    for column, button in row.buttons.items():
        assert button.geometry().right() <= row.width(), column
    assert name in row.label.toolTip()
    window.hide()


@needs_qt
def test_wizard_panel(window):
    from pymol.wizard import Wizard

    class TestWizard(Wizard):
        def get_panel(self):
            return [
                [1, 'Test Wizard', ''],
                [3, 'Pick a size', 'size'],
                [2, 'Big spheres', 'set sphere_scale, 2'],
                [2, 'Done', 'cmd.set_wizard()'],
            ]

    wizard = TestWizard(_self=cmd)
    wizard.menu['size'] = [[2, 'Size:', ''],
                           [1, '\\900small', 'set sphere_scale, 0.5']]
    panel = window.content_panel.wizard_panel

    panel.refresh()
    assert panel.isHidden()

    cmd.set_wizard(wizard)
    panel.refresh()
    assert not panel.isHidden()
    title = panel.findChild(QtWidgets.QLabel, 'wizard_title')
    assert title.text() == 'Test Wizard'
    buttons = {b.text(): b for b in panel.findChildren(QtWidgets.QAbstractButton)}
    assert set(buttons) == {'Pick a size', 'Big spheres', 'Done'}

    buttons['Big spheres'].click()
    assert cmd.get_setting_float('sphere_scale') == pytest.approx(2)

    menu = buttons['Pick a size'].menu()
    menu.aboutToShow.emit()
    [a for a in menu.actions() if a.text() == 'small'][0].trigger()
    assert cmd.get_setting_float('sphere_scale') == pytest.approx(0.5)

    buttons['Done'].click()
    cmd.sync()
    assert cmd.get_wizard() is None
    panel.refresh()
    assert panel.isHidden()


@needs_qt
def test_state_bar_with_mutagenesis(window):
    bar = window.content_panel.state_bar
    cmd.fab('ACDEF', 'pep')
    bar.refresh()
    assert bar.isHidden()  # one state

    cmd.wizard('mutagenesis')
    try:
        cmd.get_wizard().set_mode('LEU')
        cmd.select('sele', 'pep and resi 3')
        cmd.get_wizard().do_select('sele')
        count = cmd.count_states('mutation')
        assert count > 1

        cmd.frame(1)
        bar.refresh()
        assert not bar.isHidden()
        assert bar.label.text().startswith('State 1/%d' % count)
        assert 'strain' in bar.label.text()

        bar.buttons['next'].click()
        bar.refresh()
        assert bar.label.text().startswith('State 2/%d' % count)
        bar.buttons['last'].click()
        bar.refresh()
        assert bar.label.text().startswith('State %d/%d' % (count, count))
    finally:
        cmd.set_wizard()
