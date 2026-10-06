'''
Mouse mode menu with the mouse controls table (pymolx.gui.mouse_modes).
'''

import os

import pytest

from pymol import cmd

try:
    from pymol.Qt import QtWidgets
except ImportError:
    QtWidgets = None

from pymolx.gui import mouse_modes


@pytest.fixture(autouse=True)
def viewing_mode():
    yield
    cmd.mouse('three_button_viewing')


def test_all_menu_modes_exist():
    from pymol import controlling
    modes = [m for group in mouse_modes.MODE_GROUPS for m in group]
    assert len(modes) == len(set(modes))
    assert set(modes) <= set(controlling.mode_dict)


def test_action_labels():
    assert mouse_modes.action_label('rota') == 'Rota'
    assert mouse_modes.action_label('+Box') == '+Box'
    assert mouse_modes.action_label('none') == ''
    assert mouse_modes.action_label('pktb') == 'PkTB'


def test_controls_table_viewing():
    table = mouse_modes.controls_table('three_button_viewing')
    assert [table[('', c)] for c in ('L', 'M', 'R', 'Wheel')] == [
        'Rota', 'Move', 'MovZ', 'Slab']
    assert table[('Single Click', 'L')] == '+/-'  # toggle selection
    assert table[('Single Click', 'R')] == 'Menu'
    assert ('Single Click', 'Wheel') not in table


def test_controls_table_matches_bindings():
    # every drag binding of every menu mode is in its table
    from pymol import controlling
    for group in mouse_modes.MODE_GROUPS:
        for mode in group:
            table = mouse_modes.controls_table(mode)
            for button, modifier, action in controlling.mode_dict[mode]:
                row = {'none': '', 'shft': 'Shift', 'ctrl': 'Ctrl',
                       'ctsh': 'Ctrl+Shift'}.get(modifier)
                column = {'l': 'L', 'm': 'M', 'r': 'R', 'w': 'Wheel'}.get(
                    button)
                if row is not None and column is not None:
                    assert table[(row, column)] == \
                        mouse_modes.action_label(action), (mode, button)


@pytest.mark.skipif(QtWidgets is None, reason='no Qt')
def test_menu():
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    cmd.mouse('three_button_editing')
    menu = mouse_modes.MouseModeMenu(cmd)
    try:
        labels = [a.text() for a in menu.actions() if not a.isSeparator()]
        assert labels[:2] == ['3-Button Viewing', '3-Button Editing']
        assert labels[-1] == '3-Button Maestro'

        menu.aboutToShow.emit()
        checked = [a.text() for a in menu.actions() if a.isChecked()]
        assert checked == ['3-Button Editing']
        assert menu.panel.title.text() == '3-Button Editing Mouse Controls'

        menu.hovered.emit(menu.mode_actions['two_button_viewing'])
        assert menu.panel.title.text().startswith('2-Button Viewing')

        menu.mode_actions['one_button_viewing'].trigger()
        assert cmd.get('button_mode_name') == '1-Button Viewing'
    finally:
        menu.deleteLater()
