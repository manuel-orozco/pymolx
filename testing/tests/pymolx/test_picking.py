'''
What a click in the viewer does to the selection (pymolx.picking,
INCENTIVE_PARITY.md: L-33). The two added button actions ("adds" and
"subs") are pymolx's; PyMOL's own "+/-" flips what was clicked.
'''

import pytest

from pymol import cmd
from pymol import controlling
from pymolx import picking


@pytest.fixture(autouse=True)
def restore():
    yield
    picking.set_click_action(picking.DEFAULT_ACTION)


def test_actions_are_button_actions():
    # every action names a binding PyMOL knows
    for action, code in picking.CLICK_ACTIONS.items():
        assert code in controlling.but_act_code, action
    # pymolx's own two, next to the ones PyMOL came with
    assert controlling.but_act_code['adds'] == 57
    assert controlling.but_act_code['subs'] == 58
    assert controlling.but_act_code['+/-'] == 31
    assert picking.DEFAULT_ACTION in picking.CLICK_ACTIONS


def test_set_and_read_back():
    for action in sorted(picking.CLICK_ACTIONS):
        assert picking.set_click_action(action) == action
        assert picking.click_action() == action
    with pytest.raises(ValueError):
        picking.set_click_action('nonsense')
    assert picking.click_action() == 'toggle'   # unchanged by the failure


def test_reapply_after_a_mouse_mode_change():
    # picking a mouse mode rebinds every button, so the toolbar's choice
    # has to be put back
    picking.set_click_action('add')
    cmd.mouse('three_button_viewing')
    assert picking.click_action() == 'add'
    picking.reapply()
    assert picking.click_action() == 'add'


def test_labels_for_the_new_actions():
    mouse_modes = pytest.importorskip('pymolx.gui.mouse_modes')
    assert mouse_modes.action_label('adds') == 'Add'
    assert mouse_modes.action_label('subs') == 'Sub'
