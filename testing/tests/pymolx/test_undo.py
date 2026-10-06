'''
pymolx multi-level undo (INCENTIVE_PARITY.md B-05). The Incentive API
itself is covered by the upstream tests in tests/undo/; these cover what
pymolx adds: step grouping, UI settings, failed commands, the camera.
'''

import pytest

import pymol
from pymol import cmd
from pymolx import undo


@pytest.fixture(autouse=True)
def undo_enabled():
    state = (undo.stack._depth, undo.stack._step_open, undo.stack.enabled,
             undo.stack.paused, undo.stack.gui_thread)
    assert state == (0, False, False, 0, None), state
    cmd.reinitialize()
    cmd.undo_enable()
    yield
    undo.stack.defer_close = None
    undo.stack.gui_thread = None
    cmd.undo_disable()
    cmd.reinitialize()


def steps():
    return len(undo.stack.undo_stack)


def test_capability():
    assert 'multi_undo' in pymol.get_capabilities()


def test_top_level_commands_are_steps():
    cmd.fragment('ala')
    cmd.color('red')
    assert steps() == 2
    cmd.undo()
    assert cmd.count_atoms('color red') == 0
    cmd.undo()
    assert cmd.get_names() == []
    cmd.redo()
    cmd.redo()
    assert cmd.count_atoms('color red') == cmd.count_atoms()


def test_nested_commands_are_one_step():
    cmd.fab('ACD', 'pep')  # calls other commands
    assert steps() == 1
    cmd.undo()
    assert cmd.get_names() == []


def test_deferred_close_groups_steps():
    closers = []
    undo.enable_for_gui(closers.append)
    cmd.fragment('ala')
    cmd.color('red')
    assert steps() == 1
    closers.pop()()  # back in the event loop
    cmd.color('blue')
    assert steps() == 2
    cmd.undo()
    assert cmd.count_atoms('color red') == cmd.count_atoms()
    cmd.undo()
    assert cmd.get_names() == []


def test_ui_settings_are_not_steps():
    cmd.set('mouse_selection_mode', 2)
    cmd.set('seq_view', 1)
    assert steps() == 0
    cmd.set('sphere_scale', 0.5)
    assert steps() == 1


def test_failed_command_is_not_a_step():
    with pytest.raises(pymol.CmdException):
        cmd.color('red', 'nosuchobject')
    assert steps() == 0


def test_undo_pause_cm():
    with cmd.UndoPauseCM():
        cmd.fragment('ala')
    assert steps() == 0
    cmd.fragment('gly')
    assert steps() == 1


def test_suspend_undo():
    cmd.set('suspend_undo', 1)
    cmd.fragment('ala')
    assert steps() == 0
    cmd.set('suspend_undo', 0)


def test_new_step_clears_redo():
    cmd.fragment('ala')
    cmd.undo()
    assert undo.stack.can_redo()
    cmd.fragment('gly')
    assert not undo.stack.can_redo()


def test_camera_kept_origin_restored():
    cmd.pseudoatom('m1', pos=[1, 2, 3])
    cmd.origin(position=[4, 5, 6])
    cmd.turn('x', 30)
    view = cmd.get_view()
    cmd.undo()  # undoes origin
    restored = cmd.get_view()
    assert restored[:12] == pytest.approx(view[:12])
    assert restored[12:15] != pytest.approx([4, 5, 6])


def test_max_steps(monkeypatch):
    monkeypatch.setattr(undo, 'MAX_STEPS', 3)
    for i in range(5):
        cmd.pseudoatom('m%d' % i)
    assert steps() == 3


def test_loading_a_script_is_not_a_step(tmp_path):
    script = tmp_path / 'script.pml'
    script.write_text('fragment ala\ncolor red\n')
    cmd.load(str(script))
    assert steps() == 2  # one per command in the script
    cmd.undo()
    assert cmd.count_atoms('color red') == 0


def test_step_state_is_per_thread():
    import threading
    seen = []
    with undo.stack.step(False):
        thread = threading.Thread(
            target=lambda: seen.append(undo.stack._depth))
        thread.start()
        thread.join()
        assert undo.stack._depth == 1
    assert seen == [0]


def test_disabled_keeps_upstream_undo():
    cmd.undo_disable()
    cmd.fragment('ala')
    assert steps() == 0
    cmd.undo()  # upstream coordinate undo: no error
    assert cmd.get_names() == ['ala']


def test_clear_history():
    cmd.fragment('ala')
    cmd.undo()
    undo.stack.clear_history()
    assert not undo.stack.can_undo() and not undo.stack.can_redo()
