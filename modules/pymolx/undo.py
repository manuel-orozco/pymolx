'''
Multi-level undo (parity item B-05) with Incentive PyMOL's API:
cmd.undo_enable(), cmd.undo_disable(), cmd.undo(), cmd.redo() and the
"multi_undo" capability.

Implemented with session snapshots. Before the first state-changing
command (UNDOABLE) of a step, the session is saved (pickled and
compressed). Undo restores the snapshot but keeps the camera; the
"suspend_undo" setting pauses recording.

A step is one top-level command. A GUI can group further with
UndoStack.defer_close: then a step lasts until control returns to its
event loop, so a menu action, a command line entry or a script run from
the GUI is one step.

When undo is disabled, cmd.undo and cmd.redo keep the upstream behavior
(undo coordinate changes of the object being edited).
'''

import contextlib
import functools
import inspect
import pickle
import threading
import time
import zlib

import pymol
from pymol import cmd

MAX_STEPS = 50
MAX_BYTES = 512 * 2**20
SLOW_SNAPSHOT_SECONDS = 1.0

# commands which change the session (camera-only commands are left out)
UNDOABLE = '''
    alter alter_state align alphatoall angle attach bg_color bond cartoon
    cealign clean color copy create deprotect delete deselect dihedral
    disable distance dss enable extract fab fetch fit fix_chemistry flag
    fragment fuse group h_add h_fill hide intra_fit invert isodot isolevel
    isomesh isosurface label load load_brick load_cgo load_coords
    load_coordset load_map load_model load_raw load_traj madd map_double
    map_halve map_new map_set map_trim mask matrix_copy matrix_reset mclear
    mdelete mset mview order origin pair_fit protect pseudoatom ramp_new
    rebond reinitialize remove remove_picked rename replace rotate scene
    sculpt_iterate select set set_bond set_color set_dihedral set_geometry
    set_name set_object_color set_object_ttt set_raw_alignment set_symmetry
    set_title show show_as slice_new smooth sort spectrum split_chains
    split_states super symexp toggle transform_object transform_selection
    translate unbond ungroup unmask unset unset_bond update valence volume
    volume_color
'''.split()

# "set"/"unset" of these are GUI state or bookkeeping, not undo steps
UI_SETTINGS = {
    'defer_updates', 'suspend_updates', 'suspend_undo', 'text',
    'internal_gui', 'internal_gui_width', 'internal_feedback',
    'internal_prompt', 'session_file', 'mouse_selection_mode', 'seq_view',
    'movie_panel',
}


class UndoStack:
    '''
    Snapshot stacks and the recording state. One instance per cmd module.

    :ivar defer_close: None, or callable(callback) which calls callback
        once control is back in the GUI event loop (ends the step). Only
        used for commands on the GUI thread (gui_thread).

    Step state is per thread: commands from other threads are steps of
    their own and don't disturb the GUI thread's grouping.
    '''

    def __init__(self, _self):
        self.cmd = _self
        self.enabled = False
        self.undo_stack = []
        self.redo_stack = []
        self.defer_close = None
        self.gui_thread = None
        self.paused = 0
        self._lock = threading.RLock()
        self._local = threading.local()
        self._restoring = False
        self._warned_slow = False

    # step state of the current thread

    @property
    def _depth(self):
        return getattr(self._local, 'depth', 0)

    @_depth.setter
    def _depth(self, value):
        self._local.depth = value

    @property
    def _step_open(self):
        return getattr(self._local, 'step_open', False)

    @_step_open.setter
    def _step_open(self, value):
        self._local.step_open = value

    def _defer_close(self):
        if self.gui_thread == threading.get_ident():
            return self.defer_close
        return None

    # snapshots

    def _snapshot(self):
        t0 = time.time()
        session = self.cmd.get_session()
        blob = zlib.compress(pickle.dumps(session, pickle.HIGHEST_PROTOCOL), 1)
        seconds = time.time() - t0
        if seconds > SLOW_SNAPSHOT_SECONDS and not self._warned_slow:
            self._warned_slow = True
            print(' Undo: saving undo steps takes %.1f seconds for this '
                  'session; "undo_disable" turns undo off.' % seconds)
        return blob

    def _restore(self, blob):
        _self = self.cmd
        view = list(_self.get_view())
        session = pickle.loads(zlib.decompress(blob))
        self._restoring = True
        try:
            _self.set_session(session)
            # keep the camera, but restore the origin of rotation (session
            # views have a 4x4 rotation matrix, the origin is at 19:22)
            view[12:15] = session['view'][19:22]
            _self.set_view(view)
            # rebuild derived state, e.g. alignment object selections
            _self.refresh()
        finally:
            self._restoring = False

    def _trim(self):
        while len(self.undo_stack) > MAX_STEPS:
            del self.undo_stack[0]
        while (len(self.undo_stack) > 1 and
               sum(map(len, self.undo_stack + self.redo_stack)) > MAX_BYTES):
            del self.undo_stack[0]

    def _recording(self):
        return (self.enabled and not self._restoring and not self.paused and
                not self.cmd.get_setting_int('suspend_undo'))

    def _close_step(self):
        # called on the GUI thread by defer_close
        self._step_open = False

    @contextlib.contextmanager
    def step(self, changes):
        '''
        Context for running a command. Before the first command with
        changes=True in a step, the session is saved.
        '''
        pushed = False
        if changes and not self._step_open and self._recording():
            blob = self._snapshot()
            with self._lock:
                self.undo_stack.append(blob)
                self.redo_stack.clear()
                self._trim()
            self._step_open = pushed = True
            defer_close = self._defer_close()
            if defer_close is not None:
                defer_close(self._close_step)
        self._depth += 1
        try:
            yield
        except BaseException:
            if pushed:
                # failed before changing anything: not a step
                with self._lock:
                    if self.undo_stack and self.undo_stack[-1] is blob:
                        self.undo_stack.pop()
                self._step_open = False
            raise
        finally:
            self._depth -= 1
            if self._depth == 0 and self._defer_close() is None:
                self._step_open = False

    # API

    def can_undo(self):
        return self.enabled and bool(self.undo_stack)

    def can_redo(self):
        return self.enabled and bool(self.redo_stack)

    def undo(self):
        with self._lock:
            if not self.undo_stack:
                print(' Undo: nothing to undo.')
                return False
            current = self._snapshot()
            self._restore(self.undo_stack.pop())
            self.redo_stack.append(current)
            self._trim()
        return True

    def redo(self):
        with self._lock:
            if not self.redo_stack:
                print(' Redo: nothing to redo.')
                return False
            current = self._snapshot()
            self._restore(self.redo_stack.pop())
            self.undo_stack.append(current)
            self._trim()
        return True

    def enable(self):
        self.enabled = True

    def disable(self):
        self.enabled = False
        with self._lock:
            self.undo_stack.clear()
            self.redo_stack.clear()


stack = UndoStack(cmd)


def _wrap(func, changes):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with stack.step(changes):
            return func(*args, **kwargs)
    wrapper._pymolx_undo = True
    return wrapper


# "load" of these runs a script: its commands are the steps, not the load
SCRIPT_FORMATS = {'pml', 'py', 'pym', 'pyc', 'p1m', 'pwg'}


def _wrap_load(func):
    '''
    Wrapper for load: loading data is a step, running a script is not.
    '''
    signature = inspect.signature(func)

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            bound = signature.bind_partial(*args, **kwargs).arguments
        except TypeError:
            bound = {}
        fmt = bound.get('format') or ''
        if not fmt:
            from pymol.importing import filename_to_format
            fmt = filename_to_format(str(bound.get('filename', '')))[2]
        if fmt in SCRIPT_FORMATS:
            return func(*args, **kwargs)
        with stack.step(True):
            return func(*args, **kwargs)
    wrapper._pymolx_undo = True
    return wrapper


def _wrap_setting(func):
    '''
    Wrapper for set/unset: changes to UI_SETTINGS are not undo steps.
    '''
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        name = args[0] if args else kwargs.get('name', '')
        with stack.step(str(name) not in UI_SETTINGS):
            return func(*args, **kwargs)
    wrapper._pymolx_undo = True
    return wrapper


def enable_for_gui(defer_close):
    '''
    Turn undo on for an interactive session.

    :param defer_close: see UndoStack.defer_close
    '''
    stack.defer_close = defer_close
    stack.gui_thread = threading.get_ident()
    stack.enable()


def install(_self=cmd):
    '''
    Wrap the UNDOABLE commands and replace undo/redo.
    Called once at startup from pymolx._init.
    '''
    wrapped = {}
    for name in UNDOABLE:
        func = getattr(_self, name, None)
        if func is None or getattr(func, '_pymolx_undo', False):
            continue
        if name in ('set', 'unset'):
            wrapped[func] = _wrap_setting(func)
        elif name == 'load':
            wrapped[func] = _wrap_load(func)
        else:
            wrapped[func] = _wrap(func, True)
        setattr(_self, name, wrapped[func])

    # command language entries point at the same functions
    for entry in _self.keyword.values():
        if entry and entry[0] in wrapped:
            entry[0] = wrapped[entry[0]]

    upstream_undo = _self.undo
    upstream_redo = _self.redo

    def undo(*, _self=_self):
        '''
DESCRIPTION

    "undo" reverts the last change. With undo disabled ("undo_disable"),
    it restores the previous conformation of the object being edited.

USAGE

    undo

SEE ALSO

    redo, undo_enable, undo_disable
        '''
        if stack.enabled:
            return stack.undo()
        return upstream_undo(_self=_self)

    def redo(*, _self=_self):
        '''
DESCRIPTION

    "redo" reapplies the last change which was undone.

USAGE

    redo

SEE ALSO

    undo, undo_enable, undo_disable
        '''
        if stack.enabled:
            return stack.redo()
        return upstream_redo(_self=_self)

    def undo_enable(*, _self=_self):
        '''
DESCRIPTION

    "undo_enable" turns on multi-level undo of changes.

USAGE

    undo_enable

SEE ALSO

    undo, redo, undo_disable
        '''
        stack.enable()

    def undo_disable(*, _self=_self):
        '''
DESCRIPTION

    "undo_disable" turns off multi-level undo and forgets saved steps.

USAGE

    undo_disable

SEE ALSO

    undo, redo, undo_enable
        '''
        stack.disable()

    for func in (undo, redo, undo_enable, undo_disable):
        _self.extend(func.__name__, func)
        setattr(_self, func.__name__, func)

    @contextlib.contextmanager
    def UndoPauseCM(*, _self=_self):
        '''
        Context manager: commands inside don't create undo steps.
        '''
        stack.paused += 1
        try:
            yield
        finally:
            stack.paused -= 1

    _self.UndoPauseCM = UndoPauseCM

    # advertise like Incentive PyMOL
    get_capabilities = pymol.get_capabilities
    if not getattr(get_capabilities, '_pymolx_undo', False):
        def capabilities():
            return set(get_capabilities()) | {'multi_undo'}
        capabilities._pymolx_undo = True
        pymol.get_capabilities = capabilities
