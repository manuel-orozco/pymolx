'''
Unsaved changes: whether the session changed since it was last saved or
opened as a session file, so that quitting can offer to save it
(pymolx.gui.unsaved asks).

A fingerprint of the session (cmd.get_session) is taken when a session
file is saved or opened, and compared when quitting. Left out of it:
the camera (looking around isn't a change), the window, and settings
which are GUI state or bookkeeping. Mouse edits (dragged atoms,
sculpting) change it, unlike tracking commands.
'''
import functools
import hashlib
import inspect
import pickle

from pymol import cmd

SESSION_FORMATS = {'pse', 'psw', 'pze', 'pzw'}
IGNORED_KEYS = {'view', 'main'}    # camera, window
# settings of GUI state and bookkeeping (pse_binary_dump: set by the first
# cmd.get_session)
IGNORED_SETTINGS = {'pse_binary_dump'}

_baseline = None       # fingerprint when last saved or opened; None: never
_quit_handler = None   # the GUI's quit, asking first (set by the window)


def _ignored_settings():
    from pymol import setting
    from . import undo
    return {setting.index_dict[name]
            for name in undo.UI_SETTINGS | IGNORED_SETTINGS
            if name in setting.index_dict}


def fingerprint(_self=cmd):
    '''
    :return: a hash of what the session holds (None if it can't be read)
    '''
    try:
        session = _self.get_session()
    except Exception:
        return None
    ignored = _ignored_settings()
    kept = {key: value for key, value in session.items()
            if key not in IGNORED_KEYS}
    kept['settings'] = [s for s in session.get('settings') or []
                        if s and s[0] not in ignored]
    return hashlib.sha1(pickle.dumps(kept, pickle.HIGHEST_PROTOCOL)) \
        .hexdigest()


def mark_saved(_self=cmd):
    '''
    The session as it is now counts as saved
    '''
    global _baseline
    _baseline = fingerprint(_self)


def has_unsaved(_self=cmd):
    '''
    True if the session holds objects and changed since it was last saved
    or opened (never saved: changed)
    '''
    if not _self.get_names('all'):
        return False
    now = fingerprint(_self)
    return now is None or now != _baseline


def quit_asking(_self=cmd):
    '''
    Quit, in the GUI after offering to save the session (for menus; a
    script's cmd.quit() doesn't ask)
    '''
    if _quit_handler is not None:
        _quit_handler()
    else:
        _self.quit()


def _format(bound):
    fmt = str(bound.get('format') or '').lower()
    if not fmt:
        from pymol.importing import filename_to_format
        fmt = filename_to_format(str(bound.get('filename', '')))[2]
    return fmt


def _wrap(func):
    '''
    Wrapper for save and load: saving the session to a session file, or
    opening one in place of the session, makes it the saved one
    '''
    signature = inspect.signature(func)

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        try:
            bound = signature.bind_partial(*args, **kwargs).arguments
        except TypeError:
            return result
        # (partial: merged into the session, or part of it saved)
        if _format(bound) in SESSION_FORMATS and \
                not int(bound.get('partial') or 0):
            mark_saved(kwargs.get('_self', cmd))
        return result
    wrapper._pymolx_unsaved = True
    return wrapper


def install(_self=cmd):
    '''
    Wrap save and load. Called from pymolx._init, before undo.install.
    '''
    wrapped = {}
    for name in ('save', 'load'):
        func = getattr(_self, name)
        if getattr(func, '_pymolx_unsaved', False):
            continue
        wrapped[func] = _wrap(func)
        setattr(_self, name, wrapped[func])
    # command language entries point at the same functions
    for entry in _self.keyword.values():
        if entry and entry[0] in wrapped:
            entry[0] = wrapped[entry[0]]
