'''
What a click in the viewer does to the selection (parity item L-33).

PyMOL binds a plain left click to "+/-", which flips whatever was
clicked: it adds an atom that wasn't selected and takes away one that
was. Maestro instead has a selecting tool with separate "add" and "take
away" modes, which is what the toolbar offers. The two button actions
they need are pymolx's (cButModeSeleAdd and cButModeSeleSub in
layer1/ButMode.h, "adds" and "subs" in pymol.controlling); "replace",
"toggle" and "off" are PyMOL's own.

PyMOL can bind a button but not read the binding back, and changing the
mouse mode rebinds everything, so the chosen action is kept here and
"reapply" puts it back.
'''

from pymol import cmd

# what a plain single left click does, and the button action behind it
CLICK_ACTIONS = {
    'replace': 'sele',   # the clicked atom becomes the selection
    'add': 'adds',       # pymolx: added to it
    'subtract': 'subs',  # pymolx: taken from it
    'toggle': '+/-',     # PyMOL's own: added if absent, removed if there
    'off': 'none',       # clicks leave the selection alone
}
DEFAULT_ACTION = 'replace'

_action = DEFAULT_ACTION


def set_click_action(action, _self=cmd):
    '''
    Bind a plain single left click in the viewer

    :param action: a key of CLICK_ACTIONS
    :return: the action
    '''
    global _action
    if action not in CLICK_ACTIONS:
        raise ValueError('click action: one of %s'
                         % ', '.join(sorted(CLICK_ACTIONS)))
    _action = action
    _self.button('single_left', 'none', CLICK_ACTIONS[action])
    return action


def click_action():
    '''
    :return: what a plain single left click has been set to do
    '''
    return _action


def reapply(_self=cmd):
    '''
    Bind it again, after something else has rebound the buttons (picking
    a mouse mode rebinds them all)
    '''
    set_click_action(_action, _self)
