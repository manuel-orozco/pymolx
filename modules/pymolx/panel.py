'''
Rows of the object panel, as the OpenGL panel shows them, for GUIs which
draw their own panel (parity item L-15).
'''

from collections import namedtuple

from pymol import cmd, _cmd

PanelItem = namedtuple('PanelItem',
                       'name type enabled nest_level is_group is_open')
PanelItem.__doc__ = '''
One row: "type" is "all", "selection" or "object:..." like cmd.get_type.
Members of closed groups and hidden ("_") names are not listed.
'''


def get_panel_list(*, _self=cmd):
    '''
    :return: list of PanelItem, in panel order
    '''
    with _self.lockcm:
        rows = _cmd.get_panel_list(_self._COb)
    return [PanelItem(name, type_, bool(enabled), nest_level, bool(is_group),
                      bool(is_open))
            for (name, type_, enabled, nest_level, is_group, is_open) in rows]
