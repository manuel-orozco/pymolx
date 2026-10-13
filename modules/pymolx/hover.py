'''
Hover highlighting (parity item L-22, like Schrodinger Maestro): the atom,
residue, chain, ... under the mouse pointer, following
mouse_selection_mode, is drawn with lighter colors where it is shown as
lines or sticks. Same style and thickness: the core lightens the colors
while building lines and sticks (SceneHoverTint); the atoms' colors
don't change, so nothing reaches sessions, saved files or undo.

The Qt part (pymolx.gui.hover) finds the atom under the pointer with
_cmd.pick_at.
'''

from pymol import cmd, _cmd

from .sequence import SELE_MODE_KEYWORDS, uncounted


def unit_selection(model, index, mode=None, *, _self=cmd):
    '''
    The atom's unit to highlight, like a click selects it.

    :param model: object name
    :param index: atom index (1-based)
    :param mode: mouse_selection_mode (default: the current setting)
    :return: selection expression
    '''
    if mode is None:
        mode = _self.get_setting_int('mouse_selection_mode')
    atom = '(%%%s and index %d)' % (model, index)
    keyword = SELE_MODE_KEYWORDS[mode % len(SELE_MODE_KEYWORDS)]
    return '%s%s' % (keyword, atom) if keyword else atom


def unit_key(model, index, mode=None, *, _self=cmd):
    '''
    Identifies the atom's unit cheaply, so moving the pointer within the
    same residue (chain, ...) changes nothing.
    '''
    if mode is None:
        mode = _self.get_setting_int('mouse_selection_mode')
    mode %= len(SELE_MODE_KEYWORDS)
    if mode == 0:  # atoms
        return (model, index)
    fields = []
    _self.iterate('%%%s and index %d' % (model, index),
                  'fields.append((segi, chain, resi))',
                  space={'fields': fields})
    if not fields:
        return None
    segi, chain, resi = fields[0]
    return {
        1: (model, segi, chain, resi),   # residues
        2: (model, chain),               # chains
        3: (model, segi),                # segments
        4: (model,),                     # objects
        6: (model, segi, chain, resi),   # C-alphas (of the residue)
    }.get(mode) or (model, index, mode)  # molecules: per atom


def show(selection, *, _self=cmd):
    '''
    Highlight the selected atoms (of one object), replacing any previous
    highlight.

    :return: True if there are atoms
    '''
    # (the query's temporary selection is no news for the sequence viewer)
    with uncounted(_self=_self):
        indices = {}
        _self.iterate(selection,
                      'indices.setdefault(model, []).append(index)',
                      space={'indices': indices})
        if not indices:
            clear(_self=_self)
            return False
        model, atoms = next(iter(indices.items()))
        with _self.lockcm:
            _cmd.set_hover_atoms(_self._COb, model, atoms)
    return True


def clear(*, _self=cmd):
    '''
    No highlight
    '''
    with uncounted(_self=_self):
        with _self.lockcm:
            _cmd.set_hover_atoms(_self._COb, '', [])
