'''
Superposition/Alignment plugin (pymolx, parity item L-18)

Aligns many objects onto one ("extra_fit") or one selection onto another,
with any of PyMOL's alignment methods. Modeled on the Incentive PyMOL
dialog of the same name.
'''


def alignment_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('Alignment/Superposition', alignment_dialog)
