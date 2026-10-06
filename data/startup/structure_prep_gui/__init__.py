'''
Structure Preparation plugin (pymolx): fix structures with PDBFixer and
energy-minimize with OpenMM (pymolx.mm).
'''


def structure_prep_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('Structure Preparation', structure_prep_dialog)
