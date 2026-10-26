'''
Binding Pocket plugin (pymolx): the residues within a distance of a
ligand, as a selection and as a residue list (pymolx.pocket).
'''


def pocket_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('Binding Pocket', pocket_dialog)
