'''
Interface RMSD Clustering plugin (pymolx): cluster complexes, e.g. the
antibody-antigen models of an epitope prediction, by the RMSD of their
interface residues (pymolx.irmsd).
'''


def irmsd_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('Interface RMSD Clustering', irmsd_dialog)
