'''
RMSD Clustering plugin (pymolx): cluster structures, e.g. predicted
antibody-antigen complexes, by RMSD (pymolx.cluster).
'''


def cluster_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('RMSD Clustering', cluster_dialog)
