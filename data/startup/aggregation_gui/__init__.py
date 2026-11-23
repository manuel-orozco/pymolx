'''
Aggregation Prediction plugin (pymolx): which residues drive a folded
protein's aggregation and how that changes with pH, a free
re-implementation of Aggrescan4D (pymolx.aggrescan).
'''


def aggregation_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('Aggregation Prediction', aggregation_dialog)
