'''
Sequence Comparison plugin (pymolx): align the chains of related proteins
by sequence and color the residues by their similarity (pymolx.seqcompare).
'''


def seqcompare_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('Sequence Comparison', seqcompare_dialog)
