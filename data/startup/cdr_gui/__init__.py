'''
CDR Loops plugin (pymolx): number an antibody's chains with ANARCI and
show its framework regions and CDR loops under IMGT, Kabat, Chothia,
Martin or AHo numbering, from a structure or a pasted sequence
(pymolx.cdr).
'''


def cdr_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('CDR Loops', cdr_dialog)
