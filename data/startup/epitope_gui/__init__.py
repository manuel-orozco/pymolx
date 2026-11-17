'''
Epitope Prediction plugin (pymolx): an antibody and its antigen folded
together by AlphaFold 3, to see where the antibody binds. Structure
Prediction's dialog with the settings that matter for a complex: each
chain's role, many seeds, samples per seed, and optionally no MSA for the
antibody chains (pymolx.epitope).
'''


def epitope_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('Epitope Prediction', epitope_dialog)
