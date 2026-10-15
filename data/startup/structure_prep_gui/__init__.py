'''
Structure Preparation plugin (pymolx): a protein preparation workflow
like Maestro's Protein Preparation Wizard (pymolx.prep: PDBFixer, PROPKA,
OpenMM) and structure diagnostics (pymolx.diagnostics).
'''


def structure_prep_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('Structure Preparation', structure_prep_dialog)
