'''
Structure Prediction plugin (pymolx): sequences from a FASTA file, with
modified residues, glycans, ligands and disulfide bonds, written as the
input file of AlphaFold 3, RoseTTAFold3 or ESMFold2 (pymolx.prediction).
'''


def prediction_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('Structure Prediction', prediction_dialog)
