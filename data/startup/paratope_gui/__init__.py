'''
Paratope Prediction plugin (pymolx): which residues of an antibody bind
its antigen, with Parapred (sequence) or Paragraph (structure)
(pymolx.paratope).
'''


def paratope_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('Paratope Prediction', paratope_dialog)
