'''
Antibody Builder plugin (pymolx): a full IgG from two Fabs, the hinge
and the Fc, placed on a full-length template (pymolx.antibody).
'''


def antibody_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('Antibody Builder', antibody_dialog)
