'''
Interface Analysis plugin (pymolx): PISA-like analysis of the interfaces
between the chains of an object, computed locally (pymolx.interfaces).
'''


def interface_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('Interface Analysis (PISA-like)', interface_dialog)
