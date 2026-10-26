'''
Structure Comparison plugin (pymolx): align closely related proteins and
color them by local structural similarity, blue (similar) to red
(different) (pymolx.compare).
'''


def compare_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('Structure Comparison', compare_dialog)
