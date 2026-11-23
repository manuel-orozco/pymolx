'''
Modify Sequences plugin (pymolx): renumber residues, change chain IDs
and segment IDs, merge chains into one chain or into one object
(pymolx.modify).
'''


def modify_dialog():
    from . import main
    return main.show_dialog()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('Modify Sequences', modify_dialog)
