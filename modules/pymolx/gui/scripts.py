'''
"Scripts" in the menu bar: like File > Run Script..., with the file
dialog starting in the scripts folder (pymolx.paths).
'''

from pymol.Qt import QtWidgets

from .. import paths


def run_script(window):
    '''
    File > Run Script... (window.file_run), starting in the scripts
    folder. If no script is chosen, the dialogs' folder stays as it was.
    '''
    previous = getattr(window, 'initialdir', None)
    window.initialdir = paths.scripts_dir(create=True)
    window.file_run()
    if window.initialdir == paths.scripts_dir() and previous is not None:
        # cancelled (file_run moves it to the chosen script's folder)
        window.initialdir = previous


def setup(window):
    '''
    Add "Scripts" to the menu bar, before Help.
    '''
    menubar = window.menuBar()
    action = QtWidgets.QAction('Scripts', menubar)
    action.setObjectName('scripts_action')
    action.setToolTip('Run a script (like File > Run Script...)')
    action.triggered.connect(lambda _=False: run_script(window))
    help_menu = getattr(window, 'menudict', {}).get('Help')
    if help_menu is not None:
        menubar.insertAction(help_menu.menuAction(), action)
    else:
        menubar.addAction(action)
    return action
