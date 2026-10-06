'''
Main window toolbar (parity item L-03), replacing the upstream grid of
"quick buttons". Layout follows the Incentive PyMOL 3 main window:

    [Residues v] [undo] [redo] | [Zoom v] Orient Rock [Presets...]
                                  ... Builder... Scenes [Draw/Ray v] [...]

What the Zoom and overflow menus contain isn't visible in the reference
screenshot; their contents are our choice (see INCENTIVE_PARITY.md).
'''

from pymol.Qt import QtCore, QtGui, QtWidgets
from pymol.Qt.utils import WidgetMenu
from pymol import menu as pymol_menu

from .. import undo
from . import menus, theme

Qt = QtCore.Qt

SELECTION_MODES = [
    ('Atoms', 0),
    ('Residues', 1),
    ('Chains', 2),
    ('Segments', 3),
    ('Objects', 4),
    ('Molecules', 5),
    ('C-alphas', 6),
]

def icon(_self, name):
    '''
    Icon from the pymolx icon directory. If "<name>-on.svg" exists, it is
    used for the checked state of checkable buttons.
    '''
    import os
    path = _self.exp_path('%s/%s.svg' % (theme.ICON_DIR, name))
    result = QtGui.QIcon(path)
    path_on = path[:-4] + '-on.svg'
    if os.path.exists(path_on):
        result.addFile(path_on, QtCore.QSize(), QtGui.QIcon.Mode.Normal,
                       QtGui.QIcon.State.On)
    return result


class Toolbar(QtWidgets.QToolBar):
    '''
    :param window: pmg_qt.pymol_qt_gui.PyMOLQtGUI
    '''

    def __init__(self, window):
        super().__init__('Toolbar', window)
        self.window = window
        self.cmd = cmd = window.cmd

        self.setObjectName('pymolx_toolbar')
        self.setMovable(False)
        self.setIconSize(QtCore.QSize(16, 16))

        # like Incentive PyMOL, 1px dividers separate the buttons
        self._add_selection_mode()
        self.addSeparator()
        self._add_undo_redo()
        self.addSeparator()
        self.add_button('Zoom', lambda: cmd.zoom(animate=1.0),
                        menu=self._zoom_menu(),
                        popup=QtWidgets.QToolButton.ToolButtonPopupMode.MenuButtonPopup,
                        tooltip='Zoom on all objects')
        self.addSeparator()
        self.add_button('Orient', lambda: cmd.orient(animate=1.0),
                        tooltip='Orient on all objects')
        self.addSeparator()
        self._add_rock()
        self.addSeparator()
        presets = self.add_button('Presets...', menu=self._dynamic_menu(
            lambda: pymol_menu.presets(cmd, 'all')),
            tooltip='Apply a representation preset to all objects')
        presets.setProperty('dropdown', False)
        presets.setProperty('indicator', False)

        spacer = QtWidgets.QWidget(self)
        spacer.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                             QtWidgets.QSizePolicy.Policy.Preferred)
        self.addWidget(spacer)

        self.add_button('Builder...', window.open_builder_panel)
        self.addSeparator()
        self.add_button('Scenes', window.scene_panel_menu_dialog)
        self.addSeparator()
        # render dialog is constructed when the menu is first shown
        self.add_button('Draw/Ray', icon_name='camera',
                        menu=WidgetMenu(window).setSetupUi(window.render_dialog),
                        tooltip='Draw or ray trace an image')
        self.addSeparator()
        more = self.add_button('\u2022\u2022\u2022', menu=self._more_menu(),
                               tooltip='More actions')
        more.setObjectName('more_button')
        more.setProperty('dropdown', False)
        more.setProperty('indicator', False)

    def add_button(self, text, callback=None, *, icon_name=None, tooltip=None,
                   menu=None,
                   popup=QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup):
        btn = QtWidgets.QToolButton(self)
        btn.setText(text)
        if icon_name:
            btn.setIcon(icon(self.cmd, icon_name))
        btn.setToolButtonStyle(
            Qt.ToolButtonStyle.ToolButtonTextBesideIcon if text else
            Qt.ToolButtonStyle.ToolButtonIconOnly)
        if tooltip:
            btn.setToolTip(tooltip)
        if callback is not None:
            btn.clicked.connect(lambda _checked=False: callback())
        if menu is not None:
            btn.setMenu(menu)
            btn.setPopupMode(popup)
            # makes room for the arrow (split buttons have their own)
            btn.setProperty('dropdown', popup ==
                QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
        self.addWidget(btn)
        return btn

    def _on_setting(self, name, callback):
        '''
        Call callback(value) now and whenever setting "name" changes.
        '''
        index = self.cmd.setting._get_index(name)
        self.window.setting_callbacks[index].append(callback)
        callback(self.cmd.get_setting_tuple(name)[1][0])

    def _add_selection_mode(self):
        cmd = self.cmd
        menu = QtWidgets.QMenu(self)
        group = QtWidgets.QActionGroup(menu)
        actions = {}
        for label, value in SELECTION_MODES:
            action = menu.addAction(label, lambda v=value: cmd.set(
                'mouse_selection_mode', v, log=1, quiet=0))
            action.setCheckable(True)
            group.addAction(action)
            actions[value] = action

        btn = self.add_button('', menu=menu, icon_name='pointer',
                              tooltip='Selection mode for mouse clicks')
        btn.setProperty('accent', True)
        btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.selection_mode_button = btn

        def update(value):
            action = actions.get(value)
            if action is not None:
                action.setChecked(True)
                btn.setText(action.text())

        self._on_setting('mouse_selection_mode', update)

    def _add_undo_redo(self):
        '''
        Undo/redo buttons, enabled when there is something to undo/redo.
        '''
        cmd = self.cmd
        self.undo_button = self.add_button(
            '', lambda: cmd.undo(), icon_name='undo', tooltip='Undo')
        self.redo_button = self.add_button(
            '', lambda: cmd.redo(), icon_name='redo', tooltip='Redo')

        def update():
            self.undo_button.setEnabled(undo.stack.can_undo())
            self.redo_button.setEnabled(undo.stack.can_redo())

        self.undo_timer = QtCore.QTimer(self)
        self.undo_timer.timeout.connect(update)
        self.undo_timer.start(250)
        update()

    def _add_rock(self):
        cmd = self.cmd
        btn = self.add_button('Rock', tooltip='Toggle rocking about the Y axis')
        btn.setCheckable(True)
        btn.clicked.connect(lambda checked: cmd.rock(int(checked)))
        self.rock_button = btn

        def update(value):
            btn.blockSignals(True)
            btn.setChecked(bool(value))
            btn.blockSignals(False)

        self._on_setting('rock', update)

    def _dynamic_menu(self, get_data):
        '''
        Menu which is rebuilt from pymol.menu data each time it is shown.
        '''
        menu = menus.PyMenu('', self)

        @menu.aboutToShow.connect
        def _():
            menu.clear()
            menus.fill_menu(menu, get_data(), self.cmd)

        return menu

    def _zoom_menu(self):
        cmd = self.cmd
        menu = QtWidgets.QMenu(self)
        menu.addAction('All', lambda: cmd.zoom(animate=1.0))
        sele = menu.addAction('Selection',
                              lambda: cmd.zoom('sele', animate=1.0))
        menu.addAction('Center', lambda: cmd.center(animate=1.0))
        menu.addSeparator()
        menu.addAction('Reset View', cmd.reset)
        menu.aboutToShow.connect(
            lambda: sele.setEnabled('sele' in cmd.get_names('selections')))
        return menu

    def _more_menu(self):
        cmd = self.cmd
        window = self.window
        menu = QtWidgets.QMenu(self)
        menu.addAction('Get View', window.get_view)
        menu.addAction('Unpick', cmd.unpick)
        menu.addAction('Deselect', cmd.deselect)
        menu.addSeparator()
        menu.addAction('Properties...', window.open_props_dialog)
        menu.addAction('Rebuild', cmd.rebuild)
        menu.addSeparator()
        movie = QtWidgets.QMenu('Movie', menu)  # parent keeps it alive
        menu.addMenu(movie)
        for label, callback in [
            ('Rewind', cmd.rewind),
            ('Backward', cmd.backward),
            ('Stop', cmd.mstop),
            ('Play', cmd.mplay),
            ('Forward', cmd.forward),
            ('End', cmd.ending),
            (None, None),
            ('Clear Movie', cmd.mclear),
        ]:
            if label is None:
                movie.addSeparator()
            else:
                movie.addAction(label, callback)
        return menu


def setup(window):
    '''
    Add the toolbar to the main window and grow the window by its height,
    so the viewport keeps its requested size.
    '''
    toolbar = Toolbar(window)
    window.addToolBar(Qt.ToolBarArea.TopToolBarArea, toolbar)
    window.resize(window.width(),
                  window.height() + toolbar.sizeHint().height())
    return toolbar
