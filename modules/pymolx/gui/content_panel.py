'''
Content panel (parity items L-14, L-15): Qt replacement for the OpenGL
object panel, docked right of the viewer, with the toggle toolbar
(mouse mode, wizards, sequence, timeline, command output) at its bottom.

Layout follows the Incentive PyMOL 3 main window. Rows show a dot (green
when enabled), the name (selections in parentheses, as in the OpenGL
panel) and A/S/H/L/C buttons. The menus behind the buttons are the ones
the OpenGL panel opens, from pymol.menu.
'''

from pymol.Qt import QtCore, QtGui, QtWidgets
from pymol import menu as pymol_menu

from .. import panel
from . import menus
from .toolbar import icon

Qt = QtCore.Qt

PANEL_WIDTH = 300
QWIDGETSIZE_MAX = (1 << 24) - 1
REFRESH_MS = 250

COLUMNS = 'ASHLC'

# get_type -> pymol.menu function names for the A, S, H, L, C buttons, as
# dispatched by the OpenGL panel (layer3/Executive.cpp). None: no menu.
_MOL = ('mol_action', 'mol_show', 'mol_hide', 'mol_labels', 'mol_color')
MENUS = {
    'all': ('all_action', 'mol_show', 'mol_hide', 'mol_labels', 'mol_color'),
    'selection': ('sele_action',) + _MOL[1:],
    'object:molecule': _MOL,
    'object:group': ('group_action',) + _MOL[1:],
    'object:map': ('map_action', 'map_show', 'map_hide', None,
                   'general_color'),
    'object:mesh': ('mesh_action', 'mesh_show', 'mesh_hide', None,
                    'mesh_color'),
    'object:surface': ('surface_action', 'surface_show', 'surface_hide', None,
                       'mesh_color'),
    'object:measurement': ('simple_action', 'measurement_show',
                           'measurement_hide', None, 'measurement_color'),
    'object:cgo': ('simple_action', 'cgo_show', 'cgo_hide', None,
                   'general_color'),
    'object:alignment': ('simple_action', 'cgo_show', 'cgo_hide', None,
                         'general_color'),
    'object:volume': ('simple_action', 'volume_show', 'volume_hide', None,
                      'vol_color'),
    'object:slice': ('slice_action', 'slice_show', 'slice_hide', None,
                     'slice_color'),
    'object:ramp': ('ramp_action', None, None, None, 'ramp_color'),
}
_DEFAULT_MENUS = ('simple_action', None, None, None, None)


def menu_data(item, column, _self):
    '''
    pymol.menu data for a panel row's button.

    :param item: pymolx.panel.PanelItem
    :param column: one of "ASHLC"
    :return: menu data list, or None if the button has no menu
    '''
    fname = MENUS.get(item.type, _DEFAULT_MENUS)[COLUMNS.index(column)]
    if fname is None:
        return None
    name = item.name
    if item.type == 'all' and column == 'L':
        name = '(all)'
    args = (_self, name)
    if item.type == 'object:surface' and column == 'C':
        args += ('surface',)
    return getattr(pymol_menu, fname)(*args)


class ElidedLabel(QtWidgets.QLabel):
    '''
    Label which shortens its text with "..." instead of growing, so long
    object names don't push the A/S/H/L/C buttons out of the panel.
    '''

    def __init__(self, text, parent=None):
        super().__init__(text, parent)
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored,
                           QtWidgets.QSizePolicy.Policy.Preferred)

    def minimumSizeHint(self):
        return QtCore.QSize(20, super().minimumSizeHint().height())

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.setPen(self.palette().color(self.foregroundRole()))
        rect = self.contentsRect()
        text = self.fontMetrics().elidedText(
            self.text(), Qt.TextElideMode.ElideRight, rect.width())
        painter.drawText(rect, int(Qt.AlignmentFlag.AlignLeft |
                                   Qt.AlignmentFlag.AlignVCenter), text)
        painter.end()


class PanelRow(QtWidgets.QFrame):
    '''
    One object panel row.
    '''

    def __init__(self, item, content_panel):
        super().__init__(content_panel.rows)
        self.item = item
        self.content_panel = content_panel
        self.cmd = content_panel.cmd
        self.setObjectName('panel_row')
        self.setProperty('enabled_item', item.enabled)
        self.setProperty('all_row', item.type == 'all')

        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(8 + 14 * item.nest_level, 5, 6, 5)
        layout.setSpacing(4)

        if item.is_group:
            arrow = QtWidgets.QToolButton(self)
            arrow.setObjectName('group_arrow')
            arrow.setText('▾' if item.is_open else '▸')
            arrow.setToolTip('Collapse group' if item.is_open else
                             'Expand group')
            arrow.clicked.connect(self.toggle_open)
            layout.addWidget(arrow)

        dot = QtWidgets.QLabel('●', self)
        dot.setObjectName('panel_dot')
        layout.addWidget(dot)

        text = item.name
        if item.type == 'selection':
            text = '(%s)' % text
        elif item.type == 'all':
            text = 'All'
        self.label = ElidedLabel(text, self)
        self.label.setObjectName('panel_name')
        self.label.setToolTip('%s\n\nClick to enable or disable, right-click '
                              'for actions' % text)
        layout.addWidget(self.label, 1)

        self.buttons = {}
        for column in COLUMNS:
            data = menu_data(item, column, self.cmd)
            btn = QtWidgets.QToolButton(self)
            btn.setObjectName('ashlc')
            btn.setProperty('column', column)
            btn.setText(column)
            if data is None:
                btn.setEnabled(False)
            else:
                btn.setMenu(self._menu(column))
                btn.setPopupMode(
                    QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
            layout.addWidget(btn)
            self.buttons[column] = btn

    def _menu(self, column):
        menu = menus.PyMenu('', self)

        @menu.aboutToShow.connect
        def _():
            menu.clear()
            menus.fill_menu(menu, menu_data(self.item, column, self.cmd),
                            self.cmd)

        return menu

    def toggle_enabled(self):
        if self.item.enabled:
            self.cmd.disable(self.item.name)
        else:
            self.cmd.enable(self.item.name)
        self.content_panel.refresh()

    def toggle_open(self):
        self.cmd.group(self.item.name,
                       action='close' if self.item.is_open else 'open')
        self.content_panel.refresh()

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.toggle_enabled()
        elif event.button() == Qt.MouseButton.RightButton:
            menu = self.buttons['A'].menu()
            if menu is not None:
                pos = (event.globalPosition().toPoint()
                       if hasattr(event, 'globalPosition') else
                       event.globalPos())
                menu.popup(pos)
        else:
            super().mousePressEvent(event)


class WizardPanel(QtWidgets.QFrame):
    '''
    Controls of the active wizard (e.g. mutagenesis), from its get_panel():
    [1, title, ''], [2, button label, command], [3, menu label, menu tag].
    The OpenGL panel used to draw these.
    '''

    def __init__(self, cmd, parent=None):
        super().__init__(parent)
        self.cmd = cmd
        self.setObjectName('wizard_panel')
        self.setLayout(QtWidgets.QVBoxLayout())
        self.layout().setContentsMargins(6, 6, 6, 6)
        self.layout().setSpacing(3)
        self._signature = None
        self.hide()

    def _entries(self, wizard):
        try:
            return [tuple(entry[:3]) for entry in wizard.get_panel() or ()]
        except Exception as e:
            print(' Wizard panel error:', e)
            return []

    def refresh(self):
        wizard = self.cmd.get_wizard()
        entries = self._entries(wizard) if wizard is not None else []
        signature = (id(wizard), entries)
        if signature == self._signature:
            return
        self._signature = signature

        layout = self.layout()
        while layout.count():
            widget = layout.takeAt(0).widget()
            if widget is not None:
                widget.deleteLater()

        for kind, label, arg in entries:
            label = menus.strip_label(label)
            if kind == 1:
                title = QtWidgets.QLabel(label, self)
                title.setObjectName('wizard_title')
                layout.addWidget(title)
            elif kind == 2:
                button = QtWidgets.QPushButton(label, self)
                button.setObjectName('wizard_button')
                button.clicked.connect(
                    lambda _=False, c=arg: self.cmd.do(c, echo=0))
                layout.addWidget(button)
            elif kind == 3:
                # tool button: push buttons with menus ignore text padding
                button = QtWidgets.QToolButton(self)
                button.setObjectName('wizard_menu')
                button.setText(label)
                button.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                                     QtWidgets.QSizePolicy.Policy.Fixed)
                button.setMenu(self._menu(wizard, arg))
                button.setPopupMode(
                    QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
                layout.addWidget(button)
        self.setVisible(bool(entries))

    def _menu(self, wizard, tag):
        menu = menus.PyMenu('', self)

        @menu.aboutToShow.connect
        def _():
            menu.clear()
            menus.fill_menu(menu, wizard.get_menu(tag) or [], self.cmd)

        return menu


class StateBar(QtWidgets.QFrame):
    '''
    State controls (first, back, play, forward, last) and "State n/N",
    shown when there is more than one state, e.g. the rotamers of the
    mutagenesis wizard. The OpenGL panel used to have these buttons.
    '''

    def __init__(self, cmd, parent=None):
        super().__init__(parent)
        self.cmd = cmd
        self.setObjectName('state_bar')
        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(6, 2, 6, 2)
        layout.setSpacing(2)

        self.buttons = {}
        for name, tooltip, callback in [
            ('first', 'First state', cmd.rewind),
            ('prev', 'Previous state', cmd.backward),
            ('play', 'Play', None),
            ('next', 'Next state', cmd.forward),
            ('last', 'Last state', cmd.ending),
        ]:
            button = QtWidgets.QToolButton(self)
            button.setObjectName('state_button')
            button.setIcon(icon(cmd, 'state-' + name))
            button.setIconSize(QtCore.QSize(12, 12))
            button.setToolTip(tooltip)
            if callback is not None:
                button.clicked.connect(lambda _=False, f=callback: f())
            layout.addWidget(button)
            self.buttons[name] = button
        self.buttons['play'].setCheckable(True)
        self.buttons['play'].clicked.connect(self._play)

        self.label = QtWidgets.QLabel(self)
        self.label.setObjectName('state_label')
        layout.addWidget(self.label, 1)
        self.hide()

    def _play(self, checked):
        if checked:
            self.cmd.mplay()
        else:
            self.cmd.mstop()

    def refresh(self):
        cmd = self.cmd
        count = cmd.count_states('all')
        self.setVisible(count > 1)
        if count <= 1:
            return
        state = cmd.get_state()
        text = 'State %d/%d' % (state, count)
        # mutagenesis wizard: strain of the current rotamer
        scores = getattr(cmd.get_wizard(), 'bump_scores', None)
        if scores and 0 < state <= len(scores):
            text += '   strain %.1f' % scores[state - 1]
        self.label.setText(text)
        playing = bool(cmd.get_movie_playing())
        button = self.buttons['play']
        button.blockSignals(True)
        button.setChecked(playing)
        button.setIcon(icon(cmd, 'state-pause' if playing else 'state-play'))
        button.setToolTip('Stop' if playing else 'Play')
        button.blockSignals(False)


class ToggleToolbar(QtWidgets.QToolBar):
    '''
    Bottom strip: mouse mode, wizards, sequence viewer, movie timeline,
    command output.
    '''

    def __init__(self, window):
        super().__init__('Toggles', window)
        self.window = window
        self.cmd = cmd = window.cmd
        self.setObjectName('pymolx_toggles')
        self.setIconSize(QtCore.QSize(16, 16))

        self.mouse_button = self._button('', 'mouse',
                                         'Mouse mode', self._mouse_menu())
        self.mouse_button.setToolButtonStyle(
            Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.mouse_button.setProperty('dropdown', True)
        window.setting_callbacks[cmd.setting._get_index(
            'button_mode_name')].append(self._update_mouse_mode)
        self._update_mouse_mode(cmd.get('button_mode_name'))

        spacer = QtWidgets.QWidget(self)
        spacer.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                             QtWidgets.QSizePolicy.Policy.Preferred)
        self.addWidget(spacer)

        self._button('', 'wand', 'Wizards', window.menudict.get('Wizard'))

        self.seq_button = self._setting_toggle('SEQ', None, 'seq_view',
                                               'Sequence viewer')
        self.timeline_button = self._setting_toggle(
            '', 'timeline', 'movie_panel', 'Movie timeline')

        self.command_button = self._button('', 'terminal',
                                           'Show or hide the command output')
        self.command_button.setCheckable(True)
        self.command_button.setChecked(not window.browser.isHidden())
        self.command_button.toggled.connect(self.toggle_output)

    def _button(self, text, icon_name, tooltip, menu=None):
        btn = QtWidgets.QToolButton(self)
        btn.setText(text)
        if icon_name:
            btn.setIcon(icon(self.cmd, icon_name))
        btn.setToolTip(tooltip)
        if menu is not None:
            btn.setMenu(menu)
            btn.setPopupMode(
                QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
        self.addWidget(btn)
        return btn

    def _setting_toggle(self, text, icon_name, setting, tooltip):
        cmd = self.cmd
        btn = self._button(text, icon_name, tooltip)
        btn.setCheckable(True)
        btn.clicked.connect(lambda checked: cmd.set(setting, int(checked),
                                                    log=1, quiet=0))

        def update(value):
            btn.blockSignals(True)
            btn.setChecked(bool(value))
            btn.blockSignals(False)

        self.window.setting_callbacks[cmd.setting._get_index(setting)].append(
            update)
        update(cmd.get_setting_int(setting))
        return btn

    def _mouse_menu(self):
        menu = menus.PyMenu('', self)

        @menu.aboutToShow.connect
        def _():
            menu.clear()
            menus.fill_menu(menu, pymol_menu.mouse_config(self.cmd), self.cmd)

        return menu

    def _update_mouse_mode(self, name):
        self.mouse_button.setText(name or 'Mouse')

    def toggle_output(self, visible):
        '''
        Show or hide the output pane; the command line stays.
        '''
        window = self.window
        window.browser.setVisible(visible)
        dock = window.ext_window
        frame = dock.widget()
        frame.layout().activate()
        if visible:
            dock.setMaximumHeight(QWIDGETSIZE_MAX)
            window.resizeDocks([dock], [frame.sizeHint().height()],
                               Qt.Orientation.Vertical)
        else:
            # just the command line, no empty space above it
            dock.setMaximumHeight(frame.minimumSizeHint().height())


class ContentPanel(QtWidgets.QWidget):
    '''
    Object rows (refreshed from pymolx.panel) plus the toggle toolbar.
    '''

    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.cmd = window.cmd
        self.setObjectName('content_panel')
        self._items = None

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self.rows = QtWidgets.QWidget()
        self.rows.setObjectName('panel_rows')
        self.rows_layout = QtWidgets.QVBoxLayout(self.rows)
        self.rows_layout.setContentsMargins(0, 4, 0, 4)
        self.rows_layout.setSpacing(0)
        self.rows_layout.addStretch()

        scroll = QtWidgets.QScrollArea(self)
        scroll.setObjectName('panel_scroll')
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(self.rows)
        layout.addWidget(scroll, 1)

        self.wizard_panel = WizardPanel(self.cmd, self)
        layout.addWidget(self.wizard_panel)
        self.state_bar = StateBar(self.cmd, self)
        layout.addWidget(self.state_bar)

        self.toggles = ToggleToolbar(window)
        layout.addWidget(self.toggles)

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.refresh)
        self.timer.start(REFRESH_MS)
        self.refresh()

    def row_widgets(self):
        return [self.rows_layout.itemAt(i).widget()
                for i in range(self.rows_layout.count() - 1)]

    def refresh(self, force=False):
        '''
        Rebuild the rows if the panel list changed. Skipped while a menu
        is open, so its button isn't deleted under it.
        '''
        if QtWidgets.QApplication.activePopupWidget() is not None:
            return
        self.wizard_panel.refresh()
        self.state_bar.refresh()
        items = panel.get_panel_list(_self=self.cmd)
        if items == self._items and not force:
            return
        self._items = items

        for widget in self.row_widgets():
            self.rows_layout.removeWidget(widget)
            widget.deleteLater()
        for i, item in enumerate(items):
            self.rows_layout.insertWidget(i, PanelRow(item, self))


def setup(window):
    '''
    Replace the OpenGL object panel and its one-line command prompt with
    the content panel, docked right of the viewer above the command line.
    '''
    cmd = window.cmd
    cmd.set('internal_gui', 0)
    cmd.set('internal_feedback', 0)

    content = ContentPanel(window)
    dock = QtWidgets.QDockWidget('Content', window)
    dock.setObjectName('content_dock')
    dock.setTitleBarWidget(QtWidgets.QWidget())
    dock.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
    dock.setWidget(content)

    # bottom docks (command line) span the full window width
    window.setCorner(Qt.Corner.BottomLeftCorner,
                     Qt.DockWidgetArea.BottomDockWidgetArea)
    window.setCorner(Qt.Corner.BottomRightCorner,
                     Qt.DockWidgetArea.BottomDockWidgetArea)
    window.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
    window.resizeDocks([dock], [PANEL_WIDTH], Qt.Orientation.Horizontal)

    # the window was sized for the OpenGL panel (220 pixels, if shown)
    import pymol.invocation
    reserved = 220 if pymol.invocation.options.internal_gui else 0
    window.resize(window.width() + PANEL_WIDTH - reserved, window.height())
    return content


def forward_key_to_command_line(window, event):
    '''
    Typing while the viewer has focus goes to the command line (there is
    no command line in the viewer). Returns True if the key was taken.
    '''
    text = event.text()
    modifiers = event.modifiers() & (Qt.KeyboardModifier.ControlModifier |
                                     Qt.KeyboardModifier.AltModifier |
                                     Qt.KeyboardModifier.MetaModifier)
    if not text or not text.isprintable() or modifiers:
        return False
    window.lineedit.setFocus()
    window.lineedit.insert(text)
    return True
