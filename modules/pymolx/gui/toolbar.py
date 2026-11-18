'''
Main window toolbar (parity item L-03), replacing the upstream grid of
"quick buttons". The left half is Maestro-like (parity item L-33), the
rest follows the Incentive PyMOL 3 main window:

    [Select] [+/-] [R Residues v] | [Undo] [Redo] |
        [Zoom All] [Zoom Sel] | [Presets v]
        ... Builder... Scenes [Draw/Ray v] [...]

Every button of the left half is an icon (or the selection level's
letter) over its name in small capitals, drawn by CaptionButton, so each
one fills the bar's height.

The view actions without a button of their own (Orient, Center, Reset
View, Rock) are in the "..." menu. What that menu contains isn't visible
in the reference screenshot; its contents are our choice (see
INCENTIVE_PARITY.md).
'''

from pymol.Qt import QtCore, QtGui, QtWidgets
from pymol.Qt.utils import WidgetMenu
from pymol import menu as pymol_menu

from .. import picking, undo
from . import menus, theme

Qt = QtCore.Qt

ICON = 20       # px: icon buttons (undo, zoom, ...), as tall as the
                # letter of the drawn buttons next to them

# what a click in the viewer picks ("mouse_selection_mode"), with the
# short letter Maestro's selector shows for each
SELECTION_MODES = [
    ('Atoms', 0, 'A'),
    ('Residues', 1, 'R'),
    ('Chains', 2, 'C'),
    ('Segments', 3, 'S'),
    ('Objects', 4, 'O'),
    ('Molecules', 5, 'M'),
    ('C-alphas', 6, 'CA'),
]


class CaptionButton(QtWidgets.QToolButton):
    '''
    A toolbar button drawn as Maestro draws them: a letter or an icon
    above, the name in small capitals underneath, and a menu.

    (Painted rather than laid out so the whole thing is one button, with
    the frame, hover and pressed states the style gives it.)
    '''

    LETTER_SCALE = 1.8
    CAPTION_SCALE = 0.9
    GLYPH = 22          # px: the icon's box, when there is one

    def __init__(self, parent=None):
        super().__init__(parent)
        self.letter = ''
        self.caption = ''
        self.glyph = None       # QIcon, drawn in place of the letter
        self.setPopupMode(
            QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed,
                           QtWidgets.QSizePolicy.Policy.Expanding)
        self.setProperty('caption', True)

    def set_mode(self, letter, caption):
        self.letter = letter
        self.caption = caption.upper()
        self.updateGeometry()
        self.update()

    def set_glyph(self, glyph, caption):
        self.glyph = glyph
        self.letter = ''
        self.caption = caption.upper()
        self.updateGeometry()
        self.update()

    def _top_height(self):
        # the same for a letter and an icon, so the buttons line up
        height = QtGui.QFontMetrics(self._fonts()[0]).height()
        return max(height, self.GLYPH) if self.glyph is not None \
            else height

    def _fonts(self):
        letter = QtGui.QFont(self.font())
        letter.setPointSizeF(max(1.0, self.font().pointSizeF() *
                                 self.LETTER_SCALE))
        letter.setBold(True)
        caption = QtGui.QFont(self.font())
        caption.setPointSizeF(max(1.0, self.font().pointSizeF() *
                                  self.CAPTION_SCALE))
        return letter, caption

    def sizeHint(self):
        letter_font, caption_font = self._fonts()
        letter = QtGui.QFontMetrics(letter_font)
        caption = QtGui.QFontMetrics(caption_font)
        top = self.GLYPH if self.glyph is not None else \
            letter.horizontalAdvance(self.letter or 'CA')
        width = max(top + 22,
                    caption.horizontalAdvance(self.caption or 'MOLECULES')
                    + 10)
        height = self._top_height() + caption.height() + 8
        return QtCore.QSize(width, height)

    def paintEvent(self, event):
        painter = QtWidgets.QStylePainter(self)
        option = QtWidgets.QStyleOptionToolButton()
        self.initStyleOption(option)
        # the frame, the hover and the menu arrow, but not the label:
        # the letter and the caption are drawn here instead
        option.text = ''
        option.icon = QtGui.QIcon()
        painter.drawComplexControl(
            QtWidgets.QStyle.ComplexControl.CC_ToolButton, option)
        letter_font, caption_font = self._fonts()
        captions = QtGui.QFontMetrics(caption_font)
        inner = self.rect().adjusted(2, 3, -2, -3)
        top = QtCore.QRect(inner.left(), inner.top(), inner.width(),
                           self._top_height())
        bottom = QtCore.QRect(inner.left(),
                              inner.bottom() - captions.height(),
                              inner.width(), captions.height())
        on = self.isEnabled()
        if self.glyph is not None:
            self.glyph.paint(painter, top, Qt.AlignmentFlag.AlignCenter,
                             QtGui.QIcon.Mode.Normal if on else
                             QtGui.QIcon.Mode.Disabled)
        else:
            painter.setFont(letter_font)
            painter.setPen(QtGui.QColor(theme.COLORS[
                'text_bright' if on else 'text_dim']))
            painter.drawText(top, int(Qt.AlignmentFlag.AlignCenter),
                             self.letter)
        painter.setFont(caption_font)
        painter.setPen(QtGui.QColor(theme.COLORS['text' if on else
                                                 'text_dim']))
        painter.drawText(bottom, int(Qt.AlignmentFlag.AlignCenter),
                         self.caption)

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
        self.setIconSize(QtCore.QSize(ICON, ICON))

        # like Incentive PyMOL, 1px dividers separate the buttons
        self._add_select_tool()
        self._add_selection_mode()
        self.addSeparator()
        self._add_undo_redo()
        self.addSeparator()
        self._add_zoom()
        self.addSeparator()
        self._add_presets()

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
        btn.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed,
                          QtWidgets.QSizePolicy.Policy.Expanding)
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

    def caption_button(self, name, icon_name, caption, tooltip,
                       callback=None, menu=None, checkable=False):
        """
        An icon with its name in small capitals underneath, like the
        selection level and Presets, so every button fills the bar
        """
        btn = CaptionButton(self)
        btn.setObjectName(name)
        btn.set_glyph(icon(self.cmd, icon_name), caption)
        btn.setToolTip(tooltip)
        if checkable:
            btn.setCheckable(True)
        if callback is not None:
            btn.clicked.connect(lambda _checked=False: callback())
        if menu is not None:
            btn.setMenu(menu)
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
        """
        What a click in the viewer picks, shown as Maestro shows it: the
        level's letter with its name underneath, and a menu to change it
        """
        cmd = self.cmd
        menu = QtWidgets.QMenu(self)
        group = QtWidgets.QActionGroup(menu)
        modes, actions = {}, {}
        for label, value, letter in SELECTION_MODES:
            action = menu.addAction(
                '%s  (%s)' % (label, letter),
                lambda v=value: cmd.set('mouse_selection_mode', v, log=1,
                                        quiet=0))
            action.setCheckable(True)
            group.addAction(action)
            actions[value] = action
            modes[value] = (letter, label)

        btn = CaptionButton(self)
        btn.setObjectName('selection_mode_button')
        btn.setMenu(menu)
        btn.setToolTip('What a click in the viewer selects')
        btn.setProperty('accent', True)
        self.addWidget(btn)
        self.selection_mode_button = btn

        def update(value):
            action = actions.get(value)
            if action is not None:
                action.setChecked(True)
                btn.set_mode(*modes[value])

        self._on_setting('mouse_selection_mode', update)

    def _add_select_tool(self):
        '''
        The selecting tool: the arrow turns picking on, + and - make a
        click add to or take from the selection instead of replacing it
        (Maestro's). Nothing pressed: a click replaces the selection.
        '''
        self.select_button = btn = self.caption_button(
            'select_button', 'select-arrow', 'Select',
            'Click in the viewer to select', checkable=True)
        btn.setChecked(True)
        btn.toggled.connect(self._picking_changed)

        # + over -, flush, as one segmented control next to the arrow
        column = QtWidgets.QWidget(self)
        column.setObjectName('select_pm')
        column.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed,
                             QtWidgets.QSizePolicy.Policy.Expanding)
        layout = QtWidgets.QVBoxLayout(column)
        layout.setContentsMargins(0, 5, 6, 5)
        layout.setSpacing(0)
        self.add_button_plus = self._small_toggle('+', 'add',
                                                  'Add what you click to '
                                                  'the selection')
        self.add_button_minus = self._small_toggle('\u2212', 'subtract',
                                                   'Take what you click '
                                                   'out of the selection')
        layout.addWidget(self.add_button_plus)
        layout.addWidget(self.add_button_minus)
        self.addWidget(column)

        # picking a mouse mode rebinds every button, so the tool's choice
        # is put back then (and bound in the first place now)
        self._on_setting('button_mode_name',
                         lambda _value: self._picking_changed())

    def _small_toggle(self, text, action, tooltip):
        btn = QtWidgets.QToolButton(self)
        btn.setObjectName('select_%s' % action)
        btn.setText(text)
        btn.setCheckable(True)
        btn.setToolTip(tooltip)
        btn.setProperty('small', True)
        btn.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed,
                          QtWidgets.QSizePolicy.Policy.Expanding)
        btn.toggled.connect(
            lambda on, a=action: self._mode_toggled(a, on))
        return btn

    def _mode_toggled(self, action, on):
        '''
        + and - are alternatives: turning one on turns the other off
        '''
        other = self.add_button_minus if action == 'add' else \
            self.add_button_plus
        if on and other.isChecked():
            other.blockSignals(True)
            other.setChecked(False)
            other.blockSignals(False)
        self._picking_changed()

    def _picking_changed(self, *_):
        '''
        Put the three buttons together into one click action
        '''
        if not self.select_button.isChecked():
            action = 'off'
        elif self.add_button_plus.isChecked():
            action = 'add'
        elif self.add_button_minus.isChecked():
            action = 'subtract'
        else:
            action = 'replace'
        picking.set_click_action(action, self.cmd)
        for btn in (self.add_button_plus, self.add_button_minus):
            btn.setEnabled(self.select_button.isChecked())

    def _add_zoom(self):
        '''
        Zoom on everything, or on what is selected
        '''
        cmd = self.cmd
        self.zoom_all_button = self.caption_button(
            'zoom_all_button', 'zoom-all', 'Zoom All',
            'Zoom on all objects', lambda: cmd.zoom(animate=1.0))
        self.zoom_selection_button = self.caption_button(
            'zoom_selection_button', 'zoom-selection', 'Zoom Sel',
            'Zoom on the selection', lambda: cmd.zoom('sele', animate=1.0))

        def update():
            self.zoom_selection_button.setEnabled(
                'sele' in (cmd.get_names('selections') or []))

        self.zoom_timer = QtCore.QTimer(self)
        self.zoom_timer.timeout.connect(update)
        self.zoom_timer.start(400)
        update()

    def _add_presets(self):
        '''
        The representation presets, with the sliders icon and its name
        underneath, as Maestro's buttons look
        '''
        self.presets_button = btn = self.caption_button(
            'presets_button', 'sliders', 'Presets',
            'Apply a representation preset to all objects',
            menu=self._dynamic_menu(
                lambda: pymol_menu.presets(self.cmd, 'all')))
        return btn

    def _add_undo_redo(self):
        '''
        Undo/redo buttons, enabled when there is something to undo/redo.
        '''
        cmd = self.cmd
        self.undo_button = self.caption_button(
            'undo_button', 'undo', 'Undo', 'Undo', lambda: cmd.undo())
        self.redo_button = self.caption_button(
            'redo_button', 'redo', 'Redo', 'Redo', lambda: cmd.redo())

        def update():
            self.undo_button.setEnabled(undo.stack.can_undo())
            self.redo_button.setEnabled(undo.stack.can_redo())

        self.undo_timer = QtCore.QTimer(self)
        self.undo_timer.timeout.connect(update)
        self.undo_timer.start(250)
        update()

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

    def _more_menu(self):
        cmd = self.cmd
        window = self.window
        menu = QtWidgets.QMenu(self)
        # the view actions the toolbar no longer has its own buttons for
        menu.addAction('Orient', lambda: cmd.orient(animate=1.0))
        menu.addAction('Center', lambda: cmd.center(animate=1.0))
        menu.addAction('Reset View', cmd.reset)
        rock = menu.addAction('Rock')
        rock.setCheckable(True)
        rock.triggered.connect(lambda on: cmd.rock(int(on)))
        self.rock_action = rock

        def show_rock():
            rock.blockSignals(True)
            rock.setChecked(bool(cmd.get_setting_int('rock')))
            rock.blockSignals(False)

        menu.aboutToShow.connect(show_rock)
        menu.addSeparator()
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
