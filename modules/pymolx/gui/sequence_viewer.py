'''
Sequence viewer (parity item L-07): Qt replacement for the OpenGL
sequence bar, docked above the viewer. Rows and selection logic come
from pymolx.sequence; this module only draws them and handles the mouse.

Mouse, like the OpenGL viewer: click a residue to add it to the active
selection (or remove it if it is selected), drag for ranges, shift-click
to extend the last range, click a chain label (or "/A/" marker) for the
whole chain, middle-click to center, right-click for the residue or
selection menu, double-click empty space to clear the selection.

Rows are chains, or objects with their chains side by side like the
OpenGL viewer (right-click the label column, or Display > Sequence
Mode). The label column's right edge can be dragged; double-click it to
fit the labels again. Both choices are kept in the plugin preferences.

The SEQ toggle and Display > Sequence show this viewer instead of the
OpenGL one; "set seq_view, 1" does the same (the setting is put back to
0, so the OpenGL bar never shows as well).
'''

import contextlib

from pymol.Qt import QtCore, QtGui, QtWidgets
from pymol import menu as pymol_menu

from .. import sequence
from . import menus, theme

Qt = QtCore.Qt

REFRESH_MS = 150
MAX_VISIBLE_ROWS = 6
LABEL_MAX_CHARS = 14
PAD = 4        # pixels around the cell area
ROW_GAP = 3    # pixels between rows
MENU_SELE = '_seqview_menu'
DIVIDER_GRAB = 3  # pixels either side of the label column edge

PREF_BY_OBJECT = 'pymolx_seq_by_object'
PREF_LABEL_WIDTH = 'pymolx_seq_label_width'

LAYOUT_LABELS = [(False, 'One Row per Chain'),
                 (True, 'One Row per Object (Chains Side by Side)')]


class MemoryPrefs(dict):
    '''
    Preferences kept for this run only (tests, no plugin system)
    '''

    def set(self, key, value):
        self[key] = value


class PluginPrefs:
    '''
    Preferences saved with the plugin preferences (~/.pymolpluginsrc.py)
    '''

    def get(self, key, default=None):
        from pymol import plugins
        return plugins.pref_get(key, default)

    def set(self, key, value):
        from pymol import plugins
        plugins.pref_set(key, value)


def _event_pos(event):
    return (event.position().toPoint() if hasattr(event, 'position') else
            event.pos())


def _global_pos(event):
    return (event.globalPosition().toPoint()
            if hasattr(event, 'globalPosition') else event.globalPos())


class SequenceViewer(QtWidgets.QAbstractScrollArea):
    '''
    One row per chain (or per object): label at the left, residue numbers
    and chain markers above the residue cells. Cells scroll horizontally,
    all rows together.
    '''

    shownChanged = QtCore.Signal(bool)
    layoutChanged = QtCore.Signal(bool)  # by_object

    def __init__(self, cmd, parent=None, prefs=None):
        super().__init__(parent)
        self.cmd = cmd
        self.prefs = prefs if prefs is not None else MemoryPrefs()
        self.by_object = bool(self.prefs.get(PREF_BY_OBJECT, False))
        # pixels, or None to fit the labels
        self.label_width = self.prefs.get(PREF_LABEL_WIDTH, None)
        self._resizing_label = False
        self._divider_hover = False
        self.setObjectName('sequence_viewer')
        self.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.viewport().setMouseTracking(True)

        self.rows = []
        self.selected = set()      # residue keys of the active selection
        self.shown = False         # wanted by the user (SEQ toggle)
        self._counts = None
        self._objects = None
        self._colors = {}
        self._drag = None          # (row, anchor index, add)
        self._drag_last = None     # residue index the drag is at
        self._anchor = None        # (row key, index, add) for shift-click

        font = theme.console_font()
        self.setFont(font)
        self._set_metrics()

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.refresh)
        self.timer.start(REFRESH_MS)

    # geometry

    def _set_metrics(self):
        fm = QtGui.QFontMetrics(self.font())
        self.char_w = max(1, fm.horizontalAdvance('W'))
        self.line_h = fm.height()
        self.row_h = 2 * self.line_h + ROW_GAP
        self.ascent = fm.ascent()

    def _label_width(self):
        if self.label_width:
            return self._clamp_label_width(self.label_width)
        # one character of slack: char_w is rounded down
        chars = max([len(row.label) for row in self.rows] + [4])
        return (min(chars, LABEL_MAX_CHARS) + 1) * self.char_w + 2 * PAD

    def _clamp_label_width(self, width):
        # at least a few characters, and room left for the cells
        lo = 3 * self.char_w + 2 * PAD
        hi = max(lo, self.viewport().width() - 10 * self.char_w)
        return int(min(max(width, lo), hi))

    def set_label_width(self, width, save=True):
        '''
        :param width: pixels, or None to fit the labels
        '''
        self.label_width = (None if width is None else
                            self._clamp_label_width(width))
        if save:
            self.prefs.set(PREF_LABEL_WIDTH, self.label_width)
        self._update_scrollbars()
        self.viewport().update()

    def _on_divider(self, x):
        return abs(x - self._label_width()) <= DIVIDER_GRAB

    def set_by_object(self, by_object):
        '''
        One row per object, chains side by side (True), or one row per
        chain (False)
        '''
        by_object = bool(by_object)
        if by_object != self.by_object:
            self.by_object = by_object
            self.prefs.set(PREF_BY_OBJECT, by_object)
            self._anchor = None
            if self.shown:
                self.refresh(force=True)
        self.layoutChanged.emit(by_object)

    def layout_menu(self, parent=None):
        '''
        Menu with the row layouts and "Fit Label Width"
        '''
        menu = QtWidgets.QMenu(parent or self)
        group = QtWidgets.QActionGroup(menu)
        for by_object, label in LAYOUT_LABELS:
            action = menu.addAction(label)
            action.setCheckable(True)
            action.setChecked(by_object == self.by_object)
            action.triggered.connect(
                lambda _=False, b=by_object: self.set_by_object(b))
            group.addAction(action)
        menu.addSeparator()
        action = menu.addAction('Fit Label Width')
        action.setEnabled(bool(self.label_width))
        action.triggered.connect(lambda: self.set_label_width(None))
        return menu

    def content_width(self):
        ncols = max([row.ncols for row in self.rows] + [0])
        return ncols * self.char_w + 2 * PAD

    def _cells_width(self):
        return max(1, self.viewport().width() - self._label_width())

    def preferred_height(self):
        '''
        Height for all rows (at most MAX_VISIBLE_ROWS) plus the
        horizontal scroll bar if the cells don't fit.
        '''
        nrows = max(1, min(len(self.rows), MAX_VISIBLE_ROWS))
        height = nrows * self.row_h + 2 * PAD - ROW_GAP
        frame = self.width() - self.viewport().width()
        if len(self.rows) > MAX_VISIBLE_ROWS:
            frame = self.verticalScrollBar().sizeHint().width()
        label_w = self._label_width()
        if self.content_width() > self.width() - frame - label_w:
            height += self.horizontalScrollBar().sizeHint().height()
        return height

    def sizeHint(self):
        return QtCore.QSize(400, self.preferred_height())

    def _update_scrollbars(self):
        hbar = self.horizontalScrollBar()
        cells_w = self._cells_width()
        hbar.setRange(0, max(0, self.content_width() - cells_w))
        hbar.setPageStep(cells_w)
        hbar.setSingleStep(self.char_w)
        vbar = self.verticalScrollBar()
        total = len(self.rows) * self.row_h + 2 * PAD - ROW_GAP
        vbar.setRange(0, max(0, total - self.viewport().height()))
        vbar.setPageStep(self.viewport().height())
        vbar.setSingleStep(self.row_h)
        height = self.preferred_height()
        if self.height() != height:
            self.setFixedHeight(height)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._update_scrollbars()

    def _row_top(self, i):
        return PAD + i * self.row_h - self.verticalScrollBar().value()

    def _col_x(self, col):
        return (self._label_width() + PAD + col * self.char_w -
                self.horizontalScrollBar().value())

    def hit(self, pos):
        '''
        :return: (row, residue or None, column) under a viewport position;
            column is None on the label, row is None outside the rows
        '''
        y = pos.y() + self.verticalScrollBar().value() - PAD
        i = y // self.row_h if y >= 0 else -1
        if not 0 <= i < len(self.rows):
            return None, None, None
        row = self.rows[i]
        if pos.x() < self._label_width():
            return row, None, None
        col = ((pos.x() - self._label_width() - PAD +
                self.horizontalScrollBar().value()) // self.char_w)
        return row, row.residue_at(col), col

    def residue_rect(self, row, res):
        '''
        Viewport rectangle of a residue cell
        '''
        i = self.rows.index(row)
        top = self._row_top(i) + self.line_h
        return QtCore.QRect(self._col_x(res.col), top,
                            res.width * self.char_w, self.line_h)

    def scroll_to(self, row, res):
        '''
        Scroll so a residue is visible
        '''
        rect = self.residue_rect(row, res)
        hbar = self.horizontalScrollBar()
        left = self._label_width()
        if rect.left() < left:
            hbar.setValue(hbar.value() - (left - rect.left()) - PAD)
        elif rect.right() > self.viewport().width():
            hbar.setValue(hbar.value() + rect.right() -
                          self.viewport().width() + PAD)

    # data

    def set_shown(self, shown):
        '''
        Show or hide (SEQ toggle). Hidden while there are no polymers.
        '''
        shown = bool(shown)
        if shown != self.shown:
            self.shown = shown
            self.refresh(force=True)
            self.shownChanged.emit(shown)

    def has_polymer(self):
        return any(row.has_polymer for row in self.rows)

    def refresh(self, force=False):
        '''
        Poll the core's change counters; rebuild the rows when residues,
        colors or objects changed, update the highlight when selections
        changed. Skipped while hidden by the user.
        '''
        if not self.shown:
            self._set_visible(False)
            return
        cmd = self.cmd
        counts = sequence.get_change_counts(_self=cmd)
        if counts == self._counts and not force:
            return
        changed = force or self._counts is None or counts[0] != self._counts[0]
        if not changed:
            # selections or enabled objects
            changed = sequence._shown_objects(cmd) != self._objects
        if changed:
            self._objects = sequence._shown_objects(cmd)
            self.rows = sequence.get_rows(by_object=self.by_object,
                                          _self=cmd)
            self._colors.clear()
            self._drag = None
        self.selected = sequence.get_selected(_self=cmd)
        # after the queries: their temporary selections count as changes
        self._counts = sequence.get_change_counts(_self=cmd)
        self._update_scrollbars()
        self._set_visible(self.has_polymer())
        self.viewport().update()

    def _set_visible(self, visible):
        # the dock (if any) shows or hides with the viewer
        target = self.parentWidget()
        if not isinstance(target, QtWidgets.QDockWidget):
            target = self
        if target.isHidden() == visible:
            target.setVisible(visible)

    def qcolor(self, index):
        color = self._colors.get(index)
        if color is None:
            rgb = self.cmd.get_color_tuple(index) or (1.0, 1.0, 1.0)
            color = self._colors[index] = QtGui.QColor.fromRgbF(*rgb[:3])
        return color

    # drawing

    def paintEvent(self, event):
        painter = QtGui.QPainter(self.viewport())
        painter.setFont(self.font())
        width = self.viewport().width()
        painter.fillRect(self.viewport().rect(),
                         QtGui.QColor(theme.COLORS['base']))

        label_w = self._label_width()
        cw = self.char_w
        first_col = max(0, (self.horizontalScrollBar().value() - PAD) // cw)
        last_col = first_col + self._cells_width() // cw + 2
        dim = QtGui.QColor(theme.COLORS['text_dim'])
        dark = QtGui.QColor(theme.COLORS['base'])
        light = QtGui.QColor(theme.COLORS['text_bright'])

        painter.setClipRect(QtCore.QRect(label_w, 0, width - label_w,
                                         self.viewport().height()))
        for i, row in enumerate(self.rows):
            top = self._row_top(i)
            if top + self.row_h < 0 or top > self.viewport().height():
                continue
            # residue numbers and ticks
            numbered = set()
            painter.setPen(QtGui.QColor(theme.COLORS['text']))
            for col, text, _ in row.markers:
                if col + len(text) >= first_col and col <= last_col:
                    painter.drawText(self._col_x(col), top + self.ascent, text)
                    numbered.update(range(col, col + len(text) + 1))
            painter.setPen(dim)
            for col, text in row.numbers:
                if col + len(text) >= first_col and col <= last_col:
                    painter.drawText(self._col_x(col), top + self.ascent, text)
                    numbered.update(range(col, col + len(text)))
            for col in row.ticks:
                if first_col <= col <= last_col and col not in numbered:
                    x = self._col_x(col) + cw // 2
                    painter.drawLine(x, top + self.line_h - 4,
                                     x, top + self.line_h - 1)
            # residue cells
            base_y = top + self.line_h
            for res in row.residues:
                if res.col + res.width < first_col:
                    continue
                if res.col > last_col:
                    break
                color = self.qcolor(res.color)
                x = self._col_x(res.col)
                if res.key in self.selected:
                    # inverted, like the OpenGL viewer
                    painter.fillRect(x, base_y, res.width * cw, self.line_h,
                                     color)
                    painter.setPen(dark if color.lightnessF() > 0.45
                                   else light)
                else:
                    painter.setPen(color)
                painter.drawText(x, base_y + self.ascent, res.text)

        # label column
        painter.setClipping(False)
        painter.fillRect(0, 0, label_w, self.viewport().height(),
                         QtGui.QColor(theme.COLORS['window']))
        divider = 'accent' if (self._divider_hover or
                               self._resizing_label) else 'border'
        painter.setPen(QtGui.QColor(theme.COLORS[divider]))
        painter.drawLine(label_w - 1, 0, label_w - 1, self.viewport().height())
        painter.setPen(QtGui.QColor(theme.COLORS['text']))
        fm = painter.fontMetrics()
        for i, row in enumerate(self.rows):
            top = self._row_top(i) + self.line_h
            text = fm.elidedText(row.label, Qt.TextElideMode.ElideMiddle,
                                 label_w - 2 * PAD)
            painter.drawText(PAD, top + self.ascent, text)
        painter.end()

    # mouse

    def _no_undo(self):
        '''
        Selecting here is not an undo step, like picking in the viewer
        (and saves a session snapshot per click).
        '''
        pause = getattr(self.cmd, 'UndoPauseCM', None)
        return pause() if pause is not None else contextlib.nullcontext()

    def _select(self, residues, add, base=None):
        with self._no_undo():
            sequence.select_residues(residues, add, base, _self=self.cmd)
        self.selected = sequence.get_selected(_self=self.cmd)
        self.viewport().update()

    def _range(self, row, i, j):
        lo, hi = min(i, j), max(i, j)
        return row.residues[lo:hi + 1]

    def _toggle_all(self, residues):
        add = not all(r.key in self.selected for r in residues)
        self._select(residues, add)
        self._anchor = None

    def mousePressEvent(self, event):
        button = event.button()
        pos = _event_pos(event)
        if button == Qt.MouseButton.LeftButton and self._on_divider(pos.x()):
            self._resizing_label = True
            return
        row, res, col = self.hit(pos)
        if button == Qt.MouseButton.RightButton and res is None:
            self.layout_menu().popup(_global_pos(event))
            return
        if row is None:
            return
        if button == Qt.MouseButton.LeftButton:
            if res is None:
                # row label: the whole row; "/A/" marker: that chain
                residues = (row.residues if col is None else
                            row.marker_at(col))
                if residues:
                    self._toggle_all(residues)
                return
            index = row.residues.index(res)
            shift = event.modifiers() & Qt.KeyboardModifier.ShiftModifier
            if shift and self._anchor and self._anchor[0] == row.key:
                _, anchor, add = self._anchor
                self._select(self._range(row, anchor, index), add)
            else:
                add = res.key not in self.selected
                with self._no_undo():
                    base = sequence.save_base(_self=self.cmd)
                self._select([res], add, base)
                self._drag = (row, index, add)
                self._drag_last = index
                self._anchor = (row.key, index, add)
            if event.modifiers() & Qt.KeyboardModifier.ControlModifier:
                self._center([res])
        elif button == Qt.MouseButton.MiddleButton and res is not None:
            self._center([res])
        elif button == Qt.MouseButton.RightButton and res is not None:
            self._residue_menu(res, _global_pos(event))

    def mouseMoveEvent(self, event):
        pos = _event_pos(event)
        if self._resizing_label:
            self.set_label_width(pos.x(), save=False)
            return
        if self._drag is None:
            hover = self._on_divider(pos.x())
            if hover != self._divider_hover:
                self._divider_hover = hover
                if hover:
                    self.viewport().setCursor(Qt.CursorShape.SplitHCursor)
                else:
                    self.viewport().unsetCursor()
                self.viewport().update()
            return
        row, anchor, add = self._drag
        if row not in self.rows:
            self._drag = None
            return
        # stay in the row of the press, scroll at the edges
        hbar = self.horizontalScrollBar()
        if pos.x() > self.viewport().width() - self.char_w:
            hbar.setValue(hbar.value() + self.char_w)
        elif pos.x() < self._label_width() + self.char_w:
            hbar.setValue(hbar.value() - self.char_w)
        col = ((pos.x() - self._label_width() - PAD + hbar.value()) //
               self.char_w)
        index = row.nearest_index(col)
        if index != self._drag_last:
            self._drag_last = index
            self._select(self._range(row, anchor, index), add,
                         sequence.BASE_SELE)
            self._anchor = (row.key, anchor, add)

    def mouseReleaseEvent(self, event):
        if self._resizing_label:
            self._resizing_label = False
            self.set_label_width(_event_pos(event).x())
            return
        if self._drag is not None:
            self._drag = None
            self._drag_last = None
            with self._no_undo():
                sequence.clear_base(_self=self.cmd)

    def leaveEvent(self, event):
        if self._divider_hover and not self._resizing_label:
            self._divider_hover = False
            self.viewport().unsetCursor()
            self.viewport().update()
        super().leaveEvent(event)

    def mouseDoubleClickEvent(self, event):
        if event.button() != Qt.MouseButton.LeftButton:
            return
        pos = _event_pos(event)
        if self._on_divider(pos.x()):
            self.set_label_width(None)  # fit the labels
            return
        row, res, col = self.hit(pos)
        if res is None and (row is None or (
                col is not None and row.marker_at(col) is None)):
            # empty space: clear the selection, like the OpenGL viewer
            with self._no_undo():
                sequence.deselect_all(_self=self.cmd)
            self.selected = set()
            self.viewport().update()
        else:
            # second click of a fast double click toggles like a click
            self.mousePressEvent(event)

    def wheelEvent(self, event):
        delta = event.angleDelta()
        vbar = self.verticalScrollBar()
        shift = event.modifiers() & Qt.KeyboardModifier.ShiftModifier
        if delta.x() == 0 and (vbar.maximum() == 0 or shift):
            hbar = self.horizontalScrollBar()
            hbar.setValue(hbar.value() - delta.y() * 3 * self.char_w // 120)
            event.accept()
        else:
            super().wheelEvent(event)

    def viewportEvent(self, event):
        if event.type() == QtCore.QEvent.Type.ToolTip:
            pos = event.pos()
            row, res, col = self.hit(pos)
            what = 'object' if row is not None and row.chain is None \
                else 'chain'
            if self._on_divider(pos.x()):
                text = 'Drag to resize, double-click to fit the labels'
            elif res is not None:
                text = res.macro
            elif row is not None and col is None:
                text = '%s\n\nClick to select the %s, right-click for the ' \
                    'row layout' % (row.label, what)
            elif row is not None and row.marker_at(col):
                text = 'Click to select the chain'
            else:
                text = ''
            if text:
                QtWidgets.QToolTip.showText(event.globalPos(), text, self)
            else:
                QtWidgets.QToolTip.hideText()
            return True
        return super().viewportEvent(event)

    def _center(self, residues):
        with self._no_undo():
            name = sequence.residue_selection(residues, _self=self.cmd)
            self.cmd.center(name, animate=-1)
            self.cmd.delete(name)

    def _residue_menu(self, res, pos):
        '''
        Like the OpenGL viewer: the selection's menu on a selected
        residue, else the residue's menu.
        '''
        cmd = self.cmd
        menu = menus.PyMenu('', self)
        name = sequence.active_selection(_self=cmd)
        if name and res.key in self.selected:
            data = pymol_menu.pick_sele(cmd, name, name)
        else:
            # kept while the menu is open; its commands run later (cmd.do)
            with self._no_undo():
                sequence.residue_selection([res], MENU_SELE, _self=cmd)
            data = pymol_menu.seq_option(cmd, MENU_SELE, res.macro)
            menu.aboutToHide.connect(lambda: QtCore.QTimer.singleShot(
                0, lambda: cmd.do('delete ' + MENU_SELE, echo=0, log=0)))
        menus.fill_menu(menu, data, cmd)
        menu.popup(pos)
        return menu


def _add_layout_actions(menu, viewer):
    '''
    Row layout choices at the end of Display > Sequence Mode
    '''
    menu.addSeparator()
    group = QtWidgets.QActionGroup(menu)
    actions = {}
    for by_object, label in LAYOUT_LABELS:
        action = QtWidgets.QAction(label, menu)
        action.setCheckable(True)
        action.setChecked(by_object == viewer.by_object)
        action.triggered.connect(
            lambda _=False, b=by_object: viewer.set_by_object(b))
        group.addAction(action)
        menu.addAction(action)
        actions[by_object] = action
    viewer.layoutChanged.connect(
        lambda by_object: actions[by_object].setChecked(True))
    return actions


def setup(window, prefs=None):
    '''
    Dock the sequence viewer above the viewer (the content panel keeps
    the full height at the right) and connect the SEQ toggle, the
    Display > Sequence menu items and the seq_view setting to it.

    :param prefs: preferences store (default: plugin preferences)
    '''
    cmd = window.cmd
    cmd.set('seq_view', 0, quiet=1)

    viewer = SequenceViewer(cmd, prefs=prefs if prefs is not None
                            else PluginPrefs())
    dock = QtWidgets.QDockWidget('Sequence', window)
    dock.setObjectName('sequence_dock')
    dock.setTitleBarWidget(QtWidgets.QWidget())
    dock.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
    dock.setWidget(viewer)
    window.setCorner(Qt.Corner.TopRightCorner,
                     Qt.DockWidgetArea.RightDockWidgetArea)
    window.addDockWidget(Qt.DockWidgetArea.TopDockWidgetArea, dock)
    dock.hide()

    toggles = []
    content_panel = getattr(window, 'content_panel', None)
    if content_panel is not None:
        button = content_panel.toggles.seq_button
        button.clicked.connect(viewer.set_shown)
        toggles.append(button)
    menu = getattr(window, 'menudict', {}).get('Display')
    if menu is not None:
        for action in menu.actions():
            if action.text() == 'Sequence' and action.isCheckable():
                action.triggered.disconnect()
                action.triggered.connect(viewer.set_shown)
                toggles.append(action)
            elif action.text() == 'Sequence Mode' and action.menu():
                window.sequence_layout_actions = _add_layout_actions(
                    action.menu(), viewer)

    def sync(shown):
        for toggle in toggles:
            toggle.blockSignals(True)
            toggle.setChecked(shown)
            toggle.blockSignals(False)

    viewer.shownChanged.connect(sync)

    # "set seq_view, 1" shows this viewer; the setting goes back to 0 so
    # the OpenGL bar stays hidden (that reset must not hide us again)
    reset = [False]

    def seq_view_changed(value):
        if value:
            reset[0] = True
            cmd.set('seq_view', 0, quiet=1)
            viewer.set_shown(True)
        elif reset[0]:
            reset[0] = False
        else:
            viewer.set_shown(False)
        # callbacks of the toggles (registered earlier) follow the setting
        sync(viewer.shown)

    window.setting_callbacks[cmd.setting._get_index('seq_view')].append(
        seq_view_changed)
    window.sequence_dock = dock
    return viewer
