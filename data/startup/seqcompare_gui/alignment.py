'''
The alignment of the Sequence Comparison dialog, painted: the chains'
names stay at the left while their residues scroll. Two layouts:
wrapped (as many columns per line as fit, block under block) and
extended (each group of chains on one line, moved along with the slider
of AlignmentScroll below it).
'''

import math
import os

from pymol.Qt import QtCore, QtGui, QtWidgets

from pymolx import compare, seqcompare
from pymolx.gui import theme

Qt = QtCore.Qt

LIGHTEN = 0.45      # residue backgrounds: score colors toward white
PAD = 6             # pixels around the alignment
RULER_STEP = 10     # a residue number every so many columns...
RULER_STEPS = (1, 2, 5)  # ... or closer, with letters larger than normal
NUMBER_GAP = 0.6    # characters between two residue numbers, at least
MIN_NUMBER_SIZE = 7  # points: numbers above every residue get smaller
MIN_BLOCK = 10      # columns per line of the wrapped layout, at least
# letter sizes the + and - buttons step through (points, as for the
# interface font: the range of the sequence viewer)
FONT_SIZES = (6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 21, 24, 28, 32, 36)
WHEEL_STEP = 120    # angleDelta of one wheel notch
NAME_SHARE = 0.4    # of the view's width, at most, for the names
MIN_NAME = 8        # characters of a structure's name, at least
SCROLL_COLUMNS = 3  # columns per step of the slider (arrow keys, wheel)

WRAPPED, EXTENDED = 'wrapped', 'extended'
LAYOUTS = [(WRAPPED, 'Wrapped'), (EXTENDED, 'Extended')]
MESSAGE = 'Compare structures to see their alignment'


def background(score):
    '''
    "#rrggbb" for a residue's score (None: no partner)
    '''
    if score is None:
        return theme.COLORS['border']
    r, g, b = compare.ramp_color(score / 100.0)
    return '#%02x%02x%02x' % tuple(
        round(255 * (c + (1 - c) * LIGHTEN)) for c in (r, g, b))


def label(group, entry):
    return entry.name + (' (ref)' if entry is group.reference else '')


def shorten(names, limit):
    '''
    Names of at most `limit` characters. One that is longer first loses
    what all the names start with, then what all of them end with, so
    that what tells them apart stays; "…" stands for what is left out.

    :return: {name: shown}
    '''
    names = list(names)
    limit = max(3, limit)
    head = tail = 0
    if len(set(names)) > 1:
        head = len(os.path.commonprefix(names))
        tail = len(os.path.commonprefix([n[::-1] for n in names]))
    out = {}
    for name in names:
        shown = name
        if len(shown) > limit and head:
            drop = min(head, len(shown) - limit + 1)
            shown = '…' + shown[drop:]
        if len(shown) > limit and tail:
            drop = min(tail, len(shown) - limit + 1)
            shown = shown[:len(shown) - drop] + '…'
        if len(shown) > limit:
            # too long by itself: its two ends
            end = (limit - 1) // 2
            shown = shown[:limit - 1 - end] + '…' + shown[len(shown) - end:]
        out[name] = shown
    return out


def lines(result, block=None):
    '''
    The alignment's lines from top to bottom: (group, entry, start, stop)
    for a chain's residues in the columns start..stop, with entry None
    for the residue numbers above them, and None for the empty line
    between two blocks.

    :param block: columns per line (default: a group's columns on one)
    '''
    out = []
    for group in result.groups:
        n = len(group.columns)
        size = block or n or 1
        for start in range(0, n, size):
            stop = min(n, start + size)
            if out:
                out.append(None)
            out.append((group, None, start, stop))
            out += [(group, e, start, stop) for e in group.entries]
    return out


def ruler_plan(font_size, digits, normal=None):
    '''
    How close the residue numbers are for letters of a size, and how
    large: every tenth column up to the normal letter size; with larger
    letters every fifth or second column once the numbers fit at the
    normal size, and every column once they fit at a smaller one.

    :param font_size: the letters' size (points)
    :param digits: characters of the longest residue number
    :param normal: the normal letter size (default: theme.FONT_SIZE)
    :return: (columns from number to number, the numbers' size)
    '''
    normal = theme.FONT_SIZE if normal is None else normal
    if font_size <= normal:
        return RULER_STEP, font_size
    # (a letter is as wide as its column, and widths go with the size:
    # in points, a number takes `need` times its own size)
    need = digits + NUMBER_GAP
    size = min(normal, int(font_size / need))
    if size >= MIN_NUMBER_SIZE:
        return 1, size
    for step in RULER_STEPS[1:]:
        if need * normal <= step * font_size:
            return step, normal
    return RULER_STEP, normal


def ruler(group, start, stop, step=RULER_STEP):
    '''
    :return: [(column, the reference's residue number there)] every
        `step` columns from start (a number is written out in full, also
        where it reaches past stop)
    '''
    ref = group.reference
    marks = []
    for k in range(start, stop, step):
        pos = group.columns[k].get(ref.name)
        if pos is not None:
            marks.append((k, ref.residues[pos][0]))
    return marks


def cell(group, entry, column):
    '''
    :return: (letter, position in the entry's residues), or ("-", None)
        where the chain has no residue; lower case without coordinates
    '''
    pos = group.columns[column].get(entry.name)
    if pos is None:
        return '-', None
    resi, resn, letter, index = entry.residues[pos]
    return (letter if index is not None else letter.lower()), pos


def text(result, block=None):
    '''
    The alignment as plain text, laid out like the view
    '''
    width = max(len(label(g, e)) for g in result.groups for e in g.entries)
    out = []
    for line in lines(result, block):
        if line is None:
            out.append('')
            continue
        group, entry, start, stop = line
        if entry is None:
            chars = []
            for k, number in ruler(group, start, stop):
                chars += ' ' * (k - start - len(chars))
                chars[k - start:] = number
            out.append((' ' * (width + 1) + ''.join(chars)).rstrip())
        else:
            out.append('%s %s' % (label(group, entry).ljust(width), ''.join(
                cell(group, entry, k)[0] for k in range(start, stop))))
    return '\n'.join(out) + '\n'


def fit_block(columns):
    '''
    Columns per line of the wrapped layout, where so many fit: whole
    tens, so that the residue numbers stay at every tenth column
    '''
    return max(MIN_BLOCK, int(columns) // RULER_STEP * RULER_STEP)


class AlignmentView(QtWidgets.QAbstractScrollArea):
    '''
    The colored alignment. Up and down with the scroll bar at the right;
    along the chains with `offset` (pixels), which AlignmentScroll's
    slider sets.
    '''

    view_changed = QtCore.Signal()   # scrolled, letter size, layout

    def __init__(self, parent=None, layout=WRAPPED, font_size=None):
        super().__init__(parent)
        self.result = None
        self.layout_name = layout if layout in dict(LAYOUTS) else WRAPPED
        self.font_size = theme.FONT_SIZE if font_size is None else font_size
        self.block = None    # columns per line (None: all on one)
        self.lines = []
        self.labels = {}     # entry name: text shown at the left
        self.name_width = 0.0
        self.digits = 1      # characters of the longest residue number
        # the kinds of residues in color (seqcompare.KINDS); the others
        # are greyed out
        self.shown = set(seqcompare.KINDS)
        self._wheel = 0      # touchpads send fractions of a notch
        self._scores = {}
        self._kinds = {}     # entry name: {position: kind}
        self.setObjectName('seqcompare_alignment')
        # (the slider below takes this bar's place: see AlignmentScroll)
        self.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMinimumHeight(120)
        self._measure()

    # what is shown

    @property
    def extended(self):
        return self.layout_name == EXTENDED

    def set_result(self, result):
        self.result = result
        self._scores = {}
        self._kinds = {}
        self.digits = 1
        if result is not None:
            for name, scores in result.scores.items():
                self._scores[name] = {
                    pos: QtGui.QColor(background(score))
                    for pos, score in scores.items()}
                self._kinds[name] = {pos: result.kind(name, pos)
                                     for pos in scores}
            self.digits = max([1] + [
                len(r[0]) for g in result.groups
                for r in g.reference.residues])
        self._measure()
        self.relayout()
        self.horizontalScrollBar().setValue(0)
        self.verticalScrollBar().setValue(0)
        self.viewport().update()
        self.view_changed.emit()

    def set_layout(self, name):
        '''
        Wrapped or extended, at the columns that were in view
        '''
        if name == self.layout_name or name not in dict(LAYOUTS):
            return
        where = self.place()
        self.layout_name = name
        self.relayout()
        self.go_to(where)
        self.viewport().update()
        self.view_changed.emit()

    def set_shown(self, kinds):
        '''
        The kinds of residues to show in color (seqcompare.KINDS); the
        others keep their letters, greyed out
        '''
        kinds = set(kinds) & set(seqcompare.KINDS)
        if kinds != self.shown:
            self.shown = kinds
            self.viewport().update()

    def is_shown(self, entry, position):
        '''
        Whether a residue (its index in the entry's residues) is of a
        kind that is shown
        '''
        return self._kinds.get(entry.name, {}).get(
            position, 'none') in self.shown

    def text(self):
        '''
        The alignment as plain text, in the lines it is shown in
        '''
        return '' if self.result is None else text(self.result, self.block)

    # letter size

    def text_font(self):
        return theme.console_font(size=self.font_size)

    def number_font(self):
        return theme.console_font(size=self.number_size)

    def _measure(self):
        metrics = QtGui.QFontMetricsF(self.text_font())
        self.cell_width = max(1.0, metrics.horizontalAdvance('M'))
        self.line_height = max(1, int(math.ceil(metrics.height())))
        self.ascent = metrics.ascent()
        # the residue numbers: how close, and in what size
        self.ruler_step, self.number_size = ruler_plan(
            self.font_size, self.digits)
        metrics = QtGui.QFontMetricsF(self.number_font())
        self.number_width = metrics.horizontalAdvance('0')
        self.number_descent = metrics.descent()

    def can_zoom_in(self):
        return self.font_size < FONT_SIZES[-1]

    def can_zoom_out(self):
        return self.font_size > FONT_SIZES[0]

    def zoom_in(self):
        self.set_font_size(min(
            (s for s in FONT_SIZES if s > self.font_size),
            default=FONT_SIZES[-1]))

    def zoom_out(self):
        self.set_font_size(max(
            (s for s in FONT_SIZES if s < self.font_size),
            default=FONT_SIZES[0]))

    def set_font_size(self, size, x=None):
        '''
        Larger or smaller letters, keeping the column at the pixel x
        (default: the middle of the view) in its place
        '''
        size = min(FONT_SIZES[-1], max(FONT_SIZES[0], size))
        if size == self.font_size:
            return
        where = self.place()
        if x is None:
            x = (self.name_width + self.viewport().width()) / 2.0
        share = (x - self.name_width) / max(
            1.0, self.viewport().width() - self.name_width)
        column = (x - self.name_width + self.offset) / self.cell_width
        self.font_size = size
        self._measure()
        self.relayout()
        self.go_to(where)
        if self.extended:
            x = self.name_width + share * (
                self.viewport().width() - self.name_width)
            self.set_offset(column * self.cell_width - (x - self.name_width))
        self.viewport().update()
        self.view_changed.emit()

    # geometry

    @property
    def offset(self):
        '''
        Pixels the residues are moved to the left
        '''
        return self.horizontalScrollBar().value()

    def set_offset(self, offset):
        bar = self.horizontalScrollBar()
        value = min(bar.maximum(), max(0, int(round(offset))))
        if value != bar.value():
            bar.setValue(value)   # (scrollContentsBy tells the slider)
        elif value != offset:
            self.view_changed.emit()

    def scroll_range(self):
        return self.horizontalScrollBar().maximum()

    def page_width(self):
        '''
        Pixels of the view the residues have (right of the names)
        '''
        return max(1, int(self.viewport().width() - self.name_width))

    def scroll_step(self):
        return max(1, int(round(self.cell_width * SCROLL_COLUMNS)))

    def relayout(self):
        '''
        The names' width, the lines and what there is to scroll, for the
        result, the layout, the letter size and the view's size
        '''
        width = self.viewport().width()
        height = self.viewport().height()
        hbar, vbar = self.horizontalScrollBar(), self.verticalScrollBar()
        self.labels = {}
        if self.result is None or not self.result.groups:
            self.lines, self.block, self.name_width = [], None, 0.0
            hbar.setRange(0, 0)
            vbar.setRange(0, 0)
            return
        # the names: a share of the view at most, the structures' names
        # shortened to leave room for "/chain (ref)"
        entries = [(g, e) for g in self.result.groups for e in g.entries]
        room = int((NAME_SHARE * width - PAD) / self.cell_width) - max(
            len(label(g, e)) - len(e.obj) for g, e in entries)
        short = shorten([e.obj for g, e in entries], max(MIN_NAME, room))
        for g, e in entries:
            self.labels[e.name] = short[e.obj] + label(g, e)[len(e.obj):]
        self.name_width = PAD + self.cell_width * (
            1 + max(len(shown) for shown in self.labels.values()))
        columns = 0
        for group in self.result.groups:
            # (with the last residue number, which may reach past the end)
            n = len(group.columns)
            columns = max([columns, n] + [
                k + len(number) for k, number in
                ruler(group, 0, n, self.ruler_step)[-1:]])
        if self.extended:
            self.block = None
        else:
            self.block = fit_block(
                (width - self.name_width - PAD) / self.cell_width)
            columns = min(columns, self.block)
        self.lines = lines(self.result, self.block)
        content = self.name_width + columns * self.cell_width + PAD
        hbar.setRange(0, max(0, int(math.ceil(content - width))))
        hbar.setPageStep(self.page_width())
        hbar.setSingleStep(self.scroll_step())
        content = 2 * PAD + len(self.lines) * self.line_height
        vbar.setRange(0, max(0, content - height))
        vbar.setPageStep(max(1, height))
        vbar.setSingleStep(self.line_height)

    def _line_top(self, index):
        return PAD + index * self.line_height - \
            self.verticalScrollBar().value()

    def _column_x(self, start, column):
        return self.name_width + (column - start) * self.cell_width - \
            self.offset

    def place(self):
        '''
        What is at the view's top left corner: (group, column, lines
        below the residue numbers of its block), or None
        '''
        if not self.lines:
            return None
        index = int((self.verticalScrollBar().value() - PAD) //
                    self.line_height)
        index = min(len(self.lines) - 1, max(0, index))
        if self.lines[index] is None:
            index += 1
        group, entry, start, stop = self.lines[index]
        column = min(stop - 1, start + int(self.offset // self.cell_width))
        below = 0 if entry is None else 1 + group.entries.index(entry)
        return group, column, below

    def go_to(self, place):
        '''
        Bring a place() back to the view's top left corner
        '''
        if place is None:
            return
        group, column, below = place
        for index, line in enumerate(self.lines):
            if line is not None and line[0] is group and line[1] is None \
                    and line[2] <= column < line[3]:
                top = index + below
                # (the first line keeps the margin above it)
                self.verticalScrollBar().setValue(
                    top * self.line_height + (PAD if top else 0))
                # (a whole pixel that lies in this column, not before it)
                self.set_offset(math.ceil(
                    (column - line[2]) * self.cell_width))
                return

    def hit(self, pos):
        '''
        :return: (group, entry, column or None) under a point of the
            view, column None on the name; None off the chains
        '''
        index = int((pos.y() + self.verticalScrollBar().value() - PAD) //
                    self.line_height)
        if not 0 <= index < len(self.lines) or self.lines[index] is None:
            return None
        group, entry, start, stop = self.lines[index]
        if entry is None:
            return None
        if pos.x() < self.name_width:
            return group, entry, None
        column = start + int((pos.x() - self.name_width + self.offset) //
                             self.cell_width)
        if not start <= column < stop:
            return None
        return group, entry, column

    def describe(self, pos):
        '''
        Tooltip for a point of the view: the residue there, or a name
        that is too long to show in full
        '''
        found = self.hit(pos)
        if found is None:
            return ''
        group, entry, column = found
        if column is None:
            full = label(group, entry)
            return full if self.labels.get(entry.name) != full else ''
        index = cell(group, entry, column)[1]
        if index is None:
            return ''
        resi, resn, letter, atom = entry.residues[index]
        described = '%s %s %s' % (entry.name, resn, resi)
        if atom is None:
            described += '\nnot in the model (no coordinates)'
        return described

    def visible_columns(self, start, stop):
        '''
        The columns of a line (start..stop) in view: (first, end)
        '''
        first = start + int(self.offset // self.cell_width)
        end = start + int((self.offset + self.page_width()) //
                          self.cell_width) + 1
        return max(start, first), min(stop, end)

    # painting

    def draw(self, painter):
        colors = theme.COLORS
        width = self.viewport().width()
        height = self.viewport().height()
        painter.fillRect(QtCore.QRect(0, 0, width, height),
                         QtGui.QColor(colors['base']))
        painter.setFont(self.text_font())
        if not self.lines:
            painter.setPen(QtGui.QColor(colors['text_dim']))
            painter.drawText(QtCore.QRect(0, 0, width, height),
                             int(Qt.AlignmentFlag.AlignCenter), MESSAGE)
            return
        text_color = QtGui.QColor(colors['text'])
        dim = QtGui.QColor(colors['text_dim'])
        dark = QtGui.QColor(colors['base'])
        cw = self.cell_width
        names = int(math.floor(self.name_width))
        first = max(0, int((self.verticalScrollBar().value() - PAD) //
                           self.line_height))
        last = min(len(self.lines), first + height // self.line_height + 2)
        for index in range(first, last):
            line = self.lines[index]
            if line is None:
                continue
            group, entry, start, stop = line
            top = self._line_top(index)
            base = top + self.ascent
            painter.setClipRect(QtCore.QRect(names, 0, width - names,
                                             height))
            k0, k1 = self.visible_columns(start, stop)
            if entry is None:
                self._draw_numbers(painter, group, start, stop, top)
                continue
            scores = self._scores.get(entry.name, {})
            kinds = self._kinds.get(entry.name, {})
            none = QtGui.QColor(colors['border'])
            cells = [cell(group, entry, k) for k in range(k0, k1)]
            # (residues of a kind that is not shown: no color)
            plain = {pos for _, pos in cells if pos is not None and
                     kinds.get(pos, 'none') not in self.shown}
            # backgrounds: one rectangle per run of a color
            run_start, run_color = k0, None
            for k, (letter, pos) in enumerate(cells + [('', None)], k0):
                color = None if pos is None or pos in plain else \
                    scores.get(pos, none)
                if k == k1 or color != run_color:
                    if run_color is not None:
                        left = math.floor(self._column_x(start, run_start))
                        right = math.floor(self._column_x(start, k))
                        painter.fillRect(QtCore.QRect(
                            left, top, right - left, self.line_height),
                            run_color)
                    run_start, run_color = k, color
            for k, (letter, pos) in enumerate(cells, k0):
                painter.setPen(text_color if pos is None else
                               dim if pos in plain else dark)
                painter.drawText(QtCore.QPointF(
                    self._column_x(start, k), base), letter)
            painter.setClipping(False)
            painter.setPen(text_color)
            painter.drawText(QtCore.QPointF(PAD, base),
                             self.labels.get(entry.name, entry.name))
        painter.setClipping(False)
        if self.offset:
            # residues are hidden behind the names
            painter.setPen(QtGui.QColor(colors['border']))
            x = names - int(cw / 2)
            painter.drawLine(x, 0, x, height)

    def _draw_numbers(self, painter, group, start, stop, top):
        '''
        The reference's residue numbers above a block. In the letters'
        size they start at their column, as in the text; smaller ones
        (larger letters) are centered on their residue, just above it.
        '''
        small = self.number_size < self.font_size
        cw = self.cell_width
        base = top + self.ascent
        if small:
            base = top + self.line_height - self.number_descent - 1
            painter.setFont(self.number_font())
            # (a number is wider than its residue: let the first one
            # reach into the space after the names)
            left = int(math.floor(self.name_width - cw / 2))
            painter.setClipRect(QtCore.QRect(
                left, 0, self.viewport().width() - left,
                self.viewport().height()))
        painter.setPen(QtGui.QColor(theme.COLORS['text_dim']))
        k0, k1 = self.visible_columns(start, stop)
        for k, number in ruler(group, start, stop, self.ruler_step):
            if k0 - len(number) <= k <= k1:
                x = self._column_x(start, k)
                if small:
                    x += (cw - len(number) * self.number_width) / 2.0
                painter.drawText(QtCore.QPointF(x, base), number)
        if small:
            painter.setFont(self.text_font())

    def paintEvent(self, event):
        painter = QtGui.QPainter(self.viewport())
        self.draw(painter)

    # events

    def resizeEvent(self, event):
        where = self.place()
        block = self.block
        self.relayout()
        if self.block != block:
            self.go_to(where)
        self.view_changed.emit()
        super().resizeEvent(event)

    def scrollContentsBy(self, dx, dy):
        self.viewport().update()
        self.view_changed.emit()

    def wheelEvent(self, event):
        delta = event.angleDelta()
        if event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            # Ctrl + wheel: letter size, at the pointer
            self._wheel += delta.y() or delta.x()
            steps = int(self._wheel / WHEEL_STEP)
            if steps:
                self._wheel -= steps * WHEEL_STEP
                point = event.position() if hasattr(event, 'position') \
                    else event.pos()
                size = self.font_size
                for _ in range(abs(steps)):
                    larger = [s for s in FONT_SIZES if s > size]
                    smaller = [s for s in FONT_SIZES if s < size]
                    size = (min(larger, default=size) if steps > 0 else
                            max(smaller, default=size))
                self.set_font_size(size, max(point.x(), self.name_width))
            event.accept()
            return
        vbar = self.verticalScrollBar()
        shift = event.modifiers() & Qt.KeyboardModifier.ShiftModifier
        if delta.x() or vbar.maximum() == 0 or shift:
            # along the chains: nothing to scroll up and down, or Shift
            self.set_offset(self.offset - (delta.x() or delta.y()) *
                            self.scroll_step() / WHEEL_STEP)
            event.accept()
        else:
            super().wheelEvent(event)

    def viewportEvent(self, event):
        if event.type() == QtCore.QEvent.Type.ToolTip:
            described = self.describe(event.pos())
            if described:
                QtWidgets.QToolTip.showText(event.globalPos(), described,
                                            self.viewport())
            else:
                QtWidgets.QToolTip.hideText()
            return True
        return super().viewportEvent(event)

    def keyPressEvent(self, event):
        if event.matches(QtGui.QKeySequence.StandardKey.Copy):
            self.copy()
            event.accept()
            return
        super().keyPressEvent(event)

    def copy(self):
        QtWidgets.QApplication.clipboard().setText(self.text())


class AlignmentScroll(QtWidgets.QWidget):
    '''
    The row under the alignment: the layout, the slider that moves along
    chains too long to show at once, and the buttons for smaller and
    larger letters
    '''

    layout_chosen = QtCore.Signal(str)

    def __init__(self, view, parent=None):
        super().__init__(parent)
        self.view = view
        row = QtWidgets.QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        self.layout_box = QtWidgets.QComboBox(self)
        for name, title in LAYOUTS:
            self.layout_box.addItem(title, name)
        self.layout_box.setToolTip(
            'Wrapped: as many residues per line as fit, block under '
            'block.\nExtended: each chain on one line; move along it with '
            'the slider.')
        self.layout_box.setCurrentIndex(
            self.layout_box.findData(view.layout_name))
        self.layout_box.currentIndexChanged.connect(self._choose)
        row.addWidget(self.layout_box)
        self.slider = QtWidgets.QScrollBar(Qt.Orientation.Horizontal, self)
        self.slider.setObjectName('alignment_slider')
        self.slider.setToolTip('Move along the chains (also Shift + mouse '
                               'wheel over the alignment)')
        self.slider.valueChanged.connect(view.set_offset)
        row.addWidget(self.slider, 1)
        self.smaller = QtWidgets.QPushButton('−', self)
        self.smaller.setToolTip('Smaller letters: more residues in view '
                                '(also Ctrl + mouse wheel)')
        self.smaller.clicked.connect(view.zoom_out)
        self.larger = QtWidgets.QPushButton('+', self)
        self.larger.setToolTip('Larger letters: fewer residues in view '
                               '(also Ctrl + mouse wheel)')
        self.larger.clicked.connect(view.zoom_in)
        for button in (self.smaller, self.larger):
            button.setFixedWidth(32)
            button.setAutoDefault(False)
            row.addWidget(button)
        view.view_changed.connect(self.follow)
        self.follow()

    def _choose(self, index):
        name = self.layout_box.itemData(index)
        self.view.set_layout(name)
        self.layout_chosen.emit(name)

    def follow(self):
        '''
        Take the slider's place and length, and what the buttons can
        still do, from the view
        '''
        view = self.view
        slider = self.slider
        slider.blockSignals(True)
        # (Qt draws nothing of a slider without a range: with one of a
        # single pixel, an alignment that fits shows as a handle from
        # end to end, which set_offset keeps in its place)
        slider.setRange(0, max(1, view.scroll_range()))
        slider.setPageStep(view.page_width())
        slider.setSingleStep(view.scroll_step())
        slider.setValue(view.offset)
        slider.blockSignals(False)
        self.larger.setEnabled(view.can_zoom_in())
        self.smaller.setEnabled(view.can_zoom_out())
        index = self.layout_box.findData(view.layout_name)
        if index != self.layout_box.currentIndex():
            self.layout_box.blockSignals(True)
            self.layout_box.setCurrentIndex(index)
            self.layout_box.blockSignals(False)
