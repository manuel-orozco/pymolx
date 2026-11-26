'''
The proximity profile of the Interface Analysis dialog's Epitope tab: one
bar per residue of a target chain, the higher the closer the residue
comes to the other chains (pymolx.interfaces.proximity). A chain too long
to show at once is scrolled with the slider below it (ChartScroll).
'''

import math

from pymol.Qt import QtCore, QtGui, QtWidgets

from pymolx import interfaces
from pymolx.gui import theme

Qt = QtCore.Qt

# theme.COLORS, from a low score to a high one
RAMP = ['proximity_0', 'proximity_1', 'proximity_2', 'proximity_3',
        'proximity_4']
# the score axis reaches this (2.5 A) at least, and the colors end there,
# so that two profiles can be compared (the viewer's colors too:
# interfaces.color_proximity)
SCORE_TOP = interfaces.PROXIMITY_SCORE_TOP
# distances named along the color bar, A
DISTANCE_TICKS = (2.5, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0)
# a break in the residue numbers takes this many slots at most
MAX_GAP = 20
# pixels a residue gets: at least SLOT_WIDTH, so that a long chain is
# scrolled instead of sharing pixels (until the whole chain is asked for),
# and MAX_SLOT_WIDTH at most; a step wider or narrower is this factor
SLOT_WIDTH = 1.0
MAX_SLOT_WIDTH = 24.0
ZOOM_STEP = 1.5
NUMBER_STEPS = (1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000,
                5000, 10000)


def residue_number(resi):
    '''
    :return: (number, insertion code) of "100A"; (None, resi) if it has
        no number
    '''
    digits = resi.rstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz')
    try:
        return int(digits), resi[len(digits):]
    except ValueError:
        return None, resi


def slots(resis):
    '''
    Where the residues go along the axis: one slot per residue number, so
    that missing residues leave a gap (of MAX_GAP slots at most). A
    residue with an insertion code, or numbered out of order, takes the
    next slot.

    :return: (slot of each residue, number of slots)
    '''
    result = []
    slot = -1
    last = None
    for resi in resis:
        number, _ = residue_number(resi)
        if last is not None and number is not None and number > last:
            slot += min(number - last, MAX_GAP)
        else:
            slot += 1
        if number is not None:
            last = number
        result.append(slot)
    return result, slot + 1


def number_step(pixels_per_slot, spacing=56.0):
    '''
    The residue numbers to label: multiples of this, at least "spacing"
    pixels apart
    '''
    for step in NUMBER_STEPS:
        if step * pixels_per_slot >= spacing:
            return step
    return NUMBER_STEPS[-1]


def ramp_color(fraction):
    '''
    :param fraction: 0 (far) to 1 (close)
    '''
    fraction = min(1.0, max(0.0, fraction)) * (len(RAMP) - 1)
    k = min(int(fraction), len(RAMP) - 2)
    a, b = (QtGui.QColor(theme.COLORS[name]) for name in RAMP[k:k + 2])
    t = fraction - k
    return QtGui.QColor(*(round(x + (y - x) * t) for x, y in zip(
        a.getRgb()[:3], b.getRgb()[:3])))


def distance_text(distance):
    return '%g Å' % round(distance, 1)


class ProximityChart(QtWidgets.QWidget):
    '''
    Bars over the residue numbers of one chain. Hovering names a residue
    and what it is closest to; a click emits residue_clicked. The mouse
    wheel scrolls a chain that doesn't fit, with Ctrl it makes the bars
    wider or narrower.
    '''

    residue_clicked = QtCore.Signal(str, str)  # chain, resi
    view_changed = QtCore.Signal()  # scrolled, bars wider or narrower

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName('proximity_chart')
        self.setMinimumHeight(190)
        self.setMouseTracking(True)
        self.rows = []       # (resi, resn, distance, partner)
        self.chain = ''
        self.caption = ''
        self.message = ''    # shown in place of an empty profile
        self.cutoff = interfaces.CONTACT_DISTANCE
        self.slot_of = []
        self.slot_count = 0
        self.hovered = None  # index into rows
        # pixels asked for per residue (0: the whole chain in view); a
        # residue never gets fewer than fill the plot (slot_width)
        self.zoom = SLOT_WIDTH
        self.offset = 0      # pixels scrolled

    def set_profile(self, rows, chain='', caption=''):
        '''
        :param rows: one chain's residues from interfaces.proximity
        :param caption: written above the bars
        '''
        rows = list(rows)
        same = chain == self.chain and \
            [r[0] for r in rows] == [r[0] for r in self.rows]
        self.rows = rows
        self.chain = chain
        self.caption = caption
        self.slot_of, self.slot_count = slots([r[0] for r in self.rows])
        self.hovered = None
        if not same:  # another chain: from its start, at the usual width
            self.zoom = SLOT_WIDTH
            self.offset = 0
        self.offset = min(self.offset, self.scroll_range())
        self.view_changed.emit()
        self.update()

    def set_cutoff(self, distance):
        '''
        The dashed line: residues above it are closer than this, A
        '''
        self.cutoff = float(distance)
        self.update()

    def set_message(self, text):
        self.message = text
        self.update()

    # geometry

    def scores(self):
        return [interfaces.proximity_score(r[2]) for r in self.rows]

    def top_score(self):
        highest = max(self.scores() + [SCORE_TOP])
        return math.ceil(highest * 10 - 1e-9) / 10

    def plot_rect(self):
        metrics = self.fontMetrics()
        left = metrics.horizontalAdvance('0.0') + 12
        right = 12 + 14 + 6 + metrics.horizontalAdvance('2.5 Å') + 8
        top = metrics.height() + 12
        bottom = metrics.height() + 12
        return QtCore.QRectF(left, top, max(10, self.width() - left - right),
                             max(10, self.height() - top - bottom))

    def _y(self, plot, score):
        return plot.bottom() - plot.height() * score / self.top_score()

    # scrolling, and how wide the bars are

    def fit_width(self):
        '''
        The pixels a residue gets when the whole chain fills the plot
        '''
        return self.plot_rect().width() / max(1, self.slot_count)

    def slot_width(self):
        return max(self.fit_width(), self.zoom)

    def scroll_range(self):
        '''
        :return: how many pixels of the chain lie outside the plot
        '''
        hidden = self.slot_width() * self.slot_count - \
            self.plot_rect().width()
        return max(0, int(round(hidden)))

    def set_offset(self, offset):
        '''
        Scroll: the plot starts this many pixels into the chain
        '''
        offset = min(max(0, int(round(offset))), self.scroll_range())
        if offset != self.offset:
            self.offset = offset
            self.hovered = None
            self.view_changed.emit()
            self.update()

    def set_zoom(self, pixels, anchor=None):
        '''
        :param pixels: per residue; 0 for the whole chain in view
        :param anchor: the x that keeps its residue (default: the middle)
        '''
        plot = self.plot_rect()
        if anchor is None:
            anchor = plot.center().x()
        slot = (anchor - plot.left() + self.offset) / self.slot_width()
        self.zoom = min(MAX_SLOT_WIDTH, max(0.0, float(pixels)))
        offset = slot * self.slot_width() - (anchor - plot.left())
        self.offset = min(max(0, int(round(offset))), self.scroll_range())
        self.hovered = None
        self.view_changed.emit()
        self.update()

    def can_zoom_in(self):
        return bool(self.rows) and self.slot_width() < MAX_SLOT_WIDTH - 1e-9

    def can_zoom_out(self):
        return bool(self.rows) and \
            self.slot_width() > self.fit_width() + 1e-9

    def zoom_in(self):
        self.set_zoom(self.slot_width() * ZOOM_STEP)

    def zoom_out(self):
        pixels = self.slot_width() / ZOOM_STEP
        self.set_zoom(pixels if pixels > self.fit_width() else 0.0)

    def show_whole(self):
        self.set_zoom(0.0)

    def shown(self):
        '''
        :return: indices of the residues in view
        '''
        plot = self.plot_rect()
        # (rounded to whole pixels, the last bar of a whole chain may
        # start on the plot's edge)
        return [index for index, bar in enumerate(self.bar_rects())
                if bar.right() > plot.left() and bar.left() <= plot.right()]

    def bar_rects(self):
        '''
        The bars of all residues, also of those scrolled out of view
        (whole pixels, so that bars of neighbours touch without
        overlapping)
        '''
        plot = self.plot_rect()
        width = self.slot_width()
        start = plot.left() - self.offset
        gap = 1 if width >= 4 else 0  # for every bar, or for none
        per_score = plot.height() / self.top_score()
        bars = []
        for slot, row in zip(self.slot_of, self.rows):
            left = round(start + slot * width)
            right = max(left + 1, round(start + (slot + 1) * width) - gap)
            height = per_score * interfaces.proximity_score(row[2])
            bars.append(QtCore.QRectF(left, plot.bottom() - height,
                                      right - left, height))
        return bars

    def bar_rect(self, index):
        return self.bar_rects()[index]

    def residue_at(self, x):
        '''
        :return: index of the residue under this x, or None; of residues
            sharing a pixel, the closest to the other chains
        '''
        if not self.rows:
            return None
        plot = self.plot_rect()
        if not plot.left() <= x <= plot.right():
            return None
        found = None
        for index, bar in enumerate(self.bar_rects()):
            if bar.left() - 0.5 <= x < bar.right() + 0.5 and (
                    found is None or
                    self.rows[index][2] < self.rows[found][2]):
                found = index
        return found

    def describe(self, index):
        '''
        The hover text of a residue
        '''
        resi, resn, distance, partner = self.rows[index]
        return '%s %s, chain %s\n%.1f Å from %s %s, chain %s\nScore %.2f' % (
            resn, resi, self.chain, distance, partner[2], partner[1],
            partner[0], interfaces.proximity_score(distance))

    def within_cutoff(self):
        return sum(1 for r in self.rows if r[2] <= self.cutoff)

    # painting

    def number_ticks(self, plot):
        '''
        :return: [(x, text)] of the residue numbers to label, in view
        '''
        width = self.slot_width()
        step = number_step(width)
        labelled = []
        for index, row in enumerate(self.rows):
            number, code = residue_number(row[0])
            if number is not None and not code and number % step == 0:
                labelled.append(index)
        if len(labelled) < 2:
            labelled = sorted({0, len(self.rows) - 1})
        ticks = []
        for index in labelled:
            x = plot.left() - self.offset + (self.slot_of[index] + 0.5) * width
            if plot.left() <= x <= plot.right():
                ticks.append((x, self.rows[index][0]))
        return ticks

    def draw(self, painter):
        color = {name: QtGui.QColor(value)
                 for name, value in theme.COLORS.items()}
        metrics = painter.fontMetrics()
        painter.fillRect(QtCore.QRectF(0, 0, self.width(), self.height()),
                         color['base'])
        if not self.rows:
            painter.setPen(color['text_dim'])
            painter.drawText(QtCore.QRectF(0, 0, self.width(), self.height()),
                             int(Qt.AlignmentFlag.AlignCenter), self.message)
            return
        plot = self.plot_rect()
        top = self.top_score()

        painter.setPen(color['text'])
        painter.drawText(QtCore.QRectF(plot.left(), 4, self.width() -
                                       plot.left() - 8, metrics.height()),
                         int(Qt.AlignmentFlag.AlignLeft |
                             Qt.AlignmentFlag.AlignVCenter), self.caption)

        # score axis
        step = 0.1 if top <= 0.6 else 0.2
        score = 0.0
        while score <= top + 1e-9:
            y = round(self._y(plot, score)) + 0.5
            painter.setPen(color['border'])
            painter.drawLine(QtCore.QPointF(plot.left(), y),
                             QtCore.QPointF(plot.right(), y))
            painter.setPen(color['text_dim'])
            painter.drawText(
                QtCore.QRectF(0, y - metrics.height() / 2, plot.left() - 6,
                              metrics.height()),
                int(Qt.AlignmentFlag.AlignRight |
                    Qt.AlignmentFlag.AlignVCenter), '%.1f' % score)
            score += step

        # bars: the closest last, so that it shows where bars share a
        # pixel; of a scrolled chain, the part in view
        painter.save()
        painter.setClipRect(QtCore.QRectF(plot.left(), 0, plot.width() + 1,
                                          self.height()))
        bars = self.bar_rects()
        if self.hovered is not None:
            bar = bars[self.hovered]
            painter.fillRect(QtCore.QRectF(bar.left(), plot.top(),
                                           bar.width(), plot.height()),
                             color['hover'])
        scores = self.scores()
        for index in sorted(self.shown(), key=scores.__getitem__):
            painter.fillRect(bars[index],
                             ramp_color(scores[index] / SCORE_TOP))
        if self.hovered is not None:
            painter.setPen(color['text_bright'])
            painter.drawRect(bars[self.hovered].adjusted(
                -0.5, -0.5, 0.5, 0.5))
        painter.restore()

        # residue numbers
        painter.setPen(color['text_dim'])
        base = round(plot.bottom()) + 0.5
        painter.drawLine(QtCore.QPointF(plot.left(), base),
                         QtCore.QPointF(plot.right(), base))
        taken = -1e9
        for x, text in self.number_ticks(plot):
            width = metrics.horizontalAdvance(text)
            left = min(max(x - width / 2, 0), self.width() - width)
            if left < taken + 8:
                continue
            taken = left + width
            painter.drawLine(QtCore.QPointF(round(x) + 0.5, base),
                             QtCore.QPointF(round(x) + 0.5, base + 4))
            painter.drawText(QtCore.QRectF(left, base + 6, width,
                                           metrics.height()),
                             int(Qt.AlignmentFlag.AlignCenter), text)

        # the line residues in contact rise above
        cutoff_score = interfaces.proximity_score(self.cutoff)
        if cutoff_score <= top:
            pen = QtGui.QPen(color['text_bright'])
            pen.setStyle(Qt.PenStyle.DashLine)
            painter.setPen(pen)
            y = round(self._y(plot, cutoff_score)) + 0.5
            painter.drawLine(QtCore.QPointF(plot.left(), y),
                             QtCore.QPointF(plot.right(), y))

        # color bar, with the distances the scores stand for
        bar = QtCore.QRectF(plot.right() + 12, plot.top(), 14, plot.height())
        gradient = QtGui.QLinearGradient(bar.bottomLeft(), bar.topLeft())
        for k in range(len(RAMP)):
            fraction = k / (len(RAMP) - 1)
            gradient.setColorAt(min(1.0, fraction * SCORE_TOP / top),
                                ramp_color(fraction))
        painter.fillRect(bar, QtGui.QBrush(gradient))
        placed = []
        distances = [self.cutoff] + [d for d in DISTANCE_TICKS
                                     if d != self.cutoff]
        for distance in distances:
            score = interfaces.proximity_score(distance)
            if score > top + 1e-9:
                continue
            y = self._y(plot, score)
            if any(abs(y - other) < metrics.height() for other in placed):
                continue
            placed.append(y)
            painter.setPen(color['text'] if distance == self.cutoff
                           else color['text_dim'])
            painter.drawLine(QtCore.QPointF(bar.right(), round(y) + 0.5),
                             QtCore.QPointF(bar.right() + 3, round(y) + 0.5))
            painter.drawText(
                QtCore.QRectF(bar.right() + 6, y - metrics.height() / 2,
                              self.width() - bar.right() - 6,
                              metrics.height()),
                int(Qt.AlignmentFlag.AlignLeft |
                    Qt.AlignmentFlag.AlignVCenter), distance_text(distance))

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        self.draw(painter)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.offset = min(self.offset, self.scroll_range())
        self.view_changed.emit()

    def image(self, scale=2):
        '''
        The chart as a picture, "scale" times as fine as on screen
        '''
        picture = QtGui.QImage(self.size() * scale,
                               QtGui.QImage.Format.Format_ARGB32)
        painter = QtGui.QPainter(picture)
        painter.setFont(self.font())
        painter.scale(scale, scale)
        hovered, self.hovered = self.hovered, None
        try:
            self.draw(painter)
        finally:
            self.hovered = hovered
            painter.end()
        return picture

    # mouse

    @staticmethod
    def _position(event):
        return event.position() if hasattr(event, 'position') else event.pos()

    def mouseMoveEvent(self, event):
        index = self.residue_at(self._position(event).x())
        if index != self.hovered:
            self.hovered = index
            self.update()
        if index is None:
            QtWidgets.QToolTip.hideText()
        else:
            point = event.globalPosition().toPoint() if hasattr(
                event, 'globalPosition') else event.globalPos()
            QtWidgets.QToolTip.showText(point, self.describe(index), self)

    def leaveEvent(self, event):
        if self.hovered is not None:
            self.hovered = None
            self.update()
        super().leaveEvent(event)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            index = self.residue_at(self._position(event).x())
            if index is not None:
                self.residue_clicked.emit(self.chain, self.rows[index][0])
        super().mousePressEvent(event)

    def wheelEvent(self, event):
        delta = event.angleDelta()
        amount = delta.y() or delta.x()
        if not self.rows or not amount:
            super().wheelEvent(event)
        elif event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            self.set_zoom(self.slot_width() * (
                ZOOM_STEP if amount > 0 else 1 / ZOOM_STEP),
                self._position(event).x())
            event.accept()
        elif self.scroll_range():
            self.set_offset(self.offset - amount)
            event.accept()
        else:
            super().wheelEvent(event)


class ChartScroll(QtWidgets.QWidget):
    '''
    The row under a chart: the slider that moves along a chain too long
    to show at once, and the buttons for wider and narrower bars
    '''

    def __init__(self, chart, parent=None):
        super().__init__(parent)
        self.chart = chart
        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.slider = QtWidgets.QScrollBar(Qt.Orientation.Horizontal, self)
        self.slider.setObjectName('proximity_slider')
        self.slider.setToolTip('Move along the chain (also the mouse wheel '
                               'over the graph)')
        self.slider.valueChanged.connect(chart.set_offset)
        layout.addWidget(self.slider, 1)
        self.narrower = QtWidgets.QPushButton('−', self)
        self.narrower.setToolTip('Narrower bars: more residues in view '
                                 '(also Ctrl + mouse wheel)')
        self.narrower.clicked.connect(chart.zoom_out)
        self.wider = QtWidgets.QPushButton('+', self)
        self.wider.setToolTip('Wider bars: fewer residues in view '
                              '(also Ctrl + mouse wheel)')
        self.wider.clicked.connect(chart.zoom_in)
        for button in (self.narrower, self.wider):
            button.setFixedWidth(32)
            layout.addWidget(button)
        self.whole = QtWidgets.QPushButton('Whole Chain', self)
        self.whole.setToolTip('All the chain\'s residues in view')
        self.whole.clicked.connect(chart.show_whole)
        layout.addWidget(self.whole)
        chart.view_changed.connect(self.follow)
        self.follow()

    def follow(self):
        '''
        Take the slider's place and length, and what the buttons can
        still do, from the chart
        '''
        chart = self.chart
        slider = self.slider
        slider.blockSignals(True)
        # (Qt draws nothing of a slider without a range: with one of a
        # single pixel, a chain that fits shows as a handle from end to
        # end, which set_offset keeps in its place)
        slider.setRange(0, max(1, chart.scroll_range()))
        slider.setPageStep(max(1, int(chart.plot_rect().width())))
        slider.setSingleStep(max(1, int(round(chart.slot_width() * 10))))
        slider.setValue(chart.offset)
        slider.blockSignals(False)
        self.wider.setEnabled(chart.can_zoom_in())
        self.narrower.setEnabled(chart.can_zoom_out())
        self.whole.setEnabled(chart.can_zoom_out())
