'''
Structure Comparison dialog: choose the structures and the reference,
what to compare and how; results as a table per structure, a profile of
the similarity along the alignment and a color legend.
'''

from pymol import cmd
from pymol.Qt import QtCore, QtGui, QtWidgets

from pymolx import compare
from pymolx.gui import theme

Qt = QtCore.Qt

_dialog = None

METHODS = [('align (sequence-guided)', 'align'), ('super', 'super'),
           ('cealign (low identity)', 'cealign')]
METRICS = [('lDDT (local, 0-100)', 'lddt'),
           ('Distance after superposition (Å)', 'distance')]
MODES = [('Each to the reference', 'reference'),
         ('All against all (per position)', 'all')]
# profile lines, one per structure
LINE_COLORS = [theme.COLORS['plot_%d' % k] for k in range(10)]


def _combo(parent, items):
    combo = QtWidgets.QComboBox(parent)
    for label, value in items:
        combo.addItem(label, value)
    return combo


def _gradient(widget, red, blue):
    '''
    the color ramp, red (left) to blue (right)
    '''
    gradient = QtGui.QLinearGradient(0, 0, widget.width(), 0)
    for k in range(11):
        r, g, b = compare.ramp_color(k / 10)
        gradient.setColorAt(k / 10, QtGui.QColor.fromRgbF(r, g, b))
    return gradient


class Legend(QtWidgets.QWidget):
    '''
    The color bar with its scores
    '''

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(34)
        self.red, self.blue, self.unit = 50.0, 100.0, ''

    def set_range(self, red, blue, unit):
        self.red, self.blue, self.unit = red, blue, unit
        self.update()

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        bar = QtCore.QRectF(0, 0, self.width(), 14)
        painter.fillRect(bar, QtGui.QBrush(_gradient(self, self.red,
                                                     self.blue)))
        painter.setPen(QtGui.QColor(theme.COLORS['text']))
        text = QtCore.QRectF(0, 16, self.width(), self.height() - 16)
        painter.drawText(text, int(Qt.AlignmentFlag.AlignLeft),
                         '%g%s  different' % (self.red, self.unit))
        painter.drawText(text, int(Qt.AlignmentFlag.AlignRight),
                         'similar  %g%s' % (self.blue, self.unit))
        painter.drawText(text, int(Qt.AlignmentFlag.AlignHCenter),
                         'grey: no partner')


class Profile(QtWidgets.QWidget):
    '''
    Similarity along the alignment, one line per structure; the colored
    band below shows the reference's colors
    '''

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName('compare_profile')
        self.setMinimumHeight(140)
        self.setMouseTracking(True)
        self.result = None
        self.red, self.blue = 50.0, 100.0

    def set_result(self, result, red, blue):
        self.result, self.red, self.blue = result, red, blue
        self.update()

    def _frame(self):
        return QtCore.QRectF(36, 6, max(10, self.width() - 44),
                             max(10, self.height() - 34))

    def _y(self, frame, score):
        fraction = compare.similarity_fraction(score, self.result.metric,
                                               self.red, self.blue)
        return frame.bottom() - fraction * frame.height()

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.fillRect(self.rect(), QtGui.QColor(theme.COLORS['base']))
        if not self.result or not self.result.columns:
            return
        frame = self._frame()
        n = len(self.result.columns)
        text = QtGui.QColor(theme.COLORS['text_dim'])
        painter.setPen(text)
        painter.drawRect(frame)
        for score in (self.red, self.blue):
            painter.drawText(QtCore.QRectF(0, self._y(frame, score) - 8, 32,
                                           16),
                             int(Qt.AlignmentFlag.AlignRight), '%g' % score)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        for k, obj in enumerate(self.result.objects):
            points = [QtCore.QPointF(
                frame.left() + frame.width() * column / max(1, n - 1),
                self._y(frame, score))
                for column, score in self.result.profile(obj)]
            pen = QtGui.QPen(QtGui.QColor(LINE_COLORS[k % len(LINE_COLORS)]))
            pen.setWidthF(1.4 if obj != self.result.reference else 2.2)
            painter.setPen(pen)
            for a, b in zip(points, points[1:]):
                painter.drawLine(a, b)
        # the reference's colors along the alignment
        band = QtCore.QRectF(frame.left(), frame.bottom() + 6,
                             frame.width(), 10)
        ref = self.result.reference
        scores = dict(self.result.profile(ref))
        width = band.width() / n
        for column in range(n):
            score = scores.get(column)
            if score is None:
                color = QtGui.QColor(theme.COLORS['border'])
            else:
                color = QtGui.QColor.fromRgbF(*compare.ramp_color(
                    compare.similarity_fraction(score, self.result.metric,
                                                self.red, self.blue)))
            painter.fillRect(QtCore.QRectF(band.left() + column * width,
                                           band.top(), width + 0.5, 10),
                             color)

    def mouseMoveEvent(self, event):
        if not self.result or not self.result.columns:
            return
        frame = self._frame()
        pos = event.position() if hasattr(event, 'position') else event.pos()
        n = len(self.result.columns)
        column = round((pos.x() - frame.left()) / frame.width() * (n - 1))
        if not 0 <= column < n:
            return
        lines = []
        for obj, index in self.result.columns[column].items():
            chain, resi, resn, _ = self.result.residues[obj][index]
            score = self.result.scores[obj].get(index)
            lines.append('%s %s/%s %s: %s' % (
                obj, chain or '-', resn, resi,
                '-' if score is None else '%.1f' % score))
        self.setToolTip('\n'.join(lines))


class CompareDialog(QtWidgets.QDialog):

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent)
        self.cmd = _self
        self.result = None
        self.setObjectName('compare_dialog')
        self.setWindowTitle('Structure Comparison')
        self.resize(940, 680)

        layout = QtWidgets.QHBoxLayout(self)
        left = QtWidgets.QVBoxLayout()
        layout.addLayout(left, 2)
        right = QtWidgets.QVBoxLayout()
        layout.addLayout(right, 3)

        # structures
        box = QtWidgets.QGroupBox('Structures', self)
        box.setObjectName('prep_step')
        column = QtWidgets.QVBoxLayout(box)
        self.filter = QtWidgets.QLineEdit(box)
        self.filter.setPlaceholderText('filter, e.g. *mut*')
        self.filter.textChanged.connect(self._filter)
        column.addWidget(self.filter)
        self.models = QtWidgets.QListWidget(box)
        self.models.setObjectName('compare_models')
        self.models.itemChanged.connect(self._fill_reference)
        column.addWidget(self.models, 1)
        buttons = QtWidgets.QHBoxLayout()
        for label, state in (('All', True), ('None', False)):
            button = QtWidgets.QPushButton(label, box)
            button.clicked.connect(lambda _=False, s=state: self._check(s))
            buttons.addWidget(button)
        refresh = QtWidgets.QPushButton('Refresh', box)
        refresh.clicked.connect(self.refresh)
        buttons.addWidget(refresh)
        column.addLayout(buttons)
        left.addWidget(box, 1)

        # settings
        box = QtWidgets.QGroupBox('Settings', self)
        box.setObjectName('prep_step')
        grid = QtWidgets.QGridLayout(box)
        rows = []
        self.reference = QtWidgets.QComboBox(box)
        rows.append(('Reference', self.reference))
        self.selection = QtWidgets.QLineEdit(box)
        self.selection.setPlaceholderText('polymer.protein')
        self.selection.setToolTip('What to compare in every structure, e.g. '
                                  '"chain A"')
        rows.append(('Compare', self.selection))
        self.method = _combo(box, METHODS)
        rows.append(('Alignment', self.method))
        self.metric_box = _combo(box, METRICS)
        self.metric_box.setToolTip(
            'lDDT: how well each residue keeps its distances to its '
            'neighbors (15 Å), independent of the superposition: moved '
            'domains stay blue. Distance: after superposing the cores.')
        self.metric_box.currentIndexChanged.connect(self._metric_changed)
        rows.append(('Similarity', self.metric_box))
        self.mode = _combo(box, MODES)
        rows.append(('Colors', self.mode))
        scale = QtWidgets.QHBoxLayout()
        self.red = QtWidgets.QDoubleSpinBox(box)
        self.blue = QtWidgets.QDoubleSpinBox(box)
        for spin in (self.red, self.blue):
            spin.setRange(0.0, 100.0)
            spin.setDecimals(1)
            spin.valueChanged.connect(self._recolor)
        scale.addWidget(QtWidgets.QLabel('red at', box))
        scale.addWidget(self.red)
        scale.addWidget(QtWidgets.QLabel('blue at', box))
        scale.addWidget(self.blue)
        for row, (label, widget) in enumerate(rows):
            grid.addWidget(QtWidgets.QLabel(label, box), row, 0)
            grid.addWidget(widget, row, 1)
        grid.addWidget(QtWidgets.QLabel('Scale', box), len(rows), 0)
        grid.addLayout(scale, len(rows), 1)
        left.addWidget(box)

        run = QtWidgets.QHBoxLayout()
        self.status = QtWidgets.QLabel(self)
        self.status.setObjectName('prep_note')
        self.status.setWordWrap(True)
        run.addWidget(self.status, 1)
        self.run_button = QtWidgets.QPushButton('Compare', self)
        self.run_button.setDefault(True)
        self.run_button.clicked.connect(self.run)
        run.addWidget(self.run_button)
        left.addLayout(run)

        # results
        self.table = QtWidgets.QTableWidget(0, 5, self)
        self.table.setObjectName('compare_table')
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.horizontalHeader().setStretchLastSection(True)
        right.addWidget(self.table, 2)
        self.profile = Profile(self)
        right.addWidget(self.profile, 3)
        note = QtWidgets.QLabel('Similarity along the alignment (hover for '
                                'residues); band: the reference\'s colors',
                                self)
        note.setObjectName('prep_note')
        right.addWidget(note)
        self.legend = Legend(self)
        right.addWidget(self.legend)
        actions = QtWidgets.QHBoxLayout()
        actions.addWidget(QtWidgets.QLabel('Select residues worse than', self))
        self.threshold = QtWidgets.QDoubleSpinBox(self)
        self.threshold.setRange(0.0, 100.0)
        self.threshold.setDecimals(1)
        actions.addWidget(self.threshold)
        self.select_button = QtWidgets.QPushButton('Select', self)
        self.select_button.setToolTip('Selection "dissimilar"')
        self.select_button.clicked.connect(self.select_dissimilar)
        actions.addWidget(self.select_button)
        actions.addStretch(1)
        self.save_button = QtWidgets.QPushButton('Save CSV...', self)
        self.save_button.clicked.connect(self.save_csv)
        actions.addWidget(self.save_button)
        right.addLayout(actions)
        self._metric_changed()
        self._enable_results(False)

    # structures

    def refresh(self):
        checked = set(self.checked_models()) if self.models.count() else None
        self.models.blockSignals(True)
        self.models.clear()
        for name in self.cmd.get_object_list('all') or []:
            if self.cmd.get_type(name) != 'object:molecule' or \
                    name.startswith('_'):
                continue
            item = QtWidgets.QListWidgetItem(name, self.models)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked if checked is None or
                               name in checked else Qt.CheckState.Unchecked)
        self.models.blockSignals(False)
        self._filter(self.filter.text())

    def _filter(self, text):
        import fnmatch
        pattern = text.strip()
        if pattern and not any(ch in pattern for ch in '*?['):
            pattern = '*%s*' % pattern
        for row in range(self.models.count()):
            item = self.models.item(row)
            item.setHidden(bool(pattern) and not fnmatch.fnmatch(
                item.text(), pattern))
        self._fill_reference()

    def _check(self, state):
        self.models.blockSignals(True)
        for row in range(self.models.count()):
            item = self.models.item(row)
            if not item.isHidden():
                item.setCheckState(Qt.CheckState.Checked if state else
                                   Qt.CheckState.Unchecked)
        self.models.blockSignals(False)
        self._fill_reference()

    def checked_models(self):
        return [self.models.item(row).text()
                for row in range(self.models.count())
                if not self.models.item(row).isHidden() and
                self.models.item(row).checkState() == Qt.CheckState.Checked]

    def _fill_reference(self, *_):
        current = self.reference.currentText()
        models = self.checked_models()
        self.reference.clear()
        self.reference.addItems(models)
        if current in models:
            self.reference.setCurrentText(current)

    def showEvent(self, event):
        self.refresh()
        super().showEvent(event)

    # settings

    def _metric_changed(self, *_):
        metric = self.metric_box.currentData()
        red, blue = compare.RANGES[metric]
        for spin, value in ((self.red, red), (self.blue, blue)):
            spin.blockSignals(True)
            spin.setSuffix(' Å' if metric == 'distance' else '')
            spin.setValue(value)
            spin.blockSignals(False)
        self.threshold.setSuffix(' Å' if metric == 'distance' else '')
        self.threshold.setValue(80.0 if metric == 'lddt' else 2.0)
        self.legend.set_range(red, blue, ' Å' if metric == 'distance'
                              else '')

    # running

    def run(self):
        models = self.checked_models()
        if len(models) < 2:
            self.status.setText('Choose at least two structures')
            return
        QtWidgets.QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            self.result = compare.compute(
                ' '.join(models), self.reference.currentText(),
                self.selection.text().strip() or 'polymer.protein',
                self.method.currentData(), self.metric_box.currentData(),
                self.mode.currentData(), _self=self.cmd)
            compare.apply(self.result, self.red.value(), self.blue.value(),
                          _self=self.cmd)
        except Exception as e:
            self.status.setText('Comparison failed: %s' % e)
            return
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        r = self.result
        self.status.setText('%d structures on %s, %d aligned positions' % (
            len(r.objects), r.reference, len(r.columns)))
        print(r.table())
        self._fill_table()
        self.profile.set_result(r, self.red.value(), self.blue.value())
        self.legend.set_range(self.red.value(), self.blue.value(),
                              ' Å' if r.metric == 'distance' else '')
        self._enable_results(True)

    def _recolor(self, *_):
        unit = ' Å' if self.metric_box.currentData() == 'distance' else ''
        self.legend.set_range(self.red.value(), self.blue.value(), unit)
        if self.result and self.result.metric == self.metric_box.currentData():
            compare.apply(self.result, self.red.value(), self.blue.value(),
                          _self=self.cmd)
            self.profile.set_result(self.result, self.red.value(),
                                    self.blue.value())

    def _fill_table(self):
        r = self.result
        label = 'Mean lDDT' if r.metric == 'lddt' else 'Mean distance (Å)'
        self.table.setHorizontalHeaderLabels(
            ['Structure', 'Aligned', 'Identity', 'RMSD (Å)', label])
        self.table.setRowCount(len(r.objects))
        for row, obj in enumerate(r.objects):
            s = r.stats[obj]
            values = [obj + (' (reference)' if obj == r.reference else ''),
                      str(s['aligned']),
                      '%.0f%%' % s['identity'] if s['identity'] is not None
                      else '-',
                      '%.2f' % s['rmsd'] if s['rmsd'] is not None else '-',
                      '%.1f' % s['score'] if s['score'] is not None else '-']
            for col, value in enumerate(values):
                item = QtWidgets.QTableWidgetItem(value)
                if col == 0:
                    item.setForeground(QtGui.QBrush(QtGui.QColor(
                        LINE_COLORS[row % len(LINE_COLORS)])))
                self.table.setItem(row, col, item)
        self.table.resizeColumnsToContents()

    def _enable_results(self, enabled):
        for widget in (self.select_button, self.save_button,
                       self.threshold):
            widget.setEnabled(enabled)

    def select_dissimilar(self):
        if not self.result:
            return
        count = compare.select_dissimilar(self.result, self.threshold.value(),
                                          _self=self.cmd)
        self.status.setText('Selected %d residues ("dissimilar")' % count)

    def save_csv(self):
        if not self.result:
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save similarity', 'similarity.csv', 'CSV (*.csv)')
        if path:
            self.result.save_csv(path)
            self.status.setText('Saved %s' % path)


def show_dialog():
    global _dialog
    if _dialog is None:
        from pymol.gui import get_qtwindow
        _dialog = CompareDialog(get_qtwindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
