'''
Sequence Comparison dialog: choose the structures, the reference and what
to compare; results as the colored alignment, a percent identity matrix
and a summary per chain.
'''

import html

from pymol import cmd
from pymol.Qt import QtCore, QtGui, QtWidgets

from pymolx import compare, seqcompare
from pymolx.gui import theme

Qt = QtCore.Qt

_dialog = None

MODES = [('Each to the reference', 'reference'),
         ('All against all (per position)', 'all')]
BLOCK = 60          # alignment columns per line
LIGHTEN = 0.45      # residue backgrounds: score colors toward white


def _background(score):
    '''
    "#rrggbb" for a residue's score (None: no partner)
    '''
    if score is None:
        return theme.COLORS['border']
    r, g, b = compare.ramp_color(score / 100.0)
    return '#%02x%02x%02x' % tuple(
        round(255 * (c + (1 - c) * LIGHTEN)) for c in (r, g, b))


def alignment_html(result):
    '''
    The alignment, colored, in blocks of BLOCK columns with the
    reference's residue numbers above
    '''
    width = max(len(e.name) + len(' (ref)') for e in result.entries)
    dark = theme.COLORS['base']
    out = []
    for group in result.groups:
        ref = group.reference
        n = len(group.columns)
        for start in range(0, n, BLOCK):
            stop = min(n, start + BLOCK)
            # ruler: the reference's residue number every 10 columns
            ruler = [' '] * (stop - start)
            for k in range(start, stop):
                pos = group.columns[k].get(ref.name)
                if pos is not None and (k - start) % 10 == 0:
                    label = ref.residues[pos][0]
                    for j, ch in enumerate(label):
                        if k - start + j < len(ruler):
                            ruler[k - start + j] = ch
            out.append('<span style="color:%s">%s%s</span>' % (
                theme.COLORS['text_dim'], ' ' * (width + 1),
                html.escape(''.join(ruler))))
            for e in group.entries:
                label = e.name + (' (ref)' if e is ref else '')
                cells = []
                for k in range(start, stop):
                    pos = group.columns[k].get(e.name)
                    if pos is None:
                        cells.append('-')
                        continue
                    resi, resn, letter, index = e.residues[pos]
                    letter = letter if index is not None else letter.lower()
                    score = result.scores[e.name].get(pos)
                    cells.append(
                        '<span style="background-color:%s;color:%s" '
                        'title="%s %s %s">%s</span>' % (
                            _background(score), dark, e.name, resn, resi,
                            letter))
                out.append('%s %s' % (html.escape(label.ljust(width)),
                                      ''.join(cells)))
            out.append('')
    return '<pre style="font-family:%s">%s</pre>' % (
        html.escape(theme.console_font().family()), '\n'.join(out))


class Legend(QtWidgets.QWidget):

    def __init__(self, parent=None):
        super().__init__(parent)
        row = QtWidgets.QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        for score, label in seqcompare.CATEGORIES + [(None, 'no partner')]:
            swatch = QtWidgets.QLabel(self)
            swatch.setFixedSize(14, 14)
            swatch.setStyleSheet('background-color:%s' % _background(score))
            row.addWidget(swatch)
            row.addWidget(QtWidgets.QLabel(label.split(' (')[0], self))
            row.addSpacing(8)
        row.addStretch(1)


class SequenceDialog(QtWidgets.QDialog):

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent)
        self.cmd = _self
        self.result = None
        self.setObjectName('seqcompare_dialog')
        self.setWindowTitle('Sequence Comparison')
        self.resize(980, 660)

        layout = QtWidgets.QHBoxLayout(self)
        left = QtWidgets.QVBoxLayout()
        layout.addLayout(left, 2)
        right = QtWidgets.QVBoxLayout()
        layout.addLayout(right, 5)

        box = QtWidgets.QGroupBox('Structures', self)
        box.setObjectName('prep_step')
        column = QtWidgets.QVBoxLayout(box)
        self.filter = QtWidgets.QLineEdit(box)
        self.filter.setPlaceholderText('filter, e.g. *mut*')
        self.filter.textChanged.connect(self._filter)
        column.addWidget(self.filter)
        self.models = QtWidgets.QListWidget(box)
        self.models.setObjectName('seqcompare_models')
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

        box = QtWidgets.QGroupBox('Settings', self)
        box.setObjectName('prep_step')
        grid = QtWidgets.QGridLayout(box)
        self.reference = QtWidgets.QComboBox(box)
        self.selection = QtWidgets.QLineEdit(box)
        self.selection.setPlaceholderText('polymer.protein')
        self.selection.setToolTip('What to compare, e.g. "chain H"; chains '
                                  'are matched by chain ID')
        self.mode = QtWidgets.QComboBox(box)
        for label, value in MODES:
            self.mode.addItem(label, value)
        for row, (label, widget) in enumerate((
                ('Reference', self.reference), ('Compare', self.selection),
                ('Colors', self.mode))):
            grid.addWidget(QtWidgets.QLabel(label, box), row, 0)
            grid.addWidget(widget, row, 1)
        self.superpose = QtWidgets.QCheckBox(
            'Superimpose on the reference by the alignment', box)
        grid.addWidget(self.superpose, 3, 0, 1, 2)
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

        self.tabs = QtWidgets.QTabWidget(self)
        self.view = QtWidgets.QTextBrowser(self)
        self.view.setObjectName('seqcompare_alignment')
        self.view.setLineWrapMode(QtWidgets.QTextEdit.LineWrapMode.NoWrap)
        self.view.setFont(theme.console_font())
        self.tabs.addTab(self.view, 'Alignment')
        self.matrix = QtWidgets.QTableWidget(self)
        self.matrix.setObjectName('seqcompare_identity')
        self.matrix.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.tabs.addTab(self.matrix, 'Identity (%)')
        self.summary = QtWidgets.QTableWidget(0, 6, self)
        self.summary.setObjectName('seqcompare_summary')
        self.summary.setHorizontalHeaderLabels(
            ['Chain', 'Length', 'Aligned', 'Identity', 'Similarity', 'Gaps'])
        self.summary.verticalHeader().setVisible(False)
        self.summary.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.summary.horizontalHeader().setStretchLastSection(True)
        self.tabs.addTab(self.summary, 'Summary')
        right.addWidget(self.tabs, 1)
        right.addWidget(Legend(self))
        note = QtWidgets.QLabel(
            'Similar, neutral, different: BLOSUM62 above, at, below 0. Lower '
            'case: residues without coordinates. Hover a residue for its '
            'number.', self)
        note.setObjectName('prep_note')
        note.setWordWrap(True)
        right.addWidget(note)
        actions = QtWidgets.QHBoxLayout()
        actions.addStretch(1)
        self.fasta_button = QtWidgets.QPushButton('Save FASTA...', self)
        self.fasta_button.clicked.connect(self.save_fasta)
        actions.addWidget(self.fasta_button)
        self.csv_button = QtWidgets.QPushButton('Save CSV...', self)
        self.csv_button.clicked.connect(self.save_csv)
        actions.addWidget(self.csv_button)
        right.addLayout(actions)
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

    # running

    def run(self):
        models = self.checked_models()
        if len(models) < 2:
            self.status.setText('Choose at least two structures')
            return
        QtWidgets.QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            self.result = seqcompare.compute(
                ' '.join(models), self.reference.currentText(),
                self.selection.text().strip() or 'polymer.protein',
                self.mode.currentData(), _self=self.cmd)
            seqcompare.apply(self.result, self.superpose.isChecked(),
                             _self=self.cmd)
        except Exception as e:
            self.status.setText('Comparison failed: %s' % e)
            return
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        r = self.result
        self.status.setText('%d chains in %d group%s on %s' % (
            len(r.entries), len(r.groups), '' if len(r.groups) == 1 else 's',
            r.reference))
        print(r.table())
        self.view.setHtml(alignment_html(r))
        self._fill_matrix()
        self._fill_summary()
        self._enable_results(True)

    def _fill_matrix(self):
        names, values = self.result.identity_matrix()
        self.matrix.setRowCount(len(names))
        self.matrix.setColumnCount(len(names))
        self.matrix.setHorizontalHeaderLabels(names)
        self.matrix.setVerticalHeaderLabels(names)
        for i, a in enumerate(names):
            for j, b in enumerate(names):
                value = values.get((a, b))
                item = QtWidgets.QTableWidgetItem(
                    '' if value is None else '%.1f' % value)
                item.setTextAlignment(int(Qt.AlignmentFlag.AlignRight |
                                          Qt.AlignmentFlag.AlignVCenter))
                self.matrix.setItem(i, j, item)
        self.matrix.resizeColumnsToContents()

    def _fill_summary(self):
        r = self.result
        self.summary.setRowCount(len(r.entries))
        for row, e in enumerate(r.entries):
            s = r.stats[e.name]
            values = [e.name + (' (reference)' if s['reference'] else ''),
                      str(s['length']), str(s['aligned']),
                      '-' if s['identity'] is None else
                      '%.1f%%' % s['identity'],
                      '-' if s['similarity'] is None else
                      '%.1f%%' % s['similarity'],
                      '-' if s['gaps'] is None else str(s['gaps'])]
            for col, value in enumerate(values):
                self.summary.setItem(row, col,
                                     QtWidgets.QTableWidgetItem(value))
        self.summary.resizeColumnsToContents()

    def _enable_results(self, enabled):
        self.fasta_button.setEnabled(enabled)
        self.csv_button.setEnabled(enabled)

    def save_fasta(self):
        if not self.result:
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save alignment', 'alignment.fasta',
            'FASTA (*.fasta *.fa);;All files (*)')
        if path:
            with open(path, 'w') as handle:
                handle.write(self.result.fasta())
            self.status.setText('Saved %s' % path)

    def save_csv(self):
        if not self.result:
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save scores', 'sequence_similarity.csv', 'CSV (*.csv)')
        if path:
            self.result.save_csv(path)
            self.status.setText('Saved %s' % path)


def show_dialog():
    global _dialog
    if _dialog is None:
        from pymol.gui import get_qtwindow
        _dialog = SequenceDialog(get_qtwindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
