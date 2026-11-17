'''
Paratope Prediction dialog: an antibody's binding residues predicted by
Parapred (from the sequence) or Paragraph (from the structure), as a
table, a selection and colour on the model (pymolx.paratope)
'''

from pymol import cmd
from pymol.Qt import QtCore, QtWidgets

from pymolx import paratope

Qt = QtCore.Qt

_dialog = None

COLUMNS = ['Chain', 'Residue', 'Score', 'Binds', 'Region']


class ParatopeDialog(QtWidgets.QDialog):

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent)
        self.cmd = _self
        self.residues = []
        self.setObjectName('paratope_dialog')
        self.setWindowTitle('Paratope Prediction')
        self.resize(900, 620)

        layout = QtWidgets.QHBoxLayout(self)
        left = QtWidgets.QVBoxLayout()
        layout.addLayout(left, 2)
        right = QtWidgets.QVBoxLayout()
        layout.addLayout(right, 3)

        box = QtWidgets.QGroupBox('Antibody', self)
        box.setObjectName('prep_step')
        grid = QtWidgets.QGridLayout(box)
        self.source = QtWidgets.QButtonGroup(self)
        structure = QtWidgets.QRadioButton('A structure', box)
        typed = QtWidgets.QRadioButton('A sequence', box)
        structure.setChecked(True)
        self.source.addButton(structure, 0)
        self.source.addButton(typed, 1)
        grid.addWidget(structure, 0, 0)
        grid.addWidget(typed, 0, 1)
        self.source.idToggled.connect(self._source_changed)
        grid.addWidget(QtWidgets.QLabel('Object', box), 1, 0)
        self.object = QtWidgets.QComboBox(box)
        self.object.setObjectName('paratope_object')
        self.object.currentIndexChanged.connect(self._object_changed)
        grid.addWidget(self.object, 1, 1)
        grid.addWidget(QtWidgets.QLabel('Chain', box), 2, 0)
        self.chain = QtWidgets.QComboBox(box)
        grid.addWidget(self.chain, 2, 1)
        self.sequence = QtWidgets.QPlainTextEdit(box)
        self.sequence.setObjectName('paratope_sequence')
        self.sequence.setPlaceholderText('an antibody chain, e.g.\n'
                                         'EVQLVESGGG...')
        self.sequence.setFont(_console_font())
        self.sequence.setMaximumHeight(110)
        grid.addWidget(self.sequence, 3, 0, 1, 2)
        grid.setColumnStretch(1, 1)
        left.addWidget(box)

        box = QtWidgets.QGroupBox('Predictor', self)
        box.setObjectName('prep_step')
        column = QtWidgets.QVBoxLayout(box)
        self.method = QtWidgets.QComboBox(box)
        self.method.setObjectName('paratope_method')
        for label, value in paratope.METHODS:
            self.method.addItem(label, value)
        self.method.currentIndexChanged.connect(self._method_changed)
        column.addWidget(self.method)
        self.about = QtWidgets.QLabel(box)
        self.about.setObjectName('prep_note')
        self.about.setWordWrap(True)
        column.addWidget(self.about)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel('Binds at or above', box))
        self.threshold = QtWidgets.QDoubleSpinBox(box)
        self.threshold.setRange(0.05, 0.99)
        self.threshold.setSingleStep(0.05)
        self.threshold.setDecimals(2)
        self.threshold.setValue(paratope.THRESHOLD)
        self.threshold.valueChanged.connect(self._fill_table)
        row.addWidget(self.threshold)
        row.addStretch(1)
        column.addLayout(row)
        self.colour = QtWidgets.QCheckBox('Colour the model and set the '
                                          'B-factors', box)
        self.colour.setChecked(True)
        self.select = QtWidgets.QCheckBox('Selection of the paratope', box)
        self.select.setChecked(True)
        column.addWidget(self.colour)
        column.addWidget(self.select)
        left.addWidget(box)
        left.addStretch(1)

        self.run_button = QtWidgets.QPushButton('Predict', self)
        self.run_button.setDefault(True)
        self.run_button.clicked.connect(self.run)
        left.addWidget(self.run_button)

        self.table = QtWidgets.QTableWidget(0, len(COLUMNS), self)
        self.table.setObjectName('paratope_table')
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.horizontalHeader().setStretchLastSection(True)
        right.addWidget(self.table, 1)
        self.status = QtWidgets.QLabel(self)
        self.status.setObjectName('prep_note')
        self.status.setWordWrap(True)
        right.addWidget(self.status)
        copy = QtWidgets.QPushButton('Copy Table', self)
        copy.clicked.connect(self.copy)
        right.addWidget(copy)

        self._source_changed()
        self._method_changed()
        self.refresh()

    # the session

    def refresh(self):
        names = [n for n in self.cmd.get_names('objects') or []
                 if self.cmd.get_type(n) == 'object:molecule']
        current = self.object.currentText()
        self.object.blockSignals(True)
        self.object.clear()
        self.object.addItems(names)
        if current in names:
            self.object.setCurrentText(current)
        self.object.blockSignals(False)
        self._object_changed()

    def _object_changed(self, *_):
        name = self.object.currentText()
        self.chain.clear()
        self.chain.addItem('All antibody chains', '')
        for one in (self.cmd.get_chains('%' + name) if name else []) or []:
            self.chain.addItem(one, one)

    def _source_changed(self, *_):
        structure = self.source.checkedId() == 0
        self.object.setEnabled(structure)
        self.chain.setEnabled(structure)
        self.sequence.setEnabled(not structure)
        self._method_changed()

    def _method_changed(self, *_):
        method = self.method.currentData()
        missing = paratope.available(method)
        if method == 'paragraph':
            self.about.setText(
                'Paragraph (Chinery et al. 2023): a graph network over '
                'the Fv\'s structure. The chains are IMGT numbered with '
                'ANARCI first.')
        else:
            self.about.setText(
                'Parapred (Liberis et al. 2018, the PyTorch port): a '
                'network over each CDR with two residues either side, '
                'from the sequence alone.')
        structure = self.source.checkedId() == 0
        trouble = missing or ('' if structure else
                              'Paragraph needs a structure'
                              if method == 'paragraph' else '')
        self.run_button.setEnabled(not trouble)
        if trouble:
            self.status.setText(trouble)
        elif self.status.text().startswith(('Paragraph needs', 'PyTorch',
                                            'Parapred', 'ANARCI',
                                            'HMMER')):
            self.status.setText('')

    def showEvent(self, event):
        self.refresh()
        super().showEvent(event)

    # running

    def run(self):
        method = self.method.currentData()
        self.status.setText('Predicting...')
        QtWidgets.QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            if self.source.checkedId() == 0:
                name = self.object.currentText()
                if not name:
                    raise ValueError('Open a structure')
                self.residues = paratope.predict(
                    name, method, self.chain.currentData() or '',
                    _self=self.cmd)
                if self.colour.isChecked() or self.select.isChecked():
                    paratope.apply_scores(
                        name, self.residues, self.colour.isChecked(),
                        self.select.isChecked(), self.colour.isChecked(),
                        self.cmd)
            else:
                text = ''.join(c for c in self.sequence.toPlainText()
                               if c.isalpha())
                if not text:
                    raise ValueError('Paste an antibody sequence')
                self.residues = paratope.predict_sequence(text)
            if not self.residues:
                raise ValueError('No antibody variable domain found')
        except Exception as e:
            self.residues = []
            self.table.setRowCount(0)
            self.status.setText(str(e))
            return
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self._fill_table()

    def _fill_table(self, *_):
        threshold = self.threshold.value()
        rows = sorted(self.residues, key=lambda r: -r.score)
        self.table.setRowCount(len(rows))
        for r, residue in enumerate(rows):
            values = (residue.chain,
                      '%s%s' % (residue.letter, residue.resi),
                      '%.2f' % residue.score,
                      'yes' if residue.score >= threshold else '',
                      residue.region)
            for c, text in enumerate(values):
                self.table.setItem(r, c, QtWidgets.QTableWidgetItem(text))
        self.table.resizeColumnsToContents()
        if self.residues:
            binders = sum(1 for x in self.residues if x.score >= threshold)
            self.status.setText(
                '%s: %d residue(s) scored, %d at or above %.2f. A residue '
                'with no score is one the predictor did not look at.'
                % (paratope.method_label(self.method.currentData()),
                   len(self.residues), binders, threshold))

    def copy(self):
        lines = ['\t'.join(COLUMNS)]
        for r in range(self.table.rowCount()):
            lines.append('\t'.join(
                self.table.item(r, c).text()
                for c in range(self.table.columnCount())))
        QtWidgets.QApplication.clipboard().setText('\n'.join(lines))
        self.status.setText('Table copied')


def _console_font():
    from pymolx.gui import theme
    return theme.console_font()


def show_dialog():
    global _dialog
    if _dialog is None:
        from pymol.gui import get_qtwindow
        _dialog = ParatopeDialog(get_qtwindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
