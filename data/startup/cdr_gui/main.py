'''
CDR Loops dialog: a structure's chains or a pasted sequence numbered by
ANARCI, its framework regions and CDR loops under the chosen scheme, as a
table, as selections and in colour (pymolx.cdr)
'''

from pymol import cmd
from pymol.Qt import QtCore, QtWidgets

from pymolx import cdr

Qt = QtCore.Qt

_dialog = None

COLUMNS = ['Chain', 'Type', 'Loop', 'Sequence', 'Length', 'Residues']


class CDRDialog(QtWidgets.QDialog):

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent)
        self.cmd = _self
        self.found = []
        self.setObjectName('cdr_dialog')
        self.setWindowTitle('CDR Loops')
        self.resize(900, 620)

        layout = QtWidgets.QHBoxLayout(self)
        left = QtWidgets.QVBoxLayout()
        layout.addLayout(left, 2)
        right = QtWidgets.QVBoxLayout()
        layout.addLayout(right, 3)

        # what to number
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
        self.object.setObjectName('cdr_object')
        self.object.currentIndexChanged.connect(self._object_changed)
        grid.addWidget(self.object, 1, 1)
        grid.addWidget(QtWidgets.QLabel('Chain', box), 2, 0)
        self.chain = QtWidgets.QComboBox(box)
        self.chain.setObjectName('cdr_chain')
        grid.addWidget(self.chain, 2, 1)
        self.sequence = QtWidgets.QPlainTextEdit(box)
        self.sequence.setObjectName('cdr_sequence')
        self.sequence.setPlaceholderText(
            'one or more sequences, e.g.\n>Chain H\nEVQLVESGGG...')
        self.sequence.setFont(_console_font())
        self.sequence.setMaximumHeight(120)
        grid.addWidget(self.sequence, 3, 0, 1, 2)
        grid.addWidget(QtWidgets.QLabel('Numbering', box), 4, 0)
        self.scheme = QtWidgets.QComboBox(box)
        self.scheme.setObjectName('cdr_scheme')
        for label, value in cdr.SCHEMES:
            self.scheme.addItem(label, value)
        self.scheme.setToolTip(
            'The numbering, which also decides where each loop begins and '
            'ends: the schemes disagree, most of all about CDR1')
        grid.addWidget(self.scheme, 4, 1)
        grid.setColumnStretch(1, 1)
        left.addWidget(box)

        box = QtWidgets.QGroupBox('Show', self)
        box.setObjectName('prep_step')
        column = QtWidgets.QVBoxLayout(box)
        self.select = QtWidgets.QCheckBox('Selections per loop', box)
        self.select.setChecked(True)
        self.select.setToolTip('e.g. "igg_H_CDR3", and "igg_cdrs" for all '
                               'of them')
        self.colour = QtWidgets.QCheckBox('Colour the loops', box)
        self.colour.setChecked(True)
        self.colour.setToolTip('CDR1 orange, CDR2 green, CDR3 blue; the '
                               'variable domains\' framework grey')
        column.addWidget(self.select)
        column.addWidget(self.colour)
        left.addWidget(box)
        left.addStretch(1)

        self.run_button = QtWidgets.QPushButton('Find CDR Loops', self)
        self.run_button.setDefault(True)
        self.run_button.clicked.connect(self.run)
        left.addWidget(self.run_button)

        self.table = QtWidgets.QTableWidget(0, len(COLUMNS), self)
        self.table.setObjectName('cdr_table')
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
        missing = cdr.available()
        if missing:
            self.status.setText(missing)
            self.run_button.setEnabled(False)

    def _object_changed(self, *_):
        name = self.object.currentText()
        chains = self.cmd.get_chains('%' + name) if name else []
        self.chain.clear()
        self.chain.addItem('All chains', '')
        for one in chains or []:
            self.chain.addItem(one, one)

    def _source_changed(self, *_):
        structure = self.source.checkedId() == 0
        self.object.setEnabled(structure)
        self.chain.setEnabled(structure)
        self.sequence.setEnabled(not structure)

    def showEvent(self, event):
        self.refresh()
        super().showEvent(event)

    # running

    def run(self):
        self.status.setText('')
        try:
            if self.source.checkedId() == 0:
                self._from_structure()
            else:
                self._from_sequence()
        except Exception as e:
            self.status.setText(str(e))
            self.table.setRowCount(0)
            return
        self._fill_table()

    def _from_structure(self):
        name = self.object.currentText()
        if not name:
            raise ValueError('Open a structure, or give a sequence')
        self.found = cdr.find_cdrs(
            name, self.chain.currentData() or '',
            self.scheme.currentData(), select=self.select.isChecked(),
            color=self.colour.isChecked(), quiet=1, _self=self.cmd)
        if not self.found:
            raise ValueError('No antibody variable domain found in "%s"'
                             % name)

    def _from_sequence(self):
        text = self.sequence.toPlainText()
        entries = _read_sequences(text)
        if not entries:
            raise ValueError('Paste one or more sequences')
        self.found = []
        for title, sequence in entries:
            for domain in cdr.domains(sequence, self.scheme.currentData()):
                self.found.append((title, domain, []))
        if not self.found:
            raise ValueError('No antibody variable domain in what you '
                             'pasted')

    def _fill_table(self):
        rows = []
        for chain_id, domain, residues in self.found:
            by_name = dict(cdr.region_residues(domain, residues)) \
                if residues else {}
            for name, positions in domain.regions():
                if name not in cdr.CDR_NAMES:
                    continue
                sequence = ''.join(letter for _, letter in positions)
                where = by_name.get(name) or []
                rows.append((
                    chain_id,
                    '%s/%s' % (domain.chain_type, domain.species or '?'),
                    name, sequence, str(len(sequence)),
                    '%s-%s' % (where[0][1], where[-1][1]) if where else
                    '%s-%s' % (domain.label(positions[0][0]),
                               domain.label(positions[-1][0]))
                    if positions else '-'))
        self.table.setRowCount(len(rows))
        for r, row in enumerate(rows):
            for c, text in enumerate(row):
                item = QtWidgets.QTableWidgetItem(text)
                if c == 3:
                    item.setFont(_console_font())
                self.table.setItem(r, c, item)
        self.table.resizeColumnsToContents()
        kinds = {d.chain_type for _, d, _ in self.found}
        self.status.setText(
            '%d variable domain(s), %s numbering%s' % (
                len(self.found),
                cdr.scheme_label(self.scheme.currentData()),
                '; chain types ' + ', '.join(sorted(kinds)) if kinds
                else ''))

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


def _read_sequences(text):
    '''
    FASTA, or one bare sequence

    :return: [(title, sequence)]
    '''
    entries = []
    title, letters = '', []
    for line in str(text).splitlines():
        line = line.strip()
        if line.startswith('>'):
            if letters:
                entries.append((title or 'sequence', ''.join(letters)))
            title, letters = line[1:].strip(), []
        elif line:
            letters.append(''.join(c for c in line if c.isalpha()))
    if letters:
        entries.append((title or 'sequence', ''.join(letters)))
    return [(t, s) for t, s in entries if s]


def show_dialog():
    global _dialog
    if _dialog is None:
        from pymol.gui import get_qtwindow
        _dialog = CDRDialog(get_qtwindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
