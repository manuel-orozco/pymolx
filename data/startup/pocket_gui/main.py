'''
Binding Pocket dialog: choose the ligand (from the object's small
molecules and ions, or any selection) and the maximum distance; the
pocket's residues become a selection, a residue list to copy and a table
with each residue's closest contact.
'''

from pymol import cmd
from pymol.Qt import QtCore, QtGui, QtWidgets

from pymolx import pocket

Qt = QtCore.Qt

_dialog = None

NAMES = [('1-letter code (D52)', 'one'), ('3-letter code (ASP52)', 'three'),
         ('Number only (52)', 'none')]
COLUMNS = ['Chain', 'Residue', 'Number', 'Closest (Å)', 'Atom',
           'Ligand atom', 'Polar']


class PocketDialog(QtWidgets.QDialog):

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent)
        self.cmd = _self
        self.result = None
        self.setObjectName('pocket_dialog')
        self.setWindowTitle('Binding Pocket')
        self.resize(860, 600)

        layout = QtWidgets.QHBoxLayout(self)
        left = QtWidgets.QVBoxLayout()
        layout.addLayout(left, 2)
        right = QtWidgets.QVBoxLayout()
        layout.addLayout(right, 3)

        # ligand
        box = QtWidgets.QGroupBox('Ligand', self)
        box.setObjectName('prep_step')
        column = QtWidgets.QVBoxLayout(box)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel('Object', box))
        self.object = QtWidgets.QComboBox(box)
        self.object.currentIndexChanged.connect(self._fill_ligands)
        row.addWidget(self.object, 1)
        refresh = QtWidgets.QPushButton('Refresh', box)
        refresh.clicked.connect(self.refresh)
        row.addWidget(refresh)
        column.addLayout(row)
        self.ligands = QtWidgets.QListWidget(box)
        self.ligands.setObjectName('pocket_ligands')
        self.ligands.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.ligands.setToolTip('Small molecules and ions of the object; '
                                'select several for a pocket around all')
        self.ligands.itemDoubleClicked.connect(lambda _: self.run())
        column.addWidget(self.ligands, 1)
        self.custom = QtWidgets.QLineEdit(box)
        self.custom.setPlaceholderText('or a selection, e.g. resn ATP')
        self.custom.setToolTip('Any selection; when given, it is the '
                               'ligand instead of the list')
        self.custom.returnPressed.connect(self.run)
        column.addWidget(self.custom)
        left.addWidget(box, 1)

        # settings
        box = QtWidgets.QGroupBox('Settings', self)
        box.setObjectName('prep_step')
        grid = QtWidgets.QGridLayout(box)
        self.cutoff = QtWidgets.QDoubleSpinBox(box)
        self.cutoff.setRange(1.0, 15.0)
        self.cutoff.setDecimals(1)
        self.cutoff.setSingleStep(0.5)
        self.cutoff.setValue(4.0)
        self.cutoff.setSuffix(' Å')
        self.cutoff.setToolTip('Maximum distance between a ligand atom and '
                               'an atom of the residue')
        self.names = QtWidgets.QComboBox(box)
        for label, value in NAMES:
            self.names.addItem(label, value)
        self.names.currentIndexChanged.connect(self._fill_text)
        self.name = QtWidgets.QLineEdit('pocket', box)
        self.name.setToolTip('Name of the selection made')
        rows = [('Maximum distance', self.cutoff), ('Residue names',
                                                     self.names),
                ('Selection', self.name)]
        for k, (label, widget) in enumerate(rows):
            grid.addWidget(QtWidgets.QLabel(label, box), k, 0)
            grid.addWidget(widget, k, 1)
        self.hydrogens = QtWidgets.QCheckBox('Hydrogens count too', box)
        self.hydrogens.setToolTip('Off: distances between heavy atoms, so '
                                  'the pocket doesn\'t depend on whether '
                                  'hydrogens were added')
        grid.addWidget(self.hydrogens, len(rows), 0, 1, 2)
        self.waters = QtWidgets.QCheckBox('Include waters', box)
        self.waters.setToolTip('Waters near the ligand (e.g. bridging '
                               'waters) are pocket residues too')
        grid.addWidget(self.waters, len(rows) + 1, 0, 1, 2)
        left.addWidget(box)

        run = QtWidgets.QHBoxLayout()
        self.status = QtWidgets.QLabel(self)
        self.status.setObjectName('prep_note')
        self.status.setWordWrap(True)
        run.addWidget(self.status, 1)
        self.run_button = QtWidgets.QPushButton('Find Pocket', self)
        self.run_button.setDefault(True)
        self.run_button.clicked.connect(self.run)
        run.addWidget(self.run_button)
        left.addLayout(run)

        # results
        right.addWidget(QtWidgets.QLabel('Pocket residues', self))
        row = QtWidgets.QHBoxLayout()
        self.text = QtWidgets.QPlainTextEdit(self)
        self.text.setObjectName('pocket_text')
        self.text.setReadOnly(True)
        self.text.setMaximumHeight(110)
        row.addWidget(self.text, 1)
        self.copy_button = QtWidgets.QPushButton('Copy', self)
        self.copy_button.clicked.connect(
            lambda: QtWidgets.QApplication.clipboard().setText(
                self.text.toPlainText()))
        row.addWidget(self.copy_button, 0, Qt.AlignmentFlag.AlignTop)
        right.addLayout(row)
        self.table = QtWidgets.QTableWidget(0, len(COLUMNS), self)
        self.table.setObjectName('pocket_table')
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setToolTip('Double-click a residue to zoom on it')
        self.table.cellDoubleClicked.connect(self._zoom_residue)
        right.addWidget(self.table, 1)
        note = QtWidgets.QLabel('Polar: an N/O atom of the residue within '
                                '%.1f Å of an N/O of the ligand (likely '
                                'hydrogen bond or salt bridge)' %
                                pocket.POLAR_CUTOFF, self)
        note.setObjectName('prep_note')
        note.setWordWrap(True)
        right.addWidget(note)
        actions = QtWidgets.QHBoxLayout()
        self.show_button = QtWidgets.QPushButton('Show', self)
        self.show_button.setToolTip('Ligand and pocket side chains as '
                                    'sticks, zoomed')
        self.show_button.clicked.connect(self.show_pocket)
        actions.addWidget(self.show_button)
        actions.addStretch(1)
        self.save_button = QtWidgets.QPushButton('Save CSV...', self)
        self.save_button.clicked.connect(self.save_csv)
        actions.addWidget(self.save_button)
        right.addLayout(actions)
        self._enable_results(False)

    # ligand

    def refresh(self):
        current = self.object.currentText()
        names = [name for name in self.cmd.get_object_list('all') or []
                 if self.cmd.get_type(name) == 'object:molecule' and
                 not name.startswith('_')]
        self.object.blockSignals(True)
        self.object.clear()
        self.object.addItems(names)
        if current in names:
            self.object.setCurrentText(current)
        else:
            # the first object with a ligand
            for k, name in enumerate(names):
                if self.cmd.count_atoms('%%%s and not polymer and not '
                                        'solvent' % name):
                    self.object.setCurrentIndex(k)
                    break
        self.object.blockSignals(False)
        self._fill_ligands()

    def _fill_ligands(self, *_):
        selected = {item.data(Qt.ItemDataRole.UserRole)
                    for item in self.ligands.selectedItems()}
        self.ligands.clear()
        obj = self.object.currentText()
        if not obj:
            return
        found = pocket.ligands(obj, _self=self.cmd)
        for ligand in found:
            ion = ligand['atoms'] == 1
            label = '%s %s  (chain %s, %s)' % (
                ligand['resn'], ligand['resi'], ligand['chain'] or '-',
                'ion' if ion else '%d atoms' % ligand['atoms'])
            item = QtWidgets.QListWidgetItem(label, self.ligands)
            item.setData(Qt.ItemDataRole.UserRole, ligand['selection'])
            if ligand['selection'] in selected:
                item.setSelected(True)
        if not self.ligands.selectedItems() and self.ligands.count():
            # the largest: usually the drug or cofactor, not a buffer
            sizes = [x['atoms'] for x in found]
            self.ligands.item(sizes.index(max(sizes))).setSelected(True)

    def ligand_selection(self):
        custom = self.custom.text().strip()
        if custom:
            return custom
        parts = ['(%s)' % item.data(Qt.ItemDataRole.UserRole)
                 for item in self.ligands.selectedItems()]
        return ' or '.join(parts)

    def showEvent(self, event):
        self.refresh()
        super().showEvent(event)

    # running

    def run(self):
        ligand = self.ligand_selection()
        if not ligand:
            self.status.setText('Choose a ligand or type a selection')
            return
        try:
            self.result = pocket.find(
                ligand, self.cutoff.value(),
                hydrogens=self.hydrogens.isChecked(),
                waters=self.waters.isChecked(),
                name=self.name.text().strip() or 'pocket', _self=self.cmd)
        except Exception as e:
            self.status.setText('Failed: %s' % e)
            return
        r = self.result
        polar = sum(x['polar'] for x in r.residues)
        self.status.setText('%d residues within %g Å, %d with polar '
                            'contacts; selection "%s"' % (
                                len(r.residues), r.cutoff, polar,
                                r.selection))
        self._fill_text()
        self._fill_table()
        self._enable_results(bool(r.residues))
        print(' Binding pocket of %s (%g A):\n %s' % (
            ligand, r.cutoff, self.text.toPlainText().replace('\n', '\n ')))

    def _fill_text(self, *_):
        if self.result is None:
            return
        self.text.setPlainText(self.result.text(self.names.currentData()) or
                               '(no residues)')

    def _fill_table(self):
        residues = self.result.residues
        self.table.setRowCount(len(residues))
        for row, r in enumerate(residues):
            values = [r['chain'] or '-', r['resn'], r['resi'],
                      '%.2f' % r['distance'], r['atom'], r['ligand_atom'],
                      'yes' if r['polar'] else '']
            for col, value in enumerate(values):
                item = QtWidgets.QTableWidgetItem(value)
                if col == 3:
                    item.setTextAlignment(int(
                        Qt.AlignmentFlag.AlignRight |
                        Qt.AlignmentFlag.AlignVCenter))
                self.table.setItem(row, col, item)
        self.table.resizeColumnsToContents()

    def _enable_results(self, enabled):
        for widget in (self.copy_button, self.show_button, self.save_button):
            widget.setEnabled(enabled)

    def _zoom_residue(self, row, _column):
        if not self.result or row >= len(self.result.residues):
            return
        r = self.result.residues[row]
        sele = '%%%s and chain "%s" and resi %s' % (
            r['model'], r['chain'], pocket._resi(r['resi']))
        self.cmd.show('sticks', '(%s) and not hydro' % sele)
        self.cmd.zoom('(%s) or (%s)' % (sele, self.result.ligand), 2,
                      animate=0.5)

    def show_pocket(self):
        if self.result:
            pocket.show(self.result, _self=self.cmd)

    def save_csv(self):
        if not self.result:
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save pocket residues', 'pocket.csv', 'CSV (*.csv)')
        if path:
            self.result.save_csv(path)
            self.status.setText('Saved %s' % path)


def show_dialog():
    global _dialog
    if _dialog is None:
        from pymol.gui import get_qtwindow
        _dialog = PocketDialog(get_qtwindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
