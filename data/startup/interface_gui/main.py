'''
Interface Analysis dialog: PISA-like interface list, residues and bonds,
and the epitope: the residues of chosen target chains buried by the other
chains (and the partners' residues burying them).
'''

import csv

from pymol import cmd
from pymol.Qt import QtCore, QtGui, QtWidgets

from pymolx import interfaces

Qt = QtCore.Qt

_dialog = None

INTERFACE_COLUMNS = ['#', 'Chains', 'Interface area (Å²)',
                     'ΔiG (kcal/mol)', 'H-bonds', 'Salt bridges',
                     'Disulfides', 'Residues']
RESIDUE_COLUMNS = ['Chain', 'Residue', 'Number', 'ASA (Å²)',
                   'BSA (Å²)', 'ΔiG (kcal/mol)', 'Bonds']
BOND_COLUMNS = ['Type', 'Atom 1', 'Atom 2', 'Distance (Å)']


class NumericItem(QtWidgets.QTableWidgetItem):
    '''
    Sorts by its number (UserRole), not by the displayed text.
    '''

    def __lt__(self, other):
        mine = self.data(Qt.ItemDataRole.UserRole)
        theirs = other.data(Qt.ItemDataRole.UserRole)
        if mine is None or theirs is None:
            return super().__lt__(other)
        return mine < theirs


def _item(value, fmt=None):
    '''
    Table item; numbers (with fmt) sort numerically.
    '''
    item = NumericItem() if fmt is not None else QtWidgets.QTableWidgetItem()
    if fmt is not None:
        item.setData(Qt.ItemDataRole.DisplayRole, fmt % value)
        item.setData(Qt.ItemDataRole.UserRole, float(value))
        item.setTextAlignment(int(Qt.AlignmentFlag.AlignRight |
                                  Qt.AlignmentFlag.AlignVCenter))
    else:
        item.setData(Qt.ItemDataRole.DisplayRole, str(value))
    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
    return item


def _table(columns, parent):
    table = QtWidgets.QTableWidget(0, len(columns), parent)
    table.setHorizontalHeaderLabels(columns)
    table.setSelectionBehavior(
        QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
    table.setSelectionMode(
        QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
    table.verticalHeader().hide()
    table.horizontalHeader().setStretchLastSection(True)
    table.setAlternatingRowColors(True)
    return table


def _fill(table, rows):
    '''
    :param rows: lists of QTableWidgetItem
    '''
    table.setSortingEnabled(False)
    table.setRowCount(len(rows))
    for r, items in enumerate(rows):
        for c, item in enumerate(items):
            table.setItem(r, c, item)
    table.resizeColumnsToContents()
    # keep the filled order until the user clicks a header
    table.horizontalHeader().setSortIndicator(-1, Qt.SortOrder.AscendingOrder)
    table.setSortingEnabled(True)


class InterfaceDialog(QtWidgets.QDialog):

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent)
        self.cmd = _self
        self.interfaces = []
        self.analyzed_object = ''
        self.setObjectName('interface_dialog')
        self.setWindowTitle('Interface Analysis (PISA-like)')
        self.resize(820, 640)

        layout = QtWidgets.QVBoxLayout(self)

        top = QtWidgets.QHBoxLayout()
        top.addWidget(QtWidgets.QLabel('Object', self))
        self.object_combo = QtWidgets.QComboBox(self)
        self.object_combo.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed)
        top.addWidget(self.object_combo, 1)
        top.addWidget(QtWidgets.QLabel('State', self))
        self.state = QtWidgets.QSpinBox(self)
        self.state.setRange(-1, 99999)
        self.state.setValue(-1)
        self.state.setToolTip('-1 = current state')
        top.addWidget(self.state)
        self.analyze_button = QtWidgets.QPushButton('Analyze', self)
        self.analyze_button.setDefault(True)
        self.analyze_button.clicked.connect(self.analyze)
        top.addWidget(self.analyze_button)
        layout.addLayout(top)

        self.interface_table = _table(INTERFACE_COLUMNS, self)
        self.interface_table.itemSelectionChanged.connect(
            self.interface_selected)
        layout.addWidget(self.interface_table, 2)

        self.tabs = QtWidgets.QTabWidget(self)
        self.residue_table = _table(RESIDUE_COLUMNS, self)
        self.residue_table.itemDoubleClicked.connect(self.zoom_residue)
        self.bond_table = _table(BOND_COLUMNS, self)
        self.bond_table.itemDoubleClicked.connect(self.zoom_bond)
        self.tabs.addTab(self.residue_table, 'Residues')
        self.tabs.addTab(self.bond_table, 'Bonds')
        self.tabs.addTab(self._epitope_tab(), 'Epitope')
        layout.addWidget(self.tabs, 3)

        buttons = QtWidgets.QHBoxLayout()
        self.status = QtWidgets.QLabel(self)
        self.status.setObjectName('interface_status')
        buttons.addWidget(self.status, 1)
        self.show_button = QtWidgets.QPushButton('Show in viewer', self)
        self.show_button.clicked.connect(self.show_selected)
        buttons.addWidget(self.show_button)
        self.export_button = QtWidgets.QPushButton('Export CSV...', self)
        self.export_button.clicked.connect(self.export_csv)
        buttons.addWidget(self.export_button)
        layout.addLayout(buttons)

        note = QtWidgets.QLabel(interfaces.METHOD_NOTE, self)
        note.setWordWrap(True)
        note.setObjectName('interface_note')
        layout.addWidget(note)

        self._update_buttons()

    def _epitope_tab(self):
        tab = QtWidgets.QWidget(self)
        layout = QtWidgets.QVBoxLayout(tab)

        top = QtWidgets.QHBoxLayout()
        top.addWidget(QtWidgets.QLabel('Target chain(s)', tab))
        self.target_list = QtWidgets.QListWidget(tab)
        self.target_list.setObjectName('epitope_targets')
        self.target_list.setFlow(QtWidgets.QListView.Flow.LeftToRight)
        self.target_list.setWrapping(True)
        self.target_list.setMaximumHeight(52)
        self.target_list.setToolTip(
            'The chain(s) whose buried residues form the epitope (by '
            'default the chains that aren\'t antibody chains)')
        self.target_list.itemChanged.connect(self.update_epitope)
        top.addWidget(self.target_list, 1)
        top.addWidget(QtWidgets.QLabel('Min. buried area', tab))
        self.min_bsa = QtWidgets.QDoubleSpinBox(tab)
        self.min_bsa.setRange(0.0, 200.0)
        self.min_bsa.setDecimals(1)
        self.min_bsa.setSingleStep(1.0)
        self.min_bsa.setSuffix(' Å²')
        self.min_bsa.setToolTip('Residues need more buried area than this '
                                '(0: any, like PISA\'s interface residues)')
        self.min_bsa.valueChanged.connect(self.update_epitope)
        top.addWidget(self.min_bsa)
        top.addWidget(QtWidgets.QLabel('Names', tab))
        self.name_style = QtWidgets.QComboBox(tab)
        for label, value in (('1-letter (D52)', 'one'),
                             ('3-letter (ASP52)', 'three'),
                             ('Numbers only (52)', 'none')):
            self.name_style.addItem(label, value)
        self.name_style.currentIndexChanged.connect(self.update_epitope)
        top.addWidget(self.name_style)
        layout.addLayout(top)

        self.epitope_label = QtWidgets.QLabel(
            'Epitope (target residues buried by the other chains)', tab)
        layout.addWidget(self.epitope_label)
        self.epitope_text, row = self._residue_box(tab)
        layout.addLayout(row, 1)
        self.paratope_label = QtWidgets.QLabel(
            'Partner residues (paratope)', tab)
        layout.addWidget(self.paratope_label)
        self.paratope_text, row = self._residue_box(tab)
        layout.addLayout(row, 1)

        bottom = QtWidgets.QHBoxLayout()
        self.epitope_summary = QtWidgets.QLabel(tab)
        self.epitope_summary.setObjectName('interface_note')
        bottom.addWidget(self.epitope_summary, 1)
        self.select_epitope_button = QtWidgets.QPushButton(
            'Select in Viewer', tab)
        self.select_epitope_button.setToolTip(
            'Selections "epitope" and "paratope"')
        self.select_epitope_button.clicked.connect(self.select_epitope)
        bottom.addWidget(self.select_epitope_button)
        layout.addLayout(bottom)
        self.epitope = {}
        self.paratope = {}
        return tab

    def _residue_box(self, parent):
        row = QtWidgets.QHBoxLayout()
        text = QtWidgets.QPlainTextEdit(parent)
        text.setReadOnly(True)
        text.setObjectName('epitope_text')
        row.addWidget(text, 1)
        copy = QtWidgets.QPushButton('Copy', parent)
        copy.clicked.connect(lambda: QtWidgets.QApplication.clipboard()
                             .setText(text.toPlainText()))
        row.addWidget(copy, 0, Qt.AlignmentFlag.AlignTop)
        return text, row

    def targets(self):
        return [self.target_list.item(i).text()
                for i in range(self.target_list.count())
                if self.target_list.item(i).checkState() ==
                Qt.CheckState.Checked]

    def _fill_targets(self):
        obj = self.analyzed_object
        chains = sorted(self.cmd.get_chains('%%%s and polymer' % obj))
        targets = interfaces.default_targets(obj, _self=self.cmd)
        self.target_list.blockSignals(True)
        self.target_list.clear()
        for chain in chains:
            item = QtWidgets.QListWidgetItem(chain, self.target_list)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked if chain in targets
                               else Qt.CheckState.Unchecked)
        self.target_list.blockSignals(False)

    def update_epitope(self, *_):
        targets = self.targets()
        if not targets:
            self.epitope, self.paratope = {}, {}
            self.epitope_text.setPlainText('')
            self.paratope_text.setPlainText('')
            self.epitope_summary.setText('Choose the target chain(s)')
            self.select_epitope_button.setEnabled(False)
            return
        min_bsa = self.min_bsa.value()
        self.epitope = interfaces.contact_residues(
            self.interfaces, targets, min_bsa)
        self.paratope = interfaces.contact_residues(
            self.interfaces, targets, min_bsa, side='partner')
        names = self.name_style.currentData()
        self.epitope_text.setPlainText(
            interfaces.format_epitope(self.epitope, names))
        self.paratope_text.setPlainText(
            interfaces.format_epitope(self.paratope, names))
        count = sum(len(r) for r in self.epitope.values())
        area = sum(b for r in self.epitope.values() for _, _, b in r)
        self.epitope_summary.setText(
            '%d epitope residues on %s, %.0f Å² buried; %d partner '
            'residues' % (count, '+'.join(targets), area,
                          sum(len(r) for r in self.paratope.values()))
            if count else 'No residues of %s buried by the other chains' %
            '+'.join(targets))
        self.select_epitope_button.setEnabled(bool(count))

    def select_epitope(self):
        obj = self.analyzed_object
        for name, by_chain in (('epitope', self.epitope),
                               ('paratope', self.paratope)):
            parts = ['(chain "%s" and resi %s)' % (
                chain, '+'.join(interfaces._resi(resi)
                                for resi, _, _ in residues))
                for chain, residues in by_chain.items()]
            self.cmd.select(name, '%%%s and (%s)' % (
                obj, ' or '.join(parts) or 'none'), enable=0)
        self.cmd.enable('epitope')
        self.cmd.zoom('epitope', 6, animate=1)
        self.status.setText('Selected "epitope" (%d atoms) and "paratope"' %
                            self.cmd.count_atoms('epitope'))

    # session

    def refresh_objects(self):
        current = self.object_combo.currentText()
        names = [name for name in self.cmd.get_object_list()
                 if len(self.cmd.get_chains('%%%s and polymer' % name)) > 1]
        self.object_combo.clear()
        self.object_combo.addItems(names)
        if current in names:
            self.object_combo.setCurrentText(current)

    def showEvent(self, event):
        self.refresh_objects()
        super().showEvent(event)

    # analysis

    def analyze(self):
        obj = self.object_combo.currentText()
        if not obj:
            self.status.setText('No object with two or more chains')
            return
        QtWidgets.QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            self.interfaces = interfaces.analyze(
                obj, self.state.value(), _self=self.cmd)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self.analyzed_object = obj
        print(interfaces.format_table(obj, self.interfaces))

        rows = []
        for n, iface in enumerate(self.interfaces, 1):
            hb, sb, ds = iface.counts()
            rows.append([
                _item(n, '%d'), _item('%s – %s' % iface.chains),
                _item(iface.area, '%.1f'), _item(iface.dg, '%.1f'),
                _item(hb, '%d'), _item(sb, '%d'), _item(ds, '%d'),
                _item('%d + %d' % (iface.residue_count(0),
                                   iface.residue_count(1))),
            ])
            rows[-1][0].setData(Qt.ItemDataRole.UserRole + 1, n - 1)
        _fill(self.interface_table, rows)
        self.status.setText('%d interfaces in %s' % (len(self.interfaces), obj))
        self._fill_targets()
        self.update_epitope()
        if self.interfaces:
            self.interface_table.selectRow(0)
        else:
            _fill(self.residue_table, [])
            _fill(self.bond_table, [])
        self._update_buttons()

    def selected_interface(self):
        rows = self.interface_table.selectionModel().selectedRows()
        if not rows:
            return None
        item = self.interface_table.item(rows[0].row(), 0)
        return self.interfaces[item.data(Qt.ItemDataRole.UserRole + 1)]

    def interface_selected(self):
        iface = self.selected_interface()
        self._update_buttons()
        if iface is None:
            return
        residues = sorted(iface.residues,
                          key=lambda r: (iface.chains.index(r['chain']),
                                         _resi_key(r['resi'])))
        _fill(self.residue_table, [[
            _item(r['chain']), _item(r['resn']), _item(r['resi']),
            _item(r['asa'], '%.1f'), _item(r['bsa'], '%.1f'),
            _item(r['dg'], '%.2f'), _item(r['flags']),
        ] for r in residues])
        bonds = []
        for kind, items in (('H-bond', iface.hbonds),
                            ('Salt bridge', iface.salt_bridges),
                            ('Disulfide', iface.disulfides)):
            for a, b, d in items:
                row = [_item(kind), _item(a.label()), _item(b.label()),
                       _item(d, '%.2f')]
                row[0].setData(Qt.ItemDataRole.UserRole + 1, (a.index, b.index))
                bonds.append(row)
        _fill(self.bond_table, bonds)

    def _update_buttons(self):
        selected = self.selected_interface() is not None
        self.show_button.setEnabled(selected)
        self.export_button.setEnabled(bool(self.interfaces))

    # viewer

    def show_selected(self):
        iface = self.selected_interface()
        if iface is not None:
            group = interfaces.show_interface(iface, _self=self.cmd)
            self.status.setText('Shown as group "%s"' % group)

    def zoom_residue(self, item):
        row = item.row()
        chain = self.residue_table.item(row, 0).text()
        resi = self.residue_table.item(row, 2).text()
        self.cmd.zoom('%%%s and chain "%s" and resi %s' % (
            self.analyzed_object, chain, interfaces._resi(resi)), 6,
            animate=1)

    def zoom_bond(self, item):
        a, b = self.bond_table.item(item.row(), 0).data(
            Qt.ItemDataRole.UserRole + 1)
        self.cmd.zoom('%%%s and index %d+%d' % (self.analyzed_object, a, b),
                      6, animate=1)

    # export

    def export_csv(self):
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Export interfaces',
            'interfaces_%s.csv' % self.analyzed_object, 'CSV (*.csv)')
        if filename:
            write_csv(filename, self.analyzed_object, self.interfaces)
            self.status.setText('Exported to %s' % filename)


def _resi_key(resi):
    digits = ''.join(c for c in resi if c.isdigit() or c == '-')
    try:
        return (int(digits), resi)
    except ValueError:
        return (0, resi)


def write_csv(filename, obj, interface_list):
    '''
    One row per interface residue, with the interface summary repeated.
    '''
    with open(filename, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['object', 'interface', 'chains', 'area', 'dG',
                         'hbonds', 'salt_bridges', 'disulfides', 'chain',
                         'resn', 'resi', 'asa', 'bsa', 'residue_dG', 'bonds'])
        for n, iface in enumerate(interface_list, 1):
            hb, sb, ds = iface.counts()
            summary = [obj, n, '%s-%s' % iface.chains, '%.1f' % iface.area,
                       '%.2f' % iface.dg, hb, sb, ds]
            for r in iface.residues:
                writer.writerow(summary + [
                    r['chain'], r['resn'], r['resi'], '%.1f' % r['asa'],
                    '%.1f' % r['bsa'], '%.3f' % r['dg'], r['flags']])


def show_dialog():
    global _dialog
    if _dialog is None:
        from pymol.gui import get_qtwindow
        _dialog = InterfaceDialog(get_qtwindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
