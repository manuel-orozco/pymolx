'''
Interface Analysis dialog: PISA-like interface list, residues and bonds.
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
