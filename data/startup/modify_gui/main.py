'''
Modify Sequences dialog: the chains of the loaded structures in a
table; for the chains selected there (all of their residues, or the
residue numbers typed below the table): renumber the residues, change
the chain ID or the segment ID, swap two chain IDs, merge the chains
into one chain, or make one new object of them.
'''

import fnmatch

from pymol import cmd
from pymol.Qt import QtCore, QtWidgets

from pymolx import modify, sequence, undo

Qt = QtCore.Qt

_dialog = None

COLUMNS = ['Object', 'Chain', 'Segment', 'Residues', 'Numbers', 'Other',
           'Sequence']
OBJECT, CHAIN, SEGMENT, RESIDUES, NUMBERS, OTHER, SEQUENCE = range(7)

NO_ID = '(none)'

# renumbering: label, (start or offset, keep gaps), what it does
MODES = [
    ('First residue becomes', 'start',
     'The first residue of each chain gets this number and the others '
     'follow at their distance: gaps in the numbering and insertion '
     'codes (28A) are kept'),
    ('Shift every number by', 'offset',
     'This is added to every residue number (negative: subtracted); '
     'gaps and insertion codes are kept'),
    ('Number one by one from', 'count',
     'The residues are numbered one by one from this number: no gaps, '
     'no insertion codes'),
]
NUMBERINGS = [
    ('Continue after the chain before', 'continue',
     'The first chain keeps its numbers; each later chain is shifted to '
     'start after the one before (its own gaps and insertion codes '
     'kept)'),
    ('Number one by one from 1', 'consecutive',
     'All residues are numbered anew: 1, 2, 3, ...'),
    ('Keep the numbers', 'keep',
     'No number changes; fails if two of the chains share a number'),
]


def _box(title, parent):
    box = QtWidgets.QGroupBox(title, parent)
    box.setObjectName('prep_step')
    return box


class ModifyDialog(QtWidgets.QDialog):

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent)
        self.cmd = _self
        self.rows = []
        self._counts = None
        self.setObjectName('modify_dialog')
        self.setWindowTitle('Modify Sequences')
        self.resize(1040, 640)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self._chains_box(), 1)
        row = QtWidgets.QHBoxLayout()
        for box in (self._renumber_box(), self._ids_box(),
                    self._merge_box(), self._combine_box()):
            row.addWidget(box, 1)
        layout.addLayout(row)

        row = QtWidgets.QHBoxLayout()
        self.status = QtWidgets.QLabel(
            'Select chains in the table, then change them below', self)
        self.status.setObjectName('prep_note')
        self.status.setWordWrap(True)
        row.addWidget(self.status, 1)
        self.undo_button = QtWidgets.QPushButton('Undo', self)
        self.undo_button.setToolTip('Take back the last change (as Undo in '
                                    'the toolbar)')
        self.undo_button.clicked.connect(self.undo)
        row.addWidget(self.undo_button, 0, Qt.AlignmentFlag.AlignBottom)
        layout.addLayout(row)

    # widgets

    def _chains_box(self):
        box = _box('Chains', self)
        column = QtWidgets.QVBoxLayout(box)
        row = QtWidgets.QHBoxLayout()
        self.filter = QtWidgets.QLineEdit(box)
        self.filter.setPlaceholderText('filter objects, e.g. *seed*')
        self.filter.textChanged.connect(self._filter)
        row.addWidget(self.filter, 1)
        for label, state in (('All', True), ('None', False)):
            button = QtWidgets.QPushButton(label, box)
            button.clicked.connect(
                lambda _=False, s=state: self._select_all(s))
            row.addWidget(button)
        refresh = QtWidgets.QPushButton('Refresh', box)
        refresh.clicked.connect(self.refresh)
        row.addWidget(refresh)
        column.addLayout(row)

        self.table = QtWidgets.QTableWidget(0, len(COLUMNS), box)
        self.table.setObjectName('modify_table')
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setToolTip('Select the chains to change (Ctrl or Shift '
                              'for several); double-click to zoom on one')
        self.table.itemSelectionChanged.connect(self._selection_changed)
        self.table.cellDoubleClicked.connect(self._zoom)
        column.addWidget(self.table, 1)

        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel('Only residues', box))
        self.residues = QtWidgets.QLineEdit(box)
        self.residues.setObjectName('modify_residues')
        self.residues.setPlaceholderText('all; or numbers, e.g. 100-150, 200')
        self.residues.setToolTip(
            'Empty: the whole chains selected above. Residue numbers or '
            'ranges limit every change to those residues of the chains: '
            'renumber a stretch, or give part of a chain another chain ID')
        row.addWidget(self.residues, 1)
        row.addStretch(2)
        column.addLayout(row)
        return box

    def _renumber_box(self):
        box = _box('Renumber residues', self)
        grid = QtWidgets.QGridLayout(box)
        self.mode = QtWidgets.QComboBox(box)
        for k, (label, value, tip) in enumerate(MODES):
            self.mode.addItem(label, value)
            self.mode.setItemData(k, tip, Qt.ItemDataRole.ToolTipRole)
        self.number = QtWidgets.QSpinBox(box)
        self.number.setRange(-99999, 99999)
        self.number.setValue(1)
        self.hetero = QtWidgets.QCheckBox('Also ligands and waters', box)
        self.hetero.setToolTip(
            'Off: only the protein or nucleic acid residues are '
            'renumbered; ligands, ions and waters with the same chain ID '
            'keep their numbers')
        grid.addWidget(self.mode, 0, 0, 1, 2)
        grid.addWidget(QtWidgets.QLabel('Number', box), 1, 0)
        grid.addWidget(self.number, 1, 1)
        grid.addWidget(self.hetero, 2, 0, 1, 2)
        grid.setRowStretch(3, 1)
        self.renumber_button = QtWidgets.QPushButton('Renumber', box)
        self.renumber_button.clicked.connect(self.renumber)
        grid.addWidget(self.renumber_button, 4, 0, 1, 2)
        return box

    def _ids_box(self):
        box = _box('Chain and segment IDs', self)
        grid = QtWidgets.QGridLayout(box)
        self.chain_check = QtWidgets.QCheckBox('Chain ID', box)
        self.chain_check.setChecked(True)
        self.chain_id = QtWidgets.QLineEdit(box)
        self.chain_id.setToolTip(
            'The new chain ID (a PDB file holds one character; mmCIF '
            'more). If another chain has it already, the residues join '
            'that chain, as long as no number is used twice')
        self.segi_check = QtWidgets.QCheckBox('Segment ID', box)
        self.segi_id = QtWidgets.QLineEdit(box)
        self.segi_id.setToolTip(
            'The new segment ID (a PDB file holds four characters); '
            'empty with the box checked: no segment ID. Atoms are '
            'ordered by segment first, so a chain should have one')
        self.chain_id.textEdited.connect(
            lambda _: self.chain_check.setChecked(True))
        self.segi_id.textEdited.connect(
            lambda _: self.segi_check.setChecked(True))
        grid.addWidget(self.chain_check, 0, 0)
        grid.addWidget(self.chain_id, 0, 1)
        grid.addWidget(self.segi_check, 1, 0)
        grid.addWidget(self.segi_id, 1, 1)
        grid.setRowStretch(2, 1)
        row = QtWidgets.QHBoxLayout()
        self.ids_button = QtWidgets.QPushButton('Change IDs', box)
        self.ids_button.clicked.connect(self.change_ids)
        row.addWidget(self.ids_button)
        self.swap_button = QtWidgets.QPushButton('Swap Two', box)
        self.swap_button.setToolTip('With two chains selected: each gets '
                                    'the other\'s chain ID')
        self.swap_button.clicked.connect(self.swap)
        row.addWidget(self.swap_button)
        grid.addLayout(row, 3, 0, 1, 2)
        return box

    def _merge_box(self):
        box = _box('Merge into one chain', self)
        grid = QtWidgets.QGridLayout(box)
        self.merge_chain = QtWidgets.QLineEdit(box)
        self.merge_chain.setPlaceholderText('the first chain\'s')
        self.merge_chain.setToolTip('Chain ID of the merged chain')
        self.order = QtWidgets.QLineEdit(box)
        self.order.setToolTip('Chain IDs in the order the merged chain '
                              'has them (filled in from the chains '
                              'selected; change it to e.g. "L H")')
        self.numbering = QtWidgets.QComboBox(box)
        for k, (label, value, tip) in enumerate(NUMBERINGS):
            self.numbering.addItem(label, value)
            self.numbering.setItemData(k, tip, Qt.ItemDataRole.ToolTipRole)
        self.numbering.currentIndexChanged.connect(self._numbering_changed)
        self.gap = QtWidgets.QSpinBox(box)
        self.gap.setRange(0, 9999)
        self.gap.setToolTip('Numbers left out between two chains, e.g. '
                            'for a linker to be built')
        rows = [('Chain ID', self.merge_chain), ('Order', self.order),
                ('Numbers', self.numbering), ('Gap', self.gap)]
        for k, (label, widget) in enumerate(rows):
            grid.addWidget(QtWidgets.QLabel(label, box), k, 0)
            grid.addWidget(widget, k, 1)
        grid.setRowStretch(len(rows), 1)
        self.merge_button = QtWidgets.QPushButton('Merge Chains', box)
        self.merge_button.setToolTip(
            'The selected chains of an object become one chain. Only '
            'the names change: no bond is made')
        self.merge_button.clicked.connect(self.merge)
        grid.addWidget(self.merge_button, len(rows) + 1, 0, 1, 2)
        return box

    def _combine_box(self):
        box = _box('Merge into one object', self)
        grid = QtWidgets.QGridLayout(box)
        self.new_name = QtWidgets.QLineEdit('merged', box)
        self.new_name.setToolTip('Name of the new object')
        self.delete_sources = QtWidgets.QCheckBox(
            'Remove from their objects', box)
        self.delete_sources.setToolTip(
            'On: the chains are moved (an object left empty is deleted). '
            'Off: they are copied')
        note = QtWidgets.QLabel(
            'Chains of several objects with the same ID get a free one, '
            'so they stay apart', box)
        note.setObjectName('prep_note')
        note.setWordWrap(True)
        grid.addWidget(QtWidgets.QLabel('Name', box), 0, 0)
        grid.addWidget(self.new_name, 0, 1)
        grid.addWidget(self.delete_sources, 1, 0, 1, 2)
        grid.addWidget(note, 2, 0, 1, 2)
        grid.setRowStretch(3, 1)
        self.combine_button = QtWidgets.QPushButton('Make Object', box)
        self.combine_button.clicked.connect(self.combine)
        grid.addWidget(self.combine_button, 4, 0, 1, 2)
        return box

    # the table

    def refresh(self, select=None):
        '''
        Fill the table from the session.

        :param select: [(object, segi, chain)] to select {default: what
            was selected}
        '''
        if select is None:
            select = [(r['object'], r['segi'], r['chain'])
                      for r in self.selected_rows()]
        select = set(select)
        self.rows = [row for row in modify.chain_table('all', self.cmd)
                     if not row['object'].startswith('_')]
        self._counts = sequence.get_change_counts(_self=self.cmd)
        table = self.table
        table.blockSignals(True)
        table.clearSelection()
        table.setRowCount(len(self.rows))
        right = int(Qt.AlignmentFlag.AlignRight |
                    Qt.AlignmentFlag.AlignVCenter)
        for k, row in enumerate(self.rows):
            numbers = row['first'] if row['first'] == row['last'] else \
                '%s to %s' % (row['first'], row['last'])
            if row['gaps']:
                numbers += ' (%d gap%s)' % (row['gaps'],
                                            '' if row['gaps'] == 1 else 's')
            values = [row['object'], row['chain'] or NO_ID, row['segi'],
                      str(row['residues']) if row['residues'] else '',
                      numbers, str(row['others']) if row['others'] else '',
                      row['sequence']]
            for col, value in enumerate(values):
                item = QtWidgets.QTableWidgetItem(value)
                if col in (RESIDUES, OTHER):
                    item.setTextAlignment(right)
                table.setItem(k, col, item)
        table.resizeColumnsToContents()
        table.setColumnWidth(OBJECT, min(table.columnWidth(OBJECT), 320))
        table.horizontalHeaderItem(OTHER).setToolTip(
            'Ligands, ions and waters with the chain\'s ID')
        mode = QtCore.QItemSelectionModel.SelectionFlag
        for k, row in enumerate(self.rows):
            if (row['object'], row['segi'], row['chain']) in select:
                table.selectionModel().select(
                    table.model().index(k, 0), mode.Select | mode.Rows)
        table.blockSignals(False)
        self._filter()
        self._selection_changed()

    def _filter(self, *_):
        pattern = self.filter.text().strip()
        for k, row in enumerate(self.rows):
            hidden = bool(pattern) and not fnmatch.fnmatch(
                row['object'], pattern if any(c in pattern for c in '*?[')
                else '*%s*' % pattern)
            self.table.setRowHidden(k, hidden)

    def _select_all(self, state):
        self.table.clearSelection()
        if state:
            mode = QtCore.QItemSelectionModel.SelectionFlag
            for k in range(len(self.rows)):
                if not self.table.isRowHidden(k):
                    self.table.selectionModel().select(
                        self.table.model().index(k, 0),
                        mode.Select | mode.Rows)

    def selected_rows(self):
        '''
        The rows of the chains selected (and not filtered out)
        '''
        if not self.rows:
            return []
        chosen = sorted({index.row() for index in
                         self.table.selectionModel().selectedRows()})
        return [self.rows[k] for k in chosen
                if k < len(self.rows) and not self.table.isRowHidden(k)]

    def _selection_changed(self):
        rows = self.selected_rows()
        chains = list(dict.fromkeys(r['chain'] for r in rows))
        self.order.setText(' '.join(c for c in chains if c)
                           if len(rows) > 1 else '')
        self.swap_button.setEnabled(len(chains) == 2)

    def _numbering_changed(self, *_):
        self.gap.setEnabled(self.numbering.currentData() == 'continue')

    def _zoom(self, row, _column):
        if row < len(self.rows):
            self.cmd.zoom(self.rows[row]['selection'], 3, animate=0.5)

    def showEvent(self, event):
        self.refresh()
        super().showEvent(event)

    def event(self, event):
        # back from the main window: something may have been loaded,
        # deleted or altered there
        if event.type() == QtCore.QEvent.Type.WindowActivate and \
                self._counts is not None and self.isVisible() and \
                self._counts != sequence.get_change_counts(_self=self.cmd):
            self.refresh()
        return super().event(event)

    # what to change

    def selection(self, polymer=False):
        '''
        The chains selected in the table (and of them the residues
        typed), as a selection
        '''
        rows = self.selected_rows()
        if not rows:
            raise ValueError('select one or more chains in the table first')
        total = {}
        for row in self.rows:
            total[row['object']] = total.get(row['object'], 0) + 1
        chosen = {}
        for row in rows:
            chosen.setdefault(row['object'], []).append(row)
        parts = []
        for obj, mine in chosen.items():
            if len(mine) == total[obj]:
                parts.append('%' + obj)
            else:
                parts.append('(%%%s and (%s))' % (obj, ' or '.join(
                    '(segi "%s" and chain "%s")' % (r['segi'], r['chain'])
                    for r in mine)))
        selection = '(%s)' % ' or '.join(parts)
        residues = modify.resi_selection(self.residues.text())
        if residues:
            selection += ' and ' + residues
        if polymer:
            selection += ' and polymer'
        return selection

    def _run(self, change):
        '''
        Make a change (a function returning a modify.Change) as one undo
        step, and show what it did
        '''
        try:
            with undo.stack.step(True):
                result = change()
        except Exception as e:
            self.status.setText('Nothing changed: %s' % e)
            return None
        text = result.text[:1].upper() + result.text[1:]
        self.status.setText(text)
        print(' Modify Sequences: %s' % result.text)
        self.refresh(select=result.chains if result.count else None)
        return result

    def renumber(self):
        mode, number = self.mode.currentData(), self.number.value()
        return self._run(lambda: modify.renumber(
            self.selection(polymer=not self.hetero.isChecked()),
            start=None if mode == 'offset' else number,
            offset=number if mode == 'offset' else None,
            gaps=mode != 'count', _self=self.cmd))

    def change_ids(self):
        def change():
            if not (self.chain_check.isChecked() or
                    self.segi_check.isChecked()):
                raise ValueError('check Chain ID, Segment ID or both')
            return modify.set_ids(
                self.selection(),
                chain=self.chain_id.text() if self.chain_check.isChecked()
                else None,
                segi=self.segi_id.text() if self.segi_check.isChecked()
                else None, _self=self.cmd)
        return self._run(change)

    def swap(self):
        def change():
            chains = list(dict.fromkeys(
                r['chain'] for r in self.selected_rows()))
            if len(chains) != 2:
                raise ValueError('select two chains with different chain '
                                 'IDs to swap them')
            return modify.swap(self.selection(), chains[0], chains[1],
                               _self=self.cmd)
        return self._run(change)

    def merge(self):
        return self._run(lambda: modify.merge(
            self.selection(), self.merge_chain.text().strip() or None,
            self.numbering.currentData(), self.gap.value(),
            order=self.order.text().replace(',', ' ').split(),
            _self=self.cmd))

    def combine(self):
        return self._run(lambda: modify.combine(
            self.new_name.text(), self.selection(),
            self.delete_sources.isChecked(), _self=self.cmd))

    def undo(self):
        if undo.stack.can_undo():
            self.cmd.undo()
            self.status.setText('Undone')
        else:
            self.status.setText('Nothing to undo')
        self.refresh()


def show_dialog():
    global _dialog
    if _dialog is None:
        from pymol.gui import get_qtwindow
        _dialog = ModifyDialog(get_qtwindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
