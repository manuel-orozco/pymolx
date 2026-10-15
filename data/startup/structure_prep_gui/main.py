'''
Structure Preparation dialog, a workflow like Maestro's Protein
Preparation Wizard (pymolx.prep): specify the structure, preprocess,
protonate, minimize, then diagnostics (pymolx.diagnostics) in a tab of
their own.

The preparation runs in a worker thread; PyMOL is only read before and
written after it, on the GUI thread.
'''

import traceback

from pymol import cmd
from pymol.Qt import QtCore, QtWidgets

from pymolx import diagnostics, mm, prep

Qt = QtCore.Qt

_dialog = None

SOLVENTS = [('Implicit (GBn2)', 'gbn2'), ('Implicit (OBC2)', 'obc2'),
            ('Implicit (HCT)', 'hct'), ('Explicit water box', 'explicit'),
            ('Vacuum', 'vacuum')]
RESTRAINTS = [('Heavy atoms', 'heavy'), ('Backbone', 'backbone'),
              ('C-alpha', 'ca'), ('None', 'none')]
PKA_METHODS = [('PROPKA 3', 'propka'), ('Standard values', 'standard')]


class PrepThread(QtCore.QThread):
    '''
    Runs PrepJob.run (no PyMOL calls) off the GUI thread.
    '''

    def __init__(self, job, parent=None):
        super().__init__(parent)
        self.job = job
        self.error = None
        self.cancelled = False
        self.iteration = 0
        self.energy = None
        job.progress = self._progress

    def _progress(self, iteration, energy):
        self.iteration = iteration
        self.energy = energy

    def run(self):
        try:
            self.job.run()
        except prep.Cancelled:
            self.cancelled = True
        except Exception as e:
            self.error = e
            traceback.print_exc()


def _combo(parent, items):
    combo = QtWidgets.QComboBox(parent)
    for label, value in items:
        combo.addItem(label, value)
    return combo


def _set_item_enabled(combo, value, enabled, tooltip=''):
    index = combo.findData(value)
    item = combo.model().item(index)
    item.setEnabled(enabled)
    item.setToolTip('' if enabled else tooltip)


def _first_enabled(combo):
    model = combo.model()
    for row in range(combo.count()):
        if model.item(row).isEnabled():
            combo.setCurrentIndex(row)
            return


class StructurePrepDialog(QtWidgets.QDialog):

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent)
        self.cmd = _self
        self.thread = None
        self.issues = []
        self._log_lines = 0
        self.setObjectName('structure_prep_dialog')
        self.setWindowTitle('Structure Preparation')
        self.resize(600, 820)

        layout = QtWidgets.QVBoxLayout(self)
        self.tabs = QtWidgets.QTabWidget(self)
        layout.addWidget(self.tabs)
        self.tabs.addTab(self._workflow_tab(), 'Preparation Workflow')
        self.tabs.addTab(self._diagnostics_tab(), 'Diagnostics')

        self._update_compatibility()
        if not mm.available():
            self._log('OpenMM and PDBFixer are not installed:\n'
                      '  pip install openmm pdbfixer')
            self.run_button.setEnabled(False)

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self._poll)

    # workflow tab

    def _workflow_tab(self):
        tab = QtWidgets.QWidget(self)
        outer = QtWidgets.QVBoxLayout(tab)
        outer.setContentsMargins(4, 4, 4, 4)

        steps = QtWidgets.QWidget(tab)
        column = QtWidgets.QVBoxLayout(steps)
        column.setContentsMargins(0, 0, 0, 0)
        column.addWidget(self._structure_step(steps))
        column.addWidget(self._preprocess_step(steps))
        column.addWidget(self._protonation_step(steps))
        column.addWidget(self._minimization_step(steps))
        column.addWidget(self._diagnostics_step(steps))
        column.addStretch(1)
        scroll = QtWidgets.QScrollArea(tab)
        scroll.setObjectName('prep_steps')
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        scroll.setWidget(steps)

        # steps above, log below; the divider can be dragged
        split = QtWidgets.QSplitter(Qt.Orientation.Vertical, tab)
        split.setChildrenCollapsible(False)
        split.addWidget(scroll)
        bottom = QtWidgets.QWidget(split)
        below = QtWidgets.QVBoxLayout(bottom)
        below.setContentsMargins(0, 0, 0, 0)
        split.addWidget(bottom)
        outer.addWidget(split, 1)

        bar = QtWidgets.QHBoxLayout()
        self.progress = QtWidgets.QLabel(tab)
        self.progress.setObjectName('structure_prep_progress')
        bar.addWidget(self.progress, 1)
        self.cancel_button = QtWidgets.QPushButton('Cancel', tab)
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self.cancel)
        bar.addWidget(self.cancel_button)
        self.run_button = QtWidgets.QPushButton('Run', tab)
        self.run_button.setDefault(True)
        self.run_button.clicked.connect(self.run)
        bar.addWidget(self.run_button)
        below.addLayout(bar)

        self.log = QtWidgets.QPlainTextEdit(bottom)
        self.log.setObjectName('structure_prep_log')
        self.log.setReadOnly(True)
        self.log.setMinimumHeight(60)
        below.addWidget(self.log, 1)
        split.setStretchFactor(0, 4)
        split.setStretchFactor(1, 1)
        split.setSizes([620, 160])
        return tab

    @staticmethod
    def _step(title, parent, checkable=True):
        box = QtWidgets.QGroupBox(title, parent)
        box.setObjectName('prep_step')
        box.setCheckable(checkable)
        if checkable:
            box.setChecked(True)
        return box

    def _structure_step(self, parent):
        box = self._step('Specify Structure', parent, checkable=False)
        grid = QtWidgets.QGridLayout(box)
        grid.addWidget(QtWidgets.QLabel('Structure', box), 0, 0)
        self.selection = QtWidgets.QComboBox(box)
        self.selection.setEditable(True)
        self.selection.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                                     QtWidgets.QSizePolicy.Policy.Fixed)
        grid.addWidget(self.selection, 0, 1, 1, 2)
        fetch = QtWidgets.QPushButton('Get PDB...', box)
        fetch.clicked.connect(self._fetch)
        grid.addWidget(fetch, 0, 3)

        grid.addWidget(QtWidgets.QLabel('State', box), 1, 0)
        self.state = QtWidgets.QSpinBox(box)
        self.state.setRange(-1, 99999)
        self.state.setValue(-1)
        self.state.setSpecialValueText('current')
        grid.addWidget(self.state, 1, 1)
        grid.addWidget(QtWidgets.QLabel('New object', box), 2, 0)
        self.new_name = QtWidgets.QLineEdit(box)
        self.new_name.setPlaceholderText('<object>_prep')
        grid.addWidget(self.new_name, 2, 1, 1, 3)
        grid.setColumnStretch(2, 1)
        return box

    def _preprocess_step(self, parent):
        box = self.preprocess = self._step('Preprocess', parent)
        grid = QtWidgets.QGridLayout(box)
        self.delete_hydrogens = QtWidgets.QCheckBox(
            'Delete existing hydrogens', box)
        self.delete_hydrogens.setChecked(True)
        self.delete_hydrogens.setToolTip(
            'Unchecked: hydrogens with names OpenMM knows are kept (e.g. '
            'from an earlier preparation); others are added again')
        grid.addWidget(self.delete_hydrogens, 0, 0, 1, 2)
        self.cap_termini = QtWidgets.QCheckBox(
            'Cap termini (ACE/NME)', box)
        self.fill_side_chains = QtWidgets.QCheckBox(
            'Fill in missing side chains', box)
        self.fill_side_chains.setChecked(True)
        grid.addWidget(self.cap_termini, 1, 0)
        grid.addWidget(self.fill_side_chains, 1, 1)

        self.fill_loops = QtWidgets.QCheckBox(
            'Fill in missing residues', box)
        self.fill_loops.setChecked(True)
        self.fill_termini = QtWidgets.QCheckBox(
            'including chain ends', box)
        self.fill_termini.setToolTip(
            'Also build missing residues before the first and after the '
            'last modeled residue (often long, floppy tails)')
        grid.addWidget(self.fill_loops, 2, 0)
        grid.addWidget(self.fill_termini, 2, 1)

        source = QtWidgets.QHBoxLayout()
        source.addWidget(QtWidgets.QLabel('Sequence from', box))
        self.source = QtWidgets.QLineEdit(box)
        self.source.setPlaceholderText(
            'the structure (default), or a file or PDB ID')
        self.source.setToolTip(
            'Missing residues are known from the loaded file (PDB REMARK '
            '465, mmCIF); give the original file or a PDB ID otherwise')
        source.addWidget(self.source, 1)
        browse = QtWidgets.QPushButton('...', box)
        browse.setMaximumWidth(32)
        browse.clicked.connect(self._browse_source)
        source.addWidget(browse)
        grid.addLayout(source, 3, 0, 1, 2)
        self.fill_loops.toggled.connect(self.fill_termini.setEnabled)
        self.fill_loops.toggled.connect(self.source.setEnabled)

        self.replace_nonstandard = QtWidgets.QCheckBox(
            'Replace nonstandard residues', box)
        self.replace_nonstandard.setChecked(True)
        self.replace_nonstandard.setToolTip('e.g. MSE -> MET, SEP -> SER')
        grid.addWidget(self.replace_nonstandard, 4, 0)
        altloc = QtWidgets.QHBoxLayout()
        altloc.addWidget(QtWidgets.QLabel('Alternate locations', box))
        self.altloc = _combo(box, [('Highest occupancy', 'occupancy')] + [
            (letter, letter) for letter in 'ABCD'])
        self.altloc.setToolTip('Which alternate conformation to keep (a '
                               'letter: where present, else the highest '
                               'occupancy)')
        altloc.addWidget(self.altloc)
        altloc.addStretch(1)
        grid.addLayout(altloc, 4, 1)

        self.delete_heterogens = QtWidgets.QCheckBox(
            'Delete other heterogens', box)
        self.delete_heterogens.setToolTip('Ligands, ions, other non-polymer '
                                          'residues')
        grid.addWidget(self.delete_heterogens, 5, 0)
        self.disulfides = QtWidgets.QCheckBox('Create disulfide bonds', box)
        self.disulfides.setChecked(True)
        self.disulfides.setToolTip('Cysteine pairs close enough are bonded '
                                   '(unchecked: free cysteines)')
        grid.addWidget(self.disulfides, 5, 1)
        self.metals = QtWidgets.QCheckBox(
            'Treat metal sites (zero-order bonds, coordinating residues\' '
            'states, coordination held)', box)
        self.metals.setChecked(True)
        self.metals.setToolTip(
            'Atoms within coordination distance of metal ions get zero-order '
            'bonds; coordinating Cys become thiolates, His are protonated on '
            'the other nitrogen; the minimization keeps the distances')
        grid.addWidget(self.metals, 8, 0, 1, 2)
        mutate = QtWidgets.QHBoxLayout()
        mutate.addWidget(QtWidgets.QLabel('Mutate residues', box))
        self.mutations = QtWidgets.QLineEdit(box)
        self.mutations.setPlaceholderText('e.g. A:45:GLY, A:102:W')
        self.mutations.setToolTip('chain:residue:new amino acid (one- or '
                                  'three-letter code); several separated by '
                                  'commas or spaces')
        mutate.addWidget(self.mutations, 1)
        grid.addLayout(mutate, 9, 0, 1, 2)
        waters = QtWidgets.QHBoxLayout()
        self.delete_waters = QtWidgets.QCheckBox('Delete waters', box)
        waters.addWidget(self.delete_waters)
        self.water_far = QtWidgets.QCheckBox('only beyond', box)
        waters.addWidget(self.water_far)
        self.water_distance = QtWidgets.QDoubleSpinBox(box)
        self.water_distance.setRange(1.0, 50.0)
        self.water_distance.setDecimals(1)
        self.water_distance.setValue(5.0)
        self.water_distance.setSuffix(' A')
        waters.addWidget(self.water_distance)
        waters.addWidget(QtWidgets.QLabel('from het groups', box))
        waters.addStretch(1)
        grid.addLayout(waters, 6, 0, 1, 2)
        hbonds = QtWidgets.QHBoxLayout()
        self.water_hbonds = QtWidgets.QCheckBox(
            'Delete waters with fewer than', box)
        hbonds.addWidget(self.water_hbonds)
        self.water_min_hbonds = QtWidgets.QSpinBox(box)
        self.water_min_hbonds.setRange(1, 4)
        self.water_min_hbonds.setValue(2)
        hbonds.addWidget(self.water_min_hbonds)
        hbonds.addWidget(QtWidgets.QLabel('H-bonds to non-waters', box))
        hbonds.addStretch(1)
        grid.addLayout(hbonds, 7, 0, 1, 2)

        def waters_changed(*_):
            on = self.delete_waters.isChecked()
            self.water_far.setEnabled(on)
            self.water_distance.setEnabled(on and self.water_far.isChecked())
            self.water_min_hbonds.setEnabled(self.water_hbonds.isChecked())
        for check in (self.delete_waters, self.water_far, self.water_hbonds):
            check.toggled.connect(waters_changed)
        waters_changed()
        self.cap_termini.toggled.connect(self._update_compatibility)
        box.toggled.connect(self._update_compatibility)
        return box

    def _protonation_step(self, parent):
        box = self.protonate = self._step('Protonation', parent)
        grid = QtWidgets.QGridLayout(box)
        grid.addWidget(QtWidgets.QLabel('pH', box), 0, 0)
        self.ph = QtWidgets.QDoubleSpinBox(box)
        self.ph.setRange(0.0, 14.0)
        self.ph.setSingleStep(0.5)
        self.ph.setDecimals(1)
        self.ph.setValue(7.0)
        grid.addWidget(self.ph, 0, 1)
        grid.addWidget(QtWidgets.QLabel('pKa values', box), 0, 2)
        self.pka_method = _combo(box, PKA_METHODS)
        if not prep.propka_available():
            _set_item_enabled(self.pka_method, 'propka', False,
                              'PROPKA is not installed: pip install propka')
            _first_enabled(self.pka_method)
        grid.addWidget(self.pka_method, 0, 3)
        self.optimize_hbonds = QtWidgets.QCheckBox(
            'Optimize hydrogen bonds', box)
        self.optimize_hbonds.setChecked(True)
        self.optimize_hbonds.setToolTip(
            'Flip Asn, Gln and His side chains, choose His tautomers, turn '
            'OH, SH and NH3+ hydrogens and waters for hydrogen bonds')
        grid.addWidget(self.optimize_hbonds, 1, 0, 1, 4)
        self.ligand_hydrogens = QtWidgets.QCheckBox(
            'Hydrogens for ligands (PDB Chemical Component Dictionary)', box)
        self.ligand_hydrogens.setChecked(True)
        self.ligand_hydrogens.setToolTip(
            'Ligand definitions are downloaded from the PDB (needs internet)')
        grid.addWidget(self.ligand_hydrogens, 2, 0, 1, 4)
        note = QtWidgets.QLabel(
            'Asp, Glu, His, Lys states from the pKa values; flips, His '
            'tautomers and rotatable hydrogens from hydrogen bonds', box)
        note.setObjectName('prep_note')
        note.setWordWrap(True)
        grid.addWidget(note, 3, 0, 1, 4)
        grid.setColumnStretch(3, 1)
        return box

    def _minimization_step(self, parent):
        box = self.minimize = self._step('Restrained Minimization', parent)
        grid = QtWidgets.QGridLayout(box)
        grid.addWidget(QtWidgets.QLabel('Force field', box), 0, 0)
        self.forcefield = _combo(box, [
            (info['label'], name)
            for name, info in mm.FORCEFIELD_INFO.items()])
        grid.addWidget(self.forcefield, 0, 1)
        grid.addWidget(QtWidgets.QLabel('Solvent', box), 0, 2)
        self.solvent = _combo(box, SOLVENTS)
        self.solvent.setToolTip('Explicit: a water box (1 nm, 0.15 M NaCl) '
                                'during the minimization, removed afterwards')
        grid.addWidget(self.solvent, 0, 3)

        grid.addWidget(QtWidgets.QLabel('Restraints', box), 1, 0)
        self.restrain = _combo(box, RESTRAINTS)
        self.restrain.setToolTip('Harmonic restraints of the input atoms to '
                                 'their starting positions (built atoms move '
                                 'freely)')
        grid.addWidget(self.restrain, 1, 1)
        grid.addWidget(QtWidgets.QLabel('k (kJ/mol/nm²)', box), 1, 2)
        self.restraint_k = QtWidgets.QDoubleSpinBox(box)
        self.restraint_k.setRange(0.0, 1e6)
        self.restraint_k.setDecimals(0)
        self.restraint_k.setSingleStep(100.0)
        self.restraint_k.setValue(1000.0)
        grid.addWidget(self.restraint_k, 1, 3)

        grid.addWidget(QtWidgets.QLabel('Convergence', box), 2, 0)
        self.tolerance = QtWidgets.QDoubleSpinBox(box)
        self.tolerance.setRange(0.01, 1000.0)
        self.tolerance.setDecimals(2)
        self.tolerance.setValue(10.0)
        self.tolerance.setSuffix(' kJ/mol/nm')
        self.tolerance.setToolTip('Stop when the root-mean-square force is '
                                  'below this')
        grid.addWidget(self.tolerance, 2, 1)
        grid.addWidget(QtWidgets.QLabel('Max iterations', box), 2, 2)
        self.max_iterations = QtWidgets.QSpinBox(box)
        self.max_iterations.setRange(0, 1000000)
        self.max_iterations.setSingleStep(100)
        self.max_iterations.setSpecialValueText('until converged')
        grid.addWidget(self.max_iterations, 2, 3)

        grid.addWidget(QtWidgets.QLabel('Platform', box), 3, 0)
        self.platform = QtWidgets.QComboBox(box)
        self.platform.addItems(['auto'] + list(mm.PLATFORMS))
        grid.addWidget(self.platform, 3, 1)
        self.restrain.currentIndexChanged.connect(
            lambda _: self.restraint_k.setEnabled(
                self.restrain.currentData() != 'none'))
        self.forcefield.currentIndexChanged.connect(
            self._update_compatibility)
        box.toggled.connect(self._update_compatibility)
        return box

    def _diagnostics_step(self, parent):
        box = self.run_diagnostics = self._step('Diagnostics', parent)
        row = QtWidgets.QVBoxLayout(box)
        note = QtWidgets.QLabel(
            'Report missing atoms and residues, clashes, backbone '
            'outliers, cis peptides and more for the prepared structure',
            box)
        note.setObjectName('prep_note')
        note.setWordWrap(True)
        row.addWidget(note)
        return box

    def _update_compatibility(self, *_):
        '''
        CHARMM here: no ACE/NME caps, no implicit solvent
        '''
        info = mm.FORCEFIELD_INFO[self.forcefield.currentData()]
        for value in mm.IMPLICIT_SOLVENTS:
            _set_item_enabled(self.solvent, value, info['implicit'],
                              'Implicit solvent needs an AMBER force field')
        if not self.solvent.model().item(
                self.solvent.currentIndex()).isEnabled():
            self.solvent.setCurrentIndex(self.solvent.findData('explicit'))
        caps_ok = info['caps'] or not self.minimize.isChecked()
        if not caps_ok:
            self.cap_termini.setChecked(False)
        self.cap_termini.setEnabled(caps_ok)
        self.cap_termini.setToolTip(
            '' if caps_ok else '%s has no ACE/NME caps here: use an AMBER '
            'force field' % info['label'])

    def _browse_source(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Original structure', '',
            'Structures (*.pdb *.cif *.ent *.mmcif);;All files (*)')
        if filename:
            self.source.setText(filename)
            self.fill_loops.setChecked(True)

    def _fetch(self):
        code, ok = QtWidgets.QInputDialog.getText(
            self, 'Get PDB', 'PDB ID:')
        code = code.strip()
        if not ok or not code:
            return
        QtWidgets.QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            self.cmd.fetch(code)
            self.refresh()
            self.selection.setEditText(code.lower())
            self.diag_selection.setEditText(code.lower())
        except Exception as e:
            self._log('Get PDB %s failed: %s' % (code, e))
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()

    # diagnostics tab

    def _diagnostics_tab(self):
        tab = QtWidgets.QWidget(self)
        layout = QtWidgets.QVBoxLayout(tab)
        layout.setContentsMargins(4, 4, 4, 4)
        top = QtWidgets.QHBoxLayout()
        top.addWidget(QtWidgets.QLabel('Structure', tab))
        self.diag_selection = QtWidgets.QComboBox(tab)
        self.diag_selection.setEditable(True)
        self.diag_selection.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed)
        top.addWidget(self.diag_selection, 1)
        self.diag_button = QtWidgets.QPushButton('Run Diagnostics', tab)
        self.diag_button.clicked.connect(self.diagnose)
        top.addWidget(self.diag_button)
        layout.addLayout(top)

        self.diag_summary = QtWidgets.QLabel(tab)
        self.diag_summary.setObjectName('prep_note')
        self.diag_summary.setWordWrap(True)
        layout.addWidget(self.diag_summary)

        self.tree = QtWidgets.QTreeWidget(tab)
        self.tree.setObjectName('prep_diagnostics')
        self.tree.setColumnCount(3)
        self.tree.setHeaderLabels(['Issue', 'Severity', 'Details'])
        self.tree.setAlternatingRowColors(False)
        self.tree.itemClicked.connect(self._show_issue)
        header = self.tree.header()
        header.setStretchLastSection(True)
        layout.addWidget(self.tree, 1)
        hint = QtWidgets.QLabel('Click an issue to select ("sele"), show '
                                'and zoom to it', tab)
        hint.setObjectName('prep_note')
        layout.addWidget(hint)
        return tab

    def diagnose(self, selection=None, label=None):
        selection = selection or self.diag_selection.currentText() or 'all'
        QtWidgets.QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            self.issues = diagnostics.run(selection, _self=self.cmd)
        except Exception as e:
            self.issues = []
            self.diag_summary.setText('Diagnostics failed: %s' % e)
            self.tree.clear()
            return []
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self._fill_tree(label or selection)
        return self.issues

    def _fill_tree(self, selection):
        self.tree.clear()
        counts = diagnostics.summary(self.issues)
        total = sum(counts.values())
        self.diag_summary.setText('%s: %s' % (
            selection, '%d issue%s' % (total, '' if total == 1 else 's')
            if total else 'no issues found'))
        for category, description in diagnostics.CATEGORIES.items():
            found = [i for i in self.issues if i.category == category]
            if category == 'Scores' and found:
                parent = QtWidgets.QTreeWidgetItem([category, '', description])
                self.tree.addTopLevelItem(parent)
                for issue in found:
                    parent.addChild(QtWidgets.QTreeWidgetItem(
                        ['', '', issue.detail]))
                parent.setExpanded(True)
                continue
            if category in diagnostics.INFO_CATEGORIES:
                for issue in found:
                    item = QtWidgets.QTreeWidgetItem(
                        [category, issue.severity, issue.detail])
                    item.setData(0, Qt.ItemDataRole.UserRole, issue.selection)
                    self.tree.addTopLevelItem(item)
                continue
            parent = QtWidgets.QTreeWidgetItem(
                ['%s (%d)' % (category, len(found)), '', description])
            parent.setToolTip(2, description)
            self.tree.addTopLevelItem(parent)
            for issue in found:
                child = QtWidgets.QTreeWidgetItem(
                    [issue.residue, issue.severity, issue.detail])
                child.setData(0, Qt.ItemDataRole.UserRole, issue.selection)
                parent.addChild(child)
            parent.setExpanded(0 < len(found) <= 20)
        self.tree.resizeColumnToContents(0)
        self.tree.resizeColumnToContents(1)

    def _show_issue(self, item, _column=0):
        selection = item.data(0, Qt.ItemDataRole.UserRole)
        if not selection:
            return
        try:
            self.cmd.select('sele', selection)
            # (side chains too, not only the cartoon)
            self.cmd.show('sticks', 'sele')
            self.cmd.zoom('sele and present', 4, animate=0.5)
        except Exception as e:
            self.diag_summary.setText('Cannot show the issue: %s' % e)

    # session

    def refresh(self):
        names = (self.cmd.get_names('public_objects') +
                 self.cmd.get_names('public_selections'))
        for combo in (self.selection, self.diag_selection):
            text = combo.currentText()
            combo.clear()
            combo.addItems(names)
            combo.setEditText(text if text else
                              (names[0] if names else 'all'))

    def showEvent(self, event):
        self.refresh()
        super().showEvent(event)

    def _log(self, text):
        self.log.appendPlainText(text)

    # running

    def options(self):
        '''
        :return: prep.prepare_job keyword arguments from the dialog
        '''
        return dict(
            preprocess=self.preprocess.isChecked(),
            altloc=self.altloc.currentData(),
            delete_hydrogens=self.delete_hydrogens.isChecked(),
            cap_termini=self.cap_termini.isChecked(),
            fill_loops=self.fill_loops.isChecked(),
            fill_termini=self.fill_termini.isChecked(),
            source=self.source.text().strip(),
            fill_side_chains=self.fill_side_chains.isChecked(),
            replace_nonstandard=self.replace_nonstandard.isChecked(),
            delete_waters=self.delete_waters.isChecked(),
            water_distance=(self.water_distance.value()
                            if self.water_far.isChecked() else 0.0),
            water_min_hbonds=(self.water_min_hbonds.value()
                              if self.water_hbonds.isChecked() else 0),
            delete_heterogens=self.delete_heterogens.isChecked(),
            mutations=prep.parse_mutations(self.mutations.text()),
            disulfides=self.disulfides.isChecked(),
            metals=self.metals.isChecked(),
            protonate=self.protonate.isChecked(),
            ph=self.ph.value(),
            pka_method=self.pka_method.currentData(),
            ligand_hydrogens=self.ligand_hydrogens.isChecked(),
            optimize_hbonds=self.optimize_hbonds.isChecked(),
            minimize=self.minimize.isChecked(),
            forcefield=self.forcefield.currentData(),
            solvent=self.solvent.currentData(),
            restrain=self.restrain.currentData(),
            restraint_k=self.restraint_k.value(),
            tolerance=self.tolerance.value(),
            max_iterations=self.max_iterations.value(),
            platform=self.platform.currentText(),
            diagnostics=self.run_diagnostics.isChecked())

    def run(self):
        selection = self.selection.currentText().strip() or 'all'
        try:
            job = prep.prepare_job(selection, self.new_name.text().strip(),
                                   self.state.value(), _self=self.cmd,
                                   **self.options())
        except Exception as e:
            self._log('Preparation failed: %s' % e)
            return
        self._log_lines = 0
        self._log('Preparing "%s" -> "%s"' % (selection, job.name))
        self.thread = PrepThread(job, self)
        self.thread.finished.connect(self._finished)
        self._set_running(True)
        self.progress.setText('Starting...')
        self.thread.start()
        self.timer.start(200)

    def cancel(self):
        if self.thread is not None:
            self.thread.job.cancel()
            self.progress.setText('Cancelling...')

    def _flush_log(self, job):
        for line in job.log[self._log_lines:]:
            self._log(line)
        self._log_lines = len(job.log)

    def _poll(self):
        thread = self.thread
        if thread is None:
            return
        self._flush_log(thread.job)
        text = thread.job.status
        if thread.energy is not None and text.startswith('Minimizing'):
            text += ' iteration %d, %.0f kJ/mol' % (thread.iteration,
                                                   thread.energy)
        self.progress.setText(text)

    def _finished(self):
        self.timer.stop()
        thread, self.thread = self.thread, None
        self._set_running(False)
        job = thread.job
        self._flush_log(job)
        if thread.cancelled:
            self.progress.setText('Cancelled, nothing changed')
        elif thread.error is not None:
            self.progress.setText('Failed')
            self._log('Preparation failed: %s' % thread.error)
        else:
            try:
                job.apply(self.cmd)
            except Exception as e:
                self.progress.setText('Failed')
                self._log('Loading the result failed: %s' % e)
            else:
                self._flush_log(job)
                self.refresh()
                self.selection.setEditText(job.name)
                self.diag_selection.setEditText(job.name)
                self.progress.setText('Done: "%s"' % job.name)
                if job.diagnostics:
                    self.diagnose('%' + job.name, job.name)
                    self.diag_selection.setEditText(job.name)
                    counts = diagnostics.summary(self.issues)
                    self._log('Diagnostics: %s' % (', '.join(
                        '%s %d' % (c.lower(), n)
                        for c, n in counts.items() if n) or 'no issues'))
                    self.tabs.setCurrentIndex(1)
        thread.deleteLater()

    def _set_running(self, running):
        self.run_button.setEnabled(not running)
        self.cancel_button.setEnabled(running)

    def closeEvent(self, event):
        # let a running preparation finish its current step
        if self.thread is not None:
            self.thread.job.cancel()
            self.thread.wait()
        super().closeEvent(event)


def show_dialog():
    global _dialog
    if _dialog is None:
        from pymol.gui import get_qtwindow
        _dialog = StructurePrepDialog(get_qtwindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
