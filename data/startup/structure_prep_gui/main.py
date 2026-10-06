'''
Structure Preparation dialog: PDBFixer fixes and OpenMM minimization.
Minimization runs in a worker thread; PyMOL is only read before and
written after it, on the GUI thread.
'''

import traceback

from pymol import cmd
from pymol.Qt import QtCore, QtWidgets

from pymolx import mm

Qt = QtCore.Qt

_dialog = None

HETEROGENS = [('Keep water', 'keep_water'), ('Remove all', 'all'),
              ('Keep all', 'none')]
SOLVENTS = [('Implicit (GBn2)', 'gbn2'), ('Implicit (OBC2)', 'obc2'),
            ('Vacuum', 'vacuum')]
RESTRAINTS = [('Backbone', 'backbone'), ('Heavy atoms', 'heavy'),
              ('None', 'none')]


class MinimizeThread(QtCore.QThread):
    '''
    Runs MinimizationJob.run (OpenMM only) off the GUI thread.
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
        except mm.MinimizationCancelled:
            self.cancelled = True
        except Exception as e:
            self.error = e
            traceback.print_exc()


class StructurePrepDialog(QtWidgets.QDialog):

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent)
        self.cmd = _self
        self.thread = None
        self.setObjectName('structure_prep_dialog')
        self.setWindowTitle('Structure Preparation')
        self.resize(640, 720)

        layout = QtWidgets.QVBoxLayout(self)

        top = QtWidgets.QHBoxLayout()
        top.addWidget(QtWidgets.QLabel('Selection', self))
        self.selection = QtWidgets.QComboBox(self)
        self.selection.setEditable(True)
        self.selection.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                                     QtWidgets.QSizePolicy.Policy.Fixed)
        top.addWidget(self.selection, 1)
        top.addWidget(QtWidgets.QLabel('State', self))
        self.state = QtWidgets.QSpinBox(self)
        self.state.setRange(-1, 99999)
        self.state.setValue(-1)
        self.state.setToolTip('-1 = current state')
        top.addWidget(self.state)
        layout.addLayout(top)

        layout.addWidget(self._fix_group())
        layout.addWidget(self._minimize_group())

        self.log = QtWidgets.QPlainTextEdit(self)
        self.log.setObjectName('structure_prep_log')
        self.log.setReadOnly(True)
        layout.addWidget(self.log, 1)

        if not mm.available():
            self._log('OpenMM and PDBFixer are not installed:\n'
                      '  pip install openmm pdbfixer')
            self.fix_button.setEnabled(False)
            self.minimize_button.setEnabled(False)

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self._poll)

    # layout

    def _fix_group(self):
        group = QtWidgets.QGroupBox('Fix structure (PDBFixer)', self)
        form = QtWidgets.QGridLayout(group)

        self.add_atoms = QtWidgets.QCheckBox('Add missing heavy atoms', group)
        self.add_atoms.setChecked(True)
        self.replace_nonstandard = QtWidgets.QCheckBox(
            'Replace nonstandard residues', group)
        self.replace_nonstandard.setChecked(True)
        form.addWidget(self.add_atoms, 0, 0, 1, 2)
        form.addWidget(self.replace_nonstandard, 0, 2, 1, 2)

        self.add_hydrogens = QtWidgets.QCheckBox('Add hydrogens at pH', group)
        self.add_hydrogens.setChecked(True)
        self.ph = QtWidgets.QDoubleSpinBox(group)
        self.ph.setRange(0.0, 14.0)
        self.ph.setSingleStep(0.5)
        self.ph.setValue(7.0)
        form.addWidget(self.add_hydrogens, 1, 0)
        form.addWidget(self.ph, 1, 1)

        form.addWidget(QtWidgets.QLabel('Heterogens', group), 1, 2)
        self.heterogens = QtWidgets.QComboBox(group)
        for label, value in HETEROGENS:
            self.heterogens.addItem(label, value)
        form.addWidget(self.heterogens, 1, 3)

        self.add_residues = QtWidgets.QCheckBox(
            'Add missing residues, sequence from', group)
        self.source = QtWidgets.QLineEdit(group)
        self.source.setPlaceholderText('original file or PDB ID')
        browse = QtWidgets.QPushButton('...', group)
        browse.setMaximumWidth(32)
        browse.clicked.connect(self._browse_source)
        form.addWidget(self.add_residues, 2, 0, 1, 2)
        form.addWidget(self.source, 2, 2)
        form.addWidget(browse, 2, 3)

        form.addWidget(QtWidgets.QLabel('New object', group), 3, 0)
        self.new_name = QtWidgets.QLineEdit(group)
        self.new_name.setPlaceholderText('<object>_fixed')
        form.addWidget(self.new_name, 3, 1, 1, 2)
        self.fix_button = QtWidgets.QPushButton('Fix Structure', group)
        self.fix_button.clicked.connect(self.fix)
        form.addWidget(self.fix_button, 3, 3)
        return group

    def _minimize_group(self):
        group = QtWidgets.QGroupBox('Minimize (OpenMM)', self)
        form = QtWidgets.QGridLayout(group)

        form.addWidget(QtWidgets.QLabel('Force field', group), 0, 0)
        self.forcefield = QtWidgets.QComboBox(group)
        self.forcefield.addItems(list(mm.FORCEFIELDS))
        form.addWidget(self.forcefield, 0, 1)

        form.addWidget(QtWidgets.QLabel('Solvent', group), 0, 2)
        self.solvent = QtWidgets.QComboBox(group)
        for label, value in SOLVENTS:
            self.solvent.addItem(label, value)
        form.addWidget(self.solvent, 0, 3)

        form.addWidget(QtWidgets.QLabel('Restraints', group), 1, 0)
        self.restrain = QtWidgets.QComboBox(group)
        for label, value in RESTRAINTS:
            self.restrain.addItem(label, value)
        form.addWidget(self.restrain, 1, 1)

        form.addWidget(QtWidgets.QLabel('k (kJ/mol/nm²)', group), 1, 2)
        self.restraint_k = QtWidgets.QDoubleSpinBox(group)
        self.restraint_k.setRange(0.0, 1e6)
        self.restraint_k.setDecimals(0)
        self.restraint_k.setSingleStep(100.0)
        self.restraint_k.setValue(1000.0)
        form.addWidget(self.restraint_k, 1, 3)

        form.addWidget(QtWidgets.QLabel('Max iterations', group), 2, 0)
        self.max_iterations = QtWidgets.QSpinBox(group)
        self.max_iterations.setRange(0, 1000000)
        self.max_iterations.setSpecialValueText('until converged')
        form.addWidget(self.max_iterations, 2, 1)

        form.addWidget(QtWidgets.QLabel('Platform', group), 2, 2)
        self.platform = QtWidgets.QComboBox(group)
        self.platform.addItems(['auto'] + list(mm.PLATFORMS))
        form.addWidget(self.platform, 2, 3)

        buttons = QtWidgets.QHBoxLayout()
        self.progress = QtWidgets.QLabel(group)
        self.progress.setObjectName('structure_prep_progress')
        buttons.addWidget(self.progress, 1)
        self.cancel_button = QtWidgets.QPushButton('Cancel', group)
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self.cancel)
        buttons.addWidget(self.cancel_button)
        self.minimize_button = QtWidgets.QPushButton('Minimize', group)
        self.minimize_button.clicked.connect(self.minimize)
        buttons.addWidget(self.minimize_button)
        form.addLayout(buttons, 3, 0, 1, 4)

        self.solvent.currentIndexChanged.connect(self._check_solvent)
        self.forcefield.currentIndexChanged.connect(self._check_solvent)
        self._check_solvent()
        return group

    def _check_solvent(self, *_):
        # CHARMM36 here has no implicit solvent files
        implicit_ok = self.forcefield.currentText() != 'charmm36'
        model = self.solvent.model()
        for row, (_, value) in enumerate(SOLVENTS):
            item = model.item(row)
            item.setEnabled(value == 'vacuum' or implicit_ok)
        if not implicit_ok and self.solvent.currentData() != 'vacuum':
            self.solvent.setCurrentIndex(len(SOLVENTS) - 1)

    def _browse_source(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Original structure', '',
            'Structures (*.pdb *.cif *.ent *.mmcif);;All files (*)')
        if filename:
            self.source.setText(filename)
            self.add_residues.setChecked(True)

    # session

    def refresh(self):
        text = self.selection.currentText()
        names = (self.cmd.get_names('public_objects') +
                 self.cmd.get_names('public_selections'))
        self.selection.clear()
        self.selection.addItems(names)
        self.selection.setEditText(text if text else
                                   (names[0] if names else 'all'))

    def showEvent(self, event):
        self.refresh()
        super().showEvent(event)

    def _log(self, text):
        self.log.appendPlainText(text)

    # actions

    def fix(self):
        selection = self.selection.currentText() or 'all'
        QtWidgets.QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            report = mm.fix_structure(
                selection, self.new_name.text().strip(),
                add_missing_residues=int(self.add_residues.isChecked()),
                add_missing_atoms=int(self.add_atoms.isChecked()),
                replace_nonstandard=int(self.replace_nonstandard.isChecked()),
                heterogens=self.heterogens.currentData(),
                add_hydrogens=int(self.add_hydrogens.isChecked()),
                ph=self.ph.value(), source=self.source.text().strip(),
                state=self.state.value(), quiet=1, _self=self.cmd)
            self._log(mm.format_fix_report(report))
            self.refresh()
            self.selection.setEditText(report['object'])
        except Exception as e:
            self._log(' fix_structure failed: %s' % e)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()

    def minimize(self):
        selection = self.selection.currentText() or 'all'
        try:
            job = mm.prepare_minimization(
                selection,
                forcefield=self.forcefield.currentText(),
                solvent=self.solvent.currentData(),
                restrain=self.restrain.currentData(),
                restraint_k=self.restraint_k.value(),
                max_iterations=self.max_iterations.value(),
                platform=self.platform.currentText(),
                ph=self.ph.value(), state=self.state.value(), _self=self.cmd)
        except Exception as e:
            self._log(' minimize failed: %s' % e)
            return
        self.thread = MinimizeThread(job, self)
        self.thread.finished.connect(self._finished)
        self._set_running(True)
        self.progress.setText('Setting up...')
        self.thread.start()
        self.timer.start(200)

    def cancel(self):
        if self.thread is not None:
            self.thread.job.cancel()
            self.progress.setText('Cancelling...')

    def _poll(self):
        thread = self.thread
        if thread is not None and thread.energy is not None:
            self.progress.setText('Iteration %d, energy %.1f kJ/mol' % (
                thread.iteration, thread.energy))

    def _finished(self):
        self.timer.stop()
        thread, self.thread = self.thread, None
        self._set_running(False)
        if thread.cancelled:
            self.progress.setText('Cancelled, coordinates unchanged')
        elif thread.error is not None:
            self.progress.setText('Failed')
            self._log(' minimize failed: %s' % thread.error)
        else:
            report = thread.job.apply(self.cmd)
            self.progress.setText('Done')
            self._log(mm.format_minimize_report(report))
        thread.deleteLater()

    def _set_running(self, running):
        self.minimize_button.setEnabled(not running)
        self.fix_button.setEnabled(not running)
        self.cancel_button.setEnabled(running)

    def closeEvent(self, event):
        # let a running minimization finish its current step
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
