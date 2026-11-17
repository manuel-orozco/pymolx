'''
Epitope Prediction dialog: Structure Prediction's, with each chain's role,
the seeds and samples that an antibody-antigen prediction needs, and the
choice of leaving the antibody chains without an MSA (pymolx.epitope)
'''
import os

from pymol import cmd
from pymol.Qt import QtWidgets

from pmg_tk.startup.prediction_gui.main import PredictionDialog
from pymolx import epitope, prediction

_dialog = None


class EpitopeDialog(PredictionDialog):
    '''
    An antibody and its antigen folded together by AlphaFold 3
    '''

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent, _self=_self)
        self.setObjectName('epitope_dialog')
        self.setWindowTitle('Epitope Prediction')
        self.job = epitope.Job('epitope')
        self.name.setText('epitope')
        self.refresh()

    # the settings: roles, seeds, samples, the antibody's MSA

    def _settings_box(self):
        box = QtWidgets.QGroupBox('AlphaFold 3 cofolding', self)
        box.setObjectName('prep_step')
        grid = QtWidgets.QGridLayout(box)
        self.name = QtWidgets.QLineEdit('epitope', box)
        grid.addWidget(QtWidgets.QLabel('Job name', box), 0, 0)
        grid.addWidget(self.name, 0, 1)

        self.roles = QtWidgets.QTableWidget(0, 2, box)
        self.roles.setObjectName('epitope_roles')
        self.roles.setHorizontalHeaderLabels(['Chain', 'Is'])
        self.roles.verticalHeader().setVisible(False)
        self.roles.horizontalHeader().setStretchLastSection(True)
        self.roles.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.NoSelection)
        self.roles.setToolTip('What each chain is: guessed from its '
                              'sequence, change it if the guess is wrong')
        grid.addWidget(self.roles, 1, 0, 1, 2)
        grid.setRowStretch(1, 1)

        self.seed_count = QtWidgets.QSpinBox(box)
        self.seed_count.setObjectName('epitope_seeds')
        self.seed_count.setRange(1, epitope.MAX_SEEDS)
        self.seed_count.setValue(epitope.DEFAULT_SEEDS)
        self.seed_count.setToolTip(
            'One prediction per seed. Whether a seed finds the epitope is '
            'largely luck, so predict many and compare them')
        grid.addWidget(QtWidgets.QLabel('Seeds', box), 2, 0)
        grid.addWidget(self.seed_count, 2, 1)

        self.samples = QtWidgets.QSpinBox(box)
        self.samples.setObjectName('epitope_samples')
        self.samples.setRange(1, epitope.MAX_SAMPLES)
        self.samples.setValue(epitope.DEFAULT_SAMPLES)
        self.samples.setToolTip(
            'Structures each seed produces. They share the seed\'s one MSA '
            'search, so they cost little; passed to the run as %s, not '
            'written in the input file' % epitope.SAMPLES_FLAG)
        grid.addWidget(QtWidgets.QLabel('Samples per seed', box), 3, 0)
        grid.addWidget(self.samples, 3, 1)

        self.antibody_msa = QtWidgets.QCheckBox(
            'Search an MSA for the antibody chains', box)
        self.antibody_msa.setObjectName('epitope_antibody_msa')
        self.antibody_msa.setChecked(True)
        self.antibody_msa.setToolTip(
            'Unchecked, the antibody chains are given no MSA and no '
            'templates (AlphaFold 3 takes that as a custom MSA with no '
            'sequences): the slowest part of the pipeline is skipped and '
            'an antibody\'s MSA says little about it, but it does change '
            'what the model is given. The antigen keeps its MSA either way')
        grid.addWidget(self.antibody_msa, 4, 0, 1, 2)

        for signal in (self.seed_count.valueChanged,
                       self.samples.valueChanged,
                       self.antibody_msa.toggled):
            signal.connect(self.refresh)

        self.status = QtWidgets.QLabel(box)
        self.status.setObjectName('prep_note')
        self.status.setWordWrap(True)
        grid.addWidget(self.status, 5, 0, 1, 2)
        self.save_button = QtWidgets.QPushButton('Save Input File...', box)
        self.save_button.setDefault(True)
        self.save_button.clicked.connect(self.save)
        grid.addWidget(self.save_button, 6, 0, 1, 2)
        self.submit_button = QtWidgets.QPushButton('Submit Job...', box)
        self.submit_button.setObjectName('epitope_submit')
        self.submit_button.clicked.connect(self.submit_job)
        grid.addWidget(self.submit_button, 7, 0, 1, 2)
        # (the method is AlphaFold 3: this is cofolding with it)
        self.method = _OneMethod()
        return box

    def _method_changed(self, method):
        pass            # (the seeds are the seed count here)

    def _fill_roles(self):
        '''
        One row per chain, its role in a combo
        '''
        self.roles.setRowCount(len(self.job.chains))
        for row, chain in enumerate(self.job.chains):
            item = QtWidgets.QTableWidgetItem(chain.id)
            item.setFlags(item.flags() & ~QtWidgets.QTableWidgetItem(
                '').flags().ItemIsEditable)
            self.roles.setItem(row, 0, item)
            combo = self.roles.cellWidget(row, 1)
            if combo is None:
                combo = QtWidgets.QComboBox(self.roles)
                for label, value in epitope.ROLES:
                    combo.addItem(label, value)
                combo.currentIndexChanged.connect(self._roles_changed)
                self.roles.setCellWidget(row, 1, combo)
            combo.setProperty('chain', chain.id)
            index = combo.findData(self.job.role(chain.id))
            if index >= 0 and index != combo.currentIndex():
                combo.blockSignals(True)
                combo.setCurrentIndex(index)
                combo.blockSignals(False)
        self.roles.resizeColumnsToContents()

    def _roles_changed(self, *_):
        for row in range(self.roles.rowCount()):
            combo = self.roles.cellWidget(row, 1)
            if combo is not None and combo.property('chain'):
                self.job.set_role(combo.property('chain'),
                                  combo.currentData())
        self.refresh()

    # what the job holds

    def _apply(self):
        self.job.name = self.name.text().strip() or 'epitope'
        self.job.seeds = list(range(1, self.seed_count.value() + 1))
        self.job.samples = self.samples.value()
        self.job.antibody_msa = self.antibody_msa.isChecked()

    def refresh(self):
        if not hasattr(self, 'roles'):
            return super().refresh()
        if not isinstance(self.job, epitope.Job):
            job = epitope.Job(self.job.name)
            job.chains = self.job.chains
            job.ligands = self.job.ligands
            job.disulfides = self.job.disulfides
            job.detect_roles()
            self.job = job
        self._apply()
        self._fill_roles()
        super().refresh()

    def _status(self):
        lines = []
        if self.pending:
            lines.append('Disulfide bond from %s %d: click its partner'
                         % self.pending)
        if self.job.chains:
            lines.append(self.job.summary())
            lines += self.job.warnings('af3')
        if self.message:
            lines.append(self.message)
        return '\n'.join(lines) or 'Open a FASTA file to begin'

    # the file and the run

    def save(self):
        self._apply()
        try:
            text = epitope.write(self.job)
        except Exception as e:
            self.message = 'Not ready: %s' % e
            self.refresh()
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save the input file', '%s_af3.json' % self.job.name,
            'JSON (*.json)')
        if not path:
            return
        with open(path, 'w') as handle:
            handle.write(text)
        self.saved = path
        print(' Epitope Prediction: wrote %s (AlphaFold 3, %d seed(s) x '
              '%d sample(s))' % (path, len(self.job.seeds),
                                 self.job.samples))
        for line in self.job.warnings('af3'):
            print('  ' + line)
        print('  run it with %s %d' % (epitope.SAMPLES_FLAG,
                                       self.job.samples))
        self.message = 'Wrote %s' % os.path.basename(path)
        self.refresh()

    def submit_job(self):
        self._apply()
        super().submit_job()

    def _run(self, command, mode, folder):
        '''
        The run, with the samples flag after the input file's name
        '''
        script = next((part for part in command
                       if part.endswith('.sh')), '')
        clash = self.job.seed_clash(script)
        if clash:
            self.message = 'Not submitted: %s' % clash
            print(' Epitope Prediction: %s' % clash)
            self.refresh()
            return
        forwards = epitope.script_forwards_arguments(script) if script \
            else None
        if forwards:
            command = list(command) + epitope.run_arguments(self.job)
        elif forwards is False:
            print(' Epitope Prediction: the job script does not pass extra '
                  'arguments on ("$@"), so %s %d could not be added; the '
                  'run uses the script\'s own value'
                  % (epitope.SAMPLES_FLAG, self.job.samples))
        super()._run(command, mode, folder)


class _OneMethod:
    '''
    Stands in for Structure Prediction's method combo: here it is always
    AlphaFold 3
    '''

    def currentData(self):
        return 'af3'

    def findData(self, value):
        return 0 if value == 'af3' else -1

    def setCurrentIndex(self, index):
        pass

    def currentIndex(self):
        return 0


def show_dialog():
    global _dialog
    if _dialog is None:
        from pymol.gui import get_qtwindow
        _dialog = EpitopeDialog(get_qtwindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
