'''
Aggregation Prediction dialog: the aggregation propensity of a folded
protein, its pH dependence, the propensity of its fluctuations and
solubilising mutations, shown as tables and colour on the model
(pymolx.aggrescan).

The options follow the Aggrescan4D web server's own panel.
'''

import traceback

from pymol import cmd
from pymol.Qt import QtCore, QtWidgets

from pymolx import aggrescan

Qt = QtCore.Qt

_dialog = None

SCORE_COLUMNS = ['Chain', 'Residue', 'Score', 'Exposure %', 'Prone']
PH_COLUMNS = ['pH', 'Maximum', 'Average', 'Total']
MODEL_COLUMNS = ['Model', 'Average', 'Maximum', 'Total']
MUTATION_COLUMNS = ['Mutation', 'Chain', 'dA4D avg', 'dA4D total',
                    'ddG kcal/mol']

MUTATION_CHOICES = [('None', ''),
                    ('Enhance solubility (charged)', 'charged'),
                    ('Improve solubility (evolutionary)', 'evolutionary')]

ABOUT = (
    'A free re-implementation of Aggrescan4D. Above 0 a residue is '
    'aggregation prone, above %.1f it is likely to matter; a buried '
    'residue scores 0. FoldX is replaced by OpenMM and EvoEF, '
    'CABS-flex by molecular dynamics and pKa-ANI by PROPKA.'
    % aggrescan.HOTSPOT_THRESHOLD)


class AggregationThread(QtCore.QThread):
    '''
    Runs Job.run (no PyMOL calls) off the GUI thread.
    '''

    def __init__(self, job, parent=None):
        super().__init__(parent)
        self.job = job
        self.error = None
        self.cancelled = False

    def run(self):
        try:
            self.job.run()
        except aggrescan.Cancelled:
            self.cancelled = True
        except Exception as e:          # reported in the dialog
            self.error = e
            traceback.print_exc()


def _yes_no(parent, label, checked):
    '''
    The server shows each option as a Yes/No pair; a checkbox says the
    same thing in less room
    '''
    box = QtWidgets.QCheckBox(label, parent)
    box.setChecked(checked)
    return box


class AggregationDialog(QtWidgets.QDialog):

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent)
        self.cmd = _self
        self.job = None
        self.thread = None
        self.result = None
        self.setObjectName('aggregation_dialog')
        self.setWindowTitle('Aggregation Prediction')
        self.resize(980, 660)

        layout = QtWidgets.QHBoxLayout(self)
        left = QtWidgets.QVBoxLayout()
        layout.addLayout(left, 2)
        right = QtWidgets.QVBoxLayout()
        layout.addLayout(right, 3)

        box = QtWidgets.QGroupBox('Structure', self)
        box.setObjectName('prep_step')
        grid = QtWidgets.QGridLayout(box)
        grid.addWidget(QtWidgets.QLabel('Object', box), 0, 0)
        self.object = QtWidgets.QComboBox(box)
        self.object.setObjectName('aggregation_object')
        self.object.currentIndexChanged.connect(self._object_changed)
        grid.addWidget(self.object, 0, 1)
        grid.addWidget(QtWidgets.QLabel('Chain(s)', box), 1, 0)
        self.chain = QtWidgets.QComboBox(box)
        self.chain.setObjectName('aggregation_chain')
        grid.addWidget(self.chain, 1, 1)
        grid.setColumnStretch(1, 1)
        left.addWidget(box)

        box = QtWidgets.QGroupBox('Options', self)
        box.setObjectName('prep_step')
        grid = QtWidgets.QGridLayout(box)
        row = 0
        self.ph = _yes_no(box, 'pH calculations', True)
        self.ph.setObjectName('aggregation_ph')
        self.ph.setToolTip('Score at the 11 pH values from %.1f to %.1f'
                           % (aggrescan.PH_VALUES[0],
                              aggrescan.PH_VALUES[-1]))
        self.ph.toggled.connect(self._options_changed)
        grid.addWidget(self.ph, row, 0, 1, 2)
        row += 1
        self.optimize = _yes_no(box, 'Stability calculations', True)
        self.optimize.setObjectName('aggregation_optimize')
        self.optimize.setToolTip('Minimize the structure first with '
                                 'OpenMM, where A4D runs FoldX')
        self.optimize.toggled.connect(self._options_changed)
        grid.addWidget(self.optimize, row, 0, 1, 2)
        row += 1
        self.globular = _yes_no(box, 'Analysis of the globular regions',
                                False)
        self.globular.setObjectName('aggregation_globular')
        self.globular.setToolTip(
            'Leave out the parts of a predicted model with a pLDDT '
            'below %g, where A4D runs AlphaCutter'
            % aggrescan.PLDDT_CUTOFF)
        grid.addWidget(self.globular, row, 0, 1, 2)
        row += 1
        self.dynamic = _yes_no(box, 'Dynamic mode', False)
        self.dynamic.setObjectName('aggregation_dynamic')
        self.dynamic.setToolTip('Sample fluctuations with molecular '
                                'dynamics and score a set of models')
        self.dynamic.toggled.connect(self._options_changed)
        grid.addWidget(self.dynamic, row, 0, 1, 2)
        row += 1
        grid.addWidget(QtWidgets.QLabel('Distance of aggregation '
                                        'analysis', box), row, 0)
        self.distance = QtWidgets.QComboBox(box)
        self.distance.setObjectName('aggregation_distance')
        for one in aggrescan.DISTANCES:
            self.distance.addItem('%d A' % one, one)
        self.distance.setCurrentIndex(
            list(aggrescan.DISTANCES).index(aggrescan.DEFAULT_DISTANCE))
        grid.addWidget(self.distance, row, 1)
        row += 1
        self.mutate_label = QtWidgets.QLabel('Mutate residues', box)
        grid.addWidget(self.mutate_label, row, 0)
        self.mutate = QtWidgets.QLineEdit(box)
        self.mutate.setObjectName('aggregation_mutate')
        self.mutate.setPlaceholderText('L8E, A/V70R')
        self.mutate.setToolTip('Mutate these residues and analyse the '
                               'mutant. Name the chain as A/L8E when '
                               'there is more than one')
        self.mutate.textChanged.connect(self._options_changed)
        grid.addWidget(self.mutate, row, 1)
        row += 1
        grid.addWidget(QtWidgets.QLabel('Automated mutations', box), row, 0)
        self.mutations = QtWidgets.QComboBox(box)
        self.mutations.setObjectName('aggregation_mutations')
        for label, value in MUTATION_CHOICES:
            self.mutations.addItem(label, value)
        self.mutations.currentIndexChanged.connect(self._options_changed)
        grid.addWidget(self.mutations, row, 1)
        row += 1
        self.sites_label = QtWidgets.QLabel('Residues to mutate', box)
        grid.addWidget(self.sites_label, row, 0)
        self.sites = QtWidgets.QSpinBox(box)
        self.sites.setObjectName('aggregation_sites')
        self.sites.setRange(1, aggrescan.MAX_SITES)
        self.sites.setValue(aggrescan.MAX_SITES)
        grid.addWidget(self.sites, row, 1)
        row += 1
        self.pka_label = QtWidgets.QLabel('pKa values', box)
        grid.addWidget(self.pka_label, row, 0)
        self.pka = QtWidgets.QComboBox(box)
        self.pka.setObjectName('aggregation_pka')
        for label, value in aggrescan.PKA_MODES:
            self.pka.addItem(label, value)
        grid.addWidget(self.pka, row, 1)
        grid.setColumnStretch(1, 1)
        left.addWidget(box)

        self.about = QtWidgets.QLabel(ABOUT, self)
        self.about.setObjectName('prep_note')
        self.about.setWordWrap(True)
        left.addWidget(self.about)
        left.addStretch(1)

        buttons = QtWidgets.QHBoxLayout()
        self.run_button = QtWidgets.QPushButton('Run', self)
        self.run_button.setDefault(True)
        self.run_button.clicked.connect(self.run)
        buttons.addWidget(self.run_button)
        self.cancel_button = QtWidgets.QPushButton('Cancel', self)
        self.cancel_button.clicked.connect(self.cancel)
        self.cancel_button.setEnabled(False)
        buttons.addWidget(self.cancel_button)
        left.addLayout(buttons)

        self.tabs = QtWidgets.QTabWidget(self)
        self.tabs.setObjectName('aggregation_tabs')
        self.scores = self._table(SCORE_COLUMNS, 'aggregation_scores')
        self.tabs.addTab(self.scores, 'Score')
        self.ph_table = self._table(PH_COLUMNS, 'aggregation_ph_table')
        self.ph_table.itemSelectionChanged.connect(self._ph_chosen)
        self.tabs.addTab(self.ph_table, 'pH')
        self.model_table = self._table(MODEL_COLUMNS,
                                       'aggregation_model_table')
        self.tabs.addTab(self.model_table, 'Dynamic mode')
        self.mutation_table = self._table(MUTATION_COLUMNS,
                                          'aggregation_mutation_table')
        self.tabs.addTab(self.mutation_table, 'Mutations')
        right.addWidget(self.tabs, 1)

        self.status = QtWidgets.QLabel(self)
        self.status.setObjectName('prep_note')
        self.status.setWordWrap(True)
        right.addWidget(self.status)
        bar = QtWidgets.QHBoxLayout()
        copy = QtWidgets.QPushButton('Copy Table', self)
        copy.clicked.connect(self.copy)
        bar.addWidget(copy)
        self.save_button = QtWidgets.QPushButton('Save Results...', self)
        self.save_button.setToolTip('Write the scores as CSV and as PDB '
                                    'files carrying them in their '
                                    'B-factors, one per pH tested')
        self.save_button.clicked.connect(self.save)
        bar.addWidget(self.save_button)
        self.load_button = QtWidgets.QPushButton('Load Models', self)
        self.load_button.clicked.connect(self.load_models)
        self.load_button.setEnabled(False)
        bar.addWidget(self.load_button)
        right.addLayout(bar)

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self._poll)

        self._options_changed()
        self.refresh()

    def _table(self, columns, name):
        table = QtWidgets.QTableWidget(0, len(columns), self)
        table.setObjectName(name)
        table.setHorizontalHeaderLabels(columns)
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        table.horizontalHeader().setStretchLastSection(True)
        return table

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
        self.chain.addItem('All chains', '')
        for one in (self.cmd.get_chains('%' + name) if name else []) or []:
            self.chain.addItem(one, one)

    def _options_changed(self, *_):
        '''
        The server refuses the same combinations: the pH scan stands
        apart from the dynamic mode and from automated mutations
        '''
        ph = self.ph.isChecked()
        mutating = bool(self.mutations.currentData())
        self.dynamic.setEnabled(not ph and not mutating)
        self.mutations.setEnabled(not ph and not self.dynamic.isChecked())
        self.ph.setEnabled(not self.dynamic.isChecked() and not mutating)
        self.pka.setEnabled(ph)
        self.pka_label.setEnabled(ph)
        self.sites.setEnabled(mutating)
        self.sites_label.setEnabled(mutating)
        # the server needs the optimisation for any mutation, since
        # that is the step FoldX builds them in
        typed = bool(self.mutate.text().strip())
        if (mutating or typed) and not self.optimize.isChecked():
            self.optimize.setChecked(True)
        self.optimize.setEnabled(not mutating and not typed)
        self.mutate.setEnabled(not mutating)
        self.mutate_label.setEnabled(not mutating)
        missing = aggrescan.available(ph and self.pka.currentData()
                                      == 'propka',
                                      self.optimize.isChecked()
                                      or self.dynamic.isChecked()
                                      or mutating)
        self.run_button.setEnabled(not missing and not self._running())
        if missing:
            self.status.setText(missing)
        elif mutating and not aggrescan.evoef_available():
            self.status.setText('EvoEF was not found, so the mutations '
                                'come without a stability estimate.')

    def _running(self):
        return self.thread is not None and self.thread.isRunning()

    def showEvent(self, event):
        self.refresh()
        super().showEvent(event)

    # running

    def run(self):
        if self._running():
            return
        name = self.object.currentText()
        if not name:
            self.status.setText('Open a structure')
            return
        if self.job is not None:
            self.job.cleanup()      # the run before this one
        chain = self.chain.currentData() or ''
        selection = '%' + name
        if chain:
            selection += ' and chain "%s"' % chain
        try:
            self.job = aggrescan.prepare(
                selection, self.distance.currentData(),
                self.ph.isChecked(), self.pka.currentData(),
                self.optimize.isChecked(), self.dynamic.isChecked(),
                self.mutations.currentData(), self.mutate.text().strip(),
                1, name, self.cmd, globular=self.globular.isChecked(),
                max_sites=self.sites.value())
        except Exception as e:
            self.status.setText(str(e))
            return
        self.job.object = name
        for table in (self.scores, self.ph_table, self.model_table,
                      self.mutation_table):
            table.setRowCount(0)
        self.load_button.setEnabled(False)
        self.thread = AggregationThread(self.job, self)
        self.thread.finished.connect(self._finished)
        self._set_running(True)
        self.status.setText('Starting...')
        self.thread.start()
        self.timer.start(200)

    def cancel(self):
        if self.thread is not None:
            self.thread.job.cancel()
            self.status.setText('Cancelling...')

    def _set_running(self, running):
        self.run_button.setEnabled(not running)
        self.cancel_button.setEnabled(running)

    def _poll(self):
        if self.thread is not None:
            self.status.setText(self.thread.job.status or 'Working...')

    def _finished(self):
        self.timer.stop()
        thread, self.thread = self.thread, None
        self._set_running(False)
        job = thread.job
        if thread.cancelled:
            self.status.setText('Cancelled, nothing changed')
            job.cleanup()
            return
        if thread.error is not None:
            self.status.setText('Failed: %s' % thread.error)
            job.cleanup()
            return
        try:
            shown = job.apply(self.cmd)
        except Exception as e:
            self.status.setText('Showing the result failed: %s' % e)
            shown = {}
        self.result = job.chosen()
        self._fill_tables(job)
        self._say_summary(job, shown)
        # the analysed structure stays until the next run or the
        # dialog closes, so the results can still be written out
        self.load_button.setEnabled(bool(job.models))

    # the results

    def _fill_tables(self, job):
        self._fill(self.scores, [
            (r.chain, '%s%s' % (r.letter, r.resi), '%.4f' % r.score,
             '%.1f' % r.exposure, 'yes' if r.prone else '')
            for r in sorted(job.chosen().residues, key=lambda r: -r.score)])
        self._fill(self.ph_table, [
            ('%.1f' % value, '%.4f' % one.maximum, '%.4f' % one.average,
             '%.4f' % one.total) for value, one in job.ph_results])
        self._fill(self.model_table, [
            (name, '%.4f' % one.average, '%.4f' % one.maximum,
             '%.4f' % one.total) for name, _, one in job.models])
        self._fill(self.mutation_table, [
            (one.label, one.chain, '%.4f' % one.d_average,
             '%.4f' % one.d_total,
             '%+.2f' % one.ddg if one.ddg is not None else '-')
            for one in job.suggestions])
        for one, table in ((job.ph_results, self.ph_table),
                           (job.models, self.model_table),
                           (job.suggestions, self.mutation_table)):
            self.tabs.setTabEnabled(self.tabs.indexOf(table), bool(one))
        if job.ph_results:
            self.tabs.setCurrentWidget(self.ph_table)
        elif job.suggestions:
            self.tabs.setCurrentWidget(self.mutation_table)
        elif job.models:
            self.tabs.setCurrentWidget(self.model_table)
        else:
            self.tabs.setCurrentWidget(self.scores)

    def _fill(self, table, rows):
        table.setRowCount(len(rows))
        for r, values in enumerate(rows):
            for c, text in enumerate(values):
                table.setItem(r, c, QtWidgets.QTableWidgetItem(text))
        table.resizeColumnsToContents()

    def _say_summary(self, job, shown):
        result = job.chosen()
        value, used = result.metric()
        parts = ['%d residues, %d aggregation prone, %d above %.1f.'
                 % (len(result), shown.get('prone', 0),
                    shown.get('hotspots', 0),
                    aggrescan.HOTSPOT_THRESHOLD),
                 'The structure looks %s, so the %s score (%.4f) is the '
                 'one to read.'
                 % ('globular' if result.is_globular() else 'extended',
                    'maximum' if used == 'max' else 'average', value)]
        if shown.get('selection'):
            parts.append('Selection "%s".' % shown['selection'])
        if job.report.get('globular_dropped'):
            parts.append('%d low-confidence residues left out.'
                         % job.report['globular_dropped'])
        if job.globular and not job.report.get('looks_predicted', True):
            parts.append('The B-factors are not pLDDT, so nothing was '
                         'trimmed.')
        if job.ph_results:
            worst = max(job.ph_results,
                        key=lambda row: row[1].metric()[0])
            parts.append('Most aggregation prone at pH %.1f.' % worst[0])
        if job.mutated:
            parts.append('Mutated %s.'
                         % ', '.join('%s/%s' % (one.chain, one.label)
                                     for one in job.mutated))
        if job.mutations and not job.suggestions:
            parts.append('No mutation improved the predicted solubility.')
        self.status.setText(' '.join(parts))

    def _ph_chosen(self):
        '''
        Colour the structure by the pH the user clicked
        '''
        if self.job is None or not self.job.ph_results:
            return
        rows = self.ph_table.selectionModel().selectedRows()
        if not rows:
            return
        value, result = self.job.ph_results[rows[0].row()]
        self.result = result
        aggrescan.apply_scores(self.job.object, result, _self=self.cmd)
        self._fill(self.scores, [
            (r.chain, '%s%s' % (r.letter, r.resi), '%.4f' % r.score,
             '%.1f' % r.exposure, 'yes' if r.prone else '')
            for r in sorted(result.residues, key=lambda r: -r.score)])
        self.status.setText('Showing pH %.1f: maximum %.4f, average %.4f'
                            % (value, result.maximum, result.average))

    def load_models(self):
        '''
        Bring the dynamic mode's models into the session, coloured the
        same way
        '''
        if self.job is None or not self.job.models:
            return
        loaded = 0
        for name, path, result in self.job.models:
            if name == 'input':
                continue
            label = self.cmd.get_legal_name('%s_%s' % (self.job.object,
                                                       name))
            self.cmd.load(path, label)
            aggrescan.apply_scores(label, result, select=False,
                                   _self=self.cmd)
            loaded += 1
        self.status.setText('Loaded %d model(s)' % loaded)
        self.load_button.setEnabled(False)

    def _current_table(self):
        return self.tabs.currentWidget()

    def copy(self):
        table = self._current_table()
        lines = ['\t'.join(
            table.horizontalHeaderItem(c).text()
            for c in range(table.columnCount()))]
        for r in range(table.rowCount()):
            lines.append('\t'.join(table.item(r, c).text()
                                   for c in range(table.columnCount())))
        QtWidgets.QApplication.clipboard().setText('\n'.join(lines))
        self.status.setText('Table copied')

    def save(self):
        '''
        Everything the run produced, into a folder the user picks: the
        scores as CSV and as PDB files carrying them in the B-factors,
        one per pH tested
        '''
        if self.job is None or self.job.result is None:
            self.status.setText('Nothing to save yet')
            return
        folder = QtWidgets.QFileDialog.getExistingDirectory(
            self, 'Save the results into a folder')
        if not folder:
            return
        try:
            written = self.job.save(folder)
        except Exception as e:
            self.status.setText('Could not save: %s' % e)
            return
        self.status.setText('Wrote %d file(s) to %s'
                            % (len(written), folder))

    def closeEvent(self, event):
        if self.job is not None:
            self.job.cleanup()
        super().closeEvent(event)


def show_dialog():
    global _dialog
    if _dialog is None:
        from pymol.gui import get_qtwindow
        _dialog = AggregationDialog(get_qtwindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
