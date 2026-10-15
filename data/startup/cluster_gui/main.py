'''
RMSD Clustering dialog: choose the structures, what to superimpose on and
what to measure, the method and cutoff; results as a table (click a
cluster to show it) and an RMSD heatmap ordered by cluster.
'''

from pymol import cmd
from pymol.Qt import QtCore, QtGui, QtWidgets

from pymolx import cluster
from pymolx.gui import theme

Qt = QtCore.Qt

_dialog = None

ATOMS = [('C-alpha', 'ca'), ('Backbone', 'backbone'), ('Heavy atoms', 'heavy')]
BEST_BY = [('ipTM', 'iptm'), ('pTM', 'ptm'), ('Ranking score', 'ranking'),
           ('pLDDT', 'plddt')]
METHODS = [('Hierarchical, average linkage', 'average'),
           ('Hierarchical, complete linkage', 'complete'),
           ('Hierarchical, single linkage', 'single'),
           ('GROMOS', 'gromos')]


def _combo(parent, items):
    combo = QtWidgets.QComboBox(parent)
    for label, value in items:
        combo.addItem(label, value)
    return combo


class Heatmap(QtWidgets.QWidget):
    '''
    The RMSD matrix, models ordered by cluster (representative first),
    cluster boundaries drawn; hover for the values
    '''

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName('cluster_heatmap')
        self.setMouseTracking(True)
        self.setMinimumSize(160, 160)
        self.result = None
        self.order = []
        self.top = 1.0

    def set_result(self, result):
        self.result = result
        self.order = []
        for c in result.clusters:
            members = [c['representative']] + [
                m for m in c['members'] if m != c['representative']]
            self.order += [result.models.index(m) for m in members]
        self.top = max(float(result.matrix.max()), 1e-6)
        self.update()

    def sizeHint(self):
        return QtCore.QSize(320, 320)

    def _cell(self):
        n = len(self.order)
        side = min(self.width(), self.height())
        return side / n if n else 0

    @staticmethod
    def _color(fraction):
        low, mid, high = (QtGui.QColor(theme.COLORS[k]) for k in (
            'heat_low', 'heat_mid', 'heat_high'))
        a, b, t = (low, mid, fraction * 2) if fraction < 0.5 else \
            (mid, high, fraction * 2 - 1)
        return QtGui.QColor(*(round(x + (y - x) * t) for x, y in zip(
            a.getRgb()[:3], b.getRgb()[:3])))

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.fillRect(self.rect(), QtGui.QColor(theme.COLORS['base']))
        if not self.result or not self.order:
            return
        cell = self._cell()
        matrix = self.result.matrix
        for row, i in enumerate(self.order):
            for col, j in enumerate(self.order):
                painter.fillRect(QtCore.QRectF(col * cell, row * cell,
                                               cell + 0.5, cell + 0.5),
                                 self._color(matrix[i, j] / self.top))
        pen = QtGui.QPen(QtGui.QColor(theme.COLORS['text']))
        pen.setWidthF(1.0)
        painter.setPen(pen)
        start = 0
        for c in self.result.clusters:
            painter.drawRect(QtCore.QRectF(start * cell, start * cell,
                                           c['size'] * cell, c['size'] * cell))
            start += c['size']

    def mouseMoveEvent(self, event):
        cell = self._cell()
        if not cell:
            return
        pos = event.position() if hasattr(event, 'position') else event.pos()
        row, col = int(pos.y() // cell), int(pos.x() // cell)
        if 0 <= row < len(self.order) and 0 <= col < len(self.order):
            i, j = self.order[row], self.order[col]
            models = self.result.models
            self.setToolTip('%s\n%s\n%.2f A' % (
                models[i], models[j], self.result.matrix[i, j]))


class ClusterDialog(QtWidgets.QDialog):

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent)
        self.cmd = _self
        self.result = None
        self.setObjectName('cluster_dialog')
        self.setWindowTitle('RMSD Clustering')
        self.resize(900, 660)

        layout = QtWidgets.QHBoxLayout(self)
        left = QtWidgets.QVBoxLayout()
        layout.addLayout(left, 2)
        right = QtWidgets.QVBoxLayout()
        layout.addLayout(right, 3)

        # structures
        box = QtWidgets.QGroupBox('Structures', self)
        box.setObjectName('prep_step')
        column = QtWidgets.QVBoxLayout(box)
        self.filter = QtWidgets.QLineEdit(box)
        self.filter.setPlaceholderText('filter, e.g. *seed*')
        self.filter.textChanged.connect(self._filter)
        column.addWidget(self.filter)
        self.models = QtWidgets.QListWidget(box)
        self.models.setObjectName('cluster_models')
        column.addWidget(self.models, 1)
        buttons = QtWidgets.QHBoxLayout()
        for label, state in (('All', True), ('None', False)):
            button = QtWidgets.QPushButton(label, box)
            button.clicked.connect(lambda _=False, s=state: self._check(s))
            buttons.addWidget(button)
        refresh = QtWidgets.QPushButton('Refresh', box)
        refresh.clicked.connect(self.refresh)
        buttons.addWidget(refresh)
        column.addLayout(buttons)
        left.addWidget(box, 1)

        # settings
        box = QtWidgets.QGroupBox('Settings', self)
        box.setObjectName('prep_step')
        grid = QtWidgets.QGridLayout(box)
        grid.addWidget(QtWidgets.QLabel('Superimpose on', box), 0, 0)
        self.fit = QtWidgets.QLineEdit(box)
        self.fit.setPlaceholderText('antigen chains (detected)')
        grid.addWidget(self.fit, 0, 1)
        grid.addWidget(QtWidgets.QLabel('RMSD of', box), 1, 0)
        self.rmsd = QtWidgets.QLineEdit(box)
        self.rmsd.setPlaceholderText('antibody chains (detected)')
        grid.addWidget(self.rmsd, 1, 1)
        detect = QtWidgets.QPushButton('Detect antibody and antigen', box)
        detect.clicked.connect(self.detect)
        grid.addWidget(detect, 2, 1)
        grid.addWidget(QtWidgets.QLabel('Atoms', box), 3, 0)
        self.atoms = _combo(box, ATOMS)
        grid.addWidget(self.atoms, 3, 1)
        grid.addWidget(QtWidgets.QLabel('Method', box), 4, 0)
        self.method = _combo(box, METHODS)
        grid.addWidget(self.method, 4, 1)
        grid.addWidget(QtWidgets.QLabel('RMSD cutoff', box), 5, 0)
        self.cutoff = QtWidgets.QDoubleSpinBox(box)
        self.cutoff.setRange(0.1, 50.0)
        self.cutoff.setSingleStep(0.5)
        self.cutoff.setDecimals(1)
        self.cutoff.setValue(3.0)
        self.cutoff.setSuffix(' A')
        grid.addWidget(self.cutoff, 5, 1)
        grid.addWidget(QtWidgets.QLabel('Best model by', box), 6, 0)
        self.best_by = _combo(box, BEST_BY)
        self.best_by.setToolTip(
            'Scores from the predictor\'s JSON next to each model file '
            '(AlphaFold 3, Boltz, ColabFold); pLDDT from the B-factors')
        grid.addWidget(self.best_by, 6, 1)
        grid.addWidget(QtWidgets.QLabel('Scores folder', box), 7, 0)
        scores = QtWidgets.QHBoxLayout()
        self.scores = QtWidgets.QLineEdit(box)
        self.scores.setPlaceholderText('next to the loaded files')
        self.scores.setToolTip('A folder to search for score files by '
                               'object name (for models loaded before, or '
                               'from elsewhere)')
        scores.addWidget(self.scores, 1)
        browse = QtWidgets.QPushButton('...', box)
        browse.setMaximumWidth(32)
        browse.clicked.connect(self._browse_scores)
        scores.addWidget(browse)
        grid.addLayout(scores, 7, 1)
        grid.addWidget(QtWidgets.QLabel('Scores CSV', box), 8, 0)
        scores_csv = QtWidgets.QHBoxLayout()
        self.scores_csv = QtWidgets.QLineEdit(box)
        self.scores_csv.setPlaceholderText('e.g. <job>_ranking_scores.csv')
        self.scores_csv.setToolTip(
            'A row per model: a model/name/file column, or seed and sample '
            '(AlphaFold 3); score columns. A pattern (*) reads several '
            'files. Its scores come first.')
        self.scores_csv.editingFinished.connect(self._csv_changed)
        scores_csv.addWidget(self.scores_csv, 1)
        browse = QtWidgets.QPushButton('...', box)
        browse.setMaximumWidth(32)
        browse.clicked.connect(self._browse_csv)
        scores_csv.addWidget(browse)
        grid.addLayout(scores_csv, 8, 1)
        self.lower = QtWidgets.QCheckBox('Lower score is better', box)
        self.lower.setToolTip('Set from the score\'s name (e.g. a rank or an '
                              'energy); change it if needed')
        grid.addWidget(self.lower, 9, 1)
        self.best_by.currentIndexChanged.connect(self._best_by_changed)
        self.superpose = QtWidgets.QCheckBox(
            'Superimpose all on the largest cluster', box)
        self.group = QtWidgets.QCheckBox('Group by cluster', box)
        self.color = QtWidgets.QCheckBox('Color by cluster', box)
        for row, check in enumerate((self.superpose, self.group,
                                     self.color), 10):
            check.setChecked(True)
            grid.addWidget(check, row, 0, 1, 2)
        left.addWidget(box)

        run = QtWidgets.QHBoxLayout()
        self.status = QtWidgets.QLabel(self)
        self.status.setObjectName('prep_note')
        self.status.setWordWrap(True)
        run.addWidget(self.status, 1)
        self.run_button = QtWidgets.QPushButton('Cluster', self)
        self.run_button.setDefault(True)
        self.run_button.clicked.connect(self.run)
        run.addWidget(self.run_button)
        left.addLayout(run)

        # results
        self.table = QtWidgets.QTableWidget(0, 7, self)
        self.table.setObjectName('cluster_table')
        self.table.setHorizontalHeaderLabels(
            ['Cluster', 'Size', 'Mean RMSD', 'pLDDT', 'Best', 'Score',
             'Representative'])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self._show_cluster)
        right.addWidget(self.table, 2)
        actions = QtWidgets.QHBoxLayout()
        show_all = QtWidgets.QPushButton('Show All', self)
        show_all.clicked.connect(self.show_all)
        actions.addWidget(show_all)
        show_best = QtWidgets.QPushButton('Show Best Models', self)
        show_best.setToolTip('The best model of each cluster only')
        show_best.clicked.connect(self.show_best)
        actions.addWidget(show_best)
        actions.addStretch(1)
        self.save_button = QtWidgets.QPushButton('Save CSV...', self)
        self.save_button.setEnabled(False)
        self.save_button.clicked.connect(self.save_csv)
        actions.addWidget(self.save_button)
        right.addLayout(actions)
        self.heatmap = Heatmap(self)
        right.addWidget(self.heatmap, 3)
        note = QtWidgets.QLabel('RMSD matrix by cluster (hover for values)',
                                self)
        note.setObjectName('prep_note')
        right.addWidget(note)

    # structures

    def refresh(self):
        checked = set(self.checked_models()) if self.models.count() else None
        self.models.clear()
        names = [n for n in (self.cmd.get_object_list('all') or [])
                 if self.cmd.get_type(n) == 'object:molecule' and
                 not n.startswith('_')]
        for name in names:
            item = QtWidgets.QListWidgetItem(name, self.models)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked if checked is None or
                               name in checked else Qt.CheckState.Unchecked)
        self._filter(self.filter.text())

    def _filter(self, text):
        import fnmatch
        pattern = text.strip()
        if pattern and not any(ch in pattern for ch in '*?['):
            pattern = '*%s*' % pattern
        for row in range(self.models.count()):
            item = self.models.item(row)
            item.setHidden(bool(pattern) and not fnmatch.fnmatch(
                item.text(), pattern))

    def _check(self, state):
        for row in range(self.models.count()):
            item = self.models.item(row)
            if not item.isHidden():
                item.setCheckState(Qt.CheckState.Checked if state else
                                   Qt.CheckState.Unchecked)

    def checked_models(self):
        return [self.models.item(row).text()
                for row in range(self.models.count())
                if not self.models.item(row).isHidden() and
                self.models.item(row).checkState() == Qt.CheckState.Checked]

    def detect(self):
        models = self.checked_models()
        if not models:
            self.status.setText('No structures chosen')
            return
        fit, rmsd = cluster.default_selections(models[0], _self=self.cmd)
        self.fit.setText(fit)
        self.rmsd.setText(rmsd)
        self.status.setText('Antibody: %s; antigen: %s' % (rmsd, fit)
                            if fit != rmsd else
                            'No antibody chains found: all polymer atoms')

    def showEvent(self, event):
        self.refresh()
        super().showEvent(event)

    # running

    def run(self):
        models = self.checked_models()
        if len(models) < 2:
            self.status.setText('Choose at least two structures')
            return
        if not self.fit.text().strip() or not self.rmsd.text().strip():
            self.detect()
        QtWidgets.QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            self.result = cluster.compute(
                ' '.join(models), self.fit.text(), self.rmsd.text(),
                self.atoms.currentData(), self.cutoff.value(),
                self.method.currentData(),
                best_by=self.best_by.currentData(),
                scores=self.scores.text().strip(),
                scores_csv=self.scores_csv.text().strip(),
                higher=not self.lower.isChecked(), _self=self.cmd)
            cluster.apply(self.result, self.superpose.isChecked(),
                          self.group.isChecked(), self.color.isChecked(),
                          _self=self.cmd)
        except Exception as e:
            self.status.setText('Clustering failed: %s' % e)
            return
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        r = self.result
        self.status.setText(
            '%d structures, %d clusters; fit on %d atoms, RMSD of %d%s' % (
                len(r.models), len(r.clusters), r.fit_atoms, r.rmsd_atoms,
                '; CSV rows without a model: %d' % len(r.unmatched_scores)
                if r.unmatched_scores else ''))
        print(r.table())
        self._fill_table()
        self.heatmap.set_result(r)
        self.save_button.setEnabled(True)

    def _fill_table(self):
        clusters = self.result.clusters
        self.table.blockSignals(True)
        self.table.setRowCount(len(clusters))
        for row, c in enumerate(clusters):
            short = cluster.short_names(self.result.models)
            values = [str(c['cluster']), str(c['size']),
                      '%.2f A' % c['mean_rmsd'],
                      '%.1f' % c['score'] if c['score'] is not None else '-',
                      short[c['best']] if c['best'] else '-',
                      '%s %s' % (cluster.score_label(c['best_by']),
                                 cluster.format_value(c['best_score']))
                      if c['best'] else '-',
                      short[c['representative']]]
            tips = {4: c['best'], 6: c['representative']}
            for col, value in enumerate(values):
                item = QtWidgets.QTableWidgetItem(value)
                item.setToolTip(tips.get(col) or '\n'.join(c['members']))
                self.table.setItem(row, col, item)
        self.table.resizeColumnsToContents()
        self.table.blockSignals(False)

    def _show_cluster(self):
        rows = {i.row() for i in self.table.selectedIndexes()}
        if not self.result or not rows:
            return
        shown = [m for row in rows
                 for m in self.result.clusters[row]['members']]
        self.cmd.disable(' '.join(self.result.models))
        self.cmd.enable(' '.join(shown), parents=1)

    def show_best(self):
        if not self.result:
            return
        best = [c['best'] for c in self.result.clusters if c['best']]
        if not best:
            self.status.setText('No scores found: see "Scores folder"')
            return
        self.table.blockSignals(True)
        self.table.clearSelection()
        self.table.blockSignals(False)
        self.cmd.disable(' '.join(self.result.models))
        self.cmd.enable(' '.join(best), parents=1)

    def _browse_csv(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Scores CSV', '', 'CSV (*.csv);;All files (*)')
        if path:
            self.scores_csv.setText(path)
            self._csv_changed()

    def _csv_changed(self):
        '''
        The CSV's score columns in "Best model by"; its ranking score
        chosen if it has one
        '''
        text = self.scores_csv.text().strip()
        columns = []
        if text:
            try:
                table = cluster.read_scores_csvs(text, _self=self.cmd)
            except Exception as e:
                self.status.setText('Cannot read the CSV: %s' % e)
                return
            columns = sorted({k for v in table.values() for k in v})
        current = self.best_by.currentData()
        self.best_by.blockSignals(True)
        self.best_by.clear()
        for label, value in BEST_BY:
            self.best_by.addItem(label, value)
        for column in columns:
            if column not in cluster.SCORES:
                self.best_by.addItem(column, column)
        choice = 'ranking' if 'ranking' in columns else current
        index = self.best_by.findData(choice)
        self.best_by.setCurrentIndex(max(index, 0))
        self.best_by.blockSignals(False)
        self._best_by_changed()
        if columns:
            self.status.setText('CSV scores: %s' % ', '.join(
                cluster.score_label(c) for c in columns))

    def _best_by_changed(self, *_):
        self.lower.setChecked(not cluster.higher_is_better(
            self.best_by.currentData()))

    def _browse_scores(self):
        folder = QtWidgets.QFileDialog.getExistingDirectory(
            self, 'Folder with score files')
        if folder:
            self.scores.setText(folder)

    def show_all(self):
        if self.result:
            self.table.clearSelection()
            self.cmd.enable(' '.join(self.result.models), parents=1)

    def save_csv(self):
        if not self.result:
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save clusters', 'clusters.csv', 'CSV (*.csv)')
        if path:
            self.result.save_csv(path)
            self.status.setText('Saved %s' % path)


def show_dialog():
    global _dialog
    if _dialog is None:
        from pymol.gui import get_qtwindow
        _dialog = ClusterDialog(get_qtwindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
