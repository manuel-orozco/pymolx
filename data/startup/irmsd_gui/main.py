'''
Interface RMSD Clustering dialog: RMSD Clustering's, with the settings of
an interface RMSD (what counts as interface, the largest i-RMSD within a
cluster, the smallest cluster worth the name) in place of what to
superimpose and what to measure (pymolx.irmsd)
'''

from pymol import cmd
from pymol.Qt import QtWidgets

from pmg_tk.startup.cluster_gui.main import ClusterDialog, METHODS, _combo
from pymolx import irmsd

_dialog = None

ATOMS = [('Backbone (as CAPRI)', 'backbone'), ('C-alpha', 'ca'),
         ('Heavy atoms', 'heavy')]
SCOPES = [('Each pair of models', 'pair'), ('All models together', 'all')]


class InterfaceClusterDialog(ClusterDialog):
    '''
    Complexes clustered by the RMSD of their interface residues
    '''

    TITLE = 'Interface RMSD Clustering'
    SPREAD = 'Mean i-RMSD'
    HEATMAP_NOTE = 'Interface RMSD matrix by cluster (hover for values)'

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent, _self=_self)
        self.setObjectName('irmsd_dialog')

    def _settings_box(self):
        box = QtWidgets.QGroupBox('Settings', self)
        box.setObjectName('prep_step')
        grid = QtWidgets.QGridLayout(box)
        # (named as RMSD Clustering's: it superimposes on the one and
        # colors the other, and "detect" fills both)
        grid.addWidget(QtWidgets.QLabel('Antigen', box), 0, 0)
        self.fit = QtWidgets.QLineEdit(box)
        self.fit.setObjectName('irmsd_antigen')
        self.fit.setPlaceholderText('antigen chains (detected)')
        grid.addWidget(self.fit, 0, 1)
        grid.addWidget(QtWidgets.QLabel('Antibody', box), 1, 0)
        self.rmsd = QtWidgets.QLineEdit(box)
        self.rmsd.setObjectName('irmsd_antibody')
        self.rmsd.setPlaceholderText('antibody chains (detected)')
        grid.addWidget(self.rmsd, 1, 1)
        detect = QtWidgets.QPushButton('Detect antibody and antigen', box)
        detect.clicked.connect(self.detect)
        grid.addWidget(detect, 2, 1)

        grid.addWidget(QtWidgets.QLabel('Interface within', box), 3, 0)
        self.distance = QtWidgets.QDoubleSpinBox(box)
        self.distance.setObjectName('irmsd_distance')
        self.distance.setRange(2.0, 30.0)
        self.distance.setSingleStep(1.0)
        self.distance.setDecimals(1)
        self.distance.setValue(irmsd.INTERFACE)
        self.distance.setSuffix(' A')
        self.distance.setToolTip(
            'A residue is part of the interface if one of its heavy atoms '
            'is this close to the other side. 10 A is CAPRI\'s; 5 A keeps '
            'only the residues in direct contact')
        grid.addWidget(self.distance, 3, 1)
        grid.addWidget(QtWidgets.QLabel('Interface of', box), 4, 0)
        self.scope = _combo(box, SCOPES)
        self.scope.setObjectName('irmsd_scope')
        self.scope.setToolTip(
            'Each pair: two models are compared over the residues that '
            'are interface in either of them.\nAll models: over the '
            'residues that are interface in any model, the same for every '
            'pair. With models spread over the antigen this comes close '
            'to the RMSD of the whole complex')
        grid.addWidget(self.scope, 4, 1)
        grid.addWidget(QtWidgets.QLabel('Atoms', box), 5, 0)
        self.atoms = _combo(box, ATOMS)
        self.atoms.setObjectName('irmsd_atoms')
        grid.addWidget(self.atoms, 5, 1)
        grid.addWidget(QtWidgets.QLabel('Method', box), 6, 0)
        self.method = _combo(box, METHODS)
        grid.addWidget(self.method, 6, 1)

        self.cutoff_label = QtWidgets.QLabel('Maximum i-RMSD', box)
        grid.addWidget(self.cutoff_label, 7, 0)
        self.cutoff = QtWidgets.QDoubleSpinBox(box)
        self.cutoff.setObjectName('irmsd_cutoff')
        self.cutoff.setRange(0.1, 50.0)
        self.cutoff.setSingleStep(0.5)
        self.cutoff.setDecimals(1)
        self.cutoff.setValue(irmsd.CUTOFF)
        self.cutoff.setSuffix(' A')
        self.cutoff.setToolTip(
            'The largest interface RMSD at which models still belong to '
            'one cluster. CAPRI calls an agreement within %s' % ', '.join(
                '%g A %s' % (limit, name) for name, limit in irmsd.CAPRI))
        grid.addWidget(self.cutoff, 7, 1)
        grid.addWidget(QtWidgets.QLabel('Minimum cluster size', box), 8, 0)
        self.min_size = QtWidgets.QSpinBox(box)
        self.min_size.setObjectName('irmsd_min_size')
        self.min_size.setRange(1, 999)
        self.min_size.setValue(irmsd.MIN_SIZE)
        self.min_size.setToolTip(
            'Clusters of fewer models are no clusters: their models are '
            'kept together as "unclustered", in grey')
        grid.addWidget(self.min_size, 8, 1)

        row = self._score_rows(box, grid, 9)
        self._display_rows(box, grid, row)
        return box

    def _compute(self, models):
        return irmsd.compute(
            ' '.join(models), self.fit.text(), self.rmsd.text(),
            atoms=self.atoms.currentData(), cutoff=self.cutoff.value(),
            method=self.method.currentData(),
            distance=self.distance.value(),
            scope=self.scope.currentData(), min_size=self.min_size.value(),
            best_by=self.best_by.currentData(),
            scores=self.scores.text().strip(),
            scores_csv=self.scores_csv.text().strip(),
            higher=not self.lower_better.isChecked(), _self=self.cmd)

    def _summary(self, r):
        return irmsd.summary(r)


def show_dialog():
    global _dialog
    if _dialog is None:
        from pymol.gui import get_qtwindow
        _dialog = InterfaceClusterDialog(get_qtwindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
