'''
The toolbar's Labels tool: a small window that stays open next to the
viewer. Pick residues, choose what the label says (K593, LYS593, 593),
which atom it sits on and how large the letters are, and press Label.
'''

from pymol import cmd
from pymol.Qt import QtCore, QtWidgets

from .. import labels

Qt = QtCore.Qt

FORMATS = [('1-letter code (%s)', 'one'), ('3-letter code (%s)', 'three'),
           ('Number only (%s)', 'none')]
ALL = 'all'
ACTIVE = 'sele'         # PyMOL's selection of what was clicked
SIZES = (6, 72)         # points
DEFAULT_SIZE = 14       # PyMOL's "label_size"

_dialog = None


def targets(_self=cmd):
    '''
    What can be labeled: the selections, the molecules, then everything.

    :return: [(text shown, selection)]
    '''
    found = [('(%s)' % name, name)
             for name in _self.get_names('public_selections') or []]
    found += [(name, name) for name in _self.get_names('public_objects') or []
              if _self.get_type(name) == 'object:molecule']
    return found + [('All residues', ALL)]


class TargetBox(QtWidgets.QComboBox):
    '''
    The list is read again each time it opens: selections come and go
    '''

    about_to_open = QtCore.Signal()

    def showPopup(self):
        self.about_to_open.emit()
        super().showPopup()


class LabelDialog(QtWidgets.QDialog):

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent)
        self.cmd = _self
        # the residues the user asked for by name, while they exist;
        # None: the clicked selection if there is one, else everything
        self.chosen = None
        self.setObjectName('label_dialog')
        self.setWindowTitle('Labels')

        layout = QtWidgets.QVBoxLayout(self)
        form = QtWidgets.QFormLayout()
        layout.addLayout(form)

        self.target = TargetBox(self)
        self.target.setToolTip('The residues to label: a selection (what '
                               'you clicked in the viewer is "sele"), a '
                               'molecule or everything')
        self.target.about_to_open.connect(self.refresh)
        self.target.activated.connect(self._choose)
        form.addRow('Residues', self.target)

        self.names = QtWidgets.QComboBox(self)
        for text, value in FORMATS:
            self.names.addItem(text % labels.NAMES[value], value)
        self.names.setToolTip('Ligands, ions and waters have no 1-letter '
                              'code: they keep their name (NAG601)')
        form.addRow('Format', self.names)

        self.atom = QtWidgets.QComboBox(self)
        self.atom.setEditable(True)
        self.atom.addItems(labels.ATOMS)
        self.atom.setToolTip('The atom the label is put on; any atom name '
                             'can be typed. A residue without it has its '
                             'label on CA, or on its first atom')
        form.addRow('Atom', self.atom)

        self.size = QtWidgets.QSpinBox(self)
        self.size.setRange(*SIZES)
        self.size.setSuffix(' pt')
        self.size.setToolTip('Font size of all labels in the viewer')
        self.size.valueChanged.connect(self._set_size)
        form.addRow('Font size', self.size)

        row = QtWidgets.QHBoxLayout()
        self.clear_button = QtWidgets.QPushButton('Clear', self)
        self.clear_button.setToolTip('Take the labels off these residues')
        self.clear_button.clicked.connect(self.clear)
        row.addWidget(self.clear_button)
        row.addStretch(1)
        self.label_button = QtWidgets.QPushButton('Label', self)
        self.label_button.setDefault(True)
        self.label_button.setToolTip('Label these residues, in place of '
                                     'the labels they have')
        self.label_button.clicked.connect(self.apply)
        row.addWidget(self.label_button)
        layout.addLayout(row)

        self.status = QtWidgets.QLabel(self)
        self.status.setObjectName('prep_note')
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        self.refresh()

    # what to label

    def refresh(self):
        '''
        Read the selections and molecules again, and the font size
        '''
        found = targets(self.cmd)
        selections = [selection for _, selection in found]
        if self.chosen not in selections:
            self.chosen = None
        current = self.chosen or (ACTIVE if ACTIVE in selections else ALL)
        self.target.blockSignals(True)
        self.target.clear()
        for text, selection in found:
            self.target.addItem(text, selection)
        self.target.setCurrentIndex(selections.index(current))
        self.target.blockSignals(False)

        size = self.cmd.get_setting_float('label_size')
        self.size.blockSignals(True)
        # (a negative size is in Angstrom: left alone until one is set)
        self.size.setValue(round(size) if size > 0 else DEFAULT_SIZE)
        self.size.blockSignals(False)

    def _choose(self, row):
        self.chosen = self.target.itemData(row)

    def selection(self):
        return self.target.currentData()

    def event(self, event):
        # back from the viewer, where a selection may have been made
        if event.type() == QtCore.QEvent.Type.WindowActivate:
            self.refresh()
        return super().event(event)

    # actions

    def _set_size(self, value):
        self.cmd.set('label_size', value, quiet=1)

    def _where(self):
        return self.target.currentText()

    def apply(self):
        self.refresh()
        try:
            count = self.cmd.label_residues(
                self.selection(), self.names.currentData(),
                self.atom.currentText(), quiet=1)
        except Exception as e:
            self.status.setText(str(e))
            return
        if not count:
            self.status.setText('No residues to label in %s.' %
                                self._where())
            return
        self.status.setText('Labeled %d residue%s of %s.' % (
            count, '' if count == 1 else 's', self._where()))

    def clear(self):
        self.refresh()
        count = labels.clear(self.selection(), _self=self.cmd)
        self.status.setText(
            'Took the labels off %d residue%s of %s.' % (
                count, '' if count == 1 else 's', self._where())
            if count else 'No labels in %s.' % self._where())


def show(parent=None, _self=cmd):
    '''
    Open the Labels window (one per session), or bring it to the front.
    '''
    global _dialog
    if _dialog is None:
        _dialog = LabelDialog(parent, _self=_self)
    _dialog.refresh()
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
