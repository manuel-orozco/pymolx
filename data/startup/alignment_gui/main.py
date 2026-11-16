'''
Superposition/Alignment dialog
'''

from pymol import cmd
from pymol.Qt import QtCore, QtWidgets

from . import command

Qt = QtCore.Qt

_dialog = None


class AlignmentDialog(QtWidgets.QDialog):

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent)
        self.cmd = _self
        self.setObjectName('alignment_dialog')
        self.setWindowTitle('Superposition/Alignment')
        self.resize(760, 560)

        layout = QtWidgets.QVBoxLayout(self)

        link = QtWidgets.QLabel(
            'Visit <a href="{0}">{0}</a> to learn about alignment methods '
            'in PyMOL'.format(command.WIKI_URL), self)
        link.setOpenExternalLinks(True)
        layout.addWidget(link)

        # mode
        modes = QtWidgets.QHBoxLayout()
        self.many_to_one = QtWidgets.QRadioButton('Many to one', self)
        self.one_to_one = QtWidgets.QRadioButton('One to one', self)
        self.many_to_one.setChecked(True)
        modes.addWidget(self.many_to_one, 1)
        modes.addWidget(self.one_to_one, 1)
        layout.addLayout(modes)

        form = QtWidgets.QGridLayout()
        form.setColumnStretch(1, 1)
        layout.addLayout(form)

        self.method = QtWidgets.QComboBox(self)
        self.method.addItems(command.available_methods(_self))
        form.addWidget(QtWidgets.QLabel('Method', self), 0, 0)
        form.addWidget(self.method, 0, 1, 1, 2)

        self.mobile = self._selection_combo()
        self.mobile_state = self._state_spin()
        form.addWidget(QtWidgets.QLabel('Mobile Selection', self), 1, 0)
        form.addWidget(self.mobile, 1, 1)
        form.addWidget(self.mobile_state, 1, 2)
        self.mobile_chains, self.mobile_residues = self._region_row(form, 2)

        self.target = self._selection_combo()
        self.target_state = self._state_spin()
        form.addWidget(QtWidgets.QLabel('Target Selection', self), 3, 0)
        form.addWidget(self.target, 3, 1)
        form.addWidget(self.target_state, 3, 2)
        self.target_chains, self.target_residues = self._region_row(form, 4)

        self.object_name = QtWidgets.QLineEdit(self)
        self.object_name.setPlaceholderText('name of alignment object')
        form.addWidget(QtWidgets.QLabel('Create alignment object', self), 5, 0)
        form.addWidget(self.object_name, 5, 1, 1, 2)

        # outlier rejection (methods with "cycles")
        self.outlier_rejection = QtWidgets.QCheckBox('Outlier rejection', self)
        self.outlier_rejection.setChecked(True)
        layout.addWidget(self.outlier_rejection)

        self.outlier_frame = QtWidgets.QFrame(self)
        self.outlier_frame.setObjectName('outlier_frame')
        self.outlier_frame.setFrameShape(QtWidgets.QFrame.Shape.StyledPanel)
        outlier = QtWidgets.QHBoxLayout(self.outlier_frame)
        self.cycles = QtWidgets.QSpinBox(self)
        self.cycles.setRange(1, 100)
        self.cycles.setValue(5)
        self.cutoff = QtWidgets.QDoubleSpinBox(self)
        self.cutoff.setRange(0.1, 100.0)
        self.cutoff.setSingleStep(0.5)
        self.cutoff.setDecimals(1)
        self.cutoff.setValue(2.0)
        outlier.addWidget(QtWidgets.QLabel('Cycles', self))
        outlier.addWidget(self.cycles)
        outlier.addStretch(1)
        outlier.addWidget(QtWidgets.QLabel('Cutoff', self))
        outlier.addWidget(self.cutoff)
        outlier.addStretch(1)
        layout.addWidget(self.outlier_frame)

        # command preview
        layout.addWidget(QtWidgets.QLabel(
            'This will run the following command', self))
        self.preview = QtWidgets.QPlainTextEdit(self)
        self.preview.setObjectName('alignment_preview')
        self.preview.setReadOnly(True)
        layout.addWidget(self.preview, 1)

        self.ok_button = QtWidgets.QPushButton('OK', self)
        self.ok_button.setDefault(True)
        self.ok_button.clicked.connect(self.run)
        layout.addWidget(self.ok_button)

        for signal in [
            self.many_to_one.toggled,
            self.method.currentIndexChanged,
            self.mobile.currentTextChanged,
            self.target.currentTextChanged,
            self.mobile_chains.textChanged,
            self.mobile_residues.textChanged,
            self.target_chains.textChanged,
            self.target_residues.textChanged,
            self.mobile_state.valueChanged,
            self.target_state.valueChanged,
            self.object_name.textChanged,
            self.outlier_rejection.toggled,
            self.cycles.valueChanged,
            self.cutoff.valueChanged,
        ]:
            signal.connect(self.update_preview)

        self.refresh()

    def _region_row(self, form, row):
        '''
        "Chains" and "Residues" under a selection: the part of it to fit
        on. Empty means all of it.
        '''
        line = QtWidgets.QHBoxLayout()
        line.setContentsMargins(0, 0, 0, 0)
        chains = QtWidgets.QLineEdit(self)
        chains.setPlaceholderText('all chains, e.g. H or H+L')
        chains.setToolTip('Fit on these chains only; the whole object '
                          'still moves')
        residues = QtWidgets.QLineEdit(self)
        residues.setPlaceholderText('all residues, e.g. 1-120 or 10-20, 44')
        residues.setToolTip('Fit on these residues only; the whole object '
                            'still moves')
        for label, widget in (('Chains', chains), ('Residues', residues)):
            tag = QtWidgets.QLabel(label, self)
            tag.setObjectName('prep_note')
            line.addWidget(tag)
            line.addWidget(widget, 1)
        form.addWidget(QtWidgets.QLabel('', self), row, 0)
        form.addLayout(line, row, 1, 1, 2)
        return chains, residues

    def _selection_combo(self):
        combo = QtWidgets.QComboBox(self)
        combo.setEditable(True)
        combo.setInsertPolicy(QtWidgets.QComboBox.InsertPolicy.NoInsert)
        combo.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                            QtWidgets.QSizePolicy.Policy.Fixed)
        return combo

    def _state_spin(self):
        spin = QtWidgets.QSpinBox(self)
        spin.setRange(-1, 99999)
        spin.setValue(-1)
        spin.setToolTip('State: -1 = current state, 0 = all states')
        return spin

    def refresh(self):
        '''
        Fill the selection lists from the session; keeps edited text.
        '''
        _self = self.cmd
        objects = _self.get_names('public_objects')
        selections = _self.get_names('public_selections')

        def fill(combo, items, default):
            text = combo.currentText()
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(items)
            combo.setEditText(text or default)
            combo.blockSignals(False)

        fill(self.mobile, ['*'] + objects + selections, '*')
        fill(self.target, objects + selections, objects[0] if objects else '')
        self.update_preview()

    def command_parts(self):
        return command.build_command(
            self.method.currentText(),
            self.mobile.currentText(),
            self.target.currentText(),
            many_to_one=self.many_to_one.isChecked(),
            mobile_state=self.mobile_state.value(),
            target_state=self.target_state.value(),
            object_name=self.object_name.text(),
            outlier_rejection=self.outlier_rejection.isChecked(),
            cycles=self.cycles.value(),
            cutoff=self.cutoff.value(),
            mobile_chains=self.mobile_chains.text(),
            mobile_residues=self.mobile_residues.text(),
            target_chains=self.target_chains.text(),
            target_residues=self.target_residues.text(),
            _self=self.cmd)

    def update_preview(self, *_):
        method = self.method.currentText()
        if not method:
            return
        outliers = command.supports_outlier_rejection(method, self.cmd)
        self.outlier_rejection.setEnabled(outliers)
        self.outlier_frame.setEnabled(
            outliers and self.outlier_rejection.isChecked())
        try:
            parts = self.command_parts()
        except ValueError as e:      # chains or residues not understood
            self.preview.setPlainText(str(e))
            self.ok_button.setEnabled(False)
            return
        self.ok_button.setEnabled(True)
        self.preview.setPlainText(command.format_command(*parts))

    def run(self):
        try:
            parts = self.command_parts()
        except ValueError as e:
            self.preview.setPlainText(str(e))
            return
        text = command.format_command(*parts, multiline=False)
        self.cmd.do(text)
        self.accept()

    def showEvent(self, event):
        self.refresh()
        super().showEvent(event)


def show_dialog():
    global _dialog
    if _dialog is None:
        from pymol.gui import get_qtwindow
        _dialog = AlignmentDialog(get_qtwindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
