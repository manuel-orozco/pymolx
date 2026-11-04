'''
Structure Prediction dialog: a FASTA file's sequences shown chain by
chain with their numbering; clicking a residue offers the modifications,
glycans and disulfide bonds that suit it; ligands by CCD code; the whole
job written as the input file of a prediction method (pymolx.prediction).
'''

import html
import os

from pymol import cmd
from pymol.Qt import QtCore, QtGui, QtWidgets

from pymolx import prediction
from pymolx.gui import theme
from pymolx.gui.sequence_viewer import PluginPrefs

Qt = QtCore.Qt

_dialog = None

ROW = 60            # residues of a sequence row
BLOCK = 10          # ... in blocks of this many
MARKS = {'ptm': 'prediction_ptm', 'glycan': 'prediction_glycan',
         'disulfide': 'prediction_disulfide'}


def _colour(key):
    return theme.COLORS[key]


def _title(chain):
    '''
    A chain's title, unless it only names the chain ("Chain H" on H)
    '''
    title = (chain.title or '').strip()
    if prediction.title_chain_id(title) == chain.id and \
            len(title.split()) <= 2:
        return ''
    return title


class SequenceView(QtWidgets.QTextBrowser):
    '''
    The chains' sequences, numbered, each residue a link; what is on a
    residue colours it. Emits clicked(chain, residue number).
    '''
    clicked = QtCore.Signal(str, int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName('prediction_sequence')
        self.setOpenLinks(False)
        self.setFont(theme.console_font())
        self.anchorClicked.connect(self._clicked)
        self.job = None

    def _clicked(self, url):
        # (a fragment: Qt reads "A:12" as a URL scheme and lower-cases it)
        chain, _, number = url.fragment().partition(':')
        if number.isdigit():
            self.clicked.emit(chain, int(number))

    def set_job(self, job):
        self.job = job
        bar = self.verticalScrollBar().value()
        self.setHtml(self._html())
        self.verticalScrollBar().setValue(bar)

    def marks(self, chain):
        '''
        {residue number: (kind, what it is)}
        '''
        found = {}
        for number in chain.sequons():
            found[number] = ('sequon', 'N-X-S/T site')
        for number, code in chain.modifications.items():
            found[number] = ('ptm', prediction.ptm(code)[0])
        for number, name in chain.glycans.items():
            found[number] = ('glycan', name)
        for a, b in self.job.disulfides:
            for one, other in ((a, b), (b, a)):
                if one[0] == chain.id:
                    found[one[1]] = ('disulfide', 'disulfide to %s %d'
                                     % other)
        return found

    def _residue(self, chain, number, marks):
        letter = chain.sequence[number - 1]
        kind, what = marks.get(number, ('', ''))
        # (a link's own colour would make every residue look marked)
        style = 'color:%s;' % _colour('text_bright')
        if kind == 'sequon':
            style = 'color:%s;' % _colour('prediction_sequon')
        elif kind:
            style = 'background-color:%s;color:%s;' % (
                _colour(MARKS[kind]), _colour('base'))
        title = '%s %d %s%s' % (chain.id, number, letter,
                                ': ' + what if what else '')
        return ('<a href="#%s:%d" title="%s" style="text-decoration:none;'
                '%s">%s</a>' % (chain.id, number, html.escape(title), style,
                                letter))

    def _html(self):
        if self.job is None or not self.job.chains:
            return ('<p style="color:%s">Open a FASTA file to see its '
                    'sequences here.</p>' % _colour('text_dim'))
        parts = ['<style>a{text-decoration:none}</style>']
        for chain in self.job.chains:
            marks = self.marks(chain)
            parts.append(
                '<p style="color:%s">Chain %s&nbsp;&nbsp;%s<span '
                'style="color:%s">&nbsp;&nbsp;%d residues%s</span></p>'
                % (_colour('text'), chain.id,
                   html.escape(_title(chain)), _colour('text_dim'),
                   len(chain.sequence),
                   ', %d copies' % chain.copies if chain.copies > 1 else ''))
            lines = []
            for start in range(1, len(chain.sequence) + 1, ROW):
                blocks = []
                for at in range(start, min(start + ROW,
                                           len(chain.sequence) + 1), BLOCK):
                    blocks.append(''.join(
                        self._residue(chain, number, marks)
                        for number in range(at, min(at + BLOCK,
                                                    len(chain.sequence) + 1))
                    ))
                lines.append('<span style="color:%s">%5d</span> %s' % (
                    _colour('text_dim'), start, '&nbsp;'.join(blocks)))
            parts.append('<pre style="color:%s;margin:0 0 10px 0">%s</pre>'
                         % (_colour('text_bright'), '\n'.join(lines)))
        return '\n'.join(parts)


class SequenceTextDialog(QtWidgets.QDialog):
    '''
    Sequences typed or pasted in, as FASTA; what it reads is shown as it
    is typed
    '''

    def __init__(self, text='', parent=None):
        super().__init__(parent)
        self.setObjectName('prediction_text_dialog')
        self.setWindowTitle('Type or Paste Sequences')
        self.resize(640, 460)
        layout = QtWidgets.QVBoxLayout(self)
        note = QtWidgets.QLabel(
            'One title line per chain, then its sequence:\n'
            '>Chain H\nEVQLVESGGG...\n\n'
            'A title that names its chain ("Chain H", or just "H") gives '
            'that chain its ID; the rest are lettered in order. Line '
            'breaks, gaps and lower case are fine.', self)
        note.setObjectName('prep_note')
        note.setWordWrap(True)
        layout.addWidget(note)
        self.text = QtWidgets.QPlainTextEdit(self)
        self.text.setObjectName('prediction_text')
        self.text.setFont(theme.console_font())
        self.text.setPlainText(text)
        self.text.setPlaceholderText('>Chain H\nEVQLVESGGGLVQPGG...')
        self.text.textChanged.connect(self._changed)
        layout.addWidget(self.text, 1)
        self.status = QtWidgets.QLabel(self)
        self.status.setObjectName('prep_note')
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        self.buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok |
            QtWidgets.QDialogButtonBox.StandardButton.Cancel, self)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)
        self.found = []
        self._changed()

    def _changed(self):
        try:
            self.found = prediction.read_fasta(self.text.toPlainText())
        except ValueError as e:
            self.found = []
            self.status.setText(str(e)
                                if self.text.toPlainText().strip() else '')
        else:
            job = prediction.job_from_sequences(self.found)
            self.status.setText('; '.join(
                '%s: %d residues%s' % (
                    c.id, len(c.sequence),
                    '  (%s)' % c.title if c.title else '')
                for c in job.chains))
        self.buttons.button(
            QtWidgets.QDialogButtonBox.StandardButton.Ok).setEnabled(
                bool(self.found))


class SubmitDialog(QtWidgets.QDialog):
    '''
    How to run the prediction: the job script, the input file's name, and
    whether it goes through Slurm. It shows the command it will run.
    '''

    SCRIPT_PREF = 'pymolx_af3_script_%s'

    def __init__(self, input_name='', parent=None, prefs=None):
        super().__init__(parent)
        self.prefs = prefs if prefs is not None else PluginPrefs()
        self.setObjectName('prediction_submit_dialog')
        self.setWindowTitle('Submit Job')
        self.setMinimumWidth(620)
        layout = QtWidgets.QVBoxLayout(self)
        self.modes = QtWidgets.QButtonGroup(self)
        row = QtWidgets.QHBoxLayout()
        for k, (label, value) in enumerate(prediction.SUBMIT_MODES):
            button = QtWidgets.QRadioButton(label, self)
            button.setProperty('mode', value)
            self.modes.addButton(button, k)
            row.addWidget(button)
        row.addStretch(1)
        self.modes.button(0).setChecked(True)
        self.modes.idToggled.connect(self._mode_changed)
        layout.addLayout(row)

        box = QtWidgets.QGroupBox('What to run', self)
        box.setObjectName('prep_step')
        grid = QtWidgets.QGridLayout(box)
        grid.addWidget(QtWidgets.QLabel('Job script', box), 0, 0)
        self.script = QtWidgets.QLineEdit(box)
        self.script.setObjectName('prediction_script')
        self.script.setPlaceholderText(
            'e.g. .../alphafold3/scripts/run_alphafold3.sh')
        self.script.textChanged.connect(self._refresh)
        grid.addWidget(self.script, 0, 1)
        browse = QtWidgets.QPushButton('Browse...', box)
        browse.clicked.connect(self._browse)
        grid.addWidget(browse, 0, 2)
        grid.addWidget(QtWidgets.QLabel('Input file', box), 1, 0)
        self.input_name = QtWidgets.QLineEdit(input_name, box)
        self.input_name.setObjectName('prediction_input_name')
        self.input_name.setPlaceholderText('e.g. Ipilimumab_Fab_af3.json')
        self.input_name.setToolTip('Its name only: the script looks for it '
                                   'in its own input folder')
        self.input_name.textChanged.connect(self._refresh)
        grid.addWidget(self.input_name, 1, 1, 1, 2)
        self.where = QtWidgets.QLabel(box)
        self.where.setObjectName('prep_note')
        self.where.setWordWrap(True)
        grid.addWidget(self.where, 2, 0, 1, 2)
        self.copy_button = QtWidgets.QPushButton('Copy it there', box)
        self.copy_button.clicked.connect(self._copy)
        grid.addWidget(self.copy_button, 2, 2)
        grid.setColumnStretch(1, 1)
        layout.addWidget(box)

        layout.addWidget(QtWidgets.QLabel('The command:', self))
        self.command = QtWidgets.QPlainTextEdit(self)
        self.command.setObjectName('prediction_command')
        self.command.setReadOnly(True)
        self.command.setFont(theme.console_font())
        self.command.setMaximumHeight(70)
        layout.addWidget(self.command)
        self.buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel, self)
        self.submit = self.buttons.addButton(
            'Submit', QtWidgets.QDialogButtonBox.ButtonRole.AcceptRole)
        self.submit.clicked.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)
        self.saved = ''              # the file just written, to copy
        self._mode_changed()

    def mode(self):
        return self.modes.checkedButton().property('mode')

    def _mode_changed(self, *_):
        # each way of running it remembers its own script (they are
        # different scripts), but what is typed in is not thrown away
        remembered = self.prefs.get(self.SCRIPT_PREF % self.mode(), '')
        if remembered:
            self.script.setText(remembered)
        self._refresh()

    def _browse(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'The job script', os.path.dirname(self.script.text()),
            'Shell scripts (*.sh);;All files (*)')
        if path:
            self.script.setText(path)

    def folder(self):
        return prediction.script_input_dir(self.script.text().strip())

    def _copy(self):
        import shutil
        folder = self.folder()
        if self.saved and folder:
            shutil.copy(self.saved, os.path.join(
                folder, self.input_name.text().strip()))
            self._refresh()

    def _refresh(self, *_):
        folder = self.folder()
        name = self.input_name.text().strip()
        there = bool(folder and name and
                     os.path.exists(os.path.join(folder, name)))
        if not folder:
            self.where.setText('(the script does not say where its input '
                               'folder is: give the name it expects)'
                               if self.script.text().strip() else '')
        elif there:
            self.where.setText('Found in %s' % folder)
        else:
            self.where.setText('Not in %s' % folder)
        self.copy_button.setEnabled(bool(self.saved and folder and name and
                                         not there))
        try:
            command = prediction.submit_command(
                self.script.text().strip(), name, self.mode())
        except ValueError as e:
            self.command.setPlainText(str(e))
            self.submit.setEnabled(False)
            return
        self.command.setPlainText(' '.join(command))
        self.submit.setEnabled(True)

    def accept(self):
        self.prefs.set(self.SCRIPT_PREF % self.mode(),
                       self.script.text().strip())
        super().accept()


class PredictionDialog(QtWidgets.QDialog):

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent)
        self.cmd = _self
        self.job = prediction.Job()
        self.saved = ''            # the input file last written
        self.process = None        # a prediction started from here
        self.typed = ''            # the sequences typed in, to come back to
        self.prefs = PluginPrefs()  # (tests swap in MemoryPrefs)
        self.message = ''          # what the job last said, kept by refresh
        self.pending = None        # the first cysteine of a disulfide
        self.setObjectName('prediction_dialog')
        self.setWindowTitle('Structure Prediction')
        self.resize(1000, 700)

        layout = QtWidgets.QVBoxLayout(self)
        # the job
        box = QtWidgets.QGroupBox('Sequences', self)
        box.setObjectName('prep_step')
        row = QtWidgets.QHBoxLayout(box)
        open_button = QtWidgets.QPushButton('Open FASTA...', box)
        open_button.setObjectName('prediction_open')
        open_button.clicked.connect(self.open_fasta)
        row.addWidget(open_button)
        type_button = QtWidgets.QPushButton('Type or Paste...', box)
        type_button.setObjectName('prediction_type')
        type_button.setToolTip('Write the sequences here instead, as '
                               '">Chain H" and the sequence')
        type_button.clicked.connect(self.type_sequences)
        row.addWidget(type_button)
        self.file_label = QtWidgets.QLabel('no file', box)
        self.file_label.setObjectName('prep_note')
        row.addWidget(self.file_label, 1)
        row.addWidget(QtWidgets.QLabel('Copies of', box))
        self.copy_chain = QtWidgets.QComboBox(box)
        self.copy_chain.currentIndexChanged.connect(self._show_copies)
        row.addWidget(self.copy_chain)
        self.copies = QtWidgets.QSpinBox(box)
        self.copies.setRange(1, 26)
        self.copies.setToolTip('How many copies of this chain the '
                               'prediction has')
        self.copies.valueChanged.connect(self._set_copies)
        row.addWidget(self.copies)
        layout.addWidget(box)

        note = QtWidgets.QLabel(
            'Click a residue to add a modification, a glycan or a '
            'disulfide bond. Residues of an N-X-S/T site are blue.', self)
        note.setObjectName('prep_note')
        layout.addWidget(note)
        self.view = SequenceView(self)
        self.view.clicked.connect(self.residue_clicked)
        layout.addWidget(self.view, 3)

        middle = QtWidgets.QHBoxLayout()
        layout.addLayout(middle, 2)
        # what has been added
        box = QtWidgets.QGroupBox('On the sequences', self)
        box.setObjectName('prep_step')
        column = QtWidgets.QVBoxLayout(box)
        self.table = QtWidgets.QTableWidget(0, 3, box)
        self.table.setObjectName('prediction_table')
        self.table.setHorizontalHeaderLabels(['Residue', 'What', 'CCD'])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.horizontalHeader().setStretchLastSection(True)
        column.addWidget(self.table, 1)
        remove = QtWidgets.QPushButton('Remove', box)
        remove.clicked.connect(self.remove_selected)
        column.addWidget(remove)
        middle.addWidget(box, 2)

        # ligands
        box = QtWidgets.QGroupBox('Ligands', self)
        box.setObjectName('prep_step')
        column = QtWidgets.QVBoxLayout(box)
        row = QtWidgets.QHBoxLayout()
        self.ligand = QtWidgets.QLineEdit(box)
        self.ligand.setPlaceholderText('CCD code, e.g. ATP')
        self.ligand.setMaxLength(5)
        self.ligand.returnPressed.connect(self.add_ligand)
        self.ligand.setToolTip('The code of the small molecule in the PDB '
                               'Chemical Component Dictionary')
        row.addWidget(self.ligand, 1)
        self.ligand_copies = QtWidgets.QSpinBox(box)
        self.ligand_copies.setRange(1, 20)
        row.addWidget(self.ligand_copies)
        add = QtWidgets.QPushButton('Add', box)
        add.clicked.connect(self.add_ligand)
        row.addWidget(add)
        column.addLayout(row)
        self.ligand_list = QtWidgets.QListWidget(box)
        self.ligand_list.setObjectName('prediction_ligands')
        column.addWidget(self.ligand_list, 1)
        drop = QtWidgets.QPushButton('Remove', box)
        drop.clicked.connect(self.remove_ligand)
        column.addWidget(drop)
        middle.addWidget(box, 1)

        # the method
        box = QtWidgets.QGroupBox('Prediction', self)
        box.setObjectName('prep_step')
        grid = QtWidgets.QGridLayout(box)
        self.name = QtWidgets.QLineEdit('prediction', box)
        grid.addWidget(QtWidgets.QLabel('Job name', box), 0, 0)
        grid.addWidget(self.name, 0, 1)
        self.method = QtWidgets.QComboBox(box)
        for label, value in prediction.METHODS:
            self.method.addItem(label, value)
        self.method.currentIndexChanged.connect(self.refresh)
        grid.addWidget(QtWidgets.QLabel('Method', box), 1, 0)
        grid.addWidget(self.method, 1, 1)
        self.seeds = QtWidgets.QLineEdit('1', box)
        self.seeds.setToolTip('Model seeds, separated by commas (one '
                              'prediction each)')
        self.seeds.textChanged.connect(self.refresh)
        self.seed_label = QtWidgets.QLabel('Seeds', box)
        grid.addWidget(self.seed_label, 2, 0)
        grid.addWidget(self.seeds, 2, 1)
        self.status = QtWidgets.QLabel(box)
        self.status.setObjectName('prep_note')
        self.status.setWordWrap(True)
        grid.addWidget(self.status, 3, 0, 1, 2)
        grid.setRowStretch(3, 1)
        self.save_button = QtWidgets.QPushButton('Save Input File...', box)
        self.save_button.setDefault(True)
        self.save_button.clicked.connect(self.save)
        grid.addWidget(self.save_button, 4, 0, 1, 2)
        self.submit_button = QtWidgets.QPushButton('Submit Job...', box)
        self.submit_button.setObjectName('prediction_submit')
        self.submit_button.clicked.connect(self.submit_job)
        grid.addWidget(self.submit_button, 5, 0, 1, 2)
        middle.addWidget(box, 1)

        self.refresh()

    # the sequences

    def open_fasta(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open sequences', '',
            'FASTA (*.fasta *.fa *.faa *.seq *.txt);;All files (*)')
        if not path:
            return
        try:
            self.job = prediction.job_from_fasta(path)
        except Exception as e:
            self.status.setText('Could not read it: %s' % e)
            return
        self.pending = None
        self.typed = ''
        self.file_label.setText(os.path.basename(path))
        self.name.setText(self.job.name)
        self.refresh()

    def type_sequences(self):
        dialog = SequenceTextDialog(self.typed, self)
        if dialog.exec() != QtWidgets.QDialog.DialogCode.Accepted:
            return
        self.typed = dialog.text.toPlainText()
        self.job = prediction.job_from_sequences(
            dialog.found, self.name.text().strip() or 'prediction')
        self.pending = None
        self.file_label.setText('%d sequence(s) typed in'
                                % len(self.job.chains))
        self.refresh()

    def _show_copies(self, *_):
        chain_id = self.copy_chain.currentData()
        if chain_id:
            self.copies.blockSignals(True)
            self.copies.setValue(self.job.chain(chain_id).copies)
            self.copies.blockSignals(False)

    def _set_copies(self, value):
        chain_id = self.copy_chain.currentData()
        if chain_id:
            self.job.chain(chain_id).copies = value
            self.refresh()

    # residues

    def residue_clicked(self, chain_id, number):
        self.residue_menu(chain_id, number).exec(QtGui.QCursor.pos())

    def residue_menu(self, chain_id, number):
        '''
        What can be put on a residue, by what it is
        '''
        chain = self.job.chain(chain_id)
        letter = chain.letter(number)
        menu = QtWidgets.QMenu(self)
        menu.addAction('%s %d %s' % (chain_id, number, letter)).setEnabled(
            False)
        menu.addSeparator()
        if number in chain.modifications or number in chain.glycans:
            action = menu.addAction('Remove what is here')
            action.triggered.connect(
                lambda _=False: self._clear(chain, number))
            menu.addSeparator()
        for label, code in prediction.ptms_for(letter):
            action = menu.addAction('%s  (%s)' % (label, code))
            action.triggered.connect(
                lambda _=False, c=code: self._set_ptm(chain, number, c))
        kind = {'N': 'n', 'S': 'o', 'T': 'o'}.get(letter)
        if kind:
            glycans = menu.addMenu('%s-glycan' % kind.upper())
            for name in prediction.glycan_names(kind):
                action = glycans.addAction(name)
                action.triggered.connect(
                    lambda _=False, n=name: self._set_glycan(chain, number,
                                                             n))
        if letter == 'C':
            menu.addSeparator()
            if self.pending and self.pending != (chain_id, number):
                action = menu.addAction('Disulfide bond to %s %d'
                                        % self.pending)
                action.triggered.connect(
                    lambda _=False: self._add_disulfide((chain_id, number)))
            elif self.pending == (chain_id, number):
                action = menu.addAction('Cancel this disulfide bond')
                action.triggered.connect(lambda _=False: self._pend(None))
            else:
                action = menu.addAction('Disulfide bond from here...')
                action.triggered.connect(
                    lambda _=False: self._pend((chain_id, number)))
        return menu

    def _clear(self, chain, number):
        chain.modifications.pop(number, None)
        chain.glycans.pop(number, None)
        self.job.disulfides = [
            pair for pair in self.job.disulfides
            if (chain.id, number) not in pair]
        self.refresh()

    def _set_ptm(self, chain, number, code):
        self._clear(chain, number)
        chain.modifications[number] = code
        self.refresh()

    def _set_glycan(self, chain, number, name):
        self._clear(chain, number)
        chain.glycans[number] = name
        self.refresh()

    def _pend(self, where):
        self.pending = where
        self.refresh()

    def _add_disulfide(self, other):
        pair = (self.pending, other)
        self.pending = None
        if pair not in self.job.disulfides:
            self.job.disulfides.append(pair)
        self.refresh()

    def remove_selected(self):
        rows = {i.row() for i in self.table.selectedIndexes()}
        for row in sorted(rows, reverse=True):
            what = self.table.item(row, 0).data(Qt.ItemDataRole.UserRole)
            if what[0] == 'disulfide':
                self.job.disulfides.remove(what[1])
            else:
                chain = self.job.chain(what[1])
                self._clear(chain, what[2])
        self.refresh()

    # ligands

    def add_ligand(self):
        code = self.ligand.text().strip().upper()
        if not code:
            return
        self.job.ligands.append((code, self.ligand_copies.value()))
        self.ligand.clear()
        self.refresh()

    def remove_ligand(self):
        for item in self.ligand_list.selectedItems():
            self.job.ligands.pop(self.ligand_list.row(item))
        self.refresh()

    # showing it all

    def _fill_table(self):
        rows = []
        for chain in self.job.chains:
            for number, code in sorted(chain.modifications.items()):
                rows.append((('ptm', chain.id, number),
                             '%s %d %s' % (chain.id, number,
                                           chain.letter(number)),
                             prediction.ptm(code)[0], code))
            for number, name in sorted(chain.glycans.items()):
                codes, _, _ = prediction._sugars(chain, number, name)
                rows.append((('glycan', chain.id, number),
                             '%s %d %s' % (chain.id, number,
                                           chain.letter(number)),
                             '%s-glycan %s' % (
                                 prediction.glycan(name)[0].upper(), name),
                             ' '.join(codes)))
        for pair in self.job.disulfides:
            rows.append((('disulfide', pair),
                         '%s %d - %s %d' % (pair[0] + pair[1]),
                         'disulfide bond', 'SG-SG'))
        self.table.setRowCount(len(rows))
        for row, (what, where, label, code) in enumerate(rows):
            item = QtWidgets.QTableWidgetItem(where)
            item.setData(Qt.ItemDataRole.UserRole, what)
            self.table.setItem(row, 0, item)
            self.table.setItem(row, 1, QtWidgets.QTableWidgetItem(label))
            self.table.setItem(row, 2, QtWidgets.QTableWidgetItem(code))
        self.table.resizeColumnsToContents()

    def refresh(self):
        current = self.copy_chain.currentData()
        self.copy_chain.blockSignals(True)
        self.copy_chain.clear()
        for chain in self.job.chains:
            self.copy_chain.addItem('chain %s' % chain.id, chain.id)
        if current:
            index = self.copy_chain.findData(current)
            if index >= 0:
                self.copy_chain.setCurrentIndex(index)
        self.copy_chain.blockSignals(False)
        self._show_copies()
        self.copy_chain.setEnabled(bool(self.job.chains))
        self.copies.setEnabled(bool(self.job.chains))
        self.view.set_job(self.job)
        self._fill_table()
        self.ligand_list.clear()
        for code, copies in self.job.ligands:
            self.ligand_list.addItem(
                code if copies == 1 else '%s  x%d' % (code, copies))
        method = self.method.currentData()
        self.seeds.setEnabled(method == 'af3')
        self.seed_label.setEnabled(method == 'af3')
        self.status.setText(self._status())
        self.save_button.setEnabled(bool(self.job.chains))
        running = self.process is not None
        self.submit_button.setEnabled(method == 'af3' and not running)
        self.submit_button.setToolTip(
            'A prediction started here is already running'
            if running else
            'Run the prediction with a job script, here or through Slurm'
            if method == 'af3' else
            'Submitting is set up for AlphaFold 3 so far')

    def _status(self):
        '''
        What the dialog says below the method: the job, its warnings, and
        the last word from a prediction started here
        '''
        lines = []
        if self.pending:
            lines.append('Disulfide bond from %s %d: click its partner'
                         % self.pending)
        if self.job.chains:
            lines.append(self.job.summary())
            lines += self.job.warnings(self.method.currentData())
        if self.message:
            lines.append(self.message)
        return '\n'.join(lines) or 'Open a FASTA file to begin'

    # the file

    def _seeds(self):
        text = self.seeds.text().replace(',', ' ').split()
        found = [int(one) for one in text if one.lstrip('-').isdigit()]
        return found or [1]

    # running it

    def submit_job(self):
        name = os.path.basename(self.saved) if self.saved else \
            '%s_af3.json' % (self.name.text().strip() or 'prediction')
        dialog = SubmitDialog(name, self, self.prefs)
        dialog.saved = self.saved
        dialog._refresh()
        if dialog.exec() != QtWidgets.QDialog.DialogCode.Accepted:
            return
        try:
            command = prediction.submit_command(
                dialog.script.text().strip(),
                dialog.input_name.text().strip(), dialog.mode())
        except ValueError as e:
            self.status.setText('Not submitted: %s' % e)
            return
        self._run(command, dialog.mode(),
                  os.path.dirname(dialog.script.text().strip()))

    def _run(self, command, mode, folder):
        '''
        Start the job, its output going to the log
        '''
        self.process = QtCore.QProcess(self)
        self.process.setWorkingDirectory(folder or os.getcwd())
        self.process.setProcessChannelMode(
            QtCore.QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._output)
        self.process.finished.connect(
            lambda code, status: self._finished(mode, code))
        print(' Structure Prediction: %s' % ' '.join(command))
        self.message = '%s: starting...' % (
            'Slurm' if mode == 'slurm' else 'Running here')
        self.refresh()
        self.process.start(command[0], command[1:])

    def _output(self):
        if self.process is None:
            return
        text = bytes(self.process.readAllStandardOutput()).decode(
            'utf-8', 'replace')
        for line in text.splitlines():
            if line.strip():
                print('  ' + line)
                self.message = line.strip()[:200]
        self.status.setText(self._status())

    def _finished(self, mode, code):
        # (the process object goes too, or Qt tears it down at exit)
        process, self.process = self.process, None
        if process is not None:
            process.deleteLater()
        if mode == 'slurm':
            self.message = ('Submitted to Slurm: %s' % self.message
                            if code == 0 else
                            'sbatch failed (exit %d), see the output pane'
                            % code)
        else:
            self.message = ('The prediction finished' if code == 0 else
                            'The prediction stopped (exit %d), see the '
                            'output pane' % code)
        self.refresh()

    def closeEvent(self, event):
        # (a run started here belongs to this session)
        if self.process is not None:
            self.message = 'Still running; it stops when PyMOL quits'
            self.status.setText(self._status())
            self.process.setParent(None)      # (it outlives this dialog)
        super().closeEvent(event)

    def save(self):
        self.job.name = self.name.text().strip() or 'prediction'
        self.job.seeds = self._seeds()
        method = self.method.currentData()
        try:
            text = prediction.write(self.job, method)
        except Exception as e:
            self.message = 'Not ready: %s' % e
            self.refresh()
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save the input file',
            '%s_%s.json' % (self.job.name, method), 'JSON (*.json)')
        if not path:
            return
        with open(path, 'w') as handle:
            handle.write(text)
        self.saved = path
        print(' Structure Prediction: wrote %s (%s)' % (
            path, dict((v, k) for k, v in prediction.METHODS)[method]))
        for line in self.job.warnings(method):
            print('  ' + line)
        self.message = 'Wrote %s' % os.path.basename(path)
        self.refresh()


def show_dialog():
    global _dialog
    if _dialog is None:
        from pymol.gui import get_qtwindow
        _dialog = PredictionDialog(get_qtwindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
