'''
Mouse mode menu with a mouse controls table, like Incentive PyMOL: the
mode list (current mode checked) and, attached to its left, the button
assignments of the current or hovered mode.

The table is built from pymol.controlling.mode_dict, so it shows what
the mode actually binds.
'''

from pymol.Qt import QtCore, QtWidgets

Qt = QtCore.Qt

# mode menu, in groups (like Incentive PyMOL)
MODE_GROUPS = [
    ['three_button_viewing', 'three_button_editing'],
    ['three_button_lights', 'three_button_motions'],
    ['two_button_viewing', 'two_button_selecting', 'two_button_editing'],
    ['two_button_lights'],
    ['one_button_viewing'],
    ['three_button_maestro'],
]

# short action names, indexed by action code (see layer1/ButMode.cpp)
ACTION_LABELS = {
    0: 'Rota', 1: 'Move', 2: 'MovZ', 3: 'Clip', 4: 'RotZ', 5: 'ClpN',
    6: 'ClpF', 7: 'lb', 8: 'mb', 9: 'rb', 10: '+lb', 11: '+mb', 12: '+rb',
    13: 'PkAt', 14: 'PkBd', 15: 'RotF', 16: 'TorF', 17: 'MovF', 18: 'Orig',
    19: '+lBx', 20: '-lBx', 21: 'lbBx', 22: '', 23: 'Cent', 24: 'PkTB',
    25: 'Slab', 26: 'MovS', 27: 'Pk1', 28: 'MovA', 29: 'Menu', 30: 'Sele',
    31: '+/-', 32: '+Box', 33: '-Box', 34: 'MvSZ', 35: 'Clik', 36: 'RotD',
    37: 'MovD', 38: 'MvDZ', 39: 'RotO', 40: 'MovO', 41: 'MvOZ', 42: 'MvFZ',
    43: 'MvAZ', 44: 'DrgM', 45: 'RotV', 46: 'MovV', 47: 'MvVZ', 49: 'DrgO',
    50: 'IMSZ', 51: 'IMvZ', 52: 'Box', 53: 'IRtZ', 54: 'RotL', 55: 'MovL',
    56: 'MvzL',
}

COLUMNS = [('L', 'l'), ('M', 'm'), ('R', 'r'), ('Wheel', 'w')]
ROWS = [  # label, modifier, button prefix ('' = drag)
    ('', 'none', ''),
    ('Shift', 'shft', ''),
    ('Ctrl', 'ctrl', ''),
    ('Ctrl+Shift', 'ctsh', ''),
    ('Single Click', 'none', 'single_'),
    ('Double Click', 'none', 'double_'),
]
CLICK_BUTTONS = {'l': 'left', 'm': 'middle', 'r': 'right'}


def mode_names():
    from pymol import controlling
    return controlling.mode_name_dict


def action_label(action):
    from pymol import controlling
    code = controlling.but_act_code.get(action.lower())
    return ACTION_LABELS.get(code, action)


def controls_table(mode):
    '''
    :return: {(row label, column label): action label}
    '''
    from pymol import controlling
    bindings = {(button, modifier): action
                for button, modifier, action in controlling.mode_dict[mode]}
    table = {}
    for row, modifier, prefix in ROWS:
        for column, button in COLUMNS:
            if prefix:
                if button == 'w':
                    continue
                button = prefix + CLICK_BUTTONS[button]
            action = bindings.get((button, modifier))
            table[(row, column)] = action_label(action) if action else ''
    return table


class ControlsPanel(QtWidgets.QFrame):
    '''
    The mouse controls table, shown next to the mode menu.
    '''

    def __init__(self, parent=None):
        super().__init__(parent, Qt.WindowType.ToolTip |
                         Qt.WindowType.FramelessWindowHint)
        self.setObjectName('mouse_controls_panel')
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 8)
        layout.setSpacing(0)

        self.title = QtWidgets.QLabel(self)
        self.title.setObjectName('mouse_controls_title')
        self.title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.title)

        grid = QtWidgets.QGridLayout()
        grid.setContentsMargins(12, 8, 16, 0)
        grid.setHorizontalSpacing(18)
        grid.setVerticalSpacing(6)
        layout.addLayout(grid)
        layout.addStretch(1)

        header = QtWidgets.QLabel('Buttons &\nKeys', self)
        header.setObjectName('mouse_controls_keys')
        header.setAlignment(Qt.AlignmentFlag.AlignRight)
        grid.addWidget(header, 0, 0)
        for c, (column, _) in enumerate(COLUMNS, 1):
            label = QtWidgets.QLabel(column, self)
            label.setObjectName('mouse_controls_axis')
            label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            grid.addWidget(label, 0, c)

        self.cells = {}
        for r, (row, _, _) in enumerate(ROWS, 1):
            label = QtWidgets.QLabel(row, self)
            label.setObjectName('mouse_controls_axis')
            label.setAlignment(Qt.AlignmentFlag.AlignRight)
            grid.addWidget(label, r, 0)
            for c, (column, _) in enumerate(COLUMNS, 1):
                cell = QtWidgets.QLabel(self)
                cell.setObjectName('mouse_controls_cell')
                cell.setAlignment(Qt.AlignmentFlag.AlignCenter)
                grid.addWidget(cell, r, c)
                self.cells[(row, column)] = cell

    def show_mode(self, mode):
        self.title.setText('%s Mouse Controls' % mode_names().get(mode, mode))
        for key, text in controls_table(mode).items():
            self.cells[key].setText(text)


class MouseModeMenu(QtWidgets.QMenu):
    '''
    Mouse mode list; the controls table follows the current or hovered
    mode.
    '''

    def __init__(self, cmd, parent=None):
        super().__init__(parent)
        self.cmd = cmd
        self.setObjectName('mouse_mode_menu')
        self.panel = ControlsPanel(self)
        self.group = QtWidgets.QActionGroup(self)
        self.mode_actions = {}

        names = mode_names()
        for i, group in enumerate(MODE_GROUPS):
            if i:
                self.addSeparator()
            for mode in group:
                action = self.addAction(names.get(mode, mode))
                action.setCheckable(True)
                action.setData(mode)
                action.triggered.connect(
                    lambda _=False, m=mode: self.cmd.mouse(m))
                self.group.addAction(action)
                self.mode_actions[mode] = action

        self.aboutToShow.connect(self._update_checked)
        self.aboutToHide.connect(self.panel.hide)
        self.hovered.connect(self._hovered)

    def current_mode(self):
        name = self.cmd.get('button_mode_name')
        for mode, label in mode_names().items():
            if label == name:
                return mode
        return None

    def _update_checked(self):
        mode = self.current_mode()
        for m, action in self.mode_actions.items():
            action.setChecked(m == mode)
        self.panel.show_mode(mode or 'three_button_viewing')

    def _hovered(self, action):
        mode = action.data()
        if mode:
            self.panel.show_mode(mode)

    def showEvent(self, event):
        super().showEvent(event)
        # attach the table to the left of the menu, same height
        geometry = self.geometry()
        self.panel.adjustSize()
        width = self.panel.sizeHint().width()
        self.panel.setGeometry(geometry.left() - width + 1, geometry.top(),
                               width, geometry.height())
        self.panel.show()
        self.panel.raise_()
