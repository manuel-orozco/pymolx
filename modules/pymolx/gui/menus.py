'''
Qt menus built from pymol.menu data, the menu definitions the OpenGL
popups use (Presets, the content panel's A/S/H/L/C menus, mouse modes).

pymol.menu data is a list of [kind, label, command] entries:
    kind 0: separator
    kind 1: item; command is a command string, or a list (submenu)
    kind 2: title (first entry) or note line
Labels may contain PyMOL color codes: "\\RGB" with digits 0-9, e.g.
"\\900red" is red text. PyMenu draws them like the OpenGL menus do.
'''

import re

from pymol.Qt import QtCore, QtGui, QtWidgets

from . import theme

Qt = QtCore.Qt

_COLOR_CODE = re.compile(r'\\(\d\d\d)')

# must match the left padding of "QMenu::item" in dark.qss
TEXT_LEFT = 20


def strip_label(label):
    '''
    Remove PyMOL color codes and surrounding whitespace from a menu label.
    '''
    return _COLOR_CODE.sub('', label).strip()


def label_segments(label):
    '''
    Split a label with PyMOL color codes into colored runs.

    :return: list of (color, text); color is "#rrggbb" or None (default)
    '''
    label = label.strip()
    segments = []
    color = None
    pos = 0
    for match in _COLOR_CODE.finditer(label):
        if match.start() > pos:
            segments.append((color, label[pos:match.start()]))
        color = '#' + ''.join('%02x' % round(int(d) * 255 / 9)
                              for d in match.group(1))
        pos = match.end()
    if pos < len(label):
        segments.append((color, label[pos:]))
    return segments


class PyMenu(QtWidgets.QMenu):
    '''
    QMenu which draws item labels with PyMOL color codes in color, and
    a title header. The stylesheet makes the regular item text
    transparent ("QMenu#pymol_menu"); this class draws it instead.
    '''

    def __init__(self, title='', parent=None):
        super().__init__(title, parent)
        self.setObjectName('pymol_menu')

    def add_label_action(self, label, callback=None):
        action = self.addAction(strip_label(label).replace('&', '&&'))
        action.setData(label)
        if callback is not None:
            action.triggered.connect(callback)
        return action

    def add_title(self, label):
        action = self.add_label_action(label)
        action.setEnabled(False)
        action.setProperty('pymol_title', True)
        return action

    def paintEvent(self, event):
        super().paintEvent(event)
        painter = QtGui.QPainter(self)
        palette = self.palette()
        default = palette.color(QtGui.QPalette.ColorRole.Text)
        dim = palette.color(QtGui.QPalette.ColorGroup.Disabled,
                            QtGui.QPalette.ColorRole.Text)
        fonts = {False: self.font(), True: QtGui.QFont(self.font())}
        fonts[True].setBold(True)

        for action in self.actions():
            if action.isSeparator() or not action.isVisible():
                continue
            rect = self.actionGeometry(action)
            if not rect.intersects(event.rect()):
                continue

            label = action.data()
            if not isinstance(label, str):
                label = action.text().replace('&&', '&')
            is_title = bool(action.property('pymol_title'))
            if is_title:
                painter.fillRect(rect, QtGui.QColor(theme.COLORS['selection']))

            painter.setFont(fonts[is_title])
            metrics = QtGui.QFontMetrics(fonts[is_title])
            x = rect.left() + TEXT_LEFT
            plain = default if (action.isEnabled() or is_title) else dim
            for color, text in label_segments(label):
                painter.setPen(QtGui.QColor(color) if color else plain)
                painter.drawText(
                    QtCore.QRect(x, rect.top(), rect.right() - x,
                                 rect.height()),
                    int(Qt.AlignmentFlag.AlignLeft |
                        Qt.AlignmentFlag.AlignVCenter |
                        Qt.TextFlag.TextSingleLine), text)
                x += metrics.horizontalAdvance(text)
        painter.end()


def fill_menu(menu, data, _self):
    '''
    Fill a PyMenu from pymol.menu data.

    :param menu: PyMenu
    :param data: pymol.menu data
    :param _self: PyMOL cmd module the commands run in
    '''
    for i, (kind, label, command) in enumerate(entry[:3] for entry in data):
        if kind == 0:
            menu.addSeparator()
        elif kind == 2:
            if i == 0:
                menu.add_title(label)
            else:
                menu.add_label_action(label).setEnabled(False)
        elif isinstance(command, list):
            # explicit parent: PySide6 deletes menus created by
            # addMenu(title) once the Python reference is gone
            submenu = PyMenu(strip_label(label).replace('&', '&&'), menu)
            menu.addMenu(submenu).setData(label)
            fill_menu(submenu, command, _self)
        else:
            menu.add_label_action(label,
                                  lambda _=False, c=command: _self.do(c, echo=0))
    return menu
