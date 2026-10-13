'''
Hover highlighting in the 3D viewer (pymolx.hover): the atom under the
mouse pointer is picked (_cmd.pick_at) at most every INTERVAL_MS while no
button is pressed, and its unit (mouse_selection_mode) is highlighted.
Display > Hover Highlight switches it on and off.

Picking draws atoms in code colors (mostly red) into the frame buffer, so
it runs in the GL widget's paintGL, right before the scene is drawn
(before_draw): run anywhere else, Qt could show the pick image whenever
another widget repaints (a red flash).
'''

from pymol import _cmd
from pymol.Qt import QtCore, QtWidgets

from .. import hover
from ..sequence import uncounted

Qt = QtCore.Qt
INTERVAL_MS = 40
PREF_HOVER = 'pymolx_hover_highlight'


class HoverHighlighter(QtCore.QObject):
    '''
    Event filter on the OpenGL widget
    '''

    def __init__(self, window, prefs=None):
        super().__init__(window)
        self.window = window
        self.cmd = window.cmd
        self.widget = window.pymolwidget
        self.prefs = prefs
        self.enabled = True
        self._pos = None
        self._pressed = False
        self._current = None  # key of the highlighted unit
        self._pick_pending = False
        self.timer = QtCore.QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.setInterval(INTERVAL_MS)
        self.timer.timeout.connect(self.update)
        # moves without a pressed button only reach widgets that track
        self.widget.setMouseTracking(True)
        self.widget.installEventFilter(self)
        self.widget.before_draw.append(self.pick)

    def set_enabled(self, enabled, save=True):
        self.enabled = bool(enabled)
        if save and self.prefs is not None:
            self.prefs.set(PREF_HOVER, self.enabled)
        if not self.enabled:
            self.clear()

    def eventFilter(self, obj, event):
        kind = event.type()
        Type = QtCore.QEvent.Type
        if kind == Type.MouseMove:
            if event.buttons() == Qt.MouseButton.NoButton:
                self._pos = (event.position().toPoint()
                             if hasattr(event, 'position') else event.pos())
                if not self.timer.isActive():
                    self.timer.start()
        elif kind == Type.MouseButtonPress:
            # rotating, dragging, clicking: out of the way
            self._pressed = True
            self.clear()
        elif kind == Type.MouseButtonRelease:
            self._pressed = False
            self.timer.start()
        elif kind == Type.Leave:
            self._pos = None
            self.clear()
        return False

    def _active(self):
        return self.enabled and not self._pressed and \
            self._pos is not None and not self.cmd.get_setting_int('text')

    def update(self):
        '''
        Highlight what is under the pointer: pick in the next paintGL
        '''
        if not self._active():
            self._pick_pending = False
            self.clear()
            return
        self._pick_pending = True
        self.widget.update()

    def pick(self):
        '''
        before_draw hook (paintGL, context current)
        '''
        if not self._pick_pending:
            return
        self._pick_pending = False
        if not self._active():
            self.clear()
            return
        # (the queries' temporary selections are no news for the sequence
        # viewer)
        with uncounted(_self=self.cmd):
            self._update()

    def _update(self):
        cmd = self.cmd
        g = self.widget
        x = int(g.fb_scale * self._pos.x())
        y = int(g.fb_scale * (g.height() - self._pos.y()))
        with g:  # (PyMOL only draws with a valid context, even in paintGL)
            picked = _cmd.pick_at(cmd._COb, x, y)
        if picked is None:
            self.clear()
            return
        mode = cmd.get_setting_int('mouse_selection_mode')
        key = hover.unit_key(*picked, mode=mode, _self=cmd)
        if key == self._current:
            return
        self._current = key
        if not hover.show(hover.unit_selection(*picked, mode=mode, _self=cmd),
                          _self=cmd):
            self._current = None

    def clear(self):
        if self._current is not None:
            self._current = None
            hover.clear(_self=self.cmd)


def setup(window, prefs=None):
    '''
    Hover highlighting for the main window, and Display > Hover Highlight
    '''
    if prefs is None:
        from .sequence_viewer import PluginPrefs
        prefs = PluginPrefs()
    highlighter = HoverHighlighter(window, prefs)

    action = QtWidgets.QAction('Hover Highlight', window)
    action.setCheckable(True)
    action.setChecked(True)
    action.toggled.connect(highlighter.set_enabled)
    menu = getattr(window, 'menudict', {}).get('Display')
    if menu is not None:
        before = next((a for a in menu.actions()
                       if a.text() == 'Orthoscopic View'), None)
        if before is not None:
            menu.insertAction(before, action)
        else:
            menu.addAction(action)
    highlighter.action = action

    def load_prefs():
        # plugin preferences are read after the window is built
        enabled = bool(prefs.get(PREF_HOVER, True))
        highlighter.set_enabled(enabled, save=False)
        action.blockSignals(True)
        action.setChecked(enabled)
        action.blockSignals(False)

    QtCore.QTimer.singleShot(0, load_prefs)
    return highlighter
