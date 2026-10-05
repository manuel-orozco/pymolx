'''
Dark theme with a purple accent (parity item L-01).

Colors are defined once in COLORS. They feed the Qt palette (for widgets
the stylesheet doesn't cover, like dialogs) and the "$name" tokens in
data/pymolx/styles/dark.qss.
'''

import string

COLORS = {
    'window': '#2b2b2e',      # main window, toolbars
    'panel': '#333337',       # menus, docks, buttons
    'base': '#1e1e21',        # text areas: output, command line
    'text': '#dcdce0',
    'text_dim': '#9a9aa2',
    'border': '#44444a',
    'hover': '#3d3d44',
    'accent': '#9b7be0',      # active state, highlights
    'accent_text': '#b9a2f0', # prompt, selection mode
    'selection': '#5a4a8a',   # selected text
    'enabled': '#5cb85c',     # enabled object dot, "All" C button
    'button': '#4a4a52',      # A/S/H/L/C buttons
}

STYLESHEET_PATH = '$PYMOL_DATA/pymolx/styles/dark.qss'
ICON_DIR = '$PYMOL_DATA/pymolx/icons'


def stylesheet(_self):
    '''
    The theme stylesheet with color and "$icons" tokens substituted.

    :raises KeyError: if the stylesheet uses an unknown token
    '''
    with open(_self.exp_path(STYLESHEET_PATH)) as handle:
        template = string.Template(handle.read())
    # Qt stylesheet urls need forward slashes, also on Windows
    icons = _self.exp_path(ICON_DIR).replace('\\', '/')
    return template.substitute(COLORS, icons=icons)


def palette():
    '''
    QPalette matching COLORS.
    '''
    from pymol.Qt import QtGui
    c = {k: QtGui.QColor(v) for k, v in COLORS.items()}
    role = QtGui.QPalette.ColorRole
    group = QtGui.QPalette.ColorGroup

    pal = QtGui.QPalette()
    for r, color in [
        (role.Window, c['window']),
        (role.WindowText, c['text']),
        (role.Base, c['base']),
        (role.AlternateBase, c['panel']),
        (role.ToolTipBase, c['panel']),
        (role.ToolTipText, c['text']),
        (role.PlaceholderText, c['text_dim']),
        (role.Text, c['text']),
        (role.Button, c['panel']),
        (role.ButtonText, c['text']),
        (role.BrightText, c['accent_text']),
        (role.Highlight, c['selection']),
        (role.HighlightedText, c['text']),
        (role.Link, c['accent_text']),
        (role.Light, c['hover']),
        (role.Midlight, c['panel']),
        (role.Mid, c['border']),
        (role.Dark, c['base']),
        (role.Shadow, c['base']),
    ]:
        pal.setColor(r, color)

    for r in (role.WindowText, role.Text, role.ButtonText):
        pal.setColor(group.Disabled, r, c['text_dim'])

    return pal


def apply(app, window, _self):
    '''
    Apply the theme: Fusion style and palette on the application, the
    stylesheet on the main window (appended to any existing one).
    '''
    app.setStyle('Fusion')
    app.setPalette(palette())
    window.setStyleSheet(window.styleSheet() + '\n' + stylesheet(_self))


def feedback_html(text):
    '''
    Output pane HTML for feedback text, with echoed commands ("PyMOL>")
    in the accent color.
    '''
    from pymol import colorprinting
    lines = []
    for line in text.split('\n'):
        html = colorprinting.text2html(line)
        if line.startswith('PyMOL>'):
            html = '<span style="color:%s">%s</span>' % (
                COLORS['accent_text'], html)
        lines.append(html)
    return '<br>'.join(lines)
