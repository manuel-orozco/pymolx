'''
Dark theme with a purple accent (parity item L-01).

Colors are defined once in COLORS. They feed the Qt palette (for widgets
the stylesheet doesn't cover, like dialogs) and the "$name" tokens in
data/pymolx/styles/dark.qss.
'''

import string

COLORS = {
    # measured from Incentive PyMOL 3.1 screenshots: neutral greys
    'window': '#343434',      # menu bar, toolbars, content panel
    'panel': '#3e3e3e',       # inputs, buttons
    'base': '#222222',        # output pane, active toggles
    'menu': '#2e2e2e',        # popup menus
    'text': '#cecece',
    'text_bright': '#e0e0e0', # output text
    'text_dim': '#7a7a7a',
    'border': '#494949',      # 1px dividers
    'hover': '#404040',
    'accent': '#9b81fd',      # selection mode, prompt, active toggles
    'accent_text': '#9b81fd',
    'selection': '#4b3d7a',   # selected text and menu items
    'enabled': '#33ff33',     # enabled entry dot
    'button': '#575757',      # A/S/H/L/C buttons of enabled entries
    # A/S/H/L/C buttons of the "All" row
    'all_a': '#7b519e',
    'all_s': '#7e5da8',
    'all_h': '#816eb7',
    'all_l': '#837dc5',
    'all_c': '#54bd41',
}

# all GUI text (menus, toolbars, panels, dialogs, output, command line,
# sequence viewer), in points
FONT_SIZE = 11

# output pane, command line and sequence viewer font: first installed
# family wins (DejaVu Sans Mono comes with most Linux systems; the others
# are fallbacks for Windows and macOS)
CONSOLE_FONTS = ['DejaVu Sans Mono', 'Consolas', 'Cascadia Mono',
                 'Inconsolata']

STYLESHEET_PATH = '$PYMOL_DATA/pymolx/styles/dark.qss'
ICON_DIR = '$PYMOL_DATA/pymolx/icons'


def stylesheet(_self):
    '''
    The theme stylesheet with color, "$icons" and "$font_size" tokens
    substituted.

    :raises KeyError: if the stylesheet uses an unknown token
    '''
    with open(_self.exp_path(STYLESHEET_PATH)) as handle:
        template = string.Template(handle.read())
    # Qt stylesheet urls need forward slashes, also on Windows
    icons = _self.exp_path(ICON_DIR).replace('\\', '/')
    return template.substitute(COLORS, icons=icons,
                               font_size='%dpt' % FONT_SIZE)


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


# fonts offered for the sequence viewer (Display > Sequence > Font):
# menu name -> (family, bold, how to get it if not installed)
_MS_CORE_FONTS = 'sudo apt install ttf-mscorefonts-installer'
SEQUENCE_FONTS = {
    'DejaVu Sans Mono': ('DejaVu Sans Mono', False,
                         'sudo apt install fonts-dejavu-core'),
    'Courier New': ('Courier New', False, _MS_CORE_FONTS),
    'Courier New Bold': ('Courier New', True, _MS_CORE_FONTS),
    'Aptos Mono': ('Aptos Mono', False,
                   'download "Microsoft Aptos Fonts" from microsoft.com '
                   'and unzip it into ~/.local/share/fonts'),
    'Consolas': ('Consolas', False,
                 'copy consola*.ttf from Windows (C:\\Windows\\Fonts) '
                 'to ~/.local/share/fonts'),
}


def installed_families():
    from pymol.Qt import QtGui
    try:
        return set(QtGui.QFontDatabase.families())    # Qt 6
    except TypeError:
        return set(QtGui.QFontDatabase().families())  # Qt 5


def console_font(families=None, size=None, family=None, bold=False):
    '''
    QFont for the output pane, command line and sequence viewer: the
    first installed family of CONSOLE_FONTS, else the system monospace
    font. Its letters are as tall as the interface font's at the same
    size (monospace fonts run larger), so all text looks the same size.

    :param families: installed families (default: ask Qt)
    :param size: points, as for the interface font (default: FONT_SIZE)
    :param family: wanted family, if installed (default: CONSOLE_FONTS)
    :param bold: bold weight (e.g. "Courier New Bold", thin otherwise)
    '''
    from pymol.Qt import QtGui, QtWidgets
    installed = set(families) if families is not None else \
        installed_families()
    if family not in installed:
        family = next((f for f in CONSOLE_FONTS if f in installed),
                      'Monospace')
    size = size or FONT_SIZE
    font = QtGui.QFont(family)
    font.setStyleHint(QtGui.QFont.StyleHint.Monospace)
    font.setBold(bool(bold))
    # scale exactly: hinting snaps letter heights to a few pixel sizes
    # (DejaVu's capitals jump from 10 to 12 px around 11 pt)
    font.setHintingPreference(QtGui.QFont.HintingPreference.PreferNoHinting)
    ui = QtGui.QFont(QtWidgets.QApplication.font())
    # the fonts' design ratio, at a size where pixel rounding is negligible
    for f in (ui, font):
        f.setPointSizeF(100)
    font.setPointSizeF(round(size * _letter_height(ui) /
                             _letter_height(font), 1))
    return font


def _letter_height(font):
    from pymol.Qt import QtGui
    metrics = QtGui.QFontMetricsF(font)
    if hasattr(metrics, 'capHeight'):  # Qt 6
        return metrics.capHeight()
    return metrics.xHeight()


def apply(app, window, _self):
    '''
    Apply the theme: Fusion style, palette and FONT_SIZE on the
    application, the stylesheet on the main window (appended to any
    existing one).
    '''
    from pymol.Qt import QtGui
    app.setStyle('Fusion')
    app.setPalette(palette())
    font = QtGui.QFont(app.font())
    font.setPointSize(FONT_SIZE)
    app.setFont(font)
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
