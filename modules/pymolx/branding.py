'''
Product identity for pymolx.

The PyMOL license requires that Schrodinger's copyright notice is kept,
and the trademark notice requires derivative software to be "plainly
distinguished" from Schrodinger's PyMOL products. All user-visible
product naming should come from this module.
'''

NAME = 'pymolx'

# Reverse-DNS id for desktop integration (Wayland/GNOME .desktop name,
# Windows AppUserModelID)
APP_ID = 'io.github.manuel_orozco.pymolx'

HOMEPAGE = 'https://github.com/manuel-orozco/pymolx'
ISSUES_URL = HOMEPAGE + '/issues'

UPSTREAM_NAME = 'Open-Source PyMOL'

DISCLAIMER = ('Based on Open-Source PyMOL(TM). Not affiliated with or '
              'endorsed by Schr\xF6dinger, LLC.')

# Must be kept (see LICENSE)
UPSTREAM_COPYRIGHT = ('Open-Source PyMOL is Copyright (C) Schr\xF6dinger, LLC.'
                      ' All rights reserved.')

SPLASH_PNG = '$PYMOL_DATA/pymolx/splash.png'

# Show SPLASH_PNG in the viewer at startup. Off: like Incentive PyMOL, the
# viewer shows the startup text (banner) until something is loaded.
SHOW_SPLASH_IMAGE = False


def version():
    from . import __version__
    return __version__


def version_message(upstream_version):
    '''
    Product + version string for the startup banner and About dialog.

    :param upstream_version: Open-Source PyMOL version string, e.g. "3.2.0a0"
    '''
    return '%s %s (%s %s)' % (NAME, version(), UPSTREAM_NAME, upstream_version)
