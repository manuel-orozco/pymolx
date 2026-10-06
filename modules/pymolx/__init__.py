'''
pymolx: extensions to Open-Source PyMOL working towards parity with
Incentive PyMOL. The gap inventory and work tracker is INCENTIVE_PARITY.md
in the source tree; parity IDs used in this package (e.g. "P-04") refer
to entries in that file.

pymolx is not affiliated with or endorsed by Schrodinger, LLC.
'''

__version__ = '0.1.0'

from . import branding
from . import features
from . import defaults


def _init(_self):
    '''
    Called once from pymol/__init__.py after the "cmd" module is set up.
    Feature modules which replace Incentive-only stubs get imported here.
    '''
    from . import info, undo, interfaces, mm
    info.extend(_self)
    mm.extend(_self)  # before undo.install: one undo step per call
    undo.install(_self)
    _self.extend('interface_analysis', interfaces.interface_analysis)
    _self.interface_analysis = interfaces.interface_analysis


_startup_text = False


def show_startup_text(_self):
    '''
    Like Incentive PyMOL, show the startup text (banner) in the viewer
    until something is loaded or the viewer is clicked.
    '''
    global _startup_text
    _startup_text = True
    _self.set('text', 1, quiet=1)


def hide_startup_text(_self):
    '''
    End the startup text, if it is still shown.
    '''
    global _startup_text
    if _startup_text:
        _startup_text = False
        _self.set('text', 0, quiet=1)


def _ready(_self):
    '''
    Called once from pymol.adapt_to_hardware(), after the adapted state
    is stored as default. Startup changes are not undo steps.
    '''
    from . import undo
    undo.stack.clear_history()


def _started(_self):
    '''
    Called once from pymol.adapt_to_hardware(), when the PyMOL instance
    (including its command parser) is set up and before pymolrc files run.
    '''
    defaults.apply(_self)
