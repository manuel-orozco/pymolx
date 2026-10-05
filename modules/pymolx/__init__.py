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
    from . import info
    info.extend(_self)
