'''
pymolx startup defaults (parity item L-04).

Settings in DEFAULTS are applied once at startup, after hardware
adaptation and before pymolrc files run (so users can still override
them). They also become the state that "reinitialize" restores.

Only add a value after confirming it side by side with Incentive PyMOL
(record the comparison in docs/parity/), not from memory.
'''

DEFAULTS = {
    # 'setting_name': value,
}

# Display > Quality level (util.performance mode) applied before DEFAULTS:
# 0 = Maximum Quality (owner's choice), None = leave PyMOL's defaults.
# Depends on use_shaders, so it runs after the graphics are detected.
QUALITY = 0

# values replaced by apply(), so revert() can restore them
_replaced = {}


def apply(_self):
    '''
    Apply QUALITY, then DEFAULTS. Called from pymol.adapt_to_hardware()
    right before the adapted state is stored with "reinitialize store".
    '''
    if QUALITY is not None:
        from pymol import util
        util.performance(QUALITY, _self=_RecordingCmd(_self))
    for name, value in DEFAULTS.items():
        _replaced.setdefault(name, _self.get(name))
        _self.set(name, value, quiet=1)


class _RecordingCmd:
    '''
    cmd stand-in for util.performance: sets quietly and records the
    values it replaces for revert()
    '''

    def __init__(self, _self):
        self._self = _self

    def __getattr__(self, name):
        return getattr(self._self, name)

    def set(self, name, value=1, *args, **kwargs):
        _replaced.setdefault(name, self._self.get(name))
        self._self.set(name, value, quiet=1)

    def do(self, *args, **kwargs):
        # util.performance ends with "rebuild": nothing is loaded yet at
        # startup, and it would echo "PyMOL>rebuild" in the output
        pass


def revert(_self):
    '''
    Restore the values which apply() replaced. The upstream test suite
    is written against upstream defaults and calls this at startup.

    :return: True if any setting was restored
    '''
    for name, value in _replaced.items():
        _self.set(name, value, quiet=1)
    reverted = bool(_replaced)
    _replaced.clear()
    return reverted
