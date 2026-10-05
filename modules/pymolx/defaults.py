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

# values replaced by apply(), so revert() can restore them
_replaced = {}


def apply(_self):
    '''
    Apply DEFAULTS. Called from pymol.adapt_to_hardware() right before
    the adapted state is stored with "reinitialize store".
    '''
    for name, value in DEFAULTS.items():
        _replaced.setdefault(name, _self.get(name))
        _self.set(name, value, quiet=1)


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
