'''
Registry of implemented parity items.

A module which implements an item from INCENTIVE_PARITY.md calls
register() with the item's ID when it is imported, e.g.:

    from pymolx import features
    features.register('P-04')

The test harness (testing/testing.py) uses this registry to run the
upstream tests which are marked @testing.requires('incentive') once the
corresponding item is implemented (see testing/pymolx_parity.py).
'''

import re

_ID_RE = re.compile(r'^[A-Z]-\d\d$')

_implemented = set()


def register(*ids):
    '''
    Mark parity items as implemented.
    '''
    for id_ in ids:
        if not _ID_RE.match(id_):
            raise ValueError('invalid parity ID: %r' % (id_,))
        _implemented.add(id_)


def provides(id_):
    '''
    True if the parity item is implemented.
    '''
    return id_ in _implemented


def implemented():
    '''
    Sorted list of implemented parity item IDs.
    '''
    return sorted(_implemented)
