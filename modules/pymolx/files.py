'''
Where loaded objects came from: "load" (wrapped in pymolx.missing)
records each new object's file, so tools can find files next to it (e.g.
the confidence scores of predicted models, pymolx.cluster). Not kept in
sessions; renamed objects are forgotten, and a new object with the name
of a deleted one gets its file.
'''

import os

_sources = {}
# "load" runs scripts too: their objects come from the files they load
SCRIPTS = ('.py', '.pyc', '.pym', '.pml')


def remember(names, filename):
    path = os.path.abspath(os.path.expanduser(filename))
    if path.lower().endswith(SCRIPTS):
        return
    for name in names:
        _sources[name] = path


def source(name):
    '''
    :return: the file object "name" was loaded from, or None
    '''
    path = _sources.get(name)
    return path if path and os.path.exists(path) else None
