'''
Folders pymolx uses: the source checkout it was built from, and the
scripts folder (opened by the Scripts menu).
'''

import os

SCRIPTS = 'scripts'
# without a source checkout (e.g. installed from a wheel)
USER_SCRIPTS_DIR = os.path.join(
    os.environ.get('XDG_DATA_HOME') or os.path.expanduser('~/.local/share'),
    'pymolx', SCRIPTS)


def source_dir():
    '''
    :return: the source checkout pymolx was built from (recorded by
        setup.py), or None
    '''
    try:
        from ._source import SOURCE_DIR
    except ImportError:
        return None
    return SOURCE_DIR if os.path.isdir(SOURCE_DIR) else None


def scripts_dir(create=False):
    '''
    :param create: make the folder if it doesn't exist
    :return: "scripts" in the source checkout, else USER_SCRIPTS_DIR
    '''
    source = source_dir()
    path = os.path.join(source, SCRIPTS) if source else USER_SCRIPTS_DIR
    if create:
        os.makedirs(path, exist_ok=True)
    return path
