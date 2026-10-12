'''
File > Open Recent: up to 20 files, newest first.
'''

import sqlite3

import pytest

from pymol._gui import PyMOLDesktopGUI


@pytest.fixture
def gui():
    g = PyMOLDesktopGUI()
    # in memory, never the user's ~/.pymol/recent.db
    db = sqlite3.connect(':memory:')
    db.execute('CREATE TABLE recent (filename text unique, timestamp integer)')
    g._recent_filenames_db = db
    return g


def test_keeps_the_20_newest(gui):
    assert gui.RECENT_FILES_MAX == 20
    # 25 files within the same second: order still right
    for i in range(25):
        gui.recent_filenames_add('/data/f%02d.pdb' % i)
    recent = gui.recent_filenames
    assert recent == ['/data/f%02d.pdb' % i for i in range(24, 4, -1)]
    count = gui._recent_filenames_db.execute(
        'SELECT COUNT(*) FROM recent').fetchone()[0]
    assert count == 20  # not cut back to 16


def test_reopening_moves_to_top(gui):
    for i in range(5):
        gui.recent_filenames_add('/data/f%d.pdb' % i)
    gui.recent_filenames_add('/data/f1.pdb')
    assert gui.recent_filenames[:3] == [
        '/data/f1.pdb', '/data/f4.pdb', '/data/f3.pdb']
    assert len(gui.recent_filenames) == 5  # no duplicates


def test_older_databases_with_more_entries(gui):
    # a database written by another PyMOL can hold more: show 20
    for i in range(30):
        gui._recent_filenames_db.execute(
            "INSERT INTO recent VALUES (?, datetime('now', ?))",
            ('/old/f%02d.pdb' % i, '-%d minutes' % (30 - i)))
    assert len(gui.recent_filenames) == 20
    assert gui.recent_filenames[0] == '/old/f29.pdb'
