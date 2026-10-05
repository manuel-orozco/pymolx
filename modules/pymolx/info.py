'''
"pymolx" command: version, implemented parity items and optional
dependencies (parity item D-03).
'''

import importlib.util
import shutil

from pymol import cmd

# import name -> what needs it (see INCENTIVE_PARITY.md)
OPTIONAL_MODULES = {
    'numpy': 'core',
    'scipy': 'P-07 rigid-body morph',
    'rdkit': 'P-01 clean, P-02 assign_stereo, P-03 pi_interactions, C-03',
    'gemmi': 'P-09 load_mtz',
    'Bio': 'sequence tools (biopython)',
}

# executable -> what needs it
OPTIONAL_PROGRAMS = {
    'apbs': 'D-01 electrostatics (apbs_gui plugin)',
    'pdb2pqr30': 'D-01 electrostatics (apbs_gui plugin)',
    'ffmpeg': 'D-02 MP4 movie export',
}


def dependency_status():
    '''
    :return: two dicts {name: bool} for optional modules and programs
    '''
    modules = {
        name: importlib.util.find_spec(name) is not None
        for name in OPTIONAL_MODULES
    }
    programs = {
        name: shutil.which(name) is not None
        for name in OPTIONAL_PROGRAMS
    }
    return modules, programs


def pymolx(quiet=0, *, _self=cmd):
    '''
DESCRIPTION

    Show pymolx version, implemented Incentive parity items, and which
    optional dependencies are available.

USAGE

    pymolx
    '''
    from . import branding, features

    modules, programs = dependency_status()

    if not int(quiet):
        print(' ' + branding.version_message(_self.get_version()[0]))
        print(' ' + branding.DISCLAIMER)
        print(' Implemented parity items: %s' % (
            ', '.join(features.implemented()) or 'none'))
        print(' Optional dependencies:')
        for table, status in ((OPTIONAL_MODULES, modules),
                              (OPTIONAL_PROGRAMS, programs)):
            for name, purpose in table.items():
                print('  %-10s %-8s %s' % (
                    name, 'found' if status[name] else 'MISSING', purpose))

    return {
        'version': branding.version(),
        'features': features.implemented(),
        'modules': modules,
        'programs': programs,
    }


def extend(_self):
    # keyword for the command language, attribute for the Python API
    _self.extend('pymolx', pymolx)
    _self.pymolx = pymolx
