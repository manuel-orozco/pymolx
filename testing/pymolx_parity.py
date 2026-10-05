'''
pymolx: maps parity items (INCENTIVE_PARITY.md) to the upstream tests
which are marked @testing.requires('incentive') and verify them.

testing.py runs such a test (instead of skipping it) once every parity
item it is listed under is registered in pymolx.features. This keeps the
upstream test files unchanged.

Test IDs are "<file relative to tests/>" (whole file) or
"<file>::<qualified name>" of the class or method which carries the
requires('incentive') decorator.
'''

PARITY_TESTS = {
    'P-01': [  # clean
        'api/computing.py::TestComputing.testClean',
        'jira/PYMOL-2985.py',
    ],
    'P-02': [  # assign_stereo
        'api/stereochemistry.py',
    ],
    'P-03': [  # pi_interactions
        'api/querying.py::TestQuerying.testPiInteractions',
    ],
    'P-06': [  # focal_blur
        'api/experimenting.py::TestExperimenting.testFocalblur',
    ],
    'P-07': [  # morph (TestMorphing.testMorphRigimol imports epymol.rigimol)
        'api/morphing.py',
        'jira/PYMOL-174.py',
    ],
    'P-08': [  # help_setting
        'api/helping.py::TestHelping.testHelpSetting',
    ],
    'F-01': [  # MAE load
        'api/load_mae.py',
        'api/symop.py::TestBondSymOp.test_load_mae',
        'api/importing.py::TestImporting.testLoad_mimic',
        'jira/PYMOL-771.py',
        'jira/PYMOL-1498.py',
        'jira/PYMOL-1514.py',
        'jira/PYMOL-1567.py',
        'jira/PYMOL-2838.py',
        'jira/PYMOL-3221.py',
        'properties/test_simple.py',
        'properties/test_atom_props.py',
        'properties/test_copy_props.py',
        'properties/mae_io.py',
        'jira/PYMOL-3231.py',
    ],
    'F-02': [  # MAE save
        'properties/test_simple.py::TestProperties.testMAEsaveLoadSessions',
        'properties/test_atom_props.py::TestAtomProperties.testMAEsaveLoadSessionsWithAtomProperties',
        'properties/test_atom_props.py::TestAtomProperties.testMAEsaveLoadSessionsWithAtomPropertiesBinaryDump',
        'properties/mae_io.py',
        'jira/PYMOL-3231.py',
    ],
    'F-03': [  # STL export
        'api/exporting_geom.py::TestExportingGeom.testSTL',
    ],
    'F-06': [  # MOE
        'jira/PYMOL-1191.py',
    ],
    'F-07': [  # VIS
        'jira/PYMOL-1441.py',
    ],
    'F-08': [  # Phase hypotheses
        'api/epymol_ph4.py',
    ],
    'F-09': [  # pdb_header object property
        'api/importing.py::TestImporting.testPdbHeader',
    ],
    'C-01': [  # sphere_mode 10/11
        'cgos/reps.py::TestReps.test_sphere_mode_10_11',
    ],
    'C-03': [  # stereo / text_type atom properties
        'api/mmpymolx.py',
        'jira/PYMOL-317.py',
        'jira/PYMOL-1276.py',
    ],
    'C-07': [  # antialias_shader
        'settings/settings.py::TestSettings.testAA',
    ],
    'B-03': [  # bg_image_filename embedded in sessions
        'jira/PYMOL-1571.py',
    ],
    'B-04': [  # undo after remove on discrete objects
        'jira/PYMOL-1697.py::TestPYMOL1697.testUndoAfterRemoveAtomOnDiscrete',
    ],
}


def required_items(test_id):
    '''
    Parity items which must be implemented for a test to run.

    :param test_id: "<file>::<qualname>" of the decorated class or method
    :return: set of parity IDs (empty if the test isn't mapped)
    '''
    filename, _, qualname = test_id.partition('::')
    required = set()
    for item, entries in PARITY_TESTS.items():
        for entry in entries:
            entry_file, _, entry_qualname = entry.partition('::')
            if entry_file != filename:
                continue
            if (not entry_qualname or entry_qualname == qualname or
                    qualname.startswith(entry_qualname + '.')):
                required.add(item)
    return required


def unlocked(test_id, provides):
    '''
    True if the test is mapped and all its parity items are implemented.

    :param provides: callable(parity_id) -> bool, e.g. pymolx.features.provides
    '''
    required = required_items(test_id)
    return bool(required) and all(provides(item) for item in required)
