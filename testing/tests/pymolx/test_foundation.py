'''
pymolx foundation: branding, parity item registry, startup defaults,
optional dependency report, and the map from parity items to upstream
Incentive-only tests (INCENTIVE_PARITY.md: T-01, L-04, H-01..H-05, D-03)
'''

import ast
import importlib.util
import re
from pathlib import Path

import pytest

import pymol
from pymol import cmd
import pymolx
from pymolx import branding, defaults, features

TESTING_DIR = Path(__file__).resolve().parents[2]
TESTS_DIR = TESTING_DIR / 'tests'
INVENTORY = TESTING_DIR.parent / 'INCENTIVE_PARITY.md'


def _load_parity_map():
    spec = importlib.util.spec_from_file_location(
        'pymolx_parity_under_test', TESTING_DIR / 'pymolx_parity.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


parity = _load_parity_map()


def _requires_incentive(node):
    for deco in node.decorator_list:
        if (isinstance(deco, ast.Call) and
                getattr(deco.func, 'attr', None) == 'requires' and
                any(isinstance(arg, ast.Constant) and arg.value == 'incentive'
                    for arg in deco.args)):
            return True
    return False


def _definitions(path):
    '''
    Qualified names of all classes and methods in a test file, and the
    subset decorated with requires('incentive')
    '''
    tree = ast.parse(path.read_text())
    names, gated = set(), set()
    for node in tree.body:
        nodes = [(node, node.name)] if isinstance(
            node, (ast.ClassDef, ast.FunctionDef)) else []
        if isinstance(node, ast.ClassDef):
            nodes += [(child, node.name + '.' + child.name)
                      for child in node.body
                      if isinstance(child, ast.FunctionDef)]
        for child, qualname in nodes:
            names.add(qualname)
            if _requires_incentive(child):
                gated.add(qualname)
    return names, gated


def _gated_test_ids():
    for path in sorted(TESTS_DIR.rglob('*.py')):
        _, gated = _definitions(path)
        relpath = path.relative_to(TESTS_DIR).as_posix()
        for qualname in sorted(gated):
            yield relpath + '::' + qualname


# branding (H-01..H-05)

def test_version_message():
    message = pymol.get_version_message()
    assert message.startswith(branding.NAME + ' ' + pymolx.__version__)
    assert cmd.get_version()[0] in message


def test_notices():
    assert 'Not affiliated' in branding.DISCLAIMER
    assert 'Copyright' in branding.UPSTREAM_COPYRIGHT


def test_splash_installed():
    assert Path(cmd.exp_path(branding.SPLASH_PNG)).is_file()


# parity item registry

def test_features_register(monkeypatch):
    monkeypatch.setattr(features, '_implemented', set())
    assert not features.provides('P-04')
    features.register('P-04', 'F-03')
    assert features.provides('P-04')
    assert features.implemented() == ['F-03', 'P-04']


@pytest.mark.parametrize('bad_id', ['P4', 'p-04', 'P-004', 'desaturate'])
def test_features_register_invalid(bad_id):
    with pytest.raises(ValueError):
        features.register(bad_id)


# startup defaults (L-04)

def test_defaults_apply_revert(monkeypatch):
    monkeypatch.setattr(defaults, 'DEFAULTS', {'sphere_scale': 0.5})
    monkeypatch.setattr(defaults, '_replaced', {})
    original = cmd.get_setting_float('sphere_scale')
    assert original != 0.5

    defaults.apply(cmd)
    assert cmd.get_setting_float('sphere_scale') == 0.5

    assert defaults.revert(cmd)
    assert cmd.get_setting_float('sphere_scale') == original
    assert not defaults.revert(cmd)


# "pymolx" command (D-03)

def test_pymolx_command():
    result = cmd.pymolx(quiet=1)
    assert result['version'] == pymolx.__version__
    assert result['features'] == features.implemented()
    assert set(result['modules']) >= {'numpy', 'rdkit', 'gemmi', 'scipy'}
    assert result['modules']['numpy']


# parity map (T-01)

def test_parity_ids_are_in_inventory():
    inventory = INVENTORY.read_text()
    for item in parity.PARITY_TESTS:
        assert re.match(r'^[A-Z]-\d\d$', item)
        assert '| %s |' % item in inventory, item


def test_parity_map_entries_exist():
    for item, entries in parity.PARITY_TESTS.items():
        for entry in entries:
            relpath, _, qualname = entry.partition('::')
            path = TESTS_DIR / relpath
            assert path.is_file(), (item, entry)
            if qualname:
                names, _ = _definitions(path)
                assert qualname in names, (item, entry)


def test_every_incentive_test_is_mapped():
    unmapped = [test_id for test_id in _gated_test_ids()
                if not parity.required_items(test_id)]
    assert not unmapped, (
        'add these to testing/pymolx_parity.py: %s' % unmapped)


def test_required_items():
    assert parity.required_items(
        'api/computing.py::TestComputing.testClean') == {'P-01'}
    assert parity.required_items(
        'properties/test_simple.py::TestProperties.testMAEsaveLoadSessions'
    ) == {'F-01', 'F-02'}
    assert parity.required_items('api/computing.py::TestOther') == set()


def test_unlocked():
    test_id = 'properties/mae_io.py::TestMaePropertiesIO'
    assert not parity.unlocked(test_id, {'F-01'}.__contains__)
    assert parity.unlocked(test_id, {'F-01', 'F-02'}.__contains__)
    assert not parity.unlocked('api/unmapped.py::Test', lambda item: True)


# startup text in the viewer, like Incentive PyMOL

def test_startup_text():
    cmd.set('text', 0)
    pymolx.show_startup_text(cmd)
    assert cmd.get_setting_boolean('text')
    pymolx.hide_startup_text(cmd)
    assert not cmd.get_setting_boolean('text')
    # only once: later "text" changes are the user's
    cmd.set('text', 1)
    pymolx.hide_startup_text(cmd)
    assert cmd.get_setting_boolean('text')
    cmd.set('text', 0)


def test_multisample_default():
    # smooth edges in the viewer: 4x multisampling unless "-E N" says
    # otherwise (Qt and GLUT fall back to none if unavailable)
    import copy
    from pymol import invocation
    assert invocation.options.multisample == 4
    for argv, samples in [(['pymol', '-E', '0'], 0), (['pymol', '-E', '8'], 8)]:
        options = copy.copy(invocation.options)
        invocation.parse_args(argv, options=options, restricted=1)
        assert options.multisample == samples
