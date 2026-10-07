'''
Lighting presets and ambient occlusion (Display > Lighting, Display >
Occlusion): the "lighting" and "occlusion" commands and their menus.
'''

import pytest

from pymol import cmd
from pymolx import lighting, undo


@pytest.fixture(autouse=True)
def clean_session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


def current(names):
    return {name: cmd.get_setting_float(name) for name in names}


@pytest.mark.parametrize('preset', sorted(lighting.PRESETS))
def test_presets(preset):
    values = lighting.PRESETS[preset]
    assert cmd.lighting(preset) == preset
    assert current(values) == pytest.approx(values)


def test_presets_get_softer():
    p = lighting.PRESETS
    order = ['default', 'soft', 'softer', 'softest']
    for name in ('light_count', 'ambient'):
        assert [p[n][name] for n in order] == sorted(p[n][name] for n in order)
    for name in ('specular', 'direct'):
        assert [p[n][name] for n in order] == sorted(
            (p[n][name] for n in order), reverse=True)
    assert p['softest']['specular'] == 0
    # every preset sets the same settings: switching never leaves any over
    assert len({frozenset(v) for v in p.values()}) == 1


def test_default_restores_pymol_lighting():
    names = lighting.PRESETS['default']
    stock = current(names)
    cmd.lighting('softest')
    assert current(names) != stock
    cmd.lighting('default')
    assert current(names) == pytest.approx(stock)


def test_command_language_and_errors():
    cmd.do('lighting softer')
    cmd.sync()
    assert cmd.get_setting_int('light_count') == 6
    assert cmd.lighting('softes') == 'softest'  # unique prefix
    with pytest.raises(Exception):
        cmd.lighting('soft-ish')
    with pytest.raises(Exception):
        cmd.lighting('softe')  # softer or softest
    assert 'lighting' in cmd.auto_arg[0]
    assert 'occlusion' in cmd.auto_arg[0]


def test_occlusion():
    assert cmd.occlusion() is True
    assert cmd.get_setting_int('ambient_occlusion_mode') == 1
    # set after the mode, which resets it to 25
    assert cmd.get_setting_float('ambient_occlusion_scale') == 15
    cmd.occlusion('on', 30)
    assert cmd.get_setting_float('ambient_occlusion_scale') == 30
    assert cmd.occlusion('off') is False
    assert cmd.get_setting_int('ambient_occlusion_mode') == 0


def test_occlusion_darkens_surfaces(tmp_path):
    # ray traced, so it runs headless
    from PIL import Image, ImageStat

    cmd.fab('ACDEFGHIKLMNPQRSTVWY', 'pep', ss=1)
    cmd.show_as('surface')
    cmd.color('white')
    cmd.orient()

    def brightness(name):
        path = str(tmp_path / name)
        cmd.png(path, width=160, height=120, ray=1)
        image = Image.open(path).convert('RGBA')
        mask = image.getchannel('A')  # the surface, not the background
        return ImageStat.Stat(image.convert('L'), mask).mean[0]

    plain = brightness('plain.png')
    cmd.occlusion('on')
    darker = brightness('occlusion.png')
    cmd.occlusion('off')
    assert brightness('off.png') == pytest.approx(plain, abs=1)
    assert darker < plain - 3


def test_one_undo_step_each():
    cmd.undo_enable()
    try:
        before = len(undo.stack.undo_stack)
        cmd.lighting('softer')
        cmd.occlusion('on')
        assert len(undo.stack.undo_stack) == before + 2
        cmd.undo()
        assert cmd.get_setting_int('ambient_occlusion_mode') == 0
        assert cmd.get_setting_int('light_count') == 6
        cmd.undo()
        assert cmd.get_setting_int('light_count') == 2
    finally:
        cmd.undo_disable()


def test_display_menus():
    from pymol._gui import PyMOLDesktopGUI

    class Gui(PyMOLDesktopGUI):
        def __getattr__(self, name):  # handlers only the Qt window has
            return None

    menus = {label: data for _, label, data in Gui().get_menudata(cmd)}
    display = {item[1]: item for item in menus['Display']
               if item[0] == 'menu'}
    items = {item[1]: item[2] for item in display['Lighting'][2]}
    assert items == {'Default': 'lighting default',
                     'Soft': 'lighting soft',
                     'Softer': 'lighting softer',
                     'Softest (Matte)': 'lighting softest'}
    items = {item[1]: item[2] for item in display['Occlusion'][2]}
    assert items == {'On (Surfaces)': 'occlusion on', 'Off': 'occlusion off'}
    # right after Quality
    names = [item[1] for item in menus['Display'] if item[0] == 'menu']
    assert names.index('Lighting') == names.index('Quality') + 1

    # the menu commands work
    for command in ['lighting softest', 'occlusion on']:
        cmd.do(command)
    cmd.sync()
    assert cmd.get_setting_float('specular') == 0
    assert cmd.get_setting_int('ambient_occlusion_mode') == 1
