'''
Lighting presets and ambient occlusion, like ChimeraX's "lighting soft"
(Display > Lighting and Display > Occlusion).

Softer presets spread the light over more light sources, add ambient
light and dull the specular highlight. Ambient occlusion darkens buried
parts of surfaces (PyMOL computes it for surfaces only).
'''

from pymol import cmd
from pymol.shortcut import Shortcut

# all presets set the same settings, so they can be switched freely
PRESETS = {
    # PyMOL's own values
    'default': dict(light_count=2, spec_count=-1, specular=1.0,
                    shininess=55.0, ambient=0.14, direct=0.45, reflect=0.45),
    'soft': dict(light_count=4, spec_count=1, specular=0.3, shininess=20.0,
                 ambient=0.3, direct=0.35, reflect=0.35),
    'softer': dict(light_count=6, spec_count=1, specular=0.15,
                   shininess=10.0, ambient=0.4, direct=0.25, reflect=0.3),
    # no highlight at all: a matte, clay-like look
    'softest': dict(light_count=8, spec_count=1, specular=0.0,
                    shininess=10.0, ambient=0.5, direct=0.2, reflect=0.25),
}

OCCLUSION_MODE = 1        # ambient_occlusion_mode: by atom neighbors
OCCLUSION_STRENGTH = 15.0  # ambient_occlusion_scale (PyMOL's own: 25)

preset_sc = Shortcut(PRESETS)
onoff_sc = Shortcut(['on', 'off'])


def lighting(preset='default', quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    "lighting" applies a lighting preset. Softer presets spread the
    light over more light sources and dull the shine, like ChimeraX's
    "lighting soft".

USAGE

    lighting [ preset ]

ARGUMENTS

    preset = default, soft, softer or softest (no shine) {default: default}

EXAMPLE

    lighting softer
    occlusion on

SEE ALSO

    occlusion, set
    '''
    name = preset_sc.auto_err(str(preset).strip() or 'default',
                              'lighting preset')
    for setting, value in PRESETS[name].items():
        _self.set(setting, value, quiet=1)
    if not int(quiet):
        print(' Lighting: %s' % name)
    return name


def occlusion(state='on', strength=OCCLUSION_STRENGTH, quiet=1, *,
              _self=cmd):
    '''
DESCRIPTION

    "occlusion" turns ambient occlusion on or off: buried parts of
    surfaces (pockets, grooves, interfaces) get darker. PyMOL computes
    it for surfaces only, when they are built.

USAGE

    occlusion [ state [, strength ]]

ARGUMENTS

    state = on or off {default: on}

    strength = float: how dark buried parts get {default: 15}

SEE ALSO

    lighting, ambient_occlusion_mode
    '''
    state = onoff_sc.auto_err(str(state).strip() or 'on', 'occlusion state')
    if state == 'on':
        # the mode setting resets the scale: mode first
        _self.set('ambient_occlusion_mode', OCCLUSION_MODE, quiet=1)
        _self.set('ambient_occlusion_scale', float(strength), quiet=1)
    else:
        _self.set('ambient_occlusion_mode', 0, quiet=1)
    if not int(quiet):
        print(' Occlusion: %s' % state)
    return state == 'on'


def extend(_self):
    for func, sc in ((lighting, preset_sc), (occlusion, onoff_sc)):
        _self.extend(func.__name__, func)
        setattr(_self, func.__name__, func)
        _self.auto_arg[0][func.__name__] = [sc, func.__name__, '']
