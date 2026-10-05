'''
Render data/pymolx/splash.png (parity item H-01).

The molecule is ray traced by pymolx itself, the text is drawn with
Pillow. Run from the source root:

    pymol -cq tools/pymolx/make_splash.py
'''

import os
import tempfile

from PIL import Image, ImageDraw, ImageFont

from pymol import cmd
import pymolx.branding as branding

WIDTH, HEIGHT = 640, 480
MOL_WIDTH = 330
MOL_HEIGHT = 400  # leaves room for the footer text
BG_TOP = (24, 32, 46)
BG_BOTTOM = (10, 14, 22)
TEXT = (236, 240, 245)
TEXT_DIM = (150, 160, 175)
ACCENT = (86, 180, 233)

OUT = os.path.join('data', 'pymolx', 'splash.png')
MOLECULE = os.path.join('data', 'demo', 'il2.pdb')

FONT_DIRS = [
    '/usr/share/fonts/truetype/dejavu',
    '/Library/Fonts',
    'C:\\Windows\\Fonts',
]


def font(size, bold=False):
    name = 'DejaVuSans-Bold.ttf' if bold else 'DejaVuSans.ttf'
    for d in FONT_DIRS:
        path = os.path.join(d, name)
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default(size)


def render_molecule(filename):
    cmd.reinitialize()
    cmd.load(MOLECULE, 'mol')
    cmd.remove('solvent or hetatm')
    cmd.hide('everything')
    cmd.show('cartoon')
    cmd.spectrum('count', 'rainbow', 'mol and name CA')
    cmd.set('cartoon_fancy_helices')
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_trace_color', 'black')
    cmd.set('ray_opaque_background', 0)
    cmd.set('antialias', 2)
    cmd.orient('mol')
    cmd.turn('z', 90)  # long axis vertical, to fit the tall render box
    cmd.turn('y', -30)
    cmd.zoom('mol', buffer=2, complete=1)
    cmd.png(filename, width=MOL_WIDTH, height=MOL_HEIGHT, ray=1, quiet=1)


def gradient():
    img = Image.new('RGB', (WIDTH, HEIGHT))
    draw = ImageDraw.Draw(img)
    for y in range(HEIGHT):
        t = y / (HEIGHT - 1)
        draw.line([(0, y), (WIDTH, y)], fill=tuple(
            round(a + (b - a) * t) for a, b in zip(BG_TOP, BG_BOTTOM)))
    return img


def main():
    img = gradient()

    with tempfile.TemporaryDirectory() as tmpdir:
        mol_png = os.path.join(tmpdir, 'mol.png')
        render_molecule(mol_png)
        mol = Image.open(mol_png).convert('RGBA')
        img.paste(mol, (WIDTH - MOL_WIDTH - 10, 10), mol)

    draw = ImageDraw.Draw(img)
    x = 36

    draw.text((x, 120), branding.NAME, font=font(64, bold=True), fill=TEXT)
    draw.rectangle([x, 205, x + 64, 209], fill=ACCENT)
    draw.text((x, 225), 'Molecular graphics', font=font(20), fill=TEXT)
    draw.text((x, 255), 'version %s' % branding.version(), font=font(15),
              fill=TEXT_DIM)

    small = font(12)
    draw.text((x, HEIGHT - 58), branding.DISCLAIMER, font=small, fill=TEXT_DIM)
    draw.text((x, HEIGHT - 40), branding.UPSTREAM_COPYRIGHT, font=small,
              fill=TEXT_DIM)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    img.save(OUT, optimize=True)
    print(' Wrote', OUT)


main()
