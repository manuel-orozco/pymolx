'''
Build top8000.npz from the Richardson Lab's MolProbity reference data
(github.com/rlabduke/reference_data, folder Top8000, CC BY 4.0; see
README.md here).

USAGE

    python build_top8000.py path/to/reference_data/Top8000

Each table (rama8000-*, rota8000-*) becomes a dense float16 grid of
percentiles at the bin centers (missing rows: 0) plus its dimensions
(lower bound, upper bound, number of bins) per axis; all axes wrap.
'''

import os
import sys

import numpy as np

RAMA = {
    'general': 'rama8000-general-noGPIVpreP',
    'glycine': 'rama8000-gly-sym',
    'ileval': 'rama8000-ileval-nopreP',
    'prepro': 'rama8000-prepro-noGP',
    'transpro': 'rama8000-transpro',
    'cispro': 'rama8000-cispro',
}
ROTA = ['arg', 'asn', 'asp', 'cys', 'gln', 'glu', 'his', 'ile', 'leu',
        'lys', 'met', 'phetyr', 'pro', 'ser', 'thr', 'trp', 'val']


def read(path):
    dims = []
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.startswith('#'):
                parts = line[1:].split()
                if parts and parts[0].startswith('x') and \
                        parts[0].endswith(':'):
                    dims.append((float(parts[1]), float(parts[2]),
                                 int(parts[3])))
                continue
            if line.strip():
                rows.append([float(v) for v in line.split()])
    rows = np.array(rows)
    grid = np.zeros([n for _, _, n in dims], dtype=np.float32)
    index = tuple(np.floor((rows[:, k] - lo) / ((hi - lo) / n)).astype(int)
                  for k, (lo, hi, n) in enumerate(dims))
    grid[index] = rows[:, -1]
    return np.array(dims), grid


def main(top8000, out='top8000.npz'):
    arrays = {}
    for name, stem in RAMA.items():
        dims, grid = read(os.path.join(
            top8000, 'Top8000_ramachandran_pct_contour_grids',
            stem + '.data'))
        arrays['rama_' + name] = grid.astype(np.float16)
        arrays['rama_' + name + '_dims'] = dims
    for name in ROTA:
        dims, grid = read(os.path.join(
            top8000, 'Top8000_rotamer_pct_contour_grids',
            'rota8000-%s.data' % name))
        arrays['rota_' + name] = grid.astype(np.float16)
        arrays['rota_' + name + '_dims'] = dims
    np.savez_compressed(out, **arrays)


if __name__ == '__main__':
    main(*sys.argv[1:])
