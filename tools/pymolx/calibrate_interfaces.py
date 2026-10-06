'''
Calibrate pymolx.interfaces against PDBe PISA.

Downloads protein complexes (PDB format) and their PISA interface results,
runs pymolx.interfaces.analyze, and:

- fits the atomic solvation parameters (per ASP type) so that dG matches
  PISA's int_solv_en, cross-validated by entry (fit on half, test on half)
- finds the H-bond distance cutoff whose counts best match PISA's

Run from the source root (downloads are cached in the given directory):

    pymol -cq tools/pymolx/calibrate_interfaces.py -- CACHE_DIR
'''

import json
import os
import sys
import urllib.request
import xml.etree.ElementTree as ET

import numpy

from pymol import cmd
from pymolx import interfaces

# protein complexes (mostly from protein docking benchmarks), plus 1tii
ENTRIES = '''
    1avx 1ay7 1bvn 1cgi 1d6r 1dfj 1e6e 1eaw 1ewy 1f34 1fle 1gl1 1gxd 1jtg
    1mah 1oph 1ppe 1r0r 1tmq 1udi 1yvb 2b42 2mta 2o8v 2oul 2pcc 2sic 2sni
    2uuy 3sgq 4cpa 7cei 1brs 2ptc 4hhb 1tii 1a2k 1ak4 1atn 1bkd 1buh 1e96
    1f51 1fqj 1gcq 1h9d 1he1 1i4d 1j2j 1k74 1kac 1ktz 1ml0 1qa9 1s1q 1sbb
    1wq1 1xd3 1z0k 2a9k 2ajf 2btf 2hle 2hqs 3cph
'''.split()

TYPES = ['C', 'N', 'O', 'O-', 'N+', 'S']
PISA_URL = 'https://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?%s'


def pisa_interfaces(code, cache):
    path = os.path.join(cache, code + '_pisa.xml')
    if not os.path.exists(path):
        with urllib.request.urlopen(PISA_URL % code, timeout=120) as r:
            data = r.read()
        with open(path, 'wb') as handle:
            handle.write(data)
    result = {}
    for iface in ET.parse(path).getroot().iter('interface'):
        mols = iface.findall('molecule')
        if len(mols) != 2:
            continue
        if any(m.findtext('symop', '').strip().lower() != 'x,y,z'
               for m in mols):
            continue
        if any((m.findtext('class') or '').strip() != 'Protein'
               for m in mols):
            continue
        chains = tuple(sorted(m.findtext('chain_id').strip() for m in mols))
        hb = iface.find('h-bonds')
        result[chains] = {
            'area': float(iface.findtext('int_area')),
            'dg': float(iface.findtext('int_solv_en')),
            'hb': int(hb.findtext('n_bonds')) if hb is not None else 0,
        }
    return result


def collect(cache):
    rows = []
    for code in ENTRIES:
        try:
            pisa = pisa_interfaces(code, cache)
            cmd.reinitialize()
            # explicit path: reinitialize resets the fetch_path setting
            cmd.fetch(code, 'm', type='pdb', path=cache, quiet=1)
            ours = interfaces.analyze('m')
        except Exception as e:
            print(' skip %s: %s' % (code, e))
            continue
        for iface in ours:
            key = tuple(sorted(iface.chains))
            if key not in pisa:
                continue
            rows.append({
                'entry': code, 'chains': key,
                'pisa': pisa[key],
                'area': iface.area,
                'bsa': {t: iface.bsa_by_type.get(t, 0.0) for t in TYPES},
                'hb_dist': [d for _, _, d in iface.hbonds],
            })
        print(' %s: %d matched interfaces' % (code, sum(
            1 for r in rows if r['entry'] == code)))
    return rows


def design(rows):
    x = numpy.array([[r['bsa'][t] for t in TYPES] for r in rows])
    y = numpy.array([r['pisa']['dg'] for r in rows])
    return x, y


def fit(rows):
    # dG = -sum(sigma_t * bsa_t) / 1000  ->  solve for sigma (cal/mol/A^2)
    x, y = design(rows)
    sigma, *_ = numpy.linalg.lstsq(-x / 1000.0, y, rcond=None)
    return dict(zip(TYPES, sigma))


def predict(rows, sigma):
    x, _ = design(rows)
    return -x.dot([sigma[t] for t in TYPES]) / 1000.0


def stats(pred, actual):
    err = pred - actual
    ss_res = (err ** 2).sum()
    ss_tot = ((actual - actual.mean()) ** 2).sum()
    return {'mae': float(numpy.abs(err).mean()),
            'r2': float(1 - ss_res / ss_tot)}


def main(cache):
    os.makedirs(cache, exist_ok=True)
    rows = collect(cache)
    entries = sorted({r['entry'] for r in rows})
    train_entries = set(entries[0::2])
    train = [r for r in rows if r['entry'] in train_entries]
    test = [r for r in rows if r['entry'] not in train_entries]
    _, y_test = design(test)

    report = {'n_entries': len(entries), 'n_interfaces': len(rows)}
    report['eisenberg_test'] = stats(predict(test, interfaces.ASP_EISENBERG),
                                     y_test)
    sigma_train = fit(train)
    report['fitted_test'] = stats(predict(test, sigma_train), y_test)
    report['sigma_all'] = fit(rows)
    _, y_all = design(rows)
    report['fitted_all'] = stats(predict(rows, report['sigma_all']), y_all)

    areas = numpy.array([r['area'] for r in rows])
    pisa_areas = numpy.array([r['pisa']['area'] for r in rows])
    report['area'] = stats(areas, pisa_areas)
    report['area_rel_err_median'] = float(numpy.median(
        numpy.abs(areas - pisa_areas) / numpy.maximum(pisa_areas, 1)))

    pisa_hb = numpy.array([r['pisa']['hb'] for r in rows])
    report['hbond_cutoffs'] = {}
    for cutoff in numpy.arange(3.2, 3.95, 0.05):
        ours = numpy.array([sum(d <= cutoff for d in r['hb_dist'])
                            for r in rows])
        report['hbond_cutoffs']['%.2f' % cutoff] = {
            'mae': float(numpy.abs(ours - pisa_hb).mean()),
            'bias': float((ours - pisa_hb).mean())}

    print(json.dumps(report, indent=1, sort_keys=True))
    with open(os.path.join(cache, 'calibration.json'), 'w') as handle:
        json.dump(report, handle, indent=1, sort_keys=True)


main(sys.argv[-1] if len(sys.argv) > 1 else 'pisa_calibration')
