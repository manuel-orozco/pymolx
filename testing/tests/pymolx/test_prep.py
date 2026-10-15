'''
Protein preparation workflow (pymolx.prep: prepare_structure), structure
diagnostics (pymolx.diagnostics: structure_diagnostics) and the Structure
Preparation dialog. Skipped without OpenMM and PDBFixer. Uses the CPU
platform, so no GPU is needed.
'''

import math
import os

import pytest

from pymol import cmd
from pymolx import diagnostics, mm, prep, undo

pytestmark = pytest.mark.skipif(not mm.available(),
                                reason='needs openmm and pdbfixer')

SEQUENCE = 'ACDEFGHIKLMNPQRSTVWY'
DATA = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
FAST = dict(platform='CPU', diagnostics=0, quiet=1)


@pytest.fixture(autouse=True)
def session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


def peptide(name='pep', sequence=SEQUENCE):
    cmd.fab(sequence, name, ss=1)
    cmd.remove('%s and hydro' % name)


def gapped(tmp_path, cut=range(8, 12), name='loop'):
    '''
    The peptide without residues "cut", loaded from a PDB file whose
    REMARK 465 lists them (PyMOL keeps them without coordinates)
    '''
    peptide('full')
    resn = {}
    cmd.iterate('full and name CA', 'resn_[int(resi)] = resn',
                space={'resn_': resn})
    cmd.remove('full and resi %s' % '+'.join(map(str, cut)))
    lines = ['REMARK 465 MISSING RESIDUES',
             'REMARK 465   M RES C SSSEQI']
    lines += ['REMARK 465     %s A %5d' % (resn[i], i) for i in cut]
    cmd.alter('full', 'chain="A"')
    path = str(tmp_path / ('%s.pdb' % name))
    with open(path, 'w') as handle:
        handle.write('\n'.join(lines) + '\n' + cmd.get_pdbstr('full'))
    cmd.delete('full')
    cmd.load(path, name)
    return name


def residues(selection):
    found = []
    cmd.iterate('(%s) and name CA' % selection, 'found.append(resn)',
                space={'found': found})
    return found


def problems(selection):
    return [(i.category, i.residue, i.detail)
            for i in diagnostics.run(selection)
            if i.category in ('Chirality', 'Cis peptides', 'Bond lengths',
                              'Missing atoms')]


# the workflow

def test_missing_loop_filled_from_the_structure(tmp_path):
    name = gapped(tmp_path)
    assert cmd.count_atoms('%s and not present' % name) == 4
    report = cmd.prepare_structure(name, **FAST)
    assert report['object'] == 'loop_prep'
    assert report['missing_residues'] == 4
    assert residues('loop_prep') == residues(name)  # 20, in order
    # the built loop: L-amino acids, trans peptides, sensible bonds
    assert problems('loop_prep') == []
    # input atoms stay close (restrained); built ones were free
    assert report['rmsd'] < 0.5
    assert cmd.count_atoms('loop_prep and hydro') > 0
    assert not cmd.get_object_list('enabled and %s' % name)


def test_loop_left_out_and_chain_ends(tmp_path):
    name = gapped(tmp_path, cut=[1, 2, 10])
    # only the gap inside the chain by default
    cmd.prepare_structure(name, 'inner', minimize=0, protonate=0, **FAST)
    assert len(residues('inner')) == 18
    cmd.prepare_structure(name, 'whole', fill_termini=1, minimize=0,
                          protonate=0, **FAST)
    assert len(residues('whole')) == 20
    cmd.prepare_structure(name, 'nofill', fill_loops=0, minimize=0,
                          protonate=0, **FAST)
    assert len(residues('nofill')) == 17


def test_caps():
    peptide()
    cmd.prepare_structure('pep', cap_termini=1, max_iterations=200, **FAST)
    names = []
    cmd.iterate('pep_prep', 'names.append(resn)', space={'names': names})
    assert names[0] == 'ACE' and names[-1] == 'NME'
    assert cmd.count_atoms('pep_prep and name OXT') == 0
    # part of the chain, bonded to it
    assert cmd.count_atoms('pep_prep and resn ACE and name C and '
                           'neighbor (resn ALA and name N)') == 1
    assert cmd.count_atoms('pep_prep and resn NME and name N and '
                           'neighbor (resn TYR and name C)') == 1
    assert problems('pep_prep') == []


def test_existing_hydrogens_deleted_or_kept():
    peptide()
    cmd.prepare_structure('pep', 'once', max_iterations=100, **FAST)
    moved = 'once and resi 9 and name HZ1'  # an OpenMM name: known
    cmd.alter_state(1, moved, 'x = x + 0.3')
    before = cmd.get_coords(moved)[0].tolist()
    std = dict(FAST, pka_method='standard', minimize=0)
    cmd.prepare_structure('once', 'kept', delete_hydrogens=0, **std)
    assert cmd.get_coords('kept and resi 9 and name HZ1')[0].tolist() == \
        pytest.approx(before, abs=1e-3)
    cmd.prepare_structure('once', 'redone', **std)
    assert cmd.get_coords('redone and resi 9 and name HZ1')[0].tolist() != \
        pytest.approx(before, abs=1e-3)
    assert cmd.count_atoms('kept and hydro') == \
        cmd.count_atoms('redone and hydro')


def test_kept_hydrogens_with_pymol_names(capsys):
    # fab's hydrogen names (e.g. "1HB") aren't OpenMM's: added again
    cmd.fab(SEQUENCE, 'pep', ss=1)
    report = cmd.prepare_structure('pep', delete_hydrogens=0,
                                   max_iterations=100, platform='CPU',
                                   diagnostics=0)
    assert report['energy_after'] < report['energy_before']
    out = capsys.readouterr().out
    assert 'existing hydrogens kept:' in out
    assert "names OpenMM doesn't know replaced" in out
    cmd.prepare_structure('pep', 'fresh', minimize=0, platform='CPU',
                          diagnostics=0)
    assert 'existing hydrogens deleted: %d' % cmd.count_atoms(
        'pep and hydro') in capsys.readouterr().out


def test_side_chains_filled():
    peptide()
    cmd.remove('pep and resn LYS+TRP and not name N+CA+C+O+CB')
    report = cmd.prepare_structure('pep', minimize=0, **FAST)
    assert report['missing_atoms'] == 4 + 9  # Lys, Trp
    assert problems('pep_prep') == []
    cmd.prepare_structure('pep', 'kept', fill_side_chains=0, minimize=0,
                          protonate=0, **FAST)
    assert ('Missing atoms', '-/LYS 9', 'missing CD CE CG NZ') in [
        (c, r, d) for c, r, d in problems('kept')]


def variants(report, resn):
    return report['variants'].get(resn, {})


def test_protonation_follows_ph():
    peptide(sequence='GDEGHKGDEGHKG')
    std = dict(FAST, pka_method='standard')
    low = cmd.prepare_structure('pep', 'low', ph=2.0, minimize=0, **std)
    assert variants(low, 'ASP') == {'ASH': 2}
    assert variants(low, 'GLU') == {'GLH': 2}
    assert variants(low, 'HIS') == {'HIP': 2}
    high = cmd.prepare_structure('pep', 'high', ph=12.0, minimize=0,
                                 **std)
    assert variants(high, 'LYS') == {'LYN': 2}
    assert set(variants(high, 'HIS')) <= {'HID', 'HIE'}
    # His and Lys +1, Asp and Glu neutral, termini +1 and -1 (the force
    # fields have charged termini only)
    assert low['charge_estimate'] == 2 + 2 + 1 - 1
    assert cmd.count_atoms('low and resn ASP and name HD2') == 2


@pytest.mark.skipif(not prep.propka_available(), reason='needs propka')
def test_propka_pka_values_reported():
    peptide(sequence='GDEGHKGDEGHKG')
    report = cmd.prepare_structure('pep', minimize=0, **FAST)
    pkas = report['pka']
    assert {key[2] for key in pkas} == {'ASP', 'GLU', 'HIS', 'LYS'}
    asp = [v for key, v in pkas.items() if key[2] == 'ASP']
    assert all(2.0 < v < 6.0 for v in asp)
    standard = cmd.prepare_structure('pep', 'std', pka_method='standard',
                                     minimize=0, **FAST)
    assert standard.get('pka') == {}


def ligand_complex():
    '''
    A peptide, a water and a two-atom "ligand" (chain A, like PDB files)
    '''
    peptide()
    cmd.alter('pep', 'chain="A"')
    cmd.pseudoatom('w', pos=[0, 0, 30], name='O', resn='HOH', resi=100,
                   chain='A', elem='O')
    cmd.pseudoatom('l1', pos=[0, 30, 0], name='C1', resn='LIG', resi=200,
                   chain='A', elem='C')
    cmd.pseudoatom('l2', pos=[1.5, 30, 0], name='O1', resn='LIG',
                   resi=200, chain='A', elem='O')
    cmd.create('cpx', 'pep or w or l1 or l2')
    cmd.delete('pep w l1 l2')
    cmd.alter('cpx and resn HOH+LIG', 'type="HETATM"')
    cmd.bond('cpx and name C1', 'cpx and name O1', 2)
    return 'cpx'


def test_waters_and_ligands():
    name = ligand_complex()
    report = cmd.prepare_structure(name, max_iterations=100, **FAST)
    out = 'cpx_prep'
    # back in their chain, the ligand with its bond, unchanged
    assert cmd.count_atoms('%s and resn LIG and chain A' % out) == 2
    assert cmd.count_atoms('%s and resn HOH and chain A' % out) == 3
    assert cmd.count_atoms('%s and name C1 and neighbor name O1' % out) == 1
    assert cmd.get_coords('%s and name C1' % out).tolist() == [[0, 30, 0]]
    assert report['excluded_residues'] == ['LIG']

    cmd.prepare_structure(name, 'dry', delete_waters=1, minimize=0, **FAST)
    assert cmd.count_atoms('dry and resn HOH') == 0
    assert cmd.count_atoms('dry and resn LIG') == 2
    cmd.prepare_structure(name, 'apo', delete_heterogens=1, minimize=0,
                          **FAST)
    assert cmd.count_atoms('apo and resn LIG') == 0
    assert cmd.count_atoms('apo and resn HOH') == 3


def test_fix_structure_keeps_ligand_bonds():
    # (the older command, pymolx.mm, shares the bond restoring)
    name = ligand_complex()
    cmd.fix_structure(name, heterogens='none', quiet=1)
    assert cmd.count_atoms('cpx_fixed and name C1 and neighbor name O1') == 1


@pytest.mark.parametrize('forcefield, solvent', [
    ('amber19', 'obc2'), ('amber99sbildn', 'hct'),
    ('charmm36_2024', 'vacuum')])
def test_forcefields_and_solvents(forcefield, solvent):
    peptide()
    report = cmd.prepare_structure('pep', forcefield=forcefield,
                                   solvent=solvent, max_iterations=100,
                                   **FAST)
    assert report['energy_after'] < report['energy_before']


def test_explicit_solvent_is_removed():
    peptide(sequence='ACDEFGH')
    implicit = cmd.prepare_structure('pep', 'imp', minimize=0, **FAST)
    report = cmd.prepare_structure('pep', solvent='explicit',
                                   max_iterations=50, **FAST)
    assert report['atoms'] == implicit['atoms']
    assert cmd.count_atoms('pep_prep and resn HOH+NA+CL') == 0
    assert report['energy_after'] < report['energy_before']


def test_restraints_and_tolerance():
    peptide()
    heavy = cmd.prepare_structure('pep', 'h', max_iterations=200, **FAST)
    free = cmd.prepare_structure('pep', 'f', restrain='none',
                                 max_iterations=200, **FAST)
    assert heavy['restrained_atoms'] == cmd.count_atoms('pep')
    assert free['restrained_atoms'] == 0
    assert heavy['rmsd'] < free['rmsd']
    ca = cmd.prepare_structure('pep', 'c', restrain='ca',
                               max_iterations=50, **FAST)
    assert ca['restrained_atoms'] == 20


def test_options_checked():
    peptide()
    with pytest.raises(ValueError, match='no ACE/NME caps'):
        cmd.prepare_structure('pep', cap_termini=1,
                              forcefield='charmm36', solvent='vacuum',
                              **FAST)
    with pytest.raises(ValueError, match='implicit solvent'):
        cmd.prepare_structure('pep', forcefield='charmm36', **FAST)
    with pytest.raises(ValueError, match='solvent must be'):
        cmd.prepare_structure('pep', solvent='water', **FAST)
    with pytest.raises(TypeError):
        prep.prepare_job('pep', colour='red')
    assert 'pep_prep' not in cmd.get_names()


def test_cancel():
    peptide()
    job = prep.prepare_job('pep', platform='CPU')
    job.cancel()
    with pytest.raises(prep.Cancelled):
        job.run()


def test_one_undo_step():
    peptide()
    cmd.undo_enable()
    try:
        before = len(undo.stack.undo_stack)
        cmd.prepare_structure('pep', minimize=0, **FAST)
        assert len(undo.stack.undo_stack) == before + 1
        cmd.undo()
        assert 'pep_prep' not in cmd.get_names()
    finally:
        cmd.undo_disable()


def test_command_language_and_diagnostics_output(capsys):
    peptide()
    cmd.do('prepare_structure pep, ready, minimize=0, cap_termini=1, '
           'platform=CPU')
    cmd.sync()
    out = capsys.readouterr().out
    assert 'ready' in cmd.get_names()
    assert 'Protonation at pH 7.0' in out
    assert 'Diagnostics' in out and 'Hydrogens:' in out


def test_alternate_locations(tmp_path):
    peptide('pep', 'GSG')
    cmd.alter('pep', 'chain="A"')
    # Ser 2 side chain in two places, B more occupied
    cmd.create('alts', 'pep')
    cmd.alter('alts and resi 2 and name OG', 'alt="A"; q=0.3')
    cmd.create('other', 'alts and resi 2 and name OG')
    cmd.alter('other', 'alt="B"; q=0.7')
    cmd.translate([0, 0, 1.5], 'other')
    path = str(tmp_path / 'alt.pdb')
    cmd.save(path, 'alts or other')
    cmd.delete('all')
    cmd.load(path, 'two')
    assert cmd.count_atoms('two and name OG') == 2
    where = {k: cmd.get_coords('two and name OG and alt %s' % k)[0].tolist()
             for k in 'AB'}
    std = dict(FAST, minimize=0, protonate=0)
    report = cmd.prepare_structure('two', 'occ', **std)
    assert report['alternates'] == {'A/SER 2': 'B'}
    assert cmd.get_coords('occ and name OG')[0].tolist() == \
        pytest.approx(where['B'], abs=1e-3)
    cmd.prepare_structure('two', 'choice_a', altloc='A', **std)
    assert cmd.get_coords('choice_a and name OG')[0].tolist() == \
        pytest.approx(where['A'], abs=1e-3)
    assert cmd.count_atoms('choice_a and not alt ""') == 0
    with pytest.raises(ValueError, match='altloc'):
        prep.prepare_job('two', altloc='AB')


def test_selective_water_deletion():
    name = ligand_complex()
    # a second water next to the ligand, a third H-bonded to the peptide
    cmd.pseudoatom('w2', pos=[0, 32.8, 0], name='O', resn='HOH', resi=101,
                   chain='A', elem='O')
    o = cmd.get_coords('%s and resi 5 and name O' % name)[0]
    cmd.pseudoatom('w3', pos=list(map(float, o + [2.8, 0, 0])), name='O',
                   resn='HOH', resi=102, chain='A', elem='O')
    cmd.create('cpx2', '%s or w2 or w3' % name)
    cmd.alter('cpx2 and resn HOH', 'type="HETATM"')
    std = dict(FAST, minimize=0, protonate=0)

    def waters(obj):
        found = []
        cmd.iterate('%s and resn HOH' % obj, 'found.append(int(resi))',
                    space={'found': found})
        return sorted(found)

    cmd.prepare_structure('cpx2', 'near', delete_waters=1,
                          water_distance=4.0, **std)
    assert waters('near') == [101]           # 100 and 102 are far
    cmd.prepare_structure('cpx2', 'partners', water_min_hbonds=1, **std)
    assert waters('partners') == [101, 102]  # 100 has no partner
    cmd.prepare_structure('cpx2', 'dry', delete_waters=1, **std)
    assert waters('dry') == []


def ccd_reachable():
    import urllib.request
    try:
        urllib.request.urlopen(
            'https://files.rcsb.org/ligands/download/GOL.cif', timeout=10)
        return True
    except Exception:
        return False


@pytest.mark.skipif(not ccd_reachable(), reason='needs the PDB (internet)')
def test_ligand_hydrogens_from_the_ccd():
    cmd.load(os.path.join(DATA, '1oky.pdb.gz'))
    report = cmd.prepare_structure('1oky', minimize=0, optimize_hbonds=0,
                                   **FAST)
    assert report['ligand_hydrogens'] == {'GOL': 8 * 3, 'SO4': 0, 'STU': 26}
    gol = cmd.get_model('1oky_prep and resn GOL and resi 1360')
    assert sum(a.symbol == 'H' for a in gol.atom) == 8
    # bonded to their atoms
    assert cmd.count_atoms('1oky_prep and resn GOL and hydro and '
                           'not neighbor (resn GOL and not hydro)') == 0
    cmd.prepare_structure('1oky', 'bare', minimize=0, ligand_hydrogens=0,
                          optimize_hbonds=0, **FAST)
    assert cmd.count_atoms('bare and resn GOL+STU and hydro') == 0


def hydrogen_bonds(obj):
    '''
    Polar hydrogens 1.4 to 2.5 A from an N or O (not their own atom)
    '''
    import math
    pairs = cmd.find_pairs('%s and hydro and neighbor elem N+O' % obj,
                           '%s and elem N+O' % obj, cutoff=2.5)
    coords = {}
    cmd.iterate_state(1, obj, 'coords[index] = (x, y, z)',
                      space={'coords': coords})
    return sum(1 for (_, a), (_, b) in pairs
               if math.dist(coords[a], coords[b]) > 1.4)


def test_hydrogen_bond_optimization():
    cmd.load(os.path.join(DATA, '1oky.pdb.gz'))
    std = dict(FAST, minimize=0, ligand_hydrogens=0)
    report = cmd.prepare_structure('1oky', 'opt', **std)
    # (residues next to built ones may vary: PDBFixer builds randomly)
    assert set(report['flips']) >= {
        'Asn A/122 flipped', 'His A/203 flipped', 'Asn A/210 flipped',
        'Gln A/247 flipped'}
    # built residues aren't flipped (Gln 236: in the built loop)
    assert 'Gln A/236 flipped' not in report['flips']
    assert report['rotated'] > 50
    cmd.prepare_structure('1oky', 'plain', optimize_hbonds=0, **std)
    # flipped: the amide O and N swapped places
    assert cmd.get_coords('opt and resi 122 and name OD1')[0].tolist() == \
        pytest.approx(cmd.get_coords('plain and resi 122 and name ND2')[0]
                      .tolist(), abs=1e-3)
    # neutral His tautomers from their partners (Phe 141 N-H to ND1)
    assert cmd.count_atoms('opt and resi 139 and name HE2') == 1
    assert cmd.count_atoms('opt and resi 139 and name HD1') == 0
    assert hydrogen_bonds('opt') > hydrogen_bonds('plain')


def test_flip_scores():
    from pymolx import hbonds
    D, A, B = hbonds.DONOR, hbonds.ACCEPTOR, hbonds.BOTH
    assert hbonds.pair_score(D, A, 2.9) == 1
    assert hbonds.pair_score(A, D, 3.45) == pytest.approx(0.5)
    assert hbonds.pair_score(D, D, 3.0) == -1
    assert hbonds.pair_score(A, A, 3.0) == -1
    assert hbonds.pair_score(B, A, 3.0) == 1
    assert hbonds.pair_score(D, A, 2.2) == -1   # too close
    assert hbonds.pair_score(D, D, 3.4) == 0


def test_disulfides():
    cmd.load(os.path.join(DATA, '1v5a-3models.cif'), 'nmr')
    std = dict(FAST, state=1, minimize=0)
    report = cmd.prepare_structure('nmr', 'bridged', **std)
    assert report['disulfides'] == 3
    assert report['variants']['CYS'] == {'CYX': 6}

    def ss_bonds(obj):
        # (PyMOL's "neighbor" doesn't report these bonds)
        return len(cmd.get_model('%s and name SG' % obj).bond)

    assert ss_bonds('bridged') == 3
    report = cmd.prepare_structure('nmr', 'thiols', disulfides=0, **std)
    assert report['disulfides'] == 0
    assert ss_bonds('thiols') == 0
    assert cmd.count_atoms('thiols and resn CYS and name HG') == 6


def zinc_site(*bound):
    '''
    A peptide with a zinc ion 2.1 A from the given atoms' far side, e.g.
    ("6", "NE2") for His 6 bound through NE2
    '''
    import numpy as np
    peptide('pep', 'GAGAGHGAGCGAG')
    cmd.alter('pep', 'chain="A"')
    tips = []
    for resi, name in bound:
        atom = cmd.get_coords('pep and resi %s and name %s' % (resi, name))[0]
        neighbors = cmd.get_coords('pep and resi %s and neighbor name %s' % (
            resi, name))
        out = atom - neighbors.mean(axis=0)
        tips.append(atom + 2.1 * out / np.linalg.norm(out))
    where = np.mean(tips, axis=0) if len(tips) == 1 else tips[0]
    cmd.pseudoatom('zn', pos=list(map(float, where)), name='ZN', resn='ZN',
                   resi=401, chain='A', elem='Zn')
    cmd.create('site', 'pep or zn')
    cmd.alter('site and resn ZN', 'type="HETATM"; formal_charge=2')
    cmd.delete('pep zn')
    return 'site'


def test_metal_site():
    name = zinc_site(('6', 'NE2'))
    report = cmd.prepare_structure(name, max_iterations=200, **FAST)
    assert report['metal_sites'] == {'ZN A/401': 1}
    assert report['metal_fixed_atoms'] == 2
    assert report['metal_dashes'] == 'site_prep_metal'
    assert cmd.get_type('site_prep_metal') == 'object:measurement'
    out = 'site_prep'
    # His bound through NE2: protonated on ND1
    assert cmd.count_atoms('%s and resi 6 and name HD1' % out) == 1
    assert cmd.count_atoms('%s and resi 6 and name HE2' % out) == 0
    # a zero-order bond, the distance kept
    model = cmd.get_model('%s and (resn ZN or (resi 6 and name NE2))' % out)
    assert [b.order for b in model.bond] == [0]
    d = cmd.get_distance('%s and resn ZN' % out,
                         '%s and resi 6 and name NE2' % out)
    assert d == pytest.approx(2.1, abs=0.01)

    cmd.prepare_structure(name, 'plain', metals=0, minimize=0, **FAST)
    assert cmd.get_model('plain and (resn ZN or (resi 6 and name NE2))'
                         ).bond == []


def test_metal_bound_cysteine_is_a_thiolate():
    name = zinc_site(('10', 'SG'))
    report = cmd.prepare_structure(name, minimize=0, **FAST)
    assert report['variants']['CYS'] == {'CYM': 1}
    assert cmd.count_atoms('site_prep and resi 10 and name HG') == 0


def test_parse_mutations():
    assert prep.parse_mutations('A:45:GLY, A/102/w;7:Y') == [
        ('A', '45', 'GLY'), ('A', '102', 'TRP'), ('', '7', 'TYR')]
    assert prep.parse_mutations('') == []
    for bad, message in [('A:45', 'chain:residue:new'),
                         ('A:45:XYZ', 'not a standard'),
                         ('A:45A:GLY', 'insertion codes'),
                         ('45', 'chain:residue:new')]:
        with pytest.raises(ValueError, match=message):
            prep.parse_mutations(bad)


def test_mutations():
    peptide()
    report = cmd.prepare_structure('pep', mutations='5:W 6:ALA',
                                   max_iterations=200, **FAST)
    assert report['mutations'] == 2
    assert residues('pep_prep')[4:6] == ['TRP', 'ALA']
    assert cmd.count_atoms('pep_prep and resi 5 and not hydro') == 14
    assert problems('pep_prep') == []
    # side chains are built even with "fill in missing side chains" off
    cmd.prepare_structure('pep', 'mut2', mutations='5:W', minimize=0,
                          fill_side_chains=0, **FAST)
    assert cmd.count_atoms('mut2 and resi 5 and not hydro') == 14
    with pytest.raises(ValueError, match='no such residue'):
        cmd.prepare_structure('pep', mutations='99:G', **FAST)
    with pytest.raises(ValueError, match='Cannot|actually|KeyError|'
                       'no such'):
        cmd.prepare_structure('pep', mutations='B:5:G', **FAST)


# diagnostics

def categories(selection):
    return diagnostics.summary(diagnostics.run(selection))


def test_clean_peptide_has_no_issues():
    # (fab's side chains of other residues may clash)
    cmd.fab('AAAAAAAAAAAA', 'pep', ss=1)
    counts = categories('pep')
    assert sum(counts.values()) == 0, diagnostics.format_report(
        diagnostics.run('pep'))


def test_diagnostics_find_problems(tmp_path):
    name = gapped(tmp_path)
    cmd.alter(name, 'segi=""')
    issues = diagnostics.run(name)
    missing = [i for i in issues if i.category == 'Missing residues']
    assert [i.residue for i in missing] == ['A/8-11']  # no chain break too
    assert cmd.count_atoms(missing[0].selection) > 0

    peptide()
    cmd.remove('pep and resi 3 and name CG+OD1')         # missing atoms
    cmd.remove('pep and resi 15')                        # chain break
    cmd.set_dihedral('pep and resi 9 and name CA', 'pep and resi 9 and name C',
                     'pep and resi 10 and name N',
                     'pep and resi 10 and name CA', 0)   # cis peptide
    cmd.alter_state(1, 'pep and resi 12 and name CB',     # D-amino acid
                    '(x, y, z) = mirror(x, y, z)', space={
                        'mirror': mirror_through(cmd, 'pep and resi 12')})
    by = {}
    for issue in diagnostics.run('pep'):
        by.setdefault(issue.category, []).append(issue)
    assert by['Missing atoms'][0].detail == 'missing CG OD1'
    assert by['Missing residues'][0].detail.startswith('chain break')
    assert [i.residue for i in by['Cis peptides']] == ['-/LYS 9 - -/LEU 10']
    assert [i.residue for i in by['Chirality']] == ['-/ASN 12']
    assert by['Hydrogens'][0].severity == 'warning'  # none


def mirror_through(_self, residue):
    '''
    Reflection through the residue's N-CA-C plane, for alter_state
    '''
    n, ca, c = (_self.get_coords('(%s) and name %s' % (residue, name))[0]
                for name in ('N', 'CA', 'C'))
    u, v = n - ca, c - ca
    normal = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
              u[0] * v[1] - u[1] * v[0]]
    length = math.sqrt(sum(x * x for x in normal))
    normal = [x / length for x in normal]

    def mirror(x, y, z):
        d = sum((p - q) * k for p, q, k in zip((x, y, z), ca, normal))
        return tuple(p - 2 * d * k for p, k in zip((x, y, z), normal))
    return mirror


def test_diagnostics_clashes_and_bonds():
    peptide('pep', 'AAAAA')
    peptide('gly', 'GG')
    # a glycine pushed into the peptide
    target = cmd.get_coords('pep and resi 3 and name CB')[0]
    here = cmd.get_coords('gly and resi 1 and name CA')[0]
    cmd.translate(list(target - here + [0.6, 0, 0]), 'gly')
    cmd.create('both', 'pep or gly')
    cmd.alter('both and resn GLY', 'chain="B"')
    clashes = [i for i in diagnostics.run('both')
               if i.category == 'Steric clashes']
    assert clashes and clashes[0].severity == 'error'
    # a stretched bond (CB moved 1 A further out)
    ca, cb = (cmd.get_coords('pep and resi 2 and name ' + n)[0]
              for n in ('CA', 'CB'))
    shift = (cb - ca) / math.dist(cb, ca)
    cmd.alter_state(1, 'pep and resi 2 and name CB',
                    '(x, y, z) = (x + dx, y + dy, z + dz)',
                    space=dict(zip(('dx', 'dy', 'dz'), map(float, shift))))
    bonds = [i for i in diagnostics.run('pep')
             if i.category == 'Bond lengths']
    assert [i.detail.split()[0] for i in bonds] == ['CA-CB']


def test_molprobity_outliers_match_wwpdb():
    # wwPDB validation report of 1OKY (MolProbity): Ramachandran outliers
    # Ser 231, Thr 359; rotamer outliers as below
    cmd.load(os.path.join(DATA, '1oky.pdb.gz'))
    issues = diagnostics.run('1oky')
    rama = [i.residue for i in issues
            if i.category == 'Ramachandran outliers']
    rota = [i.residue for i in issues if i.category == 'Rotamer outliers']
    assert rama == ['A/SER 231', 'A/THR 359']
    assert rota == ['A/GLU 80', 'A/PHE 93', 'A/VAL 97', 'A/THR 180',
                    'A/GLU 256', 'A/LEU 289', 'A/GLU 328', 'A/GLU 343']
    scores = [i.detail for i in issues if i.category == 'Scores']
    assert 'outliers 0.7%' in scores[1] and '3.3% (8 of 245)' in scores[2]
    assert 'no hydrogens, not comparable' in scores[0]


def test_clashscore_and_cbeta():
    cmd.fab('AAAAAAAA', 'pep', ss=1)
    scores = {i.detail.split()[0]: i.detail
              for i in diagnostics.run('pep') if i.category == 'Scores'}
    assert scores['clashscore'].startswith('clashscore 0.0 ')
    assert scores['C-beta'] == 'C-beta deviations: 0'
    cb = cmd.get_coords('pep and resi 4 and name CB')[0]
    cmd.alter_state(1, 'pep and resi 4 and name CB', 'x = x + 0.5')
    found = [i for i in diagnostics.run('pep')
             if i.category == 'C-beta deviations']
    assert [i.residue for i in found] == ['-/ALA 4']
    assert cb is not None


def test_coarse_backbone_check_without_data(monkeypatch):
    monkeypatch.setattr(diagnostics, '_top8000', {})
    cmd.fab('AAAAA', 'pep', ss=1)
    cmd.set_dihedral('pep and resi 2 and name C', 'pep and resi 3 and name N',
                     'pep and resi 3 and name CA',
                     'pep and resi 3 and name C', 150)  # phi > 0
    found = [i.residue for i in diagnostics.run('pep')
             if i.category == 'Ramachandran outliers']
    assert found == ['-/ALA 3']


def test_diagnostics_command(capsys):
    cmd.fab('AAAA', 'pep', ss=1)
    issues = cmd.structure_diagnostics('pep')
    assert isinstance(issues, list)
    assert 'Steric clashes: none' in capsys.readouterr().out


# dialog

@pytest.fixture
def dialog():
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import pymol.plugins  # noqa: F401 (adds the startup plugin path)
    from pmg_tk.startup.structure_prep_gui.main import StructurePrepDialog
    d = StructurePrepDialog(None, _self=cmd)
    yield d
    d.deleteLater()


def test_dialog_runs_the_workflow(dialog, tmp_path):
    name = gapped(tmp_path)
    dialog.refresh()
    dialog.selection.setEditText(name)
    dialog.cap_termini.setChecked(True)
    dialog.platform.setCurrentText('CPU')
    dialog.max_iterations.setValue(300)
    options = dialog.options()
    assert options['cap_termini'] and options['fill_loops']
    assert options['delete_hydrogens']
    assert options['disulfides'] and options['metals']
    dialog.mutations.setText('A:9:W')
    assert dialog.options()['mutations'] == [('A', '9', 'TRP')]
    dialog.mutations.setText('')
    assert options['forcefield'] == 'amber14' and options['solvent'] == 'gbn2'
    assert options['restrain'] == 'heavy' and options['tolerance'] == 10

    dialog.run()
    assert dialog.thread is not None
    assert not dialog.run_button.isEnabled()
    dialog.thread.wait()
    dialog._finished()  # normally delivered by the event loop
    assert dialog.progress.text() == 'Done: "loop_prep"'
    assert 'loop_prep' in cmd.get_names()
    log = dialog.log.toPlainText()
    for step in ('Preprocess', 'Protonation at pH 7.0', 'Minimization',
                 'New object "loop_prep"', 'Diagnostics:'):
        assert step in log
    # diagnostics shown in their tab
    assert dialog.tabs.currentIndex() == 1
    tops = [dialog.tree.topLevelItem(i).text(0)
            for i in range(dialog.tree.topLevelItemCount())]
    assert tops[0] == 'Scores'  # its lines below it
    assert dialog.tree.topLevelItem(0).childCount() >= 3
    assert any(t.startswith('Rotamer outliers (') for t in tops)
    assert dialog.diag_summary.text().startswith('loop_prep:')
    assert dialog.run_button.isEnabled()


def test_dialog_steps_off(dialog):
    peptide()
    dialog.refresh()
    dialog.selection.setEditText('pep')
    dialog.new_name.setText('plain')
    dialog.minimize.setChecked(False)
    dialog.protonate.setChecked(False)
    dialog.run_diagnostics.setChecked(False)
    dialog.run()
    dialog.thread.wait()
    dialog._finished()
    assert cmd.count_atoms('plain and hydro') == 0
    assert dialog.tabs.currentIndex() == 0


def test_dialog_compatibility(dialog):
    dialog.cap_termini.setChecked(True)
    dialog.forcefield.setCurrentIndex(dialog.forcefield.findData('charmm36'))
    # no caps, no implicit solvent with CHARMM
    assert not dialog.cap_termini.isChecked()
    assert not dialog.cap_termini.isEnabled()
    assert dialog.solvent.currentData() == 'explicit'
    model = dialog.solvent.model()
    assert not model.item(dialog.solvent.findData('gbn2')).isEnabled()
    # caps are fine without the minimization
    dialog.minimize.setChecked(False)
    assert dialog.cap_termini.isEnabled()
    dialog.minimize.setChecked(True)
    dialog.forcefield.setCurrentIndex(dialog.forcefield.findData('amber14'))
    assert dialog.cap_termini.isEnabled()
    assert model.item(dialog.solvent.findData('gbn2')).isEnabled()


def test_dialog_diagnostics_tab(dialog):
    cmd.fab('AADAAAAAAA', 'pep', ss=1)
    cmd.remove('pep and resi 3 and name CG+OD1+OD2')
    dialog.refresh()
    dialog.diag_selection.setEditText('pep')
    dialog.diagnose()
    assert dialog.diag_summary.text() == 'pep: 1 issue'
    tree = dialog.tree
    missing = [tree.topLevelItem(i) for i in range(tree.topLevelItemCount())
               if tree.topLevelItem(i).text(0).startswith('Missing atoms')]
    assert missing[0].text(0) == 'Missing atoms (1)'
    dialog._show_issue(missing[0].child(0))
    assert cmd.count_atoms('sele') == cmd.count_atoms('pep and resi 3')
