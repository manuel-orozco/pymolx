'''
Antibody builder (pymolx.antibody: build_antibody) and the Antibody
Builder dialog. Structures come from the PDB (skipped offline): 1HZH
(full human IgG1), 1N8Z (trastuzumab Fab), 3AVE (human IgG1 Fc).
'''

import os
import urllib.request

import numpy as np
import pytest

from pymol import cmd
from pymolx import antibody, diagnostics, undo

VH = ('EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRF'
      'TISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSK')
VL = ('DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSG'
      'TDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVV')
FC = ('TCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPRE'
      'EQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTK')


@pytest.fixture(scope='module')
def pdb(tmp_path_factory):
    folder = tmp_path_factory.mktemp('pdb')
    paths = {}
    for code in ('1hzh', '1n8z', '3ave'):
        path = folder / (code + '.cif')
        try:
            urllib.request.urlretrieve(
                'https://files.rcsb.org/download/%s.cif' % code, path)
        except Exception as e:
            pytest.skip('PDB entries not available: %s' % e)
        paths[code] = str(path)
    return paths


@pytest.fixture(autouse=True)
def session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


def split_template(pdb):
    '''
    1HZH cut into Fab1 (H/L), Fab2 (K/M) and Fc (H/K, with glycans), each
    moved away
    '''
    cmd.load(pdb['1hzh'], '1hzh')
    cmd.create('fab', '1hzh and ((chain H and resi 1-225) or chain L) and '
               'polymer')
    cmd.create('fab2', '1hzh and ((chain K and resi 1-225) or chain M) and '
               'polymer')
    cmd.create('fc', '1hzh and ((chain H+K and resi 236-500 and polymer) or '
               'chain A+B)')
    for obj, angle in (('fab', 40), ('fab2', -70), ('fc', 100)):
        cmd.rotate('y', angle, obj, camera=0, origin=[0, 0, 0])
        cmd.translate([15, -20, 30], obj, camera=0)


def sequence(selection):
    return ''.join(cmd.get_fastastr(selection).split('\n')[1:])


def ca_coords(selection):
    return np.array(cmd.get_coords('(%s) and name CA and present' %
                                   selection))


def test_chain_kinds():
    ch1 = 'GGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSS'
    assert antibody.chain_kind(VH + ch1) == 'heavy'
    assert antibody.chain_kind(VL) == 'light'
    assert antibody.chain_kind(FC) == 'fc'
    assert antibody.chain_kind('GSHMKLLV' * 10) == ''
    # single chains: scFv (VH-linker-VL or VL-linker-VH), VHH
    vh, vl = VH[:120], VL[:107]
    for scfv in (vh + 'GGGGS' * 3 + vl, vl + 'GGGGS' * 3 + vh):
        assert antibody.chain_kind(scfv) == 'scfv'
    assert [d[0] for d in antibody.v_domains(vl + 'GGGGS' * 3 + vh)] == [
        'light', 'heavy']
    assert antibody.chain_kind(vh) == 'vhh'
    # (framework 4 searched after CDR2: not the WGDG of VIWGDG)
    vh_wgdg = vh[:49] + 'VIWGDG' + vh[55:]
    assert antibody.v_domains(vh_wgdg)[0][2] == antibody.v_domains(vh)[0][2]
    # the constant domains start after framework 4
    assert VH[antibody.constant_start(VH):].startswith('ASTKGP')
    assert VL[antibody.constant_start(VL):].startswith('TVAAP')
    assert antibody.constant_start(FC) == 0


def test_parse_chains():
    assert antibody.parse_chains('HL', 2) == ['H', 'L']
    assert antibody.parse_chains('AA, B', 2) == ['AA', 'B']
    assert antibody.parse_chains('', 2) == ['', '']
    with pytest.raises(ValueError, match='expected 4'):
        antibody.parse_chains('HLK', 4)


def test_superpose_rejects_outliers():
    rng = np.random.default_rng(1)
    target = rng.normal(size=(60, 3)) * 10
    angle = 0.7
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    mobile = (target - [5, 2, 1]) @ rotation   # target = mobile R^T + t
    mobile[:6] += 8.0                          # outliers
    r, t, rmsd, kept = antibody.superpose(mobile, target)
    assert kept == 54 and rmsd < 1e-6
    assert np.allclose(mobile[6:] @ r.T + t, target[6:], atol=1e-6)


def test_detection(pdb):
    cmd.load(pdb['1hzh'], '1hzh')
    cmd.load(pdb['1n8z'], '1n8z')
    cmd.load(pdb['3ave'], '3ave')
    assert antibody.detect_template('1hzh') == ('H', 'L', 'K', 'M')
    assert antibody.detect_fab('1n8z') == ('B', 'A')   # (C: the antigen)
    assert antibody.detect_fc('3ave') == ('A', 'B')
    with pytest.raises(ValueError, match='two full heavy chains'):
        antibody.detect_template('1n8z')


def test_parts_of_the_template_rebuild_it(pdb):
    split_template(pdb)
    a = antibody.assemble('fab', 'fc', fab2='fab2', template='1hzh',
                          fab1_chains='HL', fab2_chains='KM',
                          fc_chains='HK')
    assert [fit[1] < 0.01 for fit in a.fits] == [True] * 3
    name = antibody.load(a, 'igg')
    assert cmd.get_chains(name) == ['H', 'K', 'L', 'M']
    for chain in 'HLKM':
        assert sequence('%s and chain %s' % (name, chain)) == \
            sequence('1hzh and chain %s and polymer' % chain)
        model = ca_coords('%s and chain %s' % (name, chain))
        template = ca_coords('1hzh and chain %s and polymer' % chain)
        assert model.shape == template.shape
        assert np.abs(model - template).max() < 0.01
    # residues numbered from 1; parts in "custom"
    numbers = []
    cmd.iterate('%s and chain H and name CA' % name, 'numbers.append(resv)',
                space={'numbers': numbers})
    assert numbers == list(range(1, len(numbers) + 1))
    counts = {tag: cmd.count_atoms('%s and name CA and custom %s' % (
        name, tag)) for tag in antibody.PART_TAGS.values()}
    assert counts['FAB1'] == counts['FAB2'] and counts['HNGE'] >= 14
    assert counts['FC'] > 400
    # glycans went with their Fc chains
    assert cmd.count_atoms('%s and not polymer' % name) == \
        cmd.count_atoms('1hzh and chain A+B and not solvent')
    assert cmd.count_atoms('%s and not polymer and chain H' % name) and \
        cmd.count_atoms('%s and not polymer and chain K' % name)


def make_scfv(name='scfv'):
    '''
    1HZH's second Fv as one chain K: VH, then VL (numbered from 301), the
    linker missing
    '''
    cmd.create(name, '1hzh and ((chain K and resi 1-113) or (chain M and '
               'resi 1-107)) and polymer')
    # (renumbered by resi: setting resv drops insertion codes, e.g. 28A)
    cmd.alter('%s and chain M' % name, 'chain="K"; '
              'resi = str(resv + 300) + resi.lstrip("-0123456789")')
    cmd.sort(name)
    cmd.translate([0, 30, -20], name, camera=0)


def test_formats(pdb):
    split_template(pdb)
    make_scfv()
    assert antibody.detect_single('scfv', 'scfv') == 'K'
    a = antibody.assemble('fab', 'fc', 'scfv', '1hzh', arms='fab+scfv',
                          fab1_chains='HL')
    assert a.labels['fab1'] == 'Fab' and a.labels['fab2'] == 'scFv'
    assert a.tags()['fab2'] == 'SCFV2'
    # the scFv's V domains on the template's, before moving to the hinge
    assert dict((f[0], f[1]) for f in a.fits)['scFv'] < 0.01
    assert any('scFv (scfv, chain K, VH-VL)' in line and 'on the '
               'template\'s Fv' in line for line in a.log)
    name = antibody.load(a, 'bispecific')
    # on the template's Fv (where it came from)
    model = ca_coords('%s and custom SCFV2' % name)
    fv = np.vstack([ca_coords('1hzh and chain K and resi 1-113'),
                    ca_coords('1hzh and chain M and resi 1-107')])
    assert np.linalg.norm(model[:, None] - fv[None], axis=2).min(
        axis=1).max() < 0.01
    assert cmd.get_chains(name) == ['H', 'K', 'L']
    # chain K: scFv + hinge + Fc, its VL lighter
    vh = sequence('scfv and chain K and resi 1-113')
    vl = sequence('scfv and chain K and resi 301-500')
    k = sequence('%s and chain K' % name)
    assert k.startswith(vh + vl)
    assert k.endswith(sequence('fc and chain K and polymer'))
    assert sequence('%s and chain K and custom SCFV2' % name) == vh + vl
    light = a.light_residues()
    assert ('K', str(len(vh) + 1)) in light and ('K', '1') not in light
    # next to the hinge instead, with a linker to it
    b = antibody.assemble('fab', 'fc', 'scfv', '1hzh', arms='fab+scfv',
                          fab1_chains='HL', single_position='hinge',
                          arm_linkers=('', '(G4S)2'))
    assert any('toward the hinge' in line for line in b.log)
    k = antibody._sequence(b.chains['K'])
    assert k.startswith(vh + vl + 'GGGGS' * 2)
    linker = [('K', str(i + 1)) for i, r in enumerate(b.chains['K'])
              if r.part == 'linker']
    assert len(linker) == 10 and set(linker) <= b.free
    with pytest.raises(ValueError, match='fv or hinge'):
        antibody.assemble('fab', 'fc', 'scfv', '1hzh', arms='fab+scfv',
                          single_position='middle')
    # one-armed: chain K is the hinge and the Fc
    a = antibody.assemble('fab', 'fc', '', '1hzh', arms='fab+none')
    assert sorted(a.chains) == ['H', 'K', 'L']
    assert a.labels['fab2'] == 'Arm 2 (none)'
    assert {r.part for r in a.chains['K']} == {'hinge', 'fc'}
    for bad in ('fab+igm', 'none+none', 'fab+fab+fab'):
        with pytest.raises(ValueError, match='arm'):
            antibody.parse_arms(bad)
    assert antibody.parse_arms('scFv') == ['scfv', 'scfv']


def test_fusions(pdb):
    assert antibody.parse_linker('(G4S)3') == 'GGGGS' * 3
    assert antibody.parse_linker('G4S x2') == 'GGGGS' * 2
    assert antibody.parse_linker('') == ''
    assert antibody.parse_linker('PAPAP gsg') == 'PAPAPGSG'   # custom
    # spaces separate parts: the repeat is of EAAAK alone
    assert antibody.parse_linker('GSPG EAAAK x2 GS') == \
        'GSPG' + 'EAAAK' * 2 + 'GS'
    assert antibody.parse_linker('GSPG (EAAAK)2 GS') == \
        'GSPG' + 'EAAAK' * 2 + 'GS'
    assert antibody.parse_linker('(EAAAK)3') == 'EAAAK' * 3
    with pytest.raises(ValueError, match='unknown letters'):
        antibody.parse_linker('GGGGB')
    found = antibody.parse_fusions('heavy:scfv/K:(G4S)3 light1:v')
    assert [(f['site'], f['object'], f['chain'], f['linker'])
            for f in found] == [
        ('heavy1', 'scfv', 'K', 'GGGGS' * 3),
        ('heavy2', 'scfv', 'K', 'GGGGS' * 3),
        ('light1', 'v', '', antibody.DEFAULT_LINKER)]
    assert all(f['orientation'] == '' and
               f['scfv_linker'] == antibody.DEFAULT_SCFV_LINKER
               for f in found)
    with pytest.raises(ValueError, match='site'):
        antibody.parse_fusions('tail:scfv')

    split_template(pdb)
    make_scfv()
    plain = antibody.assemble('fab', 'fc', 'fab2', '1hzh')
    a = antibody.assemble('fab', 'fc', 'fab2', '1hzh',
                          fusions='heavy:scfv:(G4S)3')
    assert not a.warnings
    unit = antibody._sequence(antibody.read_chain('scfv', 'K'))
    for chain in 'HK':
        before = antibody._sequence(plain.chains[chain])
        assert antibody._sequence(a.chains[chain]) == \
            before + 'GGGGS' * 3 + unit
    tags = a.tags()
    assert tags['fusion_heavy1'] == 'FUSH1' and tags['linker'] == 'LINK'
    name = antibody.load(a, 'igg_scfv')
    assert cmd.count_atoms('%s and custom FUSH2 and name CA' % name) == \
        cmd.count_atoms('scfv and name CA')
    # the linker: no coordinates yet, free in refinement
    assert cmd.count_atoms('%s and custom LINK and present' % name) == 0
    linker = [(c, str(k)) for c in 'HK' for k, r in
              enumerate(a.chains[c], 1) if r.part == 'linker']
    assert len(linker) == 30 and set(linker) <= a.free
    # the unit apart from the rest
    unit_ca = ca_coords('%s and custom FUSH1' % name)
    rest = ca_coords('%s and not custom FUSH1+FUSH2+LINK' % name)
    d = np.linalg.norm(unit_ca[:, None] - rest[None], axis=2).min()
    assert d > 3.5
    # light chains need a Fab arm
    with pytest.raises(ValueError, match='no light chain'):
        antibody.assemble('fab', 'fc', 'scfv', '1hzh', arms='fab+scfv',
                          fusions='light2:scfv')


def _cdr_angle(assembly, part, chain):
    """
    Degrees between a fused unit's CDR direction and straight out (away
    from the rest of the model)
    """
    unit = [r for r in assembly.chains[chain] if r.part == part]
    body = np.array([r.ca for r in unit if r.ca is not None])
    rest = np.array([r.ca for residues in assembly.chains.values()
                     for r in residues
                     if r.ca is not None and r.part != part])
    return antibody._angle_between(
        antibody.cdr_direction(unit),
        body.mean(axis=0) - rest.mean(axis=0))


def test_cdrs_face_outwards(pdb):
    # CDR3 from the sequence: after the last Cys, up to framework 4
    for code, chain, loops in (('1n8z', 'B', ['WGGDGFYAMDY']),
                               ('2ghw', 'B', ['DRSYYLDY', 'RSNWPPT'])):
        cmd.load(pdb[code], code) if code in pdb else cmd.fetch(code)
        residues = antibody.read_chain(code, chain)
        sequence = antibody._sequence(residues)
        found = [sequence[slice(*antibody.cdr3_range(sequence, s, e))]
                 for _, s, e in antibody.v_domains(sequence)]
        assert found == loops
    # a Fab's own CDRs point away from its constant domains
    heavy = antibody.read_chain('1n8z', 'B')
    light = antibody.read_chain('1n8z', 'A')
    ends = [antibody.constant_start(antibody._sequence(c))
            for c in (heavy, light)]
    fv = heavy[:ends[0]] + light[:ends[1]]
    constant = [r for c, end in zip((heavy, light), ends) for r in c[end:]]
    out = np.mean([r.ca for r in fv if r.ca is not None], axis=0) - \
        np.mean([r.ca for r in constant if r.ca is not None], axis=0)
    assert antibody._angle_between(antibody.cdr_direction(fv), out) < 30

    split_template(pdb)
    make_scfv()
    # the fused unit is turned so its binding site faces the solvent
    sites = (('fusion_heavy1', 'H'), ('fusion_heavy2', 'K'))
    a = antibody.assemble('fab', 'fc', 'fab2', '1hzh',
                          fusions='heavy:scfv:(G4S)3')
    assert _cdr_angle(a, *sites[0]) < antibody.CDR_OUTWARD
    assert sum('CDRs face' in line for line in a.log) == 2
    # "any": the pose is chosen by its contacts alone, and the binding
    # sites end up anywhere (here against the antibody)
    b = antibody.assemble('fab', 'fc', 'fab2', '1hzh', fusions=[
        {'site': site[len('fusion_'):], 'object': 'scfv', 'chain': 'K',
         'linker': 'GGGGS' * 3, 'cdrs': 'any'} for site, _ in sites])
    assert not any('CDRs face' in line for line in b.log)
    for site in sites:
        assert _cdr_angle(a, *site) < _cdr_angle(b, *site) - 25
    assert min(_cdr_angle(b, *site) for site in sites) > 85
    # a VHH is turned the same way (one domain, its CDRs at one end)
    cmd.create('nano', 'scfv and chain K and resi 1-113')
    c = antibody.assemble('fab', 'fc', 'fab2', '1hzh',
                          fusions='heavy1:nano:(G4S)2')
    assert _cdr_angle(c, 'fusion_heavy1', 'H') < antibody.CDR_OUTWARD


def _span(assembly, chain):
    """
    The linker's end-to-end distance, A per residue, in a chain of an
    assembly (its residues have the placed coordinates)
    """
    linker = [k for k, r in enumerate(assembly.chains[chain])
              if r.part == 'linker']
    residues = assembly.chains[chain]
    before = max(k for k in range(linker[0]) if residues[k].ca is not None)
    after = min(k for k in range(linker[-1] + 1, len(residues))
                if residues[k].ca is not None)
    distance = np.linalg.norm(residues[after].ca - residues[before].ca)
    return distance / (after - before)


def test_extended_linker(pdb):
    split_template(pdb)
    make_scfv()
    relaxed = antibody.assemble('fab', 'fc', 'fab2', '1hzh',
                                fusions='heavy1:scfv:(G4S)3')
    extended = antibody.assemble('fab', 'fc', 'fab2', '1hzh',
                                 fusions='heavy1:scfv:(G4S)3',
                                 fusion_linkers='extended')
    # the unit is placed at the stretched linker's length, so the gap the
    # refinement has to span is nearly one residue's step per residue
    assert _span(relaxed, 'H') < 1.0
    assert antibody.EXTENDED_STEP - 0.3 < _span(extended, 'H') <= \
        antibody.CA_STEP
    assert any('extended' in line for line in extended.log)
    assert not any('extended' in line for line in relaxed.log)
    # ... and that is not a gap the builder is warned about
    assert not extended.warnings
    # per fusion, not only for all of them
    mixed = antibody.assemble('fab', 'fc', 'fab2', '1hzh', fusions=[
        {'site': 'heavy1', 'object': 'scfv', 'chain': 'K',
         'linker': 'GGGGS' * 3, 'linker_shape': 'extended'},
        {'site': 'heavy2', 'object': 'scfv', 'chain': 'K',
         'linker': 'GGGGS' * 3}])
    assert _span(mixed, 'H') > 2.5 and _span(mixed, 'K') < 1.0
    with pytest.raises(ValueError, match='relaxed or extended'):
        antibody.assemble('fab', 'fc', 'fab2', '1hzh',
                          fusions='heavy1:scfv', fusion_linkers='straight')


def test_collapsed_linker_setup(pdb):
    split_template(pdb)
    make_scfv()
    a = antibody.assemble('fab', 'fc', 'fab2', '1hzh',
                          fusions='heavy1:scfv:(G4S)3',
                          fusion_linkers='collapsed')
    # placed extended, and the dynamics told what to draw in
    assert _span(a, 'H') > 2.5
    assert any('collapsed by the dynamics' in line for line in a.log)
    assert len(a.collapse) == 1
    one = a.collapse[0]
    assert one['chain'] == 'H' and one['residues'] == 15
    chain = a.chains['H']
    assert chain[one['anchor'] - 1].part == 'fc'
    assert chain[one['linker'] - 1].part == 'linker'
    assert chain[one['unit'] - 1].part == 'fusion_heavy1'
    # the job: the unit and its linker move, pulled back to a coil's span
    pytest.importorskip('openmm')
    pytest.importorskip('pdbfixer')
    name = antibody.load(a, 'built')
    job = antibody.refine_job(a, name, 'refined', minimize=False)
    settings = job.dynamics
    assert ('H', str(one['unit'])) in settings['mobile']
    assert ('H', str(len(chain))) in settings['mobile']     # the unit's end
    assert settings['mobile'] <= job.free_residues    # none are restrained
    (anchor_atom, unit_atom, distance, _), = settings['pull']
    assert anchor_atom == ('H', str(one['anchor']), 'CA')
    assert unit_atom == (('H', str(one['unit']), 'CA'),)
    # pulled back to the span of a coil of the whole flexible junction
    flexible = one['unit'] - one['anchor'] - 1
    assert flexible >= 15                     # the linker, and the tail
    assert distance == pytest.approx(antibody.COIL_STEP * flexible ** 0.5)
    # the domain itself is held in shape while it is pulled
    assert ('H', str(one['unit']), 'CA') in settings['rigid'][0]
    assert ('H', str(len(chain)), 'O') in settings['rigid'][0]
    # the other shapes ask for no dynamics
    for shape in ('relaxed', 'extended'):
        b = antibody.assemble('fab', 'fc', 'fab2', '1hzh',
                              fusions='heavy1:scfv:(G4S)3',
                              fusion_linkers=shape)
        assert not b.collapse
        assert antibody.refine_job(b, name, 'x', minimize=False).dynamics \
            is None


def test_collapsed_linker_runs(pdb):
    pytest.importorskip('openmm')
    pytest.importorskip('pdbfixer')
    split_template(pdb)
    make_scfv()
    # a single domain: no linker of its own, so nothing but the pull can
    # strain it (this scFv's two domains are fused directly, which does)
    cmd.create('nano', 'scfv and chain K and resi 1-113')
    a = cmd.build_antibody('fab', 'fc', 'fab2', '1hzh',
                           fusions='heavy:nano:(G4S)3',
                           fusion_linkers='collapsed', name='igg', quiet=1)
    assert not a.warnings
    assert not antibody.long_peptide_bonds('igg')
    for chain in 'HK':
        ca = cmd.get_coords('igg and chain %s and custom LINK and name CA'
                            % chain)
        span = np.linalg.norm(ca[-1] - ca[0])
        # started extended (about 3.1 A per residue), now partly coiled
        assert span < 2.5 * (len(ca) - 1)
        assert span > 3.0            # not a knot either
    # the fused domain keeps its fold: a network holds it in shape while
    # it is pulled. Without it the pull peeled its sheets open (5.8 A and
    # nine strand residues lost on a real scFv), so the same build
    # without the dynamics is the measure
    cmd.build_antibody('fab', 'fc', 'fab2', '1hzh',
                       fusions='heavy:nano:(G4S)3',
                       fusion_linkers='extended', name='still', quiet=1)
    strands = {}
    for name, model in (('collapsed', 'igg'), ('still', 'still')):
        cmd.create(name + '_unit', '%s and custom FUSH1 and polymer' % model)
        cmd.dss(name + '_unit')
        found = {}
        cmd.iterate('%s_unit and name CA' % name,
                    'counts[ss] = counts.get(ss, 0) + 1',
                    space={'counts': found})
        strands[name] = found.get('S', 0)
    assert cmd.align('collapsed_unit', 'still_unit', cycles=0)[0] < 3.0
    assert strands['collapsed'] >= strands['still'] - 4
    # the dynamics didn't spoil the geometry (prep minimizes again after)
    scores = {i.detail.split()[0]: i.detail
              for i in diagnostics.run('igg') if i.category == 'Scores'}
    assert float(scores['clashscore'].split()[1]) < 2.0
    assert 'C-beta deviations: ' in scores['C-beta']
    assert int(scores['C-beta'].rsplit(' ', 1)[1]) < 20


def test_scfv_orientation(pdb):
    assert antibody.parse_orientation('VL-VH') == 'vl-vh'
    assert antibody.parse_orientation('vhvl') == 'vh-vl'
    assert antibody.parse_orientation('') == ''
    with pytest.raises(ValueError, match='vh-vl or vl-vh'):
        antibody.parse_orientation('vh-vh')
    spec = antibody.parse_fusions('heavy1:scfv/K:(G4S)2:vl-vh:(G4S)4')[0]
    assert spec['orientation'] == 'vl-vh'
    assert spec['scfv_linker'] == 'GGGGS' * 4

    split_template(pdb)
    make_scfv()          # chain K: VH (1-113) then VL (301-407)
    vh = sequence('scfv and chain K and resi 1-113')
    vl = sequence('scfv and chain K and resi 301-500')
    joined = 'GGGGS' * 3
    for site, chain in (('heavy1', 'H'), ('heavy2', 'K')):
        a = antibody.assemble('fab', 'fc', 'fab2', '1hzh',
                              fusions='%s:scfv:(G4S)2:vl-vh' % site)
        unit = antibody._sequence(
            [r for r in a.chains[chain]
             if r.part.startswith('fusion') or r.part == 'linker'])
        # outer linker, then VL, the new linker, VH (the structure's own
        # linker between the domains is left out)
        assert unit == 'GGGGS' * 2 + vl + joined + vh
        assert any('reordered to VL-VH' in line for line in a.log)
    # the same Fv: the domains keep their relative positions
    plain = antibody.assemble('fab', 'fc', 'fab2', '1hzh',
                              fusions='heavy1:scfv:(G4S)2')
    name = antibody.load(plain, 'as_is')
    other = antibody.load(a, 'reordered')
    cmd.create('unit_a', '%s and custom FUSH2 and present' % other)
    cmd.create('unit_b', '%s and custom FUSH1 and present' % name)
    assert cmd.align('unit_a', 'unit_b', cycles=0)[0] < 0.01
    # VH-VL asked for, the structure already VH-VL: unchanged
    b = antibody.assemble('fab', 'fc', 'fab2', '1hzh',
                          fusions='heavy1:scfv:(G4S)2:vh-vl')
    assert not any('reordered' in line for line in b.log)
    assert antibody._sequence(b.chains['H']) == \
        antibody._sequence(plain.chains['H'])
    # an scFv arm reordered the same way
    c = antibody.assemble('fab', 'fc', 'scfv', '1hzh', arms='fab+scfv',
                          fab1_chains='HL', arm_orientations=('', 'vl-vh'),
                          arm_scfv_linkers=('', '(G4S)3'))
    assert antibody._sequence(c.chains['K']).startswith(vl + joined + vh)
    assert any('reordered to VL-VH' in line for line in c.log)
    # a single domain (VHH): nothing to reorder, the order ignored
    cmd.create('vhh', 'scfv and chain K and resi 1-113')
    d = antibody.assemble('fab', 'fc', 'fab2', '1hzh',
                          fusions='heavy1:vhh:(G4S)2:vl-vh')
    assert not any('reordered' in line for line in d.log)
    assert antibody._sequence(
        [r for r in d.chains['H'] if r.part.startswith('fusion')]) == vh


def test_hinge_sequence(pdb):
    split_template(pdb)
    igg4 = dict(antibody.HINGES)['Human IgG4 S228P']
    a = antibody.assemble('fab', 'fc', 'fab2', '1hzh', hinge_sequence=igg4)
    name = antibody.load(a, 'igg')
    for chain in 'HK':
        hinge = sequence('%s and chain %s and custom HNGE' % (name, chain))
        assert hinge and hinge in igg4
        full = sequence('%s and chain %s' % (name, chain))
        fab = sequence('%s and chain %s and custom FAB*' % (name, chain))
        fc = sequence('%s and chain %s and custom FC' % (name, chain))
        assert full == fab + hinge + fc
    with pytest.raises(ValueError, match='unknown letters'):
        antibody.assemble('fab', 'fc', 'fab2', '1hzh',
                          hinge_sequence='EPKSX')


def test_hinge_structure(pdb):
    split_template(pdb)
    cmd.create('hinge', '1hzh and chain H+K and resi 222-245 and polymer')
    cmd.translate([40, 0, 0], 'hinge', camera=0)
    a = antibody.assemble('fab', 'fc', 'fab2', '1hzh', hinge='hinge',
                          hinge_chains='HK')
    assert any('hinge: placed' in line for line in a.log)
    name = antibody.load(a, 'igg')
    # the hinge back on the template's, the chains as the template's
    for chain in 'HK':
        assert sequence('%s and chain %s' % (name, chain)) == \
            sequence('1hzh and chain %s and polymer' % chain)
    model = ca_coords('%s and chain H and custom HNGE' % name)
    assert len(model)
    d = np.linalg.norm(model[:, None] - ca_coords('1hzh and chain H'),
                       axis=2).min(axis=1)
    assert d.max() < 0.01


def test_errors(pdb):
    split_template(pdb)
    with pytest.raises(ValueError, match='Fc'):
        antibody.assemble('fab', '', template='1hzh')
    cmd.create('full_heavy', '1hzh and chain H+L and polymer')
    with pytest.raises(ValueError, match='reaches into the Fc'):
        antibody.assemble('full_heavy', 'fc', template='1hzh',
                          fab1_chains='HL')
    with pytest.raises(ValueError, match='different IDs'):
        antibody.assemble('fab', 'fc', template='1hzh', chains='HLHM')
    with pytest.raises(ValueError, match='not an object'):
        antibody.assemble('no_such_thing', 'fc', template='1hzh')


def test_command_one_undo_step(pdb):
    split_template(pdb)
    cmd.undo_enable()
    try:
        steps = len(undo.stack.undo_stack)
        cmd.build_antibody('fab', 'fc', 'fab2', '1hzh', name='model',
                           refine=0, quiet=1)
        assert len(undo.stack.undo_stack) == steps + 1
    finally:
        cmd.undo_disable()
    # a selection keyword as name: renamed like PyMOL does
    assert 'model_' in cmd.get_names()


def test_cysteine_pairs(pdb):
    for code in ('1hzh', '1n8z', '3ave'):
        cmd.load(pdb[code], code)
    a = antibody.assemble('1n8z', '3ave', template='1hzh')
    pairs = {antibody.pair_label((p['a'], p['b'])): p
             for p in antibody.cysteine_pairs(a)}
    # domain disulfides and the first hinge one form by themselves
    assert pairs['L/23-L/88']['bonded'] and pairs['H/229-K/229']['bonded']
    assert sum(p['bonded'] for p in pairs.values()) == 13
    # heavy-light chain disulfides: close, not bonded
    for label in ('H/223-L/214', 'K/223-M/214'):
        assert not pairs[label]['bonded'] and pairs[label]['distance'] < 5
        assert set(pairs[label]['parts']) in ({'hinge', 'fab1'},
                                              {'hinge', 'fab2'})
    # the second hinge disulfide is offered though far (18.6 A in 1HZH);
    # hinge cysteines of different places aren't
    assert pairs['H/232-K/232']['distance'] > 15
    assert 'H/223-K/229' not in pairs and 'H/229-K/232' not in pairs
    # without refinement: just the bonds
    name = antibody.load(a, 'plain')
    chosen = [(('H', '223'), ('L', '214')), (('H', '229'), ('K', '229'))]
    antibody.set_disulfides(name, chosen)
    assert len(cmd.get_model('%s and name SG' % name).bond) == 2


def test_refinement(pdb):
    pytest.importorskip('openmm')
    pytest.importorskip('pdbfixer')
    cmd.load(pdb['1hzh'], '1hzh')
    cmd.load(pdb['1n8z'], '1n8z')
    cmd.load(pdb['3ave'], '3ave')
    pairs = antibody.cysteine_pairs(antibody.assemble('1n8z', '3ave',
                                                      template='1hzh'))
    chosen = ', '.join(antibody.pair_label((p['a'], p['b'])) for p in pairs
                       if p['distance'] is not None and p['distance'] < 5)
    a = cmd.build_antibody('1n8z', '3ave', template='1hzh',
                           disulfides=chosen, quiet=1)
    assert a.name == 'igg_model'
    assert 'igg_model' in cmd.get_names()
    assert not [n for n in cmd.get_names() if n.startswith('_ab')]
    assert not a.warnings
    assert not antibody.long_peptide_bonds('igg_model')
    # no residues left without coordinates, but at the chain ends
    assert cmd.count_atoms('igg_model and not present') == 0
    # the glycosylated Asn keeps one H on ND2 (none where C1 of NAG is)
    assert cmd.count_atoms('igg_model and hydro and neighbor (resn ASN and '
                           'name ND2 and neighbor resn NAG)') == 2
    # disulfides: the 13 close ones and the chosen heavy-light ones, made
    log = '\n'.join(a.log)
    assert 'made (SG-SG before minimizing): CYS H/223-CYS L/214' in log
    assert 'disulfide bonds: 15' in log
    for c1, r1, c2, r2 in (('H', '223', 'L', '214'),
                           ('K', '223', 'M', '214'),
                           ('H', '229', 'K', '229')):
        assert cmd.get_distance(
            'igg_model and chain %s and resi %s and name SG' % (c1, r1),
            'igg_model and chain %s and resi %s and name SG' % (c2, r2)) \
            < 2.3
    assert len(cmd.get_model('igg_model and name SG').bond) == 15


def test_dialog(pdb):
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    QtCore = pytest.importorskip('pymol.Qt').QtCore
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import pymol.plugins  # noqa: F401 (adds the startup plugin path)
    from pmg_tk.startup.antibody_gui.main import AntibodyDialog, PartDialog

    split_template(pdb)
    dialog = AntibodyDialog(None, _self=cmd)
    try:
        dialog.resize(900, 600)
        sch = dialog.schematic
        sch.resize(440, 470)
        t = sch._transform()
        # each box and the hinge map back to their part
        for part, _, path in sch.boxes:
            center = t.map(path.boundingRect().center()).toPoint()
            assert sch.part_at(center) == part
        assert sch.part_at(t.map(sch.hinge_area.center()).toPoint()) == \
            'hinge'
        assert sch.part_at(t.map(sch.label_at['fc']).toPoint()) == ''
        # a part's prompt recognizes its chains
        prompt = PartDialog('fab1', None, dialog, _self=cmd)
        prompt.object.setCurrentText('fab')
        prompt._accept()
        assert prompt.spec() == {'object': 'fab', 'chains': ('H', 'L')}
        prompt.deleteLater()
        prompt = PartDialog('hinge', None, dialog, _self=cmd)
        prompt._use_preset(prompt.preset.findText('Human IgG4'))
        spec = prompt.spec()
        assert spec['mode'] == 'sequence' and spec['sequence'] == \
            'ESKYGPPCPSCP'
        prompt.deleteLater()
        dialog.specs.update({
            'fab1': {'object': 'fab', 'chains': ('H', 'L')},
            'fab2': {'object': 'fab2', 'chains': ('K', 'M')},
            'fc': {'object': 'fc', 'chains': ('H', 'K')},
            'template': {'object': '1hzh', 'chains': ('H', 'L', 'K', 'M')}})
        dialog._show_specs()
        assert sch.assigned['fab1'] and not sch.assigned['hinge']
        assert dialog.part_labels['fab2'].text() == 'fab2 (heavy K, light M)'
        from pmg_tk.startup.antibody_gui.main import DisulfideDialog
        chooser = DisulfideDialog(antibody.assemble(**dialog._arguments()),
                                  None, dialog)
        rows = {chooser.table.item(k, 0).text(): k
                for k in range(chooser.table.rowCount())}
        assert chooser.table.item(rows['H/22-H/96'], 0).checkState() == \
            QtCore.Qt.CheckState.Checked
        chooser.table.item(rows['H/22-H/96'], 0).setCheckState(
            QtCore.Qt.CheckState.Unchecked)
        chooser.other.setText('H/1-K/1')
        chooser._accept()
        assert (('H', '22'), ('H', '96')) not in chooser.chosen
        assert chooser.chosen[-1] == (('H', '1'), ('K', '1'))
        dialog.disulfides = chooser.chosen
        dialog.disulfide_args = dialog._arguments()
        dialog._show_disulfides()
        assert dialog.disulfide_label.text() == '%d chosen' % len(
            chooser.chosen)
        chooser.deleteLater()
        # arm formats: the picture, the parts, the arguments
        boxes = lambda part: len([b for b in sch.boxes if b[0] == part])
        assert boxes('fab2') == 4
        make_scfv()
        dialog.format_boxes['fab2'].setCurrentIndex(
            dialog.format_boxes['fab2'].findData('scfv'))
        assert dialog.specs['fab2'] is None and boxes('fab2') == 2
        assert dialog.part_titles['fab2'].text() == 'scFv'
        assert dialog.part_titles['fab1'].text() == 'Fab'
        with pytest.raises(ValueError, match='give the scFv structure'):
            dialog._arguments()
        prompt = PartDialog('fab2', None, dialog, _self=cmd, fmt='scfv',
                            title='scFv')
        prompt.object.setCurrentText('scfv')
        prompt._use_linker(prompt.linker_preset.findData('GGGGS' * 3))
        assert prompt.spec() == {
            'object': 'scfv', 'chains': ('K',), 'linker': 'GGGGS' * 3,
            'orientation': '',
            'scfv_linker': antibody.DEFAULT_SCFV_LINKER}
        prompt.deleteLater()
        dialog.specs['fab2'] = {'object': 'scfv', 'chains': ('K',),
                                'linker': 'GGGGS' * 3}
        dialog._show_specs()
        assert dialog.position.isEnabled()
        args = dialog._arguments()
        assert args['arms'] == 'fab+scfv' and args['fab2_chains'] == 'K'
        assert args['single_position'] == 'fv'
        assert args['arm_linkers'] == ('', 'GGGGS' * 3)
        # on the Fv: drawn where a Fab's Fv is, with a linker line
        assert [p for p, _ in sch.linker_lines] == ['fab2']
        dialog.format_boxes['fab2'].setCurrentIndex(
            dialog.format_boxes['fab2'].findData('none'))
        assert boxes('fab2') == 0
        assert not dialog.part_buttons['fab2'].isEnabled()
        dialog.format_boxes['fab1'].setCurrentIndex(
            dialog.format_boxes['fab1'].findData('none'))
        assert dialog.formats['fab1'] == 'fab'     # (one arm at least)
        dialog.format_boxes['fab2'].setCurrentIndex(
            dialog.format_boxes['fab2'].findData('fab'))
        dialog.specs['fab2'] = {'object': 'fab2', 'chains': ('K', 'M')}
        dialog._show_specs()
        # fusions: dashed sites under the Fc, the prompt, the arguments
        sites = {part: t.map(center).toPoint()
                 for part, _, center in sch.sites}
        assert sorted(sites) == ['fusion_heavy1', 'fusion_heavy2']
        assert sch.part_at(sites['fusion_heavy2']) == 'fusion_heavy2'
        prompt = PartDialog('fusion_heavy1', None, dialog, _self=cmd)
        prompt.object.setCurrentText('scfv')
        prompt._use_linker(prompt.linker_preset.findData('GGGGS' * 3))
        spec = prompt.spec()
        assert spec == {'object': 'scfv', 'chains': ('K',),
                        'linker': 'GGGGS' * 3, 'both': True,
                        'kind': 'scfv', 'orientation': '',
                        'scfv_linker': antibody.DEFAULT_SCFV_LINKER,
                        'cdrs': 'out', 'linker_shape': 'relaxed'}
        # any sequence typed in (shorthand expanded), wrong letters shown
        prompt.linker.setText('GSPG EAAAK x2 GS')
        assert prompt.linker_note.text() == '16 residues: GSPGEAAAKEAAAKGS'
        assert prompt.spec()['linker'] == 'GSPGEAAAKEAAAKGS'
        prompt.linker.setText('GGGGB')
        assert 'unknown letters B' in prompt.linker_note.text()
        prompt.linker.setText('')
        assert prompt.linker_note.text() == 'no linker (joined directly)'
        # the domain order: an scFv's, with its own linker field
        assert prompt.orientation.isEnabled()
        assert not prompt.scfv_linker.isEnabled()   # (as in the structure)
        prompt.orientation.setCurrentIndex(
            prompt.orientation.findData('vl-vh'))
        assert prompt.scfv_linker.isEnabled()
        prompt.scfv_linker.setText('(G4S)4')
        assert prompt.spec()['orientation'] == 'vl-vh'
        assert prompt.spec()['scfv_linker'] == 'GGGGS' * 4
        # the linker's shape and where the binding site faces
        assert prompt.spec()['linker_shape'] == 'relaxed'
        prompt.shape.setCurrentIndex(prompt.shape.findData('extended'))
        assert prompt.spec()['linker_shape'] == 'extended'
        prompt.shape.setCurrentIndex(prompt.shape.findData('collapsed'))
        assert prompt.spec()['linker_shape'] == 'collapsed'
        prompt.cdrs.setChecked(False)
        assert prompt.spec()['cdrs'] == 'any'
        prompt.deleteLater()
        dialog.specs['fusion_heavy1'] = spec
        dialog._show_specs()
        assert boxes('fusion_heavy1') == 2 and len(sch.sites) == 1
        assert dialog.part_labels['fusions'].text() == \
            'scFv on heavy chain 1'
        assert dialog._arguments()['fusions'] == [
            {'site': 'heavy1', 'object': 'scfv', 'chain': 'K',
             'linker': 'GGGGS' * 3, 'orientation': '',
             'scfv_linker': antibody.DEFAULT_SCFV_LINKER, 'cdrs': 'out',
             'linker_shape': 'relaxed'}]
        dialog.specs['fusion_heavy1'] = None
        dialog._show_specs()
        dialog.refine.setChecked(False)
        dialog.name.setText('built')
        dialog.build()
        assert dialog.status.text().startswith('Done: "built"')
        assert sequence('built and chain H') == sequence('1hzh and chain H '
                                                         'and polymer')
        bonds = cmd.get_model('built and name SG').bond
        assert len(bonds) == len(dialog.disulfides) - 1   # (H/1, K/1: no SG)
        assert not dialog.grab().isNull()
    finally:
        dialog.deleteLater()
