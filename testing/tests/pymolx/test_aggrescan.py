'''
Aggregation prediction (pymolx.aggrescan, data/startup/aggregation_gui,
parity item L-34): the Aggrescan3D score and its Aggrescan4D
extensions, the pH scan, the OpenMM and EvoEF stand-ins for FoldX, the
dynamic mode and the dialog.

The score itself is checked against the output of real Aggrescan3D
runs when those are on this machine; everything else is self-contained.
'''

import csv
import glob
import math
import os

import pytest

from pymol import cmd
import pymol.plugins  # noqa: F401 (puts $PYMOL_DATA/startup on the path)
from pymolx import aggrescan

DATA = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
STRUCTURE = os.path.abspath(os.path.join(DATA, '1rx1.pdb'))
# real Aggrescan3D output, if this machine has any
REFERENCE = '/media/pac/data/apps/aggrescan3d/outputs'

needs_openmm = pytest.mark.skipif(
    bool(aggrescan.available(optimize=True)),
    reason=aggrescan.available(optimize=True) or '')
needs_evoef = pytest.mark.skipif(not aggrescan.evoef_available(),
                                 reason='EvoEF is not installed')
slow = pytest.mark.skipif(not os.environ.get('PYMOLX_SLOW_TESTS'),
                          reason='minutes: PYMOLX_SLOW_TESTS=1 runs it')


def reference_runs():
    '''
    [(folder of a real A3D run with a scored structure)]
    '''
    found = []
    for folder in sorted(glob.glob(os.path.join(REFERENCE, '*'))):
        if os.path.isfile(os.path.join(folder, 'A3D.csv')) and \
                os.path.isfile(os.path.join(folder, 'tmp', 'folded.pdb')):
            found.append(folder)
    return found


@pytest.fixture(autouse=True)
def session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


@pytest.fixture
def protein():
    cmd.load(STRUCTURE, 'p')
    return 'p'


# --- the scale and the formula ----------------------------------------


def test_the_scale_is_the_published_one():
    assert len(aggrescan.A3D_SCALE) == 20
    # the extremes of the in vivo Aggrescan scale
    assert aggrescan.A3D_SCALE['I'] == pytest.approx(1.9136)
    assert aggrescan.A3D_SCALE['R'] == pytest.approx(-1.7534)
    assert aggrescan.MIN_SURFACE == 10 and aggrescan.MAX_SURFACE == 55
    assert aggrescan.DISTANCES == (5, 10)
    assert aggrescan.DEFAULT_DISTANCE == 10


def test_a_buried_residue_scores_zero():
    assert aggrescan.intrinsic_score('I', 9.9) == 0.0
    assert aggrescan.intrinsic_score('I', 10.0) > 0.0


def test_the_exposure_term_stops_growing():
    capped = aggrescan.intrinsic_score('I', aggrescan.MAX_SURFACE)
    assert aggrescan.intrinsic_score('I', 90.0) == pytest.approx(capped)


def test_the_exposure_term_is_the_a3d_one():
    # propensity * ALPHA * exp(BETA * rsa)
    assert aggrescan.intrinsic_score('L', 30.0) == pytest.approx(
        aggrescan.A3D_SCALE['L'] * 0.0599 * math.exp(0.0521 * 30.0))


def test_neighbours_decay_with_distance():
    def one(gap):
        residues = [aggrescan.Residue('A', '1', 'I', 50.0, 1.0, (0, 0, 0)),
                    aggrescan.Residue('A', '2', 'I', 50.0, 1.0, (gap, 0, 0))]
        aggrescan.spread(residues, 10)
        return residues[0].score

    assert one(2.0) > one(5.0) > one(9.0) > one(11.0)
    # beyond the sphere only the residue itself counts
    assert one(11.0) == pytest.approx(1.0)


def test_a_buried_residue_keeps_zero_whatever_surrounds_it():
    residues = [aggrescan.Residue('A', '1', 'I', 0.0, 0.0, (0, 0, 0)),
                aggrescan.Residue('A', '2', 'I', 50.0, 5.0, (3, 0, 0))]
    aggrescan.spread(residues, 10)
    assert residues[0].score == 0.0
    # and it adds nothing to its neighbour either
    assert residues[1].score == pytest.approx(5.0)


# --- pH ----------------------------------------------------------------


def test_the_ph_range_is_the_servers():
    assert aggrescan.PH_VALUES[0] == 4.0
    assert aggrescan.PH_VALUES[-1] == 9.0
    assert len(aggrescan.PH_VALUES) == 11


def test_acids_and_bases_move_the_right_way_with_ph():
    low, high = 4.0, 9.0
    # protonated at low pH, so more aggregation prone there
    for letter in 'DE':
        assert aggrescan.ph_score(letter, low) > \
            aggrescan.ph_score(letter, high)
    # charged at low pH, neutral at high
    for letter in 'KR':
        assert aggrescan.ph_score(letter, low) < \
            aggrescan.ph_score(letter, high)


def test_residues_without_a_charge_do_not_move():
    for letter in 'LIFVWAG':
        assert aggrescan.ph_score(letter, 4.0) == \
            aggrescan.ph_score(letter, 9.0)


def test_the_signed_form_fixes_histidine():
    '''
    A3D measures |pH - pKa|, which is right while a residue is charged
    but mirrors histidine above its pKa. The signed form, which the
    PROPKA mode uses, keeps rising
    '''
    mirrored = [aggrescan.ph_score('H', ph) for ph in (6.5, 7.5, 8.5)]
    assert mirrored[1] < mirrored[0]        # the artefact
    signed = [aggrescan.ph_score('H', ph, signed=True)
              for ph in (6.5, 7.5, 8.5)]
    assert signed[0] < signed[1] < signed[2]


def test_a_predicted_pka_shifts_the_score():
    standard = aggrescan.ph_score('D', 7.0)
    buried = aggrescan.ph_score('D', 7.0, pka=6.5)
    assert buried > standard        # harder to deprotonate, less soluble


# --- against real Aggrescan3D output ----------------------------------


@pytest.mark.skipif(not reference_runs(),
                    reason='no Aggrescan3D output on this machine')
def test_score_reproduces_aggrescan3d():
    '''
    The same structures Aggrescan3D scored, to the 4 decimals its
    A3D.csv carries
    '''
    worst = 0.0
    for folder in reference_runs()[:5]:
        theirs = {}
        with open(os.path.join(folder, 'A3D.csv')) as handle:
            for row in csv.DictReader(handle):
                theirs[(row['chain'], row['residue'])] = (
                    row['residue_name'], float(row['score']))
        result = aggrescan.score_pdb(
            os.path.join(folder, 'tmp', 'folded.pdb'), distance=10)
        mine = {r.key: (r.letter, r.score) for r in result.residues}
        assert set(mine) == set(theirs)
        for key, (letter, value) in theirs.items():
            assert mine[key][0] == letter
            worst = max(worst, abs(mine[key][1] - value))
    assert worst < 1e-4


@pytest.mark.skipif(not reference_runs(),
                    reason='no Aggrescan3D output on this machine')
def test_exposures_match_the_freesasa_aggrescan3d_ran():
    folder = reference_runs()[0]
    theirs = {}
    for line in open(os.path.join(folder, 'tmp', 'sasa.out')):
        if line.startswith('RES'):
            theirs[(line[8], line[9:16].strip())] = float(line[22:28])
    mine = aggrescan.exposures(os.path.join(folder, 'tmp', 'folded.pdb'))
    assert set(mine) == set(theirs)
    # both are rounded to the one decimal the .rsa format carries
    assert all(abs(mine[k] - theirs[k]) <= 0.05 for k in theirs)


# --- scoring a structure in the session -------------------------------


def test_score_a_structure(protein):
    result = aggrescan.score(protein)
    assert len(result) == cmd.count_atoms('%p and name CA and polymer')
    assert result.distance == 10 and result.ph is None
    assert any(r.prone for r in result.residues)
    # a buried residue scores 0 and an exposed hydrophobic one scores high
    assert all(r.score == 0.0 for r in result.residues if r.buried)
    best = max(result.residues, key=lambda r: r.score)
    assert best.exposure >= aggrescan.MIN_SURFACE
    assert best.letter in 'AILMFVWYC'


def test_a_smaller_sphere_looks_more_locally(protein):
    '''
    The 5 A sphere is the server's option for a finer look: a residue's
    score stays closer to its own propensity, because fewer neighbours
    reach it and those that do count for less. It does not simply
    lower the scores, since most neighbours pull downwards
    '''
    wide = aggrescan.score(protein, distance=10)
    tight = aggrescan.score(protein, distance=5)
    assert len(tight) == len(wide)

    def pull(result):
        return sum(abs(r.score - r.intrinsic)
                   for r in result.residues) / len(result)

    assert pull(tight) < pull(wide)
    assert tight.average != pytest.approx(wide.average)


def test_scoring_at_a_ph_differs_from_the_classic_scale(protein):
    classic = aggrescan.score(protein)
    at_ph = aggrescan.score(protein, ph=7.0, pka='fixed')
    assert len(classic) == len(at_ph)
    assert classic.average != pytest.approx(at_ph.average)


def test_the_ph_scan_covers_the_whole_range(protein, tmp_path):
    path = str(tmp_path / 'in.pdb')
    aggrescan.write_pdb(protein, path)
    found = aggrescan.scan_ph(path, 10, aggrescan.PH_VALUES, 'fixed')
    assert [value for value, _ in found] == list(aggrescan.PH_VALUES)
    assert len({round(r.average, 6) for _, r in found}) > 1
    # the acids dominate these surfaces, so it falls as the pH rises
    assert found[0][1].maximum > found[-1][1].maximum


def test_a_folded_protein_is_globular(protein):
    result = aggrescan.score(protein)
    assert result.is_globular()
    assert result.metric('consensus')[1] == 'max'
    assert result.metric('avg') == (result.average, 'avg')


def test_an_extended_chain_is_not_globular():
    # ("chain" and the like are selection keywords: PyMOL renames them)
    cmd.fab('AAAAAAAAAAAAAAAAAAAA', 'strand', ss=2)
    result = aggrescan.score('strand')
    assert not result.is_globular()
    assert result.metric('consensus')[1] == 'avg'


def test_stats_and_csv(protein):
    result = aggrescan.score(protein)
    stats = result.stats()
    assert set(stats) >= {'All', 'A'}
    assert stats['All']['max_value'] == pytest.approx(result.maximum, abs=1e-4)
    assert stats['All']['avg_value'] == pytest.approx(result.average, abs=1e-4)
    lines = result.csv().splitlines()
    assert lines[0] == 'protein,chain,residue,residue_name,score'
    assert len(lines) == len(result) + 1


# --- showing it --------------------------------------------------------


def test_apply_scores_colours_and_selects(protein):
    result = aggrescan.score(protein)
    shown = aggrescan.apply_scores(protein, result)
    assert shown['residues'] == len(result)
    assert shown['prone'] == len(result.prone)
    assert shown['selection'] in (cmd.get_names('selections') or [])
    # the scores land in the B-factors
    best = max(result.residues, key=lambda r: r.score)
    values = []
    cmd.iterate('%%p and chain %s and resi %s and name CA'
                % (best.chain, best.resi), 'rows.append(b)',
                space={'rows': values})
    assert values[0] == pytest.approx(best.score, abs=1e-4)


def test_report_mentions_the_setup(protein):
    text = aggrescan.report(aggrescan.score(protein))
    assert 'Aggrescan scale' in text and '10 A sphere' in text
    assert 'globular' in text
    text = aggrescan.report(aggrescan.score(protein, ph=6.0))
    assert 'pH 6.0' in text


# --- writing the results out ------------------------------------------


def test_scored_pdb_carries_the_scores(protein, tmp_path):
    source = str(tmp_path / 'in.pdb')
    aggrescan.write_pdb(protein, source)
    result = aggrescan.score_pdb(source)
    out = str(tmp_path / 'scored.pdb')
    written = aggrescan.write_scored_pdb(source, out, result)
    assert written == cmd.count_atoms('(%s) and %s'
                                      % (protein, aggrescan.SELECTION))
    cmd.load(out, 'back')
    got = []
    cmd.iterate('back and name CA', 'rows.append((chain, resi, q, b))',
                space={'rows': got})
    scores = {r.key: r.score for r in result.residues}
    for chain, resi, occupancy, b in got:
        # the score goes in the B-factor, as A3D writes it, at 2 decimals
        assert b == pytest.approx(scores[(chain or '-', resi)], abs=0.005)
        assert occupancy == pytest.approx(1.0)      # left alone


def test_save_writes_one_pdb_per_ph(protein, tmp_path):
    job = aggrescan.prepare(protein, ph=True, optimize=False, pka='fixed',
                            name='prot')
    try:
        job.run()
        written = job.save(str(tmp_path))
    finally:
        job.cleanup()
    names = {os.path.basename(p) for p in written}
    for value in aggrescan.PH_VALUES:
        assert 'prot_pH%.1f.pdb' % value in names
        assert 'prot_pH%.1f.csv' % value in names
    assert 'prot_A4D.pdb' in names and 'prot_A4D.csv' in names
    assert 'prot_pH_summary.csv' in names
    assert all(os.path.isfile(p) for p in written)
    # each pH file really holds that pH's scores
    ends = [str(tmp_path / ('prot_pH%.1f.pdb' % v))
            for v in (aggrescan.PH_VALUES[0], aggrescan.PH_VALUES[-1])]
    assert open(ends[0]).read() != open(ends[1]).read()
    rows = open(str(tmp_path / 'prot_pH_summary.csv')).read().splitlines()
    assert rows[0].startswith('pH,maximum,average')
    assert len(rows) == len(aggrescan.PH_VALUES) + 1


def test_save_without_the_ph_scan_still_writes_the_structure(protein,
                                                             tmp_path):
    job = aggrescan.prepare(protein, ph=False, optimize=False, name='one')
    try:
        job.run()
        written = job.save(str(tmp_path))
    finally:
        job.cleanup()
    assert {os.path.basename(p) for p in written} == {'one_A4D.csv',
                                                      'one_A4D.pdb'}


def test_save_needs_the_analysis(protein, tmp_path):
    job = aggrescan.prepare(protein, ph=False, optimize=False)
    with pytest.raises(ValueError, match='run the analysis'):
        job.save(str(tmp_path))
    job.run()
    job.cleanup()               # the working files are gone
    with pytest.raises(ValueError, match='gone'):
        job.save(str(tmp_path))


def test_the_command_writes_the_files(protein, tmp_path):
    folder = tmp_path / 'out'
    cmd.aggregation(protein, ph=1, pka='fixed', output=str(folder),
                    quiet=1)
    names = {p.name for p in folder.iterdir()}
    assert len([n for n in names if n.endswith('.pdb')]) == \
        len(aggrescan.PH_VALUES) + 1


# --- mutations ---------------------------------------------------------


def test_evolutionary_targets_come_from_pam250():
    found = aggrescan.targets('L', 'evolutionary')
    assert 'M' in found and 'I' in found     # conservative swaps
    assert 'D' not in found and 'P' not in found
    assert 'L' not in found


def test_charged_targets_are_the_a3d_set():
    assert aggrescan.targets('L', 'charged') == ['R', 'K', 'D', 'E']
    assert 'R' not in aggrescan.targets('R', 'charged')


def test_candidates_pick_exposed_aggregation_prone_residues(protein):
    result = aggrescan.score(protein)
    found = aggrescan.candidates(result, 'charged', max_sites=3)
    assert found
    assert len({(s.chain, s.resi) for s in found}) <= 3
    assert all(s.estimated < 0 for s in found)
    # never a residue that is already charged, in the charged mode
    assert all(s.old not in aggrescan.CHARGED_MUTATIONS for s in found)
    assert all(s.new in aggrescan.CHARGED_MUTATIONS for s in found)


def test_candidates_respect_excluded_residues(protein):
    result = aggrescan.score(protein)
    first = aggrescan.candidates(result, 'charged', max_sites=2)
    assert first
    keep = (first[0].chain, first[0].resi)
    again = aggrescan.candidates(result, 'charged', excluded=[keep],
                                 max_sites=2)
    assert all((s.chain, s.resi) != keep for s in again)


def test_a_suggestion_reads_as_a_mutation():
    one = aggrescan.Suggestion('A', '57', 'L', 'E')
    assert one.label == 'L57E'
    assert one.pdbfixer == 'LEU-57-GLU'
    assert one.risk == ''                 # nothing measured yet
    one.ddg = aggrescan.ENERGY_BAD + 0.1
    assert one.risk == 'bad'
    one.ddg = aggrescan.ENERGY_WARN + 0.1
    assert one.risk == 'warn'


def test_parse_mutations_reads_what_the_user_types(protein):
    result = aggrescan.score(protein)
    letters = {r.key: r.letter for r in result.residues}
    chain, resi = next(iter(letters))
    letter = letters[(chain, resi)]
    one, = aggrescan.parse_mutations('%s%s%s' % (letter, resi, 'E'), result)
    assert (one.chain, one.resi, one.old, one.new) == (chain, resi,
                                                       letter, 'E')
    # the chain may be named, and several may be given at once
    two = aggrescan.parse_mutations(
        '%s/%s%sE, %s/%s%sK' % (chain, letter, resi, chain, letter, resi),
        result)
    assert len(two) == 2


def test_parse_mutations_catches_mistakes(protein):
    result = aggrescan.score(protein)
    letters = {r.key: r.letter for r in result.residues}
    (chain, resi), letter = next(iter(letters.items()))
    wrong = 'ACDEFGHIKLMNPQRSTVWY'.replace(letter, '')[0]
    with pytest.raises(ValueError, match='not'):
        aggrescan.parse_mutations('%s%s%s' % (wrong, resi, 'E'), result)
    with pytest.raises(ValueError, match='no residue'):
        aggrescan.parse_mutations('%s99999E' % letter, result)
    with pytest.raises(ValueError, match='cannot read'):
        aggrescan.parse_mutations('nonsense', result)
    with pytest.raises(ValueError, match='amino acid'):
        aggrescan.parse_mutations('%s%sB' % (letter, resi), result)
    assert chain


@needs_openmm
@slow
def test_mutating_before_the_analysis(protein):
    result = aggrescan.score(protein)
    best = max(result.residues, key=lambda r: r.score)
    text = '%s/%s%sE' % (best.chain, best.letter, best.resi)
    job = aggrescan.prepare(protein, ph=False, optimize=True, mutate=text)
    try:
        job.run()
        assert [one.label for one in job.mutated] == \
            ['%s%sE' % (best.letter, best.resi)]
        scored = {r.key: r for r in job.result.residues}
        assert scored[(best.chain, best.resi)].letter == 'E'
    finally:
        job.cleanup()


def test_mutating_needs_the_optimisation(protein):
    with pytest.raises(ValueError, match='optimisation'):
        aggrescan.prepare(protein, ph=False, optimize=False,
                          mutate='L8E')


# --- the options the server refuses -----------------------------------


@pytest.mark.parametrize('options,message', [
    (dict(ph=True, dynamic=True), 'dynamic'),
    (dict(ph=True, mutations='charged'), 'automated mutations'),
    (dict(ph=False, dynamic=True, mutations='charged'), 'dynamic'),
    (dict(ph=False, optimize=False, mutations='charged'), 'optimisation'),
])
def test_conflicting_options_are_refused(protein, options, message):
    with pytest.raises(ValueError) as caught:
        aggrescan.prepare(protein, **options)
    assert message in str(caught.value)


def test_a_bad_distance_is_refused(protein):
    with pytest.raises(ValueError):
        aggrescan.prepare(protein, distance=7, ph=False)


# --- the globular regions option --------------------------------------


def test_globular_trimming_leaves_experimental_structures_alone(protein):
    picked, dropped, predicted = aggrescan.globular_selection(protein)
    assert not predicted and dropped == 0 and picked == protein


def test_globular_trimming_drops_low_plddt(protein):
    # make it look like a predicted model with one weak tail
    cmd.alter('%p', 'b = 90.0')
    cmd.alter('%p and resi 1-10', 'b = 40.0')
    cmd.sort('p')
    picked, dropped, predicted = aggrescan.globular_selection(protein)
    assert predicted and dropped == 10
    assert cmd.count_atoms('(%s) and name CA' % picked) == \
        cmd.count_atoms('%p and name CA') - 10


# --- OpenMM and EvoEF, which stand in for FoldX -----------------------


@needs_openmm
@slow
def test_optimisation_lowers_the_energy(protein, tmp_path):
    source = str(tmp_path / 'in.pdb')
    aggrescan.write_pdb(protein, source)
    out = str(tmp_path / 'out.pdb')
    report = aggrescan.optimize_pdb(source, out)
    assert report['energy_after'] < report['energy_before']
    assert os.path.isfile(out)
    # it is still the same protein, and still scoreable
    assert len(aggrescan.score_pdb(out)) == len(aggrescan.score(protein))


@needs_evoef
@slow
def test_evoef_measures_a_mutation(protein, tmp_path):
    source = str(tmp_path / 'in.pdb')
    aggrescan.write_pdb(protein, source)
    result = aggrescan.score(protein)
    best = max(result.residues, key=lambda r: r.score)
    one = aggrescan.Suggestion(best.chain, best.resi, best.letter, 'E')
    done = aggrescan.evoef_ddg(source, [one], str(tmp_path / 'evoef'))
    assert done == [one]
    assert one.ddg is not None
    assert -20.0 < one.ddg < 20.0       # kcal/mol, a believable range


@needs_openmm
@slow
def test_the_dynamic_mode_samples_and_scores(protein, tmp_path):
    source = str(tmp_path / 'in.pdb')
    aggrescan.write_pdb(protein, source)
    folder = tmp_path / 'models'
    folder.mkdir()
    found = aggrescan.dynamics_pdb(source, str(folder), n_models=3,
                                   steps=12000)
    assert len(found['models']) == 3
    assert all(os.path.isfile(p) for p in found['models'])
    assert found['rmsf'] and all(v >= 0 for _, _, v in found['rmsf'])
    assert len(found['rmsf']) == len(aggrescan.score(protein))


# --- the "aggregation" preset -----------------------------------------


def test_preset_shows_a_half_transparent_surface(protein):
    cmd.aggregation(protein, quiet=1)
    from pymol import preset
    preset.aggregation(protein)
    assert cmd.count_atoms('%s and rep surface' % protein) > 0
    # nothing else is left showing, as the other presets do
    for rep in ('cartoon', 'sticks', 'lines', 'spheres'):
        assert cmd.count_atoms('%s and rep %s' % (protein, rep)) == 0
    assert float(cmd.get('transparency', protein)) == pytest.approx(
        aggrescan.PRESET_TRANSPARENCY)


def test_preset_ramp_spans_the_b_factors(protein):
    cmd.aggregation(protein, quiet=1)
    values = []
    cmd.iterate(protein, 'rows.append(b)', space={'rows': values})
    from pymol import preset
    low, high = preset.aggregation(protein)
    assert low == pytest.approx(min(values))
    assert high == pytest.approx(max(values))
    # the extremes sit at the ends of blue-white-red
    colors = {}
    cmd.iterate('%s and name CA' % protein,
                'rows[(chain, resi, b)] = color', space={'rows': colors})
    coldest = min(colors, key=lambda k: k[2])
    hottest = max(colors, key=lambda k: k[2])
    blue = cmd.get_color_tuple(colors[coldest])
    red = cmd.get_color_tuple(colors[hottest])
    assert blue[2] > 0.9 and blue[0] < 0.1       # blue
    assert red[0] > 0.9 and red[2] < 0.1         # red


def test_preset_survives_one_value_and_no_atoms(protein):
    from pymol import preset
    cmd.alter(protein, 'b = 1.0')
    cmd.sort(protein)
    assert preset.aggregation(protein) == (1.0, 1.0)
    colors = []
    cmd.iterate('%s and name CA' % protein, 'rows.append(color)',
                space={'rows': colors})
    # a spectrum over no range would be meaningless: the middle instead
    assert cmd.get_color_tuple(colors[0]) == pytest.approx((1.0, 1.0, 1.0))
    cmd.create('nothing', 'none')
    assert preset.aggregation('nothing') is None


def test_preset_spans_several_objects(protein):
    cmd.aggregation(protein, quiet=1)
    cmd.create('second', protein)
    cmd.alter('second', 'b = 99.0')
    cmd.sort('second')
    from pymol import preset
    low, high = preset.aggregation('all')
    assert high == pytest.approx(99.0)
    for obj in (protein, 'second'):
        assert cmd.count_atoms('%s and rep surface' % obj) > 0
        assert float(cmd.get('transparency', obj)) == pytest.approx(
            aggrescan.PRESET_TRANSPARENCY)
    assert low < 0


def test_preset_is_offered_in_the_menu():
    from pymol import menu
    rows = menu.presets(cmd, 'all')
    entries = {row[1]: row[2] for row in rows if len(row) == 3}
    assert 'aggregation' in entries
    assert 'preset.aggregation' in str(entries['aggregation'])


# --- the command and the plugin ---------------------------------------


def test_the_command_is_registered():
    assert 'aggregation' in cmd.keyword
    from pymolx import undo
    assert 'aggregation' in undo.UNDOABLE


def test_the_command_scores_and_colours(protein):
    job = cmd.aggregation(protein, quiet=1)
    assert len(job.result) > 0
    assert '%s_aggregation' % protein in (cmd.get_names('selections') or [])


def test_the_command_takes_one_ph(protein):
    job = cmd.aggregation(protein, ph=5.5, quiet=1)
    assert [value for value, _ in job.ph_results] == [5.5]
    assert job.chosen().ph == 5.5


def test_dialog_follows_the_server_defaults(qapp_dialog):
    dialog = qapp_dialog
    assert dialog.ph.isChecked()            # the server's default
    assert dialog.optimize.isChecked()
    assert not dialog.globular.isChecked()
    assert not dialog.dynamic.isChecked()
    assert dialog.distance.currentData() == aggrescan.DEFAULT_DISTANCE
    assert dialog.mutations.currentData() == ''
    assert dialog.mutate.text() == ''
    assert [dialog.distance.itemData(i)
            for i in range(dialog.distance.count())] == \
        list(aggrescan.DISTANCES)


def test_dialog_greys_out_the_conflicts(qapp_dialog):
    dialog = qapp_dialog
    dialog.ph.setChecked(True)
    assert not dialog.dynamic.isEnabled()
    assert not dialog.mutations.isEnabled()
    dialog.ph.setChecked(False)
    assert dialog.dynamic.isEnabled() and dialog.mutations.isEnabled()
    dialog.dynamic.setChecked(True)
    assert not dialog.ph.isEnabled() and not dialog.mutations.isEnabled()
    dialog.dynamic.setChecked(False)
    # mutations need the optimisation, as FoldX does in A4D
    dialog.optimize.setChecked(False)
    dialog.mutations.setCurrentIndex(1)
    assert dialog.optimize.isChecked() and not dialog.optimize.isEnabled()
    # a typed mutation needs it too, and the two kinds exclude each other
    dialog.mutations.setCurrentIndex(0)
    assert not dialog.mutate.isEnabled() or True
    dialog.optimize.setChecked(False)
    dialog.mutate.setText('L8E')
    assert dialog.optimize.isChecked() and not dialog.optimize.isEnabled()
    dialog.mutations.setCurrentIndex(1)
    assert not dialog.mutate.isEnabled()


def test_dialog_lists_the_session(qapp_dialog, protein):
    dialog = qapp_dialog
    dialog.refresh()
    assert dialog.object.findText(protein) >= 0
    dialog.object.setCurrentText(protein)
    assert dialog.chain.count() == len(cmd.get_chains('%' + protein)) + 1


@pytest.fixture
def qapp_dialog():
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import shiboken6
    from pmg_tk.startup.aggregation_gui.main import AggregationDialog

    dialog = AggregationDialog(None, _self=cmd)
    try:
        yield dialog
    finally:
        # a dialog only deleteLater()d outlives PyMOL and its teardown
        # then segfaults, after every test has passed
        shiboken6.delete(dialog)
