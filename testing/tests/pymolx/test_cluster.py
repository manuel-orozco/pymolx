'''
RMSD clustering (pymolx.cluster: cluster_rmsd) and the RMSD Clustering
dialog, on made-up antibody-antigen "models": two binding poses, each
model moved as a whole to somewhere else.
'''

import json
import os

import numpy as np
import pytest

from pymol import cmd
from pymolx import cluster, undo

VH = ('EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG'
      'RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS')
VL = ('DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRS'
      'GTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK')
ANTIGEN = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVK'


@pytest.fixture(autouse=True)
def session():
    from pymolx import files
    cmd.reinitialize()
    files._sources.clear()
    yield
    cmd.reinitialize()
    files._sources.clear()


def rotation(axis, degrees):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    t = np.radians(degrees)
    return np.eye(3) + np.sin(t) * k + (1 - np.cos(t)) * k @ k


def models(poses=((0, 5), (120, 3)), plddt=None):
    '''
    model_1 ... : antibody (H, L) around antigen (A), turned by each
    pose's angle (+ noise), the whole model then moved anywhere
    '''
    cmd.fab(VH, 'h', chain='H', ss=1)
    cmd.fab(VL, 'l', chain='L', ss=1)
    cmd.fab(ANTIGEN, 'a', chain='A', ss=1)
    cmd.translate([40, 0, 0], 'h')
    cmd.translate([40, 25, 0], 'l')
    center = cmd.get_coords('a').mean(axis=0)
    rng = np.random.default_rng(1)
    names = []
    for angle, count in poses:
        for _ in range(count):
            name = 'model_%d' % (len(names) + 1)
            cmd.create(name, 'h or l or a')
            ab = cmd.get_coords('%s and chain H+L' % name)
            turn = rotation([0, 0, 1], angle) @ rotation(rng.normal(size=3),
                                                         rng.normal(0, 0.7))
            cmd.load_coords((ab - center) @ turn.T + center +
                            rng.normal(0, 0.2, 3),
                            '%s and chain H+L' % name)
            whole = rotation(rng.normal(size=3), rng.uniform(0, 180))
            cmd.load_coords((cmd.get_coords(name) - center) @ whole.T +
                            rng.normal(0, 20, 3), name)
            if plddt is not None:
                cmd.alter(name, 'b = %r' % plddt[len(names)])
            names.append(name)
    cmd.delete('h l a')
    return names


def _settle(name, centre, want=3.5):
    """
    Move the antibody of a model along the line from the antigen's centre
    until its closest heavy atom is about "want" from the antigen
    """
    antibody = '%s and chain H+L' % name
    antigen = '%s and chain A' % name

    def closest():
        xa = cmd.get_coords(antigen + ' and not hydro')
        xb = cmd.get_coords(antibody + ' and not hydro')
        return float(np.linalg.norm(xa[:, None] - xb[None], axis=2).min())

    direction = cmd.get_coords(antibody).mean(axis=0) - centre
    direction = direction / np.linalg.norm(direction)
    for _ in range(300):
        gap = closest()
        if abs(gap - want) < 1.0:
            return
        step = direction * (1.0 if gap > want else -1.0)
        cmd.load_coords(cmd.get_coords(antibody) - step, antibody)


def touching_models(poses=((0, 3), (120, 3))):
    '''
    Like models(), but the antibody rests against the antigen, so the
    two have contacts: each pose puts it on a different side
    '''
    cmd.fab(VH, 'h', chain='H', ss=1)
    cmd.fab(VL, 'l', chain='L', ss=1)
    cmd.fab(ANTIGEN, 'a', chain='A', ss=1)
    cmd.create('ab', 'h or l')
    cmd.delete('h l')
    # centroids together, then apart along x until they just touch
    cmd.translate(list(cmd.get_coords('a').mean(axis=0) -
                       cmd.get_coords('ab').mean(axis=0)), 'ab', camera=0)

    def closest():
        xa = cmd.get_coords('a and not hydro')
        xb = cmd.get_coords('ab and not hydro')
        return float(np.linalg.norm(xa[:, None] - xb[None], axis=2).min())

    for _ in range(200):
        if closest() >= 3.0:
            break
        cmd.translate([1.0, 0, 0], 'ab', camera=0)
    assert closest() < 5.0, 'the test models should touch'
    centre = cmd.get_coords('a').mean(axis=0)
    rng = np.random.default_rng(3)
    names = []
    for angle, count in poses:
        for _ in range(count):
            name = 'model_%d' % (len(names) + 1)
            cmd.create(name, 'ab or a')
            coords = cmd.get_coords('%s and chain H+L' % name)
            turn = rotation([0, 0, 1], angle) @ rotation(
                rng.normal(size=3), rng.normal(0, 0.5))
            cmd.load_coords((coords - centre) @ turn.T + centre,
                            '%s and chain H+L' % name)
            _settle(name, centre)
            whole = rotation(rng.normal(size=3), rng.uniform(0, 180))
            cmd.load_coords((cmd.get_coords(name) - centre) @ whole.T +
                            rng.normal(0, 20, 3), name)
            names.append(name)
    cmd.delete('ab a')
    return names


def test_measures():
    # the same models seen three ways: the antibody's position, the
    # contacting residues, and which antigen residues are touched
    names = touching_models(((0, 3), (120, 3)))
    seen = {}
    for measure in cluster.MEASURES:
        cutoff = 0.5 if measure == 'epitope' else 3.0
        result = cluster.compute(cutoff=cutoff, measure=measure)
        seen[measure] = sorted(sorted(c['members']) for c in result.clusters)
        assert result.measure == measure
    # each finds the two poses
    for measure, clusters in seen.items():
        assert len(clusters) == 2, measure
        assert sorted(len(c) for c in clusters) == [3, 3], measure
    assert seen['ligand'] == seen['epitope'] == seen['interface']
    with pytest.raises(ValueError, match='measure must be'):
        cluster.compute(measure='nonsense')


def test_epitope_distances():
    # the distance is how much of the contacted patch is not shared
    names = touching_models(((0, 2), (120, 2)))
    result = cluster.compute(cutoff=0.5, measure='epitope')
    assert len(result.epitopes) == len(names)
    assert all(e for e in result.epitopes)          # every model touches
    assert all(chain == 'A' for e in result.epitopes for chain, _ in e)
    same = result.matrix[0, 1]
    apart = result.matrix[0, len(names) - 1]
    assert 0 <= same < apart <= 1
    # two models on the same patch share most of it
    assert same < 0.5 < apart


def test_interface_rmsd_ignores_distal_domains():
    # a model whose antibody is turned far from the antigen, with the
    # epitope untouched: ligand RMSD sees it, interface RMSD much less
    names = touching_models(((0, 3),))
    before = cluster.compute(cutoff=99.0, measure='interface').matrix
    far = names[-1]
    centre = cmd.get_coords('%s and chain A' % far).mean(axis=0)
    coords = cmd.get_coords('%s and chain H+L' % far)
    turned = rotation([0, 0, 1], 25)
    cmd.load_coords((coords - centre) @ turned.T + centre,
                    '%s and chain H+L' % far)
    ligand = cluster.compute(cutoff=99.0, measure='ligand').matrix
    interface = cluster.compute(cutoff=99.0, measure='interface').matrix
    assert ligand[0, len(names) - 1] > interface[0, len(names) - 1]
    assert before.shape == interface.shape


def test_chain_pair_scores(tmp_path):
    # AlphaFold 3's chain_pair_iptm, read for the antibody-antigen pair:
    # its own documentation suggests it for such an interface, rather
    # than the whole complex's ipTM
    data = {'iptm': 0.9, 'ptm': 0.8, 'ranking_score': 0.85,
            # chains A (antigen), H, L: A-H 0.4, A-L 0.3, H-L 0.95
            'chain_pair_iptm': [[0.7, 0.4, 0.3],
                                [0.4, 0.8, 0.95],
                                [0.3, 0.95, 0.8]],
            'chain_pair_pae_min': [[0.5, 9.0, 12.0],
                                   [9.0, 0.5, 1.2],
                                   [12.0, 1.2, 0.5]]}
    path = tmp_path / 'scores.json'
    path.write_text(json.dumps(data))
    scores = cluster.read_scores(str(path), ['A', 'H', 'L'], ['H', 'L'])
    # the best antibody-antigen pair, not the heavy-light one
    assert scores['ab_ag_iptm'] == pytest.approx(0.4)
    assert scores['ab_ag_pae'] == pytest.approx(9.0)
    assert scores['iptm'] == pytest.approx(0.9)
    # without the chains, only the whole-complex scores
    assert 'ab_ag_iptm' not in cluster.read_scores(str(path))
    # a mismatched chain list is ignored rather than guessed at
    assert 'ab_ag_iptm' not in cluster.read_scores(
        str(path), ['A', 'B'], ['B'])
    # lower is better for a PAE, higher for the rest
    assert cluster.higher_is_better('ab_ag_iptm')
    assert not cluster.higher_is_better('ab_ag_pae')
    assert cluster.higher_is_better('ranking')


def test_clusters_ranked_by_size_and_confidence():
    names = touching_models(((0, 4), (120, 2)))
    result = cluster.compute(cutoff=3.0)
    assert len(result.clusters) == 2
    big, small = result.clusters
    assert big['size'] == 4 and small['size'] == 2
    # with no chain-pair scores, the order is by size
    assert [c['cluster'] for c in result.ranked()] == [1, 2]
    # the small cluster, confident; the large one, not
    for i, model in enumerate(result.models):
        value = 0.9 if model in small['members'] else 0.2
        result.confidences[i]['ab_ag_iptm'] = value
    again = cluster.Result(result.models, result.matrix, result.labels,
                           result.scores, 1, 1, result.confidences)
    big, small = again.clusters
    assert big['mean_ab_ag_iptm'] == pytest.approx(0.2)
    assert small['mean_ab_ag_iptm'] == pytest.approx(0.9)
    # 4/4 x 0.2 = 0.20 against 2/4 x 0.9 = 0.45: the confident one leads
    assert [c['cluster'] for c in again.ranked()] == [2, 1]
    assert 'Ab-Ag ipTM' in again.table()


def test_antibody_chains_detected():
    models(((0, 1),))
    assert cluster.antibody_chains('model_1') == ['H', 'L']
    assert cluster.default_selections('model_1') == ('chain A', 'chain H+L')
    cmd.fab('ACDEFGHIK', 'pep')
    assert cluster.antibody_chains('pep') == []
    assert cluster.default_selections('pep') == ('polymer', 'polymer')


@pytest.mark.parametrize('method', cluster.METHODS)
def test_two_binding_modes(method):
    models()
    result = cluster.compute('model_*', method=method)
    assert result.labels == [0] * 5 + [1] * 3
    assert [c['size'] for c in result.clusters] == [5, 3]
    assert result.fit == 'chain A' and result.rmsd == 'chain H+L'
    # same pose: small RMSD despite the whole-model moves
    assert result.matrix[:5, :5].max() < 3.0
    assert result.matrix[0, 5] > 20
    assert result.clusters[0]['representative'] in \
        result.clusters[0]['members']


def test_command_superposes_groups_and_colors(capsys):
    names = models()
    result = cmd.cluster_rmsd('model_*')
    out = capsys.readouterr().out
    assert '8 structures, 2 clusters' in out and 'Representative' in out
    reference = result.clusters[0]['representative']
    for name in names:
        # every antigen on the reference's
        assert cmd.rms_cur('%s and chain A and name CA' % name,
                           '%s and chain A and name CA' % reference,
                           matchmaker=4) == pytest.approx(0, abs=1e-3)
    assert cmd.get_names('group_objects') == ['cluster_1', 'cluster_2']
    members = cmd.get_object_list('cluster_2')
    assert sorted(members) == ['model_6', 'model_7', 'model_8']
    colors = set()
    cmd.iterate('model_6 and chain H and elem C', 'colors.add(color)',
                space={'colors': colors})
    assert colors == {cmd.get_color_index(cluster.CLUSTER_COLORS[1])}


def test_scores_and_csv(tmp_path):
    models(((0, 3), (120, 2)), plddt=[90, 80, 70, 50, 60])
    result = cluster.compute('model_*')
    assert result.clusters[0]['score'] == pytest.approx(80)
    assert result.clusters[1]['score'] == pytest.approx(55)
    path = str(tmp_path / 'c.csv')
    result.save_csv(path)
    lines = open(path).read().splitlines()
    assert lines[0].startswith('model,cluster,representative,best,pLDDT,'
                               'ipTM,pTM,ranking,model_1')
    assert len(lines) == 6
    assert lines[4].startswith('model_4,2,')


def write_af3(tmp_path, names, iptm, ptm):
    '''
    Save models like AlphaFold 3 output (model + summary_confidences)
    '''
    import json
    paths = []
    for name, i, p in zip(names, iptm, ptm):
        stem = 'job_%s' % name.replace('model_', 'seed-1_sample-')
        path = tmp_path / ('%s_model.cif' % stem)
        cmd.save(str(path), name)
        (tmp_path / ('%s_summary_confidences.json' % stem)).write_text(
            json.dumps({'iptm': i, 'ptm': p, 'ranking_score': i * 0.8 + 0.2,
                        'chain_pair_iptm': [[0.5]]}))
        paths.append(str(path))
    return paths


def test_best_model_from_alphafold3_scores(tmp_path):
    names = models()
    iptm = [0.61, 0.83, 0.70, 0.55, 0.62, 0.40, 0.48, 0.45]
    ptm = [0.70, 0.71, 0.90, 0.60, 0.65, 0.50, 0.52, 0.66]
    paths = write_af3(tmp_path, names, iptm, ptm)
    cmd.delete('all')
    for path in paths:
        cmd.load(path)  # (pymolx remembers the files)
    result = cmd.cluster_rmsd('job_*', quiet=1)
    best = [(c['best'], c['best_by'], c['best_score'])
            for c in result.clusters]
    assert best == [
        ('job_seed-1_sample-2_model', 'iptm', pytest.approx(0.83)),
        ('job_seed-1_sample-7_model', 'iptm', pytest.approx(0.48))]
    # listed first in its group
    assert cmd.get_object_list('cluster_1')[0] == 'job_seed-1_sample-2_model'
    by_ptm = cluster.compute('job_*', best_by='ptm')
    assert by_ptm.clusters[0]['best'] == 'job_seed-1_sample-3_model'
    assert 'job_seed-1_sample-2_model (ipTM 0.830)' in result.table()


def test_scripts_are_not_sources(tmp_path):
    from pymolx import files
    models(((0, 1),))
    cmd.save(str(tmp_path / 'one_model.cif'), 'model_1')
    cmd.delete('all')
    script = tmp_path / 'run.py'
    script.write_text('from pymol import cmd\ncmd.load(%r)\n' %
                      str(tmp_path / 'one_model.cif'))
    cmd.load(str(script))
    assert files.source('one_model') == str(tmp_path / 'one_model.cif')


def test_scores_found_in_a_folder(tmp_path):
    names = models(((0, 2), (120, 2)))
    write_af3(tmp_path, names, [0.5, 0.9, 0.3, 0.2], [0.5] * 4)
    # objects named like the files, but pymolx doesn't know the files
    for name in names:
        cmd.set_name(name, 'job_' + name.replace('model_', 'seed-1_sample-')
                     + '_model')
    result = cluster.compute('job_*')
    assert [c['best'] for c in result.clusters] == [None, None]
    result = cluster.compute('job_*', scores=str(tmp_path))
    assert [c['best'] for c in result.clusters] == [
        'job_seed-1_sample-2_model', 'job_seed-1_sample-3_model']


def test_score_file_names_and_fallback(tmp_path):
    import json
    names = cluster._score_names
    assert 'job_seed-1_sample-0_summary_confidences.json' in \
        names('job_seed-1_sample-0_model')
    assert 'confidence_x_model_0.json' in names('x_model_0')       # Boltz
    assert 'x_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json' \
        in names('x_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_'
                 'seed_000')                                         # ColabFold
    path = tmp_path / 'confidence_b.json'
    path.write_text(json.dumps({'confidence_score': 0.7, 'ptm': 0.8,
                                'iptm': 0.75}))
    assert cluster.read_scores(str(path)) == {
        'iptm': 0.75, 'ptm': 0.8, 'ranking': 0.7}
    # no ipTM anywhere: pLDDT (B-factors) picks the best
    models(((0, 2),), plddt=[70, 90])
    result = cluster.compute('model_*')
    assert result.clusters[0]['best'] == 'model_2'
    assert result.clusters[0]['best_by'] == 'plddt'
    with pytest.raises(ValueError, match='best_by'):
        cluster.compute('model_*', best_by='tm')


def test_short_names():
    names = ['ab_mFolR2_seed-10_sample-2_model',
             'ab_mFolR2_seed-1_sample-0_model',
             'ab_mFolR2_seed-3_sample-4_model']
    assert cluster.short_names(names) == {
        names[0]: 'seed-10_sample-2', names[1]: 'seed-1_sample-0',
        names[2]: 'seed-3_sample-4'}
    assert cluster.short_names(['model_1', 'model_2']) == {
        'model_1': '1', 'model_2': '2'}
    assert cluster.short_names(['only']) == {'only': 'only'}


# AlphaFold 3's <job>_ranking_scores.csv (rows by seed and sample)
AF3_CSV = """seed,sample,ranking_score
1,0,0.8101264229978932
1,1,0.7947506092984028
1,2,0.8051870378874295
1,3,0.7986404397709234
1,4,0.7979878862894813
2,0,0.7921657089592417
2,1,0.7946472074278781
2,2,0.785369856135175
"""
JOB = 'mIgG2aL-TMA641-me-001-D12_mFolR2'


def af3_models(job=JOB):
    '''
    8 models named like AlphaFold 3's, seed 1 samples 0-4 in one binding
    mode, seed 2 samples 0-2 in another
    '''
    names = models()
    keys = ['seed-1_sample-%d' % i for i in range(5)] + \
        ['seed-2_sample-%d' % i for i in range(3)]
    renamed = []
    for name, key in zip(names, keys):
        new = '%s_%s_model' % (job, key)
        cmd.set_name(name, new)
        renamed.append(new)
    return renamed


def test_alphafold3_ranking_scores_csv(tmp_path):
    path = tmp_path / ('%s_ranking_scores.csv' % JOB)
    path.write_text(AF3_CSV)
    table = cluster.read_scores_csv(str(path))
    assert table['seed-1_sample-2'] == {'ranking': pytest.approx(0.80519,
                                                                 abs=1e-5)}
    names = af3_models()
    result = cmd.cluster_rmsd(' '.join(names), scores_csv=str(path),
                              best_by='ranking', quiet=1)
    best = [(c['best'], c['best_by']) for c in result.clusters]
    assert best == [('%s_seed-1_sample-0_model' % JOB, 'ranking'),
                    ('%s_seed-2_sample-1_model' % JOB, 'ranking')]
    assert result.unmatched_scores == []
    # without asking for it: ipTM is missing, the ranking score is used
    result = cluster.compute(' '.join(names), scores_csv=str(path))
    assert result.clusters[0]['best_by'] == 'ranking'


def test_csv_rows_match_whole_words_and_their_job(tmp_path):
    table = {('', 'seed-1_sample-0'): {'ranking': 1.0}}
    models_ = ['job_seed-11_sample-0_model', 'job_seed-1_sample-01_model',
               'job_seed-1_sample-0_model']
    found, unmatched = cluster.match_scores(models_, table)
    assert list(found) == ['job_seed-1_sample-0_model'] and not unmatched
    # two jobs loaded: each CSV only for its own job's models
    two = ['a_seed-1_sample-0_model', 'b_seed-1_sample-0_model']
    found, _ = cluster.match_scores(two, {('b', 'seed-1_sample-0'): {'x': 1}})
    assert list(found) == ['b_seed-1_sample-0_model']
    found, unmatched = cluster.match_scores(
        two, {('', 'seed-1_sample-0'): {'x': 1}})
    assert found == {} and unmatched == ['seed-1_sample-0']  # ambiguous


def test_csv_with_names_and_any_score(tmp_path):
    names = models()
    path = tmp_path / 'scores.csv'
    path.write_text('model;dG;iptm\n' + ''.join(
        '%s.cif;%d;%.2f\n' % (n, -10 - k, 0.5 + 0.01 * k)
        for k, n in enumerate(names)))  # ";" and file names work too
    # dG: lower is better (by its name)
    result = cluster.compute('model_*', scores_csv=str(path), best_by='dg')
    assert result.clusters[0]['best'] == 'model_5'
    assert result.clusters[0]['best_score'] == -14
    result = cluster.compute('model_*', scores_csv=str(path), best_by='dg',
                             higher=True)
    assert result.clusters[0]['best'] == 'model_1'
    # the CSV's ipTM
    result = cluster.compute('model_*', scores_csv=str(path))
    assert result.clusters[1]['best'] == 'model_8'
    with pytest.raises(ValueError, match='best_by must be one of.*dg'):
        cluster.compute('model_*', scores_csv=str(path), best_by='ddg')
    out = tmp_path / 'out.csv'
    result = cluster.compute('model_*', scores_csv=str(path), best_by='dg')
    result.save_csv(str(out))
    assert out.read_text().splitlines()[0].startswith(
        'model,cluster,representative,best,pLDDT,ipTM,pTM,ranking,dg,')


def test_clustering_again_keeps_the_models():
    names = models()
    cmd.cluster_rmsd('model_*', quiet=1)
    cmd.cluster_rmsd('model_*', cutoff=500, quiet=1)  # now one cluster
    assert sorted(cmd.get_object_list('all')) == sorted(names)
    assert cmd.get_names('group_objects') == ['cluster_1']
    # a group of the same name with other things in it stays
    cmd.fab('AA', 'mine')
    cmd.group('cluster_7', 'mine')
    cmd.cluster_rmsd('model_*', quiet=1)
    assert cmd.get_object_list('cluster_7') == ['mine']


def test_errors():
    models(((0, 1),))
    with pytest.raises(ValueError, match='at least two'):
        cluster.compute('model_*')
    models(((0, 2),))
    with pytest.raises(ValueError, match='method'):
        cluster.compute('model_*', method='ward')
    with pytest.raises(ValueError, match='fewer than 3 atoms'):
        cluster.compute('model_*', fit='chain Z')


def test_one_undo_step():
    models()
    before = {n: cmd.get_coords(n).copy() for n in cmd.get_names()}
    cmd.undo_enable()
    try:
        steps = len(undo.stack.undo_stack)
        cmd.cluster_rmsd('model_*', quiet=1)
        assert len(undo.stack.undo_stack) == steps + 1
        cmd.undo()
        assert 'cluster_1' not in cmd.get_names('all')
        for name, xyz in before.items():
            assert cmd.get_coords(name) == pytest.approx(xyz)
    finally:
        cmd.undo_disable()


def test_models_read_in_one_pass():
    names = touching_models()
    cmd.remove('model_2 and chain A and resi 5-8')
    cmd.remove('model_3 and chain H and resi 10 and not name CA')
    cmd.fab('ACDEF', 'bystander', chain='A')
    antigens = cluster.read_sides(names, 'chain A', 1, ('CA',))
    antibodies = cluster.read_sides(names, 'chain H+L')
    assert list(antigens) == names and 'bystander' not in antigens
    for name in names:
        side = antigens[name]
        # the heavy atoms, in PyMOL's order, each with its residue
        assert np.allclose(side.xyz, cmd.get_coords(
            '%s and chain A and not hydro' % name))
        residues = []
        cmd.iterate('%s and chain A and not hydro' % name,
                    'residues.append((chain, resi))',
                    space={'residues': residues})
        assert [side.residues[k] for k in side.residue] == residues
        # the atoms to compare: its alpha carbons
        assert len(side.keys) == cmd.count_atoms(
            '%s and chain A and name CA' % name)
        assert np.allclose(side.coords, cmd.get_coords(
            '%s and chain A and name CA' % name))
        assert all(key[2] == 'CA' for key in side.keys)
        # the contacts are those PyMOL finds, both ways
        assert cluster.close_residues(side, antibodies[name]) == \
            cluster.contact_residues(name, 'chain A', 'chain H+L', 1)
        assert cluster.close_residues(antibodies[name], side, 8.0) == \
            cluster.contact_residues(name, 'chain H+L', 'chain A', 1, 8.0)
    assert antibodies[names[0]].coords is antibodies[names[0]].xyz
    # nothing of the selection in a model: an empty side, no contacts
    empty = cluster.read_sides(names, 'chain Z')[names[0]]
    assert len(empty.xyz) == 0 and empty.keys == []
    assert cluster.close_residues(empty, antibodies[names[0]]) == set()
    # atoms only some models have are left out, the rest matched
    keys, xyz = cluster.matched_atoms(
        [(antigens[n].keys, antigens[n].coords) for n in names])
    assert xyz.shape == (len(names), len(keys), 3)
    assert len(keys) == len(ANTIGEN) - 4
    assert not any(resi in ('5', '6', '7', '8') for _, resi, _ in keys)
    for n, name in enumerate(names):
        row = dict(zip(antigens[name].keys, antigens[name].coords))
        assert all(np.array_equal(xyz[n, k], row[key])
                   for k, key in enumerate(keys))
    both, joined = cluster.joined(antigens[names[2]], antibodies[names[2]],
                                  {('H', '10'), ('A', '3')})
    assert sorted(both) == [('A', '3', 'CA'), ('H', '10', 'CA')]
    assert joined.shape == (2, 3)


def test_matrices_as_pair_by_pair():
    rng = np.random.default_rng(5)
    n, atoms = 23, 40
    base = rng.normal(size=(atoms, 3)) * 12 + [60, -30, 90]
    fit = base + rng.normal(size=(n, atoms, 3)) * rng.uniform(
        0, 2, size=(n, 1, 1))
    measure = rng.normal(size=(n, 17, 3)) * 20
    fit[4], measure[4] = fit[2], measure[2]     # the same model twice
    slow = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            turn, shift = cluster.kabsch(fit[j], fit[i])
            moved = measure[j] @ turn.T + shift
            slow[i, j] = slow[j, i] = np.sqrt(
                ((moved - measure[i]) ** 2).sum(axis=1).mean())
    fast = cluster.rmsd_matrix(list(fit), list(measure))
    assert np.allclose(fast, slow, rtol=0, atol=1e-10)
    assert fast[2, 4] < 1e-10 and np.array_equal(fast, fast.T)
    assert cluster.rmsd_matrix(fit[:1], measure[:1]).tolist() == [[0.0]]
    # a few models at a time gives the same
    whole = cluster.BATCH
    cluster.BATCH = 3 * 17 * 4
    try:
        assert cluster.batches(1, 11, 17) == [(1, 5), (5, 9), (9, 11)]
        assert np.allclose(cluster.rmsd_matrix(fit, measure), slow,
                           rtol=0, atol=1e-10)
    finally:
        cluster.BATCH = whole
    # a reflection is not a rotation
    mirrored = fit[0] * [1, 1, -1]
    turn = cluster.rotations((mirrored - mirrored.mean(axis=0))[None],
                             fit[0] - fit[0].mean(axis=0))[0]
    assert np.linalg.det(turn) == pytest.approx(1)
    assert np.allclose(turn, cluster.kabsch(mirrored, fit[0])[0])

    sets = [set(rng.choice(30, size=rng.integers(0, 12), replace=False))
            for _ in range(n)]
    sets[3], sets[6], sets[7] = set(), set(), set(sets[1])
    slow = np.array([[0.0 if i == j else
                      1.0 - len(a & b) / len(a | b) if a | b else 1.0
                      for j, b in enumerate(sets)]
                     for i, a in enumerate(sets)])
    assert np.array_equal(cluster.jaccard_matrix(sets), slow)

    def gromos(matrix, cutoff):     # (as Daura et al. describe it)
        left, labels, number = set(range(len(matrix))), {}, 0
        while left:
            members = sorted(left)
            counts = [sum(matrix[i, j] <= cutoff for j in members)
                      for i in members]
            center = members[counts.index(max(counts))]
            for j in members:
                if matrix[center, j] <= cutoff:
                    labels[j] = number
                    left.discard(j)
            number += 1
        return [labels[i] for i in range(len(matrix))]

    for cutoff in (0.0, 0.8, 1.5, 2.5, 100.0):
        assert cluster.gromos(fast, cutoff) == gromos(fast, cutoff)
    assert cluster.gromos(np.zeros((0, 0)), 1.0) == []


def test_many_models_take_seconds():
    '''
    The time grows with the number of models, not with its square (as it
    did when PyMOL was asked about one model after the other)
    '''
    import time
    names = touching_models(((0, 2), (120, 2)))
    for copy in range(60):
        cmd.create('copy_%02d' % copy, names[copy % 4])
    start = time.time()
    result = cluster.compute('model_* copy_*', 'chain A', 'chain H+L')
    cluster.apply(result)
    elapsed = time.time() - start
    assert sorted(c['size'] for c in result.clusters) == [32, 32]
    # (a tenth of this here; a minute before)
    assert elapsed < 10
    # moved as coordinates: the objects keep no matrix of their own
    assert cmd.get_object_matrix('copy_07') == (
        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0)
    spread = np.ptp([cmd.get_coords('%s and chain A and name CA' % m)
                     for m in result.models], axis=0)
    assert spread.max() < 1e-2      # all on one antigen
    colors = set()
    cmd.iterate('chain H+L and elem C', 'colors.add(color)',
                space={'colors': colors})
    assert colors == {cmd.get_color_index(c)
                      for c in cluster.CLUSTER_COLORS[:2]}


def test_dialog(tmp_path):
    dialog_tmp = tmp_path
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import pymol.plugins  # noqa: F401 (adds the startup plugin path)
    from pmg_tk.startup.cluster_gui.main import ClusterDialog

    names = models()
    cmd.fab('ACD', 'other')
    dialog = ClusterDialog(None, _self=cmd)
    try:
        dialog.refresh()
        dialog.filter.setText('model')
        assert dialog.checked_models() == names
        dialog.detect()
        assert dialog.fit.text() == 'chain A'
        assert dialog.rmsd.text() == 'chain H+L'
        dialog.run()
        assert dialog.table.rowCount() == 2
        assert dialog.table.item(0, 1).text() == '5'
        assert dialog.save_button.isEnabled()
        assert dialog.table.item(0, 4).text() == '-'  # (no scores)
        assert dialog.table.item(0, 6).text() == \
            dialog.result.clusters[0]['representative'].split('_')[1]
        # a cluster shown alone
        dialog.table.selectRow(1)
        assert set(cmd.get_object_list('enabled and model_*')) == \
            {'model_6', 'model_7', 'model_8'}
        dialog.show_all()
        assert len(cmd.get_object_list('enabled and model_*')) == 8
        assert dialog.heatmap.order[:5] != [] and \
            len(dialog.heatmap.order) == 8
        dialog.heatmap.resize(160, 160)
        image = dialog.heatmap.grab()
        assert not image.isNull()
        # a scores CSV: its columns to choose from, "lower is better" by name
        path = os.path.join(str(dialog_tmp), 'scores.csv')
        with open(path, 'w') as handle:
            handle.write('model,ranking_score,energy\n' + ''.join(
                '%s,%.2f,%d\n' % (n, 0.5 + 0.01 * k, k)
                for k, n in enumerate(names)))
        dialog.scores_csv.setText(path)
        dialog._csv_changed()
        assert dialog.best_by.currentData() == 'ranking'
        assert dialog.best_by.findData('energy') >= 0
        assert not dialog.lower_better.isChecked()
        dialog.best_by.setCurrentIndex(dialog.best_by.findData('energy'))
        assert dialog.lower_better.isChecked()
        dialog.run()
        assert dialog.table.item(0, 5).text() == 'energy 0.000'
    finally:
        dialog.deleteLater()
