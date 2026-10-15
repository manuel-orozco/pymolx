'''
RMSD clustering (pymolx.cluster: cluster_rmsd) and the RMSD Clustering
dialog, on made-up antibody-antigen "models": two binding poses, each
model moved as a whole to somewhere else.
'''

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
        assert not dialog.lower.isChecked()
        dialog.best_by.setCurrentIndex(dialog.best_by.findData('energy'))
        assert dialog.lower.isChecked()
        dialog.run()
        assert dialog.table.item(0, 5).text() == 'energy 0.000'
    finally:
        dialog.deleteLater()
