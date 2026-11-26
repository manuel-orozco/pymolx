'''
Interface RMSD clustering (pymolx.irmsd: cluster_irmsd) and the Interface
RMSD Clustering dialog, on made-up antibody-antigen "models": the
antibody resting against one side of the antigen or another, each model
then moved as a whole to somewhere else.
'''

import os

import numpy as np
import pytest

from pymol import cmd
from pymolx import cluster, irmsd, undo

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


def gap(first, second):
    xa = cmd.get_coords('(%s) and not hydro' % first)
    xb = cmd.get_coords('(%s) and not hydro' % second)
    return float(np.linalg.norm(xa[:, None] - xb[None], axis=2).min())


def settle(name, centre, want=3.5):
    """
    Move the antibody of a model along the line from the antigen's centre
    until its closest heavy atom is about "want" from the antigen
    """
    antibody = '%s and chain H+L' % name
    direction = cmd.get_coords(antibody).mean(axis=0) - centre
    direction = direction / np.linalg.norm(direction)
    for _ in range(300):
        apart = gap('%s and chain A' % name, antibody)
        if abs(apart - want) < 1.0:
            return
        step = direction * (1.0 if apart > want else -1.0)
        cmd.load_coords(cmd.get_coords(antibody) - step, antibody)


def touching_models(poses=((0, 3), (120, 3))):
    '''
    model_1 ...: antibody (H, L) resting against the antigen (A), on the
    side each pose's angle gives (with a little noise), the whole model
    then moved anywhere
    '''
    cmd.fab(VH, 'vh', chain='H', ss=1)
    cmd.fab(VL, 'vl', chain='L', ss=1)
    cmd.fab(ANTIGEN, 'ag', chain='A', ss=1)
    cmd.create('ab', 'vh or vl')
    cmd.delete('vh vl')
    cmd.translate(list(cmd.get_coords('ag').mean(axis=0) -
                       cmd.get_coords('ab').mean(axis=0)), 'ab', camera=0)
    for _ in range(200):
        if gap('ag', 'ab') >= 3.0:
            break
        cmd.translate([1.0, 0, 0], 'ab', camera=0)
    centre = cmd.get_coords('ag').mean(axis=0)
    rng = np.random.default_rng(3)
    names = []
    for angle, count in poses:
        for _ in range(count):
            name = 'model_%d' % (len(names) + 1)
            cmd.create(name, 'ab or ag')
            coords = cmd.get_coords('%s and chain H+L' % name)
            turn = rotation([0, 0, 1], angle) @ rotation(
                rng.normal(size=3), rng.normal(0, 0.5))
            cmd.load_coords((coords - centre) @ turn.T + centre,
                            '%s and chain H+L' % name)
            settle(name, centre)
            whole = rotation(rng.normal(size=3), rng.uniform(0, 180))
            cmd.load_coords((cmd.get_coords(name) - centre) @ whole.T +
                            rng.normal(0, 20, 3), name)
            names.append(name)
    cmd.delete('ab ag')
    return names


def members(result):
    return sorted(sorted(c['members']) for c in result.clusters
                  if not c.get('unclustered'))


@pytest.mark.parametrize('method', cluster.METHODS)
def test_two_binding_modes(method):
    names = touching_models(((0, 3), (120, 3)))
    result = irmsd.compute('model_*', method=method)
    assert members(result) == [names[:3], names[3:]]
    assert result.fit == 'chain A' and result.rmsd == 'chain H+L'
    matrix = result.matrix
    assert matrix == pytest.approx(matrix.T)
    assert np.diag(matrix) == pytest.approx(0)
    # the same pose: within CAPRI's "acceptable"; another: far beyond
    assert matrix[0, 1] < irmsd.CUTOFF < matrix[0, 5]
    assert all(c['mean_rmsd'] < irmsd.CUTOFF for c in result.clusters)
    assert result.skipped == []


def test_interface_is_both_sides_within_the_distance():
    names = touching_models(((0, 2),))
    found = irmsd.interface_residues(names[0], 'chain A', 'chain H+L')
    assert {chain for chain, _ in found} >= {'A'}
    assert {chain for chain, _ in found} & {'H', 'L'}
    close = irmsd.interface_residues(names[0], 'chain A', 'chain H+L',
                                     distance=5.0)
    assert close and close < found
    # a wider interface is compared over more atoms
    wide = irmsd.compute('model_*', distance=10.0)
    narrow = irmsd.compute('model_*', distance=5.0)
    assert narrow.atoms_compared[1] < wide.atoms_compared[1]
    assert wide.interfaces[0] == found


def test_only_the_interface_counts():
    # the ends of the antibody that reach past the antigen pushed 5 A
    # aside in one model: still out of the interface's reach, so the
    # interface RMSD is what it was, while the antibody's own RMSD is not
    names = touching_models(((0, 2),))
    before = irmsd.compute('model_*')
    near = set()
    for name in names:
        near |= cluster.contact_residues(name, 'chain H+L', 'chain A', 1,
                                         irmsd.INTERFACE + 6.0)
    everything = cluster.contact_residues(names[1], 'chain H+L',
                                          'chain H+L', 1, 1.0)
    far = everything - near
    assert len(far) > 20, 'the test antibody should reach past the antigen'
    target = '%s and (%s)' % (names[1], cluster._residues_selection(far))
    ligand = cluster.compute('model_*', cutoff=99.0).matrix[0, 1]
    cmd.load_coords(cmd.get_coords(target) + [0.0, 0.0, 5.0], target)
    after = irmsd.compute('model_*')
    assert after.interfaces == before.interfaces
    assert after.matrix[0, 1] == pytest.approx(before.matrix[0, 1])
    assert cluster.compute('model_*', cutoff=99.0).matrix[0, 1] > \
        ligand + 1.0


def test_scopes():
    names = touching_models(((0, 2), (120, 2)))
    pair = irmsd.compute('model_*', scope='pair')
    shared = irmsd.compute('model_*', scope='all')
    assert pair.scope == 'pair' and shared.scope == 'all'
    # every pair over the same atoms, and at least as many
    assert shared.atoms_compared[0] == shared.atoms_compared[1]
    assert pair.atoms_compared[1] <= shared.atoms_compared[1]
    assert pair.atoms_compared[0] < shared.atoms_compared[0]
    assert members(pair) == members(shared) == [names[:2], names[2:]]
    # "all" is RMSD Clustering's measure=interface, given its settings
    theirs = cluster.compute('model_*', measure='interface', cutoff=99.0)
    ours = irmsd.compute('model_*', scope='all', atoms='ca',
                         distance=cluster.CONTACT)
    assert ours.matrix == pytest.approx(theirs.matrix, abs=1e-6)


def test_the_cutoff_decides():
    names = touching_models(((0, 3), (120, 3)))
    assert len(irmsd.compute('model_*', cutoff=0.01).clusters) == len(names)
    assert len(irmsd.compute('model_*', cutoff=500).clusters) == 1
    assert irmsd.capri_class(0.8) == 'high'
    assert irmsd.capri_class(1.5) == 'medium'
    assert irmsd.capri_class(3.9) == 'acceptable'
    assert irmsd.capri_class(4.1) == ''


def test_set_aside():
    labels, pool = cluster.set_aside([0, 0, 0, 1, 1, 2, 3], 2)
    assert labels == [0, 0, 0, 1, 1, 2, 2] and pool == 2
    assert cluster.set_aside([0, 0, 1], 1) == ([0, 0, 1], None)
    # nothing is large enough: all in the pool
    assert cluster.set_aside([0, 1, 2], 2) == ([0, 0, 0], 0)


def test_small_clusters_are_unclustered(tmp_path, capsys):
    names = touching_models(((0, 3), (120, 3), (240, 1)))
    everything = irmsd.compute('model_*')
    assert [c['size'] for c in everything.clusters] == [3, 3, 1]
    assert not any(c.get('unclustered') for c in everything.clusters)

    out = tmp_path / 'clusters.csv'
    result = cmd.cluster_irmsd('model_*', min_size=2, file=str(out))
    assert members(result) == [names[:3], names[3:6]]
    alone = result.clusters[-1]
    assert alone['unclustered'] and alone['members'] == names[6:]
    assert alone['representative'] is None and alone['best'] is None
    text = capsys.readouterr().out
    assert '1 unclustered (in clusters of fewer than 2)' in text
    assert [line.split()[0] for line in text.splitlines()
            if line.split()[1:2] in (['3'], ['1'])] == ['1', '2', '-']
    # a group of its own, in grey
    assert cmd.get_object_list('(cluster_unclustered)') == names[6:]
    assert cmd.get_names('group_objects') == [
        'cluster_1', 'cluster_2', 'cluster_unclustered']
    colors = set()
    cmd.iterate('%s and chain H and elem C' % names[6],
                'colors.add(color)', space={'colors': colors})
    assert colors == {cmd.get_color_index(cluster.UNCLUSTERED_COLOR)}
    rows = out.read_text().splitlines()
    assert rows[0].startswith('model,cluster,representative,best,')
    assert [row.split(',')[1] for row in rows[1:]] == list('1112220')
    # clustered again without a least size: the group goes
    cmd.cluster_irmsd('model_*', quiet=1)
    assert cmd.get_names('group_objects') == [
        'cluster_1', 'cluster_2', 'cluster_3']
    assert sorted(cmd.get_object_list('all')) == sorted(names)
    # nothing large enough: nothing fails, nothing is a cluster
    nothing = cmd.cluster_irmsd('model_*', min_size=5, quiet=1)
    assert [c.get('unclustered') for c in nothing.clusters] == [True]
    assert cmd.get_names('group_objects') == ['cluster_unclustered']


def test_models_that_do_not_touch_are_left_out(capsys):
    names = touching_models(((0, 3),))
    loose = names[-1]
    cmd.translate([200.0, 0, 0], '%s and chain H+L' % loose, camera=0)
    result = cmd.cluster_irmsd('model_*')
    assert result.skipped == [loose]
    assert result.models == names[:2]
    assert 'left out, their sides do not touch: %s' % loose in \
        capsys.readouterr().out
    # (it is not moved or grouped either)
    assert loose not in cmd.get_object_list('(cluster_1)')
    cmd.translate([200.0, 0, 0], '%s and chain H+L' % names[0], camera=0)
    with pytest.raises(ValueError, match='no interfaces to compare'):
        irmsd.compute('model_*')


def test_command_superposes_groups_and_colors(capsys):
    names = touching_models(((0, 3), (120, 2)))
    result = cmd.cluster_irmsd('model_*')
    text = capsys.readouterr().out
    assert 'cluster_irmsd: 5 structures, 2 clusters; i-RMSD over' in text
    assert 'backbone atoms of residues within 10 A' in text
    assert 'average, cutoff 4.0 A' in text
    assert cmd.get_names('group_objects') == ['cluster_1', 'cluster_2']
    # all on the antigen of the largest cluster's representative
    reference = cmd.get_coords('%s and chain A' %
                               result.clusters[0]['representative'])
    for name in names:
        assert cmd.get_coords('%s and chain A' % name) == \
            pytest.approx(reference, abs=1e-2)
    colors = {}
    for name in names:
        cmd.iterate('%s and chain H and name CA' % name,
                    '_ = colors.setdefault(model, color)',
                    space={'colors': colors})
    assert len({colors[n] for n in names[:3]}) == 1
    assert colors[names[0]] != colors[names[3]]


def test_one_undo_step():
    touching_models(((0, 2), (120, 2)))
    assert 'cluster_irmsd' in undo.UNDOABLE
    before = {n: cmd.get_coords(n).copy() for n in cmd.get_names()}
    cmd.undo_enable()
    try:
        steps = len(undo.stack.undo_stack)
        cmd.cluster_irmsd('model_*', quiet=1)
        assert len(undo.stack.undo_stack) == steps + 1
        cmd.undo()
        assert 'cluster_1' not in cmd.get_names('all')
        for name, xyz in before.items():
            assert cmd.get_coords(name) == pytest.approx(xyz)
    finally:
        cmd.undo_disable()


def test_errors():
    touching_models(((0, 2),))
    for wrong, message in (({'scope': 'some'}, 'scope must be'),
                           ({'cutoff': 0}, 'cutoff'),
                           ({'distance': -1}, 'distance'),
                           ({'min_size': 0}, 'min_size'),
                           ({'method': 'ward'}, 'method must be'),
                           ({'atoms': 'all'}, 'atoms must be'),
                           ({'objects': 'model_1'}, 'at least two'),
                           ({'antigen': 'chain Z'}, 'no interfaces'),
                           ({'antigen': 'chain A', 'antibody': 'chain A'},
                            'two sides')):
        with pytest.raises(ValueError, match=message):
            irmsd.compute(**dict({'objects': 'model_*'}, **wrong))


def test_matrix_as_pair_by_pair():
    rng = np.random.default_rng(1)
    n, atoms = 31, 120
    base = rng.normal(size=(atoms, 3)) * 15 + [80, -40, 120]
    xyz = base + rng.normal(size=(n, atoms, 3)) * rng.uniform(
        0.0, 3.0, size=(n, 1, 1))
    xyz[7] = xyz[3]                         # the same structure twice
    xyz[9] = xyz[3] + 1e-5 * rng.normal(size=(atoms, 3))
    masks = rng.random((n, atoms)) < 0.25
    masks[:, :5] = False                    # atoms in no interface
    slow = np.zeros((n, n))
    counts = []
    for i in range(n):
        for j in range(i + 1, n):
            use = masks[i] | masks[j]
            counts.append(int(use.sum()))
            slow[i, j] = slow[j, i] = irmsd.irmsd(xyz[i][use], xyz[j][use])
    fast, compared = irmsd.irmsd_matrix(xyz, masks, 'pair')
    assert np.allclose(fast, slow, rtol=0, atol=1e-9)
    assert compared == (min(counts), max(counts))
    assert np.array_equal(fast, fast.T) and not fast.diagonal().any()
    # (nearly) the same structure: measured, not a difference of sums
    assert fast[3, 7] < 1e-10
    assert fast[3, 9] == pytest.approx(slow[3, 9], rel=1e-6)
    assert 0 < fast[3, 9] < irmsd.EXACT_BELOW
    shared = masks.any(axis=0)
    slow = np.array([[irmsd.irmsd(xyz[i][shared], xyz[j][shared])
                      if i != j else 0.0 for j in range(n)]
                     for i in range(n)])
    fast, compared = irmsd.irmsd_matrix(xyz, masks, 'all')
    assert np.allclose(fast, slow, rtol=0, atol=1e-10)
    assert compared == (shared.sum(), shared.sum())
    assert irmsd.irmsd_matrix(xyz[:1], masks[:1])[0].tolist() == [[0.0]]
    assert irmsd.irmsd_matrix(xyz[:2], masks[:2])[0][0, 1] == \
        pytest.approx(irmsd.irmsd(xyz[0][masks[0] | masks[1]],
                                  xyz[1][masks[0] | masks[1]]))


def test_dialog(tmp_path):
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import shiboken6
    import pymol.plugins  # noqa: F401 (adds the startup plugin path)
    from pmg_tk.startup.irmsd_gui.main import InterfaceClusterDialog

    names = touching_models(((0, 3), (120, 3), (240, 1)))
    cmd.fab('ACD', 'other')
    dialog = InterfaceClusterDialog(None, _self=cmd)
    try:
        assert dialog.windowTitle() == 'Interface RMSD Clustering'
        # CAPRI's settings to start from
        assert dialog.cutoff.value() == irmsd.CUTOFF == 4.0
        assert dialog.distance.value() == irmsd.INTERFACE == 10.0
        assert dialog.atoms.currentData() == 'backbone'
        assert dialog.scope.currentData() == 'pair'
        assert dialog.min_size.value() == 1
        assert dialog.table.horizontalHeaderItem(2).text() == 'Mean i-RMSD'
        dialog.refresh()
        dialog.filter.setText('model')
        assert dialog.checked_models() == names
        dialog.run()
        assert dialog.fit.text() == 'chain A'       # (detected)
        assert dialog.rmsd.text() == 'chain H+L'
        assert dialog.table.rowCount() == 3
        assert [dialog.table.item(row, 1).text() for row in range(3)] == \
            ['3', '3', '1']
        assert '7 structures, 3 clusters' in dialog.status.text()
        assert dialog.save_button.isEnabled()
        # the smallest cluster worth the name
        dialog.min_size.setValue(2)
        dialog.run()
        assert dialog.table.rowCount() == 3
        assert dialog.table.item(2, 0).text() == '-'
        assert dialog.table.item(2, 2).text() == '-'
        assert dialog.table.item(2, 6).text() == 'unclustered'
        assert '2 clusters, 1 unclustered' in dialog.status.text()
        dialog.table.selectRow(2)
        assert cmd.get_object_list('enabled and model_*') == names[6:]
        dialog.show_best()                           # (no scores: says so)
        assert 'No scores found' in dialog.status.text()
        dialog.show_all()
        assert len(cmd.get_object_list('enabled and model_*')) == 7
        assert len(dialog.heatmap.order) == 7
        dialog.heatmap.resize(160, 160)
        assert not dialog.heatmap.grab().isNull()
        # the settings reach the clustering
        dialog.min_size.setValue(1)
        dialog.cutoff.setValue(50.0)
        dialog.distance.setValue(5.0)
        dialog.scope.setCurrentIndex(dialog.scope.findData('all'))
        dialog.atoms.setCurrentIndex(dialog.atoms.findData('ca'))
        dialog.run()
        assert dialog.table.rowCount() == 1
        assert dialog.result.scope == 'all'
        assert dialog.result.distance == 5.0 and dialog.result.atoms == 'ca'
        path = os.path.join(str(tmp_path), 'out.csv')
        dialog.result.save_csv(path)
        assert len(open(path).read().splitlines()) == 8
        # models that cannot be compared: said, not raised
        dialog.fit.setText('chain Z')
        dialog.run()
        assert dialog.status.text().startswith('Clustering failed')
    finally:
        shiboken6.delete(dialog)
