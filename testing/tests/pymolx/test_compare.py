'''
Structure comparison (pymolx.compare: compare_structures) and the
Structure Comparison dialog, on copies of 1rx1 (DHFR): one with a moved
loop, one with its C-terminal part turned, all moved as a whole.
'''

import math
import os

import numpy as np
import pytest

from pymol import cmd
from pymolx import compare, undo

DATA = os.path.join(os.path.dirname(__file__), '..', '..', 'data')


@pytest.fixture(autouse=True)
def session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


def structures():
    for name in ('wt', 'loopmut', 'hinge'):
        cmd.load(os.path.join(DATA, '1rx1.pdb'), name)
    cmd.remove('hydro or solvent or not polymer')
    cmd.alter_state(1, 'loopmut and resi 16-22', 'y = y + 4')
    center = cmd.get_coords('hinge and resi 110 and name CA')[0]
    t = math.radians(25)
    turn = np.array([[math.cos(t), -math.sin(t), 0],
                     [math.sin(t), math.cos(t), 0], [0, 0, 1]])
    xyz = cmd.get_coords('hinge and resi 110-159')
    cmd.load_coords((xyz - center) @ turn.T + center,
                    'hinge and resi 110-159')
    cmd.translate([20, 5, 0], 'loopmut')
    cmd.turn('y', 40)  # (camera only)
    cmd.rotate('x', 70, 'hinge', camera=0)


def mean(result, obj, selection):
    indices = []
    cmd.iterate('%s and name CA and (%s)' % (obj, selection),
                'indices.append(index)', space={'indices': indices})
    values = [result.scores[obj][i] for i in indices
              if i in result.scores[obj]]
    return sum(values) / len(values)


def test_lddt_shows_local_changes():
    structures()
    result = compare.compute('wt loopmut hinge', 'wt')
    assert len(result.columns) == 159
    assert mean(result, 'loopmut', 'resi 16-22') < 70
    assert mean(result, 'loopmut', 'resi 40-100') > 95
    assert mean(result, 'hinge', 'resi 40-100') > 95
    # the reference: the mean of the others
    assert mean(result, 'wt', 'resi 16-22') < mean(result, 'wt',
                                                   'resi 40-100')
    stats = result.stats
    assert stats['loopmut']['identity'] == 100
    assert stats['loopmut']['aligned'] == 159
    assert 0.5 < stats['loopmut']['rmsd'] < 1.5
    assert stats['hinge']['rmsd'] > 2
    # superimposed on the reference's core
    assert cmd.rms_cur('loopmut and resi 40-100 and name CA',
                       'wt and resi 40-100 and name CA',
                       matchmaker=4) < 0.01


def test_distance_after_superposition():
    structures()
    result = compare.compute('wt loopmut hinge', 'wt', metric='distance')
    assert mean(result, 'loopmut', 'resi 16-22') == pytest.approx(4, abs=0.3)
    assert mean(result, 'loopmut', 'resi 40-100') < 0.1
    assert mean(result, 'hinge', 'resi 120-159') > 4


def test_lddt_function():
    xyz = np.random.default_rng(0).normal(0, 6, (30, 3))
    assert compare.lddt(xyz, xyz, xyz) == pytest.approx(100)
    moved = xyz.copy()
    moved[0] += [3, 0, 0]
    scores = compare.lddt(xyz, moved, xyz)
    assert scores[0] < 60 and scores[1:].min() > scores[0]
    # neighbors the model lacks count as lost
    partial = compare.lddt(xyz[:20], xyz[:20], xyz)
    assert partial.max() <= 100 and partial.mean() < 100


def test_all_against_all_and_partners_missing():
    structures()
    cmd.remove('loopmut and resi 140-159')
    cmd.remove('hinge and resi 140-159')
    result = compare.compute('wt loopmut hinge', 'wt', mode='all')
    # one score per alignment position, the same in every structure
    for column in result.columns:
        values = {round(result.scores[o][i], 6)
                  for o, i in column.items() if i in result.scores[o]}
        assert len(values) <= 1
    # the reference's residues without partners: no score
    indices = []
    cmd.iterate('wt and resi 145 and name CA', 'indices.append(index)',
                space={'indices': indices})
    assert indices[0] not in result.scores['wt']


def test_command_colors_legend_and_csv(tmp_path, capsys):
    structures()
    path = str(tmp_path / 'sim.csv')
    result = cmd.compare_structures('wt loopmut hinge', 'wt', file=path)
    assert 'Mean lDDT' in capsys.readouterr().out
    blue = cmd.get_color_index('cmp_100')
    colors = {}
    cmd.iterate('loopmut and name CA', 'colors[int(resi)] = color',
                space={'colors': colors})
    assert colors[60] == blue
    assert colors[19] != blue
    assert cmd.get_type(result.legend) == 'object:ramp'
    lines = open(path).read().splitlines()
    assert lines[0] == 'object,chain,resi,resn,column,lddt'
    assert len(lines) == 1 + 3 * 159
    # the ramp: red, yellow, green, blue
    assert compare.ramp_color(0) == (1, 0, 0)
    assert compare.ramp_color(1) == (0, 0, 1)
    assert compare.ramp_color(1 / 3) == pytest.approx((1, 1, 0))


def test_select_dissimilar_and_errors():
    structures()
    result = compare.compute('wt loopmut', 'wt')
    count = compare.select_dissimilar(result, 70)
    assert count > 0
    assert cmd.count_atoms('dissimilar and loopmut and resi 16-22') > 0
    assert cmd.count_atoms('dissimilar and loopmut and resi 60-80') == 0
    with pytest.raises(ValueError, match='at least two'):
        compare.compute('wt')
    with pytest.raises(ValueError, match='metric'):
        compare.compute('wt loopmut', metric='tm')
    with pytest.raises(ValueError, match='reference'):
        compare.compute('wt loopmut', reference='other')


@pytest.mark.parametrize('method', ['super', 'cealign'])
def test_other_alignment_methods(method):
    structures()
    result = compare.compute('wt loopmut', 'wt', method=method)
    assert len(result.columns) > 140
    assert mean(result, 'loopmut', 'resi 40-100') > 90


def test_one_undo_step():
    structures()
    before = cmd.get_coords('loopmut').copy()
    cmd.undo_enable()
    try:
        steps = len(undo.stack.undo_stack)
        cmd.compare_structures('wt loopmut', 'wt', quiet=1)
        assert len(undo.stack.undo_stack) == steps + 1
        cmd.undo()
        assert cmd.get_coords('loopmut') == pytest.approx(before)
    finally:
        cmd.undo_disable()


def test_dialog(tmp_path):
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import pymol.plugins  # noqa: F401 (adds the startup plugin path)
    from pmg_tk.startup.compare_gui.main import CompareDialog

    structures()
    dialog = CompareDialog(None, _self=cmd)
    try:
        dialog.refresh()
        assert dialog.checked_models() == ['wt', 'loopmut', 'hinge']
        assert dialog.reference.currentText() == 'wt'
        dialog.reference.setCurrentText('loopmut')
        dialog.run()
        assert dialog.result.reference == 'loopmut'
        assert dialog.table.rowCount() == 3
        assert dialog.table.item(0, 0).text() == 'loopmut (reference)'
        assert dialog.table.horizontalHeaderItem(4).text() == 'Mean lDDT'
        assert dialog.select_button.isEnabled()
        dialog.select_dissimilar()
        assert cmd.count_atoms('dissimilar') > 0
        # the distance: scale and table follow
        dialog.metric_box.setCurrentIndex(dialog.metric_box.findData('distance'))
        assert dialog.red.value() == 3 and dialog.blue.value() == 0
        dialog.run()
        assert dialog.table.horizontalHeaderItem(4).text() == \
            'Mean distance (Å)'
        dialog.profile.resize(300, 160)
        assert not dialog.profile.grab().isNull()
    finally:
        dialog.deleteLater()
