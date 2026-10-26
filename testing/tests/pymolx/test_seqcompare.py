'''
Sequence comparison (pymolx.seqcompare: compare_sequences) and the
Sequence Comparison dialog.
'''

import os

import pytest

from pymol import cmd
from pymolx import seqcompare, undo

DATA = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
REF = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ'


@pytest.fixture(autouse=True)
def session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


def variants():
    cmd.fab(REF, 'ref', chain='A')
    cmd.fab('MKTAYIAKQRQLSFVKSHDSRQLEERLGLIEVQ', 'var', chain='A')  # I12L, F19D
    cmd.fab('MKTAYIAKQRQISFGGGVKSHFSRQLEERLGLIEVQ', 'ins', chain='A')
    cmd.fab('MKTAYIAKQRQISFVKLEERLGLIEVQ', 'del', chain='A')
    cmd.fab('YIAKQRQISFVKSHFSRQLEERLG', 'short', chain='A')  # both ends


def aligned(result):
    return dict(zip(*(iter(result.fasta().replace('\n', ' ').split()),) * 2))


def test_alignment_with_insertions_and_deletions():
    variants()
    result = seqcompare.compute('ref var ins del short', 'ref')
    rows = {}
    name = None
    for line in result.fasta().splitlines():
        if line.startswith('>'):
            name = line[1:]
            rows[name] = ''
        else:
            rows[name] += line
    assert rows['ref/A'] == 'MKTAYIAKQRQISF---VKSHFSRQLEERLGLIEVQ'
    assert rows['ins/A'] == 'MKTAYIAKQRQISFGGGVKSHFSRQLEERLGLIEVQ'
    assert rows['del/A'].replace('-', '') == 'MKTAYIAKQRQISFVKLEERLGLIEVQ'
    assert rows['del/A'].count('-') == 9
    assert rows['short/A'] == '----YIAKQRQISF---VKSHFSRQLEERLG-----'
    stats = result.stats
    assert stats['var/A']['identity'] == pytest.approx(100 * 31 / 33)
    assert stats['var/A']['similarity'] == pytest.approx(100 * 32 / 33)
    assert stats['ins/A']['gaps'] == 3 and stats['del/A']['gaps'] == 6
    assert stats['short/A']['gaps'] == 0  # (missing ends don't count)


def test_scores():
    variants()
    result = seqcompare.compute('ref var del short', 'ref')
    var = result.scores['var/A']
    assert var[10] == seqcompare.IDENTICAL
    assert var[11] == seqcompare.SIMILAR     # I -> L
    assert var[18] == seqcompare.DIFFERENT   # F -> D
    ref = result.scores['ref/A']
    # the reference: mean over the others; a deletion inside counts as
    # different, a missing end doesn't
    assert ref[18] == pytest.approx((0 + 0 + 100) / 3)  # var, del, short
    assert ref[0] == pytest.approx((100 + 100) / 2)     # short lacks it
    assert seqcompare.pair_score('K', 'R') == seqcompare.SIMILAR
    assert seqcompare.pair_score('A', 'S') == seqcompare.SIMILAR
    assert seqcompare.pair_score('A', 'T') == seqcompare.NEUTRAL
    assert seqcompare.pair_score('W', 'G') == seqcompare.DIFFERENT
    assert seqcompare.pair_score('X', 'X') == seqcompare.DIFFERENT


def test_all_against_all():
    variants()
    result = seqcompare.compute('ref var ins', 'ref', mode='all')
    group = result.groups[0]
    for k, column in enumerate(group.columns):
        values = {round(result.scores[name][pos], 6)
                  for name, pos in column.items()
                  if pos in result.scores[name]}
        assert len(values) <= 1


def test_chains_matched_by_id():
    heavy = 'EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVA'
    light = 'DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIY'
    for name, h, l in (('ab1', heavy, light),
                       ('ab2', heavy.replace('NIKD', 'TFSS'), light)):
        cmd.fab(h, name + 'h', chain='H')
        cmd.fab(l, name + 'l', chain='L')
        cmd.create(name, '%sh or %sl' % (name, name))
        cmd.delete('%sh %sl' % (name, name))
    result = seqcompare.compute('ab1 ab2', 'ab1')
    assert [[e.name for e in g.entries] for g in result.groups] == [
        ['ab1/H', 'ab2/H'], ['ab1/L', 'ab2/L']]
    assert result.stats['ab2/L']['identity'] == 100
    assert result.stats['ab2/H']['identity'] == pytest.approx(
        100 * (len(heavy) - 4) / len(heavy))
    names, matrix = result.identity_matrix()
    assert names == ['ab1/H', 'ab2/H', 'ab1/L', 'ab2/L']
    assert matrix[('ab1/L', 'ab2/L')] == 100
    assert ('ab1/H', 'ab1/L') not in matrix
    only_h = seqcompare.compute('ab1 ab2', 'ab1', selection='chain H')
    assert len(only_h.groups) == 1


def test_residues_without_coordinates():
    cmd.load(os.path.join(DATA, '1oky.pdb.gz'), 'kinase')
    cmd.create('trimmed', 'kinase and present')  # no unmodeled residues
    result = seqcompare.compute('kinase trimmed', 'kinase')
    fasta = result.fasta()
    assert 'mdgtaaeprpgagslqhaqpp' in fasta.replace('\n', '')
    # unmodeled residues: in the alignment, without color; phosphoserine
    # (SEP) as serine
    assert result.stats['trimmed/A']['identity'] == 100
    assert seqcompare.ONE_LETTER['SEP'] == 'S'


def test_command_colors_alignment_object_and_files(tmp_path, capsys):
    variants()
    path = tmp_path / 'aln.fasta'
    result = cmd.compare_sequences('ref var', 'ref', file=str(path))
    assert 'Similarity' in capsys.readouterr().out
    assert path.read_text().startswith('>ref/A\n')
    colors = {}
    cmd.iterate('var and name CA', 'colors[int(resi)] = color',
                space={'colors': colors})
    assert colors[1] == cmd.get_color_index('cmp_100')     # identical
    assert colors[19] == cmd.get_color_index('cmp_000')    # F -> D
    assert colors[12] == cmd.get_color_index('cmp_067')    # I -> L
    raw = cmd.get_raw_alignment(result.alignment)
    assert len(raw) == 33
    csv_path = tmp_path / 's.csv'
    result.save_csv(str(csv_path))
    lines = csv_path.read_text().splitlines()
    assert lines[0] == 'chain,resi,resn,column,score' and len(lines) == 67


def test_superpose():
    variants()
    cmd.translate([30, 0, 0], 'var')
    cmd.compare_sequences('ref var', 'ref', superpose=1, quiet=1)
    assert cmd.rms_cur('var and name CA', 'ref and name CA',
                       matchmaker=4) < 0.5


def test_errors():
    variants()
    with pytest.raises(ValueError, match='at least two'):
        seqcompare.compute('ref')
    with pytest.raises(ValueError, match='mode'):
        seqcompare.compute('ref var', mode='pairs')
    with pytest.raises(ValueError, match='reference'):
        seqcompare.compute('ref var', reference='zzz')
    cmd.alter('var', 'chain="B"')
    cmd.fab('AAAA', 'two', chain='C')
    cmd.create('multi', 'ref or two')
    with pytest.raises(ValueError, match='no chains to compare'):
        seqcompare.compute('multi var', 'multi')


def test_one_undo_step():
    variants()
    cmd.undo_enable()
    try:
        steps = len(undo.stack.undo_stack)
        cmd.compare_sequences('ref var', 'ref', quiet=1)
        assert len(undo.stack.undo_stack) == steps + 1
    finally:
        cmd.undo_disable()


def test_dialog():
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import pymol.plugins  # noqa: F401 (adds the startup plugin path)
    from pmg_tk.startup.seqcompare_gui.main import SequenceDialog

    variants()
    dialog = SequenceDialog(None, _self=cmd)
    try:
        dialog.refresh()
        assert dialog.reference.currentText() == 'ref'
        dialog.run()
        text = dialog.view.toPlainText()
        assert 'ref/A (ref)' in text and 'MKTAYIAKQRQISF' in text
        assert dialog.matrix.rowCount() == 5
        assert dialog.matrix.item(0, 0).text() == '100.0'
        assert dialog.summary.item(1, 3).text().endswith('%')
        assert dialog.fasta_button.isEnabled()
        dialog.view.resize(400, 300)
        assert not dialog.view.grab().isNull()
    finally:
        dialog.deleteLater()
