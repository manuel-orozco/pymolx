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


def _alignment_module():
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import pymol.plugins  # noqa: F401 (adds the startup plugin path)
    from pmg_tk.startup.seqcompare_gui import alignment
    return alignment


def _long_chains():
    '''
    Two chains of 120 residues, numbered from 101, the second with two
    substitutions and without residues 40-44 of the first
    '''
    ref = (REF * 4)[:120]
    other = ref[:9] + 'W' + ref[10:39] + ref[44:89] + 'W' + ref[90:]
    cmd.fab(ref, 'long', chain='A', resi=101)
    cmd.fab(other, 'other', chain='A', resi=101)
    return seqcompare.compute('long other', 'long')


def test_alignment_lines():
    alignment = _alignment_module()
    result = _long_chains()
    group = result.groups[0]
    assert len(group.columns) == 120
    # extended: residue numbers, then a line per chain
    lines = alignment.lines(result)
    assert [(e.name if e else None, a, b) for _, e, a, b in lines] == [
        (None, 0, 120), ('long/A', 0, 120), ('other/A', 0, 120)]
    # wrapped: blocks of 50 columns with an empty line between them
    lines = alignment.lines(result, 50)
    assert len(lines) == 3 * 3 + 2 and lines[3] is None
    assert lines[4][2:] == (50, 100) and lines[-1][2:] == (100, 120)
    assert alignment.ruler(group, 0, 30) == [(0, '101'), (10, '111'),
                                             (20, '121')]
    assert alignment.ruler(group, 50, 70) == [(50, '151'), (60, '161')]
    assert alignment.cell(group, group.entries[1], 41) == ('-', None)
    assert alignment.cell(group, group.entries[1], 44) == (
        group.entries[0].sequence[44], 39)
    assert alignment.fit_block(67.9) == 60
    assert alignment.fit_block(3) == alignment.MIN_BLOCK
    text = alignment.text(result, 50).splitlines()
    assert text[0] == ' ' * 13 + '101       111       121       131' \
        '       141'
    assert text[1] == 'long/A (ref) ' + group.entries[0].sequence[:50]
    assert text[2].startswith('other/A      ') and '-----' in text[2]
    assert text[3] == '' and len(text) == 11
    # on one line: the same residues, the numbers every tenth column
    whole = alignment.text(result).splitlines()
    assert len(whole) == 3
    assert whole[1].split()[-1] == group.entries[0].sequence
    assert whole[0].split() == [str(n) for n in range(101, 221, 10)]
    # a number at a block's end is written out in full
    short = alignment.text(result, 10).splitlines()
    assert short[-3].split() == ['211']


def test_alignment_names():
    alignment = _alignment_module()
    names = ['fold_complex_seed_%d_sample_0_model' % n for n in (17, 23)]
    # short enough: as they are
    assert alignment.shorten(names, 40) == {n: n for n in names}
    # first without what they all start with
    assert list(alignment.shorten(names, 20).values()) == [
        '…d_17_sample_0_model', '…d_23_sample_0_model']
    # then without what they all end with: the numbers stay
    assert list(alignment.shorten(names, 8).values()) == [
        '…17_sam…', '…23_sam…']
    assert list(alignment.shorten(names, 4).values()) == ['…17…', '…23…']
    # nothing in common to leave out: the two ends
    assert alignment.shorten(['abcdefghijkl', 'zyx'], 7) == {
        'abcdefghijkl': 'abc…jkl', 'zyx': 'zyx'}
    assert alignment.shorten(['abcdefghijkl'], 6) == {
        'abcdefghijkl': 'abc…kl'}
    assert all(len(v) <= 5 for v in alignment.shorten(
        names + ['other_name_altogether'], 5).values())


def _sized(view, width, height):
    # (a widget that is never shown lays itself out when it is grabbed)
    view.resize(width, height)
    assert not view.grab().isNull()


def test_alignment_numbers_follow_the_letter_size():
    alignment = _alignment_module()
    import shiboken6
    plan = alignment.ruler_plan
    # up to the normal size: every tenth column, in the letters' size
    assert [plan(size, 3, 11) for size in (6, 9, 11)] == [
        (10, 6), (10, 9), (10, 11)]
    # larger letters: every fifth, every second, then every column, for
    # which the numbers get smaller
    assert [plan(size, 3, 11)[0] for size in alignment.FONT_SIZES
            if size > 11] == [5, 5, 5, 5, 2, 2, 1, 1, 1]
    assert plan(18, 3, 11) == (5, 11) and plan(21, 3, 11) == (2, 11)
    assert plan(28, 3, 11) == (1, 7) and plan(36, 3, 11) == (1, 10)
    # short numbers get there sooner, long ones later
    assert [plan(size, 2, 11)[0] for size in (12, 16, 21)] == [5, 2, 1]
    assert [plan(size, 4, 11)[0] for size in (12, 24, 28, 36)] == [
        5, 5, 2, 1]
    # never closer than the numbers are wide
    for digits in (1, 2, 3, 4, 5, 6):
        for size in alignment.FONT_SIZES:
            step, number = plan(size, digits, 11)
            assert number <= min(size, 11)
            if size > 11 and step < 10:
                assert (digits + alignment.NUMBER_GAP) * number <= \
                    step * size

    result = _long_chains()   # (numbered 101 to 220)
    group = result.groups[0]
    assert alignment.ruler(group, 0, 12, 5) == [(0, '101'), (5, '106'),
                                                (10, '111')]
    assert [n for _, n in alignment.ruler(group, 20, 23, 1)] == [
        '121', '122', '123']
    view = alignment.AlignmentView(layout=alignment.EXTENDED)
    try:
        _sized(view, 700, 200)
        view.set_result(result)
        assert view.digits == 3
        assert (view.ruler_step, view.number_size) == (10, 11)
        steps = []
        while view.can_zoom_in():
            view.zoom_in()
            steps.append(view.ruler_step)
            assert not view.grab().isNull()
            # the numbers fit between each other as drawn
            assert 3 * view.number_width < view.ruler_step * view.cell_width
        assert steps == [5, 5, 5, 5, 2, 2, 1, 1, 1]
        assert view.number_size < view.font_size
        view.set_font_size(11)
        assert (view.ruler_step, view.number_size) == (10, 11)
        # the text keeps a number every tenth column
        view.set_font_size(36)
        assert view.text().splitlines()[0].split()[:2] == ['101', '111']
    finally:
        shiboken6.delete(view)


def test_alignment_view():
    alignment = _alignment_module()
    from pymol.Qt import QtCore, QtGui, QtWidgets
    import shiboken6
    result = _long_chains()
    view = alignment.AlignmentView()
    row = alignment.AlignmentScroll(view)
    try:
        _sized(view, 700, 300)
        assert not view.lines and not view.text()
        assert row.slider.maximum() == 1  # (drawn from end to end)
        view.set_result(result)
        # wrapped: as many tens of columns as fit, nothing to move along
        fit = (view.viewport().width() - view.name_width - alignment.PAD) \
            / view.cell_width
        assert view.block == int(fit) // 10 * 10 and 20 <= view.block < 120
        assert view.scroll_range() == 0 and view.offset == 0
        assert row.slider.maximum() == 1
        assert row.slider.pageStep() == view.page_width()
        assert view.text() == alignment.text(result, view.block)
        blocks = view.block
        _sized(view, 420, 300)
        assert view.block < blocks
        _sized(view, 420, 100)
        assert view.verticalScrollBar().maximum() > 0
        _sized(view, 700, 300)
        assert view.block == blocks

        # extended: one line per chain, moved along with the slider
        changes = []
        view.view_changed.connect(lambda: changes.append(1))
        row.layout_box.setCurrentIndex(
            row.layout_box.findData(alignment.EXTENDED))
        assert view.extended and view.block is None and changes
        assert len(view.lines) == 3
        # (the last residue number, 211 at column 110, ends within it)
        expected = view.name_width + 120 * view.cell_width + \
            alignment.PAD - view.viewport().width()
        assert 0 <= view.scroll_range() - expected < 1  # (whole pixels)
        assert row.slider.maximum() == view.scroll_range()
        row.slider.setValue(row.slider.maximum())
        assert view.offset == view.scroll_range()
        first, end = view.visible_columns(0, 120)
        assert first > 0 and end == 120
        view.set_offset(10 ** 6)
        assert view.offset == view.scroll_range()
        view.set_offset(-5)
        assert view.offset == 0 and row.slider.value() == 0

        # what is under the pointer
        cw, lh = view.cell_width, view.line_height
        point = QtCore.QPoint(int(view.name_width + 9.5 * cw),
                              int(alignment.PAD + 2.5 * lh))
        group, entry, column = view.hit(point)
        assert (entry.name, column) == ('other/A', 9)
        assert view.describe(point) == 'other/A TRP 110'
        gap = QtCore.QPoint(int(view.name_width + 41.5 * cw), point.y())
        assert view.hit(gap)[2] == 41 and view.describe(gap) == ''
        name = QtCore.QPoint(8, point.y())
        assert view.hit(name)[2] is None and view.describe(name) == ''
        assert view.hit(QtCore.QPoint(point.x(), 4 * lh + 20)) is None

        # larger letters: the column in the middle stays there
        view.set_offset(view.scroll_range() // 2)
        middle = (view.name_width + view.viewport().width()) / 2.0
        column = (middle - view.name_width + view.offset) / cw
        row.larger.click()
        assert view.font_size == 12 and view.cell_width > cw
        middle = (view.name_width + view.viewport().width()) / 2.0
        after = (middle - view.name_width + view.offset) / view.cell_width
        assert abs(after - column) < 0.5
        assert row.slider.maximum() == view.scroll_range()
        while view.can_zoom_in():
            view.zoom_in()
        assert view.font_size == alignment.FONT_SIZES[-1]
        assert not row.larger.isEnabled() and row.smaller.isEnabled()
        while view.can_zoom_out():
            row.smaller.click()
        assert view.font_size == alignment.FONT_SIZES[0]
        assert not row.smaller.isEnabled()
        view.set_font_size(11)

        # the mouse wheel: along the chains (no lines to scroll), and
        # with Ctrl the letter size
        def wheel(delta, modifiers=QtCore.Qt.KeyboardModifier.NoModifier):
            where = QtCore.QPointF(view.name_width + 100, 30)
            QtWidgets.QApplication.sendEvent(
                view.viewport(), QtGui.QWheelEvent(
                    where, QtCore.QPointF(view.viewport().mapToGlobal(
                        where.toPoint())),
                    QtCore.QPoint(0, 0), QtCore.QPoint(0, delta),
                    QtCore.Qt.MouseButton.NoButton, modifiers,
                    QtCore.Qt.ScrollPhase.NoScrollPhase, False))
        view.set_offset(0)
        wheel(-120)
        assert view.offset == view.scroll_step() == row.slider.value()
        wheel(120, QtCore.Qt.KeyboardModifier.ControlModifier)
        assert view.font_size == 12
        wheel(-120, QtCore.Qt.KeyboardModifier.ControlModifier)
        assert view.font_size == 11

        # back to wrapped at the columns that were in view (in a view
        # too low for all the blocks)
        _sized(view, 420, 100)
        view.set_offset(70 * view.cell_width)
        assert view.place()[1] == 70
        view.set_layout(alignment.WRAPPED)
        assert row.layout_box.currentData() == alignment.WRAPPED
        assert view.offset == 0 and view.scroll_range() == 0
        assert view.verticalScrollBar().value() > 0
        group, column, below = view.place()
        assert column <= 70 < column + view.block and below == 0
        view.set_layout(alignment.EXTENDED)
        assert view.place()[1] == column
        _sized(view, 700, 300)

        # a name too long for its share of the view is shortened
        long_name = 'a_structure_with_a_name_that_goes_on_and_on_' * 2
        cmd.set_name('other', long_name)
        view.set_result(seqcompare.compute('long ' + long_name, 'long'))
        shown = view.labels[long_name + '/A']
        assert shown.startswith('a_structure') and '…' in shown
        assert shown.endswith('_on_/A') and len(shown) < 40
        assert view.labels['long/A'] == 'long/A (ref)'
        assert view.name_width <= alignment.NAME_SHARE * \
            view.viewport().width() + view.cell_width
        assert view.describe(QtCore.QPoint(8, point.y())) == long_name + '/A'
        assert not view.grab().isNull()
    finally:
        shiboken6.delete(row)
        shiboken6.delete(view)


def test_kinds_of_residues():
    variants()
    assert [seqcompare.kind(s) for s in (
        100.0, seqcompare.SIMILAR, 80.0, seqcompare.NEUTRAL, 50.0, 10.0,
        0.0, None)] == ['identical', 'similar', 'similar', 'neutral',
                        'neutral', 'different', 'different', 'none']
    result = seqcompare.compute('ref var del short', 'ref')
    # a chain compared to the reference: by its own pair
    assert result.kind('var/A', 10) == 'identical'
    assert result.kind('var/A', 11) == 'similar'       # I -> L
    assert result.kind('var/A', 18) == 'different'     # F -> D
    assert result.lowest['var/A'] == result.scores['var/A']
    # the reference: by the least alike of the others, whatever the mean
    assert result.scores['ref/A'][18] == pytest.approx(100 / 3)
    assert result.kind('ref/A', 18) == 'different'     # var, del differ
    assert result.kind('ref/A', 11) == 'similar'
    assert result.kind('ref/A', 0) == 'identical'      # (short lacks it)
    assert result.kind('ref/A', 5) == 'identical'
    assert result.kind('nothing/A', 0) == 'none'
    # a residue the reference has nothing for: no partner
    with_insert = seqcompare.compute('ref ins', 'ref')
    assert with_insert.kind('ins/A', 14) == 'none'     # (the first G)
    assert with_insert.kind('ins/A', 13) == 'identical'
    # all against all: every residue of a column the same kind
    everyone = seqcompare.compute('ref var ins', 'ref', mode='all')
    assert {everyone.kind(e.name, 11) for e in everyone.entries} == {
        'similar'}
    assert everyone.scores['ref/A'][11] == pytest.approx(
        (2 * seqcompare.SIMILAR + 100) / 3)
    assert {everyone.kind(e.name, 3) for e in everyone.entries} == {
        'identical'}


def test_alignment_filter():
    alignment = _alignment_module()
    from pymol.Qt import QtGui
    from pymolx.gui import theme
    import shiboken6
    from pmg_tk.startup.seqcompare_gui import main
    variants()
    dialog = main.SequenceDialog(None, _self=cmd)
    try:
        dialog.refresh()
        dialog.run()
        view, legend = dialog.view, dialog.legend
        view.set_layout(alignment.EXTENDED)
        view.resize(700, 300)
        assert set(legend.boxes) == set(seqcompare.KINDS) == view.shown
        entries = {e.name: e for e in dialog.result.entries}

        def pixel(line, column):
            # (a corner of the residue's cell, clear of its letter)
            image = view.grab().toImage()
            return image.pixelColor(
                int(view.name_width + column * view.cell_width + 1),
                int(alignment.PAD + line * view.line_height + 1)).name()

        base = QtGui.QColor(theme.COLORS['base']).name()
        blue = alignment.background(seqcompare.IDENTICAL)
        red = alignment.background(seqcompare.DIFFERENT)
        # lines: numbers, ref, var, ins, del, short; var's F19D at column
        # 21 (after the three residues ins has more)
        assert pixel(2, 21) == red and pixel(2, 5) == blue
        assert pixel(1, 5) == blue and pixel(1, 21) not in (base, blue)
        legend.boxes['different'].click()
        assert view.shown == set(seqcompare.KINDS) - {'different'}
        assert not view.is_shown(entries['var/A'], 18)
        assert view.is_shown(entries['var/A'], 5)
        assert pixel(2, 21) == base and pixel(2, 5) == blue
        # the reference there too: one of the others differs
        assert not view.is_shown(entries['ref/A'], 18)
        assert pixel(1, 21) == base and pixel(1, 5) == blue
        # only what differs
        legend.none_button.click()
        assert view.shown == set() and pixel(2, 5) == base
        legend.boxes['different'].click()
        assert legend.shown() == view.shown == {'different'}
        assert pixel(2, 21) == red and pixel(2, 5) == base
        assert pixel(1, 21) not in (base, blue)
        # the residues ins has more: no partner, and that kind alone
        grey = QtGui.QColor(theme.COLORS['border']).name()
        assert pixel(3, 15) == base
        legend.boxes['none'].click()
        assert pixel(3, 15) == grey
        legend.all_button.click()
        assert view.shown == set(seqcompare.KINDS)
        assert all(box.isChecked() for box in legend.boxes.values())
        assert pixel(2, 21) == red and pixel(2, 5) == blue
        # the text is the whole alignment, whatever is shown
        text = view.text()
        legend.none_button.click()
        assert view.text() == text
        # a new comparison keeps the filter
        dialog.run()
        assert dialog.view.shown == set()
    finally:
        shiboken6.delete(dialog)


def test_dialog():
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    alignment = _alignment_module()
    import shiboken6
    from pmg_tk.startup.seqcompare_gui import main
    from pymolx.gui.sequence_viewer import MemoryPrefs

    variants()
    prefs = MemoryPrefs()
    dialog = main.SequenceDialog(None, _self=cmd, prefs=prefs)
    try:
        dialog.refresh()
        assert dialog.reference.currentText() == 'ref'
        assert not dialog.copy_button.isEnabled()
        dialog.run()
        # (never shown: the view has no size of its own yet)
        dialog.view.resize(600, 300)
        assert not dialog.view.grab().isNull()
        text = dialog.view.text()
        assert 'ref/A (ref)' in text and 'MKTAYIAKQRQISF' in text
        assert dialog.matrix.rowCount() == 5
        assert dialog.matrix.item(0, 0).text() == '100.0'
        assert dialog.summary.item(1, 3).text().endswith('%')
        assert dialog.fasta_button.isEnabled()
        # the layout is remembered; Copy Alignment gives the text shown
        assert dialog.view.layout_name == alignment.WRAPPED
        dialog.view_scroll.layout_box.setCurrentIndex(1)
        assert dialog.view.extended
        assert prefs[main.PREF_LAYOUT] == alignment.EXTENDED
        dialog.copy_button.click()
        copied = QtWidgets.QApplication.clipboard().text()
        assert copied == dialog.view.text() and len(copied.splitlines()) == 6
        again = main.SequenceDialog(None, _self=cmd, prefs=prefs)
        try:
            assert again.view.extended
            assert again.view_scroll.layout_box.currentData() == \
                alignment.EXTENDED
        finally:
            shiboken6.delete(again)
    finally:
        shiboken6.delete(dialog)
