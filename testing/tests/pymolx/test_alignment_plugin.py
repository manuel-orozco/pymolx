'''
Superposition/Alignment plugin (data/startup/alignment_gui, parity item
L-18): command building for each method and mode, running it, and the
dialog (headless, Qt offscreen platform).
'''

import os

import numpy as np
import pytest

from pymol import cmd
import pymol.plugins  # puts $PYMOL_DATA/startup on the plugin path
from pmg_tk.startup.alignment_gui import command

DEMO = cmd.exp_path('$PYMOL_DATA/demo/il2.pdb')


@pytest.fixture(autouse=True)
def clean_session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


def build(method='align', **kwargs):
    return command.format_command(*command.build_command(
        method, kwargs.pop('mobile', '*'), kwargs.pop('target', 'm1'),
        _self=cmd, **kwargs))


def test_methods():
    assert command.available_methods(cmd) == [
        'align', 'super', 'cealign', 'usalign', 'fit']


def test_many_to_one_matches_incentive():
    # text of the Incentive PyMOL dialog for the same input
    assert build(target='fold_lilrb3_d1d4_his_4cmut_model_0') == (
        'extra_fit *, fold_lilrb3_d1d4_his_4cmut_model_0, \\\n'
        '    method=align, \\\n'
        '    cycles=5, \\\n'
        '    cutoff=2.0, \\\n'
        '    mobile_state=-1, \\\n'
        '    target_state=-1')


def test_one_line():
    name, args, options = command.build_command(
        'super', 'm2', 'm1', many_to_one=False, _self=cmd)
    assert command.format_command(name, args, options, multiline=False) == (
        'super m2, m1, cycles=5, cutoff=2.0, mobile_state=-1, '
        'target_state=-1')


def test_options_per_method():
    # no outlier rejection for cealign/usalign; cealign takes (target, mobile)
    assert build('cealign', mobile='m2', many_to_one=False) == (
        'cealign mobile=m2, \\\n'
        '    target=m1, \\\n'
        '    mobile_state=-1, \\\n'
        '    target_state=-1')
    assert 'cycles' not in build('usalign')
    # outlier rejection off: one pass
    assert 'cycles=0' in build('align', outlier_rejection=False)
    assert 'cutoff' not in build('align', outlier_rejection=False)
    # alignment object
    assert 'object=aln' in build('align', object_name=' aln ')
    assert 'object' not in build('align', object_name='')


def test_parse_residues():
    assert command.parse_residues('10-120, 150-160') == '10-120+150-160'
    assert command.parse_residues('42') == '42'
    assert command.parse_residues('10:20') == '10-20'       # (also a range)
    assert command.parse_residues('5 7 9') == '5+7+9'
    assert command.parse_residues('27A-30') == '27A-30'     # insertion code
    for empty in ('', '   ', 'all', '*'):
        assert command.parse_residues(empty) == ''
    for bad in ('1 to 120', '10..20', '-5', 'x'):
        with pytest.raises(ValueError, match='numbers and ranges'):
            command.parse_residues(bad)


def test_parse_chains():
    assert command.parse_chains('A') == 'A'
    assert command.parse_chains('H, L') == 'H+L'
    assert command.parse_chains('H+L') == 'H+L'
    assert command.parse_chains('H L') == 'H+L'
    for empty in ('', 'all', '*'):
        assert command.parse_chains(empty) == ''
    for bad in ('H-L', 'chain H', 'ABCDE'):
        with pytest.raises(ValueError, match='IDs such as'):
            command.parse_chains(bad)


def test_region_selection():
    assert command.region_selection('m1') == 'm1'
    assert command.region_selection('m1', 'H', '1-120') == \
        '(m1) and chain H and resi 1-120'
    assert command.region_selection('*', '', '10-20,30') == \
        '(*) and resi 10-20+30'
    assert command.region_selection('m1', 'all', 'all') == 'm1'


def test_region_in_command():
    # the part fitted on; the whole object still moves
    text = build(many_to_one=False, mobile='m2', target='m1',
                 mobile_chains='H', mobile_residues='1-120',
                 target_chains='A', target_residues='5-124')
    assert text.startswith('align (m2) and chain H and resi 1-120, '
                           '(m1) and chain A and resi 5-124')
    # the other methods take them too
    assert 'mobile=(m2) and chain H' in build(
        'cealign', many_to_one=False, mobile='m2', mobile_chains='H')
    assert build(mobile='*', mobile_chains='H').startswith(
        'extra_fit (*) and chain H, m1')
    # nothing given: unchanged
    assert build(many_to_one=False).startswith('align *, m1')
    # a selection with a comma stays one argument
    assert build(many_to_one=False, mobile='a, b',
                 mobile_chains='H').startswith('align ((a, b) and chain H),')
    assert build(many_to_one=False, mobile='a, b').startswith('align (a, b),')


def test_quote_selection():
    assert command.quote_selection('') == 'all'
    assert command.quote_selection('resi 1, 2') == '(resi 1, 2)'
    assert command.quote_selection('(a, b)') == '(a, b)'
    assert command.quote_selection('m1 and name CA') == 'm1 and name CA'


@pytest.mark.parametrize('method', ['align', 'super', 'cealign', 'usalign'])
@pytest.mark.parametrize('many_to_one', [True, False])
def test_run(method, many_to_one):
    cmd.load(DEMO, 'm1')
    cmd.copy('m2', 'm1')
    cmd.rotate('x', 60, 'm2')
    cmd.translate([10, 0, 0], 'm2')
    mobile = '*' if many_to_one else 'm2'
    text = command.format_command(*command.build_command(
        method, mobile, 'm1', many_to_one=many_to_one, _self=cmd),
        multiline=False)
    cmd.do(text)
    cmd.sync()
    assert cmd.rms_cur('m2 and name CA', 'm1 and name CA') == \
        pytest.approx(0, abs=0.01)


# dialog

try:
    from pymol.Qt import QtWidgets
except ImportError:
    QtWidgets = None


@pytest.mark.skipif(QtWidgets is None, reason='no Qt')
def test_dialog():
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    from pmg_tk.startup.alignment_gui.main import AlignmentDialog

    cmd.load(DEMO, 'm1')
    cmd.copy('m2', 'm1')
    dialog = AlignmentDialog(None, _self=cmd)
    try:
        assert dialog.mobile.currentText() == '*'
        assert dialog.target.currentText() == 'm1'
        assert dialog.preview.toPlainText().startswith('extra_fit *, m1')

        dialog.method.setCurrentText('cealign')
        assert not dialog.outlier_rejection.isEnabled()
        assert 'cycles' not in dialog.preview.toPlainText()

        dialog.method.setCurrentText('align')
        dialog.outlier_rejection.setChecked(False)
        assert not dialog.outlier_frame.isEnabled()
        assert 'cycles=0' in dialog.preview.toPlainText()

        dialog.one_to_one.setChecked(True)
        dialog.mobile.setEditText('m2')
        assert dialog.preview.toPlainText().startswith('align m2, m1')

        cmd.rotate('x', 60, 'm2')
        dialog.outlier_rejection.setChecked(True)
        dialog.run()
        cmd.sync()
        assert cmd.rms_cur('m2 and name CA', 'm1 and name CA') == \
            pytest.approx(0, abs=0.01)

        # a stretch of residues: fitted on those, the whole object moves
        cmd.rotate('y', 45, 'm2')
        cmd.translate([12, 4, -6], 'm2')
        away = cmd.get_coords('m2 and resi 60-70')
        dialog.mobile_residues.setText('10-40')
        dialog.target_residues.setText('10-40')
        assert dialog.preview.toPlainText().startswith(
            'align (m2) and resi 10-40, (m1) and resi 10-40')
        dialog.run()
        cmd.sync()
        assert cmd.rms_cur('m2 and resi 10-40 and name CA',
                           'm1 and resi 10-40 and name CA') == \
            pytest.approx(0, abs=0.01)
        # the residues outside it came along (the object is not cut up)
        assert np.linalg.norm(
            cmd.get_coords('m2 and resi 60-70') - away, axis=1).mean() > 1.0
        assert cmd.count_atoms('m2') == cmd.count_atoms('m1')

        # what it can't read is said, and OK waits for it to be fixed
        dialog.mobile_residues.setText('10 to 40')
        assert 'numbers and ranges' in dialog.preview.toPlainText()
        assert not dialog.ok_button.isEnabled()
        dialog.mobile_chains.setText('H-L')
        dialog.mobile_residues.setText('')
        assert 'IDs such as' in dialog.preview.toPlainText()
        dialog.mobile_chains.setText('')
        assert dialog.ok_button.isEnabled()
        # cleared: the mobile side is the whole object again
        assert dialog.preview.toPlainText().startswith(
            'align m2, (m1) and resi 10-40')
    finally:
        dialog.deleteLater()
