'''
Superposition/Alignment plugin (data/startup/alignment_gui, parity item
L-18): command building for each method and mode, running it, and the
dialog (headless, Qt offscreen platform).
'''

import os

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
    finally:
        dialog.deleteLater()
