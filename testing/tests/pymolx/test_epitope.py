'''
Epitope Prediction (data/startup/epitope_gui, pymolx.epitope, parity item
L-30): antibody-antigen cofolding with AlphaFold 3 - chain roles, seeds
and samples, the antibody's MSA, and the dialog.
'''

import json
import os

import pytest

from pymol import cmd
import pymol.plugins  # noqa: F401 (puts $PYMOL_DATA/startup on the path)
from pymolx import epitope, prediction

# a heavy chain through CH1, a light chain, and something that is neither
HEAVY = ('EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYAD'
         'SVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSASTK'
         'GPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLS')
LIGHT = ('DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRF'
         'SGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQ'
         'LKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADY')
ANTIGEN = ('MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFE'
           'DHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPER')
TEXT = '>Chain H\n%s\n>Chain L\n%s\n>Chain A\n%s\n' % (HEAVY, LIGHT, ANTIGEN)


@pytest.fixture
def job():
    return epitope.job_from_text(TEXT, 'ab_ag')


def test_roles(job):
    assert {c.id: job.role(c.id) for c in job.chains} == {
        'H': 'heavy', 'L': 'light', 'A': 'antigen'}
    assert [c.id for c in job.antibody_chains()] == ['H', 'L']
    assert [c.id for c in job.antigen_chains()] == ['A']
    # a lone heavy variable domain looks like a VHH on its own ...
    alone = epitope.job_from_text('>Chain H\n%s\n>Chain A\n%s\n'
                                  % (HEAVY[:120], ANTIGEN))
    assert alone.role('H') == 'single'
    # ... but beside a light chain it is the heavy chain of an Fv
    pair = epitope.job_from_text('>H\n%s\n>L\n%s\n>A\n%s\n'
                                 % (HEAVY[:120], LIGHT, ANTIGEN))
    assert pair.role('H') == 'heavy' and pair.role('L') == 'light'
    assert epitope.detect_role(ANTIGEN) == 'antigen'
    assert epitope.role_label('heavy') == 'Antibody heavy chain'
    with pytest.raises(ValueError, match='role'):
        job.set_role('H', 'nonsense')


def test_defaults_and_checks(job):
    assert len(job.seeds) == epitope.DEFAULT_SEEDS >= 10
    assert job.samples == epitope.DEFAULT_SAMPLES
    assert job.predictions() == len(job.seeds) * job.samples
    assert job.antibody_msa is True            # the user's call to turn off
    job.check()
    assert not job.warnings('af3')
    # an antibody with nothing to bind, and an antigen with no antibody
    job.set_role('A', 'single')
    with pytest.raises(ValueError, match='no antigen'):
        job.check()
    for chain in job.chains:
        job.set_role(chain.id, 'antigen')
    with pytest.raises(ValueError, match='no antibody'):
        job.check()


def test_warnings(job):
    job.seeds = [1, 2]
    assert any('1,000 seeds' in w for w in job.warnings('af3'))
    job.seeds = list(range(1, 11))
    job.set_role('L', 'antigen')
    assert any('without a light one' in w for w in job.warnings('af3'))
    job.set_role('L', 'light')
    job.antibody_msa = False
    assert any('no MSA for the antibody' in w for w in job.warnings('af3'))
    job.samples = 0
    with pytest.raises(ValueError, match='samples per seed'):
        job.check()


def test_alphafold3_input(job):
    found = json.loads(epitope.write(job))
    assert found['dialect'] == 'alphafold3'
    assert found['modelSeeds'] == list(range(1, epitope.DEFAULT_SEEDS + 1))
    proteins = {e['protein']['id']: e['protein']
                for e in found['sequences'] if 'protein' in e}
    assert sorted(proteins) == ['A', 'H', 'L']
    # with the MSA searched, nothing is said about it
    for protein in proteins.values():
        assert 'unpairedMsa' not in protein and 'templates' not in protein
    # without it, the antibody chains get an empty one; the antigen keeps
    # the full pipeline
    job.antibody_msa = False
    proteins = {e['protein']['id']: e['protein'] for e
                in json.loads(epitope.write(job))['sequences']
                if 'protein' in e}
    for chain_id in ('H', 'L'):
        assert proteins[chain_id]['unpairedMsa'] == ''
        assert proteins[chain_id]['pairedMsa'] == ''
        assert proteins[chain_id]['templates'] == []
    assert 'unpairedMsa' not in proteins['A']
    # it is AlphaFold 3 cofolding, not the other methods
    with pytest.raises(ValueError, match='AlphaFold 3'):
        epitope.write(job, 'rf3')


def test_what_prediction_offers_still_works(job):
    # modifications, glycans, ligands and disulfides come from the
    # Structure Prediction core
    site = job.chain('A').sequons()
    if site:
        job.chain('A').glycans[site[0]] = 'G0F'
    job.disulfides = [(('H', 22), ('H', 96))]
    job.ligands = [('ATP', 1)]
    found = json.loads(epitope.write(job))
    bonds = found['bondedAtomPairs']
    assert [['H', 22, 'SG'], ['H', 96, 'SG']] in bonds
    ligands = [e['ligand'] for e in found['sequences'] if 'ligand' in e]
    assert any(l['ccdCodes'] == ['ATP'] for l in ligands)
    if site:
        assert any(len(l['ccdCodes']) > 1 for l in ligands)


def test_samples_are_a_run_flag(job, tmp_path):
    job.samples = 8
    assert epitope.run_arguments(job) == ['--num_diffusion_samples', '8']
    # ... and not in the input file, which has no such field
    assert 'num_diffusion' not in epitope.write(job)
    # a job script only passes it on if it forwards its arguments
    plain = tmp_path / 'plain.sh'
    plain.write_text('#!/bin/bash\ndocker run af3 --json_path=$1\n')
    assert epitope.script_forwards_arguments(str(plain)) is False
    forwards = tmp_path / 'forwards.sh'
    forwards.write_text('#!/bin/bash\ndocker run af3 --json_path=$1 "$@"\n')
    assert epitope.script_forwards_arguments(str(forwards)) is True
    assert epitope.script_forwards_arguments(str(tmp_path / 'gone.sh')) is None


def test_seed_clash(job, tmp_path):
    # AlphaFold 3 refuses --num_seeds together with a file of many seeds:
    # "Input must have one rng seed to set multiple seeds"
    expands = tmp_path / 'expands.sh'
    expands.write_text('#!/bin/bash\nrun_alphafold.py --json_path=$1 \\\n'
                       '  --num_seeds=10 --num_diffusion_samples=10\n')
    assert epitope.script_num_seeds(str(expands)) == 10
    assert len(job.seeds) > 1
    message = job.seed_clash(str(expands))
    assert '--num_seeds=10' in message and 'exactly one seed' in message
    # one seed in the file: the script may expand it
    job.seeds = [1]
    assert job.seed_clash(str(expands)) == ''
    # a script that leaves the seeds alone is fine either way
    plain = tmp_path / 'plain.sh'
    plain.write_text('#!/bin/bash\nrun_alphafold.py --json_path=$1\n')
    assert epitope.script_num_seeds(str(plain)) is None
    job.seeds = list(range(1, 11))
    assert job.seed_clash(str(plain)) == ''
    assert job.seed_clash('') == ''
    assert epitope.script_num_seeds(str(tmp_path / 'gone.sh')) is None


def test_save(job, tmp_path):
    path = str(tmp_path / 'epitope_af3.json')
    assert epitope.save(job, path) == path
    assert json.load(open(path))['name'] == 'ab_ag'


def test_dialog(tmp_path):
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import shiboken6
    from pmg_tk.startup.epitope_gui.main import EpitopeDialog

    dialog = EpitopeDialog(None, _self=cmd)
    try:
        assert dialog.windowTitle() == 'Epitope Prediction'
        assert dialog.seed_count.value() == epitope.DEFAULT_SEEDS
        assert dialog.samples.value() == epitope.DEFAULT_SAMPLES
        assert dialog.antibody_msa.isChecked()

        dialog.job = epitope.job_from_text(TEXT, 'ab_ag')
        dialog.refresh()
        # a row per chain, its role guessed
        assert dialog.roles.rowCount() == 3
        shown = {dialog.roles.item(r, 0).text():
                 dialog.roles.cellWidget(r, 1).currentData()
                 for r in range(3)}
        assert shown == {'H': 'heavy', 'L': 'light', 'A': 'antigen'}
        assert '10 seed(s) x 10 sample(s) = 100' in dialog.status.text()

        # changing a role in the table reaches the job
        combo = dialog.roles.cellWidget(2, 1)
        combo.setCurrentIndex(combo.findData('single'))
        assert dialog.job.role('A') == 'single'
        combo.setCurrentIndex(combo.findData('antigen'))

        # the settings reach the file
        dialog.seed_count.setValue(25)
        dialog.samples.setValue(3)
        dialog.antibody_msa.setChecked(False)
        dialog._apply()
        found = json.loads(epitope.write(dialog.job))
        assert len(found['modelSeeds']) == 25
        proteins = {e['protein']['id']: e['protein']
                    for e in found['sequences'] if 'protein' in e}
        assert proteins['H']['unpairedMsa'] == ''
        assert 'unpairedMsa' not in proteins['A']
        assert dialog.job.samples == 3
        assert not dialog.job.antibody_msa
    finally:
        dialog.close()
        shiboken6.delete(dialog)
        QtWidgets.QApplication.processEvents()


def test_prediction_dialog_still_works():
    # the shared dialog was split so this one can replace its settings box
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import shiboken6
    from pmg_tk.startup.prediction_gui.main import PredictionDialog

    dialog = PredictionDialog(None, _self=cmd)
    try:
        assert dialog.method.currentData() == 'af3'
        assert dialog.seeds.isEnabled()
        dialog.method.setCurrentIndex(dialog.method.findData('esmfold2'))
        assert not dialog.seeds.isEnabled()
        dialog.job = prediction.job_from_text('>Chain A\nACDEFGHIK\n')
        dialog.refresh()
        assert dialog.save_button.isEnabled()
    finally:
        dialog.close()
        shiboken6.delete(dialog)
        QtWidgets.QApplication.processEvents()
