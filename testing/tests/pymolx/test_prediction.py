'''
Structure prediction input (pymolx.prediction) and the Structure
Prediction dialog: a FASTA file's sequences with modified residues,
glycans, ligands and disulfides, written for AlphaFold 3, RoseTTAFold3
or ESMFold2.
'''

import json
import os

import pytest

from pymol import cmd
from pymolx import prediction

HEAVY = ('EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYAD'
         'SVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSNSTAE')
LIGHT = ('DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRF'
         'SGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK')


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / 'ab.fasta'
    path.write_text('>heavy chain\n%s\n>light chain\n%s\n' % (
        '\n'.join(HEAVY[k:k + 60] for k in range(0, len(HEAVY), 60)),
        LIGHT))
    return str(path)


@pytest.fixture
def job(fasta):
    return prediction.job_from_fasta(fasta)


def test_read_fasta(fasta):
    found = prediction.read_fasta(open(fasta).read())
    assert [title for title, _ in found] == ['heavy chain', 'light chain']
    assert found[0][1] == HEAVY          # the wrapped lines joined
    # blank lines, gaps and lower case
    assert prediction.read_fasta('>a\nac-de\n\nfg*\n') == [('a', 'ACDEFG')]
    # a bare sequence, with no title line, is taken as one
    assert prediction.read_fasta('ACDEF') == [('', 'ACDEF')]
    with pytest.raises(ValueError, match='no sequences'):
        prediction.read_fasta('   ')
    with pytest.raises(ValueError, match='not amino acids'):
        prediction.read_fasta('>a\nACD123\n')


def test_chains_and_sequons(job):
    assert [c.id for c in job.chains] == ['A', 'B']
    assert job.chain('A').title == 'heavy chain'
    assert job.chain('A').sequence == HEAVY
    # N-X-S/T, X not proline
    assert job.chain('A').sequons() == [HEAVY.index('NST') + 1]
    # N-P-S does not count, N-A-T does
    assert prediction.Chain('C', 'ANPSANAT').sequons() == [6]
    with pytest.raises(ValueError, match='no chain'):
        job.chain('Z')


def test_library():
    assert ('Phosphoserine', 'SEP') in prediction.ptms_for('S')
    assert prediction.ptms_for('S') == [('Phosphoserine', 'SEP')]
    assert ('N6-acetyllysine', 'ALY') in prediction.ptms_for('K')
    assert prediction.ptm('TPO') == ('Phosphothreonine', 'T')
    names = prediction.glycan_names('n')
    for one in ('G0', 'G0F', 'G1', 'G1F', 'G2F', 'Man5', 'Core (Man3)',
                'NAG-NAG', 'NAG'):
        assert one in names
    assert prediction.glycan('NAG')[1] == [('NAG', -1, '')]
    assert prediction.glycan('NAG-NAG')[1] == [('NAG', -1, ''),
                                               ('NAG', 0, 'O4')]
    assert 'Tn antigen (GalNAc)' in prediction.glycan_names('o')
    assert 'G0F' not in prediction.glycan_names('o')
    # G0F: the core, two antennae and the core fucose
    kind, sugars = prediction.glycan('G0F')
    assert kind == 'n'
    assert [code for code, _, _ in sugars] == [
        'NAG', 'NAG', 'BMA', 'MAN', 'MAN', 'NAG', 'NAG', 'FUC']
    assert sugars[-1] == ('FUC', 0, 'O6')      # on the first NAG
    assert prediction.glycan('G1F')[1][7] == ('GAL', 5, 'O4')
    with pytest.raises(ValueError, match='unknown glycan'):
        prediction.glycan('G9')


def test_chain_ids_from_titles():
    for title, expected in (('Chain H', 'H'), ('chain_B', 'B'),
                            ('chain 1', '1'), ('H', 'H'),
                            ('heavy chain', ''), ('Chain Heavy', ''),
                            ('sp|P01857|IGHG1', '')):
        assert prediction.title_chain_id(title) == expected


def test_job_from_text():
    job = prediction.job_from_text(
        '>Chain H\nEVQLVE\n>Chain L\nDIQMTQ\n>an antigen\nMKTAYIA\n')
    assert [(c.id, c.sequence) for c in job.chains] == [
        ('H', 'EVQLVE'), ('L', 'DIQMTQ'), ('A', 'MKTAYIA')]
    assert job.chains[2].title == 'an antigen'
    # the same chain named twice: the second is lettered instead
    twice = prediction.job_from_text('>Chain H\nEVQL\n>Chain H\nDIQM\n')
    assert [c.id for c in twice.chains] == ['H', 'A']
    with pytest.raises(ValueError, match='not amino acids'):
        prediction.job_from_text('>Chain H\nEVQL123\n')


def test_what_each_residue_takes():
    # phosphorylation only on serine, threonine and tyrosine
    phospho = {code: residues for label, code, residues in prediction.PTMS
               if 'hospho' in label}
    assert phospho == {'SEP': 'S', 'TPO': 'T', 'PTR': 'Y'}
    for letter in 'ACDEFGHIKLMNPQRVW':
        assert not [code for _, code in prediction.ptms_for(letter)
                    if code in phospho]
    # deamidation, on asparagine and glutamine only, written as the
    # acid the residue becomes
    deamidation = {code: residues for label, code, residues
                   in prediction.PTMS if 'Deamidation' in label}
    assert deamidation == {'ASP': 'N', 'GLU': 'Q'}
    assert ('Deamidation (Asn to Asp)', 'ASP') in prediction.ptms_for('N')
    assert ('Deamidation (Gln to Glu)', 'GLU') in prediction.ptms_for('Q')
    for letter in 'ACDEFGHIKLMPRSTVWY':
        assert not [code for _, code in prediction.ptms_for(letter)
                    if code in deamidation]
    # glycans: N-linked on asparagine, O-linked on serine and threonine
    assert sorted(prediction.GLYCOSYLATION_ATOMS) == ['N', 'S', 'T']
    chain = prediction.Chain('A', 'NSTYKC')
    for number, letter in enumerate(chain.sequence, 1):
        for name in prediction.glycan_names('n'):
            chain.glycans = {number: name}
            if letter == 'N':
                chain.check()
            else:
                with pytest.raises(ValueError):
                    chain.check()
        for name in prediction.glycan_names('o'):
            chain.glycans = {number: name}
            if letter in 'ST':
                chain.check()
            else:
                with pytest.raises(ValueError):
                    chain.check()


def test_deamidation_is_written(job):
    chain = job.chain('A')
    site = chain.sequons()[0]              # an asparagine
    chain.modifications[site] = 'ASP'
    job.check()
    found = json.loads(prediction.write(job, 'af3'))
    protein = [e['protein'] for e in found['sequences']
               if 'protein' in e][0]
    assert protein['modifications'] == [{'ptmType': 'ASP',
                                         'ptmPosition': site}]
    entry = json.loads(prediction.write(job, 'rf3'))[0]
    seq = [c for c in entry['components'] if 'seq' in c][0]['seq']
    assert seq.startswith(chain.sequence[:site - 1] + '(ASP)')
    # a deamidated residue cannot also carry a glycan
    chain.glycans[site] = 'NAG'
    with pytest.raises(ValueError, match='both a modification and a '
                       'glycan'):
        job.check()


def test_checks(job):
    chain = job.chain('A')
    chain.modifications[5] = 'SEP'           # residue 5 is V
    with pytest.raises(ValueError, match='replaces S'):
        job.check()
    del chain.modifications[5]
    chain.glycans[2] = 'G0F'                 # residue 2 is V
    with pytest.raises(ValueError, match='glycans go on N'):
        job.check()
    del chain.glycans[2]
    chain.glycans[7] = 'G0F'                 # residue 7 is S: O-linked only
    with pytest.raises(ValueError, match='N-linked'):
        job.check()
    del chain.glycans[7]
    job.disulfides.append((('A', 1), ('A', 96)))
    with pytest.raises(ValueError, match='not a cysteine'):
        job.check()
    job.disulfides.clear()
    job.ligands.append(('not a code', 1))
    with pytest.raises(ValueError, match='CCD code'):
        job.check()


def decorate(job):
    chain = job.chain('A')
    chain.glycans[chain.sequons()[0]] = 'G0F'
    chain.modifications[7] = 'SEP'
    job.chain('B').glycans[9] = 'T antigen (core 1)'
    job.disulfides.append((('A', 22), ('A', 96)))
    job.ligands.append(('atp', 2))
    job.check()
    return job


def test_alphafold3(job):
    found = json.loads(prediction.write(decorate(job), 'af3'))
    assert found['dialect'] == 'alphafold3'
    assert found['version'] == prediction.AF3_VERSION
    assert found['modelSeeds'] == [1]
    proteins = {e['protein']['id']: e['protein']
                for e in found['sequences'] if 'protein' in e}
    assert sorted(proteins) == ['A', 'B']
    assert proteins['A']['sequence'] == HEAVY      # unmodified
    assert proteins['A']['modifications'] == [
        {'ptmType': 'SEP', 'ptmPosition': 7}]      # no CCD_ prefix
    assert 'modifications' not in proteins['B']
    ligands = [e['ligand'] for e in found['sequences'] if 'ligand' in e]
    glycan = [one for one in ligands if one['ccdCodes'][0] == 'NAG'][0]
    assert glycan['ccdCodes'] == ['NAG', 'NAG', 'BMA', 'MAN', 'MAN',
                                  'NAG', 'NAG', 'FUC']
    atp = [one for one in ligands if one['ccdCodes'] == ['ATP']][0]
    assert len(atp['id']) == 2 and atp['id'][0] != 'A'    # two copies
    # the glycan hangs off the asparagine, and its sugars off each other
    bonds = found['bondedAtomPairs']
    site = job.chain('A').sequons()[0]
    assert [['A', site, 'ND2'], [glycan['id'], 1, 'C1']] in bonds
    assert [[glycan['id'], 1, 'O4'], [glycan['id'], 2, 'C1']] in bonds
    assert [[glycan['id'], 3, 'O6'], [glycan['id'], 5, 'C1']] in bonds
    assert [[glycan['id'], 1, 'O6'], [glycan['id'], 8, 'C1']] in bonds
    o_glycan = [one for one in ligands if one['ccdCodes'][0] == 'NGA'][0]
    assert [['B', 9, 'OG'], [o_glycan['id'], 1, 'C1']] in bonds
    # the disulfide, as an SG-SG pair (its reader takes those, whatever
    # docs/input.md says about bonds between polymer entities)
    assert [['A', 22, 'SG'], ['A', 96, 'SG']] in bonds
    assert not job.warnings('af3')


def test_rosettafold3(job):
    found = json.loads(prediction.write(decorate(job), 'rf3'))
    assert isinstance(found, list) and len(found) == 1
    entry = found[0]
    chains = {c['chain_id']: c for c in entry['components'] if 'seq' in c}
    # modified residues are written into the sequence
    assert chains['A']['seq'] == HEAVY[:6] + '(SEP)' + HEAVY[7:]
    assert chains['B']['seq'] == LIGHT
    sugars = [c for c in entry['components'] if 'ccd_code' in c]
    assert [c['ccd_code'] for c in sugars].count('NAG') == 4
    assert [c['ccd_code'] for c in sugars].count('ATP') == 2
    site = job.chain('A').sequons()[0]
    first = sugars[0]['chain_id']
    assert ['A/ASN/%d/ND2' % site, '%s/NAG/1/C1' % first] in entry['bonds']
    # RoseTTAFold3 takes bonds freely, so the disulfide is written
    assert ['A/CYS/22/SG', 'A/CYS/96/SG'] in entry['bonds']
    assert not job.warnings('rf3')


def test_esmfold2(job):
    found = json.loads(prediction.write(decorate(job), 'esmfold2'))
    assert found['model'] == prediction.ESMFOLD2_MODEL
    proteins = [e['protein'] for e in found['sequences'] if 'protein' in e]
    assert proteins[0]['modifications'] == [{'position': 7, 'ccd': 'SEP'}]
    ligands = [e['ligand'] for e in found['sequences'] if 'ligand' in e]
    assert ['ATP'] in [one['ccd'] for one in ligands]
    # it is the SDK's shape, not a documented file: say so
    assert any('SDK' in line for line in job.warnings('esmfold2'))


def test_disulfides_follow_the_copies(job):
    job.chain('A').copies = 2
    job.disulfides.append((('A', 22), ('A', 96)))       # inside one chain
    job.disulfides.append((('A', 22), ('B', 23)))       # between two
    pairs = prediction.disulfide_pairs(job)
    ids = [(a[0], b[0]) for a, b in pairs]
    # the inner bond is made in both copies, the outer one once
    assert ids == [('A', 'A'), ('C', 'C'), ('A', 'B')]
    found = json.loads(prediction.write(job, 'af3'))
    assert [['C', 22, 'SG'], ['C', 96, 'SG']] in found['bondedAtomPairs']
    entry = json.loads(prediction.write(job, 'rf3'))[0]
    assert ['C/CYS/22/SG', 'C/CYS/96/SG'] in entry['bonds']


def test_copies_get_their_own_chains(job):
    job.chain('A').copies = 2
    job.chain('A').glycans[job.chain('A').sequons()[0]] = 'NAG'
    found = json.loads(prediction.write(job, 'af3'))
    proteins = [e['protein'] for e in found['sequences'] if 'protein' in e]
    assert [p['id'] for p in proteins] == ['A', 'C', 'B']
    assert proteins[0]['sequence'] == proteins[1]['sequence']
    # each copy carries its own glycan, on its own chain
    ligands = [e['ligand']['id'] for e in found['sequences']
               if 'ligand' in e]
    assert len(ligands) == 2 and len(set(ligands)) == 2
    assert len({bond[0][0] for bond in found['bondedAtomPairs']}) == 2


def test_save(job, tmp_path):
    path = str(tmp_path / 'job.json')
    prediction.save(decorate(job), 'af3', path)
    assert json.load(open(path))['name'] == 'ab'
    with pytest.raises(ValueError, match='method'):
        prediction.write(job, 'esmfold9')
    empty = prediction.Job()
    with pytest.raises(ValueError, match='no sequences'):
        prediction.write(empty, 'af3')


@pytest.fixture(scope='module')
def gui():
    '''
    The dialog the GUI tests share: building several of them in one PyMOL
    leaves it crashing on quit, and one is enough
    '''
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import pymol.plugins  # noqa: F401 (adds the startup plugin path)
    from pmg_tk.startup.prediction_gui.main import PredictionDialog
    from pymolx.gui.sequence_viewer import MemoryPrefs
    dialog = PredictionDialog(None, _self=cmd)
    dialog.prefs = MemoryPrefs()           # (not the real preferences)
    yield dialog
    import shiboken6
    dialog.close()
    shiboken6.delete(dialog)          # (see test_submit_dialog)
    QtWidgets.QApplication.processEvents()


@pytest.fixture
def fresh(gui, fasta):
    '''
    The shared dialog, with the test FASTA loaded and nothing on it
    '''
    gui.job = prediction.job_from_fasta(fasta)
    gui.saved = ''
    gui.typed = ''
    gui.message = ''
    gui.pending = None
    gui.prefs.clear()
    gui.method.setCurrentIndex(gui.method.findData('af3'))
    gui.refresh()
    return gui


SCRIPT = """#!/bin/bash
LOCAL_DATA_DIR="/data"
AF3_DIR="$LOCAL_DATA_DIR/alphafold3"
docker run -it \\
  --volume $AF3_DIR/inputs_outputs:/root/af_input_output \\
  --volume "$AF3_DIR/weights:/root/models" \\
  alphafold3 python /root/repo/run_alphafold.py --json_path=/root/af_input_output/$1
"""


def test_script_input_dir(tmp_path):
    path = tmp_path / 'run_alphafold3.sh'
    path.write_text(SCRIPT)
    # the folder the script mounts as its input/output, variables resolved
    assert prediction.script_input_dir(str(path)) == \
        '/data/alphafold3/inputs_outputs'
    assert prediction.shell_variables(SCRIPT)['AF3_DIR'] == \
        '/data/alphafold3'
    # a script that says nothing, and one that is not there
    other = tmp_path / 'quiet.sh'
    other.write_text('echo hello\n')
    assert prediction.script_input_dir(str(other)) == ''
    assert prediction.script_input_dir(str(tmp_path / 'nope.sh')) == ''


def test_submit_command(tmp_path):
    path = str(tmp_path / 'run.sh')
    open(path, 'w').write(SCRIPT)
    # a single run goes through script(1), which gives docker a terminal
    assert prediction.submit_command(path, 'ab_af3.json') == [
        'script', '-qec', 'bash %s ab_af3.json' % path, '/dev/null']
    assert prediction.submit_command(path, 'ab_af3.json', pty=False) == [
        'bash', path, 'ab_af3.json']
    assert prediction.submit_command(path, 'ab_af3.json', 'slurm') == [
        'sbatch', path, 'ab_af3.json']
    # a name with spaces is quoted for the shell script(1) runs
    assert "'my job.json'" in prediction.submit_command(
        path, 'my job.json')[2]
    with pytest.raises(ValueError, match='job script'):
        prediction.submit_command('', 'ab.json')
    with pytest.raises(ValueError, match='no such job script'):
        prediction.submit_command(str(tmp_path / 'nope.sh'), 'ab.json')
    with pytest.raises(ValueError, match="input file's name"):
        prediction.submit_command(path, '')
    with pytest.raises(ValueError, match='name only'):
        prediction.submit_command(path, '/somewhere/ab.json')
    with pytest.raises(ValueError, match='single or slurm'):
        prediction.submit_command(path, 'ab.json', 'cloud')


def test_submit_dialog(fresh, tmp_path):
    import shiboken6
    from pmg_tk.startup.prediction_gui.main import SubmitDialog
    from pymol.Qt import QtWidgets

    folder = tmp_path / 'alphafold3' / 'inputs_outputs'
    folder.mkdir(parents=True)
    script = tmp_path / 'run.sh'
    script.write_text(SCRIPT.replace('"/data"', '"%s"' % tmp_path))
    saved = tmp_path / 'ab_af3.json'
    saved.write_text('{}')

    # submitting is for AlphaFold 3 so far
    assert fresh.submit_button.isEnabled()
    fresh.method.setCurrentIndex(fresh.method.findData('rf3'))
    assert not fresh.submit_button.isEnabled()
    fresh.method.setCurrentIndex(fresh.method.findData('af3'))

    sub = SubmitDialog('ab_af3.json', None, fresh.prefs)
    try:
        sub.saved = str(saved)
        sub.script.setText(str(script))
        # it knows where the script looks, and the file is not there yet
        assert sub.folder() == str(folder)
        assert sub.where.text().startswith('Not in')
        assert sub.copy_button.isEnabled()
        sub._copy()
        assert (folder / 'ab_af3.json').exists()
        assert sub.where.text().startswith('Found in')
        assert not sub.copy_button.isEnabled()
        # the command it shows, for each way of running it
        assert sub.command.toPlainText() == \
            'script -qec %s /dev/null' % ('bash %s ab_af3.json' % script)
        assert sub.submit.isEnabled()
        sub.modes.button(1).setChecked(True)      # Slurm
        assert sub.mode() == 'slurm'
        assert sub.command.toPlainText() == 'sbatch %s ab_af3.json' % script
        # a missing name is refused, with the reason shown
        sub.input_name.clear()
        assert not sub.submit.isEnabled()
        assert "input file's name" in sub.command.toPlainText()
        sub.input_name.setText('ab_af3.json')
        # the script is remembered per way of running it
        sub.accept()
        assert fresh.prefs[SubmitDialog.SCRIPT_PREF % 'slurm'] == str(script)
        assert SubmitDialog.SCRIPT_PREF % 'single' not in fresh.prefs
    finally:
        # destroyed now, not queued: a dialog still alive when PyMOL
        # quits crashes it (offscreen, after the unittest phase)
        sub.close()
        shiboken6.delete(sub)
        QtWidgets.QApplication.processEvents()


def test_submit_runs_the_command(fresh, tmp_path):
    from pymol.Qt import QtWidgets

    # a stand-in for the job script: it reports what it was given
    script = tmp_path / 'run.sh'
    script.write_text('#!/bin/bash\necho "predicting $1 in $PWD"\n')
    fresh._run(prediction.submit_command(str(script), 'ab_af3.json'),
               'single', str(tmp_path))
    assert fresh.process is not None
    assert not fresh.submit_button.isEnabled()      # one at a time
    assert fresh.process.waitForFinished(20000)
    for _ in range(20):
        QtWidgets.QApplication.processEvents()
        if fresh.process is None:
            break
    assert fresh.process is None
    assert 'finished' in fresh.message
    assert fresh.submit_button.isEnabled()


def test_dialog(fresh, fasta):
    QtCore = pytest.importorskip('pymol.Qt').QtCore
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets

    dialog = fresh
    if True:
        # sequences typed in instead of opened
        from pmg_tk.startup.prediction_gui.main import SequenceTextDialog
        typing = SequenceTextDialog('', dialog)
        ok = QtWidgets.QDialogButtonBox.StandardButton.Ok
        assert not typing.buttons.button(ok).isEnabled()
        typing.text.setPlainText('>Chain H\nEVQL123\n')
        assert 'not amino acids' in typing.status.text()
        assert not typing.buttons.button(ok).isEnabled()
        typing.text.setPlainText('>Chain H\nEVQLVESGGG\n>Chain L\nDIQMT\n')
        assert typing.buttons.button(ok).isEnabled()
        assert 'H: 10 residues' in typing.status.text()
        typing.deleteLater()
        # a title that only names its chain is not repeated in the heading
        named = prediction.job_from_text('>Chain H\nEVQLVE\n')
        dialog.job = named
        dialog.refresh()
        assert 'Chain H&nbsp;&nbsp;<span' in dialog.view._html()
        dialog.job = prediction.job_from_text('>Chain H kappa\nEVQLVE\n')
        dialog.refresh()
        assert 'Chain H kappa' in dialog.view._html()
        dialog.job = prediction.job_from_fasta(fasta)
        dialog.refresh()
        chain = dialog.job.chain('A')
        site = chain.sequons()[0]
        # the sequence is drawn with a link per residue, the site marked
        text = dialog.view._html()
        assert 'href="#A:%d"' % site in text
        assert 'N-X-S/T site' in text
        # the menu offers what suits the residue
        def labels(chain_id, number):
            return [a.text() for a in
                    dialog.residue_menu(chain_id, number).actions()]
        # a click sends the chain as it is, upper case (Qt lower-cases
        # a URL scheme, so the link is a fragment)
        seen = []
        dialog.view.clicked.disconnect()      # (its own slot opens a menu)
        dialog.view.clicked.connect(lambda c, n: seen.append((c, n)))
        dialog.view._clicked(QtCore.QUrl('#A:%d' % site))
        assert seen == [('A', site)]
        assert 'Phosphoserine  (SEP)' in labels('A', 7)       # a serine
        assert not [one for one in labels('A', 5)             # a valine
                    if 'Phospho' in one]
        assert 'N-glycan' in [a.text() for a in
                              dialog.residue_menu('A', site).actions()]
        assert 'Disulfide bond from here...' in labels('A', 22)

        # the glycan submenu follows the residue
        def glycans(chain_id, number):
            found = [a for a in dialog.residue_menu(chain_id,
                                                    number).actions()
                     if a.menu()]
            return (found[0].text(),
                    [a.text() for a in found[0].menu().actions()]) \
                if found else ('', [])

        title, names = glycans('A', site)
        assert title == 'N-glycan' and 'NAG-NAG' in names
        title, names = glycans('A', HEAVY.index('S') + 1)
        assert title == 'O-glycan' and 'G0F' not in names
        assert glycans('A', 5) == ('', [])        # a valine takes none
        # adding and removing
        dialog._set_glycan(chain, site, 'G0F')
        dialog._set_ptm(chain, 7, 'SEP')
        dialog._pend(('A', 22))
        assert 'click its partner' in dialog.status.text()
        dialog._add_disulfide(('A', 96))
        dialog.ligand.setText('atp')
        dialog.add_ligand()
        assert dialog.table.rowCount() == 3
        assert dialog.ligand_list.count() == 1
        # the disulfide reaches the file, AlphaFold 3's included
        dialog.job.seeds = [1]
        found = json.loads(prediction.write(dialog.job, 'af3'))
        assert [['A', 22, 'SG'], ['A', 96, 'SG']] in \
            found['bondedAtomPairs']
        dialog.method.setCurrentIndex(dialog.method.findData('rf3'))
        assert not dialog.seeds.isEnabled()       # seeds are AlphaFold's
        dialog.method.setCurrentIndex(dialog.method.findData('esmfold2'))
        assert 'SDK' in dialog.status.text()
        # the table removes what it names
        dialog.table.selectRow(0)
        dialog.remove_selected()
        assert dialog.table.rowCount() == 2
        assert not dialog.grab().isNull()
