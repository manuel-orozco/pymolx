'''
Paratope prediction (pymolx.paratope, data/startup/paratope_gui, parity
item L-32): Parapred on a sequence or a structure, Paragraph and
ParaSurf on a structure, and the dialog.

ParaSurf takes minutes without a GPU, so the run itself only happens
with PYMOLX_SLOW_TESTS=1; what surrounds it is always tested.
'''

import os

import pytest

from pymol import cmd
import pymol.plugins  # noqa: F401 (puts $PYMOL_DATA/startup on the path)
from pymolx import paratope

VH = ('EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG'
      'RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS')
VL = ('DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR'
      'SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK')

needs_parapred = pytest.mark.skipif(
    bool(paratope.available('parapred')),
    reason=paratope.available('parapred') or '')
needs_paragraph = pytest.mark.skipif(
    bool(paratope.available('paragraph')),
    reason=paratope.available('paragraph') or '')
needs_parasurf = pytest.mark.skipif(
    bool(paratope.available('parasurf')),
    reason=paratope.available('parasurf') or '')
slow = pytest.mark.skipif(not os.environ.get('PYMOLX_SLOW_TESTS'),
                          reason='minutes without a GPU: '
                                 'PYMOLX_SLOW_TESTS=1 runs it')


@pytest.fixture(autouse=True)
def session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


def test_methods_listed():
    assert [v for _, v in paratope.METHODS] == ['parapred', 'paragraph',
                                                'parasurf']
    assert paratope.STRUCTURE_METHODS == ('paragraph', 'parasurf')
    assert paratope.method_label('paragraph').startswith('Paragraph')
    assert paratope.method_label('parasurf').startswith('ParaSurf')
    assert 'one of' in paratope.available('nonsense')
    assert 0 < paratope.THRESHOLD < 1


@needs_parapred
def test_extended_cdrs():
    from pymolx import cdr
    domain, = cdr.domains(VH, 'chothia')
    blocks = paratope.extended_cdrs(domain)
    assert [name for name, _, _ in blocks] == list(cdr.CDR_NAMES)
    loops = domain.cdrs()
    for name, letters, places in blocks:
        # each CDR with up to FLANK residues either side
        assert loops[name] in letters
        assert len(letters) <= len(loops[name]) + 2 * paratope.FLANK
        assert len(places) == len(letters)
        assert len(letters) <= paratope.MAX_LENGTH


@needs_parapred
def test_parapred_on_a_sequence():
    residues = paratope.predict_sequence(VH, 'VH')
    assert residues
    assert all(0.0 <= r.score <= 1.0 for r in residues)
    assert {r.region for r in residues} == {'CDR1', 'CDR2', 'CDR3'}
    assert all(r.chain == 'VH' for r in residues)
    # the residues it scored are the ones it was given, in the sequence
    for r in residues:
        assert VH[int(r.resi) - 1] == r.letter
    # the likeliest binders are in the loops, with high scores
    best = max(residues, key=lambda r: r.score)
    assert best.binds and best.region in ('CDR1', 'CDR2', 'CDR3')
    assert sum(1 for r in residues if r.binds) >= 5
    # a light chain works too, and is scored separately
    assert paratope.predict_sequence(VL, 'VL')


@needs_parapred
def test_parapred_on_a_structure():
    cmd.fab(VH, 'ab', chain='H', ss=1)
    residues = cmd.predict_paratope('ab', method='parapred', quiet=1)
    assert residues and all(r.chain == 'H' for r in residues)
    # the scores reached the B-factors and the selection
    assert 'ab_paratope' in cmd.get_names('selections')
    binders = [r for r in residues if r.binds]
    assert cmd.count_atoms('ab_paratope and name CA') == len(binders)
    rows = []
    cmd.iterate('ab and name CA and b > 0.9', 'rows.append(resi)',
                space={'rows': rows})
    assert rows == [r.resi for r in residues if r.score > 0.9]
    text = paratope.report(residues, 'parapred')
    assert 'Parapred' in text and 'CDR3' in text


@needs_paragraph
def test_paragraph_on_a_structure(tmp_path):
    cmd.fab(VH, 'h', chain='H', ss=1)
    cmd.fab(VL, 'l', chain='L', ss=1)
    cmd.create('fv', 'h or l')
    cmd.delete('h l')
    assert paratope.antibody_pairs('fv') == [('H', 'L')]
    residues = paratope.paragraph('fv')
    assert residues
    assert {r.chain for r in residues} <= {'H', 'L'}
    assert all(0.0 <= r.score <= 1.0 for r in residues)
    # every score maps back onto a residue the structure has
    for r in residues:
        assert cmd.count_atoms('fv and chain %s and resi %s and name CA'
                               % (r.chain, r.resi)) == 1
    # and no residue is scored twice
    keys = [(r.chain, r.resi) for r in residues]
    assert len(keys) == len(set(keys))


@needs_paragraph
def test_paragraph_scores_each_copy_alike():
    # two copies of one Fv: renumbering them to IMGT one residue at a
    # time merged residues, and the copies then scored differently
    cmd.fab(VH, 'h', chain='H', ss=1)
    cmd.fab(VL, 'l', chain='L', ss=1)
    cmd.create('two', 'h or l')
    cmd.create('copy2', 'two')
    cmd.alter('copy2', 'chain = {"H": "K", "L": "M"}.get(chain, chain)')
    cmd.create('two', 'two or copy2')
    cmd.delete('h l copy2')
    cmd.sort('two')
    pairs = paratope.antibody_pairs('two')
    assert sorted(pairs) == [('H', 'L'), ('K', 'M')]
    residues = paratope.paragraph('two')
    counts = {}
    for r in residues:
        counts[r.chain] = counts.get(r.chain, 0) + 1
    assert counts.get('H') == counts.get('K')
    assert counts.get('L') == counts.get('M')


@needs_parasurf
def test_parasurf_is_found():
    # the clone, its weights and DMS, which it computes the surface with
    home = paratope.parasurf_home()
    assert os.path.isdir(os.path.join(home, 'ParaSurf'))
    assert os.path.isfile(paratope.parasurf_weights('paired'))
    assert os.access(paratope.parasurf_dms(), os.X_OK)
    # a kind of Fv it has no weights for falls back on the paired ones
    for kind in ('heavy', 'light', 'nonsense'):
        assert os.path.isfile(paratope.parasurf_weights(kind))
    assert 0.1 <= paratope.PARASURF_MESH <= 1.0
    # the GPU when there is one, which is the difference between
    # seconds and minutes
    assert paratope.parasurf_device() in ('cuda', 'cpu')


@needs_parapred        # ANARCI, which finds the variable domains
def test_parasurf_writes_one_pdb_per_fv():
    import shutil
    from pymolx import cdr
    cmd.fab(VH, 'h', chain='A', ss=1)
    cmd.fab(VL, 'l', chain='B', ss=1)
    cmd.create('fv', 'h or l')
    cmd.delete('h l')
    cmd.alter('fv', 'b = 0')
    folder, jobs, mapping = paratope.parasurf_inputs('fv')
    try:
        assert jobs == [('fv_0', 'paired')]
        path = os.path.join(folder, 'fv_0.pdb')
        lines = [l for l in open(path) if l.startswith(('ATOM', 'HETATM'))]
        assert lines
        # ParaSurf wants two chains, named H and L, without hydrogens or
        # het atoms, and the atoms numbered in order
        assert {l[21] for l in lines} == {'H', 'L'}
        assert not [l for l in lines if l.startswith('HETATM')]
        assert not [l for l in lines if l[12:16].strip().startswith('H')]
        assert [int(l[6:11]) for l in lines] == list(range(1, len(lines) + 1))
        # the chains it writes map back onto the structure's own
        assert {where[0] for where in mapping.values()} == {'A', 'B'}
        assert {key[1] for key in mapping} == {'H', 'L'}
        for (_, letter, resi), (chain, back, region) in mapping.items():
            assert back == resi
            assert chain == ('A' if letter == 'H' else 'B')
            assert cmd.count_atoms('fv and chain %s and resi %s and name CA'
                                   % (chain, resi)) == 1
        # and ANARCI's regions come with them, for the table
        assert {where[2] for where in mapping.values()} >= set(cdr.CDR_NAMES)
        # and the scores come back on those chains
        rows = [('fv_0', letter, resi, 'ALA', 0.9)
                for _, letter, resi in list(mapping)[:3]]
        rows.append(('fv_0', 'H', '9999', 'ALA', 0.9))   # not in the Fv
        residues = paratope.parasurf_residues(rows, mapping)
        assert len(residues) == 3
        # each one keeps the region ANARCI put it in
        for r in residues:
            assert r.region == mapping[('fv_0',
                                        'H' if r.chain == 'A' else 'L',
                                        r.resi)][2]
        assert {r.chain for r in residues} <= {'A', 'B'}
        assert all(r.score == 0.9 for r in residues)
    finally:
        shutil.rmtree(folder, ignore_errors=True)


@slow
@needs_parasurf
def test_parasurf_on_a_structure():
    # one chain, the coarsest surface: the whole run still takes
    # minutes on a CPU (a folded Fv is quicker than this built one)
    cmd.fab(VH, 'vh', chain='H', ss=1)
    residues = paratope.parasurf('vh', mesh_dense=0.1)
    assert residues
    assert {r.chain for r in residues} == {'H'}
    assert all(0.0 <= r.score <= 1.0 for r in residues)
    assert {r.region for r in residues} >= set(cdr_names())
    # every score lands on a residue of the structure, and only once
    keys = [(r.chain, r.resi) for r in residues]
    assert len(keys) == len(set(keys))
    for r in residues:
        assert cmd.count_atoms('vh and chain %s and resi %s and name CA'
                               % (r.chain, r.resi)) == 1
    # the scores go on the model like the other predictors'
    paratope.apply_scores('vh', residues, _self=cmd)
    assert 'vh_paratope' in cmd.get_names('selections')


def cdr_names():
    from pymolx import cdr
    return cdr.CDR_NAMES


@needs_parapred
def test_dialog():
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import shiboken6
    from pmg_tk.startup.paratope_gui.main import ParatopeDialog

    cmd.fab(VH, 'ab', chain='H', ss=1)
    dialog = ParatopeDialog(None, _self=cmd)
    try:
        dialog.refresh()
        dialog.object.setCurrentText('ab')
        dialog.method.setCurrentIndex(dialog.method.findData('parapred'))
        dialog.run()
        assert dialog.table.rowCount() == len(dialog.residues)
        assert 'Parapred' in dialog.status.text()
        # ranked by score, the best first
        scores = [float(dialog.table.item(r, 2).text())
                  for r in range(dialog.table.rowCount())]
        assert scores == sorted(scores, reverse=True)
        # the threshold column follows the spin box without predicting again
        dialog.threshold.setValue(0.99)
        few = sum(1 for r in range(dialog.table.rowCount())
                  if dialog.table.item(r, 3).text() == 'yes')
        dialog.threshold.setValue(0.10)
        many = sum(1 for r in range(dialog.table.rowCount())
                   if dialog.table.item(r, 3).text() == 'yes')
        assert few < many

        # Paragraph needs a structure: a pasted sequence is refused
        dialog.source.button(1).setChecked(True)
        dialog.method.setCurrentIndex(dialog.method.findData('paragraph'))
        assert not dialog.run_button.isEnabled()
        assert 'needs a structure' in dialog.status.text()
        dialog.method.setCurrentIndex(dialog.method.findData('parapred'))
        assert dialog.run_button.isEnabled()
        dialog.sequence.setPlainText(VL)
        dialog.run()
        assert dialog.residues and dialog.table.rowCount()
    finally:
        dialog.close()
        shiboken6.delete(dialog)
        QtWidgets.QApplication.processEvents()


@needs_parasurf
def test_dialog_offers_parasurf():
    QtWidgets = pytest.importorskip('pymol.Qt').QtWidgets
    if QtWidgets.QApplication.instance() is None:
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        QtWidgets.QApplication(['pymolx-test'])
    import shiboken6
    from pmg_tk.startup.paratope_gui.main import ParatopeDialog

    cmd.fab(VH, 'ab', chain='H', ss=1)
    dialog = ParatopeDialog(None, _self=cmd)
    try:
        dialog.refresh()
        dialog.object.setCurrentText('ab')
        # the surface detail only belongs to ParaSurf
        dialog.method.setCurrentIndex(dialog.method.findData('paragraph'))
        assert dialog.mesh_row.isHidden()
        dialog.method.setCurrentIndex(dialog.method.findData('parasurf'))
        assert not dialog.mesh_row.isHidden()
        assert dialog.mesh.value() == paratope.PARASURF_MESH
        assert 'ParaSurf' in dialog.about.text()
        assert dialog.run_button.isEnabled()
        # it reads a structure, not a pasted sequence
        dialog.source.button(1).setChecked(True)
        assert not dialog.run_button.isEnabled()
        assert dialog.status.text() == 'ParaSurf needs a structure'
        dialog.source.button(0).setChecked(True)
        assert dialog.run_button.isEnabled()
        # nothing is running, so no worker and no leftover folder
        assert dialog.thread is None
    finally:
        dialog.close()
        shiboken6.delete(dialog)
        QtWidgets.QApplication.processEvents()
