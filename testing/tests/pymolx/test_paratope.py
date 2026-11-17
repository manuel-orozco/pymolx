'''
Paratope prediction (pymolx.paratope, data/startup/paratope_gui, parity
item L-32): Parapred on a sequence or a structure, Paragraph on a
structure, and the dialog.
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


@pytest.fixture(autouse=True)
def session():
    cmd.reinitialize()
    yield
    cmd.reinitialize()


def test_methods_listed():
    assert [v for _, v in paratope.METHODS] == ['parapred', 'paragraph']
    assert paratope.method_label('paragraph').startswith('Paragraph')
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
