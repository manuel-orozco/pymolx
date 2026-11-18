'''
Paratope prediction (parity item L-32): which residues of an antibody
bind its antigen, from two published predictors.

- Parapred (Liberis et al. 2018), the PyTorch port from Alchemab: a
  convolutional and recurrent network over the sequence of each CDR with
  two residues either side ("extended CDRs"). Needs a sequence only.
- Paragraph (Chinery et al. 2023, OPIG): an equivariant graph network
  over the structure's Fv. Needs a structure, IMGT numbered, which is
  done here with ANARCI (pymolx.cdr) so any numbering can be given.
- ParaSurf (Papadopoulos et al. 2025): a 3D convolutional network with
  transformer blocks over the molecular surface of the Fv, which DMS
  computes and pdb2pqr gives force field features for. Needs a
  structure and the cloned repository with its weights (see
  parasurf_home), and takes minutes on a CPU.

All three score residues from 0 to 1. Parapred sees the extended CDRs
and Paragraph the CDR vicinity, so a residue without a score is one the
predictor did not consider, which is not the same as one it judged not
to bind; ParaSurf scores every residue of the Fv.
'''
import csv
import os
import re
import shutil
import sys
import tempfile

from pymol import cmd

from . import cdr

METHODS = (('Parapred (sequence)', 'parapred'),
           ('Paragraph (structure)', 'paragraph'),
           ('ParaSurf (structure)', 'parasurf'))
DEFAULT_METHOD = 'parapred'
# the methods that read a structure rather than a sequence
STRUCTURE_METHODS = ('paragraph', 'parasurf')
# Parapred reads each CDR with this many residues either side
FLANK = 2
# the longest input its encoder takes
MAX_LENGTH = 40
# a residue counts as paratope above this (the papers' usual threshold)
THRESHOLD = 0.67
PARATOPE_COLOR = 'red'
REST_COLOR = 'grey70'
# ParaSurf: points per square angstrom of surface DMS computes. Its own
# default; 0.1 is three times faster and a little coarser
PARASURF_MESH = 0.3
# its weights, by the kind of Fv they were trained on
PARASURF_WEIGHTS = (('paired', 'entire'), ('heavy', 'heavy'),
                    ('light', 'light'))
# its network, as blind_predict.py configures it
PARASURF_GRID = 41
PARASURF_CHANNELS = 22
PARASURF_BATCH = 64


def method_label(method):
    return dict((v, k) for k, v in METHODS).get(method, method)


def available(method=DEFAULT_METHOD):
    '''
    :return: '' if the method can be used, else why not
    '''
    missing = cdr.available()
    if missing:
        return missing
    try:
        import torch  # noqa: F401
    except ImportError:
        return ('PyTorch is not installed: "pip install torch" in the '
                'environment PyMOL runs in')
    if method == 'parapred':
        try:
            import parapred  # noqa: F401
        except ImportError:
            return ('Parapred is not installed: "pip install '
                    'git+https://github.com/alchemab/parapred-pytorch"')
        if not _parapred_weights():
            return ('Parapred\'s weights are missing: put '
                    'parapred_pytorch.h5 in the package\'s weights folder '
                    '(its wheel leaves it out)')
    elif method == 'paragraph':
        try:
            import Paragraph  # noqa: F401
        except ImportError:
            return ('Paragraph is not installed: "pip install '
                    'git+https://github.com/oxpig/Paragraph"')
    elif method == 'parasurf':
        return _parasurf_missing()
    else:
        return 'method: one of %s' % ', '.join(v for _, v in METHODS)
    return ''


def _parasurf_missing():
    '''
    :return: '' if ParaSurf can be run, else what is missing
    '''
    home = parasurf_home()
    if not home:
        return ('ParaSurf is not installed: clone '
                'https://github.com/aggelos-michael-papadopoulos/ParaSurf'
                ' next to the PyMOL environment, or set PARASURF_HOME')
    if not parasurf_weights():
        return ('ParaSurf\'s weights are missing: put the ".pth" files '
                'from its README in "%s", or set PARASURF_WEIGHTS'
                % os.path.join(home, 'weights'))
    try:
        from openbabel import pybel  # noqa: F401
    except ImportError:
        return ('Open Babel is not installed: "pip install '
                'openbabel-wheel" in the environment PyMOL runs in')
    if not parasurf_dms():
        return ('DMS is missing (ParaSurf computes the surface with it): '
                'build it with "make" in "%s" and put the binary on the '
                'PATH' % os.path.join(home, 'dms'))
    return ''


def parasurf_home():
    '''
    The cloned ParaSurf repository: $PARASURF_HOME, one already on the
    path, or a "ParaSurf" folder beside the Python environment or in the
    home directory

    :return: its path, or ''
    '''
    places = [os.environ.get('PARASURF_HOME', '')]
    try:
        # its package has no __init__.py, so only __path__ is reliable
        import ParaSurf
        for where in list(getattr(ParaSurf, '__path__', None) or []):
            places.append(os.path.dirname(os.path.abspath(where)))
    except ImportError:
        pass
    places += [os.path.join(os.path.dirname(sys.prefix), 'ParaSurf'),
               os.path.expanduser('~/ParaSurf')]
    for place in places:
        # the repository holds the package of the same name
        if place and os.path.isdir(os.path.join(place, 'ParaSurf')):
            return place
    return ''


def parasurf_weights(kind='paired'):
    '''
    The weights file for an Fv of this kind ("paired", "heavy" or
    "light"), falling back on the paired ones

    :return: its path, or ''
    '''
    where = os.environ.get('PARASURF_WEIGHTS', '')
    if where and os.path.isfile(where):
        return where
    if not (where and os.path.isdir(where)):
        home = parasurf_home()
        where = os.path.join(home, 'weights') if home else ''
    if not (where and os.path.isdir(where)):
        return ''
    files = sorted(f for f in os.listdir(where) if f.endswith('.pth'))
    wanted = dict(PARASURF_WEIGHTS).get(kind, 'entire')
    for mark in (wanted, 'entire'):
        for name in files:
            if mark in name.lower():
                return os.path.join(where, name)
    return os.path.join(where, files[0]) if files else ''


def parasurf_dms():
    '''
    The DMS binary ParaSurf shells out to: on the PATH, or in the
    environment's bin directory (where "make install" can put it)

    :return: its path, or ''
    '''
    found = shutil.which('dms')
    if found:
        return found
    beside = os.path.join(sys.prefix, 'bin', 'dms')
    return beside if os.path.exists(beside) else ''


def _parapred_weights():
    try:
        import parapred
    except ImportError:
        return ''
    path = os.path.join(os.path.dirname(parapred.__file__), 'weights',
                        'parapred_pytorch.h5')
    return path if os.path.exists(path) else ''


class Residue:
    '''
    One scored residue

    :ivar chain: the chain it is in (a structure's, or the sequence's name)
    :ivar resi: its number in the structure, or its 1-based place in the
        sequence
    :ivar letter: one-letter code
    :ivar score: 0 to 1
    :ivar region: the loop it belongs to (CDR1 ...), or ''
    '''

    def __init__(self, chain, resi, letter, score, region=''):
        self.chain = chain
        self.resi = str(resi)
        self.letter = letter
        self.score = float(score)
        self.region = region

    @property
    def binds(self):
        return self.score >= THRESHOLD

    def __repr__(self):
        return '<%s %s%s %.2f>' % (self.letter, self.chain, self.resi,
                                   self.score)


def extended_cdrs(domain, residues=()):
    '''
    Each CDR with FLANK residues either side, as Parapred reads them

    :param domain: a cdr.Domain
    :param residues: [(chain, resi)] of the domain, or () for a sequence
    :return: [(region name, sequence, [(chain, resi) or index])]
    '''
    regions = domain.regions()
    places = list(residues) if residues else \
        list(range(1, len(domain.positions) + 1))
    at, spans = 0, {}
    for name, positions in regions:
        spans[name] = (at, at + len(positions))
        at += len(positions)
    found = []
    for name in cdr.CDR_NAMES:
        start, end = spans[name]
        first = max(0, start - FLANK)
        last = min(len(places), end + FLANK)
        letters = ''.join(
            letter for _, letter in domain.positions[first:last])
        # (the places must stay in step with the letters)
        found.append((name, letters[:MAX_LENGTH],
                      places[first:last][:MAX_LENGTH]))
    return found


def parapred(sequences):
    '''
    Parapred's score for every residue of the sequences given

    :param sequences: [str], each an extended CDR
    :return: [[float]] in the same order
    '''
    import torch
    from parapred.model import Parapred, clean_output
    from parapred.cnn import generate_mask
    from parapred.preprocessing import encode_batch
    kept = [s for s in sequences if s]
    if not kept:
        return [[] for _ in sequences]
    # the model wants them longest first
    order = sorted(range(len(kept)), key=lambda i: -len(kept[i]))
    encoded, lengths = encode_batch([kept[i] for i in order],
                                    max_length=MAX_LENGTH)
    mask = generate_mask(encoded, lengths)
    model = Parapred()
    model.load_state_dict(torch.load(_parapred_weights(),
                                     map_location='cpu'))
    model.eval()
    with torch.no_grad():
        out = model(encoded, mask, lengths)
    scores = [None] * len(kept)
    for k, i in enumerate(order):
        scores[i] = [float(v) for v in clean_output(out[k], lengths[k])]
    found, at = [], 0
    for one in sequences:
        if one:
            found.append(scores[at])
            at += 1
        else:
            found.append([])
    return found


def predict_sequence(sequence, name='sequence', scheme='chothia'):
    '''
    Parapred on a sequence: its CDRs are found first (ANARCI)

    :return: [Residue] (resi: the 1-based place in the sequence)
    '''
    found = []
    for domain in cdr.domains(sequence, scheme):
        if not domain.is_antibody:
            continue
        blocks = extended_cdrs(domain)
        scores = parapred([letters for _, letters, _ in blocks])
        for (region, letters, places), values in zip(blocks, scores):
            for letter, place, score in zip(letters, places, values):
                found.append(Residue(name, domain.start + place, letter,
                                     score, region))
    return found


def predict_structure(obj, chain='', scheme='chothia', state=1, _self=cmd):
    '''
    Parapred on a structure: the same, mapped onto its residues

    :return: [Residue]
    '''
    found = []
    for chain_id, domain, residues in cdr.annotate(obj, chain, scheme,
                                                   state, _self):
        if not domain.is_antibody:
            continue
        blocks = extended_cdrs(domain, residues)
        scores = parapred([letters for _, letters, _ in blocks])
        for (region, letters, places), values in zip(blocks, scores):
            for letter, place, score in zip(letters, places, values):
                found.append(Residue(place[0], place[1], letter, score,
                                     region))
    return found


def _imgt_pdb(obj, pairs, folder, state, _self):
    """
    Write each Fv of an antibody as its own IMGT numbered PDB, as
    Paragraph wants them: its chains named H and L

    :param pairs: [(heavy chain, light chain or '')]
    :return: {(file name, written chain, IMGT number): (chain, resi)}
    """
    mapping = {}
    for k, (heavy, light) in enumerate(pairs):
        stem = 'fv_%d' % k
        name = _self.get_unused_name('_paratope')
        try:
            selections, here = [], {}
            for chain_id, letter in ((heavy, 'H'), (light, 'L')):
                if not chain_id:
                    continue
                for one, domain, residues in cdr.annotate(
                        obj, chain_id, 'imgt', state, _self):
                    if not domain.is_antibody:
                        continue
                    for (position, _letter), (_, resi) in zip(
                            domain.positions, residues):
                        number = domain.label(position)
                        here[(letter, number)] = (one, resi)
                        selections.append('(chain "%s" and resi %s)'
                                          % (one, resi))
            if not selections:
                continue
            _self.create(name, '(%%%s) and polymer.protein and (%s)'
                         % (obj, ' or '.join(selections)), state, 1)
            # renumber in one pass: residue by residue, an IMGT number
            # can collide with an original one not yet renumbered, and
            # the two residues then merge
            renumber = {}
            for (letter, number), (one, resi) in here.items():
                value, insertion = re.match(r'(-?\d+)([A-Za-z]?)',
                                            number).groups()
                renumber[(one, resi)] = (letter, value + insertion)
            _self.alter(name, 'chain, resi = NEW.get((chain, resi), '
                              '(chain, resi))', space={'NEW': renumber})
            _self.sort(name)
            _self.save(os.path.join(folder, stem + '.pdb'), '%' + name, 1)
            for key, where in here.items():
                mapping[(stem,) + key] = where
        finally:
            _self.delete(name)
    return mapping


def paragraph(obj, pairs=(), state=1, _self=cmd):
    """
    Paragraph on a structure: each Fv is written IMGT numbered, scored,
    and the scores mapped back onto the structure's residues

    :param pairs: [(heavy chain, light chain)]; found if not given
    :return: [Residue]
    """
    import subprocess
    import sys
    pairs = list(pairs) or antibody_pairs(obj, state, _self)
    if not pairs:
        raise ValueError('no antibody chain found in "%s"' % obj)
    found = []
    with tempfile.TemporaryDirectory() as folder:
        mapping = _imgt_pdb(obj, pairs, folder, state, _self)
        if not mapping:
            raise ValueError('no variable domain to score in "%s"' % obj)
        # Paragraph takes one run for paired Fvs and one for lone heavy
        # chains (--heavy picks the weights for them)
        groups = {True: [], False: []}
        for k, (_, light) in enumerate(pairs):
            if any(key[0] == 'fv_%d' % k for key in mapping):
                groups[bool(light)].append('fv_%d' % k)
        for paired, stems in groups.items():
            if not stems:
                continue
            found += _run_paragraph(folder, stems, paired, mapping)
    return found


def _run_paragraph(folder, stems, paired, mapping):
    import subprocess
    import sys
    key = os.path.join(folder, 'key_%s.csv' % ('paired' if paired
                                               else 'heavy'))
    with open(key, 'w') as handle:
        for stem in stems:
            handle.write('%s,H,%s\n' % (stem, 'L' if paired else ''))
    out = os.path.join(folder, 'out_%s.csv' % ('paired' if paired
                                               else 'heavy'))
    command = [os.path.join(os.path.dirname(sys.executable), 'Paragraph'),
               '--pdb_H_L_csv', key, '--pdb_folder_path', folder,
               '--out_path', out]
    if not paired:
        command.append('--heavy')
    done = subprocess.run(command, capture_output=True, text=True)
    if not os.path.exists(out):
        raise ValueError('Paragraph failed: %s'
                         % (done.stderr or done.stdout or '?')[-300:])
    found = []
    with open(out) as handle:
        for row in csv.DictReader(handle):
            where = mapping.get((row['pdb'], row['chain_id'], row['IMGT']))
            if where is None:
                continue
            found.append(Residue(where[0], where[1],
                                 _one_letter(row['AA']), row['pred']))
    return found


def parasurf_inputs(obj, pairs=(), state=1, _self=cmd):
    """
    Write each Fv as its own PDB, as ParaSurf reads them: the chains
    named H and L, no hydrogens, no het atoms, one conformer, the
    residues numbered as the structure has them

    :return: (folder, [(file stem, "paired"/"heavy"/"light")],
        {(stem, written chain, resi): (chain, resi, region)}); the
        caller removes the folder
    """
    pairs = list(pairs) or antibody_pairs(obj, state, _self)
    if not pairs:
        raise ValueError('no antibody chain found in "%s"' % obj)
    folder = tempfile.mkdtemp(prefix='parasurf_')
    jobs, mapping = [], {}
    for k, (heavy, light) in enumerate(pairs):
        stem = 'fv_%d' % k
        name = _self.get_unused_name('_paratope')
        try:
            selections, renames, here = [], {}, {}
            for chain_id, letter in ((heavy, 'H'), (light, 'L')):
                if not chain_id:
                    continue
                for one, domain, residues in cdr.annotate(
                        obj, chain_id, 'imgt', state, _self):
                    if not domain.is_antibody:
                        continue
                    renames[one] = letter
                    regions = {}
                    for region, places in cdr.region_residues(domain,
                                                              residues):
                        for place in places:
                            regions[place] = region
                    for _, resi in residues:
                        here[(letter, str(resi))] = (
                            one, str(resi), regions.get((one, resi), ''))
                        selections.append('(chain "%s" and resi %s)'
                                          % (one, resi))
            if not selections:
                continue
            _self.create(name, '(%%%s) and polymer.protein and not hydro '
                               'and (alt ""+A) and (%s)'
                         % (obj, ' or '.join(selections)), state, 1)
            # one pass, so two chains can swap their names
            _self.alter(name, 'chain = NEW.get(chain, chain)',
                        space={'NEW': renames})
            _self.sort(name)
            _self.save(os.path.join(folder, stem + '.pdb'), '%' + name, 1)
        finally:
            _self.delete(name)
        kinds = {('H', 'L'): 'paired', ('H',): 'heavy', ('L',): 'light'}
        kind = kinds.get(tuple(sorted(set(renames.values()))), 'paired')
        jobs.append((stem, kind))
        for key, where in here.items():
            mapping[(stem,) + key] = where
    return folder, jobs, mapping


def parasurf_scores(folder, jobs, mesh_dense=PARASURF_MESH, device='',
                    progress=None):
    """
    Run ParaSurf on the files parasurf_inputs wrote. Makes no PyMOL
    calls, so it can run off the GUI thread, and takes minutes on a CPU

    :param progress: called with (file stem, number done, total)
    :return: [(stem, chain, resi, residue name, score)]
    """
    import numpy as np
    import torch
    home = parasurf_home()
    if not home:
        raise ValueError(_parasurf_missing())
    if home not in sys.path:
        sys.path.insert(0, home)
    dms = parasurf_dms()
    if dms:          # ParaSurf runs it through the shell
        os.environ['PATH'] = (os.path.dirname(dms) + os.pathsep +
                              os.environ.get('PATH', ''))
    from ParaSurf.preprocess.clean_dataset import clean_dataset
    from ParaSurf.train.protein import Protein_pred
    from ParaSurf.train.network import Network
    from ParaSurf.train.utils import receptor_info
    if not device:
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    # hydrogens and het atoms out (already), atom numbers back in order
    clean_dataset(folder)
    found = []
    for done, (stem, kind) in enumerate(jobs):
        if progress is not None:
            progress(stem, done, len(jobs))
        path = os.path.join(folder, stem + '.pdb')
        protein = Protein_pred(path, save_path=folder,
                               mesh_dense=float(mesh_dense))
        surface = [f for f in os.listdir(protein.save_path)
                   if 'surfpoints' in f]
        if not surface:
            raise ValueError('ParaSurf: DMS wrote no surface for ' + stem)
        atoms = []
        with open(os.path.join(protein.save_path, surface[0])) as handle:
            for place, line in enumerate(handle):
                parts = line.split()
                if len(parts) > 6 and parts[6] == 'A':
                    atoms.append(place)       # a receptor atom, not a
                                              # point of its surface
        network = Network(parasurf_weights(kind), gridSize=PARASURF_GRID,
                          feature_channels=PARASURF_CHANNELS,
                          device=device)
        scores = network.get_lig_scores(
            protein, batch_size=PARASURF_BATCH, add_forcefields=True,
            add_atom_radius_features=True)
        best = receptor_info(path, np.array([scores[i] for i in atoms]))[1]
        names = _pdb_residue_names(path)
        for res_id, data in best.items():
            parts = res_id.split('_')         # number_chain[_insertion]
            if len(parts) < 2:
                continue
            resi = parts[0] + (parts[2] if len(parts) > 2 else '')
            found.append((stem, parts[1], resi, names.get(res_id, ''),
                          float(data['scores'])))
    if progress is not None:
        progress('', len(jobs), len(jobs))
    return found


def parasurf_residues(rows, mapping):
    """
    The scores of parasurf_scores on the structure's own chains

    :return: [Residue]
    """
    found = []
    for stem, chain, resi, resn, score in rows:
        where = mapping.get((stem, chain, resi))
        if where is None:
            continue
        found.append(Residue(where[0], where[1], _one_letter(resn), score,
                             where[2] if len(where) > 2 else ''))
    return found


def _pdb_residue_names(path):
    """
    :return: {ParaSurf's residue id: residue name} of a PDB file
    """
    names = {}
    with open(path) as handle:
        for line in handle:
            if not line.startswith('ATOM'):
                continue
            res_id = '%s_%s' % (line[22:26].strip(), line[21])
            if line[26].strip():
                res_id += '_' + line[26].strip()
            names.setdefault(res_id, line[17:20].strip())
    return names


def parasurf(obj, pairs=(), mesh_dense=PARASURF_MESH, device='', state=1,
             _self=cmd):
    """
    ParaSurf on a structure: each Fv is written out, its surface
    computed and scored, and the scores mapped back onto the residues

    :return: [Residue]
    """
    folder, jobs, mapping = parasurf_inputs(obj, pairs, state, _self)
    try:
        if not jobs:
            raise ValueError('no variable domain to score in "%s"' % obj)
        rows = parasurf_scores(folder, jobs, mesh_dense, device)
    finally:
        shutil.rmtree(folder, ignore_errors=True)
    return parasurf_residues(rows, mapping)


def _one_letter(resn):
    from . import seqcompare
    return seqcompare.ONE_LETTER.get(str(resn).upper(), 'X')


def antibody_pairs(obj, state=1, _self=cmd):
    '''
    The antibody chains of a structure, paired heavy with light

    :return: [(heavy chain, light chain or '')]
    '''
    heavy, light = [], []
    for chain_id, domain, _ in cdr.annotate(obj, '', 'imgt', state, _self):
        if not domain.is_antibody:
            continue
        (light if domain.kind == 'L' else heavy).append(chain_id)
    pairs = []
    for k, one in enumerate(heavy):
        pairs.append((one, light[k] if k < len(light) else ''))
    for extra in light[len(heavy):]:
        pairs.append(('', extra))
    return pairs


def predict(obj_or_sequence, method=DEFAULT_METHOD, chain='',
            name='sequence', scheme='chothia', state=1, _self=cmd):
    '''
    Paratope residues of an antibody

    :param obj_or_sequence: an object in the session, or a sequence
    :return: [Residue]
    '''
    method = str(method).strip().lower()
    missing = available(method)
    if missing:
        raise ValueError(missing)
    text = str(obj_or_sequence).strip()
    is_object = text in (_self.get_names('objects') or [])
    if method in STRUCTURE_METHODS:
        if not is_object:
            raise ValueError('%s needs a structure: give an object'
                             % method_label(method).split(' (')[0])
        if method == 'parasurf':
            return parasurf(text, (), PARASURF_MESH, '', state, _self)
        return paragraph(text, (), state, _self)
    if is_object:
        return predict_structure(text, chain, scheme, state, _self)
    return predict_sequence(text, name, scheme)


def apply_scores(obj, residues, color=True, select=True, b_factors=True,
                 _self=cmd):
    '''
    Show the scores on a structure: in the B-factors, as a selection of
    the residues above THRESHOLD, and coloured by score
    '''
    by_residue = {(r.chain, r.resi): r.score for r in residues}
    if b_factors:
        _self.alter('%' + obj, 'b = SCORES.get((chain, resi), 0.0)',
                    space={'SCORES': by_residue})
        _self.sort(obj)
    binding = [r for r in residues if r.binds]
    if select and binding:
        _self.select(_self.get_legal_name('%s_paratope' % obj),
                     '(%%%s) and (%s)' % (obj, _residues_selection(
                         [(r.chain, r.resi) for r in binding])), enable=0)
    if color:
        _self.color(REST_COLOR, '%' + obj)
        if binding:
            _self.color(PARATOPE_COLOR,
                        '(%%%s) and (%s)' % (obj, _residues_selection(
                            [(r.chain, r.resi) for r in binding])))
    return len(binding)


def _residues_selection(keys):
    by_chain = {}
    for chain, resi in keys:
        by_chain.setdefault(chain, []).append(str(resi))
    return ' or '.join('(chain "%s" and resi %s)' % (chain, '+'.join(res))
                       for chain, res in sorted(by_chain.items())) or 'none'


def report(residues, method=DEFAULT_METHOD, limit=0):
    '''
    The scored residues as a table, the likeliest binders first
    '''
    lines = [' Paratope (%s), %d residue(s) scored, %d at or above %.2f'
             % (method_label(method), len(residues),
                sum(1 for r in residues if r.binds), THRESHOLD)]
    if not residues:
        return lines[0]
    lines.append(' %-6s %-7s %-5s %-6s %s' % ('Chain', 'Residue', 'Score',
                                              'Binds', 'Region'))
    ranked = sorted(residues, key=lambda r: -r.score)
    for r in ranked[:limit] if limit else ranked:
        lines.append(' %-6s %-7s %-5.2f %-6s %s' % (
            r.chain, '%s%s' % (r.letter, r.resi), r.score,
            'yes' if r.binds else '', r.region))
    return '\n'.join(lines)


def predict_paratope(selection, method=DEFAULT_METHOD, chain='',
                     color=1, select=1, quiet=0, state=1, *, _self=cmd):
    '''
DESCRIPTION

    "predict_paratope" marks the residues of an antibody that are
    likely to bind its antigen, with Parapred (sequence), Paragraph
    (structure) or ParaSurf (structure).

USAGE

    predict_paratope object [, method [, chain ]]

ARGUMENTS

    object = str: the structure (or a sequence for Parapred)

    method = parapred, paragraph or parasurf {default: parapred}

    chain = str: one chain {default: every antibody chain}

    color = 0/1: colour the predicted paratope red, the rest grey, and
    put the scores in the B-factors {default: 1}

    select = 0/1: a selection "<object>_paratope" of the residues above
    the threshold {default: 1}

NOTES

    Parapred and Paragraph do not score the whole chain: a residue with
    no score is one they did not look at, not one they judged not to
    bind. ParaSurf scores every residue of the Fv, but takes minutes
    without a GPU.

EXAMPLE

    predict_paratope 1hzh, method=paragraph

SEE ALSO

    find_cdrs, epitope, interface_analysis
    '''
    method = str(method).strip().lower()
    residues = predict(selection, method, chain, 'sequence', 'chothia',
                       int(state), _self)
    if str(selection).strip() in (_self.get_names('objects') or []):
        apply_scores(str(selection).strip(), residues, bool(int(color)),
                     bool(int(select)), bool(int(color)), _self)
    if not int(quiet):
        print(report(residues, method))
    return residues


def extend(_self):
    _self.extend('predict_paratope', predict_paratope)
    _self.predict_paratope = predict_paratope
