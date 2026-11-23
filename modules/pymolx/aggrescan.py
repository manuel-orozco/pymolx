'''
Aggregation prediction (parity item L-34): which residues of a folded
protein drive its self-association, and how that changes with pH.

This is a free re-implementation of Aggrescan4D (A4D, Barcenas et al.,
NAR 2024; Zalewski et al., Protein Science 2024) built on the published
Aggrescan3D standalone code (A3D, MIT licensed, bitbucket.org/lcbio/
aggrescan3d). A4D itself needs FoldX, which is not free; everything
here runs on open source tools:

- the score and its constants are ported from A3D's aggrescan_3d.py,
- solvent exposure comes from freesasa with NACCESS radii, the same
  program and radii A3D ships,
- structure optimisation and mutant building use OpenMM and PDBFixer
  where A4D runs FoldX RepairPDB and BuildModel (pymolx.mm),
- the stability effect of a mutation comes from EvoEF, a free force
  field in FoldX's mould which takes the same mutant-file syntax. The
  difference of two OpenMM potential energies will not do: a mutant
  and its wild type differ in composition, and that swamps the
  stability effect,
- pH shifted pKa values come from PROPKA instead of pKa-ANI,
- dynamic mode samples with OpenMM implicit solvent molecular dynamics
  instead of CABS-flex, and clusters the frames (pymolx.cluster),
- the globular regions option uses the pLDDT a predicted model carries
  in its B-factors, where A4D runs AlphaCutter.

Only freesasa is needed to score a structure; the rest are needed for
the options that use them, and are reported as missing if they are.

Two deliberate differences from the A4D papers, which the published A3D
code does not implement: tyrosine is kept constant with pH (A4D varies
it, and its parameters are not published), and arginine varies with pH
as A3D's code has it (A4D freezes it; set FREEZE_ARGININE to follow the
paper instead).

The score of a residue is its intrinsic aggregation propensity, damped
by how buried it is, plus the same quantity summed over its neighbours
within a sphere and decayed with distance. Above 0 a residue is
aggregation prone; above 0.5 it is likely to matter.
'''
import math
import os
import shutil
import sys
import tempfile

from pymol import cmd

# --- the A3D scoring constants (aggrescan_3d.py) ----------------------

# the in vivo derived Aggrescan propensity scale (data/matrices/
# aggrescan.mat), Conchillo-Sole et al. 2007
A3D_SCALE = {
    'I': 1.9136, 'F': 1.8456, 'V': 1.6856, 'L': 1.4716, 'Y': 1.2506,
    'W': 1.1286, 'M': 1.0016, 'C': 0.6956, 'A': 0.0556, 'T': -0.0674,
    'S': -0.2024, 'P': -0.2424, 'G': -0.4434, 'K': -1.6164, 'H': -0.9414,
    'Q': -1.1394, 'R': -1.7534, 'N': -1.2104, 'E': -1.7304, 'D': -1.7014,
}
# relative solvent exposure (%) below which a residue counts as buried
# and scores 0, and above which the exposure term stops growing
MIN_SURFACE = 10.0
MAX_SURFACE = 55.0
# A3D reads exposure back from freesasa's .rsa text, which carries one
# decimal, and a residue rounding to 10.0 counts as exposed. Rounding
# the same way reproduces its scores; without it the residues sitting
# on the boundary differ
RSA_DECIMALS = 1
# exposure term: ALPHA * exp(BETA * rsa)
ALPHA = 0.0599
BETA = 0.0521
# neighbour term: GAMMA * exp(-DELTA / radius * distance)
GAMMA = 1.2915
DELTA = 2.56
# the sphere radius, in Angstrom, the A4D server offers
DISTANCES = (5, 10)
DEFAULT_DISTANCE = 10
# a residue at or above this score is aggregation prone; the second is
# the paper's "likely to matter" mark
APR_THRESHOLD = 0.0
HOTSPOT_THRESHOLD = 0.5

# --- the pH dependent lipophilicity scale (A3D, SolupHred) ------------

# residues whose propensity does not change with pH, on the scale the
# pH calculation uses (which is not the scale above: see the papers'
# Figure 2a, the two correlate but their ranges differ)
PH_FIXED_SCORES = {
    'A': 0.66, 'N': -0.84, 'C': 1.66, 'Q': -0.87, 'G': 0.0, 'I': 2.75,
    'L': 1.77, 'M': 1.30, 'F': 3.99, 'P': 1.69, 'S': -0.99, 'T': 0.12,
    'W': 3.29, 'Y': 1.33, 'V': 1.45,
}
# the ionizable ones: pKa, a weight on the neutral term, and the
# partition coefficients of the neutral and ionized forms
PH_IONIZABLE = {
    'R': {'pka': 12.5, 'weight': 0.857546950885663,
          'neutral': 10 ** -3.66, 'ionized': 10 ** -7.38, 'acid': False},
    'D': {'pka': 3.5, 'weight': 0.622376926398327,
          'neutral': 10 ** -3.18, 'ionized': 10 ** -8.54, 'acid': True},
    'E': {'pka': 4.2, 'weight': 0.940369170156725,
          'neutral': 10 ** -3.79, 'ionized': 10 ** -6.2, 'acid': True},
    'H': {'pka': 6.6, 'weight': 0.812212007706225,
          'neutral': 10 ** -4.67, 'ionized': 10 ** -5.97, 'acid': False},
    'K': {'pka': 10.5, 'weight': 0.9093172437,
          'neutral': 10 ** -2.19, 'ionized': 10 ** -6.81, 'acid': False},
}
# puts the pH scale on the same footing as the one above
PH_OFFSET = -3.13
# the pH values the A4D server reports: 4.0 to 9.0 in steps of 0.5
PH_VALUES = tuple(round(4.0 + 0.5 * i, 1) for i in range(11))
DEFAULT_PH = 7.0
# where the pKa of an ionizable residue comes from
PKA_MODES = (('3D pKa (PROPKA)', 'propka'), ('Standard pKa', 'fixed'))
DEFAULT_PKA = 'propka'
# A4D holds arginine constant (its 3D pKa cannot be predicted well);
# the published A3D code lets it vary, which is what this follows
FREEZE_ARGININE = False
# a predicted pKa this far from the pH saturates the formula anyway
MAX_PKA_GAP = 12.0

# --- which metric describes a structure (NAR 2024, SI 1) --------------

# a residue is "packed" with this many alpha carbons within this radius,
# and a structure is globular when this fraction of them are
COMPACT_NEIGHBOURS = 10
COMPACT_RADIUS = 10.0
COMPACT_FRACTION = 0.6
METRICS = (('Maximum score', 'max'), ('Average score', 'avg'),
           ('Consensus', 'consensus'))
DEFAULT_METRIC = 'consensus'

# --- automated mutations ----------------------------------------------

# A3D mutates aggregation prone residues to these; the first four are
# the charged set the A4D server offers
CHARGED_MUTATIONS = ('R', 'K', 'D', 'E')
# a residue needs at least this score to be worth mutating (A3D)
MUTATION_THRESHOLD = -0.2
# how many positions, and how many variants of each, the server keeps
MAX_SITES = 6
MAX_PER_SITE = 2
MAX_SUGGESTIONS = 12
# the evolutionary mode keeps substitutions scoring at least this in
# PAM250 and suggests at most this many per position (NAR 2024)
PAM_THRESHOLD = 0
MAX_PER_SITE_EVOLUTIONARY = 3
MUTATION_MODES = (('Charged mutations', 'charged'),
                  ('Evolutionary conserved mutations', 'evolutionary'))

# --- display -----------------------------------------------------------

# the A4D colours: blue soluble, white neutral, red aggregation prone.
# The ends are fixed so two structures can be compared by eye
COLOR_PALETTE = 'blue_white_red'
COLOR_RANGE = 2.0
# below this pLDDT a residue of a predicted model is left out when the
# globular regions option is on
PLDDT_CUTOFF = 70.0
# how see-through the "aggregation" preset's surface is
PRESET_TRANSPARENCY = 0.5

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    # the ones A3D also reads
    'MSE': 'M', 'HID': 'H', 'HIE': 'H', 'HIP': 'H', 'CYX': 'C',
    'ASH': 'D', 'GLH': 'E', 'LYN': 'K',
}
ONE_TO_THREE = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
                'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY',
                'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
                'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
                'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}
AA_NAMES = {'A': 'alanine', 'R': 'arginine', 'N': 'asparagine',
            'D': 'aspartic acid', 'C': 'cysteine', 'E': 'glutamic acid',
            'Q': 'glutamine', 'G': 'glycine', 'H': 'histidine',
            'I': 'isoleucine', 'L': 'leucine', 'K': 'lysine',
            'M': 'methionine', 'F': 'phenylalanine', 'P': 'proline',
            'S': 'serine', 'T': 'threonine', 'W': 'tryptophan',
            'Y': 'tyrosine', 'V': 'valine'}


def aggrescan_home():
    '''
    The folder holding the Aggrescan3D checkout this is ported from and
    its data: $AGGRESCAN_HOME, or an "aggrescan" folder beside the
    Python environment, beside its checkout, or in the home directory

    :return: its path, or ''
    '''
    places = [os.environ.get('AGGRESCAN_HOME', '')]
    beside = os.path.dirname(sys.prefix)
    for root in (beside, os.path.dirname(beside), os.path.expanduser('~')):
        places += [os.path.join(root, 'aggrescan'),
                   os.path.join(root, 'Aggrescan')]
    for place in places:
        if place and os.path.isdir(place):
            return place
    return ''


def available(ph=False, optimize=False):
    '''
    :param ph: the pH calculation is wanted (PROPKA)
    :param optimize: the structure is to be minimized (OpenMM)
    :return: '' if it can be run, else what is missing
    '''
    try:
        import freesasa  # noqa: F401
    except ImportError:
        return ('freesasa is not installed: "pip install freesasa" in '
                'the environment PyMOL runs in')
    if optimize:
        from . import mm
        if not mm.available():
            return ('Structure optimisation needs OpenMM and PDBFixer: '
                    'pip install openmm pdbfixer')
    if ph and not propka_available():
        return ('PROPKA is not installed: "pip install propka", or '
                'choose the standard pKa values')
    return ''


def propka_available():
    try:
        import propka.run  # noqa: F401
        return True
    except ImportError:
        return False


class Residue:
    '''
    One scored residue

    :ivar chain: its chain ("-" when the file has none, as in A3D)
    :ivar resi: its number, with an insertion code if it has one
    :ivar letter: one-letter code
    :ivar exposure: relative solvent accessibility, in percent
    :ivar intrinsic: its own propensity after the exposure term
    :ivar score: the A4D score, the one to report
    :ivar coords: its alpha carbon
    '''

    def __init__(self, chain, resi, letter, exposure, intrinsic, coords):
        self.chain = chain
        self.resi = str(resi)
        self.letter = letter
        self.exposure = float(exposure)
        self.intrinsic = float(intrinsic)
        self.score = 0.0
        self.coords = coords

    @property
    def key(self):
        return (self.chain, self.resi)

    @property
    def buried(self):
        return self.intrinsic == 0.0

    @property
    def prone(self):
        return self.score > APR_THRESHOLD

    @property
    def hotspot(self):
        return self.score >= HOTSPOT_THRESHOLD

    def __repr__(self):
        return '<%s %s%s %.4f>' % (self.letter, self.chain, self.resi,
                                   self.score)


class Result:
    '''
    The scores of one structure

    :ivar residues: [Residue] in the order they appear in the file
    :ivar distance: the sphere radius used
    :ivar ph: the pH it was scored at, or None for the classic scale
    :ivar name: what it was scored from
    '''

    def __init__(self, residues, distance=DEFAULT_DISTANCE, ph=None,
                 name=''):
        self.residues = list(residues)
        self.distance = distance
        self.ph = ph
        self.name = name

    def __len__(self):
        return len(self.residues)

    @property
    def scores(self):
        return [r.score for r in self.residues]

    @property
    def total(self):
        return sum(self.scores)

    @property
    def average(self):
        return self.total / len(self.residues) if self.residues else 0.0

    @property
    def maximum(self):
        return max(self.scores) if self.residues else 0.0

    @property
    def minimum(self):
        return min(self.scores) if self.residues else 0.0

    @property
    def prone(self):
        return [r for r in self.residues if r.prone]

    def metric(self, which=DEFAULT_METRIC):
        '''
        The single number that describes the structure

        :param which: "max", "avg", or "consensus" to pick by shape
        :return: (value, which was used)
        '''
        if which == 'consensus':
            which = 'max' if self.is_globular() else 'avg'
        return (self.maximum if which == 'max' else self.average), which

    def is_globular(self):
        '''
        The NAR 2024 compactness criterion: a structure is globular when
        more than COMPACT_FRACTION of its residues have at least
        COMPACT_NEIGHBOURS alpha carbons within COMPACT_RADIUS
        '''
        import numpy
        if len(self.residues) <= COMPACT_NEIGHBOURS:
            return False
        xyz = numpy.array([r.coords for r in self.residues])
        packed = 0
        for one in xyz:
            close = numpy.linalg.norm(xyz - one, axis=1) <= COMPACT_RADIUS
            # the residue itself is in there, so look for one more
            if close.sum() - 1 >= COMPACT_NEIGHBOURS:
                packed += 1
        return packed / len(self.residues) > COMPACT_FRACTION

    def by_chain(self):
        '''
        :return: {chain: [Residue]} in file order
        '''
        found = {}
        for residue in self.residues:
            found.setdefault(residue.chain, []).append(residue)
        return found

    def stats(self):
        '''
        The summary A3D writes as A3D_summary.json

        :return: {chain or "All": {max_value, avg_value, ...}}
        '''
        def one(residues):
            scores = [r.score for r in residues] or [0.0]
            return {'max_value': round(max(scores), 4),
                    'avg_value': round(sum(scores) / len(scores), 4),
                    'min_value': round(min(scores), 4),
                    'total_value': round(sum(scores), 4)}
        found = {chain: one(residues)
                 for chain, residues in self.by_chain().items()}
        found['All'] = one(self.residues)
        return found

    def csv(self):
        '''
        The per residue table, in A3D's A3D.csv column order
        '''
        lines = ['protein,chain,residue,residue_name,score']
        for r in self.residues:
            lines.append('%s,%s,%s,%s,%.4f' % (self.name, r.chain, r.resi,
                                               r.letter, round(r.score, 4)))
        return '\n'.join(lines) + '\n'


def ph_score(letter, ph, pka=None, signed=False):
    '''
    The pH dependent propensity of one residue

    :param letter: its one-letter code
    :param ph: the pH
    :param pka: its pKa, or None for the standard value
    :param signed: measure the pH from the pKa in the direction that
        charges the residue, rather than A3D's absolute difference.
        A3D's form is right while a residue is charged, which holds for
        the standard pKa values over the server's pH range, but not
        necessarily for a predicted one
    :return: the propensity, on the pH scale
    '''
    if letter in PH_FIXED_SCORES:
        return PH_FIXED_SCORES[letter]
    data = PH_IONIZABLE.get(letter)
    if data is None:
        return 0.0
    if letter == 'R' and FREEZE_ARGININE:
        ph = DEFAULT_PH
        pka = None
    if pka is None:
        pka = data['pka']
    if signed:
        # positive while the residue carries its charge
        gap = (ph - pka) if data['acid'] else (pka - ph)
    else:
        gap = abs(pka - ph)
    gap = max(-MAX_PKA_GAP, min(MAX_PKA_GAP, gap))
    ratio = 10.0 ** gap
    neutral = math.log10(data['neutral'] + data['ionized'] * ratio)
    both = math.log10(1.0 + ratio)
    return neutral * data['weight'] - both - PH_OFFSET


def intrinsic_score(letter, exposure, ph=None, pka=None, signed=False):
    '''
    A residue's own contribution: its propensity damped by how buried
    it is. Zero below MIN_SURFACE, and the exposure term stops growing
    above MAX_SURFACE
    '''
    if exposure < MIN_SURFACE:
        return 0.0
    exposure = min(exposure, MAX_SURFACE)
    if ph is None:
        propensity = A3D_SCALE.get(letter)
        if propensity is None:
            return 0.0
    else:
        propensity = ph_score(letter, ph, pka, signed)
    return propensity * ALPHA * math.exp(BETA * exposure)


def spread(residues, distance=DEFAULT_DISTANCE):
    '''
    Add each residue's neighbours to its score, decayed with distance.
    A buried residue keeps a score of 0, as in A3D

    :param residues: [Residue] with their intrinsic scores set
    :param distance: the sphere radius, in Angstrom
    '''
    import numpy
    if not residues:
        return
    xyz = numpy.array([r.coords for r in residues], dtype=float)
    own = numpy.array([r.intrinsic for r in residues], dtype=float)
    decay = -DELTA / float(distance)
    for index, residue in enumerate(residues):
        if residue.intrinsic == 0.0:
            residue.score = 0.0
            continue
        gaps = numpy.linalg.norm(xyz - xyz[index], axis=1)
        near = (gaps < distance) & (gaps != 0.0)
        residue.score = float(
            residue.intrinsic
            + (own[near] * GAMMA * numpy.exp(decay * gaps[near])).sum())


def exposures(pdb_path):
    '''
    Relative solvent accessibility of every residue of a PDB file, the
    way A3D measures it: freesasa, NACCESS radii, Lee-Richards with 100
    slices

    :return: {(chain, resi): percent}
    '''
    import freesasa
    classifier = freesasa.Classifier.getStandardClassifier('naccess')
    structure = freesasa.Structure(pdb_path, classifier)
    result = freesasa.calc(structure, freesasa.Parameters(
        {'algorithm': freesasa.LeeRichards, 'n-slices': 100}))
    found = {}
    for chain, residues in result.residueAreas().items():
        for resi, area in residues.items():
            found[(chain.strip() or '-', resi.strip())] = round(
                area.relativeTotal * 100.0, RSA_DECIMALS)
    return found


def alpha_carbons(pdb_path):
    '''
    The alpha carbons of a PDB file, in file order, the way A3D reads
    them (first model, ATOM records, no alternates but the first)

    :return: [(chain, resi, three-letter name, (x, y, z))]
    '''
    found, seen = [], set()
    with open(pdb_path) as handle:
        for line in handle:
            if line.startswith('ENDMDL'):
                break
            if not line.startswith('ATOM'):
                continue
            if line[12:16].strip() != 'CA':
                continue
            if line[16] not in (' ', 'A'):
                continue
            chain = line[21].strip() or '-'
            resi = (line[22:26].strip() + line[26].strip())
            if (chain, resi) in seen:
                continue
            seen.add((chain, resi))
            found.append((chain, resi, line[17:20].strip(),
                          (float(line[30:38]), float(line[38:46]),
                           float(line[46:54]))))
    return found


def score_pdb(pdb_path, distance=DEFAULT_DISTANCE, ph=None, pkas=None,
              name='', signed=False):
    '''
    Score a PDB file. Makes no PyMOL calls, so a worker thread can run it

    :param pdb_path: a file of protein ATOM records without hydrogens
    :param distance: the sphere radius, 5 or 10 Angstrom
    :param ph: a pH to use the pH dependent scale at, or None for the
        classic Aggrescan scale
    :param pkas: {(chain, resi): pKa} from PROPKA, or None
    :return: Result
    '''
    rsa = exposures(pdb_path)
    residues = []
    for chain, resi, resn, coords in alpha_carbons(pdb_path):
        letter = THREE_TO_ONE.get(resn.upper())
        if letter is None:
            continue
        exposure = rsa.get((chain, resi))
        if exposure is None:
            continue
        pka = (pkas or {}).get((chain, resi))
        residues.append(Residue(
            chain, resi, letter, exposure,
            intrinsic_score(letter, exposure, ph, pka, signed), coords))
    spread(residues, distance)
    return Result(residues, distance, ph,
                  name or os.path.basename(pdb_path).rsplit('.', 1)[0])


# --- reading a structure out of PyMOL ---------------------------------

# what A3D analyses: protein, no hydrogens, first alternate only, and
# nothing pymolx.missing added without coordinates
SELECTION = 'polymer.protein and present and not hydro and (alt ""+A)'


def write_pdb(selection, path, state=1, _self=cmd):
    '''
    Write the part of a selection the analysis reads

    :return: how many atoms were written
    '''
    picked = '(%s) and %s' % (selection, SELECTION)
    count = _self.count_atoms(picked, state=state)
    if not count:
        raise ValueError('no protein atoms in "%s"' % selection)
    _self.save(path, picked, state)
    return count


def pkas(pdb_path):
    '''
    PROPKA's pKa for every titratable residue of a PDB file. This takes
    the place of A4D's pKa-ANI

    :return: {(chain, resi): pKa}
    '''
    import logging

    import propka.run
    logger = logging.getLogger('propka')
    level = logger.level
    logger.setLevel(logging.CRITICAL)
    try:
        molecule = propka.run.single(pdb_path, optargs=['--quiet'],
                                     write_pka=False)
    finally:
        logger.setLevel(level)
    conformation = molecule.conformations[molecule.conformation_names[0]]
    found = {}
    for group in conformation.groups:
        # only the residues whose propensity moves with pH; the N+ and
        # C- terminal groups would collide with their residue's key
        if group.residue_type not in PH_IONIZABLE:
            continue
        atom = group.atom
        chain = (atom.chain_id or '').strip() or '-'
        resi = '%s%s' % (atom.res_num,
                         (getattr(atom, 'icode', '') or '').strip())
        found[(chain, resi)] = group.pka_value
    return found


def globular_selection(selection, cutoff=PLDDT_CUTOFF, _self=cmd):
    '''
    The confidently folded part of an AlphaFold model, which takes the
    place of A4D's "analysis of the globular regions". A4D runs
    AlphaCutter, which reads DSSP and a contact map; this uses the
    pLDDT the model carries in its B-factors, since the long loops
    AlphaCutter removes are the ones predicted with low confidence.

    Only useful for predicted models: an experimental structure's
    B-factors are not pLDDT, and are left alone

    :return: (selection, residues dropped, whether it looked like a
        predicted model)
    '''
    values = []
    _self.iterate('(%s) and name CA' % selection,
                  'rows.append((chain, resi, b))', space={'rows': values})
    if not values:
        return selection, 0, False
    scores = [b for _, _, b in values]
    # pLDDT runs 0 to 100 and a model is mostly confident
    looks_predicted = (min(scores) >= 0.0 and max(scores) <= 100.0
                       and max(scores) > 1.0
                       and sum(b > 50 for b in scores) > 0.5 * len(scores))
    if not looks_predicted:
        return selection, 0, False
    weak = [(chain, resi) for chain, resi, b in values if b < cutoff]
    if not weak or len(weak) >= len(values):
        return selection, 0, True
    return ('(%s) and not (%s)' % (selection, _residues_selection(
        [(c or '-', r) for c, r in weak])), len(weak), True)


def score(selection='polymer', distance=DEFAULT_DISTANCE, ph=None,
          pka=DEFAULT_PKA, state=1, name='', _self=cmd):
    '''
    Score a structure in the session

    :param selection: what to score
    :param distance: the sphere radius, 5 or 10 Angstrom
    :param ph: a pH, or None for the classic Aggrescan scale
    :param pka: "propka" for 3D pKa values, "fixed" for standard ones
    :return: Result
    '''
    folder = tempfile.mkdtemp(prefix='aggrescan_')
    try:
        path = os.path.join(folder, 'input.pdb')
        write_pdb(selection, path, state, _self)
        values = pkas(path) if (ph is not None and pka == 'propka'
                                and propka_available()) else None
        return score_pdb(path, distance, ph, values,
                         name or str(selection), signed=values is not None)
    finally:
        shutil.rmtree(folder, ignore_errors=True)


def scan_ph(pdb_path, distance=DEFAULT_DISTANCE, values=PH_VALUES,
            pka=DEFAULT_PKA, name='', progress=None):
    '''
    Score one structure at each of several pH values. The pKa values
    only depend on the structure, so they are predicted once

    :param progress: callable(done, total) for a worker thread
    :return: [(pH, Result)]
    '''
    predicted = pkas(pdb_path) if (pka == 'propka' and propka_available()) \
        else None
    found = []
    for index, one in enumerate(values):
        if progress is not None:
            progress(index, len(values))
        found.append((one, score_pdb(pdb_path, distance, one, predicted,
                                     name, signed=predicted is not None)))
    if progress is not None:
        progress(len(values), len(values))
    return found


# --- OpenMM: optimisation, dynamics and mutants -----------------------
#
# A4D uses FoldX to optimise the input and to build mutants, and
# CABS-flex to sample fluctuations. Neither is free, so this does both
# with OpenMM: a minimization in implicit solvent in place of FoldX's
# RepairPDB, PDBFixer's side chain building in place of BuildModel, and
# a short Langevin run clustered into models in place of CABS-flex.
# Energies are therefore OpenMM potential energies, not FoldX dG.

# kJ/mol per kcal/mol
KJ_PER_KCAL = 4.184
# how the A4D server marks a mutation's energetic effect, kcal/mol
ENERGY_WARN = 0.5
ENERGY_BAD = 1.0
# defaults for the molecular dynamics of the dynamic mode
MD_MODELS = 12
MD_TEMPERATURE = 300.0
MD_STEPS = 50000
MD_EQUILIBRATION = 5000
MD_TIMESTEP = 0.002
MD_FRAMES = 100


def _openmm():
    from . import mm
    return mm._require()


def _prepared(pdb_path, mutations=(), chain='', ph=DEFAULT_PH):
    '''
    Read a PDB file into OpenMM, optionally mutating it first

    :param mutations: ["ALA-57-GLY"] as PDBFixer spells them, or
        {chain: [...]} to mutate several chains at once
    :param chain: the chain they are in, for the plain list
    :return: (app, topology, positions, system)
    '''
    _, app, _, PDBFixer = _openmm()
    from . import mm
    with open(pdb_path) as handle:
        fixer = PDBFixer(pdbfile=handle)
    if mutations:
        if isinstance(mutations, dict):
            for one, items in mutations.items():
                fixer.applyMutations(list(items), one)
        else:
            fixer.applyMutations(list(mutations), chain)
    # the force field has no templates for anything else
    extra = [r for r in fixer.topology.residues()
             if r.name not in mm.STANDARD_RESIDUES]
    if extra:
        modeller = app.Modeller(fixer.topology, fixer.positions)
        modeller.delete(extra)
        fixer.topology, fixer.positions = (modeller.topology,
                                           modeller.positions)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    forcefield = app.ForceField('amber14-all.xml', 'implicit/gbn2.xml')
    system = forcefield.createSystem(
        fixer.topology, nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds, hydrogenMass=None)
    return app, fixer.topology, fixer.positions, system


def _platform(openmm, system, positions, integrator, platform='auto'):
    from . import mm
    names = mm.PLATFORMS if platform == 'auto' else (platform,)
    errors = []
    for name in names:
        try:
            chosen = openmm.Platform.getPlatformByName(name)
            context = openmm.Context(system, integrator, chosen)
            context.setPositions(positions)
            return context, name
        except Exception as e:
            errors.append('%s: %s' % (name, e))
    raise RuntimeError('no OpenMM platform worked: ' + '; '.join(errors))


def _save(app, topology, positions, path, heavy_only=True):
    '''
    Write a topology to a PDB file, by default without hydrogens so the
    exposure calculation reads it the way A3D's input is read
    '''
    if heavy_only:
        modeller = app.Modeller(topology, positions)
        modeller.delete([a for a in modeller.topology.atoms()
                         if a.element is not None
                         and a.element.symbol == 'H'])
        topology, positions = modeller.topology, modeller.positions
    with open(path, 'w') as handle:
        app.PDBFile.writeFile(topology, positions, handle, keepIds=True)


def _minimized(pdb_path, out_path='', mutations=(), chain='',
               fixed_outside=None, tolerance=10.0, max_iterations=0,
               platform='auto', ph=DEFAULT_PH):
    '''
    Build a structure (optionally mutating it) and minimize it

    :param fixed_outside: (chain, resi) to let only that residue move,
        everything else held where it is
    :param out_path: where to write the result, without hydrogens
    :return: {energy_before, energy_after (kcal/mol), platform, atoms}
    '''
    openmm, app, unit, _ = _openmm()
    app, topology, positions, system = _prepared(pdb_path, mutations, chain,
                                                 ph)
    if fixed_outside is not None:
        want = (str(fixed_outside[0]), str(fixed_outside[1]))
        moving = 0
        for atom in topology.atoms():
            residue = atom.residue
            key = (residue.chain.id.strip() or '-',
                   residue.id + residue.insertionCode.strip())
            if key == want:
                moving += 1
            else:
                system.setParticleMass(atom.index, 0.0)
        if not moving:
            raise ValueError('residue %s%s is not in the structure'
                             % fixed_outside)
    integrator = openmm.VerletIntegrator(0.001)
    context, used = _platform(openmm, system, positions, integrator, platform)

    def energy():
        return context.getState(getEnergy=True).getPotentialEnergy() \
            .value_in_unit(unit.kilocalorie_per_mole)

    before = energy()
    openmm.LocalEnergyMinimizer.minimize(
        context, tolerance * unit.kilojoule_per_mole / unit.nanometer,
        int(max_iterations))
    after = energy()
    if out_path:
        _save(app, topology,
              context.getState(getPositions=True).getPositions(), out_path)
    return {'energy_before': before, 'energy_after': after,
            'platform': used, 'atoms': topology.getNumAtoms()}


def optimize_pdb(pdb_path, out_path, tolerance=10.0, max_iterations=0,
                 platform='auto', mutations=(), chain='', ph=DEFAULT_PH):
    '''
    Minimize a structure, the way A4D's "stability calculations" option
    optimises its input with FoldX. Makes no PyMOL calls

    :param out_path: where the minimized structure is written, without
        hydrogens
    :return: {energy_before, energy_after (kcal/mol), platform, atoms}
    '''
    return _minimized(pdb_path, out_path, mutations, chain, None, tolerance,
                      max_iterations, platform, ph)


def dynamics_pdb(pdb_path, folder, n_models=MD_MODELS, steps=MD_STEPS,
                 temperature=MD_TEMPERATURE, platform='auto',
                 progress=None, cancelled=None):
    '''
    Sample fluctuations with implicit solvent molecular dynamics and
    reduce them to a few representative models, in place of A4D's
    CABS-flex. Makes no PyMOL calls

    :param folder: where the models are written
    :param n_models: how many representatives to keep
    :param progress: callable(done, total)
    :param cancelled: callable() -> True to stop
    :return: {'models': [path], 'rmsf': [(chain, resi, Angstrom)],
        'frames': int, 'platform': str}
    '''
    import numpy
    openmm, app, unit, _ = _openmm()
    app, topology, positions, system = _prepared(pdb_path)
    integrator = openmm.LangevinMiddleIntegrator(
        temperature * unit.kelvin, 1.0 / unit.picosecond,
        MD_TIMESTEP * unit.picoseconds)
    context, used = _platform(openmm, system, positions, integrator, platform)

    openmm.LocalEnergyMinimizer.minimize(
        context, 10.0 * unit.kilojoule_per_mole / unit.nanometer)
    context.setVelocitiesToTemperature(temperature * unit.kelvin)
    steps = max(int(steps), MD_FRAMES)
    every = max(1, (steps - MD_EQUILIBRATION) // MD_FRAMES)
    context.getIntegrator().step(min(MD_EQUILIBRATION, steps))

    alphas = [a.index for a in topology.atoms() if a.name == 'CA']
    frames, taken = [], 0
    while taken < steps - MD_EQUILIBRATION:
        if cancelled is not None and cancelled():
            raise Cancelled()
        context.getIntegrator().step(every)
        taken += every
        state = context.getState(getPositions=True)
        frames.append(state.getPositions(asNumpy=True)
                      .value_in_unit(unit.angstrom).copy())
        if progress is not None:
            progress(taken, steps - MD_EQUILIBRATION)

    # keep one model per cluster of the trajectory, plus RMSF
    from . import cluster
    carbons = [frame[alphas] for frame in frames]
    matrix = cluster.rmsd_matrix(carbons, carbons)
    picked = _medoids(matrix, min(n_models, len(frames)))
    names = []
    for order, index in enumerate(picked):
        path = os.path.join(folder, 'model_%d.pdb' % order)
        _save(app, topology,
              unit.Quantity(frames[index], unit.angstrom), path)
        names.append(path)

    stacked = numpy.array(carbons)
    reference = stacked.mean(axis=0)
    moved = []
    for frame in stacked:
        rotation, translation = cluster.kabsch(frame, reference)
        moved.append(frame @ rotation.T + translation)
    moved = numpy.array(moved)
    swing = numpy.sqrt(((moved - moved.mean(axis=0)) ** 2)
                       .sum(axis=2).mean(axis=0))
    rmsf = []
    for atom, value in zip([a for a in topology.atoms() if a.name == 'CA'],
                           swing):
        residue = atom.residue
        rmsf.append((residue.chain.id.strip() or '-',
                     residue.id + residue.insertionCode.strip(),
                     float(value)))
    return {'models': names, 'rmsf': rmsf, 'frames': len(frames),
            'platform': used}


def _medoids(matrix, count):
    '''
    Cluster a distance matrix into "count" groups and return the index
    of the most central member of each, largest group first
    '''
    import numpy
    n = len(matrix)
    if count >= n:
        return list(range(n))
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    tree = linkage(squareform(matrix, checks=False), method='average')
    labels = fcluster(tree, t=count, criterion='maxclust')
    groups = {}
    for index, label in enumerate(labels):
        groups.setdefault(label, []).append(index)
    found = []
    for label in sorted(groups, key=lambda k: -len(groups[k])):
        members = groups[label]
        inside = matrix[numpy.ix_(members, members)]
        found.append(members[int(numpy.argmin(inside.sum(axis=1)))])
    return found


class Cancelled(Exception):
    pass


# --- automated mutations ----------------------------------------------


class Suggestion:
    '''
    One point mutation and what it did

    :ivar chain, resi, old, new: where and what
    :ivar ddg: the change in folding stability, kcal/mol, from EvoEF
        (a free force field in FoldX's mould). Positive means less
        stable. None when EvoEF is not installed
    :ivar d_average, d_total, d_maximum: the change in the A4D score
        summaries; negative means more soluble
    :ivar estimated: the predicted change before the mutant was built
    '''

    def __init__(self, chain, resi, old, new):
        self.chain = chain
        self.resi = str(resi)
        self.old = old
        self.new = new
        self.ddg = None
        self.d_average = None
        self.d_total = None
        self.d_maximum = None
        self.estimated = None
        self.error = ''

    @property
    def label(self):
        return '%s%s%s' % (self.old, self.resi, self.new)

    @property
    def pdbfixer(self):
        return '%s-%s-%s' % (ONE_TO_THREE[self.old], self.resi,
                             ONE_TO_THREE[self.new])

    @property
    def solubilizing(self):
        return self.d_average is not None and self.d_average < 0

    @property
    def risk(self):
        '''
        How the A4D server marks the energetic effect: '', 'warn'
        above ENERGY_WARN kcal/mol, 'bad' above ENERGY_BAD
        '''
        if self.ddg is None:
            return ''
        if self.ddg > ENERGY_BAD:
            return 'bad'
        if self.ddg > ENERGY_WARN:
            return 'warn'
        return ''

    def __repr__(self):
        return '<%s %s dA4D %s>' % (self.chain, self.label,
                                    '%.4f' % self.d_average
                                    if self.d_average is not None else '?')


def evoef_binary():
    '''
    The EvoEF executable, which stands in for FoldX when the stability
    effect of a mutation is wanted: $EVOEF, one on the PATH, or one in
    an "evoef" folder beside the Python environment, beside its
    checkout, or in the home directory

    :return: its path, or ''
    '''
    where = os.environ.get('EVOEF', '')
    if where and os.path.isfile(where):
        return where
    found = shutil.which('EvoEF') or shutil.which('evoef')
    if found:
        return found
    beside = os.path.dirname(sys.prefix)
    for root in (beside, os.path.dirname(beside), os.path.expanduser('~')):
        for folder in ('evoef', 'EvoEF'):
            for inside in ('repo', '', 'bin'):
                path = os.path.join(root, folder, inside, 'EvoEF')
                if os.path.isfile(path) and os.access(path, os.X_OK):
                    return path
    return ''


def evoef_available():
    return bool(evoef_binary())


def _evoef(binary, folder, *arguments):
    '''
    Run EvoEF in a directory of its own. It writes its output files to
    the working directory under the bare basename of --pdb, so it is
    given one it may fill

    :return: its standard output
    '''
    import subprocess
    done = subprocess.run([binary] + list(arguments), cwd=folder,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True)
    if done.returncode != 0:
        raise RuntimeError('EvoEF failed: %s'
                           % (done.stderr or done.stdout)[-400:])
    return done.stdout


def _evoef_total(text):
    '''
    The "Total = <value>" line of a ComputeStability run
    '''
    import re
    found = re.search(r'^Total\s+=\s+(-?\d+\.?\d*)', text, re.M)
    if not found:
        raise RuntimeError('EvoEF printed no total')
    return float(found.group(1))


def evoef_ddg(pdb_path, suggestions, folder=''):
    '''
    The folding stability change of each mutation, with EvoEF in the
    place A4D uses FoldX. The structure is repaired once, every mutant
    is built in one pass, and each is compared with the wild type
    companion EvoEF writes beside it, which went through the same side
    chain optimisation

    Fills Suggestion.ddg and returns the suggestions it could measure

    :param suggestions: [Suggestion]; insertion codes are skipped,
        EvoEF's mutant file has nowhere to put them
    '''
    binary = evoef_binary()
    if not binary:
        return []
    own = not folder
    folder = folder or tempfile.mkdtemp(prefix='aggrescan_evoef_')
    os.makedirs(folder, exist_ok=True)
    try:
        shutil.copyfile(pdb_path, os.path.join(folder, 'input.pdb'))
        _evoef(binary, folder, '--command=RepairStructure',
               '--pdb=input.pdb')
        repaired = 'input_Repair.pdb'
        if not os.path.isfile(os.path.join(folder, repaired)):
            raise RuntimeError('EvoEF wrote no repaired structure')
        # one line per mutation, in FoldX's own syntax
        wanted = [one for one in suggestions if one.resi.isdigit()]
        if not wanted:
            return []
        with open(os.path.join(folder, 'list.txt'), 'w') as handle:
            for one in wanted:
                handle.write('%s%s%s%s;\n' % (
                    one.old, '_' if one.chain == '-' else one.chain,
                    one.resi, one.new))
        _evoef(binary, folder, '--command=BuildMutant',
               '--pdb=%s' % repaired, '--mutant_file=list.txt')
        stem = repaired[:-4]
        done = []
        for index, one in enumerate(wanted, 1):
            mutant = os.path.join(folder, '%s_Model_%04d.pdb' % (stem, index))
            wild = os.path.join(folder, '%s_Model_%04d_WT.pdb'
                                % (stem, index))
            if not (os.path.isfile(mutant) and os.path.isfile(wild)):
                continue
            after = _evoef_total(_evoef(
                binary, folder, '--command=ComputeStability',
                '--pdb=%s' % os.path.basename(mutant)))
            before = _evoef_total(_evoef(
                binary, folder, '--command=ComputeStability',
                '--pdb=%s' % os.path.basename(wild)))
            one.ddg = after - before
            done.append(one)
        return done
    finally:
        if own:
            shutil.rmtree(folder, ignore_errors=True)


def parse_mutations(text, result=None):
    '''
    Read mutations the way the A4D server's "mutate residues" page
    collects them: "L8E", "A/L8E", "A L8E", separated by commas or
    whitespace. The chain may be left out when the structure has one

    :param result: a Result to check the mutations against, so a wrong
        starting residue is caught here rather than deep in PDBFixer
    :return: [Suggestion]
    :raises ValueError: on anything unreadable or not in the structure
    '''
    import re
    chains = sorted({r.chain for r in result.residues}) if result else []
    letters = {r.key: r.letter for r in result.residues} if result else {}
    found = []
    for piece in re.split(r'[,;\s]+', str(text).strip()):
        if not piece:
            continue
        chain = ''
        token = piece
        if '/' in token:
            chain, token = token.split('/', 1)
        match = re.match(r'^([A-Za-z])(\d+[A-Za-z]?)([A-Za-z])$', token)
        if not match:
            raise ValueError('cannot read the mutation "%s": write it '
                             'like L8E, or A/L8E to name the chain'
                             % piece)
        old, resi, new = (match.group(1).upper(), match.group(2),
                          match.group(3).upper())
        for one in (old, new):
            if one not in ONE_TO_THREE:
                raise ValueError('"%s" is not an amino acid, in "%s"'
                                 % (one, piece))
        if not chain:
            if len(chains) == 1:
                chain = chains[0]
            elif not chains:
                raise ValueError('name the chain of "%s" (as A/%s)'
                                 % (piece, token))
            else:
                raise ValueError('"%s" needs a chain: the structure has '
                                 '%s' % (piece, ', '.join(chains)))
        if letters:
            there = letters.get((chain, resi))
            if there is None:
                raise ValueError('there is no residue %s in chain %s'
                                 % (resi, chain))
            if there != old:
                raise ValueError('residue %s of chain %s is %s, not %s'
                                 % (resi, chain, AA_NAMES[there],
                                    AA_NAMES[old]))
        found.append(Suggestion(chain, resi, old, new))
    if not found:
        raise ValueError('no mutations given')
    return found


def pam250():
    '''
    The PAM250 log-odds matrix the evolutionary mode scores
    substitutions with

    :return: {(from, to): value}
    '''
    from Bio.Align import substitution_matrices
    matrix = substitution_matrices.load('PAM250')
    found = {}
    for one in matrix.alphabet:
        for other in matrix.alphabet:
            if one in ONE_TO_THREE and other in ONE_TO_THREE:
                found[(one, other)] = float(matrix[one, other])
    return found


def targets(letter, mode='charged'):
    '''
    What a residue may be mutated into

    :param mode: "charged" for A3D's charged set, "evolutionary" for
        the substitutions PAM250 calls conserved or neutral
    :return: [one-letter code]
    '''
    if mode == 'evolutionary':
        matrix = pam250()
        found = [other for other in sorted(ONE_TO_THREE)
                 if other != letter
                 and matrix.get((letter, other), -99) >= PAM_THRESHOLD]
        return found
    return [one for one in CHARGED_MUTATIONS if one != letter]


def _swapped(result, chain, resi, new):
    '''
    Rescore a structure with one residue's identity changed but its
    shape kept: a cheap estimate of what a mutation would do, used to
    choose which mutants are worth building

    :return: Result
    '''
    residues = []
    for one in result.residues:
        if one.key == (chain, str(resi)):
            # a built side chain would change its exposure too, and its
            # pKa is taken as the standard one
            intrinsic = intrinsic_score(new, one.exposure, result.ph)
            copy = Residue(one.chain, one.resi, new, one.exposure,
                           intrinsic, one.coords)
        else:
            copy = Residue(one.chain, one.resi, one.letter, one.exposure,
                           one.intrinsic, one.coords)
        residues.append(copy)
    spread(residues, result.distance)
    return Result(residues, result.distance, result.ph, result.name)


def candidates(result, mode='charged', excluded=(), max_sites=MAX_SITES,
               per_site=0):
    '''
    The mutations worth trying, most promising first: A3D picks the
    most aggregation prone residues, and this then keeps the ones a
    cheap rescore says would help

    :param excluded: [(chain, resi)] the user protects
    :param max_sites: how many positions to mutate at most
    :param per_site: how many variants of each to keep (0: the mode's
        default)
    :return: [Suggestion] with "estimated" set
    '''
    if per_site <= 0:
        per_site = (MAX_PER_SITE_EVOLUTIONARY if mode == 'evolutionary'
                    else MAX_PER_SITE)
    protected = {(str(c), str(r)) for c, r in excluded}
    base = result.average
    sites = []
    for one in sorted(result.residues, key=lambda r: -r.score):
        if one.score <= MUTATION_THRESHOLD or one.score == 0.0:
            continue
        if one.key in protected:
            continue
        # A3D leaves residues that are already a charged target alone
        if mode == 'charged' and one.letter in CHARGED_MUTATIONS:
            continue
        sites.append(one)
        if len(sites) >= max_sites:
            break
    found = []
    for one in sites:
        tried = []
        for new in targets(one.letter, mode):
            guess = _swapped(result, one.chain, one.resi, new).average - base
            if guess >= 0:
                continue    # not solubilizing even before it is built
            item = Suggestion(one.chain, one.resi, one.letter, new)
            item.estimated = guess
            tried.append(item)
        tried.sort(key=lambda s: s.estimated)
        found.extend(tried[:per_site])
    found.sort(key=lambda s: s.estimated)
    return found


def evaluate(pdb_path, suggestions, distance=DEFAULT_DISTANCE, ph=None,
             pka=DEFAULT_PKA, base=None, limit=MAX_SUGGESTIONS,
             platform='auto', folder='', progress=None, cancelled=None):
    '''
    Build, minimize and rescore each mutant. Makes no PyMOL calls

    :param base: the wild type Result, minimized the same way
    :param limit: how many mutants to build at most
    :return: [Suggestion] with the measured changes filled in, the most
        solubilizing first
    '''
    own = not folder
    folder = folder or tempfile.mkdtemp(prefix='aggrescan_mut_')
    os.makedirs(folder, exist_ok=True)
    try:
        optimize_pdb(pdb_path, os.path.join(folder, 'wild.pdb'),
                     platform=platform)
        if base is None:
            base = score_pdb(os.path.join(folder, 'wild.pdb'), distance, ph,
                             pkas(os.path.join(folder, 'wild.pdb'))
                             if ph is not None and pka == 'propka'
                             and propka_available() else None,
                             signed=ph is not None and pka == 'propka')
        picked = list(suggestions)[:limit]
        for index, item in enumerate(picked):
            if cancelled is not None and cancelled():
                raise Cancelled()
            if progress is not None:
                progress(index, len(picked))
            out = os.path.join(folder, '%s%s.pdb' % (item.chain, item.label))
            try:
                optimize_pdb(pdb_path, out, platform=platform,
                             mutations=[item.pdbfixer], chain=item.chain)
                values = pkas(out) if (ph is not None and pka == 'propka'
                                       and propka_available()) else None
                after = score_pdb(out, distance, ph, values,
                                  signed=values is not None)
            except Exception as e:
                item.error = str(e)
                continue
            item.d_average = after.average - base.average
            item.d_total = after.total - base.total
            item.d_maximum = after.maximum - base.maximum
        if progress is not None:
            progress(len(picked), len(picked))
        # "only solubilizing mutations are considered" (NAR 2024): the
        # cheap estimate chose these, the rebuilt mutant decides
        done = [s for s in picked
                if s.d_average is not None and s.d_average < 0]
        if done and evoef_available():
            try:
                evoef_ddg(pdb_path, done, os.path.join(folder, 'evoef'))
            except Exception:
                pass    # the scores still stand without the stability
        # the server ranks by solubilization, with the stable ones first
        done.sort(key=lambda s: (s.risk != '', s.d_average))
        return done
    finally:
        if own:
            shutil.rmtree(folder, ignore_errors=True)


# --- one whole analysis, for a worker thread --------------------------


class Job:
    '''
    One aggregation analysis, following the A4D pipeline: optimise the
    input, score it, and then whichever of the pH scan, the dynamic
    mode and the automated mutations were asked for.

    Created by prepare() (which reads PyMOL), then run() (which makes
    no PyMOL calls, so a worker thread can own it), then apply().
    '''

    def __init__(self):
        self.selection = 'polymer'
        self.object = ''
        self.name = ''
        self.state = 1
        self.distance = DEFAULT_DISTANCE
        self.ph = True
        self.ph_values = PH_VALUES
        self.pka = DEFAULT_PKA
        self.optimize = True
        self.dynamic = False
        self.mutations = ''
        self.mutate = ''
        self.globular = False
        self.max_sites = MAX_SITES
        self.max_mutants = MAX_SUGGESTIONS
        self.excluded = []
        self.metric = DEFAULT_METRIC
        self.platform = 'auto'
        self.n_models = MD_MODELS
        self.md_steps = MD_STEPS
        self.folder = ''
        self.input_path = ''
        self.scored_path = ''
        # results
        self.result = None
        self.ph_results = []
        self.models = []
        self.rmsf = []
        self.suggestions = []
        self.mutated = []
        self.report = {}
        self.status = ''
        self.fraction = 0.0
        self.progress = None
        self._cancel = None

    def cancel(self):
        if self._cancel is not None:
            self._cancel.set()

    def cancelled(self):
        return self._cancel is not None and self._cancel.is_set()

    def _say(self, text, fraction=None):
        self.status = text
        if fraction is not None:
            self.fraction = fraction
        if self.progress is not None:
            self.progress(text, fraction)

    def check_options(self):
        '''
        :raises ValueError: if the options cannot be combined, the way
            the A4D server refuses them
        '''
        if self.distance not in DISTANCES:
            raise ValueError('distance must be one of: %s'
                             % ', '.join(str(d) for d in DISTANCES))
        if self.pka not in [v for _, v in PKA_MODES]:
            raise ValueError('pka must be one of: %s'
                             % ', '.join(v for _, v in PKA_MODES))
        if self.mutations and self.mutations not in \
                [v for _, v in MUTATION_MODES]:
            raise ValueError('mutations must be one of: %s'
                             % ', '.join(v for _, v in MUTATION_MODES))
        if self.ph and self.dynamic:
            raise ValueError('the pH calculation cannot be combined with '
                             'the dynamic mode')
        if self.ph and self.mutations:
            raise ValueError('the pH calculation cannot be combined with '
                             'automated mutations')
        if self.dynamic and self.mutations:
            raise ValueError('the dynamic mode cannot be combined with '
                             'automated mutations')
        if self.mutations and not self.optimize:
            raise ValueError('automated mutations need the structure '
                             'optimisation (as FoldX does in A4D)')
        if self.mutate and not self.optimize:
            raise ValueError('mutating residues needs the structure '
                             'optimisation (as FoldX does in A4D)')
        missing = available(self.ph and self.pka == 'propka',
                            self.optimize or self.dynamic or
                            bool(self.mutations))
        if missing:
            raise ValueError(missing)

    def run(self):
        '''
        Do the analysis. Fills result, ph_results, models, suggestions
        and report

        :raises Cancelled: if cancel() was called
        '''
        import threading
        if self._cancel is None:
            self._cancel = threading.Event()
        self.report = {'distance': self.distance, 'optimized': self.optimize}
        self.scored_path = self.input_path

        wanted = {}
        if self.mutate:
            # the residues the user asked for, checked against the
            # structure before anything is built
            self.mutated = parse_mutations(
                self.mutate, score_pdb(self.input_path, self.distance))
            for one in self.mutated:
                wanted.setdefault(one.chain, []).append(one.pdbfixer)
            self.report['mutated'] = [one.label for one in self.mutated]

        if self.optimize:
            self._say('Building the mutant...' if wanted else
                      'Optimising the structure with OpenMM...', 0.0)
            out = os.path.join(self.folder, 'optimized.pdb')
            energies = optimize_pdb(self.input_path, out,
                                    platform=self.platform,
                                    mutations=wanted or ())
            self.scored_path = out
            self.report['optimization'] = energies
            self._check()

        self._say('Scoring the structure...', 0.2)
        self.result = score_pdb(self.scored_path, self.distance, None,
                                name=self.name)
        self.report['residues'] = len(self.result)
        self.report['globular'] = self.result.is_globular()
        self._check()

        if self.ph:
            def step(done, total):
                self._say('Scoring at pH %s...'
                          % (self.ph_values[min(done, total - 1)]),
                          0.3 + 0.6 * done / max(total, 1))
                self._check()
            self._say('Predicting pKa values...', 0.3)
            self.ph_results = scan_ph(self.scored_path, self.distance,
                                      self.ph_values, self.pka, self.name,
                                      step)
        if self.dynamic:
            self._run_dynamic()
        if self.mutations:
            self._run_mutations()
        self._say('Done.', 1.0)
        return self.report

    def _check(self):
        if self.cancelled():
            raise Cancelled()

    def _run_dynamic(self):
        self._say('Sampling fluctuations with molecular dynamics...', 0.3)

        def step(done, total):
            self._say('Molecular dynamics: %d%%' % (100 * done / total),
                      0.3 + 0.4 * done / max(total, 1))
            self._check()

        models = os.path.join(self.folder, 'models')
        os.makedirs(models, exist_ok=True)
        found = dynamics_pdb(self.scored_path, models, self.n_models,
                             self.md_steps, platform=self.platform,
                             progress=step, cancelled=self.cancelled)
        self.rmsf = found['rmsf']
        self.report['dynamic'] = {'frames': found['frames'],
                                  'platform': found['platform']}
        scored = [('input', self.scored_path, self.result)]
        for index, path in enumerate(found['models']):
            self._say('Scoring model %d of %d...'
                      % (index + 1, len(found['models'])),
                      0.7 + 0.3 * index / max(len(found['models']), 1))
            self._check()
            scored.append(('model_%d' % index, path,
                           score_pdb(path, self.distance, None,
                                     name='model_%d' % index)))
        # the most aggregation prone first, as A3D ranks them
        scored.sort(key=lambda row: -row[2].average)
        self.models = scored

    def _run_mutations(self):
        self._say('Choosing residues to mutate...', 0.4)
        picked = candidates(self.result, self.mutations, self.excluded,
                            self.max_sites)
        self.report['mutation_candidates'] = len(picked)
        if not picked:
            self.suggestions = []
            return

        def step(done, total):
            self._say('Building and scoring mutant %d of %d...'
                      % (min(done + 1, total), total),
                      0.45 + 0.5 * done / max(total, 1))
            self._check()

        self.suggestions = evaluate(
            self.scored_path, picked, self.distance, None, self.pka,
            self.result, self.max_mutants, self.platform,
            os.path.join(self.folder, 'mutants'), step, self.cancelled)

    def chosen(self):
        '''
        The result the structure should be coloured by: the pH closest
        to neutral when the pH scan ran, else the static one
        '''
        if self.ph_results:
            return min(self.ph_results,
                       key=lambda row: abs(row[0] - DEFAULT_PH))[1]
        return self.result

    def apply(self, _self=cmd, color=True, select=True, models=False):
        '''
        Show the result in the session: scores in the B-factors, the
        aggregation prone residues as a selection, and the A4D colours
        '''
        if self.result is None:
            return {}
        done = apply_scores(self.object, self.chosen(), color, select,
                            _self=_self)
        if models and self.models:
            for name, path, _ in self.models:
                if name == 'input':
                    continue
                _self.load(path, _self.get_legal_name(
                    '%s_%s' % (self.object, name)))
        return done

    def save(self, folder, prefix=''):
        '''
        Write the results as files: the scores as CSV, the structure
        carrying them as PDB, and, when the pH scan ran, one PDB and
        one row per pH tested. The scores live in the B-factors, as
        A3D writes them

        :param folder: an existing directory to write into
        :param prefix: the stem of each name {default: the object's}
        :return: [path] in the order written
        :raises ValueError: if the analysis has been cleaned up
        '''
        if self.result is None:
            raise ValueError('nothing to save: run the analysis first')
        if not (self.scored_path and os.path.isfile(self.scored_path)):
            raise ValueError('the analysed structure is gone: run the '
                             'analysis again before saving')
        prefix = prefix or self.name or 'aggregation'
        written = []

        def put(name, text):
            path = os.path.join(folder, '%s_%s' % (prefix, name))
            with open(path, 'w') as handle:
                handle.write(text)
            written.append(path)
            return path

        def structure(name, result):
            path = os.path.join(folder, '%s_%s.pdb' % (prefix, name))
            write_scored_pdb(self.scored_path, path, result)
            written.append(path)

        put('A4D.csv', self.result.csv())
        structure('A4D', self.result)

        if self.ph_results:
            # one file per pH tested, named by it
            for value, result in self.ph_results:
                structure('pH%.1f' % value, result)
                put('pH%.1f.csv' % value, result.csv())
            rows = ['pH,maximum,average,minimum,total,metric,value']
            for value, result in self.ph_results:
                number, used = result.metric(self.metric)
                rows.append('%.1f,%.4f,%.4f,%.4f,%.4f,%s,%.4f'
                            % (value, result.maximum, result.average,
                               result.minimum, result.total, used, number))
            put('pH_summary.csv', '\n'.join(rows) + '\n')

        if self.models:
            rows = ['model,average,maximum,minimum,total']
            for name, path, result in self.models:
                rows.append('%s,%.4f,%.4f,%.4f,%.4f'
                            % (name, result.average, result.maximum,
                               result.minimum, result.total))
                if os.path.isfile(path):
                    out = os.path.join(folder, '%s_%s.pdb' % (prefix, name))
                    write_scored_pdb(path, out, result)
                    written.append(out)
            put('models.csv', '\n'.join(rows) + '\n')

        if self.suggestions:
            rows = ['mutation,chain,dA4D_average,dA4D_total,'
                    'dA4D_maximum,ddG_kcal_per_mol']
            for one in self.suggestions:
                rows.append('%s,%s,%.4f,%.4f,%.4f,%s'
                            % (one.label, one.chain, one.d_average,
                               one.d_total, one.d_maximum,
                               '%.2f' % one.ddg
                               if one.ddg is not None else ''))
            put('mutations.csv', '\n'.join(rows) + '\n')
        return written

    def cleanup(self):
        if self.folder and os.path.isdir(self.folder):
            shutil.rmtree(self.folder, ignore_errors=True)
        self.folder = ''


def prepare(selection='polymer', distance=DEFAULT_DISTANCE, ph=True,
            pka=DEFAULT_PKA, optimize=True, dynamic=False, mutations='',
            mutate='', state=1, name='', _self=cmd, **options):
    '''
    Read a structure out of PyMOL into a Job ready to run()

    :return: Job
    '''
    job = Job()
    job.selection = selection
    objects = _self.get_object_list('(%s)' % selection) or []
    if not objects:
        raise ValueError('no atoms in selection "%s"' % selection)
    job.object = objects[0]
    job.name = name or job.object
    job.state = int(state) if int(state) > 0 else _self.get_state()
    job.distance = int(distance)
    job.ph = bool(ph)
    job.pka = pka
    job.optimize = bool(optimize)
    job.dynamic = bool(dynamic)
    job.mutations = mutations or ''
    job.mutate = mutate or ''
    for key, value in options.items():
        if not hasattr(job, key):
            raise ValueError('unknown option "%s"' % key)
        setattr(job, key, value)
    job.check_options()
    job.folder = tempfile.mkdtemp(prefix='aggrescan_')
    job.input_path = os.path.join(job.folder, 'input.pdb')
    picked = selection
    if job.globular:
        picked, dropped, predicted = globular_selection(
            '(%s) and %s' % (selection, SELECTION), _self=_self)
        job.report['globular_dropped'] = dropped
        job.report['looks_predicted'] = predicted
    write_pdb(picked, job.input_path, job.state, _self)
    return job


# --- showing and reporting --------------------------------------------


def preset(selection='(all)', _self=cmd):
    '''
    The "aggregation" preset (Presets > aggregation, or
    preset.aggregation): a half transparent surface coloured by the
    B-factors, which is where this module leaves its scores. Blue is
    soluble and red aggregation prone, spread over the range the
    structure actually spans, so the weakest and strongest residues of
    this one structure are always at the ends of the ramp.

    That differs from apply_scores, which holds the ramp at a fixed
    range so two structures can be compared by eye. Run "aggregation"
    first: the preset colours whatever is in the B-factors, which in
    an untouched crystal structure is the crystallographic B

    :return: (lowest, highest) B-factor the ramp was stretched over,
        or None if the selection has no atoms
    '''
    from pymol import preset as presets
    name, objects, _ = presets._prepare(selection, _self=_self)
    try:
        values = []
        _self.iterate(name, 'rows.append(b)', space={'rows': values})
        if not values:
            return None
        low, high = min(values), max(values)
        _self.show('surface', name)
        # on the object, not the selection: an atom level transparency
        # would quietly win over a later "set transparency, ..., obj"
        _self.set('transparency', PRESET_TRANSPARENCY, objects)
        if high - low < 1e-6:
            # one value throughout: the middle of the ramp, not a
            # spectrum over nothing
            _self.color('white', name)
        else:
            _self.spectrum('b', COLOR_PALETTE, name, low, high)
        return low, high
    finally:
        _self.delete(name)


def apply_scores(obj, result, color=True, select=True, b_factors=True,
                 _self=cmd):
    '''
    Show a Result on a structure: the scores in the B-factors, the
    aggregation prone residues as a selection, and the A4D colours
    (blue soluble, white neutral, red aggregation prone) over a fixed
    range so two structures can be compared by eye

    :return: {residues, prone, hotspots, selection}
    '''
    scores = {r.key: r.score for r in result.residues}
    if b_factors:
        _self.alter('%' + obj,
                    'b = SCORES.get((chain or "-", resi), 0.0)',
                    space={'SCORES': scores})
        _self.sort(obj)
    if color:
        _self.spectrum('b', COLOR_PALETTE, '%' + obj,
                       -COLOR_RANGE, COLOR_RANGE)
    prone = result.prone
    name = ''
    if select and prone:
        name = _self.get_legal_name('%s_aggregation' % obj)
        _self.select(name, '(%%%s) and (%s)' % (obj, _residues_selection(
            [r.key for r in prone])), enable=0)
    return {'residues': len(result), 'prone': len(prone),
            'hotspots': sum(1 for r in result.residues if r.hotspot),
            'selection': name}


def _residues_selection(keys):
    by_chain = {}
    for chain, resi in keys:
        by_chain.setdefault(chain, []).append(str(resi))
    return ' or '.join(
        '(chain "%s" and resi %s)' % ('' if chain == '-' else chain,
                                      '+'.join(res))
        for chain, res in sorted(by_chain.items())) or 'none'


def write_scored_pdb(source_pdb, out_path, result):
    '''
    A copy of a structure carrying its scores, the way A3D writes
    output.pdb: each atom's B-factor is its residue's score at
    "%6.2f", and the occupancy is left alone. An atom whose residue
    was not scored (a buried one is still scored, at 0) gets 0.00

    :param source_pdb: the structure the scores were computed on
    :return: how many atoms were given a score
    '''
    scores = {r.key: r.score for r in result.residues}
    written = 0
    with open(source_pdb) as handle, open(out_path, 'w') as out:
        for line in handle:
            if line.startswith(('ATOM  ', 'HETATM')):
                chain = line[21].strip() or '-'
                resi = line[22:26].strip() + line[26].strip()
                score = scores.get((chain, resi))
                if score is not None:
                    written += 1
                # the columns are fixed width, so pad a short line first
                line = '%-80s' % line.rstrip('\n')
                line = '%s%6.2f%s\n' % (line[:60], score or 0.0, line[66:])
            out.write(line)
    return written


def report(result, limit=15, metric=DEFAULT_METRIC):
    '''
    A Result as a short table, the most aggregation prone first
    '''
    value, used = result.metric(metric)
    where = 'pH %.1f' % result.ph if result.ph is not None else \
        'Aggrescan scale'
    lines = [' Aggregation (%s, %g A sphere): %d residues, %d prone, '
             '%d above %.1f'
             % (where, result.distance, len(result), len(result.prone),
                sum(1 for r in result.residues if r.hotspot),
                HOTSPOT_THRESHOLD),
             '  average %.4f, maximum %.4f, minimum %.4f, total %.4f'
             % (result.average, result.maximum, result.minimum,
                result.total),
             '  %s score %.4f (%s structure)'
             % ('maximum' if used == 'max' else 'average', value,
                'globular' if result.is_globular() else 'extended')]
    if not result.residues:
        return '\n'.join(lines)
    lines.append(' %-6s %-8s %-8s %-8s %s'
                 % ('Chain', 'Residue', 'Score', 'Exposure', 'Prone'))
    ranked = sorted(result.residues, key=lambda r: -r.score)
    for one in (ranked[:limit] if limit else ranked):
        lines.append(' %-6s %-8s %-8.4f %-8.1f %s'
                     % (one.chain, '%s%s' % (one.letter, one.resi),
                        one.score, one.exposure,
                        'yes' if one.prone else ''))
    if limit and len(ranked) > limit:
        lines.append('  ... %d more' % (len(ranked) - limit))
    return '\n'.join(lines)


def format_ph(ph_results, metric=DEFAULT_METRIC):
    '''
    The pH scan as a table, with the metric the shape of the structure
    calls for
    '''
    if not ph_results:
        return ' no pH results'
    _, first = ph_results[0]
    _, used = first.metric(metric)
    lines = [' pH scan (%s, %s score recommended):'
             % ('globular' if first.is_globular() else 'extended',
                'maximum' if used == 'max' else 'average'),
             ' %-6s %-10s %-10s %-10s' % ('pH', 'Maximum', 'Average',
                                          'Total')]
    for value, result in ph_results:
        mark = ' <' if result.metric(metric)[0] == max(
            r.metric(metric)[0] for _, r in ph_results) else ''
        lines.append(' %-6.1f %-10.4f %-10.4f %-10.4f%s'
                     % (value, result.maximum, result.average,
                        result.total, mark))
    lines.append('  (< marks the pH with the highest aggregation '
                 'propensity)')
    return '\n'.join(lines)


def format_suggestions(suggestions, mode='charged'):
    '''
    The automated mutations as a table, the best first
    '''
    label = dict(MUTATION_MODES).get(mode, mode)
    if not suggestions:
        return (' no %s improved the predicted solubility'
                % label.lower())
    lines = [' %s, %d suggested:' % (label, len(suggestions)),
             ' %-8s %-7s %-11s %-11s %s'
             % ('Mutation', 'Chain', 'dA4D avg', 'dA4D total',
                'ddG (kcal/mol)')]
    for one in suggestions:
        flag = {'warn': ' !', 'bad': ' !!'}.get(one.risk, '')
        lines.append(' %-8s %-7s %-11.4f %-11.4f %s%s'
                     % (one.label, one.chain, one.d_average, one.d_total,
                        '%+.2f' % one.ddg if one.ddg is not None else '-',
                        flag))
    if any(one.ddg is not None for one in suggestions):
        lines.append('  ddG from EvoEF, positive means less stable: '
                     '! above %.1f, !! above %.1f kcal/mol'
                     % (ENERGY_WARN, ENERGY_BAD))
    else:
        lines.append('  no stability estimate: EvoEF is not installed')
    return '\n'.join(lines)


def format_models(models):
    '''
    The dynamic mode's models, the most aggregation prone first
    '''
    if not models:
        return ' no models'
    lines = [' Dynamic mode, %d models (most aggregation prone first):'
             % len(models),
             ' %-10s %-10s %-10s %-10s' % ('Model', 'Average', 'Maximum',
                                           'Total')]
    for name, _, result in models:
        lines.append(' %-10s %-10.4f %-10.4f %-10.4f'
                     % (name, result.average, result.maximum,
                        result.total))
    return '\n'.join(lines)


def aggregation(selection='polymer', distance=DEFAULT_DISTANCE, ph=0,
                optimize=0, dynamic=0, mutations='', mutate='',
                pka=DEFAULT_PKA, output='', color=1, select=1, quiet=0,
                state=1, *, _self=cmd):
    '''
DESCRIPTION

    "aggregation" predicts which residues of a folded protein drive
    its aggregation, with a free re-implementation of Aggrescan4D
    built on the published Aggrescan3D code. The score goes into the
    B-factors and colours the structure blue (soluble) through white
    to red (aggregation prone).

USAGE

    aggregation [ selection [, distance [, ph [, optimize
        [, dynamic [, mutations ]]]]]]

ARGUMENTS

    selection = str: the structure to score {default: polymer}

    distance = 5 or 10: radius of the sphere centred on each alpha
    carbon {default: 10}

    ph = 0/1 or float: 1 scores the 11 pH values from 4.0 to 9.0, a
    number scores that one pH, 0 uses the classic Aggrescan scale
    {default: 0}

    optimize = 0/1: minimize the structure first with OpenMM, which
    takes the place of A4D's FoldX step {default: 0}

    dynamic = 0/1: sample fluctuations with molecular dynamics and
    score a set of models, in place of A4D's CABS-flex {default: 0}

    mutations = charged or evolutionary: suggest solubilising point
    mutations {default: none}

    mutate = str: mutate these residues first and analyse the mutant,
    e.g. "L8E" or "A/L8E, B/V70R" {default: none}

    pka = propka or fixed: where the pKa of an ionizable residue comes
    from when scoring at a pH {default: propka}

    output = str: a directory to write the results into, as CSV and as
    PDB files carrying the scores in their B-factors. With the pH scan
    that is one PDB per pH tested {default: none}

NOTES

    The command scores the coordinates as they are; the Aggregation
    prediction plugin follows the A4D server and optimises first.

    A residue above 0 is aggregation prone and one above 0.5 is likely
    to matter. A buried residue (under 10% exposed) scores 0.

    As on the A4D server, the pH scan cannot be combined with the
    dynamic mode or with automated mutations.

EXAMPLE

    aggregation 1ubq
    aggregation 1ubq, ph=5.5
    aggregation polymer, optimize=1, mutations=charged
    aggregation 1ubq, optimize=1, mutate=L8E
    aggregation 1ubq, ph=1, output=/tmp/1ubq_aggregation

SEE ALSO

    predict_paratope, interface_analysis, minimize
    '''
    one_ph = None
    scan = False
    text = str(ph).strip().lower()
    if text in ('1', 'on', 'yes', 'true'):
        scan = True
    elif text not in ('', '0', 'off', 'no', 'false', 'none'):
        one_ph = float(text)

    job = prepare(selection, int(distance), scan, str(pka),
                  bool(int(optimize)), bool(int(dynamic)),
                  str(mutations).strip().lower(), str(mutate),
                  int(state), '', _self)
    try:
        job.run()
        if one_ph is not None:
            job.ph_results = scan_ph(job.scored_path, job.distance,
                                     (one_ph,), job.pka, job.name)
        saved = []
        if output:
            folder = _self.exp_path(str(output))
            os.makedirs(folder, exist_ok=True)
            saved = job.save(folder)
        shown = job.apply(_self, bool(int(color)), bool(int(select)))
        if not int(quiet):
            print(report(job.chosen()))
            if job.ph_results and len(job.ph_results) > 1:
                print(format_ph(job.ph_results))
            if job.models:
                print(format_models(job.models))
            if job.mutated:
                print(' mutated: %s'
                      % ', '.join('%s/%s' % (one.chain, one.label)
                                  for one in job.mutated))
            if job.mutations:
                print(format_suggestions(job.suggestions, job.mutations))
            if shown.get('selection'):
                print(' selection "%s": the aggregation prone residues'
                      % shown['selection'])
            if saved:
                print(' wrote %d file(s) to %s'
                      % (len(saved), os.path.dirname(saved[0])))
        return job
    finally:
        job.cleanup()


def extend(_self):
    _self.extend('aggregation', aggregation)
    _self.aggregation = aggregation
