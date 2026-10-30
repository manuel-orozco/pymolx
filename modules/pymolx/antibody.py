'''
Antibody builder (command "build_antibody", Antibody Builder dialog): a
full IgG model from its parts (two Fabs, the hinge, the Fc) placed on a
full-length IgG template.

1. parts: each Fab (heavy + light chain), the Fc (two chains), optionally
   the hinge (a structure, or a sequence); the template: a full IgG
   (two heavy chains, two light chains), by default PDB 1HZH (human
   IgG1). Chains are recognized by their sequence (variable domain
   framework motifs; the Fc by its CH2/CH3 motifs).
2. placement: every chain is aligned by sequence to its template chain
   (BLOSUM62, pymolx.seqcompare); each part is superposed on the matched
   C-alpha atoms (outliers rejected like "align"): the Fabs on their
   constant domains (the elbow angle varies between antibodies), the Fc
   on whichever of its two chain orders fits better.
3. stitching: each heavy chain is the Fab heavy chain + hinge + an Fc
   chain, joined where their template positions meet (overlapping
   residues are taken from the Fab, then the Fc; unaligned tails at the
   junction, e.g. tags, are left out). The hinge is the user's (its
   structure superposed on the template's hinge and the joined ends), or
   the template's own residues. Hinge residues without coordinates take
   the template's (whole residue if the same amino acid, else backbone
   and CB); the rest are left without coordinates.
   Arms can also be an scFv or a VHH (superposed on the template arm's
   variable domains, staying there or moved next to the hinge; a linker
   to the hinge optional) or none (one-armed).
   Units (scFv, VHH) can be appended to chain C-termini through a linker:
   placed pointing away from the rest of the model, the linker built.
4. refinement (optional, pymolx.prep): residues without coordinates are
   built (PDBFixer), hydrogens added, and the model minimized (OpenMM)
   with restraints on all but the hinge and the residues next to each
   junction.

The model's chains (default H, L for Fab1 and K, M for Fab2, as in 1HZH)
are numbered from 1; the atom property "custom" tells each residue's
part (FAB1, FAB2 or SCFV1, VHH2, ...; HNGE, FC; FUSH1, FUSL2, ... for
fused units, LINK for their linkers; select with e.g. "custom HNGE"),
also shown by color. (Not "segi": PyMOL orders atoms by segment first,
which would put the hinge after the Fc.)
'''

import os
import re

import numpy as np

from pymol import cmd

from . import seqcompare

PARTS = ('fab1', 'fab2', 'hinge', 'fc')    # (fab1, fab2: arms 1 and 2)
# arm formats: Fab (heavy + light chain), scFv (VH and VL in one chain),
# VHH (single domain), none (one-armed)
ARM_FORMATS = ('fab', 'scfv', 'vhh', 'none')
FORMAT_LABELS = {'fab': 'Fab', 'scfv': 'scFv', 'vhh': 'VHH', 'none': 'none'}
FORMAT_TAGS = {'fab': 'FAB', 'scfv': 'SCFV', 'vhh': 'VHH'}
PART_LABELS = {'fab1': 'Fab1', 'fab2': 'Fab2', 'hinge': 'Hinge',
               'fc': 'Fc'}
PART_TAGS = {'fab1': 'FAB1', 'fab2': 'FAB2', 'hinge': 'HNGE', 'fc': 'FC'}
# heavy chain colors; light chains a lighter shade (LIGHT_TINT)
PART_COLORS = {'fab1': '#e0883a', 'fab2': '#4a90d9', 'hinge': '#b267cf',
               'fc': '#35a98a', 'fusion': '#d6b23a'}
LIGHT_TINT = 0.45
# C-terminal fusions (scFv, VHH appended through a linker): where
FUSION_SITES = {'heavy1': 'heavy chain 1', 'heavy2': 'heavy chain 2',
                'light1': 'light chain 1', 'light2': 'light chain 2'}
FUSION_SITE_ALIASES = {'heavy': ('heavy1', 'heavy2'),
                       'light': ('light1', 'light2')}
DEFAULT_LINKER = 'GGGGSGGGGS'   # (G4S)2
# fused scFvs: the domain order of the structure ('') or VH-VL / VL-VH
# (reordered: the structure's own linker left out, the domains joined by
# DEFAULT_SCFV_LINKER, which the refinement builds)
ORIENTATIONS = ('', 'vh-vl', 'vl-vh')
DEFAULT_SCFV_LINKER = 'GGGGS' * 3   # (G4S)3
# where scFv and VHH arms go: on the template arm's Fv, or next to the hinge
SINGLE_POSITIONS = ('fv', 'hinge')
LINKERS = [('(G4S)2, flexible', 'GGGGSGGGGS'),
           ('(G4S)3, flexible', 'GGGGS' * 3),
           ('(G4S)4, flexible', 'GGGGS' * 4),
           ('(EAAAK)3, rigid', 'EAAAK' * 3), ('None (direct)', '')]
# (not 5DK3, IgG4: one CH2 domain is turned, so a usual Fc fitted on it
# lands on a Fab)
TEMPLATES = [('1hzh', 'Human IgG1 (b12), PDB 1HZH'),
             ('1igt', 'Mouse IgG2a, PDB 1IGT'),
             ('1igy', 'Mouse IgG1, PDB 1IGY')]
DEFAULT_TEMPLATE = '1hzh'
# hinge sequences (upper and core hinge)
HINGES = [('Human IgG1', 'EPKSCDKTHTCPPCP'),
          ('Human IgG2', 'ERKCCVECPPCP'),
          ('Human IgG3', 'ELKTPLGDTTHTCPRCP' + 'EPKSCDTPPPCPRCP' * 3),
          ('Human IgG4', 'ESKYGPPCPSCP'),
          ('Human IgG4 S228P', 'ESKYGPPCPPCP'),
          ('Mouse IgG1', 'VPRDCGCKPCICT'),
          ('Mouse IgG2a', 'EPRGPTIKPCPPCKCP')]

BACKBONE = ('N', 'CA', 'C', 'O')
# variable domain: framework 1 Cys, then the framework 2 Trp (like
# pymolx.cluster's motif, also with His: WYQH of b12's light chain)
V_MOTIF = r'C.{9,16}W[VIFLYM].[QKREH]'
FR4_WINDOW = (45, 110)  # framework 4 motif: residues after framework 2
FC_MOTIFS = (r'C[VIL]V[VIL]D[VI]S', r'LTC[LM][VI]')
MIN_PAIRS = 20        # matched C-alpha atoms needed to place a Fab or Fc
HINGE_WINDOW = 15     # template residues beyond the junction searched
FREE_FLANK = 2        # residues next to a junction moved in refinement
CA_STEP = 3.8         # A between consecutive C-alpha atoms
STRETCHED = 3.3       # A per residue: a span this long is nearly stretched
RELAXED = 2.8         # A per residue: a relaxed span
MIN_CUTOFF = 0.5      # A: superposition outliers are at least this far
REACH = 3.8           # A per residue: a joined end can reach this far
REACH_SLACK = 0.3     # A
OVERLAP = 2.0         # A: atoms of different parts this close overlap
OVERLAP_WARN = 10     # ... a warning from this many atom pairs
SINGLE_CLEARANCE = 3.0  # A: an scFv/VHH moves toward the hinge while it
SINGLE_OVERLAPS = 10    # has at most this many atom pairs this close
DISULFIDE_SEARCH = 8.0  # A: SG-SG of cysteine pairs offered as disulfides
DISULFIDE_BOND = 2.5    # A: SG-SG of disulfides made by distance
GLYCAN_CONTACT = 2.0  # A: het groups bonded to a chain go with it
FUSION_GAP = 6.0      # A: a fused unit's N-terminus from the chain's end
EXTENDED_STEP = 3.2   # A per residue of an extended linker
LINKER_SHAPES = ('relaxed', 'extended', 'collapsed')
COIL_STEP = 3.8       # A: a random coil's end-to-end is this times sqrt(N)
PULL_PS = 40          # ps of dynamics drawing a fused unit in ...
PULL_SPREAD = 1       # residues of it the pull acts on (the network
                      # passes the force on to the rest)
SETTLE_PS = 20       # ... and of letting it settle afterwards
FUSION_CLEARANCE = 6.5  # A: C-alpha atoms closer to others count as overlap
FUSION_POSES = 40     # best poses by C-alpha atoms checked with all atoms
CDR_OUTWARD = 60      # degrees: CDRs count as facing the solvent
# where a fused unit's binding site is turned: away from the antibody
# ('out', as far as clashes allow) or wherever it fits best ('any')
CDR_CHOICES = ('out', 'any')
GLYCOSYLATION_CONTACT = 3.0  # A: sugar C1 to Asn ND2 / Ser, Thr O
CLASH = 2.2           # A: template atoms this close to the parts' clash

def _three(letter):
    from .prep import ONE_LETTER
    return ONE_LETTER.get(letter, 'ALA')


class Residue:
    '''
    One residue: name and heavy atoms {name: (element, xyz, b)} (none: no
    coordinates); part (PARTS); its place in the template chain
    '''

    def __init__(self, resn, atoms=None, part=''):
        self.resn = resn
        self.atoms = atoms or {}
        self.part = part
        self.template = None    # template residue index
        self.threaded = False   # coordinates from the template
        self.light = False      # light chain or VL domain (lighter color)

    @property
    def letter(self):
        return seqcompare.ONE_LETTER.get(self.resn, 'X')

    @property
    def ca(self):
        atom = self.atoms.get('CA')
        return None if atom is None else atom[1]

    def moved(self, rotation, shift):
        result = Residue(self.resn, {
            name: (elem, xyz @ rotation.T + shift, b)
            for name, (elem, xyz, b) in self.atoms.items()}, self.part)
        result.template = self.template
        result.light = self.light
        return result


def _sequence(residues):
    return ''.join(r.letter for r in residues)


def _object_state(obj, state, _self):
    state = int(state)
    if state < 1:
        state = _self.get_state()
    return max(1, min(state, _self.count_states('%' + obj) or 1))


def read_chain(obj, chain, state=-1, _self=cmd):
    '''
    :return: [Residue] of a protein chain in order, with residues
        without coordinates (no atoms); first alternate location; no
        terminal OXT
    '''
    sele = '%%%s and chain "%s" and polymer.protein' % (obj, chain)
    order = {}
    _self.iterate(sele, '_ = order.setdefault(resi, resn)',
                  space={'order': order})
    residues = {resi: Residue(resn) for resi, resn in order.items()}
    rows = []
    _self.iterate_state(_object_state(obj, state, _self),
                        sele + ' and not hydro',
                        'rows.append((resi, name, elem, x, y, z, b))',
                        space={'rows': rows})
    for resi, name, elem, x, y, z, b in rows:
        atoms = residues[resi].atoms
        # (OXT: a chain end may become a junction; refinement adds it)
        if name not in atoms and name != 'OXT':
            atoms[name] = (elem, np.array((x, y, z)), b)
    return [residues[resi] for resi in order]


def attached_groups(obj, chain, state=-1, _self=cmd):
    '''
    Het groups bonded to a protein chain (glycans, covalent ligands): the
    connected het atoms with one bonded to the chain, within
    GLYCAN_CONTACT of it, or a sugar C1 at a glycosylation site.

    :return: [Residue]
    '''
    het = '(%%%s and not polymer and not solvent and not hydro)' % obj
    protein = '(%%%s and chain "%s" and polymer.protein)' % (obj, chain)
    temp = _self.get_unused_name('_ab_het')
    # bonded, or close; or a sugar's C1 at a glycosylation site (some
    # entries model the link long and without a bond, e.g. 1HZH: 2.6 A)
    _self.select(temp, '%s and ((neighbor %s) or (%s within %g of %s) or '
                 '(name C1 within %g of (%s and ((resn ASN and name ND2) or '
                 '(resn SER and name OG) or (resn THR and name OG1)))))' % (
                     het, protein, het, GLYCAN_CONTACT, protein,
                     GLYCOSYLATION_CONTACT, protein))
    count = _self.count_atoms(temp)
    while count:
        _self.select(temp, '%s and ((%s) extend 1)' % (het, temp))
        new = _self.count_atoms(temp)
        if new == count:
            break
        count = new
    rows = []
    if count:
        _self.iterate_state(_object_state(obj, state, _self), temp,
                            'rows.append((chain, resi, resn, name, elem, '
                            'x, y, z, b))', space={'rows': rows})
    _self.delete(temp)
    groups = {}
    for c, resi, resn, name, elem, x, y, z, b in rows:
        residue = groups.setdefault((c, resi), Residue(resn))
        residue.atoms.setdefault(name, (elem, np.array((x, y, z)), b))
    return list(groups.values())


# chains


def v_domains(sequence):
    '''
    Variable domains in order: [(kind, start, end)], kind "heavy" or
    "light" by the framework 4 motif (WGxG / FGxG, some VH LGxG followed
    by TVS), searched FR4_WINDOW residues after framework 2 (not in CDR2:
    e.g. VIWGDG); start about 25 residues before the framework 1 Cys,
    end 7 residues after framework 4
    '''
    found = []
    position = 0
    while True:
        first = re.search(V_MOTIF, sequence[position:])
        if not first:
            break
        cys = position + first.start()
        rest = position + first.end()
        low, high = FR4_WINDOW
        fr4 = re.search(r'[WFL]G.G', sequence[rest + low:rest + high])
        if not fr4:
            position = rest
            continue
        at = rest + low + fr4.start()
        kind = 'heavy' if sequence[at] == 'W' or (
            sequence[at] == 'L' and 'TVS' in sequence[at:at + 14]) \
            else 'light'
        end = min(len(sequence), at + 4 + 7)
        found.append((kind, max(position, cys - 25), end))
        position = end
    return found


def ch1_end(sequence):
    '''
    Index of the last CH1 residue of a heavy chain (before the hinge:
    ...KVDKKV in human, ...KVDKKI in mouse IgG), or None
    '''
    start = constant_start(sequence)
    found = re.search(r'KVDK[KR][VIA]', sequence[start:start + 130])
    return start + found.end() - 1 if found else None


def chain_kind(sequence):
    '''
    "heavy" or "light" (antibody variable domain followed by a constant
    one, by its framework 4 motif), "scfv" (two variable domains), "vhh"
    (one heavy-type variable domain alone), "fc" (CH2/CH3 without a
    variable domain) or ""
    '''
    domains = v_domains(sequence)
    if len(domains) >= 2:
        return 'scfv'
    if domains:
        kind, _, end = domains[0]
        if kind == 'heavy' and len(sequence) - end < 40:
            return 'vhh'
        return kind
    if any(re.search(m, sequence) for m in FC_MOTIFS):
        return 'fc'
    return ''


def constant_start(sequence):
    '''
    Index of the first residue after the variable domain (0 if none)
    '''
    domains = v_domains(sequence)
    return domains[0][2] if domains else 0


def protein_chains(obj, _self=cmd):
    '''
    {chain: (sequence, kind)} in chain order
    '''
    result = {}
    for chain in _self.get_chains('%%%s and polymer.protein' % obj) or []:
        residues = read_chain(obj, chain, _self=_self)
        sequence = _sequence(residues)
        result[chain] = (sequence, chain_kind(sequence))
    return result


def _contacts(obj, a, b, _self):
    return _self.count_atoms(
        '(%%%s and chain "%s" and name CA) within 8 of (%%%s and chain "%s" '
        'and name CA)' % (obj, a, obj, b))


def _pair(obj, heavy, light, _self):
    '''
    Heavy chains with the light chain each binds most: [(heavy, light)]
    '''
    pairs = []
    free = list(light)
    for h in heavy:
        if not free:
            break
        best = max(free, key=lambda c: _contacts(obj, h, c, _self))
        pairs.append((h, best))
        free.remove(best)
    return pairs


def detect_fab(obj, _self=cmd):
    '''
    :return: (heavy chain, light chain) or ("", "")
    '''
    chains = protein_chains(obj, _self)
    heavy = [c for c, (_, kind) in chains.items() if kind == 'heavy']
    light = [c for c, (_, kind) in chains.items() if kind == 'light']
    pairs = _pair(obj, heavy, light, _self)
    if pairs:
        return pairs[0]
    return (heavy or [''])[0], (light or [''])[0]


def detect_single(obj, fmt, _self=cmd):
    '''
    :return: the chain of an scFv or VHH ("" if none)
    '''
    chains = protein_chains(obj, _self)
    found = [c for c, (_, kind) in chains.items() if kind == fmt]
    if not found and fmt == 'vhh':
        found = [c for c, (_, kind) in chains.items() if kind == 'heavy']
    return found[0] if found else ''


def detect_fc(obj, _self=cmd):
    '''
    :return: (chain, chain): the Fc chains (a single chain twice)
    '''
    chains = protein_chains(obj, _self)
    found = [c for c, (_, kind) in chains.items() if kind == 'fc'] or \
        [c for c, (_, kind) in chains.items() if kind == 'heavy']
    if not found:
        found = list(chains)
    if not found:
        return ('', '')
    return (found[0], found[1] if len(found) > 1 else found[0])


def detect_hinge(obj, _self=cmd):
    chains = list(protein_chains(obj, _self))
    if not chains:
        return ('', '')
    return (chains[0], chains[1] if len(chains) > 1 else chains[0])


def detect_template(obj, _self=cmd):
    '''
    :return: (heavy 1, light 1, heavy 2, light 2) of a full IgG
    '''
    chains = protein_chains(obj, _self)
    heavy = [c for c, (seq, kind) in chains.items()
             if kind == 'heavy' and len(seq) >= 300]
    light = [c for c, (_, kind) in chains.items() if kind == 'light']
    if len(heavy) < 2 or len(light) < 2:
        raise ValueError(
            'the template "%s" needs two full heavy chains and two light '
            'chains (found heavy: %s; light: %s)' % (
                obj, ', '.join(heavy) or 'none', ', '.join(light) or 'none'))
    pairs = _pair(obj, heavy[:2], light, _self)
    return pairs[0][0], pairs[0][1], pairs[1][0], pairs[1][1]


def parse_chains(text, count):
    '''
    "HL", "H,L" or "H L" -> ["H", "L"] ("" -> [""] * count)
    '''
    text = str(text or '').strip()
    if not text:
        return [''] * count
    if ',' in text or ' ' in text:
        parts = [p.strip() for p in re.split(r'[,\s]+', text) if p.strip()]
    else:
        parts = list(text)
    if len(parts) != count:
        raise ValueError('expected %d chains in "%s"' % (count, text))
    return parts


def resolve(spec, _self=cmd):
    '''
    An object name, a file (loaded) or a PDB ID (fetched)

    :return: object name
    '''
    spec = str(spec).strip()
    if not spec:
        return ''
    names = _self.get_names('objects')
    if spec in names:
        return spec
    path = _self.exp_path(spec)
    if os.path.exists(path):
        before = set(names)
        _self.load(path, quiet=1)
        new = [n for n in _self.get_names('objects') if n not in before]
        if not new:
            raise ValueError('nothing loaded from %s' % spec)
        return new[0]
    if re.fullmatch(r'\d[0-9a-zA-Z]{3}', spec):
        if spec.lower() in names:
            return spec.lower()
        _self.fetch(spec.lower(), quiet=1)
        if spec.lower() in _self.get_names('objects'):
            return spec.lower()
    raise ValueError('"%s" is not an object, a file or a PDB ID' % spec)


# superposition


def _kabsch(mobile, target):
    '''
    rotation, shift with target ~ mobile @ rotation.T + shift
    '''
    cm, ct = mobile.mean(axis=0), target.mean(axis=0)
    u, _, vt = np.linalg.svd((mobile - cm).T @ (target - ct))
    d = np.sign(np.linalg.det(u @ vt))
    rotation = (u @ np.diag((1.0, 1.0, d)) @ vt).T
    return rotation, ct - cm @ rotation.T


def superpose(mobile, target, cycles=5, cutoff=2.0):
    '''
    Least-squares fit with outlier rejection like PyMOL's "align" (atoms
    beyond cutoff times the RMSD, and beyond MIN_CUTOFF, left out, up to
    "cycles" times)

    :return: rotation, shift, RMSD of the kept atoms, number kept
    '''
    mobile = np.asarray(mobile, float)
    target = np.asarray(target, float)
    keep = np.ones(len(mobile), bool)
    for cycle in range(cycles + 1):
        rotation, shift = _kabsch(mobile[keep], target[keep])
        d = np.linalg.norm(mobile @ rotation.T + shift - target, axis=1)
        rmsd = float(np.sqrt(np.mean(d[keep] ** 2)))
        new = keep & (d <= max(cutoff * rmsd, MIN_CUTOFF))
        if cycle == cycles or new.sum() == keep.sum() or new.sum() < 3:
            break
        keep = new
    return rotation, shift, rmsd, int(keep.sum())


def _map(template, residues, start=0, end=None):
    '''
    Sequence alignment of residues to template[start:end]

    :return: {residue index: template index}
    '''
    end = len(template) if end is None else end
    mapping = seqcompare._align_pair(_sequence(template[start:end]),
                                     _sequence(residues))
    return {i: t + start for t, i in mapping.items()}


def _ca_pairs(residues, mapping, template, first=0):
    pairs = [(residues[i].ca, template[t].ca)
             for i, t in mapping.items()
             if i >= first and residues[i].ca is not None and
             template[t].ca is not None]
    return [p[0] for p in pairs], [p[1] for p in pairs]


# assembly


class Assembly:
    '''
    The stitched model: chains {chain: [Residue]}, het groups {chain:
    [Residue]}, residues to move in refinement, fits and messages
    '''

    def __init__(self):
        self.name = ''
        self.chains = {}
        self.groups = {}
        self.heavy = []         # output heavy chains, light chains
        self.light = []
        self.free = set()       # (chain, resi)
        self.fits = []          # (label, RMSD, atoms)
        self.log = []
        self.warnings = []
        self.inputs = []        # objects used
        self.disulfides = None  # chosen [((chain, resi), (chain, resi))]
                                # or None: by distance
        self.formats = {'fab1': 'fab', 'fab2': 'fab'}
        self.labels = dict(PART_LABELS, linker='linker')
        self.fusion_sites = []  # sites with a fusion (FUSION_SITES)
        self.collapse = []      # fusions whose linker the dynamics
                                # collapses (refine_job)

    def tags(self):
        '''
        {part: "custom" tag}: FAB1, SCFV2, VHH1, ...; HNGE, FC
        '''
        tags = dict(PART_TAGS)
        for k, part in enumerate(('fab1', 'fab2')):
            fmt = self.formats[part]
            if fmt != 'none':
                tags[part] = FORMAT_TAGS[fmt] + str(k + 1)
        for site in self.fusion_sites:
            # (heavy1 -> FUSH1)
            tags['fusion_' + site] = 'FUS' + site[0].upper() + site[-1]
        tags['linker'] = 'LINK'
        return tags

    def light_residues(self):
        return {(chain, str(k)) for chain, residues in self.chains.items()
                for k, r in enumerate(residues, 1) if r.light}

    def say(self, text):
        self.log.append(text)

    def warn(self, text):
        self.warnings.append(text)
        self.log.append('Warning: ' + text)

    def parts(self):
        '''
        {(chain, resi): part}
        '''
        result = {}
        for chain, residues in self.chains.items():
            for k, r in enumerate(residues, 1):
                result[(chain, str(k))] = r.part
        for chain, groups in self.groups.items():
            for k, r in enumerate(groups, 1001):
                result[(chain, str(k))] = r.part
        return result

    def missing(self):
        return sum(1 for residues in self.chains.values()
                   for r in residues if not r.atoms)


def _thread(template_residue, resn, part):
    '''
    The template residue's coordinates for residue "resn": all atoms if
    the same amino acid, else the backbone and CB
    '''
    atoms = template_residue.atoms
    if template_residue.resn == resn:
        names = list(atoms)
    else:
        names = list(BACKBONE) + ([] if resn == 'GLY' else ['CB'])
    residue = Residue(resn, {n: atoms[n] for n in names if n in atoms}, part)
    residue.template = None
    residue.threaded = bool(residue.atoms)
    return residue


def _place_fab(label, part, obj, heavy, light, template, t_heavy, t_light,
               state, assembly, glycans, _self, obstacles=()):
    '''
    The Fab superposed on the template's Fab: on the constant domains, or
    on the whole Fab when its elbow angle would put it into the obstacles
    (coordinates of parts placed already) that way

    :return: heavy chain residues, light chain residues (template indices
        set), het groups of each
    '''
    h = read_chain(obj, heavy, state, _self)
    l = read_chain(obj, light, state, _self)
    if not h or not l:
        raise ValueError('%s: no protein chain "%s" or "%s" in "%s"' % (
            label, heavy, light, obj))
    th = read_chain(template, t_heavy, -1, _self)
    tl = read_chain(template, t_light, -1, _self)
    map_h = _map(th, h)
    map_l = _map(tl, l)
    mobile, target = [], []
    for residues, mapping, tres in ((h, map_h, th), (l, map_l, tl)):
        m, t = _ca_pairs(residues, mapping, tres,
                         constant_start(_sequence(residues)))
        mobile += m
        target += t
    if len(mobile) < MIN_PAIRS:
        raise ValueError(
            '%s: only %d residues of its constant domains match the '
            'template; check the heavy/light chains or use a closer '
            'template' % (label, len(mobile)))
    rotation, shift, rmsd, kept = superpose(mobile, target)
    what, total = 'constant-domain', len(mobile)
    if len(obstacles):
        from scipy.spatial import cKDTree
        tree = cKDTree(np.asarray(obstacles))
        xyz = np.array([xyz for r in h + l for _, xyz, _ in
                        r.atoms.values()])
        n = _overlaps(xyz @ rotation.T + shift, tree, SINGLE_CLEARANCE)
        if n > SINGLE_OVERLAPS:
            whole_m, whole_t = [], []
            for residues, mapping, tres in ((h, map_h, th), (l, map_l, tl)):
                m, t = _ca_pairs(residues, mapping, tres)
                whole_m += m
                whole_t += t
            fit = superpose(whole_m, whole_t)
            n_whole = _overlaps(xyz @ fit[0].T + fit[1], tree,
                                SINGLE_CLEARANCE)
            if n_whole < n:
                rotation, shift, rmsd, kept = fit
                what, total = 'Fab', len(whole_m)
                assembly.say('%s: fitted on the whole Fab (on its constant '
                             'domains, its elbow angle would put it into '
                             'the other parts)' % label)
    assembly.fits.append((label, rmsd, kept))
    assembly.say('%s (%s, chains %s/%s): fitted on %d of %d %s C-alpha '
                 'atoms, RMSD %.2f A' % (label, obj, heavy, light, kept,
                                         total, what, rmsd))
    result = []
    for residues, mapping in ((h, map_h), (l, map_l)):
        moved = []
        for i, r in enumerate(residues):
            r = r.moved(rotation, shift)
            r.part = part
            r.template = mapping.get(i)
            r.light = residues is l
            moved.append(r)
        result.append(moved)
    groups = []
    for chain in (heavy, light):
        found = []
        if glycans:
            for g in attached_groups(obj, chain, state, _self):
                g = g.moved(rotation, shift)
                g.part = part
                found.append(g)
        groups.append(found)
    return result[0], result[1], groups[0], groups[1]


def _ca_centroid(residues):
    xyz = [r.ca for r in residues if r.ca is not None]
    return np.mean(xyz, axis=0) if xyz else None


def _overlaps(xyz, tree, cutoff):
    return sum(len(found) for found in tree.query_ball_point(xyz, cutoff))


def _place_single(label, part, fmt, obj, chain, template, t_heavy, t_light,
                  state, assembly, glycans, _self, obstacles=(),
                  position='fv', orientation='',
                  scfv_linker=DEFAULT_SCFV_LINKER):
    '''
    An scFv or VHH: its variable domains superposed on the template arm's
    (position "fv": where a Fab's Fv is), or then moved toward the hinge
    (position "hinge"), into the place of the template's constant domains
    (CH1, CL) as the Fv of an scFv-Fc sits next to the hinge, or as far as
    it goes without overlapping the obstacles (coordinates of the parts
    placed already); residues after its last variable domain (tags) are
    left out

    :return: residues (template indices not set), het groups
    '''
    residues = read_chain(obj, chain, state, _self)
    if not residues:
        raise ValueError('%s: no protein chain "%s" in "%s"' % (
            label, chain, obj))
    domains = v_domains(_sequence(residues))
    if domains and domains[-1][2] < len(residues):
        tail = residues[domains[-1][2]:]
        assembly.say('%s: %d residues after its last variable domain left '
                     'out (%s)' % (label, len(tail), _sequence(tail)))
        residues = residues[:domains[-1][2]]
    kinds = [kind for kind, _, _ in domains]
    if fmt == 'scfv' and not ('heavy' in kinds and 'light' in kinds):
        raise ValueError('%s: chain %s of "%s" needs a VH and a VL domain '
                         '(found: %s)' % (label, chain, obj,
                                          ', '.join(kinds) or 'none'))
    if fmt == 'vhh' and not domains:
        raise ValueError('%s: no variable domain in chain %s of "%s"' % (
            label, chain, obj))
    if fmt == 'scfv' and orientation:
        residues, order = reorder_scfv(residues, domains, orientation,
                                       scfv_linker)
        if order.get('to'):
            assembly.say('%s: %s in the structure, reordered to %s (its own '
                         'linker %s left out, the domains joined by %s, '
                         'built in the refinement)' % (
                             label, order['from'], order['to'],
                             order['left_out'] or 'none',
                             order['linker'] or 'nothing'))
            domains = v_domains(_sequence(residues))
    th = read_chain(template, t_heavy, -1, _self)
    tl = read_chain(template, t_light, -1, _self)
    h_end = constant_start(_sequence(th))
    l_end = constant_start(_sequence(tl))
    c_end = ch1_end(_sequence(th)) or h_end + 100
    targets = {'heavy': th[:h_end], 'light': tl[:l_end]}
    mobile, target = [], []
    used = []
    for kind, start, end in domains[:2] if fmt == 'scfv' else domains[:1]:
        mapping = _map(targets[kind], residues[start:end])
        m, t = _ca_pairs(residues[start:end], mapping, targets[kind])
        mobile += m
        target += t
        used.append(kind)
        if kind == 'light':
            for r in residues[start:end]:
                r.light = True
    if len(mobile) < MIN_PAIRS:
        raise ValueError('%s: only %d residues of its variable domains '
                         'match the template\'s' % (label, len(mobile)))
    rotation, shift, rmsd, kept = superpose(mobile, target)
    # from the template's V domains to its constant ones (of the kinds
    # used: VH and VL, or VH alone)
    v_parts = [targets[kind] for kind in used]
    c_parts = [th[h_end:c_end + 1] if kind == 'heavy' else tl[l_end:]
               for kind in used]
    move = _ca_centroid([r for p in c_parts for r in p]) - \
        _ca_centroid([r for p in v_parts for r in p])
    full = float(np.linalg.norm(move))
    fraction, note = 1.0, ''
    if position == 'fv':
        fraction, move = 0.0, move * 0.0
    elif len(obstacles):
        from scipy.spatial import cKDTree
        tree = cKDTree(np.asarray(obstacles))
        xyz = np.array([xyz for r in residues
                        for _, xyz, _ in r.atoms.values()])
        xyz = xyz @ rotation.T + shift
        best = None
        for step in range(10, -1, -1):
            n = _overlaps(xyz + move * step / 10.0, tree, SINGLE_CLEARANCE)
            if best is None or n < best[1]:
                best = (step / 10.0, n)
            if n <= SINGLE_OVERLAPS:
                break
        fraction = best[0]
        if fraction < 1.0:
            note = ' (of %.0f: further would overlap the other parts%s)' % (
                full, '; the other arm of the template, or another '
                'template, may suit it better' if fraction < 0.5 else '')
        move = move * fraction
    assembly.fits.append((label, rmsd, kept))
    where = 'on the template\'s Fv' if position == 'fv' else \
        'then moved %.0f A toward the hinge%s' % (full * fraction, note)
    assembly.say('%s (%s, chain %s, %s): fitted on %d of %d variable-domain '
                 'C-alpha atoms, RMSD %.2f A, %s' % (
                     label, obj, chain, scfv_order(domains) if fmt == 'scfv'
                     else 'VHH',
                     kept, len(mobile), rmsd, where))
    placed = []
    for r in residues:
        light, own = r.light, r.part == 'linker'   # (a reordered scFv's)
        r = r.moved(rotation, shift + move)
        if not own:
            r.part = part
        r.light = light
        placed.append(r)
    groups = []
    if glycans:
        for g in attached_groups(obj, chain, state, _self):
            g = g.moved(rotation, shift + move)
            g.part = part
            groups.append(g)
    return placed, groups


def _place_fc(obj, chains, template, t_heavy, state, assembly, glycans,
              _self):
    '''
    The Fc superposed on the template's (both chain orders tried)

    :return: [residues for heavy 1, for heavy 2], [groups, groups]
    '''
    fc = [read_chain(obj, c, state, _self) for c in chains]
    th = [read_chain(template, c, -1, _self) for c in t_heavy]
    if not fc[0] or not fc[1]:
        raise ValueError('Fc: no protein chains %s in "%s"' % (
            '/'.join(chains), obj))
    best = None
    orders = [(0, 1), (1, 0)] if chains[0] != chains[1] else [(0, 1)]
    for order in orders:
        maps = [_map(th[k], fc[order[k]]) for k in range(2)]
        mobile, target = [], []
        for k in range(2):
            m, t = _ca_pairs(fc[order[k]], maps[k], th[k])
            mobile += m
            target += t
        if len(mobile) < MIN_PAIRS:
            continue
        fit = superpose(mobile, target)
        if best is None or fit[2] < best[0][2]:
            best = (fit, order, maps, len(mobile))
    if best is None:
        raise ValueError('Fc: too few residues match the template\'s Fc')
    (rotation, shift, rmsd, kept), order, maps, total = best
    assembly.fits.append(('Fc', rmsd, kept))
    assembly.say('Fc (%s, chains %s): fitted on %d of %d C-alpha atoms, '
                 'RMSD %.2f A' % (obj, '/'.join(chains[k] for k in order),
                                  kept, total, rmsd))
    result, groups = [], []
    for k in range(2):
        moved = []
        for i, r in enumerate(fc[order[k]]):
            r = r.moved(rotation, shift)
            r.part = 'fc'
            r.template = maps[k].get(i)
            moved.append(r)
        result.append(moved)
        found = []
        if glycans:
            for g in attached_groups(obj, chains[order[k]], state, _self):
                g = g.moved(rotation, shift)
                g.part = 'fc'
                found.append(g)
        groups.append(found)
    return result, groups


def _hinge_residues(k, fab, fc, template, p1, p2, hinge, assembly, label):
    '''
    The hinge between template positions p1 (Fab end) and p2 (Fc start):
    the user's residues (structure, or sequence) or the template's

    :param hinge: None, ('sequence', text) or ('structure', [Residue])
    '''
    if hinge is None:
        residues = []
        for t in range(p1 + 1, p2):
            r = Residue(template[t].resn, dict(template[t].atoms), 'hinge')
            r.template = t
            r.threaded = bool(r.atoms)
            residues.append(r)
        return residues
    kind, value = hinge
    if kind == 'sequence':
        residues = [Residue(_three(c), part='hinge') for c in value]
    else:
        residues = [Residue(r.resn, dict(r.atoms), 'hinge') for r in value]
    start = max(0, p1 - HINGE_WINDOW)
    end = min(len(template), p2 + HINGE_WINDOW)
    mapping = _map(template, residues, start, end)
    for i, r in enumerate(residues):
        r.template = mapping.get(i)
    # residues overlapping the Fab or the Fc are theirs
    first, last = 0, len(residues)
    for i, r in enumerate(residues):
        if r.template is not None and r.template <= p1:
            first = i + 1
    for i in range(len(residues) - 1, -1, -1):
        r = residues[i]
        if r.template is not None and r.template >= p2:
            last = i
    if first >= last and any(r.template is not None for r in residues):
        assembly.say('%s: the hinge given overlaps the Fab and the Fc '
                     'entirely; none of it added' % label)
        return []
    kept = residues[first:last]
    if kind == 'structure' and any(r.atoms for r in kept):
        # superposed on the template's hinge and the joined ends
        fab_at = {r.template: r for r in fab if r.template is not None}
        fc_at = {r.template: r for r in fc if r.template is not None}
        mobile, target = [], []
        for r in residues:
            t = r.template
            if t is None or r.ca is None:
                continue
            anchor = fab_at.get(t) if t <= p1 else \
                fc_at.get(t) if t >= p2 else template[t]
            if anchor is not None and anchor.ca is not None:
                mobile.append(r.ca)
                target.append(anchor.ca)
        if len(mobile) >= 3:
            rotation, shift, rmsd, n = superpose(mobile, target, cycles=0)
            assembly.say('%s hinge: placed on %d C-alpha atoms, RMSD %.2f '
                         'A' % (label, n, rmsd))
            kept = [r.moved(rotation, shift) for r in kept]
            for r in kept:
                r.part = 'hinge'
        else:
            assembly.warn('%s hinge: it doesn\'t overlap the template\'s '
                          'hinge or the joined ends; its sequence is used, '
                          'with the template\'s coordinates' % label)
            for r in kept:
                r.atoms = {}
    return kept


def _reaches(a, b, steps):
    '''
    C-alpha atoms a and b can be "steps" residues apart (over more than 3
    residues: without being nearly stretched)
    '''
    limit = REACH * steps + REACH_SLACK if steps <= 3 else STRETCHED * steps
    return np.linalg.norm(a - b) <= limit


def _stitch_heavy(k, label, fab, fc, template, hinge, assembly, p1=None):
    '''
    Fab heavy chain (or scFv, VHH, nothing) + hinge + Fc chain, by template
    positions; p1: the template position the arm ends at, for arms not
    aligned to the template (scFv, VHH, none: the end of CH1)
    '''
    if p1 is None:
        fab_mapped = [i for i, r in enumerate(fab)
                      if r.template is not None]
        if not fab_mapped:
            raise ValueError('%s: its heavy chain doesn\'t match the '
                             'template' % label)
        last = fab_mapped[-1]
        p1 = fab[last].template
        tail = fab[last + 1:]
        if tail:
            assembly.say('%s: %d residues after the heavy chain\'s last '
                         'template-aligned one left out (%s)' % (
                             label, len(tail), _sequence(tail)))
        fab = fab[:last + 1]
    fc_after = [i for i, r in enumerate(fc)
                if r.template is not None and r.template > p1]
    if not fc_after:
        raise ValueError('%s: its heavy chain reaches into the Fc; give a '
                         'Fab (heavy chain up to the hinge)' % label)
    first = fc_after[0]
    p2 = fc[first].template
    if first:
        assembly.say('Fc chain for %s: %d residues before the junction '
                     'left out (%s)' % (label, first,
                                        _sequence(fc[:first])))
    fc = fc[first:]
    middle = _hinge_residues(k, fab, fc, template, p1, p2, hinge, assembly,
                             label)
    chain = fab + middle + fc
    # the span between the Fab's last and the Fc's first residues with
    # coordinates: missing ones get the template's
    last_fab = max((i for i, r in enumerate(fab) if r.atoms), default=-1)
    first_fc = len(fab) + len(middle) + min(
        (i for i, r in enumerate(fc) if r.atoms), default=len(fc))
    threaded = 0
    for i in range(last_fab + 1, first_fc):
        r = chain[i]
        if r.atoms or r.template is None or not template[r.template].atoms:
            continue
        new = _thread(template[r.template], r.resn, r.part)
        if new.atoms:
            new.template = r.template
            chain[i] = new
            threaded += 1
    # template coordinates the arm's end can't reach (an arm placed
    # differently from the template's, e.g. an scFv) are left to build
    end = chain[last_fab].ca if last_fab >= 0 else None
    for i in range(last_fab + 1, first_fc):
        r = chain[i]
        if end is None or r.ca is None:
            continue
        if not r.threaded or _reaches(r.ca, end, i - last_fab):
            break
        r.atoms = {}
        r.threaded = False
    # ... nor the Fc's first residue with its own coordinates
    end = chain[first_fc].ca if first_fc < len(chain) else None
    for i in range(first_fc - 1, last_fab, -1):
        r = chain[i]
        if end is None or r.ca is None:
            continue
        if not r.threaded or _reaches(r.ca, end, first_fc - i):
            break
        r.atoms = {}
        r.threaded = False
    span = chain[last_fab + 1:first_fc]
    assembly.say('%s heavy chain: %d arm + %d hinge + %d Fc residues; '
                 'junction %s (%d residues, %d with template coordinates, '
                 '%d to build)' % (
                     label, len(fab), len(middle), len(fc),
                     _sequence(span) or '-', len(span),
                     sum(1 for r in span if r.threaded),
                     sum(1 for r in span if not r.atoms)))
    free = range(max(0, last_fab + 1 - FREE_FLANK),
                 min(len(chain), first_fc + FREE_FLANK))
    return chain, free


def _unthread_clashes(assembly, cutoff=CLASH):
    '''
    Residues with the template's coordinates which clash with the parts'
    own atoms lose them (refinement builds them instead): the template's
    hinge needn't fit the parts' ends
    '''
    own, own_ids, threaded = [], [], []
    for chain, residues in assembly.chains.items():
        for k, r in enumerate(residues):
            for name, (elem, xyz, _) in r.atoms.items():
                if r.threaded:
                    threaded.append((chain, k, name, elem, xyz))
                else:
                    own.append(xyz)
                    own_ids.append((chain, k, name, elem))
    for chain, groups in assembly.groups.items():
        for r in groups:
            for name, (elem, xyz, _) in r.atoms.items():
                own.append(xyz)
                own_ids.append((chain, None, name, elem))
    if not own or not threaded:
        return
    own = np.array(own)
    clashing = set()
    for chain, k, name, elem, xyz in threaded:
        if (chain, k) in clashing:
            continue
        d = np.linalg.norm(own - xyz, axis=1)
        for j in np.nonzero(d < cutoff)[0]:
            c, kk, other, other_elem = own_ids[j]
            if c == chain and kk is not None and (
                    (kk == k + 1 and name == 'C' and other == 'N') or
                    (kk == k - 1 and name == 'N' and other == 'C')):
                continue  # the peptide bond
            if elem.upper() == 'S' and other_elem.upper() == 'S':
                continue  # a disulfide
            clashing.add((chain, k))
            break
    for chain, k in clashing:
        r = assembly.chains[chain][k]
        r.atoms = {}
        r.threaded = False
    if clashing:
        assembly.say('Template residues clashing with the parts, built '
                     'instead: %s' % ', '.join(
                         '%s/%d' % (c, k + 1) for c, k in sorted(clashing)))


def _report_overlaps(assembly):
    '''
    Warn where different parts overlap (atoms closer than OVERLAP), e.g.
    an arm placed where the template's other parts are
    '''
    from scipy.spatial import cKDTree
    atoms, xyz = [], []
    for chain, residues in assembly.chains.items():
        for k, r in enumerate(residues):
            for name, (elem, coord, _) in r.atoms.items():
                atoms.append((chain, k, r.part, elem.upper()))
                xyz.append(coord)
    if not xyz:
        return
    counts = {}
    for i, j in cKDTree(np.array(xyz)).query_pairs(OVERLAP):
        a, b = atoms[i], atoms[j]
        if a[2] == b[2] or (a[0] == b[0] and abs(a[1] - b[1]) <= 1) or \
                a[3] == b[3] == 'S':
            continue
        key = tuple(sorted((a[2], b[2])))
        counts[key] = counts.get(key, 0) + 1
    for (a, b), n in sorted(counts.items()):
        text = '%s and %s overlap: %d atom pairs closer than %g A' % (
            assembly.labels[a], assembly.labels[b], n, OVERLAP)
        if n < OVERLAP_WARN:
            assembly.say(text + ' (refinement resolves that)')
        elif a.startswith('fusion') or b.startswith('fusion'):
            assembly.warn(text + ' (refinement may not resolve it; a longer '
                          'linker or another site may help)')
        else:
            assembly.warn(text + ' (refinement may not resolve it; another '
                          'template may fit better)')


def _more_linker(distance, steps):
    '''
    "; about N more residues (a linker) would relax it"
    '''
    more = int(np.ceil(distance / RELAXED)) - steps
    return '; about %d more residues (e.g. a linker) would relax it' % \
        more if more > 0 else ''


def _check_gaps(chain_id, residues, assembly):
    '''
    Warn where residues without coordinates can't bridge their gap, or
    consecutive residues are far apart
    '''
    previous = None
    for k, r in enumerate(residues):
        if r.ca is None:
            continue
        if previous is not None:
            j, ca = previous
            d = float(np.linalg.norm(r.ca - ca))
            steps = k - j
            if steps == 1 and d > 4.3:
                assembly.say('Chain %s: residues %d-%d are %.1f A apart '
                             '(refinement joins them)' % (
                                 chain_id, j + 1, k + 1, d))
            elif steps > 1 and d > CA_STEP * steps + 0.3:
                assembly.warn('chain %s: %d residues can\'t span the %.1f A '
                              'between residues %d and %d%s' % (
                                  chain_id, steps - 1, d, j + 1, k + 1,
                                  _more_linker(d, steps)))
            elif steps > 3 and d > STRETCHED * steps:
                assembly.warn('chain %s: the %d residues between residues '
                              '%d and %d span %.1f A, nearly stretched%s' % (
                                  chain_id, steps - 1, j + 1, k + 1, d,
                                  _more_linker(d, steps)))
        previous = (k, r.ca)


def parse_arms(text):
    '''
    "fab+scfv" (or "fab,scfv", "Fab/scFv") -> ["fab", "scfv"]
    '''
    parts = [p.strip().lower() for p in re.split(r'[+,/\s]+',
                                                 str(text or 'fab+fab'))
             if p.strip()]
    if len(parts) == 1:
        parts *= 2
    if len(parts) != 2 or not set(parts) <= set(ARM_FORMATS):
        raise ValueError('arms: two of %s, e.g. fab+scfv' %
                         ', '.join(ARM_FORMATS))
    if parts == ['none', 'none']:
        raise ValueError('arms: at least one arm')
    return parts


def arm_labels(formats):
    '''
    {"fab1": "Fab1", "fab2": "Fab2"}, or "Fab", "scFv" when they differ
    '''
    same = formats[0] == formats[1]
    return {'fab%d' % (k + 1): 'Arm %d (none)' % (k + 1) if f == 'none'
            else FORMAT_LABELS[f] + ('%d' % (k + 1) if same else '')
            for k, f in enumerate(formats)}


def assemble(fab1, fc, fab2='', template=DEFAULT_TEMPLATE, hinge='',
             hinge_sequence='', fab1_chains='', fab2_chains='',
             fc_chains='', hinge_chains='', template_chains='',
             chains='HLKM', glycans=True, state=-1, arms='fab+fab',
             fusions='', single_position='fv', arm_linkers=('', ''),
             arm_orientations=('', ''),
             arm_scfv_linkers=(DEFAULT_SCFV_LINKER,) * 2,
             fusion_linkers='relaxed',
             _self=cmd):
    '''
    Place and stitch the parts (see the module's description); fab1,
    fab2: the arms' structures, of the formats "arms"; fusions: units
    appended to chain C-termini (parse_fusions), their linkers "relaxed"
    (the unit near the chain's end), "extended" (as far out as the linker
    reaches) or "collapsed" (extended, then drawn in by dynamics) unless
    the fusion says so itself; single_position:
    scFv and VHH arms "fv" (on the template's Fv) or "hinge" (next to it);
    arm_linkers: residues between an scFv or VHH arm and the hinge

    :return: Assembly (not loaded)
    '''
    formats = parse_arms(arms)
    fusions = parse_fusions(fusions)
    if single_position not in SINGLE_POSITIONS:
        raise ValueError('scFv/VHH position: fv or hinge')
    if fusion_linkers not in LINKER_SHAPES:
        raise ValueError('fusion linkers: relaxed or extended')
    for fusion in fusions:
        fusion.setdefault('linker_shape', fusion_linkers)
    out = parse_chains(chains, 4)
    if len(set(out)) != 4:
        raise ValueError('the model\'s four chains need different IDs')
    assembly = Assembly()
    template = resolve(template or DEFAULT_TEMPLATE, _self)
    t_chains = parse_chains(template_chains, 4)
    if not all(t_chains):
        t_chains = detect_template(template, _self)
    assembly.say('Template %s: Fab1 chains %s/%s, Fab2 chains %s/%s' % (
        template, *t_chains))
    labels = arm_labels(formats)
    assembly.formats = dict(zip(('fab1', 'fab2'), formats))
    assembly.labels.update(labels)
    objects = [resolve(fab1, _self) if formats[0] != 'none' else '',
               resolve(fab2, _self) if formats[1] != 'none' else '']
    if formats[1] != 'none' and not objects[1] and formats[0] == formats[1]:
        objects[1] = objects[0]     # (the same structure on both arms)
        fab2_chains = fab2_chains or fab1_chains
    fc = resolve(fc, _self)
    if not fc:
        raise ValueError('give the Fc')
    fabs = []
    for k, (fmt, obj, spec) in enumerate(zip(formats, objects,
                                             (fab1_chains, fab2_chains))):
        label = labels['fab%d' % (k + 1)]
        if fmt == 'none':
            fabs.append((label, fmt, '', ()))
            continue
        if not obj:
            raise ValueError('give the %s structure' % label)
        if fmt == 'fab':
            found = parse_chains(spec, 2)
            if not all(found):
                found = detect_fab(obj, _self)
            if not all(found):
                raise ValueError('%s: no heavy and light chain found in '
                                 '"%s"; give them' % (label, obj))
        else:
            found = parse_chains(spec, 1)
            if not found[0]:
                found = [detect_single(obj, fmt, _self)]
            if not found[0]:
                raise ValueError('%s: no %s chain found in "%s"; give it' % (
                    label, FORMAT_LABELS[fmt], obj))
        fabs.append((label, fmt, obj, tuple(found)))
    fc_spec = parse_chains(fc_chains, 2)
    if not all(fc_spec):
        fc_spec = detect_fc(fc, _self)
    hinge_parts = [None, None]
    hinge = resolve(hinge, _self)
    sequence = re.sub(r'[^A-Za-z]', '', str(hinge_sequence)).upper()
    if hinge:
        h_spec = parse_chains(hinge_chains, 2) if hinge_chains else \
            detect_hinge(hinge, _self)
        for k in range(2):
            hinge_parts[k] = ('structure',
                              read_chain(hinge, h_spec[k], state, _self))
        assembly.inputs.append(hinge)
    elif sequence:
        bad = set(sequence) - set('ACDEFGHIKLMNPQRSTVWY')
        if bad:
            raise ValueError('hinge sequence: unknown letters %s' %
                             ''.join(sorted(bad)))
        hinge_parts = [('sequence', sequence)] * 2

    t_res = {c: read_chain(template, c, -1, _self) for c in t_chains}
    fc_res, fc_groups = _place_fc(fc, fc_spec, template,
                                  (t_chains[0], t_chains[2]), state,
                                  assembly, glycans, _self)

    def coordinates(*chains):
        return [xyz for residues in chains for r in residues
                for _, xyz, _ in r.atoms.values()]

    # Fabs first: scFv and VHH arms keep clear of them and the Fc
    placed = [([], [], [], []), ([], [], [], [])]
    for k, (label, fmt, obj, found) in enumerate(fabs):
        if fmt == 'fab':
            obstacles = coordinates(*fc_res, *placed[0][:2])
            placed[k] = _place_fab(label, 'fab%d' % (k + 1), obj, found[0],
                                   found[1], template, t_chains[2 * k],
                                   t_chains[2 * k + 1], state, assembly,
                                   glycans, _self, obstacles)
    for k, (label, fmt, obj, found) in enumerate(fabs):
        if fmt in ('scfv', 'vhh'):
            obstacles = coordinates(*fc_res, *placed[0][:2],
                                    *placed[1][:2])
            residues, groups = _place_single(
                label, 'fab%d' % (k + 1), fmt, obj, found[0], template,
                t_chains[2 * k], t_chains[2 * k + 1], state, assembly,
                glycans, _self, obstacles, single_position,
                parse_orientation(arm_orientations[k]
                                  if arm_orientations else ''),
                parse_linker(arm_scfv_linkers[k] if arm_scfv_linkers
                             else DEFAULT_SCFV_LINKER))
            linker = parse_linker(arm_linkers[k] if arm_linkers else '')
            residues = residues + [Residue(_three(c), part='linker')
                                   for c in linker]
            if linker:
                assembly.say('%s: linker %s to the hinge' % (label, linker))
            placed[k] = (residues, [], groups, [])
    assembly.heavy = [out[0], out[2]]
    assembly.light = [out[2 * k + 1] for k in range(2)
                      if formats[k] == 'fab']
    for k in range(2):
        label = labels['fab%d' % (k + 1)]
        h, l, h_groups, l_groups = placed[k]
        p1 = None
        if formats[k] != 'fab':
            p1 = ch1_end(_sequence(t_res[t_chains[2 * k]]))
            if p1 is None:
                raise ValueError('the end of CH1 (KVDKKV) wasn\'t found in '
                                 'template chain %s' % t_chains[2 * k])
        chain, free = _stitch_heavy(k, label, h, fc_res[k],
                                    t_res[t_chains[2 * k]], hinge_parts[k],
                                    assembly, p1)
        assembly.chains[out[2 * k]] = chain
        assembly.free |= {(out[2 * k], str(i + 1)) for i in free}
        # a reordered scFv's own linker, inside the chain, moves too
        for i, r in enumerate(chain):
            if r.part == 'linker':
                assembly.free |= {
                    (out[2 * k], str(j + 1)) for j in
                    range(max(0, i - FREE_FLANK),
                          min(len(chain), i + 1 + FREE_FLANK))}
        assembly.groups[out[2 * k]] = h_groups + fc_groups[k]
        if l:
            assembly.chains[out[2 * k + 1]] = l
            assembly.groups[out[2 * k + 1]] = l_groups
    sites = {'heavy1': out[0], 'light1': out[1], 'heavy2': out[2],
             'light2': out[3]}
    for fusion in sorted(fusions, key=lambda f: list(FUSION_SITES).index(
            f['site'])):
        chain_id = sites[fusion['site']]
        k = int(fusion['site'][-1])
        if fusion['site'].startswith('light'):
            if chain_id not in assembly.chains:
                raise ValueError('fusion on light chain %d: arm %d has no '
                                 'light chain (%s)' % (k, k, formats[k - 1]))
            owner = assembly.chains[chain_id] + [
                r for r in assembly.chains[sites['heavy%d' % k]]
                if r.part == 'fab%d' % k]
        else:
            owner = [r for c in (out[0], out[2])
                     for r in assembly.chains[c] if r.part == 'fc']
        _place_fusion(assembly, fusion, chain_id, owner, state, _self)
    _unthread_clashes(assembly)
    for chain_id, residues in assembly.chains.items():
        _check_gaps(chain_id, residues, assembly)
    _report_overlaps(assembly)
    assembly.inputs += [template] + [f[2] for f in fabs if f[2]] + [fc]
    assembly.inputs = list(dict.fromkeys(assembly.inputs))
    return assembly


# C-terminal fusions


def parse_linker(text):
    '''
    "GGGGSGGGGS", "(G4S)2", "(GGGGS)3", "GGGGS x 2" -> the sequence;
    parts can be mixed ("GSPG (EAAAK)2 GS"), spaces separate them (so
    "GSPG EAAAK x2" repeats EAAAK alone)
    '''
    text = str(text or '').upper()
    text = re.sub(r'G(\d)S', lambda m: 'G' * int(m.group(1)) + 'S', text)
    text = re.sub(r'\(([A-Z]+)\)\s*[X*]?\s*(\d+)',
                  lambda m: m.group(1) * int(m.group(2)), text)
    text = re.sub(r'([A-Z]+)\s*[X*]\s*(\d+)',
                  lambda m: m.group(1) * int(m.group(2)), text)
    text = text.replace(' ', '')
    bad = set(text) - set('ACDEFGHIKLMNPQRSTVWY')
    if bad:
        raise ValueError('linker "%s": unknown letters %s' % (
            text, ''.join(sorted(bad))))
    return text


def parse_orientation(text):
    '''
    "VH-VL", "vlvh", "VL/VH" or "" (the structure's order) -> "vh-vl",
    "vl-vh" or ""
    '''
    text = re.sub(r'[^a-z]', '', str(text or '').lower())
    if text in ('', 'asis', 'asgiven', 'structure'):
        return ''
    if text in ('vhvl', 'hl'):
        return 'vh-vl'
    if text in ('vlvh', 'lh'):
        return 'vl-vh'
    raise ValueError('scFv orientation: vh-vl or vl-vh')


def parse_fusions(text):
    '''
    "heavy:2ghw/B:(G4S)2 light1:vhh_model heavy1:2ghw/B:(G4S)2:vl-vh"
    (site:object[/chain][:linker[:orientation[:scfv_linker]]]; sites
    heavy1, heavy2, light1, light2, or heavy / light for both; linker
    default (G4S)2; an scFv's orientation vh-vl or vl-vh, default: the
    structure's, its domains then joined by scfv_linker, default (G4S)3)

    :return: [{'site', 'object', 'chain', 'linker', 'orientation',
        'scfv_linker'}]
    '''
    if isinstance(text, (list, tuple)):
        return [dict(f) for f in text]
    found = []
    for token in re.split(r'[\s|;]+', str(text or '').strip()):
        if not token:
            continue
        parts = token.split(':')
        if not 2 <= len(parts) <= 5:
            raise ValueError('fusion "%s": write site:object[/chain]'
                             '[:linker[:orientation[:scfv_linker]]], e.g. '
                             'heavy:2ghw/B:(G4S)2:vl-vh' % token)
        site = parts[0].lower()
        sites = FUSION_SITE_ALIASES.get(site, (site,))
        if not set(sites) <= set(FUSION_SITES):
            raise ValueError('fusion "%s": the site is one of %s, heavy, '
                             'light' % (token, ', '.join(FUSION_SITES)))
        obj, _, chain = parts[1].partition('/')
        linker = parse_linker(parts[2]) if len(parts) >= 3 else \
            DEFAULT_LINKER
        orientation = parse_orientation(parts[3] if len(parts) >= 4 else '')
        scfv_linker = parse_linker(parts[4]) if len(parts) >= 5 else \
            DEFAULT_SCFV_LINKER
        for one in sites:
            found.append({'site': one, 'object': obj, 'chain': chain,
                          'linker': linker, 'orientation': orientation,
                          'scfv_linker': scfv_linker, 'cdrs': 'out'})
    return found


def _rotation_between(a, b):
    '''
    Rotation matrix turning direction a onto direction b
    '''
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-8:
        if c > 0:
            return np.eye(3)
        # (half turn about any axis perpendicular to a)
        axis = np.cross(a, [1.0, 0, 0] if abs(a[0]) < 0.9 else [0, 1.0, 0])
        return _axis_rotation(axis, np.pi)
    k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + k + k @ k / (1 + c)


def _axis_rotation(axis, angle):
    axis = axis / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


def _fusion_unit(fusion, state, _self):
    '''
    The fused unit's residues (variable domains only: tags before and after
    left out) and its format ("scfv" or "vhh")
    '''
    obj = resolve(fusion['object'], _self)
    if not obj:
        raise ValueError('fusion on %s: give its structure' %
                         FUSION_SITES[fusion['site']])
    chain = fusion.get('chain') or detect_single(obj, 'scfv', _self) or \
        detect_single(obj, 'vhh', _self)
    if not chain:
        raise ValueError('fusion on %s: no scFv or VHH chain in "%s"' % (
            FUSION_SITES[fusion['site']], obj))
    residues = read_chain(obj, chain, state, _self)
    domains = v_domains(_sequence(residues))
    if not domains:
        raise ValueError('fusion on %s: no variable domain in chain %s of '
                         '"%s"' % (FUSION_SITES[fusion['site']], chain, obj))
    fmt = 'scfv' if len(domains) >= 2 else 'vhh'
    for kind, start, end in domains:
        for r in residues[start:end]:
            r.light = kind == 'light'
    residues, order = reorder_scfv(
        residues, domains, fusion.get('orientation', '') if fmt == 'scfv'
        else '', fusion.get('scfv_linker', DEFAULT_SCFV_LINKER))
    return obj, chain, fmt, residues, order


def cdr3_range(sequence, start, end):
    '''
    The CDR3 of a variable domain (v_domains range): after the last
    cysteine before the framework 4 motif, up to it

    :return: (first, last + 1) in the sequence
    '''
    found = list(re.finditer(r'[WFL]G.G', sequence[start:end]))
    fr4 = start + found[-1].start() if found else max(start, end - 11)
    cys = sequence.rfind('C', start, fr4)
    return (cys + 3, fr4) if cys > start else (max(start, fr4 - 10), fr4)


def cdr_direction(residues):
    '''
    Where a variable domain or an Fv (one or two domains) has its binding
    site: the direction from its center to the center of its CDR3 loops,
    the face the other CDRs are on too

    :return: unit vector, or None without coordinates
    '''
    sequence = _sequence(residues)
    loops, body = [], []
    for _, start, end in v_domains(sequence)[:2]:
        first, last = cdr3_range(sequence, start, end)
        loops += [r.ca for r in residues[first:last] if r.ca is not None]
        body += [r.ca for r in residues[start:end] if r.ca is not None]
    if not loops or not body:
        return None
    direction = np.mean(loops, axis=0) - np.mean(body, axis=0)
    norm = np.linalg.norm(direction)
    return direction / norm if norm > 1e-6 else None


def _angle_between(a, b):
    '''
    Degrees between two vectors
    '''
    cosine = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def scfv_order(domains):
    '''
    "VH-VL" or "VL-VH" of an scFv's domains (v_domains)
    '''
    return '-'.join('VH' if kind == 'heavy' else 'VL'
                    for kind, _, _ in domains[:2])


def reorder_scfv(residues, domains, orientation, linker=DEFAULT_SCFV_LINKER):
    '''
    An scFv's domains in the order "orientation" ("vh-vl", "vl-vh", or ""
    for the structure's): the same Fv, its two domains joined the other way
    round. Their own linker can't span that, so it is left out and the
    domains are joined by "linker" (no coordinates: the refinement builds
    it). Residues before the first domain and after the last (tags) go.

    :return: residues, {'from', and when reordered 'to', 'left_out',
        'linker'}
    '''
    order = scfv_order(domains)
    if len(domains) < 2 or not orientation or \
            orientation == order.lower():
        return residues[domains[0][1]:domains[-1][2]], {'from': order}
    (s1, e1), (s2, e2) = [(_domain_start(residues, start, end), end)
                          for _, start, end in domains[:2]]
    joined = [Residue(_three(c), part='linker') for c in linker]
    return residues[s2:e2] + joined + residues[s1:e1], {
        'from': order, 'to': orientation.upper(),
        'left_out': _sequence(residues[e1:s2]),
        'linker': ''.join(r.letter for r in joined)}


def _domain_start(residues, start, end):
    '''
    The first residue of a variable domain: 22 before its framework 1 Cys
    (Cys22 of VH, Cys23 of VL in Kabat numbering), and not a glycine of a
    linker before it
    '''
    found = re.search(V_MOTIF, _sequence(residues[start:end]))
    first = start + max(0, found.start() - 22) if found else start
    while first < end and residues[first].resn == 'GLY':
        first += 1
    return first


def _place_fusion(assembly, fusion, chain_id, owner, state, _self):
    '''
    Append a unit (scFv, VHH) through a linker to chain chain_id: its
    N-terminus FUSION_GAP from the chain's last residue, pointing away from
    the chain's part ("owner" residues: the Fc for a heavy chain, the Fab
    for a light chain); of the directions within 75 degrees and turns
    about them, the one overlapping the model least (C-alpha atoms)
    '''
    from scipy.spatial import cKDTree
    site = fusion['site']
    obj, unit_chain, fmt, unit, order = _fusion_unit(fusion, state, _self)
    chain = assembly.chains[chain_id]
    last = max((k for k, r in enumerate(chain) if r.ca is not None),
               default=None)
    if last is None:
        raise ValueError('fusion on %s: chain %s has no coordinates' % (
            FUSION_SITES[site], chain_id))
    anchor = chain[last].ca
    between = len(chain) - 1 - last + len(fusion['linker'])
    centre = _ca_centroid(owner)
    outward = anchor - centre
    outward /= np.linalg.norm(outward)
    first = next(k for k, r in enumerate(unit) if r.ca is not None)
    ca = np.array([r.ca for r in unit if r.ca is not None])
    n_term = unit[first].ca
    body = ca.mean(axis=0) - n_term          # N-terminus to the unit's body
    # how far out the unit sits: within reach of a relaxed linker, or at
    # the length of an extended one (which it then has to stretch to)
    shape = fusion.get('linker_shape', 'relaxed')
    extended = shape in ('extended', 'collapsed')
    gap = EXTENDED_STEP * (between + 1) if extended else \
        min(FUSION_GAP, max(2.0, 0.8 * CA_STEP * (between + 1)))
    # (not the chain's last residues: the unit is joined to them)
    joined = {id(r) for r in chain[max(0, last - 3):last + 1]}
    others = np.array([r.ca for residues in assembly.chains.values()
                       for r in residues
                       if r.ca is not None and id(r) not in joined])
    tree = cKDTree(others)
    atoms = np.array([xyz for residues in assembly.chains.values()
                      for r in residues if id(r) not in joined
                      for _, xyz, _ in r.atoms.values()])
    atom_tree = cKDTree(atoms)
    unit_xyz = np.array([xyz for r in unit for _, xyz, _ in
                         r.atoms.values()])
    # directions around "outward"
    e1 = np.cross(outward, [1.0, 0, 0] if abs(outward[0]) < 0.9
                  else [0, 1.0, 0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(outward, e1)
    # the binding site, to be turned toward the solvent (away from the
    # rest of the model) rather than against the antibody
    cdrs = cdr_direction(unit) if fusion.get('cdrs', 'out') == 'out' \
        else None
    middle = others.mean(axis=0)

    def cdr_angle(rotation, shift):
        if cdrs is None:
            return 0.0
        centre = ca.mean(axis=0) @ rotation.T + shift
        return _angle_between(rotation @ cdrs, centre - middle)

    # C-alpha atoms first (fast), then the best poses with all atoms. The
    # N-terminus stays put whatever the rotation, so the unit can be laid
    # out along the direction it reaches (its body), or turned to face its
    # CDRs that way; both are tried and the clash-free ones compared.
    poses = []
    for tilt in np.radians((0, 25, 50, 75)):
        for azimuth in (np.radians(np.arange(0, 360, 45)) if tilt else (0,)):
            direction = np.cos(tilt) * outward + np.sin(tilt) * (
                np.cos(azimuth) * e1 + np.sin(azimuth) * e2)
            aligns = [_rotation_between(body, direction)]
            if cdrs is not None:
                aligns.append(_rotation_between(cdrs, direction))
            for align in aligns:
                for spin in np.radians(np.arange(0, 360, 30)):
                    rotation = _axis_rotation(direction, spin) @ align
                    shift = anchor + gap * direction - rotation @ n_term
                    n = _overlaps(ca @ rotation.T + shift, tree,
                                  FUSION_CLEARANCE)
                    poses.append(((n, cdr_angle(rotation, shift), tilt),
                                  rotation, shift))
    poses.sort(key=lambda pose: pose[0])
    best = None
    for (n, angle, tilt), rotation, shift in poses[:FUSION_POSES]:
        xyz = unit_xyz @ rotation.T + shift
        # real clashes first, then contacts, then where the CDRs face
        clashes = _overlaps(xyz, atom_tree, OVERLAP)
        close = _overlaps(xyz, atom_tree, SINGLE_CLEARANCE)
        score = (clashes, close, angle, n, tilt)
        if best is None or score < best[0]:
            best = (score, rotation, shift)
    (_, overlaps, angle, _, tilt), rotation, shift = best
    part = 'fusion_' + site
    linker = []
    for letter in fusion['linker']:
        r = Residue(_three(letter), part='linker')
        linker.append(r)
    placed = []
    for r in unit:
        moved = r.moved(rotation, shift)     # (keeps "light")
        if r.part != 'linker':               # (a reordered scFv's own)
            moved.part = part
        placed.append(moved)
    start = len(chain)
    chain.extend(linker + placed)
    # the linkers and their neighbors move in refinement
    free = set(range(max(0, last + 1 - FREE_FLANK),
                     min(len(chain), start + len(linker) + FREE_FLANK)))
    for k, r in enumerate(placed, start + len(linker)):
        if r.part == 'linker':
            free |= set(range(max(0, k - FREE_FLANK),
                              min(len(chain), k + 1 + FREE_FLANK)))
    assembly.free |= {(chain_id, str(k + 1)) for k in free}
    if shape == 'collapsed':
        # the dynamics pulls the unit's first residue back to the anchor
        assembly.collapse.append({
            'chain': chain_id, 'anchor': last + 1, 'linker': start + 1,
            'unit': start + len(linker) + 1, 'residues': len(linker)})
    if site not in assembly.fusion_sites:
        assembly.fusion_sites.append(site)
    assembly.labels[part] = '%s on %s' % (FORMAT_LABELS[fmt],
                                          FUSION_SITES[site])
    assembly.inputs.append(obj)
    if order.get('to'):
        assembly.say('%s (%s, chain %s): %s in the structure, reordered to '
                     '%s (its own linker %s left out, the domains joined by '
                     '%s, built in the refinement)' % (
                         FORMAT_LABELS[fmt], obj, unit_chain, order['from'],
                         order['to'], order['left_out'] or 'none',
                         order['linker'] or 'nothing'))
    assembly.say('%s (%s, chain %s%s) appended to the C-terminus of %s '
                 '(chain %s) through %s: pointing %.0f degrees from straight '
                 'out%s' % (FORMAT_LABELS[fmt], obj, unit_chain,
                            ', ' + order.get('to', order['from'])
                            if fmt == 'scfv' else '',
                            FUSION_SITES[site], chain_id,
                            ('%s, %s (%.0f A out)' % (
                                fusion['linker'],
                                'extended then collapsed by the dynamics'
                                if shape == 'collapsed' else 'extended',
                                gap)) if extended
                            else fusion['linker'] or 'no linker',
                            np.degrees(tilt),
                            ', %d atom pairs within %g A of the rest' % (
                                overlaps, SINGLE_CLEARANCE)
                            if overlaps else ''))
    if cdrs is not None:
        assembly.say('  its CDRs face %.0f degrees from straight out (%s)'
                     % (angle, 'into the solvent' if angle <= CDR_OUTWARD
                        else 'toward the antibody: every pose with them '
                        'further out would clash more'))


# disulfide bonds


def cysteine_pairs(assembly, cutoff=DISULFIDE_SEARCH):
    '''
    Candidate disulfide bonds of the model: cysteine pairs with SG atoms
    within cutoff, and pairs of cysteines at the same place in the hinge
    region (moved in refinement) of the two heavy chains (distance
    unknown when one is still to be built)

    :return: [{'a': (chain, resi), 'b': (chain, resi), 'distance': SG-SG
        A or None, 'parts': (part, part), 'bonded': forms by distance}],
        nearest first
    '''
    cys = []
    for chain, residues in assembly.chains.items():
        for k, r in enumerate(residues, 1):
            if r.resn == 'CYS':
                sg = r.atoms.get('SG')
                cys.append(((chain, str(k)), r.part,
                            None if sg is None else sg[1],
                            k if r.template is None else r.template))
    found = []
    for i, (a, part_a, sg_a, at_a) in enumerate(cys):
        for b, part_b, sg_b, at_b in cys[i + 1:]:
            # (hinge: the same place in both heavy chains)
            hinge = a[0] != b[0] and a in assembly.free and \
                b in assembly.free and at_a == at_b
            d = None
            if sg_a is not None and sg_b is not None:
                d = float(np.linalg.norm(sg_a - sg_b))
            if not hinge and (d is None or d > cutoff):
                continue
            found.append({'a': a, 'b': b, 'distance': d,
                          'parts': (part_a, part_b),
                          'bonded': d is not None and d <= DISULFIDE_BOND})
    # (between chains first: those are the ones to choose)
    found.sort(key=lambda x: (x['a'][0] == x['b'][0], x['distance'] is None,
                              x['distance'] or 0.0))
    return found


def pair_label(pair):
    return '%s/%s-%s/%s' % (pair[0] + pair[1])


def set_disulfides(name, pairs, _self=cmd):
    '''
    Without refinement: SG-SG bonds of the object made exactly the given
    pairs (bonds only; the atoms don't move)
    '''
    from .interfaces import _resi
    _self.unbond('%%%s and name SG' % name, '%%%s and name SG' % name)
    for (ca, ra), (cb, rb) in pairs:
        a = '%%%s and chain "%s" and resi %s and name SG' % (name, ca,
                                                              _resi(ra))
        b = '%%%s and chain "%s" and resi %s and name SG' % (name, cb,
                                                              _resi(rb))
        if _self.count_atoms(a) == 1 and _self.count_atoms(b) == 1:
            _self.bond(a, b, quiet=1)


# loading


def _atom_name(name, elem):
    if len(name) < 4 and len(elem) == 1:
        return ' %-3s' % name
    return '%-4s' % name[:4]


def pdb_text(assembly):
    '''
    PDB text of the residues with coordinates
    '''
    lines = []
    serial = 0

    def write(record, chain, resv, r):
        nonlocal serial
        for name, (elem, xyz, b) in r.atoms.items():
            serial += 1
            lines.append(
                '%-6s%5d %s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      '
                '    %2s' % (record, serial % 100000, _atom_name(name, elem),
                             r.resn[:3], chain, resv, xyz[0], xyz[1], xyz[2],
                             1.0, b, elem.upper()))

    for chain, residues in assembly.chains.items():
        for k, r in enumerate(residues, 1):
            write('ATOM', chain, k, r)
        lines.append('TER')
    for chain, groups in assembly.groups.items():
        for k, r in enumerate(groups, 1001):
            write('HETATM', chain, k, r)
    lines.append('END')
    return '\n'.join(lines) + '\n'


def object_name(name, _self=cmd):
    '''
    A legal object name; selection keywords (e.g. "model") get "_" as
    PyMOL gives them (found with a one-atom probe: the check isn't
    exposed)
    '''
    name = _self.get_legal_name(str(name).strip() or 'igg_model')
    names = _self.get_names('all')
    if name in names:
        return name
    _self.load_raw('HETATM    1  C   UNK     1       0.000   0.000   0.000'
                   '  1.00  0.00           C\nEND\n', 'pdb', name, quiet=1)
    new = [n for n in _self.get_names('all') if n not in names]
    for n in new:
        _self.delete(n)
    return new[0] if len(new) == 1 else name


def load(assembly, name, _self=cmd):
    '''
    The assembly as object "name" (residues without coordinates as CA
    atoms without coordinates, like PyMOL's mmCIF reader)

    :return: the object's name
    '''
    from . import missing
    _self.delete(name)
    before = set(_self.get_names('objects'))
    _self.load_raw(pdb_text(assembly), 'pdb', name, quiet=1)
    # (PyMOL renames objects named like selection keywords, e.g. "model")
    new = [n for n in _self.get_names('objects') if n not in before]
    name = new[0] if len(new) == 1 else name
    absent = [(chain, str(k), k, r.resn)
              for chain, residues in assembly.chains.items()
              for k, r in enumerate(residues, 1) if not r.atoms]
    if absent:
        missing.add_missing_residues(name, absent, _self=_self)
    set_parts(name, assembly.parts(), _self, assembly.tags())
    assembly.name = name
    return name


def set_parts(name, parts, _self=cmd, tags=PART_TAGS):
    _self.alter('%' + name,
                'custom = TAGS.get(parts.get((chain, resi)), custom)',
                space={'parts': parts, 'TAGS': tags})


def _residues_selection(keys):
    '''
    "(chain H and resi 1-120+130) or ..." for {(chain, resi)} (numbers)
    '''
    by_chain = {}
    for chain, resi in keys:
        by_chain.setdefault(chain, []).append(int(resi))
    parts = []
    for chain, numbers in sorted(by_chain.items()):
        numbers.sort()
        ranges = []
        start = previous = numbers[0]
        for n in numbers[1:] + [None]:
            if n is not None and n == previous + 1:
                previous = n
                continue
            ranges.append('%d-%d' % (start, previous) if previous > start
                          else '%d' % start)
            if n is not None:
                start = previous = n
        parts.append('(chain "%s" and resi %s)' % (chain, '+'.join(ranges)))
    return ' or '.join(parts) or 'none'


def _lighter(hex_color, amount):
    rgb = [int(hex_color[k:k + 2], 16) / 255.0 for k in (1, 3, 5)]
    return [c + (1.0 - c) * amount for c in rgb]


def color_parts(name, assembly, _self=cmd):
    '''
    Each part in its color, light chains and VL domains lighter; het
    groups by element with the part's carbons
    '''
    tags = assembly.tags()
    light = _residues_selection(assembly.light_residues())
    for key, hex_color in PART_COLORS.items():
        _self.set_color('ab_' + key, _lighter(hex_color, 0.0))
        _self.set_color('ab_%s_light' % key, _lighter(hex_color,
                                                      LIGHT_TINT))
    for part, tag in tags.items():
        key = 'fusion' if part.startswith('fusion') else \
            'hinge' if part == 'linker' else part
        _self.color('ab_' + key, '%%%s and custom %s' % (name, tag))
        _self.color('ab_%s_light' % key, '%%%s and custom %s and (%s)' % (
            name, tag, light))
    _self.color('atomic', '%%%s and not polymer and not elem C' % name)


def show_model(name, _self=cmd):
    _self.show_as('cartoon', '%' + name)
    _self.show('sticks', '%%%s and not polymer and not solvent' % name)
    _self.zoom('%' + name, animate=0.5)


# refinement


def refine_job(assembly, source, name, _self=cmd, **options):
    '''
    pymolx.prep job for the loaded assembly "source" into "name": build
    missing residues and side chains, hydrogens, minimization with the
    hinge and the junctions free

    :return: prep.PrepJob (run() it, then finish())
    '''
    from . import prep
    settings = dict(fill_loops=True, fill_termini=False,
                    fill_side_chains=True, delete_waters=True,
                    pka_method='standard', ligand_hydrogens=False,
                    diagnostics=False)
    settings.update(options)
    settings['free_residues'] = set(assembly.free)
    if assembly.disulfides is not None:
        settings['disulfide_pairs'] = list(assembly.disulfides)
    if assembly.collapse:
        # the fused units move with their linkers, drawn in from their
        # extended start by the dynamics and then left to settle
        mobile, pull, rigid = set(), [], []
        for one in assembly.collapse:
            chain = one['chain']
            end = len(assembly.chains[chain])
            # everything between the anchor and the chain's end moves:
            # the linker, the unit, and any residues of the chain's own
            # tail before them (built, so also stretched out at first)
            mobile |= {(chain, str(k))
                       for k in range(one['anchor'] + 1, end + 1)}
            # pulled by its first residues, so the force is shared, and
            # held in shape so the pull can't unpick its sheets
            pull.append(((chain, str(one['anchor']), 'CA'),
                         tuple((chain, str(k), 'CA') for k in
                               range(one['unit'], one['unit'] + PULL_SPREAD)),
                         COIL_STEP * max(1, one['unit'] - one['anchor']
                                         - 1) ** 0.5, prep.PULL_K))
            rigid.append([(chain, str(k), name)
                          for k in range(one['unit'], end + 1)
                          for name in BACKBONE])
        settings['free_residues'] = set(settings['free_residues']) | mobile
        settings['dynamics'] = {
            'mobile': mobile, 'pull': pull, 'rigid': rigid,
            'pull_ps': PULL_PS, 'settle_ps': SETTLE_PS,
            'temperature': 300.0, 'label': 'Collapsing the linkers'}
    return prep.prepare_job('%' + source, name, _self=_self, **settings)


def finish(assembly, job, source, _self=cmd):
    '''
    Load the refined model, delete the unrefined one, parts and colors
    '''
    job.apply(_self)
    assembly.name = job.name
    assembly.log += job.log
    flipped = inverted_residues(source, job.name, _self)
    if source != job.name:
        _self.delete(source)
    present(assembly, _self)
    if flipped:
        assembly.warn('refinement inverted residues (now D): %s' % ', '.join(
            '%s/%s' % key for key in flipped))
    stretched = long_peptide_bonds(job.name, _self)
    if stretched:
        assembly.warn('peptide bonds still longer than 1.5 A: %s' % ', '.join(
            '%s %s-%s %.1f A' % s for s in stretched[:8]))


def present(assembly, _self=cmd):
    name = assembly.name
    set_parts(name, assembly.parts(), _self, assembly.tags())
    color_parts(name, assembly, _self)
    show_model(name, _self)
    _self.disable(' '.join(o for o in assembly.inputs if o != name))


def _handedness(name, _self):
    from .diagnostics import _chirality
    atoms = {}
    _self.iterate_state(1, '%%%s and polymer.protein and name N+CA+C+CB' %
                        name, '_ = atoms.setdefault((chain, resi), {})'
                        '.setdefault(name, {"coord": (x, y, z)})',
                        space={'atoms': atoms})
    return {key: _chirality(named) > -0.5 for key, named in atoms.items()
            if len(named) == 4}


def inverted_residues(before, after, _self=cmd):
    '''
    :return: [(chain, resi)] L in "before" and D in "after"
    '''
    old = _handedness(before, _self)
    new = _handedness(after, _self)
    return sorted((key for key, left in old.items()
                   if left and new.get(key) is False),
                  key=lambda k: (k[0], int(k[1])))


def long_peptide_bonds(name, _self=cmd, limit=1.5):
    '''
    :return: [(chain, resi, next resi, C-N distance)] longer than limit
    '''
    atoms = {}
    _self.iterate_state(1, '%%%s and polymer.protein and name C+N' % name,
                        'atoms[(chain, resv, name)] = (x, y, z)',
                        space={'atoms': atoms})
    found = []
    for (chain, resv, atom), xyz in atoms.items():
        if atom != 'C':
            continue
        n = atoms.get((chain, resv + 1, 'N'))
        if n is not None:
            d = float(np.linalg.norm(np.subtract(xyz, n)))
            if d > limit:
                found.append((chain, resv, resv + 1, d))
    return sorted(found)


def build_antibody(fab1, fc, fab2='', template=DEFAULT_TEMPLATE, hinge='',
                   hinge_sequence='', name='igg_model', fab1_chains='',
                   fab2_chains='', fc_chains='', hinge_chains='',
                   template_chains='', chains='HLKM', glycans=1, refine=1,
                   disulfides='', arms='fab+fab', fusions='',
                   scfv_position='fv', linker1='', linker2='',
                   orientation1='', orientation2='',
                   fusion_linkers='relaxed',
                   forcefield='amber14', platform='auto', state=-1, quiet=0,
                   *, _self=cmd):
    '''
DESCRIPTION

    "build_antibody" builds a full IgG from two Fabs, the Fc and the
    hinge, placed on a full-length IgG template: each part is superposed
    on the template by sequence-matched C-alpha atoms, the heavy chains
    are joined through the hinge, and (refine=1) residues without
    coordinates are built and the model minimized with the hinge free.

USAGE

    build_antibody fab1, fc [, fab2 [, template [, hinge
        [, hinge_sequence [, name ]]]]]

ARGUMENTS

    fab1, fab2 = str: object, file or PDB ID of each arm: a Fab (heavy +
    light chain), an scFv or a VHH (see arms); fab2 empty: the same as
    fab1 (arms of the same format)

    fc = str: object, file or PDB ID of the Fc (two chains)

    template = str: full IgG (object, file or PDB ID) {default: 1hzh}

    hinge = str: object, file or PDB ID of the hinge (one chain per heavy
    chain, or one for both) {default: the template's hinge}

    hinge_sequence = str: hinge sequence instead (e.g. EPKSCDKTHTCPPCP)

    name = str: the new object {default: igg_model}

    fab1_chains, fab2_chains = str: heavy and light chain of a Fab, e.g.
    "HL", or the chain of an scFv or VHH {default: recognized from the
    sequences}

    fc_chains, hinge_chains = str: two chains, e.g. "AB"

    template_chains = str: heavy 1, light 1, heavy 2, light 2

    chains = str: the model's chains, same order {default: HLKM}

    glycans = 0/1: keep glycans bonded to the chains {default: 1}

    refine = 0/1: build missing residues and minimize (OpenMM)
    {default: 1}

    arms = str: the arms' formats: fab, scfv (VH and VL in one chain,
    placed next to the hinge), vhh (single domain) or none (one-armed),
    e.g. "fab+scfv" {default: fab+fab}

    scfv_position = fv or hinge: scFv and VHH arms on the template arm's
    Fv, or next to the hinge (in place of CH1/CL, as in scFv-Fc)
    {default: fv}

    linker1, linker2 = str: residues between an scFv or VHH arm and the
    hinge, e.g. (G4S)3 {default: none}

    orientation1, orientation2 = vh-vl or vl-vh: the domain order of an
    scFv arm, reordered if the structure's differs (the same Fv, its own
    linker replaced by (G4S)3) {default: the structure's}

    fusion_linkers = relaxed, extended or collapsed: whether a fused unit
    sits near the chain's end (its linker then coiled), as far out as the
    linker reaches when stretched, or starts there and is drawn in by a
    short dynamics run (OpenMM), which lets the linker find its own coil
    {default: relaxed}

    fusions = str: units (scFv, VHH) appended to chain C-termini through a
    linker: site:object[/chain][:linker], several separated by spaces;
    sites heavy1, heavy2, light1, light2, or heavy / light for both;
    linker e.g. (G4S)3 {default linker: (G4S)2}

    disulfides = str: exactly these disulfide bonds, in the model's
    numbering, e.g. "H/223-L/214, H/229-K/229" (chosen cysteines apart
    are brought together by the refinement), or "none" {default: those
    close enough}

EXAMPLE

    build_antibody 1n8z, 3ave, template=1hzh
    build_antibody fab_model, fc_model, hinge_sequence=ESKYGPPCPPCP
    build_antibody 1n8z, 3ave, scfv_model, arms=fab+scfv
    build_antibody 1n8z, 3ave, fusions=heavy:2ghw/B:(G4S)3

SEE ALSO

    prepare_structure, align
    '''
    from .prep import parse_disulfides
    pairs = parse_disulfides(disulfides)
    assembly = assemble(fab1, fc, fab2, template, hinge, hinge_sequence,
                        fab1_chains, fab2_chains, fc_chains, hinge_chains,
                        template_chains, chains, bool(int(glycans)), state,
                        arms, fusions, str(scfv_position).strip().lower(),
                        (linker1, linker2), (orientation1, orientation2),
                        fusion_linkers=str(fusion_linkers).strip().lower(),
                        _self=_self)
    assembly.disulfides = pairs
    name = object_name(name, _self)
    if int(refine):
        source = _self.get_unused_name('_ab_assembly')
        load(assembly, source, _self)
        try:
            job = refine_job(assembly, source, name, _self,
                             forcefield=forcefield, platform=platform)
            job.run()
            finish(assembly, job, source, _self)
        finally:
            _self.delete(source)
    else:
        load(assembly, name, _self)
        if pairs is not None:
            set_disulfides(assembly.name, pairs, _self)
        present(assembly, _self)
    if not int(quiet):
        print(' build_antibody: "%s"' % name)
        for line in assembly.log:
            print('  ' + line)
    return assembly


def extend(_self):
    _self.extend('build_antibody', build_antibody)
    _self.build_antibody = build_antibody
