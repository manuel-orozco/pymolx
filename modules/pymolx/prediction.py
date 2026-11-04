'''
Structure prediction input (command "prediction_job", Structure Prediction
dialog): protein sequences from a FASTA file, with modified residues,
glycans, ligands and disulfide bonds, written as the input file of a
structure prediction method.

The methods' formats, as their documentation defines them:

- AlphaFold 3 (docs/input.md of google-deepmind/alphafold3): one JSON
  object, "sequences" of {"protein": {"id", "sequence", "modifications":
  [{"ptmType": CCD code, "ptmPosition": 1-based}]}} and {"ligand": {"id",
  "ccdCodes": [...]}}, "bondedAtomPairs" of [[chain, residue, atom], ...],
  "dialect": "alphafold3", "version" (1 to 4, 4 being the current one).
  A glycan is one ligand entity holding its sugars, bonded to the protein
  and to each other; a disulfide is an SG-SG pair. (docs/input.md says
  bonds "between or within polymer entities" are not supported, but the
  reader itself, folding_input.Input.from_json, only checks that the
  chains exist, the residues are in range, the atoms are not a SMILES
  ligand's and the bonds are unique, so SG-SG pairs go through.)
- RoseTTAFold3 (models/rf3/README.md of RosettaCommons/foundry): a JSON
  list of {"name", "components": [{"seq", "chain_id"} | {"ccd_code"} |
  {"smiles"}], "bonds": ["chain/residue name/number/atom", ...]}. Modified
  residues are written in the sequence itself as "(CCD)"; bonds are free,
  so disulfides can be given.
- ESMFold2 (Biohub): its documented interface is a Python SDK, whose
  StructurePredictionInput holds ProteinInput / LigandInput chains and
  Modification(position, ccd). The file written here mirrors those
  objects; check it against the SDK version in use.
'''

import json
import os
import re
import shlex

ONE_LETTER = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN',
    'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
    'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP',
    'Y': 'TYR', 'V': 'VAL',
}
AMINO_ACIDS = set(ONE_LETTER) | set('UOBZX')

METHODS = (('AlphaFold 3', 'af3'), ('RoseTTAFold3', 'rf3'),
           ('ESMFold2', 'esmfold2'))
ESMFOLD2_MODEL = 'esmfold2-fast-2026-05'
AF3_VERSION = 4

# modified residues: label, CCD code, the residues it replaces
PTMS = [
    ('Phosphoserine', 'SEP', 'S'),
    ('Phosphothreonine', 'TPO', 'T'),
    ('Phosphotyrosine', 'PTR', 'Y'),
    ('N6-acetyllysine', 'ALY', 'K'),
    ('N6-methyllysine', 'MLZ', 'K'),
    ('N6,N6-dimethyllysine', 'MLY', 'K'),
    ('N6,N6,N6-trimethyllysine', 'M3L', 'K'),
    ('4-hydroxyproline', 'HYP', 'P'),
    # deamidation: the residue becomes the acid, so its CCD code is the
    # product's (isoaspartate, IAS, is left out: it links through its
    # side chain, C-gamma, not as an ordinary peptide residue)
    ('Deamidation (Asn to Asp)', 'ASP', 'N'),
    ('Deamidation (Gln to Glu)', 'GLU', 'Q'),
    ('Selenomethionine', 'MSE', 'M'),
    ('S-hydroxycysteine', 'CSO', 'C'),
    ('Citrulline', 'CIR', 'R'),
]

# glycans: each sugar is (CCD code, the sugar it hangs off (0-based, -1 =
# the protein), the parent's oxygen). The first sugar's bond to the
# protein uses the residue's own atom (GLYCOSYLATION_ATOMS).
GLYCANS = [
    ('G0', 'n', [('NAG', -1, ''), ('NAG', 0, 'O4'), ('BMA', 1, 'O4'),
                 ('MAN', 2, 'O3'), ('MAN', 2, 'O6'), ('NAG', 3, 'O2'),
                 ('NAG', 4, 'O2')]),
    ('G0F', 'n', [('NAG', -1, ''), ('NAG', 0, 'O4'), ('BMA', 1, 'O4'),
                  ('MAN', 2, 'O3'), ('MAN', 2, 'O6'), ('NAG', 3, 'O2'),
                  ('NAG', 4, 'O2'), ('FUC', 0, 'O6')]),
    ('G1', 'n', [('NAG', -1, ''), ('NAG', 0, 'O4'), ('BMA', 1, 'O4'),
                 ('MAN', 2, 'O3'), ('MAN', 2, 'O6'), ('NAG', 3, 'O2'),
                 ('NAG', 4, 'O2'), ('GAL', 5, 'O4')]),
    ('G1F', 'n', [('NAG', -1, ''), ('NAG', 0, 'O4'), ('BMA', 1, 'O4'),
                  ('MAN', 2, 'O3'), ('MAN', 2, 'O6'), ('NAG', 3, 'O2'),
                  ('NAG', 4, 'O2'), ('GAL', 5, 'O4'), ('FUC', 0, 'O6')]),
    ('G2', 'n', [('NAG', -1, ''), ('NAG', 0, 'O4'), ('BMA', 1, 'O4'),
                 ('MAN', 2, 'O3'), ('MAN', 2, 'O6'), ('NAG', 3, 'O2'),
                 ('NAG', 4, 'O2'), ('GAL', 5, 'O4'), ('GAL', 6, 'O4')]),
    ('G2F', 'n', [('NAG', -1, ''), ('NAG', 0, 'O4'), ('BMA', 1, 'O4'),
                  ('MAN', 2, 'O3'), ('MAN', 2, 'O6'), ('NAG', 3, 'O2'),
                  ('NAG', 4, 'O2'), ('GAL', 5, 'O4'), ('GAL', 6, 'O4'),
                  ('FUC', 0, 'O6')]),
    ('Man5', 'n', [('NAG', -1, ''), ('NAG', 0, 'O4'), ('BMA', 1, 'O4'),
                   ('MAN', 2, 'O3'), ('MAN', 2, 'O6'), ('MAN', 3, 'O2'),
                   ('MAN', 4, 'O3'), ('MAN', 4, 'O6')]),
    ('Core (Man3)', 'n', [('NAG', -1, ''), ('NAG', 0, 'O4'),
                          ('BMA', 1, 'O4'), ('MAN', 2, 'O3'),
                          ('MAN', 2, 'O6')]),
    ('NAG-NAG', 'n', [('NAG', -1, ''), ('NAG', 0, 'O4')]),
    ('NAG', 'n', [('NAG', -1, '')]),
    ('Tn antigen (GalNAc)', 'o', [('NGA', -1, '')]),
    ('T antigen (core 1)', 'o', [('NGA', -1, ''), ('GAL', 0, 'O3')]),
    ('Core 2', 'o', [('NGA', -1, ''), ('GAL', 0, 'O3'), ('NAG', 0, 'O6')]),
    ('Sialyl-T', 'o', [('NGA', -1, ''), ('GAL', 0, 'O3'),
                       ('SIA', 1, 'O3')]),
]
# where a glycan joins its residue
GLYCOSYLATION_ATOMS = {'N': 'ND2', 'S': 'OG', 'T': 'OG1'}
N_SEQUON = re.compile(r'N[^P][ST]')


def glycan(name):
    for label, kind, sugars in GLYCANS:
        if label == name:
            return kind, sugars
    raise ValueError('unknown glycan "%s"' % name)


def glycan_names(kind=''):
    return [label for label, this, _ in GLYCANS if not kind or this == kind]


def ptm(code):
    for label, this, residues in PTMS:
        if this == code:
            return label, residues
    raise ValueError('unknown modification "%s"' % code)


def ptms_for(letter):
    '''
    The modifications that replace a residue, as [(label, CCD code)]
    '''
    return [(label, code) for label, code, residues in PTMS
            if letter in residues]


def read_fasta(text):
    '''
    Sequences of a FASTA file.

    :return: [(title, sequence)], the sequence upper case without gaps
    '''
    found = []
    title, letters = '', []
    for line in str(text).splitlines():
        line = line.strip()
        if line.startswith('>'):
            if letters:
                found.append((title, ''.join(letters)))
            title, letters = line[1:].strip(), []
        elif line and not line.startswith(';'):
            letters.append(re.sub(r'[\s\-*.]', '', line).upper())
    if letters:
        found.append((title, ''.join(letters)))
    if not found:
        raise ValueError('no sequences in the FASTA text')
    for title, sequence in found:
        unknown = set(sequence) - AMINO_ACIDS
        if unknown:
            raise ValueError('sequence "%s": not amino acids: %s' % (
                title or '?', ''.join(sorted(unknown))))
        if not sequence:
            raise ValueError('sequence "%s" is empty' % (title or '?'))
    return found


def more_ids(taken=()):
    '''
    Chain IDs as they are needed: A, B, ... then AA, AB, ..., leaving out
    the ones taken
    '''
    import string
    used = set(taken)
    letters = string.ascii_uppercase
    for one in list(letters) + [a + b for a in letters for b in letters]:
        if one not in used:
            yield one
    raise ValueError('too many chains')


def chain_ids(count, taken=()):
    '''
    The next "count" chain IDs
    '''
    spare = more_ids(taken)
    return [next(spare) for _ in range(count)]


class Chain:
    '''
    One protein chain: its sequence, and what is on it

    modifications: {residue number (1-based): CCD code}
    glycans: {residue number: glycan name}
    '''

    def __init__(self, chain_id, sequence, title=''):
        self.id = chain_id
        self.sequence = sequence
        self.title = title
        self.copies = 1
        self.modifications = {}
        self.glycans = {}

    def letter(self, number):
        return self.sequence[number - 1]

    def residue(self, number):
        '''
        "S12" for the sequence's own residue, the CCD code if modified
        '''
        code = self.modifications.get(number)
        return code if code else self.letter(number)

    def sequons(self):
        '''
        N-X-S/T sites (X not proline), 1-based
        '''
        return [m.start() + 1 for m in N_SEQUON.finditer(self.sequence)]

    def check(self):
        for number, code in sorted(self.modifications.items()):
            label, residues = ptm(code)
            if self.letter(number) not in residues:
                raise ValueError(
                    '%s %d is %s: %s replaces %s' % (
                        self.id, number, self.letter(number), label,
                        ' or '.join(residues)))
        for number, name in sorted(self.glycans.items()):
            kind, _ = glycan(name)
            letter = self.letter(number)
            if letter not in GLYCOSYLATION_ATOMS:
                raise ValueError('%s %d is %s: glycans go on N (N-linked) '
                                 'or S/T (O-linked)' % (self.id, number,
                                                        letter))
            if (kind == 'n') != (letter == 'N'):
                raise ValueError('%s %d: %s is %s-linked, %s is not' % (
                    self.id, number, name, kind.upper(), letter))
            if number in self.modifications:
                raise ValueError('%s %d has both a modification and a '
                                 'glycan' % (self.id, number))


class Job:
    '''
    What to predict: chains, ligands, disulfide bonds and the method
    '''

    def __init__(self, name='prediction'):
        self.name = name
        self.chains = []
        self.ligands = []     # [(CCD code, copies)]
        self.disulfides = []  # [((chain, number), (chain, number))]
        self.seeds = [1]

    def chain(self, chain_id):
        for one in self.chains:
            if one.id == chain_id:
                return one
        raise ValueError('no chain "%s"' % chain_id)

    def check(self):
        if not self.chains:
            raise ValueError('no sequences: open a FASTA file')
        for one in self.chains:
            one.check()
        for a, b in self.disulfides:
            if a == b:
                raise ValueError('a disulfide bond needs two cysteines')
            for chain_id, number in (a, b):
                chain = self.chain(chain_id)
                if chain.letter(number) != 'C':
                    raise ValueError('%s %d is %s, not a cysteine' % (
                        chain_id, number, chain.letter(number)))
        for code, copies in self.ligands:
            if not re.fullmatch(r'[A-Za-z0-9]{1,5}', code):
                raise ValueError('ligand "%s": give a CCD code, e.g. ATP'
                                 % code)
            if copies < 1:
                raise ValueError('ligand "%s": at least one copy' % code)

    def warnings(self, method):
        found = []
        if method == 'esmfold2':
            found.append(
                'ESMFold2\'s documented input is its Python SDK, not a file: '
                'this JSON mirrors the SDK\'s objects, check it against your '
                'SDK version')
            if self.glycan_count():
                found.append('ESMFold2: glycans are written as ligands '
                             'without their bonds, which the SDK does not '
                             'document')
        return found

    def glycan_count(self):
        return sum(len(c.glycans) for c in self.chains)

    def summary(self):
        return '%d chain(s), %d modification(s), %d glycan(s), %d ' \
            'ligand(s), %d disulfide(s)' % (
                sum(c.copies for c in self.chains),
                sum(len(c.modifications) for c in self.chains),
                self.glycan_count(), sum(n for _, n in self.ligands),
                len(self.disulfides))


def _sugars(chain, number, name):
    '''
    A glycan's sugars and the bonds inside it

    :return: [CCD codes], [(from sugar, atom, to sugar)] 1-based, and the
        protein's own atom for the first bond
    '''
    kind, tree = glycan(name)
    codes = [code for code, _, _ in tree]
    bonds = [(parent + 1, atom, k + 1)
             for k, (_, parent, atom) in enumerate(tree) if parent >= 0]
    return codes, bonds, GLYCOSYLATION_ATOMS[chain.letter(number)]


def _copies(job):
    '''
    Every chain of the job, its copies included: [(chain, chain id)]
    '''
    taken = [c.id for c in job.chains]
    extra = chain_ids(sum(c.copies - 1 for c in job.chains), taken)
    found, spare = [], list(extra)
    for chain in job.chains:
        found.append((chain, chain.id))
        for _ in range(chain.copies - 1):
            found.append((chain, spare.pop(0)))
    return found


def disulfide_pairs(job):
    '''
    The disulfide bonds with the chains' copies filled in: a bond inside
    one chain is made in each of its copies, a bond between two chains
    only between the first of each (which copy pairs with which is not
    something the job says)

    :return: [((chain id, number), (chain id, number))]
    '''
    ids = {}
    for chain, chain_id in _copies(job):
        ids.setdefault(chain.id, []).append(chain_id)
    found = []
    for (chain_a, number_a), (chain_b, number_b) in job.disulfides:
        if chain_a == chain_b:
            found += [((one, number_a), (one, number_b))
                      for one in ids.get(chain_a, [chain_a])]
        else:
            found.append(((ids.get(chain_a, [chain_a])[0], number_a),
                          (ids.get(chain_b, [chain_b])[0], number_b)))
    return found


def alphafold3(job):
    '''
    :return: the AlphaFold 3 input as a dict
    '''
    sequences, bonds = [], []
    spare = more_ids([c.id for c in job.chains] +
                     [i for _, i in _copies(job)])
    for chain, chain_id in _copies(job):
        protein = {'id': chain_id, 'sequence': chain.sequence}
        if chain.modifications:
            protein['modifications'] = [
                {'ptmType': code, 'ptmPosition': number}
                for number, code in sorted(chain.modifications.items())]
        sequences.append({'protein': protein})
        for number, name in sorted(chain.glycans.items()):
            codes, inside, atom = _sugars(chain, number, name)
            glycan_id = next(spare)
            sequences.append({'ligand': {'id': glycan_id,
                                         'ccdCodes': codes}})
            bonds.append([[chain_id, number, atom], [glycan_id, 1, 'C1']])
            for parent, parent_atom, child in inside:
                bonds.append([[glycan_id, parent, parent_atom],
                              [glycan_id, child, 'C1']])
    for code, copies in job.ligands:
        ids = [next(spare) for _ in range(copies)]
        sequences.append({'ligand': {'id': ids if len(ids) > 1 else ids[0],
                                     'ccdCodes': [code.upper()]}})
    for (chain_a, number_a), (chain_b, number_b) in disulfide_pairs(job):
        bonds.append([[chain_a, number_a, 'SG'], [chain_b, number_b, 'SG']])
    found = {'name': job.name, 'modelSeeds': list(job.seeds),
             'sequences': sequences}
    if bonds:
        found['bondedAtomPairs'] = bonds
    found['dialect'] = 'alphafold3'
    found['version'] = AF3_VERSION
    return found


def rosettafold3(job):
    '''
    :return: the RoseTTAFold3 input as a list of one job
    '''
    components, bonds = [], []
    spare = more_ids([i for _, i in _copies(job)])
    for chain, chain_id in _copies(job):
        # modified residues are written into the sequence itself
        letters = []
        for number, letter in enumerate(chain.sequence, 1):
            code = chain.modifications.get(number)
            letters.append('(%s)' % code if code else letter)
        components.append({'seq': ''.join(letters), 'chain_id': chain_id})
        for number, name in sorted(chain.glycans.items()):
            codes, inside, atom = _sugars(chain, number, name)
            ids = [next(spare) for _ in codes]
            for code, sugar_id in zip(codes, ids):
                components.append({'ccd_code': code, 'chain_id': sugar_id})
            name_of = ONE_LETTER[chain.letter(number)]
            bonds.append(['%s/%s/%d/%s' % (chain_id, name_of, number, atom),
                          '%s/%s/1/C1' % (ids[0], codes[0])])
            for parent, parent_atom, child in inside:
                bonds.append([
                    '%s/%s/1/%s' % (ids[parent - 1], codes[parent - 1],
                                    parent_atom),
                    '%s/%s/1/C1' % (ids[child - 1], codes[child - 1])])
    for code, copies in job.ligands:
        for _ in range(copies):
            components.append({'ccd_code': code.upper(),
                               'chain_id': next(spare)})
    for (chain_a, number_a), (chain_b, number_b) in disulfide_pairs(job):
        bonds.append(['%s/CYS/%d/SG' % (chain_a, number_a),
                      '%s/CYS/%d/SG' % (chain_b, number_b)])
    found = {'name': job.name, 'components': components}
    if bonds:
        found['bonds'] = bonds
    return [found]


def esmfold2(job):
    '''
    :return: the ESMFold2 input as a dict, shaped after its Python SDK
    '''
    sequences = []
    spare = more_ids([i for _, i in _copies(job)])
    for chain, chain_id in _copies(job):
        protein = {'id': chain_id, 'sequence': chain.sequence}
        if chain.modifications:
            protein['modifications'] = [
                {'position': number, 'ccd': code}
                for number, code in sorted(chain.modifications.items())]
        sequences.append({'protein': protein})
        for number, name in sorted(chain.glycans.items()):
            codes, _, _ = _sugars(chain, number, name)
            sequences.append({'ligand': {'id': next(spare),
                                         'ccd': codes}})
    for code, copies in job.ligands:
        for _ in range(copies):
            sequences.append({'ligand': {'id': next(spare),
                                         'ccd': [code.upper()]}})
    return {'model': ESMFOLD2_MODEL, 'name': job.name,
            'sequences': sequences}


WRITERS = {'af3': alphafold3, 'rf3': rosettafold3, 'esmfold2': esmfold2}

# running a prediction: the job script takes the input file's name
SUBMIT_MODES = (('Single run (this machine)', 'single'),
                ('Slurm (sbatch)', 'slurm'))
INPUT_MOUNT = '/root/af_input_output'    # AlphaFold 3's input/output folder


def shell_variables(text):
    '''
    The plain NAME="value" assignments of a shell script, resolved
    against each other (so AF3_DIR="$LOCAL_DATA_DIR/apps/af3" comes out
    whole)
    '''
    found = {}
    for line in text.splitlines():
        match = re.match(r'\s*([A-Za-z_][A-Za-z_0-9]*)=(.*)$', line)
        if not match:
            continue
        name, value = match.group(1), match.group(2).strip()
        if value[:1] in '"\'' and value[-1:] == value[:1]:
            value = value[1:-1]
        elif ' ' in value or '$(' in value or '`' in value:
            continue                      # (not a plain assignment)
        found[name] = expand_variables(value, found)
    return found


def expand_variables(text, variables):
    '''
    $NAME and ${NAME} replaced by what the script set them to
    '''
    return re.sub(r'\$\{?([A-Za-z_][A-Za-z_0-9]*)\}?',
                  lambda m: variables.get(m.group(1), m.group(0)), text)


def script_input_dir(path, mount=INPUT_MOUNT):
    '''
    The folder a job script mounts as the prediction's input and output
    folder, which is where it looks for the input file ('' if the script
    does not say)
    '''
    try:
        with open(path) as handle:
            text = handle.read()
    except OSError:
        return ''
    variables = shell_variables(text)
    for match in re.finditer(r'--volume[=\s]+"?([^":\s]+):([^":\s]+)"?',
                             text):
        if match.group(2) == mount:
            return os.path.normpath(expand_variables(match.group(1),
                                                     variables))
    return ''


def submit_command(script, input_name, mode='single', pty=True):
    '''
    What to run to predict the structure: the job script with the input
    file's name, through sbatch for Slurm. A single run goes through
    script(1), which gives it a terminal, because the AlphaFold 3 script
    runs "docker run -it", and docker refuses a TTY without one.

    :return: the command as a list
    '''
    script = str(script).strip()
    input_name = str(input_name).strip()
    if not script:
        raise ValueError('give the job script to run')
    if not os.path.exists(script):
        raise ValueError('no such job script: %s' % script)
    if not input_name:
        raise ValueError('give the input file\'s name')
    if os.sep in input_name:
        raise ValueError('the input file\'s name only (the script looks '
                         'for it in its own input folder)')
    if mode not in dict((v, k) for k, v in SUBMIT_MODES):
        raise ValueError('mode must be single or slurm')
    if mode == 'slurm':
        return ['sbatch', script, input_name]
    inner = ' '.join(shlex.quote(one) for one in
                     ['bash', script, input_name])
    return ['script', '-qec', inner, '/dev/null'] if pty else \
        ['bash', script, input_name]


def write(job, method):
    '''
    :return: the input file's text
    '''
    if method not in WRITERS:
        raise ValueError('method must be one of: ' + ', '.join(WRITERS))
    job.check()
    return json.dumps(WRITERS[method](job), indent=2) + '\n'


def save(job, method, path):
    text = write(job, method)
    with open(path, 'w') as handle:
        handle.write(text)
    return path


TITLE_CHAIN = re.compile(r'(?:^|\b)chain[\s_:|-]*([A-Za-z0-9]{1,2})\b',
                         re.IGNORECASE)


def title_chain_id(title):
    '''
    The chain a FASTA title names, as in ">Chain H" or ">H" ('' if it
    names none)
    '''
    title = str(title).strip()
    if re.fullmatch(r'[A-Za-z0-9]{1,2}', title):
        return title.upper()
    found = TITLE_CHAIN.search(title)
    return found.group(1).upper() if found else ''


def job_from_sequences(found, name='prediction'):
    '''
    A Job whose chains are the sequences read (read_fasta). A title that
    names its chain (">Chain H") gives that chain its ID; the rest are
    lettered in order, leaving out the ones already taken.
    '''
    job = Job(name)
    ids = [title_chain_id(title) for title, _ in found]
    for k, one in enumerate(ids):
        if one and one in ids[:k]:
            ids[k] = ''          # (named twice: letter them instead)
    spare = more_ids([one for one in ids if one])
    for one, (title, sequence) in zip(ids, found):
        job.chains.append(Chain(one or next(spare), sequence, title))
    return job


def job_from_text(text, name='prediction'):
    '''
    A Job from FASTA text, as typed or pasted in
    '''
    return job_from_sequences(read_fasta(text), name)


def job_from_fasta(path, name=''):
    '''
    A Job with the FASTA file's sequences as its chains
    '''
    with open(path) as handle:
        found = read_fasta(handle.read())
    return job_from_sequences(
        found, name or os.path.splitext(os.path.basename(path))[0])
