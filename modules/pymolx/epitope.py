'''
Epitope prediction (parity item L-30): an antibody and its antigen folded
together by AlphaFold 3, to see where the antibody binds.

It is pymolx.prediction with the parts that matter for an antibody-antigen
complex:

- each chain has a role (antibody heavy, light, single chain, or antigen),
  guessed from its sequence (pymolx.antibody.chain_kind);
- many seeds by default. The AlphaFold 3 paper (Abramson et al. 2024)
  scored antibody-antigen interfaces from "1,000 model seeds ... each
  with five diffusion samples", against "the typical 5 seeds" elsewhere,
  and its Fig. 5a shows the quality still rising at 1,000. Seeds are the
  lever there; the samples per seed stayed at AlphaFold 3's own 5;
- the number of diffusion samples each seed produces, which is a flag of
  the run rather than part of the input file;
- optionally no MSA and no templates for the antibody chains. This is
  the combination AlphaFold 3's docs/input.md lists as "both unpairedMsa
  and pairedMsa fields are set to an empty string", which it says runs
  "completely MSA-free" (the two fields have to be set or unset
  together, which is why both are written). An antibody's MSA is mostly
  its germline relatives, and leaving it out skips the slowest part of
  the data pipeline. It is off by default, since input.md calls building
  both MSAs "the recommended option" and this changes what the model is
  given.

The input file is AlphaFold 3's (pymolx.prediction.alphafold3) with the
per-chain MSA fields added.
'''

from . import antibody
from . import prediction

# what a chain is in the complex
ROLES = (('Antibody heavy chain', 'heavy'),
         ('Antibody light chain', 'light'),
         ('Antibody, one chain (scFv, VHH)', 'single'),
         ('Antigen', 'antigen'))
ANTIBODY_ROLES = ('heavy', 'light', 'single')
DEFAULT_ROLE = 'antigen'

# a seed's chance of finding the epitope is modest, so predict many
DEFAULT_SEEDS = 10
MAX_SEEDS = 200
# diffusion samples per seed (AlphaFold 3's --num_diffusion_samples,
# whose own default is 5): they share the seed's MSA search, so more of
# them costs little and widens the spread of poses to compare
DEFAULT_SAMPLES = 10
MAX_SAMPLES = 50
SAMPLES_FLAG = '--num_diffusion_samples'

# what AlphaFold 3 takes as "no MSA" / "no templates" for a chain
NO_MSA = {'unpairedMsa': '', 'pairedMsa': '', 'templates': []}

_KINDS = {'heavy': 'heavy', 'light': 'light', 'scfv': 'single',
          'vhh': 'single', 'fc': 'heavy'}


def role_label(role):
    return dict((v, k) for k, v in ROLES).get(role, role)


def detect_role(sequence):
    '''
    The role a sequence looks like: an antibody chain by its variable
    domain (or an Fc by its constant ones), otherwise the antigen
    '''
    return _KINDS.get(antibody.chain_kind(str(sequence).upper()),
                      DEFAULT_ROLE)


class Job(prediction.Job):
    '''
    An antibody and its antigen to fold together

    roles: {chain id: role}; chains without one are antigens
    '''

    def __init__(self, name='epitope'):
        super().__init__(name)
        self.roles = {}
        self.seeds = list(range(1, DEFAULT_SEEDS + 1))
        self.samples = DEFAULT_SAMPLES
        self.antibody_msa = True    # False: no MSA for the antibody chains

    def role(self, chain_id):
        return self.roles.get(chain_id, DEFAULT_ROLE)

    def set_role(self, chain_id, role):
        if role not in dict((v, k) for k, v in ROLES):
            raise ValueError('role: one of %s'
                             % ', '.join(v for _, v in ROLES))
        self.roles[chain_id] = role

    def detect_roles(self):
        '''
        Guess every chain's role from its sequence, then from the company
        it keeps: a lone heavy variable domain looks like a VHH, but
        beside a light chain it is the heavy chain of an Fv
        '''
        kinds = {}
        for chain in self.chains:
            kinds[chain.id] = antibody.chain_kind(chain.sequence.upper())
            self.roles[chain.id] = _KINDS.get(kinds[chain.id], DEFAULT_ROLE)
        if any(k == 'light' for k in kinds.values()):
            for chain_id, kind in kinds.items():
                if kind == 'vhh':
                    self.roles[chain_id] = 'heavy'

    def chains_of(self, *roles):
        return [c for c in self.chains if self.role(c.id) in roles]

    def antibody_chains(self):
        return self.chains_of(*ANTIBODY_ROLES)

    def antigen_chains(self):
        return self.chains_of('antigen')

    def check(self):
        super().check()
        if not self.antibody_chains():
            raise ValueError('no antibody chain: give the heavy and light '
                             'chains, or an scFv or VHH')
        if not self.antigen_chains():
            raise ValueError('no antigen chain: the antibody has nothing '
                             'to bind')
        if not self.seeds:
            raise ValueError('no seeds: at least one prediction')
        if not 1 <= int(self.samples) <= MAX_SAMPLES:
            raise ValueError('samples per seed: 1 to %d' % MAX_SAMPLES)

    def seed_clash(self, script):
        '''
        What to say when a job script sets --num_seeds and the input file
        holds more than one seed: AlphaFold 3 refuses that pair

        :return: the message, or ''
        '''
        theirs = script_num_seeds(script) if script else None
        if theirs is None or len(self.seeds) <= 1:
            return ''
        return ('the job script sets --num_seeds=%d, and AlphaFold 3 then '
                'wants exactly one seed in the input file, not the %d this '
                'one has ("Input must have one rng seed to set multiple '
                'seeds"). Take --num_seeds out of the script and let the '
                'Seeds box decide, or set Seeds to 1 and let the script '
                'expand it' % (theirs, len(self.seeds)))

    def warnings(self, method='af3'):
        found = super().warnings(method)
        if len(self.seeds) < 5:
            found.append('%d seed(s): the AlphaFold 3 paper ranked its '
                         'antibody-antigen predictions across 1,000 seeds '
                         '(5 samples each) and found their quality still '
                         'rising at that point, so a handful says little'
                         % len(self.seeds))
        heavy = self.chains_of('heavy')
        light = self.chains_of('light')
        if (heavy and not light) or (light and not heavy):
            found.append('a heavy chain without a light one (or the other '
                         'way round): an Fv needs both, unless this is a '
                         'VHH or an scFv ("one chain")')
        if not self.antibody_msa:
            found.append('no MSA for the antibody chains: faster, and an '
                         'antibody\'s MSA says little about it, but it '
                         'does change what the model is given')
        return found

    def predictions(self):
        return len(self.seeds) * int(self.samples)

    def summary(self):
        return '%s; %d seed(s) x %d sample(s) = %d prediction(s)' % (
            ', '.join('%s %s' % (role_label(self.role(c.id)).lower(), c.id)
                      for c in self.chains) or 'no chains',
            len(self.seeds), int(self.samples), self.predictions())


def alphafold3(job):
    '''
    :return: the AlphaFold 3 input as a dict, the antibody chains with no
        MSA and no templates if the job says so
    '''
    found = prediction.alphafold3(job)
    if job.antibody_msa:
        return found
    antibody_ids = {chain_id for chain, chain_id in prediction._copies(job)
                    if job.role(chain.id) in ANTIBODY_ROLES}
    for entry in found['sequences']:
        protein = entry.get('protein')
        if protein is not None and protein.get('id') in antibody_ids:
            protein.update(NO_MSA)
    return found


def write(job, method='af3'):
    '''
    The input file's text (AlphaFold 3 only: this is cofolding with it)
    '''
    if method != 'af3':
        raise ValueError('epitope prediction is AlphaFold 3 cofolding')
    job.check()
    import json
    return json.dumps(alphafold3(job), indent=2) + '\n'


def save(job, path):
    text = write(job)
    with open(path, 'w') as handle:
        handle.write(text)
    return path


def run_arguments(job):
    '''
    What to add to the AlphaFold 3 run beyond the input file: the number
    of diffusion samples (not part of the input file)

    :return: [flag, value]
    '''
    return [SAMPLES_FLAG, str(int(job.samples))]


def script_num_seeds(path):
    '''
    The --num_seeds a job script sets, if it does. AlphaFold 3 expands a
    single seed into that many, and refuses an input file that already
    holds more than one ("Input must have one rng seed to set multiple
    seeds"), so the two ways of asking for seeds can't be combined.

    :return: the number, or None (not set, or the script can't be read)
    '''
    import re
    try:
        with open(path) as handle:
            text = handle.read()
    except OSError:
        return None
    match = re.search(r'--num_seeds[=\s]+(\d+)', text)
    return int(match.group(1)) if match else None


def script_forwards_arguments(path):
    '''
    Whether a job script passes extra arguments on to AlphaFold 3 (it
    needs "$@" or "$2"): without that, the samples flag can't reach it

    :return: True, False, or None when the script can't be read
    '''
    try:
        with open(path) as handle:
            text = handle.read()
    except OSError:
        return None
    return '"$@"' in text or '$@' in text or '$2' in text


def job_from_sequences(found, name='epitope'):
    '''
    An epitope job from read FASTA entries, roles guessed
    '''
    job = prediction.job_from_sequences(found, name)
    result = Job(name)
    result.chains = job.chains
    result.detect_roles()
    return result


def job_from_text(text, name='epitope'):
    return job_from_sequences(prediction.read_fasta(text), name)


def job_from_fasta(path, name=''):
    import os
    found = prediction.read_fasta(open(path).read())
    return job_from_sequences(
        found, name or os.path.splitext(os.path.basename(path))[0])
