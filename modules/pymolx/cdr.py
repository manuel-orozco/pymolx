'''
CDR loops of an antibody (parity item L-31): number a sequence or a
structure's chain with ANARCI and mark the framework regions and the
three complementarity-determining regions under the scheme asked for
(IMGT, Kabat, Chothia, Martin, AHo).

ANARCI gives the numbering; which positions are a CDR is the scheme's own
convention, kept in REGIONS below. A range is matched on the residue
number alone, so insertions (Kabat 100A, 100B ...) fall in the loop they
belong to.

ANARCI needs HMMER. The PyPI "hmmer" wheel puts a Python shim on the path
whose name reaches hmmscan's output header, where Biopython's parser
expects "hmmscan" and otherwise returns no alignment coordinates, so the
real binaries inside the wheel are used instead (hmmer_path).
'''
import os
import re

from pymol import cmd

from . import seqcompare

# numbering schemes ANARCI can apply, in menu order
SCHEMES = (('IMGT', 'imgt'), ('Kabat', 'kabat'), ('Chothia', 'chothia'),
           ('Martin (enhanced Chothia)', 'martin'), ('AHo', 'aho'))
DEFAULT_SCHEME = 'imgt'

# where each scheme puts the loops: {scheme: {'H' or 'L': ((first, last),
# ...)}} for CDR1, CDR2 and CDR3, by residue number, ends included.
# 'L' covers both kappa and lambda chains.
REGIONS = {
    'imgt': {'H': ((27, 38), (56, 65), (105, 117)),
             'L': ((27, 38), (56, 65), (105, 117))},
    'kabat': {'H': ((31, 35), (50, 65), (95, 102)),
              'L': ((24, 34), (50, 56), (89, 97))},
    'chothia': {'H': ((26, 32), (52, 56), (95, 102)),
                'L': ((24, 34), (50, 56), (89, 97))},
    'martin': {'H': ((26, 32), (52, 56), (95, 102)),
               'L': ((24, 34), (50, 56), (89, 97))},
    'aho': {'H': ((27, 42), (57, 76), (109, 138)),
            'L': ((27, 42), (57, 76), (109, 138))},
}
CDR_NAMES = ('CDR1', 'CDR2', 'CDR3')
# ANARCI's chain types: H, and the light ones
LIGHT_TYPES = ('K', 'L')
# ... and the rest, which are T-cell receptor chains
TCR_TYPES = ('A', 'B', 'G', 'D')

_hmmer_path = None


def hmmer_path():
    '''
    The folder with the real HMMER binaries: inside the "hmmer" wheel if
    it is installed (its console scripts are shims whose name confuses
    the parser ANARCI uses), else '' for whatever is on the path

    :return: a folder, or ''
    '''
    global _hmmer_path
    if _hmmer_path is None:
        _hmmer_path = ''
        try:
            import hmmer
            folder = os.path.join(os.path.dirname(hmmer.__file__), 'data',
                                  'bin')
            if os.path.exists(os.path.join(folder, 'hmmscan')):
                _hmmer_path = folder
        except ImportError:
            pass
    return _hmmer_path


def available():
    '''
    :return: '' if ANARCI can be used, else why not
    '''
    try:
        import anarci  # noqa: F401
    except ImportError:
        return ('ANARCI is not installed: "pip install anarci hmmer" in '
                'the environment PyMOL runs in')
    if not hmmer_path():
        from shutil import which
        if not which('hmmscan'):
            return ('HMMER is not installed (ANARCI needs hmmscan): '
                    '"pip install hmmer"')
    return ''


def scheme_label(scheme):
    return dict((v, k) for k, v in SCHEMES).get(scheme, scheme)


class Domain:
    '''
    One numbered variable domain of a sequence

    :ivar positions: [((number, insertion code), one-letter)] of the
        residues the numbering covers, in order
    :ivar chain_type: ANARCI's: H, K, L (or a T-cell receptor's)
    :ivar start, end: where it sits in the sequence given (0-based, end
        included)
    '''

    def __init__(self, positions, chain_type, species, start, end,
                 scheme=DEFAULT_SCHEME):
        self.positions = positions
        self.chain_type = chain_type
        self.species = species
        self.start = start
        self.end = end
        self.scheme = scheme

    @property
    def kind(self):
        '''
        'H' for a heavy chain, 'L' for a light one (the key of REGIONS)
        '''
        return 'L' if self.chain_type in LIGHT_TYPES else 'H'

    @property
    def is_antibody(self):
        return self.chain_type not in TCR_TYPES

    def sequence(self):
        return ''.join(letter for _, letter in self.positions)

    def regions(self):
        '''
        The domain split into framework and CDR regions

        :return: [(name, [((number, insertion), letter)])] in order, the
            names FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4
        '''
        ranges = REGIONS[self.scheme][self.kind]
        found = [('FR1', []), ('CDR1', []), ('FR2', []), ('CDR2', []),
                 ('FR3', []), ('CDR3', []), ('FR4', [])]
        for position, letter in self.positions:
            number = position[0]
            for k, (first, last) in enumerate(ranges):
                if first <= number <= last:
                    found[2 * k + 1][1].append((position, letter))
                    break
            else:
                # the framework before the next loop, or after the last
                after = sum(1 for first, _ in ranges if number > first)
                found[2 * min(after, 3)][1].append((position, letter))
        return found

    def cdrs(self):
        '''
        :return: {name: sequence} for CDR1, CDR2, CDR3
        '''
        return {name: ''.join(letter for _, letter in residues)
                for name, residues in self.regions()
                if name.startswith('CDR')}

    def label(self, position):
        '''
        A position as text: "27" or "100A"
        '''
        number, insertion = position
        return '%d%s' % (number, insertion.strip())


def _numbering_to_positions(numbering):
    '''
    ANARCI's numbering (with its gaps) as [(position, letter)]
    '''
    return [(position, letter) for position, letter in numbering
            if letter != '-']


def domains(sequence, scheme=DEFAULT_SCHEME, _self=cmd):
    '''
    The variable domains ANARCI finds in a sequence (an scFv has two)

    :return: [Domain]
    :raises ValueError: if ANARCI can't be used, or the scheme is unknown
    '''
    scheme = str(scheme).strip().lower()
    if scheme not in REGIONS:
        raise ValueError('scheme: one of %s'
                         % ', '.join(v for _, v in SCHEMES))
    missing = available()
    if missing:
        raise ValueError(missing)
    import anarci
    sequence = re.sub(r'[^A-Za-z]', '', str(sequence)).upper()
    if not sequence:
        raise ValueError('no sequence')
    _, numbered, details = anarci.run_anarci(
        [('query', sequence)], scheme=scheme, hmmerpath=hmmer_path())[:3]
    found = []
    if not numbered or not numbered[0]:
        return found
    for (numbering, start, end), detail in zip(numbered[0], details[0]):
        found.append(Domain(_numbering_to_positions(numbering),
                            detail.get('chain_type', ''),
                            detail.get('species', ''), start, end, scheme))
    return found


def chain_sequence(obj, chain, state=1, _self=cmd):
    '''
    A chain's one-letter sequence and its residue identifiers

    :return: (sequence, [(chain, resi)])
    '''
    rows = []
    _self.iterate('%%%s and chain "%s" and polymer.protein and name CA'
                  % (obj, chain),
                  'rows.append((resi, resn))', space={'rows': rows})
    letters = []
    residues = []
    for resi, resn in rows:
        letters.append(seqcompare.ONE_LETTER.get(resn, 'X'))
        residues.append((chain, resi))
    return ''.join(letters), residues


def annotate(obj, chain='', scheme=DEFAULT_SCHEME, state=1, _self=cmd):
    '''
    Number the protein chains of a structure and find their CDR loops

    :return: [(chain, Domain, [(chain, resi)] of the domain's residues)]
    '''
    chains = [chain] if str(chain).strip() else \
        (_self.get_chains('%' + obj) or [])
    found = []
    for one in chains:
        sequence, residues = chain_sequence(obj, one, state, _self)
        if len(sequence) < 50:
            continue
        for domain in domains(sequence, scheme, _self):
            found.append((one, domain,
                          residues[domain.start:domain.end + 1]))
    return found


def region_residues(domain, residues):
    '''
    The structure's residues of each region of a numbered domain

    :param residues: [(chain, resi)] of the domain, as annotate() gives
    :return: [(name, [(chain, resi)])]
    '''
    out = []
    at = 0
    for name, positions in domain.regions():
        out.append((name, residues[at:at + len(positions)]))
        at += len(positions)
    return out


# colors of the three loops, and of the framework
CDR_COLORS = ('tv_orange', 'tv_green', 'tv_blue')
FRAMEWORK_COLOR = 'grey70'


def selection_name(obj, chain, name, _self=cmd):
    return _self.get_legal_name('%s_%s_%s' % (obj, chain, name))


def find_cdrs(obj, chain='', scheme=DEFAULT_SCHEME, select=1, color=1,
              state=1, quiet=0, *, _self=cmd):
    '''
DESCRIPTION

    "find_cdrs" numbers the antibody chains of a structure with ANARCI
    and marks their CDR loops under a numbering scheme.

USAGE

    find_cdrs object [, chain [, scheme ]]

ARGUMENTS

    object = str: the structure

    chain = str: one chain {default: every protein chain}

    scheme = imgt, kabat, chothia, martin or aho {default: imgt}

    select = 0/1: make a selection per loop, e.g. "igg_H_CDR3", and one
    of all of them, "<object>_cdrs" {default: 1}

    color = 0/1: color the loops (CDR1 orange, CDR2 green, CDR3 blue)
    and the variable domains' framework grey {default: 1}

EXAMPLE

    find_cdrs 1hzh, scheme=kabat

SEE ALSO

    build_antibody, epitope
    '''
    scheme = str(scheme).strip().lower()
    found = annotate(obj, chain, scheme, int(state), _self)
    select, color = int(select), int(color)
    names = []
    for chain_id, domain, residues in found:
        if not domain.is_antibody:
            continue
        for k, (name, where) in enumerate(region_residues(domain,
                                                          residues)):
            if not where:
                continue
            selection = _residues_selection(where)
            if color:
                which = CDR_NAMES.index(name) if name in CDR_NAMES else None
                _self.color(CDR_COLORS[which] if which is not None
                            else FRAMEWORK_COLOR,
                            '(%%%s) and (%s)' % (obj, selection))
            if select and name in CDR_NAMES:
                full = selection_name(obj, chain_id, name, _self)
                _self.select(full, '(%%%s) and (%s)' % (obj, selection),
                             enable=0)
                names.append(full)
    if select and names:
        _self.select(_self.get_legal_name('%s_cdrs' % obj),
                     ' or '.join('(%s)' % n for n in names), enable=0)
    if not int(quiet):
        print(report(found, scheme))
    return found


def report(found, scheme=DEFAULT_SCHEME):
    '''
    The loops as a table

    :param found: what annotate() returned
    '''
    lines = [' CDR loops (%s numbering)' % scheme_label(scheme)]
    if not found:
        return lines[0] + '\n no antibody variable domain found'
    width = max(len(d.cdrs().get(n, '')) for _, d, _ in found
                for n in CDR_NAMES)
    lines.append(' %-6s %-7s %-6s %-*s %s' % (
        'Chain', 'Type', 'Loop', max(width, 8), 'Sequence', 'Residues'))
    for chain_id, domain, residues in found:
        kind = '%s/%s' % (domain.chain_type, domain.species or '?')
        for name, where in region_residues(domain, residues):
            if name not in CDR_NAMES:
                continue
            sequence = ''.join(
                letter for position, letter in domain.regions()[
                    [n for n, _ in domain.regions()].index(name)][1])
            span = '%s-%s' % (where[0][1], where[-1][1]) if where else '-'
            lines.append(' %-6s %-7s %-6s %-*s %s' % (
                chain_id, kind, name, max(width, 8), sequence, span))
    return '\n'.join(lines)


def _residues_selection(keys):
    '''
    "(chain H and resi 27+28) or ..." for [(chain, resi)]
    '''
    by_chain = {}
    for chain, resi in keys:
        by_chain.setdefault(chain, []).append(str(resi))
    return ' or '.join('(chain "%s" and resi %s)' % (chain, '+'.join(res))
                       for chain, res in sorted(by_chain.items())) or 'none'


def extend(_self):
    _self.extend('find_cdrs', find_cdrs)
    _self.find_cdrs = find_cdrs
