'''
Modify Sequences (commands "renumber_residues", "set_chain_id",
"set_segment_id", "swap_chains", "merge_chains" and "combine_chains",
Modify Sequences dialog): the residue numbers, chain IDs and segment
IDs of loaded structures, and chains joined into one chain or into one
object.

PyMOL knows a residue by its segment ID, its chain ID and its number
(with the insertion code, "28A"). Changing these with "alter" has
traps: "resv += 5" drops insertion codes, numbering residues one at a
time merges two of them when a new number meets one not yet changed,
and atoms are left out of order until "sort". So every change here is
first worked out as a map {(segi, chain, resi): (segi, chain, resi)} of
whole residues, checked (two residues must not end up as one: nothing
is changed then, and the message names them), and made in one pass per
object, which is sorted afterwards.
'''

import re
import string

from pymol import cmd

from .interfaces import _resi
from .sequence import POLYMER_FLAGS, one_letter

# given to chains of a combined object whose ID another chain has
CHAIN_IDS = string.ascii_uppercase + string.ascii_lowercase + string.digits

# how merge_chains numbers the residues
NUMBERINGS = {
    'continue': 'each chain goes on where the one before ended',
    'consecutive': 'all numbered anew, without gaps',
    'keep': 'as they are',
}
DEFAULT_NUMBERING = 'continue'

PDB_CHAIN = 1       # characters a PDB file holds
PDB_SEGI = 4

_NUMBER = re.compile(r'(-?\d+)(.*)$')
_BAD_ID = re.compile(r'[\s"\'(),]')
_RANGE = re.compile(r'(-?\d+[A-Za-z]?)(?:[-:](-?\d+[A-Za-z]?))?$')


class Change:
    '''
    What a change did

    :ivar count: number of residues changed
    :ivar chains: [(object, segi, chain)] where those residues are now
    :ivar text: for the user
    '''

    def __init__(self, count=0, chains=(), text=''):
        self.count = count
        self.chains = list(chains)
        self.text = text


def split_resi(resi):
    '''
    "28A" -> (28, "A")
    '''
    match = _NUMBER.match(resi)
    return int(match.group(1)), match.group(2)


def resi_selection(text):
    '''
    "100-150, 200A" as a selection ("resi 100-150+200A"; negative numbers
    escaped); "" for no text
    '''
    text = re.sub(r'(?<=[0-9A-Za-z])\s*(-|:|\bto\b)\s*', '-',
                  str(text).strip())
    parts = []
    for part in re.split(r'[\s,+]+', text):
        if not part:
            continue
        match = _RANGE.match(part)
        if not match:
            raise ValueError('residues: "%s" is neither a number nor a '
                             'range like 100-150' % part)
        parts.append('-'.join(_resi(x) for x in match.groups() if x))
    return 'resi ' + '+'.join(parts) if parts else ''


def check_id(value, what='chain ID'):
    value = '' if value is None else str(value).strip()
    if _BAD_ID.search(value):
        raise ValueError('a %s cannot have spaces, quotes, commas or '
                         'brackets' % what)
    return value


def chain_label(obj, segi, chain):
    label = 'chain %s' % chain if chain else 'the chain without an ID'
    if segi:
        label += ' (segment %s)' % segi
    return '%s of %s' % (label, obj)


def chain_selection(obj, segi, chain):
    return '%%%s and segi "%s" and chain "%s"' % (obj, segi, chain)


def _residue_label(key):
    segi, chain, resi = key
    label = 'residue %s of ' % resi
    label += 'chain %s' % chain if chain else 'the chain without an ID'
    if segi:
        label += ' (segment %s)' % segi
    return label


def _span(numbers):
    '''
    "1 to 214" for residue numbers ("28A")
    '''
    numbers = list(numbers)
    if not numbers:
        return ''
    low = min(numbers, key=lambda r: split_resi(r))
    high = max(numbers, key=lambda r: split_resi(r))
    return low if low == high else '%s to %s' % (low, high)


def _some(parts, limit=3):
    parts = list(parts)
    text = '; '.join(parts[:limit])
    if len(parts) > limit:
        text += '; and %d more' % (len(parts) - limit)
    return text


def residues(selection='all', _self=cmd):
    '''
    The residues of a selection (whole residues), in the order of their
    atoms.

    :return: {object: [(segi, chain, resi)]}
    '''
    found = {}

    def add(model, segi, chain, resi):
        found.setdefault(model, {}).setdefault((segi, chain, resi))

    _self.iterate('byres (%s)' % selection, 'add(model, segi, chain, resi)',
                  space={'add': add})
    return {model: list(keys) for model, keys in found.items()}


def polymer_residues(selection='all', _self=cmd):
    '''
    The residues of a selection that are part of a polymer (protein or
    nucleic acid), rather than ligands, ions or waters.

    :return: {(object, segi, chain, resi)}
    '''
    found = set()
    _self.iterate('byres (%s)' % selection,
                  '_ = flags & _polymer and '
                  'found.add((model, segi, chain, resi))',
                  space={'found': found, '_polymer': POLYMER_FLAGS})
    return found


def _by_chain(keys):
    '''
    {(segi, chain): [key]}, chains and residues in their order
    '''
    chains = {}
    for key in keys:
        chains.setdefault(key[:2], []).append(key)
    return chains


def chain_table(selection='all', _self=cmd):
    '''
    The chains of a selection: one row for each segment ID and chain ID
    of an object, in the order of their atoms.

    :return: [{'object', 'segi', 'chain', 'residues' (of the polymer),
        'others' (ligands, ions, waters), 'first', 'last', 'gaps'
        (breaks in the polymer's numbering), 'sequence', 'selection'}]
    '''
    rows = {}

    def add(model, segi, chain, resi, resv, resn, oneletter, flags):
        row = rows.get((model, segi, chain))
        if row is None:
            row = rows[(model, segi, chain)] = {
                'object': model, 'segi': segi, 'chain': chain,
                'polymer': {}, 'other': {}}
        if flags & POLYMER_FLAGS:
            if resi not in row['polymer']:
                row['polymer'][resi] = (resv, one_letter(resn, oneletter))
        else:
            row['other'].setdefault(resi, (resv, ''))

    _self.iterate(selection,
                  'add(model, segi, chain, resi, resv, resn, oneletter, '
                  'flags)', space={'add': add})
    table = []
    for (model, segi, chain), row in rows.items():
        polymer, other = row.pop('polymer'), row.pop('other')
        # a residue with polymer and other atoms is a polymer residue
        other = {k: v for k, v in other.items() if k not in polymer}
        numbered = polymer or other
        values = [v[0] for v in polymer.values()]
        row.update({
            'residues': len(polymer), 'others': len(other),
            'first': next(iter(numbered)),
            'last': next(reversed(numbered)),
            'gaps': sum(1 for a, b in zip(values, values[1:]) if b - a > 1),
            'sequence': ''.join(v[1] for v in polymer.values()),
            'selection': chain_selection(model, segi, chain)})
        table.append(row)
    return table


def apply(changes, _self=cmd):
    '''
    Make the changes, each object in one pass, after checking that no
    two residues of an object end up with the same segment ID, chain ID
    and number.

    :param changes: {object: {(segi, chain, resi): (segi, chain, resi)}}
    :return: number of residues changed
    :raises ValueError: when residues would collide; nothing is changed
    '''
    changes = {obj: {old: new for old, new in moved.items() if old != new}
               for obj, moved in changes.items()}
    changes = {obj: moved for obj, moved in changes.items() if moved}
    for obj, moved in changes.items():
        taken = {}
        for old in residues('%' + obj, _self).get(obj, ()):
            new = moved.get(old, old)
            if new not in taken:
                taken[new] = old
            elif taken[new] == new or old == new:
                # one of the two stays where it is
                comes = old if taken[new] == new else taken[new]
                raise ValueError(
                    '%s: %s would become %s, which is there already' % (
                        obj, _residue_label(comes), _residue_label(new)))
            else:
                raise ValueError('%s: %s and %s would both become %s' % (
                    obj, _residue_label(taken[new]), _residue_label(old),
                    _residue_label(new)))
    for obj, moved in changes.items():
        _self.alter('%' + obj, 'segi, chain, resi = _moved.get('
                    '(segi, chain, resi), (segi, chain, resi))',
                    space={'_moved': moved})
        _self.sort(obj)
    if changes:
        _self.rebuild()
    return sum(len(moved) for moved in changes.values())


def _done(changes, _self):
    '''
    apply(), and where the changed residues are now
    '''
    count = apply(changes, _self)
    chains = {}
    for obj, moved in changes.items():
        for old, new in moved.items():
            if old != new:
                chains.setdefault((obj,) + new[:2])
    return count, list(chains)


def _nothing(selection):
    return ValueError('no atoms in the selection "%s"' % selection)


def renumber(selection='all', start=None, offset=None, gaps=True, _self=cmd):
    '''
    Renumber the residues of a selection, each chain by itself.

    :param offset: add this to every number (insertion codes kept)
    :param start: the number of each chain's first residue: the others
        follow at their distance (gaps=True: gaps and insertion codes
        kept), or one by one (gaps=False: start, start + 1, ...; no
        gaps, no insertion codes)
    :return: Change
    '''
    if (start is None) == (offset is None):
        raise ValueError('give either the start or the offset')
    found = residues(selection, _self)
    if not found:
        raise _nothing(selection)
    changes, parts = {}, []
    for obj, keys in found.items():
        moved = changes[obj] = {}
        for (segi, chain), members in _by_chain(keys).items():
            if start is not None and not gaps:
                for k, key in enumerate(members):
                    moved[key] = (segi, chain, str(int(start) + k))
            else:
                shift = int(offset) if offset is not None else \
                    int(start) - split_resi(members[0][2])[0]
                for key in members:
                    number, code = split_resi(key[2])
                    moved[key] = (segi, chain, '%d%s' % (number + shift,
                                                         code))
            if any(moved[key] != key for key in members):
                parts.append('%s, %s is now %s' % (
                    chain_label(obj, segi, chain),
                    _span(key[2] for key in members),
                    _span(moved[key][2] for key in members)))
    count, chains = _done(changes, _self)
    text = '%d residues renumbered: %s' % (count, _some(parts)) if count \
        else 'the residues have these numbers already'
    return Change(count, chains, text)


def set_ids(selection, chain=None, segi=None, _self=cmd):
    '''
    Give the residues of a selection another chain ID, segment ID or
    both (None: left as it is).

    :return: Change
    '''
    if chain is None and segi is None:
        raise ValueError('give a chain ID or a segment ID')
    chain = None if chain is None else check_id(chain)
    segi = None if segi is None else check_id(segi, 'segment ID')
    found = residues(selection, _self)
    if not found:
        raise _nothing(selection)
    changes = {}
    for obj, keys in found.items():
        changes[obj] = {key: (key[0] if segi is None else segi,
                              key[1] if chain is None else chain, key[2])
                        for key in keys}
    count, chains = _done(changes, _self)
    what = []
    if chain is not None:
        what.append('chain %s' % chain if chain else 'no chain ID')
    if segi is not None:
        what.append('segment %s' % segi if segi else 'no segment ID')
    text = '%d residues now have %s' % (count, ' and '.join(what)) if count \
        else 'the residues have these IDs already'
    return Change(count, chains, text)


def swap(selection, first, second, _self=cmd):
    '''
    Swap two chain IDs among the residues of a selection.

    :return: Change
    '''
    first, second = check_id(first), check_id(second)
    if first == second:
        raise ValueError('the two chain IDs are the same')
    other = {first: second, second: first}
    changes = {}
    for obj, keys in residues(selection, _self).items():
        changes[obj] = {key: (key[0], other[key[1]], key[2])
                        for key in keys if key[1] in other}
    if not any(changes.values()):
        raise ValueError('no chain %s or %s in the selection "%s"' % (
            first, second, selection))
    count, chains = _done(changes, _self)
    return Change(count, chains, 'chains %s and %s swapped (%d residues)' % (
        first or '(no ID)', second or '(no ID)', count))


def merge(selection, chain=None, numbering=DEFAULT_NUMBERING, gap=0,
          start=1, order=(), _self=cmd):
    '''
    Make one chain of the chains of a selection, in each object by
    itself. The residues all get the first chain's segment ID (atoms
    are ordered by segment first, so two segment IDs in a chain
    scramble its residues). The polymers come first, and the chains'
    ligands, ions and waters are numbered after them: in PDB entries
    a chain's waters follow its last residue, where the next chain is
    to go.

    :param chain: the chain's ID {default: the first chain's}
    :param numbering: NUMBERINGS. "continue": the first chain keeps its
        numbers and each later one is shifted to start "gap" + 1 after
        the one before (its own gaps and insertion codes kept), then
        each chain's ligands and waters likewise, without the gap;
        "consecutive": start, start + 1, ...; "keep": no number changes
    :param order: chain IDs in the order wanted {default: as they come
        in the object}
    :return: Change
    '''
    numbering = str(numbering).strip().lower()
    if numbering not in NUMBERINGS:
        raise ValueError('numbering must be one of: ' + ', '.join(NUMBERINGS))
    chain = None if chain is None else check_id(chain)
    gap, start = int(gap), int(start)
    if gap < 0:
        raise ValueError('the gap cannot be negative')
    order = [check_id(c) for c in order]
    found = residues(selection, _self)
    if not found:
        raise _nothing(selection)
    polymer = polymer_residues(selection, _self)
    changes, parts = {}, []
    for obj, keys in found.items():
        chains = _by_chain(keys)
        if len(chains) < 2:
            continue
        names = list(chains)
        names.sort(key=lambda b: order.index(b[1]) if b[1] in order
                   else len(order))
        segi = names[0][0]
        target = names[0][1] if chain is None else chain
        # the polymer of each chain, then what else each chain has
        blocks = [(kind, [key for key in chains[name]
                          if ((obj,) + key in polymer) == kind])
                  for kind in (True, False) for name in names]
        moved = changes[obj] = {}
        spans = {True: [], False: []}
        last, number = None, start
        for kind, members in blocks:
            if not members:
                continue
            numbers = [split_resi(key[2])[0] for key in members]
            shift = 0
            if numbering == 'continue' and last is not None:
                shift = last + 1 + (gap if kind else 0) - min(numbers)
            last = max(numbers) + shift
            for key in members:
                if numbering == 'consecutive':
                    resi = str(number)
                    number += 1
                else:
                    resi = '%d%s' % ((split_resi(key[2])[0] + shift,) +
                                     split_resi(key[2])[1:])
                moved[key] = (segi, target, resi)
                spans[kind].append(resi)
        text = '%s of %s are now chain %s' % (
            ' + '.join(dict.fromkeys(name[1] or '(no ID)' for name in names)),
            obj,
            target or '(no ID)')
        if spans[True]:
            text += ', ' + _span(spans[True])
        if spans[False]:
            text += ', its ligands and waters ' + _span(spans[False])
        parts.append(text)
    if not changes:
        raise ValueError('merging takes two or more chains of one object')
    count, chains = _done(changes, _self)
    return Change(count, chains, 'chains %s' % _some(parts))


def combine(name, selection, delete=False, _self=cmd):
    '''
    Make a new object of the residues of a selection, which may come
    from several objects. A chain ID that a chain of an earlier object
    has is replaced by a free one, so that chains stay apart.

    :param delete: remove the residues from the objects they came from
        (and those objects, when nothing is left of them)
    :return: Change (.renamed: [(object, old chain, new chain)])
    '''
    name = str(name).strip()
    if not name:
        raise ValueError('give the new object a name')
    name = _self.get_legal_name(name)
    if name in _self.get_names('all'):
        raise ValueError('there is a "%s" already' % name)
    # the selection as it is now: the new object may match it later
    source = _self.get_unused_name('_combine_sele')
    temporary, renamed = [], []
    try:
        if not _self.select(source, 'byres (%s)' % selection, enable=0):
            raise _nothing(selection)
        found = residues(source, _self)
        chains = {obj: list(dict.fromkeys(key[1] for key in keys))
                  for obj, keys in found.items()}
        free = [c for c in CHAIN_IDS
                if not any(c in ids for ids in chains.values())]
        used = set()
        for obj in found:
            copy = _self.get_unused_name('_combine')
            _self.create(copy, '%%%s and %s' % (obj, source), zoom=0)
            temporary.append(copy)
            new = {}
            for chain in chains[obj]:
                if chain in used:
                    if not free:
                        raise ValueError('no chain IDs left for the '
                                         'chains of %s' % obj)
                    new[chain] = free.pop(0)
                    renamed.append((obj, chain, new[chain]))
                used.add(new.get(chain, chain))
            if new:
                _self.alter('%' + copy, 'chain = _new.get(chain, chain)',
                            space={'_new': new})
                _self.sort(copy)
        _self.create(name, ' or '.join('%' + t for t in temporary), zoom=0)
        if delete:
            _self.remove(source)
            for obj in found:
                if not _self.count_atoms('%' + obj):
                    _self.delete(obj)
    finally:
        for copy in temporary:
            _self.delete(copy)
        _self.delete(source)
    count = sum(len(keys) for keys in found.values())
    text = '%d residues of %d %s are now the object "%s"' % (
        count, len(found), 'object' if len(found) == 1 else 'objects', name)
    if renamed:
        text += '; ' + _some('chain %s of %s is chain %s there' % (
            old or '(no ID)', obj, new) for obj, old, new in renamed)
    change = Change(count, [(name, segi, chain) for segi, chain in _by_chain(
        residues('%' + name, _self).get(name, ()))], text)
    change.name = name
    change.renamed = renamed
    return change


def _report(command, change, quiet):
    if not int(quiet):
        print(' %s: %s' % (command, change.text))
        for what, limit, ids in (('chain', PDB_CHAIN, {c[2] for c in
                                                        change.chains}),
                                 ('segment', PDB_SEGI, {c[1] for c in
                                                        change.chains})):
            long = sorted(i for i in ids if len(i) > limit)
            if long:
                print(' %s: a PDB file holds %d character%s of a %s ID '
                      '(%s); mmCIF keeps all' % (
                          command, limit, '' if limit == 1 else 's', what,
                          ', '.join(long)))
    return change.count


def _number(value, what):
    try:
        return int(value)
    except (TypeError, ValueError):
        raise ValueError('%s must be a whole number' % what) from None


def renumber_residues(selection='all', start='', offset='', gaps=1, quiet=0,
                      *, _self=cmd):
    '''
DESCRIPTION

    "renumber_residues" gives the residues of a selection new numbers,
    each chain by itself: shifted by an offset, or so that the chain's
    first residue gets the number "start".

USAGE

    renumber_residues [ selection [, start [, offset [, gaps ]]]]

ARGUMENTS

    selection = str: the residues (whole residues) {default: all}

    start = int: number of the first residue of each chain

    offset = int: number added to every residue number, instead of
    "start"

    gaps = 0/1: with "start": 1 keeps the gaps in the numbering and the
    insertion codes (every residue is shifted by the same amount), 0
    numbers the residues one by one: start, start + 1, ... {default: 1}

EXAMPLE

    renumber_residues chain A, offset=100
    renumber_residues chain H and polymer, start=1, gaps=0
    renumber_residues chain A and resi 100-150, 200

NOTES

    Residues outside the selection keep their numbers. If a residue
    would get the number of another one of its chain, nothing is
    changed and the two are named.

    Ligands and waters with the chain's ID are renumbered with it,
    unless the selection leaves them out ("... and polymer").

SEE ALSO

    set_chain_id, set_segment_id, merge_chains, alter
    '''
    start = None if str(start).strip() == '' else _number(start, 'start')
    offset = None if str(offset).strip() == '' else _number(offset, 'offset')
    change = renumber(selection, start, offset, int(gaps), _self=_self)
    return _report('renumber_residues', change, quiet)


def set_chain_id(selection, chain='', quiet=0, *, _self=cmd):
    '''
DESCRIPTION

    "set_chain_id" gives the residues of a selection another chain ID.

USAGE

    set_chain_id selection, chain

ARGUMENTS

    selection = str: the residues (whole residues)

    chain = str: the new chain ID {default: none, an empty chain ID}

EXAMPLE

    set_chain_id chain A, H
    set_chain_id 1abc and solvent, W

NOTES

    When the chain ID is that of another chain, the residues join that
    chain, as long as none of their numbers is in use there (see
    "merge_chains" to renumber them on the way).

SEE ALSO

    swap_chains, set_segment_id, merge_chains, renumber_residues
    '''
    change = set_ids(selection, chain=chain, _self=_self)
    return _report('set_chain_id', change, quiet)


def set_segment_id(selection, segi='', quiet=0, *, _self=cmd):
    '''
DESCRIPTION

    "set_segment_id" gives the residues of a selection another segment
    ID.

USAGE

    set_segment_id selection, segi

ARGUMENTS

    selection = str: the residues (whole residues)

    segi = str: the new segment ID {default: none, an empty segment ID}

EXAMPLE

    set_segment_id chain A, PROA
    set_segment_id all

NOTES

    PyMOL orders atoms by segment ID before chain ID: a chain with two
    segment IDs has its residues out of order (also in saved files).

SEE ALSO

    set_chain_id, renumber_residues
    '''
    change = set_ids(selection, segi=segi, _self=_self)
    return _report('set_segment_id', change, quiet)


def swap_chains(chain1, chain2, selection='all', quiet=0, *, _self=cmd):
    '''
DESCRIPTION

    "swap_chains" swaps two chain IDs.

USAGE

    swap_chains chain1, chain2 [, selection ]

ARGUMENTS

    chain1, chain2 = str: the two chain IDs

    selection = str: where {default: all}

EXAMPLE

    swap_chains H, L, 1abc

SEE ALSO

    set_chain_id
    '''
    change = swap(selection, chain1, chain2, _self=_self)
    return _report('swap_chains', change, quiet)


def merge_chains(selection, chain='', numbering=DEFAULT_NUMBERING, gap=0,
                 start=1, order='', quiet=0, *, _self=cmd):
    '''
DESCRIPTION

    "merge_chains" makes one chain of the chains of a selection (in
    each object by itself).

USAGE

    merge_chains selection [, chain [, numbering [, gap [, start
        [, order ]]]]]

ARGUMENTS

    selection = str: residues of two or more chains

    chain = str: ID of the merged chain {default: the first chain's}

    numbering = continue: the first chain keeps its numbers, and each
    later chain is shifted to start after the one before / consecutive:
    all residues are numbered anew, start, start + 1, ... / keep: the
    numbers stay (fails if two chains share one) {default: continue}.
    Ligands, ions and waters of the chains are numbered after the last
    residue of the merged chain, unless the numbers are kept.

    gap = int: with "continue", numbers left out between two chains
    {default: 0}

    start = int: with "consecutive", the first number {default: 1}

    order = str: chain IDs in the order wanted, e.g. "L H" {default: as
    they come in the object}

EXAMPLE

    merge_chains chain H+L, A
    merge_chains 1abc and chain A+B+C, numbering=consecutive
    merge_chains chain H+L, order=L H, gap=100

NOTES

    Only the names change: no bond is made between the chains. The
    merged chain gets the first chain's segment ID.

SEE ALSO

    combine_chains, set_chain_id, renumber_residues
    '''
    change = merge(selection, None if str(chain).strip() == '' else chain,
                   numbering, _number(gap, 'gap'), _number(start, 'start'),
                   str(order).replace(',', ' ').split(), _self=_self)
    return _report('merge_chains', change, quiet)


def combine_chains(name, selection, delete=0, quiet=0, *, _self=cmd):
    '''
DESCRIPTION

    "combine_chains" makes one new object of chains (or any residues)
    of one or more objects.

USAGE

    combine_chains name, selection [, delete ]

ARGUMENTS

    name = str: name of the new object

    selection = str: the residues (whole residues)

    delete = 0/1: remove them from the objects they came from
    {default: 0}

EXAMPLE

    combine_chains complex, (fab and chain H+L) or (antigen and chain A)

NOTES

    A chain ID which a chain of an earlier object has already is
    replaced by a free one (reported), so that two chains "A" do not
    become one. Use "merge_chains" on the new object to join chains.

SEE ALSO

    merge_chains, create, extract
    '''
    change = combine(name, selection, int(delete), _self=_self)
    return _report('combine_chains', change, quiet)


COMMANDS = (renumber_residues, set_chain_id, set_segment_id, swap_chains,
            merge_chains, combine_chains)


def extend(_self):
    for func in COMMANDS:
        _self.extend(func.__name__, func)
        setattr(_self, func.__name__, func)
