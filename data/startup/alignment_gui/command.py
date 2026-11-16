'''
Builds the alignment command which the dialog shows and runs (no Qt).
'''

import inspect
import re

# alignment methods, in menu order (those missing from the build are left out)
METHODS = ['align', 'super', 'cealign', 'usalign', 'fit']

WIKI_URL = 'http://pymolwiki.org/index.php/Align'


def available_methods(_self):
    return [name for name in METHODS if name in _self.keyword]


def method_parameters(method, _self):
    '''
    :return: parameter names of an alignment command, in order
    '''
    func = inspect.unwrap(_self.keyword[method][0])
    return list(inspect.signature(func).parameters)


def supports_outlier_rejection(method, _self):
    return 'cycles' in method_parameters(method, _self)


def _one_group(selection):
    '''
    True if brackets enclose the whole selection, e.g. "(a, b)" but not
    "(a, b) and chain H"
    '''
    if not (selection.startswith('(') and selection.endswith(')')):
        return False
    depth = 0
    for k, letter in enumerate(selection):
        depth += (letter == '(') - (letter == ')')
        if depth == 0 and k < len(selection) - 1:
            return False
    return depth == 0


def quote_selection(selection):
    '''
    Selection as a command argument: commas would split it.
    '''
    selection = selection.strip() or 'all'
    if ',' in selection and not _one_group(selection):
        selection = '(%s)' % selection
    return selection


# a residue: number and optional insertion code; a range: two of them
_RESIDUE = r'\d+[A-Za-z]?'
_RESIDUE_TOKEN = re.compile(r'^(%s)(?:[-:](%s))?$' % (_RESIDUE, _RESIDUE))
_CHAIN_TOKEN = re.compile(r'^[A-Za-z0-9]{1,4}$')
_SEPARATORS = re.compile(r'[\s,+]+')
ALL = ('', 'all', '*')


def parse_residues(text):
    '''
    Residue numbers and ranges as typed ("10-120, 150-160", "42", "10:20")
    as a "resi" selection's argument ("10-120+150-160"); '' for all

    :raises ValueError: for anything else
    '''
    tokens = [t for t in _SEPARATORS.split(str(text).strip()) if t]
    if not tokens or tokens[0].lower() in ALL:
        return ''
    found = []
    for token in tokens:
        match = _RESIDUE_TOKEN.match(token)
        if not match:
            raise ValueError('residues: numbers and ranges, e.g. 10-120, 150 '
                             '(not "%s")' % token)
        found.append(match.group(1) + ('-' + match.group(2)
                                       if match.group(2) else ''))
    return '+'.join(found)


def parse_chains(text):
    '''
    Chain IDs as typed ("A", "A+B", "H, L") as a "chain" selection's
    argument ("H+L"); '' for all

    :raises ValueError: for anything else
    '''
    tokens = [t for t in _SEPARATORS.split(str(text).strip()) if t]
    if not tokens or tokens[0].lower() in ALL:
        return ''
    for token in tokens:
        if not _CHAIN_TOKEN.match(token):
            raise ValueError('chains: IDs such as A or H+L (not "%s")'
                             % token)
    return '+'.join(tokens)


def region_selection(selection, chains='', residues=''):
    '''
    The part of a selection to fit: its chains and residues (all if not
    given). The whole object moves whatever part is fitted.

    :raises ValueError: if chains or residues can't be read
    '''
    selection = selection.strip() or 'all'
    chains = parse_chains(chains)
    residues = parse_residues(residues)
    if not chains and not residues:
        return selection
    parts = ['(%s)' % selection]
    if chains:
        parts.append('chain ' + chains)
    if residues:
        parts.append('resi ' + residues)
    return ' and '.join(parts)


def build_command(method, mobile, target, *, many_to_one=True,
                  mobile_state=-1, target_state=-1, object_name='',
                  outlier_rejection=True, cycles=5, cutoff=2.0,
                  mobile_chains='', mobile_residues='',
                  target_chains='', target_residues='', _self):
    '''
    :param mobile_chains, mobile_residues, target_chains, target_residues:
        the part of each selection to fit on ('' for all of it); the whole
        object still moves

    :return: (command name, positional arguments, [(keyword, value)])
    '''
    params = method_parameters(method, _self)
    mobile = quote_selection(region_selection(mobile, mobile_chains,
                                              mobile_residues))
    target = quote_selection(region_selection(target, target_chains,
                                              target_residues))
    options = []

    if many_to_one:
        name, args = 'extra_fit', [mobile, target]
        options.append(('method', method))
    elif params[:2] == ['mobile', 'target']:
        name, args = method, [mobile, target]
    else:
        # e.g. cealign takes (target, mobile)
        name, args = method, []
        options += [('mobile', mobile), ('target', target)]

    if 'cycles' in params:
        if outlier_rejection:
            options += [('cycles', int(cycles)), ('cutoff', float(cutoff))]
        else:
            options.append(('cycles', 0))

    object_name = object_name.strip()
    if object_name and 'object' in params:
        options.append(('object', object_name))

    if 'mobile_state' in params:
        options += [('mobile_state', int(mobile_state)),
                    ('target_state', int(target_state))]

    return name, args, options


def format_command(name, args, options, multiline=True):
    '''
    Command line text; multi-line with "\\" continuations for display.
    '''
    parts = []
    if args:
        parts.append(', '.join(args))
    parts += ['%s=%s' % (key, value) for key, value in options]
    if not multiline:
        return name + ' ' + ', '.join(parts)
    lines = [name + ' ' + parts[0]] + ['    ' + part for part in parts[1:]]
    return ', \\\n'.join(lines)
