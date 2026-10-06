'''
Builds the alignment command which the dialog shows and runs (no Qt).
'''

import inspect

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


def quote_selection(selection):
    '''
    Selection as a command argument: commas would split it.
    '''
    selection = selection.strip() or 'all'
    if ',' in selection and not selection.startswith('('):
        selection = '(%s)' % selection
    return selection


def build_command(method, mobile, target, *, many_to_one=True,
                  mobile_state=-1, target_state=-1, object_name='',
                  outlier_rejection=True, cycles=5, cutoff=2.0, _self):
    '''
    :return: (command name, positional arguments, [(keyword, value)])
    '''
    params = method_parameters(method, _self)
    mobile = quote_selection(mobile)
    target = quote_selection(target)
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
