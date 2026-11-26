'''
Residue labels (the toolbar's Labels tool): one label per residue, its
one- or three-letter name and number ("K593", "LYS593"), on the atom
asked for.

PyMOL's own "label" writes one label per atom of the selection, so a
label per residue needs the right atoms selected first and the text
spelled out as an expression.
'''

from pymol import cmd

# name -> example, in the order the dialog offers them
NAMES = {'one': 'K593', 'three': 'LYS593', 'none': '593'}

# what PyMOL's "oneletter" gives a residue it has no letter for, and
# what it gives water
NO_LETTER = ('', 'X', '?')
WATER_LETTER = 'O'
NUCLEOTIDES = {'A', 'C', 'G', 'U', 'T', 'I',
               'DA', 'DC', 'DG', 'DT', 'DU', 'DI'}

# the atoms offered for a label; any other name can be asked for
ATOMS = ('CA', 'CB', 'N', 'C', 'O')
# where the label of a residue without the atom asked for goes: the
# first of these it has (nucleotides have a C1'), else its first atom
# (ligands, ions, water)
FALLBACK = ('CA', "C1'", 'C1*')


def residue_text(resn, resi, oneletter='', names='one'):
    '''
    The label of a residue: "K593" (one), "LYS593" (three) or "593"
    (none). A nucleotide's letter is its base ("A12"). A residue without
    a one-letter code (ligands, ions, water) keeps its name with "one"
    too, and a name ending in a digit is set apart from the number
    ("SO4 502").

    :param oneletter: PyMOL's "oneletter" of one of the residue's atoms
    '''
    if names == 'none':
        return resi
    if names == 'one':
        if oneletter not in NO_LETTER and oneletter != WATER_LETTER:
            return oneletter + resi
        if resn in NUCLEOTIDES:
            return resn[-1] + resi
    return resn + (' ' if resn[-1:].isdigit() else '') + resi


def find(selection='all', names='one', atom='CA', *, _self=cmd):
    '''
    The labels of the residues with an atom in the selection, each on
    its atom "atom" (or see FALLBACK). Atoms without coordinates
    (residues missing from the model) get none.

    :return: {(model, index): text}
    '''
    if names not in NAMES:
        raise ValueError('names must be one of: ' + ', '.join(NAMES))
    wanted = [atom.strip() or FALLBACK[0]]
    wanted += [name for name in FALLBACK if name not in wanted]
    ranks = {name: rank for rank, name in enumerate(wanted)}
    other = len(wanted)
    residues = {}

    def see(model, segi, chain, resi, resn, name, index, oneletter):
        rank = ranks.get(name, other)
        entry = residues.get((model, segi, chain, resi))
        if entry is None:
            residues[(model, segi, chain, resi)] = [
                rank, index, resn, oneletter]
        elif rank < entry[0]:
            entry[:2] = rank, index

    _self.iterate('(byres (%s)) and present' % selection,
                  '_ = see(model, segi, chain, resi, resn, name, index, '
                  'oneletter)', space={'see': see})
    return {(key[0], index): residue_text(resn, key[3], oneletter, names)
            for key, (_, index, resn, oneletter) in residues.items()}


def label_residues(selection='all', names='one', atom='CA', size=0,
                   quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    "label_residues" gives each residue of the selection one label, its
    name and number, in place of any labels its atoms had.

USAGE

    label_residues [ selection [, names [, atom [, size ]]]]

ARGUMENTS

    selection = str: atoms of the residues to label {default: all}

    names = one (K593), three (LYS593) or none (593) {default: one}

    atom = str: name of the atom the label is put on. A residue without
    it has its label on CA, a nucleotide on C1', anything else on its
    first atom {default: CA}

    size = float: font size of all labels (the "label_size" setting);
    0 leaves it as it is {default: 0}

NOTES

    Residues without a one-letter code (ligands, ions, water) are
    labeled with their name: "NAG601".

EXAMPLE

    label_residues chain A and resi 590-600
    label_residues sele, three, CB, size=18

SEE ALSO

    label, hide
    '''
    size, quiet = float(size), int(quiet)
    new = find(selection, names, atom, _self=_self)
    residues = '(byres (%s))' % selection
    if new:
        # one label per residue: what its other atoms carried goes
        _self.alter(residues, 'label = new.get((model, index), "")',
                    space={'new': new})
        _self.show('labels', residues)
        _self.rebuild(residues, 'labels')
    if size > 0:
        _self.set('label_size', size, quiet=1)
    if not quiet:
        print(' Labels: %d residue%s labeled.' % (
            len(new), '' if len(new) == 1 else 's'))
    return len(new)


def clear(selection='all', *, _self=cmd):
    '''
    Take the labels off the residues of the selection.

    :return: number of residues that had one
    '''
    residues = set()
    _self.iterate('(byres (%s))' % selection,
                  '_ = residues.add((model, segi, chain, resi)) '
                  'if label else None', space={'residues': residues})
    if residues:
        _self.label('(byres (%s))' % selection, '')
    return len(residues)


def extend(_self):
    _self.extend('label_residues', label_residues)
    _self.label_residues = label_residues
