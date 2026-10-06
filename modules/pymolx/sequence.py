'''
Sequence viewer data (parity item L-07): one row per chain of the
enabled molecular objects, residue cells with their colors, the active
selection, and selecting residues like the OpenGL sequence viewer
(layer3/Seeker.cpp) does. Qt-free; the widget is
pymolx.gui.sequence_viewer.
'''

import bisect

from pymol import cmd, _cmd

# atom flags (layer2/AtomInfo.h)
FLAG_PROTEIN = 0x00000040
FLAG_NUCLEIC = 0x00000080
FLAG_POLYMER = 0x08000000
FLAG_GUIDE = 0x80000000
# polymer residues: terminal nucleotides can be classed "organic" but
# still have the nucleic bit
POLYMER_FLAGS = FLAG_POLYMER | FLAG_PROTEIN | FLAG_NUCLEIC

# mouse_selection_mode -> selection keyword (SelModeKW, layer1/Scene.cpp)
SELE_MODE_KEYWORDS = ['', 'byresi', 'bychain', 'bysegi', 'byobject',
                      'bymol', 'bca.']

# seq_view_format values shown here; others fall back to residue codes
FORMAT_CODES = 0
FORMAT_NAMES = 1
FORMAT_CHAINS = 3

NUMBER_SPACING = 10  # residue numbers above every 10th residue
TICK_SPACING = 5     # small ticks in between

TEMP_SELE = '_seqview'
BASE_SELE = '_seqview_base'


class Residue:
    '''
    One cell of a row.

    :ivar code: one-letter code ("?" if unknown; polymer residues)
    :ivar text: what the cell shows (one-letter code, residue name, ...)
    :ivar col: first column of the cell, in characters
    :ivar color: PyMOL color index (guide atom, else a carbon, else the
        first atom, like the OpenGL viewer)
    :ivar indices: atom indices (1-based, for select_list mode "index")
    '''
    __slots__ = ('model', 'segi', 'chain', 'resi', 'resv', 'resn', 'polymer',
                 'code', 'text', 'col', 'color', 'indices', '_color_rank')

    def __init__(self, model, segi, chain, resi, resv, resn, polymer,
                 code='?'):
        self.model = model
        self.segi = segi
        self.chain = chain
        self.resi = resi
        self.resv = resv
        self.resn = resn
        self.polymer = polymer
        self.code = code
        self.text = ''
        self.col = 0
        self.color = 0
        self.indices = []
        self._color_rank = -1

    @property
    def key(self):
        return (self.model, self.segi, self.chain, self.resi)

    @property
    def width(self):
        return len(self.text)

    @property
    def macro(self):
        '''
        PyMOL macro notation, e.g. "/1aon//A/ALA`10"
        '''
        return '/%s/%s/%s/%s`%s' % (self.model, self.segi, self.chain,
                                    self.resn, self.resi)


class SeqRow:
    '''
    One chain of one object, or a whole object with its chains side by
    side (like the OpenGL viewer), each chain after a "/A/" marker.

    :ivar chain: None for a whole-object row
    :ivar residues: Residue list, in display order
    :ivar parts: (marker text, residues) per chain; marker text is "" for
        a single-chain row
    :ivar markers: (col, text, residues) of the chain markers, shown in
        the number line above blank cells
    :ivar numbers: (col, text) residue numbers shown above the cells
    :ivar ticks: columns with a tick mark (every TICK_SPACING residues)
    :ivar ncols: width in characters
    '''

    def __init__(self, model, segi, chain):
        self.model = model
        self.segi = segi
        self.chain = chain
        self.residues = []
        self.parts = []
        self.markers = []
        self.numbers = []
        self.ticks = []
        self.ncols = 0
        self._cols = []

    @classmethod
    def merge(cls, rows):
        '''
        One row for an object from its chain rows (in order)
        '''
        merged = cls(rows[0].model, None, None)
        segi = None
        for row in rows:
            if row.segi and row.segi not in (row.chain, segi):
                marker = '/%s/%s/' % (row.segi, row.chain)
            else:
                marker = '/%s/' % row.chain
            segi = row.segi
            merged.parts.append((marker, row.residues))
            merged.residues.extend(row.residues)
        return merged

    @property
    def key(self):
        return (self.model, self.segi, self.chain)

    @property
    def label(self):
        '''
        "object/chain", with the segment if there is one besides the
        chain; just "object" for a whole-object row
        '''
        parts = [self.model]
        if self.segi and self.segi != self.chain:
            parts.append(self.segi)
        if self.chain:
            parts.append(self.chain)
        return '/'.join(parts)

    @property
    def has_polymer(self):
        return any(r.polymer for r in self.residues)

    def residue_at(self, col):
        '''
        :return: Residue whose cell covers column col, or None
        '''
        i = bisect.bisect_right(self._cols, col) - 1
        if i >= 0:
            res = self.residues[i]
            if col < res.col + res.width:
                return res
        return None

    def marker_at(self, col):
        '''
        :return: residues of the chain whose marker covers column col, or
            None
        '''
        for start, text, residues in self.markers:
            if start <= col < start + len(text):
                return residues
        return None

    def nearest_index(self, col):
        '''
        :return: index of the residue at column col, or of the nearest
            one (for drags past the ends or over gaps)
        '''
        i = bisect.bisect_right(self._cols, col) - 1
        if i < 0:
            return 0
        if i + 1 < len(self.residues):
            res = self.residues[i]
            if col >= res.col + res.width and (
                    self.residues[i + 1].col - col <
                    col - (res.col + res.width - 1)):
                return i + 1
        return i

    def layout(self, fmt=FORMAT_CODES):
        '''
        Set the cell texts and columns, chain markers, residue numbers and
        ticks. Non-polymer residues (ligands, ions) show their name with a
        blank column on each side. Ligands and ions of a chain come after
        its polymer residues.
        '''
        if not self.parts:
            self.parts = [('', self.residues)]
        self.residues = []
        self.markers = []
        self.numbers = []
        self.ticks = []
        col = 0
        for marker, residues in self.parts:
            residues.sort(key=lambda res: not res.polymer)
            if marker:
                if col:
                    col += 1
                self.markers.append((col, marker, residues))
                col += len(marker) + 1  # blank cells under "/A/ "
            prev_wide = False
            first_col = col
            for res in residues:
                if fmt == FORMAT_NAMES:
                    res.text = res.resn or '?'
                    wide = True
                elif fmt == FORMAT_CHAINS:
                    res.text = res.chain or '-'
                    wide = False
                elif res.polymer:
                    res.text = res.code
                    wide = False
                else:
                    res.text = res.resn or '?'
                    wide = True
                if col > first_col and (wide or prev_wide):
                    col += 1
                res.col = col
                col += res.width
                prev_wide = wide
            self._number(residues)
            self.residues.extend(residues)
        self.ncols = col
        self._cols = [res.col for res in self.residues]

    def _number(self, residues):
        # numbers above every NUMBER_SPACING-th polymer residue (and the
        # first one of the chain), never overlapping
        free = 0
        if self.numbers:
            col, text = self.numbers[-1]
            free = col + len(text) + 1
        first = True
        for res in residues:
            if not res.polymer:
                continue
            if res.resv % TICK_SPACING == 0:
                self.ticks.append(res.col)
            if (first or res.resv % NUMBER_SPACING == 0) and res.col >= free:
                self.numbers.append((res.col, res.resi))
                free = res.col + len(res.resi) + 1
            first = False


_NUCLEIC_CODES = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'U',
                  'I': 'I', 'N': 'N', 'DA': 'A', 'DC': 'C', 'DG': 'G',
                  'DT': 'T', 'DU': 'U', 'DI': 'I', 'DN': 'N'}


def one_letter(resn, oneletter=''):
    '''
    One-letter code of a polymer residue ("?" if unknown)

    :param oneletter: iterate's "oneletter" (amino acids only)
    '''
    return _NUCLEIC_CODES.get(resn) or oneletter or '?'


def get_change_counts(*, _self=cmd):
    '''
    :return: (changed, dirty) counters of the core: "changed" grows when
        residues, colors or objects may have changed, "dirty" when
        selections may have changed
    '''
    return _cmd.get_seq_change_counts(_self._COb)


def _shown_objects(_self):
    '''
    Enabled molecular objects, in panel order, without hidden ones and
    those excluded with "set seq_view, off, object".
    '''
    index = _self.setting._get_index('seq_view')
    names = []
    for name in _self.get_names('public_objects', enabled_only=1):
        if _self.get_type(name) != 'object:molecule':
            continue
        settings = _self.get_object_settings(name) or ()
        if any(s[0] == index and not s[2] for s in settings):
            continue
        names.append(name)
    return names


def get_rows(fmt=None, by_object=False, *, _self=cmd):
    '''
    Rows of the sequence viewer: one per chain (and segment) of each
    enabled molecular object. Polymer residues come first, in atom
    order, then ligands and ions. Water is left out.

    :param fmt: seq_view_format value (default: the current setting)
    :param by_object: one row per object, chains side by side (like the
        OpenGL viewer)
    :return: list of SeqRow
    '''
    if fmt is None:
        fmt = _self.get_setting_int('seq_view_format')
    objects = _shown_objects(_self)
    if not objects:
        return []

    rows = {}
    residues = {}

    def add(model, segi, chain, resi, resv, resn, color, flags, elem, index,
            oneletter):
        key = (model, segi, chain, resi)
        res = residues.get(key)
        if res is None:
            rkey = key[:3]
            row = rows.get(rkey)
            if row is None:
                row = rows[rkey] = SeqRow(model, segi, chain)
            res = residues[key] = Residue(model, segi, chain, resi, resv,
                                          resn, bool(flags & POLYMER_FLAGS),
                                          one_letter(resn, oneletter))
            row.residues.append(res)
        res.indices.append(index)
        # guide atom > carbon > first atom
        rank = 2 if flags & FLAG_GUIDE else 1 if elem == 'C' else 0
        if rank > res._color_rank:
            res.color = color
            res._color_rank = rank

    selection = '(%s) and not solvent' % ' '.join('%' + n for n in objects)
    _self.iterate(selection, 'add(model, segi, chain, resi, resv, resn, '
                  'color, flags, elem, index, oneletter)',
                  space={'add': add})

    order = {name: i for i, name in enumerate(objects)}
    result = sorted(rows.values(), key=lambda row: order[row.model])
    if by_object:
        chains = {}
        for row in result:
            chains.setdefault(row.model, []).append(row)
        result = [SeqRow.merge(rows) for rows in chains.values()]
    for row in result:
        row.layout(fmt)
    return result


def active_selection(create=False, *, _self=cmd):
    '''
    Name of the active selection (the last enabled one), like
    ExecutiveGetActiveSeleName.

    :param create: if there is none, return a new name ("sele", or
        "selNN" with auto_number_selections) instead of None
    '''
    names = _self.get_names('public_selections', enabled_only=1)
    if names:
        return names[-1]
    if not create:
        return None
    if _self.get_setting_int('auto_number_selections'):
        counter = _self.get_setting_int('sel_counter') + 1
        _self.set('sel_counter', counter, quiet=1)
        return 'sel%02d' % counter
    return 'sele'


def get_selected(name=None, *, _self=cmd):
    '''
    Keys (model, segi, chain, resi) of residues with atoms in a selection.

    :param name: selection name (default: the active selection)
    '''
    if name is None:
        name = active_selection(_self=_self)
        if name is None:
            return set()
    elif name not in _self.get_names('selections'):
        return set()
    keys = set()
    _self.iterate(name, 'add((model, segi, chain, resi))',
                  space={'add': keys.add})
    return keys


def residue_selection(residues, name=TEMP_SELE, *, _self=cmd):
    '''
    Select the atoms of residues, by atom index (no selection syntax, so
    any chain, segment or residue names work).

    :return: name
    '''
    by_model = {}
    for res in residues:
        by_model.setdefault(res.model, []).extend(res.indices)
    _self.select(name, 'none', enable=0)
    part = name + '_obj'
    for model, indices in by_model.items():
        _self.select_list(part, model, indices, mode='index')
        _self.select(name, '?%s or ?%s' % (name, part), enable=0)
    _self.delete(part)
    return name


def select_residues(residues, add=True, base=None, *, _self=cmd):
    '''
    Add residues to the active selection (or remove them), expanded by
    mouse_selection_mode, like clicking residues in the OpenGL viewer.
    Wizards which take selections (mutagenesis) get the result.

    :param residues: Residue list
    :param add: add (True) or remove (False)
    :param base: selection to add to or remove from instead of the
        active one (e.g. its state when a drag started)
    :return: name of the selection
    '''
    name = active_selection(create=True, _self=_self)
    if base is None:
        base = name
    mode = SELE_MODE_KEYWORDS[_self.get_setting_int('mouse_selection_mode')
                              % len(SELE_MODE_KEYWORDS)]
    residue_selection(residues, _self=_self)
    if add:
        expr = '(%s(?%s)) or %s(?%s)' % (mode, base, mode, TEMP_SELE)
    else:
        expr = '(%s(?%s)) and not %s(?%s)' % (mode, base, mode, TEMP_SELE)
    enable = 1 if _self.get_setting_int('auto_show_selections') else -1
    _self.select(name, expr, enable=enable)
    _self.delete(TEMP_SELE)

    wizard = _self.get_wizard()
    if wizard is not None and hasattr(wizard, 'do_select'):
        wizard.do_select(name)
    return name


def save_base(*, _self=cmd):
    '''
    Copy the active selection to BASE_SELE (empty if there is none), the
    starting point of a drag.
    '''
    name = active_selection(_self=_self)
    _self.select(BASE_SELE, '?%s' % name if name else 'none', enable=0)
    return BASE_SELE


def clear_base(*, _self=cmd):
    _self.delete(BASE_SELE)


def deselect_all(*, _self=cmd):
    '''
    Empty the active selection (double click on empty space).
    '''
    name = active_selection(_self=_self)
    if name is not None:
        _self.select(name, 'none', enable=1)
    return name
