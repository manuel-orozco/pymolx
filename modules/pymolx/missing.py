'''
Residues missing from the model (no coordinates), for the sequence viewer
(parity item L-07).

mmCIF files: PyMOL itself adds missing polymer residues (from
_entity_poly_seq) as CA atoms without coordinates. PDB files: PyMOL
ignores REMARK 465 (missing residues), so "load" is wrapped to read it
and add the same kind of atoms (_cmd.add_missing_ca, layer4/Cmd.cpp).
Being atoms, they can be selected and are saved in sessions.
'''

import functools
import inspect
import os

from pymol import cmd, _cmd


def parse_remark465(lines):
    '''
    Missing residues from PDB REMARK 465 records.

    :param lines: iterable of PDB lines
    :return: list of (chain, resi, resv, resn), in file order, without
        repeats (NMR entries list them per model)
    '''
    result = []
    seen = set()
    in_table = False
    for line in lines:
        if not line.startswith('REMARK 465'):
            if in_table:
                break
            continue
        if not in_table:
            in_table = 'RES C SSSEQI' in line
            continue
        resn = line[15:18].strip()
        chain = line[19:20].strip()
        number = line[21:26].strip()
        icode = line[26:27].strip()
        try:
            resv = int(number)
        except ValueError:
            continue
        if not resn:
            continue
        resi = number + icode
        key = (chain, resi)
        if key not in seen:
            seen.add(key)
            result.append((chain, resi, resv, resn))
    return result


def _pdb_lines(filename, zipped):
    if zipped == 'gz':
        import gzip
        return gzip.open(filename, 'rt', errors='replace')
    if zipped == 'bz2':
        import bz2
        return bz2.open(filename, 'rt', errors='replace')
    return open(filename, errors='replace')


def add_missing_residues(name, missing, *, _self=cmd):
    '''
    Add residues missing from the model to an object as CA atoms without
    coordinates (what PyMOL's mmCIF reader does), in residue order. Each
    takes the segment of its chain; residues that already have atoms are
    skipped.

    :param missing: list of (chain, resi, resv, resn)
    :return: number of atoms added
    '''
    segis = {}
    existing = set()

    def note(segi, chain, resi):
        segis.setdefault(chain, segi)
        existing.add((chain, resi))

    _self.iterate('%' + name, 'note(segi, chain, resi)',
                  space={'note': note})
    residues = []
    for chain, resi, resv, resn in missing:
        if (chain, resi) in existing:
            continue
        inscode = resi[len(str(resv)):]
        residues.append((segis.get(chain, ''), chain, resv, inscode, resn))
    if not residues:
        return 0
    with _self.lockcm:
        return _cmd.add_missing_ca(_self._COb, name, residues)


def _record(filename, object_name, fmt, before, _self):
    '''
    After a PDB file was loaded: add its REMARK 465 residues to the
    object(s) it created or loaded into.
    '''
    from pymol.importing import filename_to_format
    noext, _, guessed, zipped = filename_to_format(filename)
    if (fmt or guessed) != 'pdb' or not os.path.isfile(filename):
        return
    with _pdb_lines(filename, zipped) as handle:
        missing = parse_remark465(handle)
    if not missing:
        return
    names = set(_self.get_names('objects'))
    targets = names - before
    if not targets:
        name = (object_name or '').strip() or noext
        targets = {name} & names
    for name in targets:
        if _self.get_type(name) == 'object:molecule':
            add_missing_residues(name, missing, _self=_self)


def _wrap_load(func, _self):
    signature = inspect.signature(func)

    @functools.wraps(func)
    def load(*args, **kwargs):
        try:
            bound = signature.bind_partial(*args, **kwargs).arguments
        except TypeError:
            bound = {}
        filename = str(bound.get('filename', ''))
        before = set(_self.get_names('objects'))
        result = func(*args, **kwargs)
        try:
            _record(filename, bound.get('object', ''),
                    str(bound.get('format', '') or ''), before, _self)
        except Exception as e:  # never break loading
            print(' pymolx: could not read missing residues:', e)
        return result

    load._pymolx_missing = True
    return load


def install(_self=cmd):
    '''
    Wrap "load" (API and command language). Called from pymolx._init
    before undo.install, so a load stays one undo step.
    '''
    upstream = _self.load
    if getattr(upstream, '_pymolx_missing', False):
        return
    wrapped = _wrap_load(upstream, _self)
    _self.load = wrapped
    for entry in _self.keyword.values():
        if entry and entry[0] is upstream:
            entry[0] = wrapped
