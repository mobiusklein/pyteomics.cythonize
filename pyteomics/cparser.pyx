# cython: embedsignature=True
# cython: profile=True
#   Copyright 2016 Joshua Klein, Lev Levitsky
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

cimport cython
from cpython.ref cimport PyObject, Py_INCREF
from cpython.dict cimport PyDict_GetItem, PyDict_Next, PyDict_SetItem
from cpython.int cimport PyInt_AsLong, PyInt_Check, PyInt_FromLong
from cpython.float cimport PyFloat_AsDouble
from cpython.list cimport (PyList_GET_ITEM, PyList_GetItem, PyList_SetItem, 
                           PyList_SET_ITEM, PyList_Append, PyList_Insert,
                           PyList_Size)
from cpython.tuple cimport PyTuple_GetItem, PyTuple_GET_ITEM
from cpython.sequence cimport PySequence_GetItem
from cpython.exc cimport PyErr_Occurred
from cpython.object cimport PyObject_CallMethodObjArgs, PyObject_Not, PyObject_IsTrue

import re
from collections import deque
import itertools as it
from pyteomics.auxiliary import PyteomicsError, memoize, BasicComposition

cdef:
    list std_amino_acids, std_labels
    str std_nterm, std_cterm
    object _modX_sequence, _modX_group, _modX_split


std_amino_acids = ['Q','W','E','R','T','Y','I','P','A','S',
                   'D','F','G','H','K','L','C','V','N','M']
"""modX labels for the 20 standard amino acids."""

std_nterm = 'H-'
"""modX label for the unmodified N-terminus."""

std_cterm = '-OH'
"""modX label for the unmodified C-terminus."""

std_labels = std_amino_acids + [std_nterm, std_cterm]
"""modX labels for the standard amino acids and unmodified termini."""

_modX_sequence = re.compile(r'^([^-]+-)?((?:[a-z]*[A-Z])+)(-[^-]+)?$')
_modX_group = re.compile(r'[a-z]*[A-Z]')
_modX_split = re.compile(r'([a-z]*)([A-Z])')


@cython.nonecheck(False)
cdef bint is_term_mod(str label):
    """Check if `label` corresponds to a terminal modification.

    Parameters
    ----------
    label : str

    Returns
    -------
    out : bool
    """
    cdef:
        bint result
    result = label.startswith("-")
    if result:
        return result
    else:
        result = label.endswith("-")
        return result

cdef inline object match_modX(str label):
    """Check if `label` is a valid 'modX' label.

    Parameters
    ----------
    label : str

    Returns
    -------
    out : re.match or None
    """
    return PyObject_CallMethodObjArgs(_modX_split, "match", <PyObject*>label, NULL)

cdef inline int is_modX(str label):
    """Check if `label` is a valid 'modX' label.

    Parameters
    ----------
    label : str

    Returns
    -------
    out : bool
    """
    return bool(match_modX(label))

cpdef int length(object sequence, 
                 int show_unmodified_termini=0, int split=0,
                 int allow_unknown_modifications=0, object labels=std_labels):
    """Calculate the number of amino acid residues in a polypeptide
    written in modX notation.

    Parameters
    ----------
    sequence : str or list or dict
        A string with a polypeptide sequence, a list with a parsed sequence or
        a dict of amino acid composition.
    labels : list, optional
        A list of allowed labels for amino acids and terminal modifications.

    Examples
    --------
    >>> length('PEPTIDE')
    7
    >>> length('H-PEPTIDE-OH')
    7
    """
    if not sequence: return 0

    if isinstance(sequence, str) or isinstance(sequence, list):
        if isinstance(sequence, str):
            parsed_sequence = parse(sequence, show_unmodified_termini, 
                                    split, allow_unknown_modifications, labels)
        else:
            parsed_sequence = sequence
        num_term_groups = 0
        if is_term_mod(parsed_sequence[0]):
            num_term_groups += 1
        if is_term_mod(parsed_sequence[-1]):
            num_term_groups += 1
        return len(parsed_sequence) - num_term_groups
    elif isinstance(sequence, dict):
        return sum([amount for aa, amount in sequence.items()
                    if not is_term_mod(aa)])

    raise PyteomicsError('Unsupported type of sequence.')


cdef inline list interpolate_labels(object labels):
    cdef list _labels
    if isinstance(labels, list):
        _labels = <list>labels
    else:
        _labels = list(labels)
    _labels.extend((std_cterm, std_nterm))
    return _labels


cpdef inline tuple _split_label(str label):
    cdef:
        str mod, X
        tuple temp

    
    temp = <tuple>PyObject_CallMethodObjArgs(match_modX(label), "groups", NULL)
    if <PyObject*>temp == NULL:
        raise PyteomicsError('Cannot split a non-modX label: %s' % label)
    mod = <str>PyTuple_GET_ITEM(<tuple>temp, 0)
    X = <str>PyTuple_GET_ITEM(<tuple>temp, 1)
    if not mod:
        return (X,)
    else:
        return temp


cdef object _modX_split_findall = _modX_split.findall


cdef list _parse_modx_sequence_split(str sequence):
    cdef:
        list parsed_sequence
        tuple parts
        Py_ssize_t i, n
        object g

    parsed_sequence = []
    parts = tuple(_modX_split_findall(sequence))
    n = len(parts)

    for i in range(n):
        g = <object>PyTuple_GET_ITEM(parts, i)

        if g[0]:
            parsed_sequence.append(g)
        else:
            parsed_sequence.append((g[1],))
    return parsed_sequence
    

cpdef list parse(str sequence, bint show_unmodified_termini=0, bint split=0,
                 bint allow_unknown_modifications=0, object labels=std_labels):
    """Parse a sequence string written in modX notation into a list of
    labels or (if `split` argument is :py:const:`True`) into a list of
    tuples representing amino acid residues and their modifications.

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.
    show_unmodified_termini : bool, optional
        If :py:const:`True` then the unmodified N- and C-termini are explicitly
        shown in the returned list. Default value is :py:const:`False`.
    split : bool, optional
        If :py:const:`True` then the result will be a list of tuples with 1 to 4
        elements: terminal modification, modification, residue. Default value is
        :py:const:`False`.
    allow_unknown_modifications : bool, optional
        If :py:const:`True` then do not raise an exception when an unknown
        modification of a known amino acid residue is found in the sequence.
        This also includes terminal groups.
        Default value is :py:const:`False`.

        .. note::
            Since version 2.5, this parameter has effect only if `labels`
            are provided.
    labels : container, optional
        A container of allowed labels for amino acids,
        modifications and terminal modifications.
        If not provided, no checks will be done.
        Separate labels for modifications (such as 'p' or 'ox')
        can be supplied, which means they are applicable to all residues.

        .. warning::
            If `show_unmodified_termini` is set to :py:const:`True`, standard
            terminal groups need to be present in `labels`.

        .. warning::
            Avoid using sequences with only one terminal group, as they are
            ambiguous. If you provide one, `labels` (or :py:const:`std_labels`)
            will be used to resolve the ambiguity.

    Returns
    -------
    out : list
        List of tuples with labels of modifications and amino acid residues.

    Examples
    --------
    >>> parse('PEPTIDE', split=True)
    [('P',), ('E',), ('P',), ('T',), ('I',), ('D',), ('E',)]
    >>> parse('H-PEPTIDE')
    ['P', 'E', 'P', 'T', 'I', 'D', 'E']
    >>> parse('PEPTIDE', show_unmodified_termini=True)
    ['H-', 'P', 'E', 'P', 'T', 'I', 'D', 'E', '-OH']
    >>> parse('TEpSToxM', labels=std_labels + ['pS', 'oxM'])
    ['T', 'E', 'pS', 'T', 'oxM']
    >>> parse('zPEPzTIDzE', True, True, labels=std_labels+['z'])
    [('H-', 'z', 'P'), ('E',), ('P',), ('z', 'T'), ('I',), ('D',), ('z', 'E', '-OH')]
    """

    cdef:
        object match
        tuple temp
        str n, body, c, nterm, cterm, mod, X
        list parsed_sequence
        tuple tgroup
        str sgroup
        list _labels
        set slabels
        Py_ssize_t i
        int is_mod, is_residue, in_labels

    _labels = interpolate_labels(labels)
    
    match = PyObject_CallMethodObjArgs(_modX_sequence, "match", <PyObject*>sequence, NULL)
    if match is None:    
        raise PyteomicsError('Not a valid modX sequence: ' + sequence)
    temp = <tuple>PyObject_CallMethodObjArgs(match, "groups", NULL)
    n = <str>PyTuple_GET_ITEM(temp, 0)
    body = <str>PyTuple_GET_ITEM(temp, 1)
    c = <str>PyTuple_GET_ITEM(temp, 2)
    # slabels = set(_labels)
    # labels help save the day when only one terminal group is given
    if c is None and n is not None:
        # we can try to resolve the ambiguity
        if n != std_nterm and n not in _labels:
            # n is the body then
            c = '-' + body
            body = n[:-1]
            n = None

    # Actual parsing
    if split:
        # parsed_sequence = [g if g[0] else (g[1],) for g in re.findall(
        #     _modX_split, body)]
        parsed_sequence = _parse_modx_sequence_split(body)
    else:
        parsed_sequence = <list>PyObject_CallMethodObjArgs(_modX_group, "findall", <PyObject*>body, NULL)
    nterm, cterm = (n or std_nterm), (c or std_cterm)

    if not allow_unknown_modifications:
        if nterm is not None and nterm not in _labels:
            raise PyteomicsError('Unknown label: {}'.format(nterm))
        if cterm is not None and cterm not in _labels:
            raise PyteomicsError('Unknown label: {}'.format(cterm))

    if split:
        i = 0
        for i in range(PyList_Size(parsed_sequence)):
            tgroup = <tuple>PyList_GET_ITEM(parsed_sequence, i)
            if len(tgroup) == 2:
                mod = <str>PyTuple_GET_ITEM(tgroup, 0)
                X = <str>PyTuple_GET_ITEM(tgroup, 1)
            else:
                mod = ''
                X = <str>PyTuple_GET_ITEM(tgroup, 0)
            if ((not mod) and X not in _labels) or\
                not ((mod+X in _labels) or\
                (X in _labels and (mod in _labels or allow_unknown_modifications))):
                raise PyteomicsError('Unknown label: {}'.format(tgroup))

    else:
        i = 0
        for i in range(PyList_Size(parsed_sequence)):
            sgroup = <str>PyList_GET_ITEM(parsed_sequence, i)

            match = PyObject_CallMethodObjArgs(_modX_split, "match", <PyObject*>sgroup, NULL)
            temp = <tuple>PyObject_CallMethodObjArgs(match, "groups", NULL)
            mod = <str>PyTuple_GET_ITEM(temp, 0)
            X = <str>PyTuple_GET_ITEM(temp, 1)
            is_mod = PyObject_IsTrue(mod)
            in_labels = X not in _labels
            if ((not is_mod) and in_labels) or\
                not ((mod+X in _labels) or\
                (X in _labels and (mod in _labels or allow_unknown_modifications))):
                raise PyteomicsError('Unknown label: {}'.format(sgroup))

    # Append terminal labels
    if show_unmodified_termini or nterm != std_nterm:
        if split:
            temp = <tuple>PyList_GET_ITEM(parsed_sequence, 0)
            temp = (<str>nterm, ) + temp
            Py_INCREF(temp)
            PyList_SET_ITEM(<list>parsed_sequence, 0, tuple(temp))
        else:
            PyList_Insert(parsed_sequence, 0, nterm)
    if show_unmodified_termini or cterm != std_cterm:
        if split:
            temp = <tuple>PyList_GET_ITEM(parsed_sequence, len(parsed_sequence) - 1)
            temp = temp + (cterm,)
            Py_INCREF(temp)
            PyList_SET_ITEM(parsed_sequence, len(parsed_sequence) - 1, temp)
        else:
            PyList_Append(parsed_sequence, cterm)


    return parsed_sequence


cdef str tostring(object parsed_sequence, bint show_unmodified_termini=True):
    """Create a string from a parsed sequence.

    Parameters
    ----------
    parsed_sequence : iterable
        Expected to be in one of the formats returned by
        :py:func:`parse`, i.e. list of labels or list of tuples.
    show_unmodified_termini : bool, optional
        Defines the behavior towards standard terminal groups in the input.
        :py:const:`True` means that they will be preserved if present (default).
        :py:const:`False` means that they will be removed. Standard terminal
        groups will not be added if not shown in `parsed_sequence`,
        regardless of this setting.

    Returns
    -------
    sequence : str
    """
    cdef:
        list labels, group_l
        object group, remove_fn
        bint is_term
        Py_ssize_t i, n
    n = len(parsed_sequence)
    labels = []
    for i in range(n):
        group = parsed_sequence[i]
        if isinstance(group, str):
            is_term = group == std_cterm or group == std_nterm
            if is_term or show_unmodified_termini:
                labels.append(group)
        else: # treat `group` as a tuple
            group_l = list(group)
            remove_fn = group_l.remove
            if not show_unmodified_termini:
                if std_cterm in group_l:
                    remove_fn(std_cterm)
                if std_nterm in group_l:
                    remove_fn(std_nterm)
            labels.append(''.join(group_l))
    return ''.join(labels)


cpdef dict amino_acid_composition(object sequence, bint show_unmodified_termini=0, bint term_aa=0,
                                  bint allow_unknown_modifications=0, object labels=std_labels):
    """Calculate amino acid composition of a polypeptide.

    Parameters
    ----------
    sequence : str or list
        The sequence of a polypeptide or a list with a parsed sequence.
    show_unmodified_termini : bool, optional
        If :py:const:`True` then the unmodified N- and C-terminus are explicitly
        shown in the returned dict. Default value is :py:const:`False`.
    term_aa : bool, optional
        If :py:const:`True` then the terminal amino acid residues are
        artificially modified with `nterm` or `cterm` modification.
        Default value is :py:const:`False`.
    allow_unknown_modifications : bool, optional
        If :py:const:`True` then do not raise an exception when an unknown
        modification of a known amino acid residue is found in the sequence.
        Default value is :py:const:`False`.
    labels : list, optional
        A list of allowed labels for amino acids and terminal modifications.

    Returns
    -------
    out : dict
        A dictionary of amino acid composition.

    Examples
    --------
    >>> amino_acid_composition('PEPTIDE') == \
    {'I': 1, 'P': 2, 'E': 2, 'T': 1, 'D': 1}
    True
    >>> amino_acid_composition('PEPTDE', term_aa=True) == \
    {'ctermE': 1, 'E': 1, 'D': 1, 'P': 1, 'T': 1, 'ntermP': 1}
    True
    >>> amino_acid_composition('PEPpTIDE', labels=std_labels+['pT']) == \
    {'I': 1, 'P': 2, 'E': 2, 'D': 1, 'pT': 1}
    True
    """
    cdef:
        list parsed_sequence
        int nterm_aa_position, cterm_aa_position, aa_count
        dict aa_dict
        tuple aa
        list _labels
        object temp
        PyObject* pvalue
        Py_ssize_t i

    _labels = interpolate_labels(labels)

    if isinstance(sequence, str):
        parsed_sequence = parse(<str>sequence, show_unmodified_termini, split=False,
                                allow_unknown_modifications=allow_unknown_modifications,
                                labels=_labels)
    elif isinstance(sequence, list):
        if sequence and isinstance(<object>PyList_GetItem(sequence, 0), tuple):
            parsed_sequence = parse(tostring(sequence, True),
                show_unmodified_termini,
                split=False,
                allow_unknown_modifications=allow_unknown_modifications,
                labels=_labels)
        else:
            parsed_sequence = sequence
    else:
        raise PyteomicsError('Unsupported type of a sequence.'
                'Must be str or list, not %s' % type(sequence))

    aa_dict = {}
    # Process terminal amino acids.
    if term_aa:
        if is_term_mod(<object>PyList_GetItem(parsed_sequence, 0)):
            nterm_aa_position = 1
        else:
            nterm_aa_position = 0
        if is_term_mod(parsed_sequence[-1]):
            cterm_aa_position = PyList_Size(parsed_sequence) - 2
        else:
            cterm_aa_position = PyList_Size(parsed_sequence) - 1
        if PyList_Size(parsed_sequence) > 1:
            # temp = PyInt_FromLong(cterm_aa_position)
            # temp = PyObject_CallMethodObjArgs(parsed_sequence, "pop", <PyObject*>temp, NULL)
            temp = parsed_sequence.pop(cterm_aa_position)
            PyDict_SetItem(aa_dict, 'cterm' + temp, 1)
        # temp = PyInt_FromLong(nterm_aa_position)
        # temp = PyObject_CallMethodObjArgs(parsed_sequence, "pop", <PyObject*>temp, NULL)
        temp = parsed_sequence.pop(nterm_aa_position)
        PyDict_SetItem(aa_dict, 'nterm' + temp, 1)
    # Process core amino acids.
    i = 0
    for i in range(PyList_Size(parsed_sequence)):
        aa = <tuple>PyList_GET_ITEM(parsed_sequence, i)
        pvalue = PyDict_GetItem(aa_dict, aa)
        if pvalue == NULL:
            aa_count = 0
        else:
            aa_count = PyInt_AsLong(<object>pvalue)
        aa_count += 1
        PyDict_SetItem(aa_dict, aa, aa_count)

    return aa_dict
