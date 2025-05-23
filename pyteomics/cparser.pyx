# cython: embedsignature=True
# cython: profile=False

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

from libc.stdlib cimport malloc, free

cimport cython
from cpython.ref cimport PyObject, Py_INCREF
from cpython.iterator cimport PyIter_Next
from cpython.dict cimport PyDict_GetItem, PyDict_Next, PyDict_SetItem, PyDict_Keys, PyDict_Values
from pyteomics.ccompat cimport PyInt_AsLong, PyInt_Check, PyInt_FromLong
from cpython.float cimport PyFloat_AsDouble
from cpython.list cimport (PyList_GET_ITEM, PyList_GetItem, PyList_SetItem,
                           PyList_SET_ITEM, PyList_Append, PyList_Insert,
                           PyList_Size, PyList_New)
from cpython.tuple cimport (
    PyTuple_GetItem, PyTuple_GET_ITEM, PyTuple_SET_ITEM, PyTuple_New, PyTuple_Size)
from cpython.sequence cimport PySequence_GetItem, PySequence_Fast, PySequence_Fast_GET_ITEM, PySequence_Fast_GET_SIZE
from cpython.exc cimport PyErr_Occurred
from cpython.object cimport PyObject_CallMethodObjArgs, PyObject_Not, PyObject_IsTrue

from pyteomics.ccompat cimport PyStr_AsUTF8AndSize, PyStr_FromStringAndSize

cdef extern from "<ctype.h>" nogil:
    int isupper(int argument)


import re
from collections import deque
import itertools as it
from pyteomics.auxiliary import PyteomicsError, memoize, BasicComposition

cdef:
    list std_amino_acids, std_labels
    str std_nterm, std_cterm
    object _modX_sequence, _modX_group, _modX_split

DEF TOSTRING_STACK_BUFFER_SIZE = 150

std_amino_acids = ['Q','W','E','R','T','Y','I','P','A','S',
                   'D','F','G','H','K','L','C','V','N','M']
"""modX labels for the 20 standard amino acids."""

std_nterm = 'H-'
"""modX label for the unmodified N-terminus."""

std_cterm = '-OH'
"""modX label for the unmodified C-terminus."""

std_labels = std_amino_acids + [std_nterm, std_cterm]
"""modX labels for the standard amino acids and unmodified termini."""

_modX_sequence = re.compile(r'^([^-]+-)?((?:[^A-Z-]*[A-Z])+)(-[^-]+)?$')
_modX_group = re.compile(r'[^A-Z-]*[A-Z]')
_modX_split = re.compile(r'([^A-Z-]*)([A-Z])')


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
    if labels is None:
        _labels = list(std_labels)
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
                 bint allow_unknown_modifications=0, object labels=None):
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

    if labels is None:
        allow_unknown_modifications = True
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
    # if split:
        # parsed_sequence = [g if g[0] else (g[1],) for g in re.findall(
        #     _modX_split, body)]
    #     parsed_sequence = _parse_modx_sequence_split(body)
    # else:
    #     parsed_sequence = <list>PyObject_CallMethodObjArgs(_modX_group, "findall", <PyObject*>body, NULL)
    parsed_sequence = _tokenize_modx_sequence(body, split)

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


cdef list _tokenize_modx_sequence(str sequence, bint split=True):
    cdef:
        char* csequence
        char c
        Py_ssize_t i, n, parsed_n, start, end
        list tokens
        str nterm, cterm

    start = 0
    end = 0
    tokens = []
    i = 0
    parsed_n = 0
    csequence = PyStr_AsUTF8AndSize(sequence, &n)
    nterm = None
    cterm = None

    while i < n:
        c = csequence[i]
        # We've reached an amino acid. Process the current modX block
        if isupper(c):
            end = i + 1
            if not split:
                tokens.append(PyStr_FromStringAndSize(&csequence[start], end - start))
            else:
                if end - start > 1:
                    mod = PyStr_FromStringAndSize(&csequence[start], (end - start) - 1)
                    x = PyStr_FromStringAndSize(&csequence[end - 1], 1)
                    tokens.append((mod, x))
                else:
                    x = PyStr_FromStringAndSize(&csequence[end - 1], 1)
                    tokens.append((x, ))
            parsed_n += 1
            start = i + 1
        # We've reached a terminal modification boundary
        elif c == '-':
            # If we haven't parsed any other tokens
            if parsed_n == 0:
                end = i + 1
                if split:
                    nterm = PyStr_FromStringAndSize(&csequence[start], end - start)
                else:
                    tokens.append(nterm)
                start = i + 1

        i += 1
    end = n
    # The end of the last position was different from the end of the string so we must have
    # a C-terminal
    if start != end:
        cterm = PyStr_FromStringAndSize(&csequence[start], end - start)
        if not split:
            tokens.append(cterm)
        else:
            tokens[-1] = tokens[-1] + (cterm, )
    if split and nterm:
        tokens[0] = (nterm, ) + tokens[0]
    return tokens


cpdef str tostring(object parsed_sequence, bint show_unmodified_termini=True):
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
        object group
        str group_j
        Py_ssize_t i, n, acc, sz
        Py_ssize_t j, m, k
        object sequence_fast
        SequencePosition pos
        str result

        char* cstr
        char[TOSTRING_STACK_BUFFER_SIZE] prealloc_buffer
        char* string_buffer
        bint needs_free
    needs_free = False
    sz = 0
    acc = 0
    sequence_fast = PySequence_Fast(parsed_sequence, "tostring requires an iterable sequence!")
    n = PySequence_Fast_GET_SIZE(sequence_fast)
    for i in range(n):
        group = <object>PySequence_Fast_GET_ITEM(sequence_fast, i)
        if isinstance(group, tuple):
            m = PyTuple_Size(group)
            for j in range(m):
                group_j = <str>PyTuple_GetItem(group, j)
                PyStr_AsUTF8AndSize(group_j, &sz)
                acc += sz
        elif isinstance(group, str):
            PyStr_AsUTF8AndSize(group, &sz)
            acc += sz
        elif isinstance(group, SequencePosition):
            pos = <SequencePosition>group
            if pos.terminal is not None:
                cstr = PyStr_AsUTF8AndSize(pos.terminal, &sz)
                acc += sz
            if pos.modification is not None:
                cstr = PyStr_AsUTF8AndSize(pos.modification, &sz)
                acc += sz
            if pos.amino_acid is not None:
                cstr = PyStr_AsUTF8AndSize(pos.amino_acid, &sz)
                acc += sz
        else:
            raise TypeError(type(group))

    if acc < TOSTRING_STACK_BUFFER_SIZE:
        string_buffer = prealloc_buffer
    else:
        string_buffer = <char*>malloc(sizeof(char*) * acc)
        needs_free = True
    acc = 0
    for i in range(n):
        group = <object>PySequence_Fast_GET_ITEM(sequence_fast, i)
        if isinstance(group, tuple):
            m = PyTuple_Size(group)
            for j in range(m):
                group_j = <str>PyTuple_GetItem(group, j)
                if (group_j == std_cterm or group_j == std_nterm) and not show_unmodified_termini:
                    continue

                cstr = PyStr_AsUTF8AndSize(group_j, &sz)
                for k in range(sz):
                    string_buffer[acc] = cstr[k]
                    acc += 1

        elif isinstance(group, str):
            group_j = <str>group
            if (group_j == std_cterm or group_j == std_nterm) and not show_unmodified_termini:
                continue

            PyStr_AsUTF8AndSize(group, &sz)
            cstr = PyStr_AsUTF8AndSize(group, &sz)
            for k in range(sz):
                string_buffer[acc] = cstr[k]
                acc += 1

        elif isinstance(group, SequencePosition):
            pos = <SequencePosition>group
            if pos.position_type == PositionType.n_term:
                if pos.terminal is not None:
                    if show_unmodified_termini or pos.terminal != std_nterm:
                        cstr = PyStr_AsUTF8AndSize(pos.terminal, &sz)
                        for j in range(sz):
                            string_buffer[acc] = cstr[j]
                            acc += 1

            if pos.modification is not None:
                cstr = PyStr_AsUTF8AndSize(pos.modification, &sz)
                for j in range(sz):
                    string_buffer[acc] = cstr[j]
                    acc += 1

            if pos.amino_acid is not None:
                cstr = PyStr_AsUTF8AndSize(pos.amino_acid, &sz)
                for j in range(sz):
                    string_buffer[acc] = cstr[j]
                    acc += 1

            if pos.position_type == PositionType.c_term:
                if pos.terminal is not None:
                    if show_unmodified_termini or pos.terminal != std_cterm:
                        cstr = PyStr_AsUTF8AndSize(pos.terminal, &sz)
                        for j in range(sz):
                            string_buffer[acc] = cstr[j]
                            acc += 1

    result = PyStr_FromStringAndSize(string_buffer, acc)
    if needs_free:
        free(string_buffer)
    return result



cpdef dict amino_acid_composition(object sequence, bint show_unmodified_termini=0, bint term_aa=0,
                                  bint allow_unknown_modifications=0, object labels=None):
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

    if labels is None:
        allow_unknown_modifications = True

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


ctypedef fused tuple_or_list:
    tuple
    list


cpdef enum PositionType:
    n_term = 1
    c_term = 2
    internal = 3
    universal = 4


cdef inline list copy_list(list source):
    cdef:
        size_t i, n
        list result
        object val
    n = PyList_Size(source)
    result = PyList_New(n)
    for i in range(n):
        val = <object>PyList_GET_ITEM(source, i)
        Py_INCREF(val)
        PyList_SET_ITEM(result, i, val)
    return result


@cython.final
@cython.freelist(5000)
cdef class SequencePosition(object):
    cdef:
        public str modification
        public PositionType position_type
        public str amino_acid
        public str terminal

    def __init__(self, amino_acid, modification=None, terminal=None, position_type=None):
        if position_type is None:
            position_type = PositionType.internal
        self.amino_acid = amino_acid
        self.modification = modification
        self.terminal = terminal
        self.position_type = position_type

    @staticmethod
    cdef SequencePosition _create(str amino_acid, str modification, str terminal, PositionType position_type):
        cdef SequencePosition self = SequencePosition.__new__(SequencePosition)
        self.amino_acid = amino_acid
        self.modification = modification
        self.terminal = terminal
        self.position_type = position_type
        return self

    def __getitem__(self, i):
        if isinstance(i, slice):
            return self.as_tuple()[i]
        return self.get_index(i)

    cpdef str get_index(self, long i):
        cdef:
            long size, n_mods
        size = self.get_size()
        n_mods = self.modification is not None
        if i >= size:
            raise IndexError(i)
        elif i < 0:
            if size + i < 0:
                raise IndexError(i)
            i = size + i
        if self.position_type == PositionType.n_term:
            if i == 0:
                return self.terminal
            elif (i - 1) == n_mods:
                return self.amino_acid
            else:
                return self.modification
        elif self.position_type == PositionType.internal:
            if i == n_mods:
                return self.amino_acid
            else:
                return self.modification
        else:
            if i - n_mods == 0:
                return self.amino_acid
            elif i - n_mods == 1:
                return self.terminal
            else:
                return self.modification

    cdef size_t get_size(self):
        cdef:
            size_t n_mods

        n_mods = self.modification is not None
        if self.position_type == PositionType.internal:
            return n_mods + 1
        else:
            return n_mods + 2

    cdef bint is_modified(self):
        return self.modification is not None

    cpdef add_modification(self, str mod):
        self.modification = mod

    def __repr__(self):
        return str(self.as_tuple())

    def __len__(self):
        return self.get_size()

    cpdef tuple as_tuple(self):
        cdef:
            size_t size, i, n_mods, j
            tuple result
            str mod
        size = self.get_size()
        result = PyTuple_New(size)
        i = 0
        if self.position_type == PositionType.n_term:
            Py_INCREF(self.terminal)
            PyTuple_SET_ITEM(result, i, self.terminal)
            i += 1
        if self.is_modified():
            mod = self.modification
            Py_INCREF(mod)
            PyTuple_SET_ITEM(result, i, mod)
            i += 1
        Py_INCREF(self.amino_acid)
        PyTuple_SET_ITEM(result, i, self.amino_acid)
        i += 1
        if self.position_type == PositionType.c_term:
            Py_INCREF(self.terminal)
            PyTuple_SET_ITEM(result, i, self.terminal)
            i += 1
        return result

    cdef SequencePosition copy(self):
        return SequencePosition._create(
            self.amino_acid, (self.modification), self.terminal, self.position_type)

    @classmethod
    def from_tuple(cls, tuple source):
        cdef:
            size_t size, i
            str entry
            str last_entry
            str terminal
            str amino_acid
            str modification
            PositionType position_type

        position_type = PositionType.internal
        modification = None
        size = PyTuple_Size(source)
        last_entry = None
        terminal = None
        for i in range(size):
            entry = <str>PyTuple_GET_ITEM(source, i)
            if entry.endswith("-"):
                if position_type != PositionType.internal:
                    raise ValueError(
                        "Position already has position type %s, second type provided %s" % (
                            position_type, PositionType.n_term))
                position_type = PositionType.n_term
                terminal = entry
            elif entry.startswith("-"):
                if position_type != PositionType.internal:
                    raise ValueError(
                        "Position already has position type %s, second type provided %s" % (
                            position_type, PositionType.c_term))
                position_type = PositionType.c_term
                terminal = entry
            else:
                if last_entry is not None:
                    modification = last_entry
                last_entry = entry
        amino_acid = last_entry
        return cls(amino_acid, modification, terminal, position_type)

    @classmethod
    def parse(cls, str sequence, **kwargs):
        cdef:
            size_t n, i
            list tokens, result
            object converter
        converter = cls.from_tuple
        tokens = parse(sequence, show_unmodified_termini=True, split=True, **kwargs)
        n = PyList_Size(tokens)
        result = []
        for i in range(n):
            result.append(converter(<tuple>PyList_GET_ITEM(tokens, i)))
        return result

    @classmethod
    def from_iterable(cls, object iterable):
        cdef:
            size_t i, n
            list result
            object converter

        converter = cls.from_tuple
        result = []
        for obj in iterable:
            if isinstance(obj, SequencePosition):
                result.append((<SequencePosition>obj).copy())
            else:
                result.append(converter(<tuple?>obj))
        return result

    cdef bint _eq(self, SequencePosition other):
        if self.amino_acid != other.amino_acid:
            return False
        elif self.modification != other.modification:
            return False
        elif self.terminal != other.terminal:
            return False
        return True

    cdef bint _eq_tuple(self, tuple other):
        return self.as_tuple() == other

    def __eq__(self, other):
        if other is None:
            return False
        if isinstance(other, SequencePosition):
            return self._eq(<SequencePosition>other)
        elif isinstance(other, tuple):
            return self._eq_tuple(<tuple>other)
        else:
            return self._eq_tuple(tuple(other))

    def __ne__(self, other):
        return not (self == other)


cpdef str _sequence_tostring(list sequence, bint show_unmodified_termini=True):
    cdef:
        size_t i, n, acc, j
        Py_ssize_t sz
        char* cstr
        char[TOSTRING_STACK_BUFFER_SIZE] prealloc_buffer
        char* string_buffer
        bint needs_free
        SequencePosition pos
        str result

    needs_free = False
    sz = 0
    n = PyList_Size(sequence)
    acc = 0
    for i in range(n):
        pos = <SequencePosition?>PyList_GET_ITEM(sequence, i)
        if pos.terminal is not None:
            cstr = PyStr_AsUTF8AndSize(pos.terminal, &sz)
            acc += sz
        if pos.modification is not None:
            cstr = PyStr_AsUTF8AndSize(pos.modification, &sz)
            acc += sz
        if pos.amino_acid is not None:
            cstr = PyStr_AsUTF8AndSize(pos.amino_acid, &sz)
            acc += sz
    if acc < TOSTRING_STACK_BUFFER_SIZE:
        string_buffer = prealloc_buffer
    else:
        string_buffer = <char*>malloc(sizeof(char*) * acc)
        needs_free = True
    acc = 0
    for i in range(n):
        pos = <SequencePosition?>PyList_GET_ITEM(sequence, i)
        if pos.position_type == PositionType.n_term:
            if pos.terminal is not None:
                cstr = PyStr_AsUTF8AndSize(pos.terminal, &sz)
                for j in range(sz):
                    string_buffer[acc] = cstr[j]
                    acc += 1

        if pos.modification is not None:
            cstr = PyStr_AsUTF8AndSize(pos.modification, &sz)
            for j in range(sz):
                string_buffer[acc] = cstr[j]
                acc += 1

        if pos.amino_acid is not None:
            cstr = PyStr_AsUTF8AndSize(pos.amino_acid, &sz)
            for j in range(sz):
                string_buffer[acc] = cstr[j]
                acc += 1

        if pos.position_type == PositionType.c_term:
            if pos.terminal is not None:
                cstr = PyStr_AsUTF8AndSize(pos.terminal, &sz)
                for j in range(sz):
                    string_buffer[acc] = cstr[j]
                    acc += 1

    result = PyStr_FromStringAndSize(string_buffer, acc)
    if needs_free:
        free(string_buffer)
    return result


cpdef list copy_sequence(list sequence):
    cdef:
        list result
        size_t i, n
        SequencePosition position
    n = PyList_Size(sequence)
    result = PyList_New(n)
    for i in range(n):
        position = <SequencePosition>PyList_GET_ITEM(sequence, i)
        position = position.copy()
        Py_INCREF(position)
        PyList_SET_ITEM(result, i, position)
    return result


@cython.freelist(100)
cdef class IsoformGenerator(object):
    cdef:
        public dict variable_mods
        public dict fixed_mods
        public list labels
        public long max_mods
        public bint override
        public bint show_unmodified_termini

        public bint persist

        # Temporary State
        public list sequence
        public list original_sequence
        public list last_indices

    def __init__(self, variable_mods=None, fixed_mods=None, labels=None, max_mods=None, override=False,
                 show_unmodified_termini=False, persist=True):
        if variable_mods is None:
            variable_mods = {}
        if fixed_mods is None:
            fixed_mods = {}
        if labels is None:
            labels = list(std_labels)
        if max_mods is None:
            max_mods = -1
        self.variable_mods = self.coerce_modification_dict(variable_mods)
        self.fixed_mods = self.coerce_modification_dict(fixed_mods)
        self.labels = labels
        self.max_mods = PyInt_AsLong(max_mods)
        self.override = override
        self.show_unmodified_termini = show_unmodified_termini
        self.persist = persist

        self.sequence = None
        self.original_sequence = None
        self.last_indices = []

    cpdef dict coerce_modification_dict(self, dict source):
        cdef:
            str label
            ModificationRule rule
            object targets
            dict result
        result = {}
        for label, targets in source.items():
            rule = ModificationRule(label, list(targets))
            PyDict_SetItem(result, label, rule)
        return result

    cpdef list initialize_from_str(self, str sequence):
        labels = self.labels + list(self.fixed_mods)
        self.sequence = SequencePosition.parse(sequence, labels=labels)
        return self.sequence

    cpdef list initialize_from_iterable(self, object sequence):
        self.sequence = SequencePosition.from_iterable(sequence)
        return self.sequence

    cpdef list apply_fixed_mods(self, list sequence):
        cdef:
            size_t i, n
            Py_ssize_t j
            PyObject* pkey
            PyObject* pvalue
            SequencePosition position
            ModificationRule rule

        n = PyList_Size(sequence)
        for i in range(n):
            position = <SequencePosition>PyList_GET_ITEM(sequence, i)
            j = 0
            while PyDict_Next(self.fixed_mods, &j, &pkey, &pvalue):
                rule = <ModificationRule>pvalue
                if rule.test(position):
                    if not position.is_modified():
                        position.add_modification(rule.modification)
        return sequence

    cpdef dict build_site_map(self, list sequence, bint include_unmodified=True):
        cdef:
            size_t i, n
            Py_ssize_t j
            PyObject* pkey
            PyObject* pvalue
            PyObject* psites
            list sites
            SequencePosition position
            ModificationRule rule
            dict site_map

        site_map = dict()
        n = PyList_Size(sequence)
        for i in range(n):
            position = <SequencePosition>PyList_GET_ITEM(sequence, i)
            j = 0
            while PyDict_Next(self.variable_mods, &j, &pkey, &pvalue):
                rule = <ModificationRule>pvalue
                if rule.test(position):
                    if not position.is_modified():
                        psites = PyDict_GetItem(site_map, i)
                        if psites == NULL:
                            sites = []
                            PyDict_SetItem(site_map, i, sites)
                        else:
                            sites = <list>psites
                        position = position.copy()
                        position.modification = rule.modification
                        sites.append(position)
        if include_unmodified:
            j = 0
            while PyDict_Next(site_map, &j, &pkey, &pvalue):
                position = <SequencePosition>PyList_GET_ITEM(sequence, PyInt_AsLong(<object>pkey))
                position = position.copy()
                (<list>pvalue).append(position)
        return site_map

    @cython.boundscheck(False)
    def generate(self, sequence_spec):
        cdef:
            list sequence
            dict variable_mod_sites
            list indices, assignments
            tuple assignment
            tuple subindices, subassignments
            object temp1, temp2
            object assignments_product
            size_t i, n, j, k, n_mods
            long* pindices

        if isinstance(sequence_spec, str):
            self.initialize_from_str(sequence_spec)
        else:
            self.initialize_from_iterable(sequence_spec)

        self.apply_fixed_mods(self.sequence)
        variable_mod_sites = self.build_site_map(self.sequence)
        indices = PyDict_Keys(variable_mod_sites)
        assignments = PyDict_Values(variable_mod_sites)
        k = PyList_Size(indices)

        if self.persist:
            self.original_sequence = list(self.sequence)

        if self.max_mods < 0 or self.max_mods > k:
            self.last_indices = indices
            assignments_product = it.product(*assignments)
            pindices = <long*>malloc(sizeof(long) * k)
            for i in range(k):
                pindices[i] = PyInt_AsLong(<object>PyList_GET_ITEM(indices, i))
            for a in assignments_product:
                assignment = <tuple>a
                yield self._apply_variable_assignment(self.sequence, pindices, assignment, k)
            free(pindices)
        else:
            self._reset_sequence()
            for n_mods in range(0, self.max_mods):
                for subindices in it.combinations(indices, n_mods):
                    self._reset_sequence()
                    self.last_indices = list(subindices)
                    k = PyTuple_Size(subindices)
                    pindices = <long*>malloc(sizeof(long) * k)
                    subassignments = PyTuple_New(k)
                    for i in range(k):
                        temp1 = <object>PyTuple_GET_ITEM(subindices, i)
                        pindices[i] = PyInt_AsLong(temp1)
                        temp2 = <object>PyDict_GetItem(variable_mod_sites, temp1)
                        Py_INCREF(temp2)
                        PyTuple_SET_ITEM(subassignments, i, temp2)
                    assignments_product = it.product(*subassignments)
                    for a in assignments_product:
                        yield self._apply_variable_assignment(self.sequence, pindices, <tuple>a, k)
                    free(pindices)
        self.sequence = None

    @cython.boundscheck(False)
    cdef list _apply_variable_assignment(self, list sequence, long* indices, tuple_or_list assignment, size_t n):
        cdef:
            size_t i
            long ix
            SequencePosition position
            str mod
        if self.persist:
            sequence = list(sequence)
        else:
            self._reset_sequence()
        for i in range(n):
            if tuple_or_list is tuple:
                position = <SequencePosition>PyTuple_GET_ITEM(assignment, i)
            else:
                position = <SequencePosition>PyList_GET_ITEM(assignment, i)
            ix = indices[i]
            Py_INCREF(position)
            PyList_SET_ITEM(sequence, ix, position)
        return sequence

    cdef void _reset_sequence(self):
        cdef:
            size_t i, n
            long ix
            SequencePosition position

        n = PyList_Size(self.last_indices)
        for i in range(n):
            ix = PyInt_AsLong(<object>PyList_GET_ITEM(self.last_indices, i))
            position = <SequencePosition>PyList_GET_ITEM(self.original_sequence, ix)
            Py_INCREF(position)
            PyList_SET_ITEM(self.sequence, ix, position)

    def __call__(self, sequence_spec):
        return self.generate(sequence_spec)


cdef bint _is_nterm_mod(str target):
    return target.startswith("nterm")


cdef bint _is_cterm_mod(str target):
    return target.startswith("cterm")


@cython.freelist(100)
cdef class ModificationRule(object):
    cdef:
        public str modification
        public list targets

    def __init__(self, str label, targets):
        self.modification = label
        self.targets = [ModificationTarget._from_string(x) if isinstance(x, str)
                        else ModificationTarget._create(None, PositionType.universal)
                        for x in targets]

    cpdef bint test(self, SequencePosition position):
        cdef:
            size_t i, n
            ModificationTarget t
        n = PyList_Size(self.targets)
        for i in range(n):
            t = <ModificationTarget>PyList_GET_ITEM(self.targets, i)
            if t.test(position):
                return True
        return False

    def __call__(self, position):
        return self.test(position)

    def __repr__(self):
        t = "{self.__class__.__name__}({self.modification!r}, {self.targets})"
        return t.format(self=self)


@cython.freelist(1000)
cdef class ModificationTarget(object):
    cdef:
        public str target
        public PositionType position_type

    def __init__(self, target, position_type=PositionType.internal):
        self.target = target
        self.position_type = position_type

    cpdef bint test(self, SequencePosition position):
        if self.position_type == PositionType.universal:
            return True
        if self.position_type != PositionType.internal:
            if position.position_type != self.position_type:
                return False
        if self.target is not None:
            if self.target != position.amino_acid:
                return False
        return True

    def __call__(self, position):
        return self.test(position)

    def __repr__(self):
        return "{self.__class__.__name__}({self.target!r}, {self.position_type!r})".format(self=self)

    @staticmethod
    cdef ModificationTarget _create(str target, PositionType position_type):
        cdef ModificationTarget inst = ModificationTarget.__new__(ModificationTarget)
        inst.target = target
        inst.position_type = position_type
        return inst

    @staticmethod
    def from_string(str string):
        return ModificationTarget._from_string(string)

    @staticmethod
    cdef ModificationTarget _from_string(str string):
        cdef:
            PositionType position_type
            size_t n
            str target
            ModificationTarget inst
        n = len(string)
        target = None
        position_type = PositionType.internal
        if _is_cterm_mod(string):
            position_type = PositionType.c_term
            if n > 5:
                if n > 6:
                    raise ValueError("Terminal Modifications with Target must be a single Amino Acid, got %s" % (string, ))
                target = string[-1]
        elif _is_nterm_mod(string):
            position_type = PositionType.n_term
            if n > 5:
                if n > 6:
                    raise ValueError("Terminal Modifications with Target must be a single Amino Acid, got %s" % (string, ))
                target = string[-1]
        else:
            if n > 1:
                raise ValueError("Target must be a single Amino Acid, got %s" % (string, ))
            target = string
        inst = ModificationTarget._create(target, position_type)
        return inst


def isoforms(str sequence, **kwargs):
    cdef:
        IsoformGenerator isoform_generator
        str _format
        object gen

    variable_mods = kwargs.get("variable_mods", None)
    fixed_mods = kwargs.get("fixed_mods", None)
    labels = kwargs.get("labels", None)
    max_mods = kwargs.get("max_mods", None)
    override = kwargs.get("override", False)
    show_unmodified_termini = kwargs.get("show_unmodified_termini", False)
    _format = kwargs.get("format", "str")

    isoform_generator = IsoformGenerator(
        variable_mods=variable_mods,
        fixed_mods=fixed_mods,
        labels=labels,
        max_mods=max_mods,
        override=override,
        show_unmodified_termini=show_unmodified_termini)
    gen = isoform_generator.generate(sequence)
    if _format == "str":
        for form in gen:
            yield _sequence_tostring(
                form, show_unmodified_termini=show_unmodified_termini)
    elif _format == "split":
        for form in gen:
            yield form
    else:
        raise PyteomicsError('Unsupported value of "format": {}'.format(_format))
