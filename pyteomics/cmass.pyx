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

import re

cimport cython
from cpython.ref cimport PyObject
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem, PyDict_Next, PyDict_Keys, PyDict_Update, PyDict_DelItem
from cpython.int cimport PyInt_AsLong, PyInt_Check, PyInt_FromLong
from cpython.tuple cimport PyTuple_GetItem, PyTuple_GET_ITEM
from cpython.list cimport PyList_GET_ITEM
from cpython.float cimport PyFloat_AsDouble
from cpython.sequence cimport PySequence_GetItem
from cpython.exc cimport PyErr_Occurred

from pyteomics.auxiliary import PyteomicsError, _nist_mass
from pyteomics.mass import (
    std_aa_mass,
    std_ion_comp,
    std_aa_comp)

from collections import defaultdict
from itertools import chain

from pyteomics.parser import amino_acid_composition as pamino_acid_composition
from pyteomics import cparser
# from pyteomics cimport cparser
from pyteomics.cparser cimport parse, amino_acid_composition, _split_label

cdef dict nist_mass = _nist_mass

cdef:
    dict _std_aa_mass = std_aa_mass
    dict _std_ion_comp = {k: CComposition(v) for k, v in std_ion_comp.items()}
    dict _std_aa_comp = {k: CComposition(v) for k, v in std_aa_comp.items()}

std_ion_comp = _std_ion_comp
std_aa_comp = _std_aa_comp


def __get_constants():
    return nist_mass, std_aa_mass, std_ion_comp, std_aa_comp


cdef inline double get_mass(dict mass_data, object key):
    '''
    Internal method to do neutral monoisotopic mass look ups

    Parameters
    ----------
    mass_data: dict
        The dictionary of symbol-to-mass mappings
    key: object
        The symbol to look up
    '''
    cdef:
        PyObject* interm
        double mass

    interim = PyDict_GetItem(mass_data, key)
    if interim == NULL:
        raise KeyError(key)
    interim = PyDict_GetItem(<dict>interim, 0)
    if interim == NULL:
        raise KeyError(0)
    mass = PyFloat_AsDouble(<object>PyTuple_GetItem(<tuple>interim, 0))
    return mass


cpdef double fast_mass(str sequence, str ion_type=None, int charge=0,
                       dict mass_data=_nist_mass, dict aa_mass=_std_aa_mass,
                       dict ion_comp=_std_ion_comp):
    """Calculate monoisotopic mass of an ion using the fast
    algorithm. May be used only if amino acid residues are presented in
    one-letter code.

    Parameters
    ----------
    sequence : str
        A polypeptide sequence string.
    ion_type : str, optional
        If specified, then the polypeptide is considered to be
        in a form of corresponding ion. Do not forget to
        specify the charge state!
    charge : int, optional
        If not 0 then m/z is calculated: the mass is increased
        by the corresponding number of proton masses and divided
        by z.
    mass_data : dict, optional
        A dict with the masses of chemical elements (the default
        value is :py:data:`nist_mass`).
    aa_mass : dict, optional
        A dict with the monoisotopic mass of amino acid residues
        (default is std_aa_mass);
    ion_comp : dict, optional
        A dict with the relative elemental compositions of peptide ion
        fragments (default is :py:data:`std_ion_comp`).

    Returns
    -------
    mass : float
        Monoisotopic mass or m/z of a peptide molecule/ion.
    """
    cdef:
        CComposition icomp
        double mass = 0
        int i, num
        Py_ssize_t pos
        object a
        PyObject* pkey
        PyObject* pvalue

    for i in range(len(sequence)):
        a = PySequence_GetItem(sequence, i)
        pvalue = PyDict_GetItem(aa_mass, a)
        if pvalue == NULL:
            raise PyteomicsError('No mass data for residue: %r' % (a,))
        mass += PyFloat_AsDouble(<object>pvalue)
    pvalue = PyErr_Occurred()
    if pvalue != NULL:
        raise (<object>pvalue)("An error occurred in cmass.fast_mass")
    mass += get_mass(mass_data, 'H') * 2 + get_mass(mass_data, 'O')

    if ion_type:
        try:
            icomp = ion_comp[ion_type]
        except KeyError:
            raise PyteomicsError('Unknown ion type: {}'.format(ion_type))
        pos = 0
        while(PyDict_Next(icomp, &pos, &pkey, &pvalue)):
            mass += get_mass(mass_data, <object>pkey) * PyFloat_AsDouble(<object>pvalue)
        pvalue = PyErr_Occurred()
        if pvalue != NULL:
            raise (<object>pvalue)("An error occurred in cmass.fast_mass")

    if charge:
        mass = (mass + get_mass(mass_data, 'H+') * charge) / charge

    return mass


cpdef double fast_mass2(str sequence, str ion_type=None, int charge=0,
                        dict mass_data=_nist_mass, dict aa_mass=_std_aa_mass,
                        dict ion_comp=_std_ion_comp):
    """Calculate monoisotopic mass of an ion using the fast
    algorithm. *modX* notation is fully supported.

    Parameters
    ----------
    sequence : str
        A polypeptide sequence string.
    ion_type : str, optional
        If specified, then the polypeptide is considered to be
        in a form of corresponding ion. Do not forget to
        specify the charge state!
    charge : int, optional
        If not 0 then m/z is calculated: the mass is increased
        by the corresponding number of proton masses and divided
        by z.
    mass_data : dict, optional
        A dict with the masses of chemical elements (the default
        value is :py:data:`nist_mass`).
    aa_mass : dict, optional
        A dict with the monoisotopic mass of amino acid residues
        (default is std_aa_mass);
    ion_comp : dict, optional
        A dict with the relative elemental compositions of peptide ion
        fragments (default is :py:data:`std_ion_comp`).

    Returns
    -------
    mass : float
        Monoisotopic mass or m/z of a peptide molecule/ion.
    """
    cdef:
        dict comp
        str aa, mod, X, element
        int num
        double mass, interim
        tuple temp
        CComposition icomp
        Py_ssize_t pos
        PyObject* pkey
        PyObject* pvalue
        PyObject* ptemp

    ptemp = PyDict_GetItem(aa_mass, 'H-')
    if ptemp == NULL:
        PyDict_SetItem(aa_mass, 'H-', get_mass(mass_data, "H"))
    ptemp = PyDict_GetItem(aa_mass, '-OH')
    if ptemp == NULL:
        PyDict_SetItem(aa_mass, '-OH', get_mass(mass_data, "H") + get_mass(mass_data, "O"))
    try:
        comp = (cparser.amino_acid_composition(
        # comp = (pamino_acid_composition(
            sequence,
            show_unmodified_termini=1,
            allow_unknown_modifications=1,
            labels=list(aa_mass)))
    except PyteomicsError:
        raise PyteomicsError('Mass not specified for label(s): {}'.format(
            ', '.join(set(cparser.parse(sequence)).difference(aa_mass))))
    mass = 0.
    pos = 0
    while(PyDict_Next(comp, &pos, &pkey, &pvalue)):
        aa = <str>pkey
        num = PyInt_AsLong(<object>pvalue)
        if aa in aa_mass:
            ptemp = PyDict_GetItem(aa_mass, aa)
            mass += PyFloat_AsDouble(<object>ptemp) * num
        else:
            temp = _split_label(aa)
            mod = <str>PyTuple_GET_ITEM(temp, 0)
            X = <str>PyTuple_GET_ITEM(temp, 1)
            ptemp = PyDict_GetItem(aa_mass, mod)
            if ptemp == NULL:
                raise KeyError("An error occurred in cmass.fast_mass: %s not found in aa_mass" % mod)
            interim = PyFloat_AsDouble(<object>ptemp)
            ptemp = PyDict_GetItem(aa_mass, X)
            if ptemp == NULL:
                raise KeyError("An error occurred in cmass.fast_mass: %s not found in aa_mass" % X)
            interim += PyFloat_AsDouble(<object>ptemp)
            mass += interim * num

    if ion_type:
        try:
            icomp = ion_comp[ion_type]
        except KeyError:
            raise PyteomicsError('Unknown ion type: {}'.format(ion_type))

        pos = 0
        while(PyDict_Next(icomp, &pos, &pkey, &pvalue)):
            mass += get_mass(mass_data, <object>pkey) * PyFloat_AsDouble(<object>pvalue)
        pvalue = PyErr_Occurred()
        if pvalue != NULL:
            raise (<object>pvalue)("An error occurred in cmass.fast_mass")

    if charge:
        mass = (mass + get_mass(mass_data, 'H+') * charge) / charge

    return mass


# Forward Declaration
cdef: 
    str _atom = r'([A-Z][a-z+]*)(?:\[(\d+)\])?([+-]?\d+)?'
    str _formula = r'^({})*$'.format(_atom)
    str _isotope_string = r'^([A-Z][a-z+]*)(?:\[(\d+)\])?$'

    object isotope_pattern = re.compile(_isotope_string)
    object formula_pattern = re.compile(_formula)


@cython.boundscheck(False)
cdef str _parse_isotope_string(str label, int* isotope_num):
    '''Parses an isotope string and extracts the element name and isotope number.
    The element name is returned, but the isotope number is returned by indirection.

    Parameters
    ----------
    label: str
        An isotope string of the format "<element><[isotope_number]>"
    isotope_num: int_ptr
        A pointer to an integer which will contain the isotope number

    Returns
    -------
    str: Element Name
    '''
    cdef:
        # int isotope_num = 0
        int i = 0
        int in_bracket = False
        # str element_name
        str current
        list name_parts = []
        list num_parts = []
        #Isotope result
    for i in range(len(label)):
        current = label[i]
        if in_bracket:
            if current == "]":
                break
            num_parts.append(current)
        elif current == "[":
            in_bracket = True
        else:
            name_parts.append(current)
    element_name = (''.join(name_parts))
    if len(num_parts) > 0:
        isotope_num[0] = (int(''.join(num_parts)))
    else:
        isotope_num[0] = 0
    return element_name


cdef str _make_isotope_string(str element_name, int isotope_num):
    """Form a string label for an isotope."""
    cdef:
        tuple parts
    if isotope_num == 0:
        return element_name
    else:
        parts = (element_name, isotope_num)
        return '%s[%d]' % parts


def marshal_ccomposition(state):
    return CComposition(state)


cdef class CComposition(dict):
    """
    A Composition object stores a chemical composition of a
    substance. Basically it is a dict object, in which keys are the names
    of chemical elements and values contain integer numbers of
    corresponding atoms in a substance.

    The main improvement over dict is that Composition objects allow
    addition and subtraction.

    If ``formula`` is not specified, the constructor will look at the first
    positional argument and try to build the object from it. Without
    positional arguments, a Composition will be constructed directly from
    keyword arguments.

    Parameters
    ----------
    formula : str, optional
        A string with a chemical formula. All elements must be present in
        `mass_data`.
    mass_data : dict, optional
        A dict with the masses of chemical elements (the default
        value is :py:data:`nist_mass`). It is used for formulae parsing only.
    """

    def _from_parsed_sequence(self, parsed_sequence, aa_comp):
        self.clear()
        comp = defaultdict(int)
        for aa in parsed_sequence:
            if aa in aa_comp:
                for elem, cnt in aa_comp[aa].items():
                    comp[elem] += cnt
            else:
                try:
                    mod, aa = cparser._split_label(aa)
                    for elem, cnt in chain(
                            aa_comp[mod].items(), aa_comp[aa].items()):
                        comp[elem] += cnt

                except (PyteomicsError, KeyError):
                    raise PyteomicsError(
                            'No information for %s in `aa_comp`' % aa)
        self._from_dict(comp)

    def _from_split_sequence(self, split_sequence, aa_comp):
        self.clear()
        comp = defaultdict(int)
        for group in split_sequence:
            i = 0
            while i < len(group):
                for j in range(len(group)+1, -1, -1):
                    try:
                        label = ''.join(group[i:j])
                        for elem, cnt in aa_comp[label].items():
                            comp[elem] += cnt
                    except KeyError:
                        continue
                    else:
                        i = j
                        break
                if j == 0:
                    raise PyteomicsError("Invalid group starting from "
                            "position %d: %s" % (i+1, group))
        self._from_dict(comp)

    def _from_sequence(self, sequence, aa_comp):
        parsed_sequence = cparser.parse(
            sequence,
            labels=aa_comp,
            show_unmodified_termini=True)
        self._from_parsed_sequence(parsed_sequence, aa_comp)

    def __str__(self):   # pragma: no cover
        return 'Composition({})'.format(dict.__repr__(self))

    def __repr__(self):  # pragma: no cover
        return str(self)

    def __iadd__(CComposition self, other):
        cdef:
            str elem
            long cnt
            PyObject *pkey
            PyObject *pvalue
            Py_ssize_t ppos = 0

        while(PyDict_Next(other, &ppos, &pkey, &pvalue)):
            elem = <str>pkey
            cnt = self.getitem(elem)
            self.setitem(elem, cnt + PyInt_AsLong(<object>pvalue))

        self._mass_args = None
        return self


    def __add__(self, other):
        cdef:
            str elem
            long cnt
            CComposition result
            PyObject *pkey
            PyObject *pvalue
            Py_ssize_t ppos = 0
        if not isinstance(self, CComposition):
            other, self = self, other
        result = CComposition(self)
        while(PyDict_Next(other, &ppos, &pkey, &pvalue)):
            elem = <str>pkey
            cnt = result.getitem(elem)
            cnt += PyInt_AsLong(<object>pvalue)
            result.setitem(elem, cnt)

        return result

    def __isub__(self, other):
        cdef:
            str elem
            long cnt
            PyObject *pkey
            PyObject *pvalue
            Py_ssize_t ppos = 0

        while(PyDict_Next(other, &ppos, &pkey, &pvalue)):
            elem = <str>pkey
            cnt = self.getitem(elem)
            self.setitem(elem, cnt - PyInt_AsLong(<object>pvalue))

        self._mass_args = None
        return self

    def __sub__(self, other):
        cdef:
            str elem
            long cnt
            CComposition result
            PyObject *pkey
            PyObject *pvalue
            Py_ssize_t ppos = 0
        if not isinstance(self, CComposition):
            self = CComposition(self)
        result = CComposition(self)
        while(PyDict_Next(other, &ppos, &pkey, &pvalue)):
            elem = <str>pkey
            cnt = result.getitem(elem)
            cnt -= PyInt_AsLong(<object>pvalue)
            result.setitem(elem, cnt)

        return result

    def __reduce__(self):
        return marshal_ccomposition, (dict(self),)

    def __getstate__(self):
        return dict(self)

    def __setstate__(self, d):
        self._from_dict(d)
        self._mass = None
        self._mass_args = None


    def __mul__(self, other):
        cdef:
            CComposition prod = CComposition()
            int rep, v
            str k

        if isinstance(other, CComposition):
            self, other = other, self
        
        if not isinstance(other, int):
            raise PyteomicsError(
                'Cannot multiply Composition by non-integer',
                other)
        rep = other
        for k, v in self.items():
            prod.setitem(k, v * rep)
        return prod

    def __richcmp__(self, other, int code):
        if code == 2:
            if not isinstance(other, dict):
                return False
            self_items = set([i for i in self.items() if i[1]])
            other_items = set([i for i in other.items() if i[1]])
            return self_items == other_items
        else:
            return NotImplemented

    def __neg__(self):
        return self * -1

    # Override the default behavior, if a key is not present
    # do not initialize it to 0.
    def __missing__(self, str key):
        return 0

    def __setitem__(self, str key, object value):
        cdef long int_value = PyInt_AsLong(round(value))
        if int_value:  # Will not occur on 0 as 0 is falsey AND an integer
            self.setitem(key, int_value)
        elif key in self:
            del self[key]
        self._mass_args = None

    def copy(self):
        return self.__class__(self)

    cdef inline long getitem(self, str elem):
        cdef:
            PyObject* resobj
            long count
        resobj = PyDict_GetItem(self, elem)
        if (resobj == NULL):
            return 0
        count = PyInt_AsLong(<object>resobj)
        return count

    cdef inline void setitem(self, str elem, long val):
        PyDict_SetItem(self, elem, val)
        self._mass_args = None

    cpdef CComposition clone(self):
        '''Create a copy of this instance

        Returns
        -------
        CComposition
        '''
        return CComposition(self)

    def update(self, *args, **kwargs):
        dict.update(self, *args, **kwargs)
        self._mass_args = None

    @cython.boundscheck(False)
    cpdef _from_formula(self, str formula, dict mass_data):
        cdef:
            str elem, isotope, number
        if '(' in formula:
            self._from_formula_parens(formula, mass_data)
        elif not formula_pattern.match(formula):
            raise PyteomicsError('Invalid formula: ' + formula)
        else:
            for elem, isotope, number in re.findall(_atom, formula):
                if not elem in mass_data:
                    raise PyteomicsError('Unknown chemical element: ' + elem)
                self[_make_isotope_string(elem, int(isotope) if isotope else 0)
                        ] += int(number) if number else 1

    @cython.boundscheck(True)
    def _from_formula_parens(self, formula, mass_data):
        # Parsing a formula backwards.
        prev_chem_symbol_start = len(formula)
        i = len(formula) - 1

        seek_mode = 0
        parse_stack = ""
        resolve_stack = []
        group_coef = None

        while i >= 0:
            if seek_mode < 1:
                if (formula[i] == ")"):
                    seek_mode += 1
                    if i + 1 == prev_chem_symbol_start:
                        group_coef = 1
                    elif formula[i + 1].isdigit():
                        group_coef = int(formula[i + 1:prev_chem_symbol_start])
                    i -= 1
                    continue
                # Read backwards until a non-number character is met.
                if (formula[i].isdigit() or formula[i] == '-'):
                    i -= 1
                    continue

                else:
                    # If the number of atoms is omitted then it is 1.
                    if i + 1 == prev_chem_symbol_start:
                        num_atoms = 1
                    else:
                        try:
                            num_atoms = int(formula[i + 1:prev_chem_symbol_start])
                        except ValueError:
                            raise PyteomicsError(
                                'Badly-formed number of atoms: %s' % formula)

                    # Read isotope number if specified, else it is undefined (=0).
                    if formula[i] == ']':
                        brace_pos = formula.rfind('[', 0, i)
                        if brace_pos == -1:
                            raise PyteomicsError(
                                'Badly-formed isotope number: %s' % formula)
                        try:
                            isotope_num = int(formula[brace_pos + 1:i])
                        except ValueError:
                            raise PyteomicsError(
                                'Badly-formed isotope number: %s' % formula)
                        i = brace_pos - 1
                    else:
                        isotope_num = 0

                    # Match the element name to the mass_data.
                    element_found = False
                    # Sort the keys from longest to shortest to workaround
                    # the overlapping keys issue
                    for element_name in sorted(mass_data, key=len, reverse=True):
                        if formula.endswith(element_name, 0, i + 1):
                            isotope_string = _make_isotope_string(
                                element_name, isotope_num)
                            self[isotope_string] += num_atoms
                            i -= len(element_name)
                            prev_chem_symbol_start = i + 1
                            element_found = True
                            break

                    if not element_found:
                        raise PyteomicsError(
                            'Unknown chemical element in the formula: %s' % formula)
            else:
                ch = formula[i]
                parse_stack += ch
                i -= 1
                if(ch == "("):
                    seek_mode -= 1
                    if seek_mode == 0:

                        resolve_stack.append(Composition(
                                             # Omit the last character, then reverse the parse
                                             # stack string.
                                             formula=parse_stack[:-1][::-1],
                                             mass_data=mass_data)
                                             * group_coef)
                        prev_chem_symbol_start = i + 1
                        seek_mode = False
                        parse_stack = ""
                elif(formula[i] == ")"):
                    seek_mode += 1
                else:
                    # continue to accumulate tokens
                    pass

        # Unspool the resolve stack, adding together the chunks
        # at this level. __add__ operates immutably, so must manually
        # loop through each chunk.
        for chunk in resolve_stack:
            for elem, cnt in chunk.items():
                self[elem] += cnt

    cpdef _from_dict(self, comp):
        '''
        Directly overwrite this object's keys with the values in
        `comp` without checking their type.
        '''
        PyDict_Update(self, comp)


    cpdef double mass(self, int average=False, charge=None, dict mass_data=nist_mass, ion_type=None) except -1:
        '''
        Calculate the mass or m/z of a Composition.
        '''
        cdef object mdid
        mdid = id(mass_data)
        if self._mass_args is not None and average is self._mass_args[0]\
                and charge == self._mass_args[1] and mdid == self._mass_args[2]\
                and ion_type == self._mass_args[3]:
            return self._mass
        else:
            self._mass_args = (average, charge, mdid, ion_type)
            self._mass = _calculate_mass(composition=self, average=average, charge=charge, mass_data=mass_data, ion_type=ion_type)
            return self._mass

    def __init__(self, *args, **kwargs):
        dict.__init__(self)
        cdef:
            dict mass_data, aa_comp
            str kwa
            set kw_sources, kw_given
        aa_comp=kwargs.get('aa_comp', _std_aa_comp)
        mass_data=kwargs.get('mass_data')
        if mass_data is None:
            mass_data = nist_mass

        kw_sources = {'formula', 'sequence', 'parsed_sequence',
                'split_sequence'}
        kw_given = kw_sources.intersection(kwargs)
        if len(kw_given) > 1:
            raise PyteomicsError('Only one of {} can be specified!\n'
                    'Given: {}'.format(', '.join(kw_sources),
                        ', '.join(kw_given)))
        elif kw_given:
            kwa = kw_given.pop()
            getattr(self, '_from_' + kwa)(kwargs[kwa],
                    mass_data if kwa == 'formula' else aa_comp)

        # can't build from kwargs
        elif args:
            if isinstance(args[0], dict):
                self._from_dict(args[0])
            elif isinstance(args[0], str):
                try:
                    self._from_sequence(args[0], aa_comp)
                except PyteomicsError:
                    try:
                        self._from_formula(args[0], mass_data)
                    except PyteomicsError:
                        raise PyteomicsError(
                                'Could not create a Composition object from '
                                'string: "{}": not a valid sequence or '
                                'formula'.format(args[0]))
            else:
                try:
                    self._from_sequence(cparser.tostring(args[0], True),
                            aa_comp)
                except:
                    raise PyteomicsError('Could not create a Composition object'
                            ' from `{}`. A Composition object must be '
                            'specified by sequence, parsed or split sequence,'
                            ' formula or dict.'.format(args[0]))
        else:
            self._from_dict(kwargs)

        self._mass = None
        self._mass_args = None

Composition = CComposition

def calculate_mass(composition=None, average=False, charge=None, mass_data=None, ion_type=None, **kwargs):
    """Calculates the monoisotopic mass of a polypeptide defined by a
    sequence string, parsed sequence, chemical formula or
    Composition object.

    One or none of the following keyword arguments is required:
    **formula**, **sequence**, **parsed_sequence**, **split_sequence**
    or **composition**.
    All arguments given are used to create a :py:class:`Composition` object,
    unless an existing one is passed as a keyword argument.

    Note that if a sequence string is supplied and terminal groups are not
    explicitly shown, then the mass is calculated for a polypeptide with
    standard terminal groups (NH2- and -OH).

    .. warning::

        Be careful when supplying a list with a parsed sequence. It must be
        obtained with enabled `show_unmodified_termini` option.

    Parameters
    ----------
    formula : str, optional
        A string with a chemical formula.
    sequence : str, optional
        A polypeptide sequence string in modX notation.
    parsed_sequence : list of str, optional
        A polypeptide sequence parsed into a list of amino acids.
    composition : Composition, optional
        A Composition object with the elemental composition of a substance.
    aa_comp : dict, optional
        A dict with the elemental composition of the amino acids (the
        default value is std_aa_comp).
    average : bool, optional
        If :py:const:`True` then the average mass is calculated. Note that mass
        is not averaged for elements with specified isotopes. Default is
        :py:const:`False`.
    charge : int, optional
        If not 0 then m/z is calculated: the mass is increased
        by the corresponding number of proton masses and divided
        by `charge`.
    mass_data : dict, optional
        A dict with the masses of the chemical elements (the default
        value is :py:data:`nist_mass`).
    ion_comp : dict, optional
        A dict with the relative elemental compositions of peptide ion
        fragments (default is :py:data:`std_ion_comp`).
    ion_type : str, optional
        If specified, then the polypeptide is considered to be in the form
        of the corresponding ion. Do not forget to specify the charge state!

    Returns
    -------
    mass : float
    """
    if composition is None:
        composition = CComposition(mass_data=mass_data, **kwargs)
    return composition.mass(average=average, charge=charge, mass_data=mass_data, ion_type=ion_type)

    

@cython.wraparound(False)
@cython.boundscheck(False)
cdef double _calculate_mass(CComposition composition,
                                   int average=False, charge=None, mass_data=None,
                                   ion_type=None) except -1:
    """Calculates the monoisotopic mass of a CComposition object.

    Parameters
    ----------
    composition : CComposition
        A Composition object with the elemental composition of a substance. Exclusive with `formula`
    average : bool, optional
        If :py:const:`True` then the average mass is calculated. Note that mass
        is not averaged for elements with specified isotopes. Default is
        :py:const:`False`.
    charge : int, optional
        If not 0 then m/z is calculated: the mass is increased
        by the corresponding number of proton masses and divided
        by z.
    mass_data : dict, optional
        A dict with the masses of the chemical elements (the default
        value is :py:data:`nist_mass`).

    Returns
    -------
        mass : float
    """
    cdef:
        int old_charge, isotope_num, isotope, quantity
        double mass, isotope_mass, isotope_frequency
        long _charge
        str isotope_string, element_name
        dict mass_provider
        CComposition ion_type_comp
        list key_list
        PyObject* interm
        Py_ssize_t iter_pos = 0

    if mass_data is None:
        mass_provider = nist_mass
    else:
        mass_provider = mass_data

    # Get charge.
    if charge is None:
        charge = composition.getitem('H+')
    else:
        if charge != 0 and composition.getitem('H+') != 0:
            raise PyteomicsError("Charge is specified both by the number of protons and parameters")
    _charge = PyInt_AsLong(charge)
    old_charge = composition.getitem('H+')
    composition.setitem('H+', charge)

    # Calculate mass.
    mass = 0.0
    key_list = PyDict_Keys(composition)
    for iter_pos in range(len(key_list)):
        isotope_string = <str>PyList_GET_ITEM(key_list, iter_pos)
        # element_name, isotope_num = _parse_isotope_string(isotope_string)
        element_name = _parse_isotope_string(isotope_string, &isotope_num)

        # Calculate average mass if required and the isotope number is
        # not specified.
        if (not isotope_num) and average:
            for isotope in mass_provider[element_name]:
                if isotope != 0:
                    quantity = <int>composition.getitem(element_name)
                    isotope_mass = <double>mass_provider[element_name][isotope][0]
                    isotope_frequency = <double>mass_provider[element_name][isotope][1]

                    mass += quantity * isotope_mass * isotope_frequency
        else:
            interim = PyDict_GetItem(mass_provider, element_name)
            interim = PyDict_GetItem(<dict>interim, isotope_num)
            isotope_mass = PyFloat_AsDouble(<object>PyTuple_GetItem(<tuple>interim, 0))

            mass += (composition.getitem(isotope_string) * isotope_mass)

    if ion_type is not None:
        interm = PyDict_GetItem(_std_ion_comp, ion_type)
        if interm == NULL:
            raise KeyError("Unknown ion_type: {}".format(ion_type))
        ion_type_comp = <CComposition>interm
        key_list = PyDict_Keys(ion_type_comp)
        for iter_pos in range(len(key_list)):
            isotope_string = <str>PyList_GET_ITEM(key_list, iter_pos)
            element_name = _parse_isotope_string(isotope_string, &isotope_num)

            # Calculate average mass if required and the isotope number is
            # not specified.
            if (not isotope_num) and average:
                for isotope in mass_provider[element_name]:
                    if isotope != 0:
                        quantity = ion_type_comp.getitem(element_name)
                        isotope_mass = <double>mass_provider[element_name][isotope][0]
                        isotope_frequency = <double>mass_provider[element_name][isotope][1]

                        mass += quantity * isotope_mass * isotope_frequency
            else:
                interim = PyDict_GetItem(mass_provider, element_name)
                interim = PyDict_GetItem(<dict>interim, isotope_num)
                isotope_mass = PyFloat_AsDouble(<object>PyTuple_GetItem(<tuple>interim, 0))

                mass += (ion_type_comp.getitem(isotope_string) * isotope_mass)


    # Calculate m/z if required.
    if _charge != 0:
        mass /= abs(_charge)


    if old_charge != 0:
        composition.setitem('H+', old_charge)
    else:
        PyDict_DelItem(composition, "H+")

    return mass

Composition = CComposition
