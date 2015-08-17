cimport cython
from cpython.ref cimport PyObject
from cpython.dict cimport PyDict_GetItem, PyDict_Next, PyDict_SetItem
from cpython.int cimport PyInt_AsLong, PyInt_Check
from cpython.float cimport PyFloat_AsDouble
from cpython.tuple cimport PyTuple_GetItem, PyTuple_GET_ITEM
from cpython.sequence cimport PySequence_GetItem
from cpython.exc cimport PyErr_Occurred

from pyteomics.auxiliary import PyteomicsError, _nist_mass
from pyteomics.mass import std_aa_mass as _std_aa_mass, std_ion_comp as _std_ion_comp, std_aa_comp as _std_aa_comp

import cparser
from cparser import parse, amino_acid_composition, _split_label
cimport cparser
from cparser cimport parse, amino_acid_composition, _split_label

cdef:
    dict nist_mass = _nist_mass
    dict std_aa_mass = _std_aa_mass
    dict std_ion_comp = _std_ion_comp
    dict std_aa_comp = _std_aa_comp


cdef inline double get_mass(dict mass_data, object key):
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
                       dict mass_data=nist_mass, dict aa_mass=std_aa_mass,
                       dict ion_comp=std_ion_comp):
    cdef:
        dict icomp
        double mass = 0
        int i, num
        Py_ssize_t pos
        str a
        PyObject* pkey
        PyObject* pvalue

    for i in range(len(sequence)):
        a = PySequence_GetItem(sequence, i)
        pvalue = PyDict_GetItem(aa_mass, a)
        if pvalue == NULL:
            raise PyteomicsError('No mass data for residue: ' + a)
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
                        dict mass_data=nist_mass, dict aa_mass=std_aa_mass,
                        dict ion_comp=std_ion_comp):
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
        comp = amino_acid_composition(sequence,
                show_unmodified_termini=1,
                allow_unknown_modifications=1,
                labels=list(aa_mass))
    except PyteomicsError:
        raise PyteomicsError('Mass not specified for label(s): {}'.format(
            ', '.join(set(parse(sequence)).difference(aa_mass))))

    mass = 0.
    pos = 0
    while(PyDict_Next(comp, &pos, &pkey, &pvalue)):
        aa = <str>pkey
        num = <int>pvalue
        if aa in aa_mass:
            ptemp = PyDict_GetItem(aa_mass, aa)
            mass += PyFloat_AsDouble(<object>ptemp) * num
        else:
            temp = _split_label(aa)
            mod = <str>PyTuple_GET_ITEM(temp, 0)
            X = <str>PyTuple_GET_ITEM(temp, 1)
            ptemp = PyDict_GetItem(aa_mass, mod)
            if ptemp is NULL:
                raise (<object>ptemp)("An error occurred in cmass.fast_mass: %s not found in aa_mass" % mod)
            interim = PyFloat_AsDouble(<object>ptemp)
            ptemp = PyDict_GetItem(aa_mass, X)
            if ptemp is NULL:
                raise (<object>ptemp)("An error occurred in cmass.fast_mass: %s not found in aa_mass" % X)
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
