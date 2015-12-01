

cpdef double fast_mass(str sequence, str ion_type=*, int charge=*,
                       dict mass_data=*, dict aa_mass=*,
                       dict ion_comp=*)

cpdef double fast_mass2(str sequence, str ion_type=*, int charge=*,
                        dict mass_data=*, dict aa_mass=*,
                        dict ion_comp=*)

cdef class CComposition(dict):
    cdef object _mass
    cdef tuple _mass_args
    cpdef CComposition clone(self)
    cpdef double mass(self, int average=?, charge=?, dict mass_data=?, ion_type=?) except -1
    cpdef _from_formula(self, str formula, dict mass_data)
    cpdef _from_dict(self, comp)
    cdef long getitem(self, str elem)
    cdef void setitem(self, str elem, long val)

cdef: 
    str _atom
    str _formula
    str _isotope_string

    object isotope_pattern
    object formula_pattern

    cdef inline str _parse_isotope_string(str label, int* isotope_num)
    cdef inline str _make_isotope_string(str element_name, int isotope_num)

