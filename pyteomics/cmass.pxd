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

    cdef double _calculate_mass(CComposition composition,
                                int average=*, charge=*, mass_data=*,
                                ion_type=*) except -1
    cdef str _parse_isotope_string(str label, int* isotope_num)
    cdef str _make_isotope_string(str element_name, int isotope_num)

