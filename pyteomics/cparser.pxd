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



cpdef list parse(str sequence, bint show_unmodified_termini=*, bint split=*,
                 bint allow_unknown_modifications=*, object  labels=*)

cpdef dict amino_acid_composition(object sequence, bint show_unmodified_termini=*, bint term_aa=*,
                                  bint allow_unknown_modifications=*, object labels=*)

cpdef tuple _split_label(str label)

cpdef str tostring(object parsed_sequence, bint show_unmodified_termini=*)
