cpdef list parse(str sequence, int show_unmodified_termini=*, int split=*,
                 int allow_unknown_modifications=*, object  labels=*)

cpdef dict amino_acid_composition(object sequence, int show_unmodified_termini=*, int term_aa=*,
                                  int allow_unknown_modifications=*, object labels=*)

cpdef tuple _split_label(str label)

