from pyteomics.parser import std_labels
from pyteomics import cparser

ref = [('H-', 'z', 'P'), ('E',), ('P',), ('z', 'T'), ('I',), ('D',), ('z', 'E', '-OH')]
out = cparser.amino_acid_composition('zPEPzTIDzE', labels=std_labels + ['z'])
print ref
print out
