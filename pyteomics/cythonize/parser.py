from pyteomics import parser

try:
    from . import cparser
    parse = cparser.parse
    tostring = cparser.tostring
    amino_acid_composition = cparser.amino_acid_composition
    length = cparser.length
except ImportError:
    parse = parser.parse
    tostring = parser.tostring
    amino_acid_composition = parser.amino_acid_composition
    length = parser.length


std_labels = parser.std_labels
cleave = parser.cleave
expasy_rules = parser.expasy_rules
isoforms = parser.isoforms
num_sites = parser.num_sites
valid = parser.valid
fast_valid = parser.fast_valid
