
echo "amino_acid_composition"
echo "cparser"
python -m timeit -s "from pyteomics import cparser" "cparser.amino_acid_composition('PEPpTIDE', allow_unknown_modifications=True, show_unmodified_termini=True)"
echo "parser"
python -m timeit -s "from pyteomics import parser" "parser.amino_acid_composition('PEPpTIDE', allow_unknown_modifications=True, show_unmodified_termini=True)"
