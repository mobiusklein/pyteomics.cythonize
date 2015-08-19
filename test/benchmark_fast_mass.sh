#!/usr/bin/env bash
echo "fast_mass"
echo cmass
python -m timeit -s "from pyteomics import cmass" "cmass.fast_mass('PEPTIDE', ion_type='y')"
echo mass
python -m timeit -s "from pyteomics import mass" "mass.fast_mass('PEPTIDE', ion_type='y')"
