#!/usr/bin/env bash
echo "fast_mass2"
echo "cmass"
python -m timeit -s "from pyteomics import cmass" "cmass.fast_mass2('PEPTIDE', ion_type='y')"
echo "mass"
python -m timeit -s "from pyteomics import mass" "mass.fast_mass2('PEPTIDE', ion_type='y')"
