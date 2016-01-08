from pyteomics import mass

try:
    from . import cmass
    fast_mass = cmass.fast_mass
    fast_mass2 = cmass.fast_mass2
    Composition = cmass.CComposition
    calculate_mass = cmass.calculate_mass
    nist_mass, std_aa_mass, std_ion_comp, std_aa_comp = cmass.__get_constants()

except ImportError:
    from pyteomics.mass import (
        nist_mass, std_aa_mass, std_ion_comp, std_aa_comp)

    fast_mass = mass.fast_mass
    fast_mass2 = mass.fast_mass2
    Composition = mass.Composition
    calculate_mass = mass.calculate_mass
