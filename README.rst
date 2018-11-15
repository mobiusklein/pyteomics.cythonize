This package re-implements several of :title-reference:`Pyteomics` functions
in C using :title-reference:`Cython` and the Python-C API. Currently, only 
commonly used functions in `pyteomics.mass` and `pyteomics.parser` are implemented,
providing faster sequence manipulation and mass calculations. Every effort has been 
made to make the user-facing interfaces identical to their pure Python counterparts.

These functions are also exposed in the package's C-API so that other C-Extensions can make
use of them.

This package also re-implements the `pyteomics.mass.Composition` type in C using Python's
`dict` as a base.


API
---

This package provides two modules, :title-reference:`pyteomics.cmass` and :title-reference:`pyteomics.cparser`,
which mimic a subset of the APIs of :title-reference:`pyteomics.mass` and :title-reference:`pyteomics.parser`
respectively. For example:

.. code:: python

    from pyteomics import cmass, mass

    assert cmass.fast_mass("PEPTIDE") == mass.fast_mass("PEPTIDE")

