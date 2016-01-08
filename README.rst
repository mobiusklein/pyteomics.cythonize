This package re-implements several of :title-reference:`Pyteomics` functions
in C using :title-reference:`Cython` and the Python-C API. Currently, only 
commonly used functions in `pyteomics.mass` and `pyteomics.parser` are implemented,
providing faster sequence manipulation and mass calculations. Every effort has been 
made to make the user-facing interfaces identical to their pure Python counterparts.

These functions are also exposed in the package's C-API so that other C-Extensions can make
use of them.

This package also re-implements the `pyteomics.mass.Composition` type in C using Python's
`dict` as a base.
