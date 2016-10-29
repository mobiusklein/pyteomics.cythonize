pyteomics-cythonize package
===========================

This package re-implements several of :title-reference:`Pyteomics` functions
in C using :title-reference:`Cython` and the Python-C API. Currently, only 
commonly used functions in `pyteomics.mass` and `pyteomics.parser` are implemented,
providing faster sequence manipulation and mass calculations. Every effort has been 
made to make the user-facing interfaces identical to their pure Python counterparts.

These functions are also exposed in the package's C-API so that other C-Extensions can make
use of them.

This package also re-implements the :class:`pyteomics.mass.Composition` type in C using Python's
:class:`dict` as a base.

cmass module
----------------------

.. automodule:: cmass
    :members:
    :undoc-members:
    :show-inheritance:

    .. autofunction:: cmass.fast_mass

    .. autofunction:: cmass.fast_mass2

    .. autofunction:: cmass.calculate_mass

    .. autoclass:: CComposition

        .. py:method:: CComposition.mass(self, int average=False, charge=None, dict mass_data=nist_mass, ion_type=None)

            Calculate the mass or m/z of a Composition.

        .. py:method:: clone(self)

            Copy this instance


cparser module
------------------------

.. automodule:: cparser
    :members:
    :undoc-members:
    :show-inheritance:
