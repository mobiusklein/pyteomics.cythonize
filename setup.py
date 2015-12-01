#!/usr/bin/env python

'''
setup.py file for pyteomics
'''

from setuptools import setup, Extension, find_packages

try:
    from Cython.Build import cythonize
    extensions = cythonize([
        Extension(name="pyteomics.cythonize.cparser", sources=["pyteomics/cythonize/cparser.pyx"]),
        Extension(name="pyteomics.cythonize.cmass", sources=["pyteomics/cythonize/cmass.pyx"])
        ])
except ImportError:
    extensions = ([
        Extension(name="pyteomics.cparser", sources=["pyteomics/cythonize/cparser.c"]),
        Extension(name="pyteomics.cmass", sources=["pyteomics/cythonize/cmass.c"])
        ])


setup(
    name='pyteomics.cythonize',
    version=0.1,
    packages=find_packages(),
    zip_safe=False,
    ext_modules=extensions,
    namespace_packages=["pyteomics", "pyteomics.cythonize"]
    )
