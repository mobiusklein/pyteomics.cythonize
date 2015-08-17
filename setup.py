#!/usr/bin/env python

'''
setup.py file for pyteomics
'''

from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize


extensions = cythonize([
    Extension(name="pyteomics.cparser", sources=["pyteomics/cparser.pyx"]),
    Extension(name="pyteomics.cmass", sources=["pyteomics/cmass.pyx"])
    ])

setup(
    name='pyteomics.cythonize',
    version=0.1,
    packages=find_packages(),
    zip_safe=False,
    ext_modules=extensions
    )
