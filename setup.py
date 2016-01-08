#!/usr/bin/env python
from setuptools import setup, Extension, find_packages

try:
    from Cython.Build import cythonize
    extensions = cythonize([
        Extension(name="pyteomics.cythonize.cparser", sources=["pyteomics/cythonize/cparser.pyx"]),
        Extension(name="pyteomics.cythonize.cmass", sources=["pyteomics/cythonize/cmass.pyx"])
        ])
except ImportError:
    extensions = ([
        Extension(name="pyteomics.cythonize.cparser", sources=["pyteomics/cythonize/cparser.c"]),
        Extension(name="pyteomics.cythonize.cmass", sources=["pyteomics/cythonize/cmass.c"])
        ])


setup(
    name='pyteomics.cythonize',
    description='An Cython-accelerated version of common pyteomics functions',
    long_description=open("README.rst").read(),
    version="0.1.0",
    packages=find_packages(),
    zip_safe=False,
    install_requires=['pyteomics'],
    ext_modules=extensions,
    maintainer='Joshua Klein',
    maintainer_email="jaklein@bu.edu",
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Education',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries'
    ],
    namespace_packages=["pyteomics", "pyteomics.cythonize"],
    license='License :: OSI Approved :: Apache Software License'
    )
