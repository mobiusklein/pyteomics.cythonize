#!/usr/bin/env python
from setuptools import setup, Extension, find_namespace_packages

try:
    from Cython.Build import cythonize
    extensions = cythonize([
        Extension(name="pyteomics.cparser", sources=[
                  "pyteomics/cparser.pyx"]),
        Extension(name="pyteomics.cmass",
                  sources=["pyteomics/cmass.pyx"])
    ])
except ImportError:
    extensions = ([
        Extension(name="pyteomics.cparser",
                  sources=["pyteomics/cparser.c"]),
        Extension(name="pyteomics.cmass",
                  sources=["pyteomics/cmass.c"])
    ])


setup(
    name="pyteomics.cythonize",
    description="An Cython-accelerated version of common pyteomics functions",
    long_description=open("README.rst").read(),
    version="0.2.8",
    packages=['pyteomics'],
    zip_safe=False,
    package_dir={"pyteomics": '.'},
    install_requires=["pyteomics"],
    include_package_data=True,
    ext_modules=extensions,
    maintainer="Joshua Klein",
    maintainer_email="jaklein@bu.edu",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Topic :: Education",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries",
    ],
    license="License :: OSI Approved :: Apache Software License",
    url="https://github.com/mobiusklein/pyteomics.cythonize",
)
