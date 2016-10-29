import os


def get_include():
    """Retrieve the path to compiled C extensions' source files to make linking simple.

    This module contains two variants of the algorithm reimplimented using C and the Python-C API.
    """
    return os.path.dirname(__file__)
