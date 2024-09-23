""" Miscellaneous helper functions 

This is not yet a proper library, so a proper description will be created later.

Author: Lorenzo Rota
"""

import itertools


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks"""

    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)