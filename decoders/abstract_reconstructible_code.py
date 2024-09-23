""" Custom multivariate polynomial interpolation library

This is not yet a proper library, so a proper description will be created later

Author: Lorenzo Rota.
"""


class AbstractReconstructibleCode:
    """Abstract class defining the main functionality of a reconstructible code"""

    def __init__(self, query_complexity, length,alphabet):
        self.query_complexity = query_complexity
        self.length = length
        self.alphabet = alphabet


    def build_query_set(self,position):
        raise NotImplementedError


    def reconstruct(self, position, symbols):
        raise NotImplementedError