""" Custom multivariate polynomial interpolation library

This is not yet a proper library, so a proper description will be created later

Author: Lorenzo Rota.
"""


class AbstractPolynomialEncoder:
    """Abstract class defining the main functionality of a reconstructible code"""

    def __init__(self, query_complexity, length, alphabet):
        self.query_complexity = query_complexity
        self.length = length
        self.alphabet = alphabet

    def encode_sequence(self, message_seq):
        raise NotImplementedError

    def compute_codeword(self):
        raise NotImplementedError
