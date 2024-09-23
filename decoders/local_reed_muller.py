""" Custom multivariate polynomial interpolation library

This is not yet a proper library, so a proper description will be created later

Author: Lorenzo Rota.
"""

from sage.all import *
from helpers.math import *
from helpers.misc import *
from decoders.abstract_reconstructible_code import AbstractReconstructibleCode
import random
import logging


class ReedMullerLocalDecoder(AbstractReconstructibleCode):
    """A container that handles local decoding for Reed-Muller codes
    through polynomial interpolation.
    
    The multivariate polynomial that encodes a message is determined by the 
    parameters m, d, q, and allows for k = (m+d) choose d elements in q to be 
    encoded. Typically a value of q = 256 would be chosen such that each
    ascii character can be mapped to a field element in the polynomial ring.
    
    Example: Pick m=3,d=17,q=256 to encode a message of 1140 bytes.

    A local decoder

    Attributes:
        m (int): The number of indeterminates.
        d (int): The maximum degree of each monomial.
        q (int): The finite field order.
        k (int): The total number of monomials.
        query_complexity (int): The number of queries used in the local reconstruction.
        length (int): Length of the codeword.
        alphabet (list): The list of symbols from GF(q).
        interpolating_set (list): The collection of sets of k points (m-tuples)
            that map to the message.
        reconstruction_set (list): The collection of reconstructed symbols from the interpolating set.
        msg (str): The message to be encoded or which has been retrieved
        polynomial (sage polynomial): The polynomial that stores the message.

    """

    def __init__(self, m, d, q):
        self.m = m
        self.d = d
        self.q = q
        self.k = binomial(m+d,d)
        self.query_complexity = q-1
        self.length = q**m
        self.alphabet = list()
        self.interpolating_set = None
        self.reconstruction_set = None
        self.msg = ''
        # construct alphabet
        for elem in GF(q):
            self.alphabet.append(elem)
        logging.debug('initializing a ReedMullerLocalDecoder')
        logging.debug('parameters:')
        logging.debug('-----------------------------------')
        logging.debug('+ m: {}'.format(str(m)))
        logging.debug('+ d: {}'.format(str(d)))
        logging.debug('+ q: {}'.format(str(q)))
        logging.debug('+ k: {} symbols or {} bytes'.format(str(self.k), str(self.k // int(log(256) / log(q)))))
        logging.debug('-----------------------------------')


    def build_query_set(self, position):
        """Builds a query set, which in the case of Reed-Muller is a line to which
        the polynomial is restricted.

        """
        
        m, d, q = self.m, self.d, self.q

        position_vec = vector(GF(q), position)

        logging.debug('construct query set for point: {}'.format(str(position_vec)))

        for attempt in range(100):
            rand_dir = vector(GF(q), GF_tuple_from_decimal(q, m, random.randint(1, self.length)))
            # check that rand_dir is not embedded in any hyperplane
            if rand_dir[-1] != to_GF(q, 0):
                break
        if attempt == 99:
            raise Exception('Unable to find appropriate query set')

        query_set = [position_vec + (t * rand_dir) for t in GF(q)][1:]

        logging.debug('constructed query set: {}'.format(str(query_set)))
        return query_set


    def reconstruct(self, position, symbols):
        """Locally reconstructs a symbol given the evaluations from a query set
        """

        m, d, q = self.m, self.d, self.q

        evals = vector(GF(q), symbols)
        C = codes.ReedSolomonCode(GF(q), Integer(q-1), Integer(d+1))
        decoder = C.decoder("BerlekampWelch")
        recon_poly = decoder.decode_to_message(evals)
        recon_symbol = recon_poly(0)

        logging.debug('reconstructed symbol: {}'.format(str(recon_symbol)))
        return recon_symbol