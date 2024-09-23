""" Custom multivariate polynomial interpolation library

This is not yet a proper library, so a proper description will be created later

Author: Lorenzo Rota.
"""

from sage.all import *
from helpers.math import *
from helpers.misc import *
import logging


class ReedMullerEncoder:
    """A container that handles message encoding for Reed-Muller codes
    through polynomial interpolation.
    
    The multivariate polynomial that encodes a message is determined by the 
    parameters m, d, q, and allows for k = (m+d) choose d elements in q to be 
    encoded. Typically a value of q = 256 would be chosen such that each
    ascii character can be mapped to a field element in the polynomial ring.
    
    Example: Pick m=3,d=17,q=256 to encode a message of 1140 bytes.

    Attributes:
        m (int): The number of indeterminates.
        d (int): The maximum degree of each monomial.
        q (int): The finite field order.
        k (int): The total number of monomials.
        poly_ring: a d-degree m-variate polynomial ring with elements in GF(q).
        variables (tuple): The tuple of variables.
        monomials (vector): A monomial basis.
        indices_dict (dict): A dictionary of multi-indices (vector) that represent
            each variable exponent in the monomial basis.
        interpolating_set (list): The collection of sets of k points (m-tuples)
            that map to the message.
        vandermonde (matrix): The vandermonde matrix associated with
            an interpolating set.
        coefs (vector): The list of coefficients that determine the unique
            interpolating polynomial.
        msg (str): The message to be encoded or which has been retrieved
        polynomial (sage polynomial): The polynomial that stores the message.

    """

    def __init__(self, m, d, q):
        self.m = m
        self.d = d
        self.q = q
        self.k = binomial(m+d,d)
        self.poly_ring, self.variables = GF(q)[','.join('x%s'%i for i in range(m))].objgens()
        self.monomials, self.indices_dict = self.construct_monomial_basis()
        self.interpolating_set = None
        self.vandermonde = None
        self.coefs = None
        self.msg = ''
        self.polynomial = None
        logging.debug('initializing a ReedMullerEncoder')
        logging.debug('parameters:')
        logging.debug('-----------------------------------')
        logging.debug('+ m: {}'.format(str(m)))
        logging.debug('+ d: {}'.format(str(d)))
        logging.debug('+ q: {}'.format(str(q)))
        logging.debug('+ k: {} symbols or {} bytes'.format(str(self.k), str(self.k // int(log(256) / log(q)))))
        logging.debug('-----------------------------------')


    def construct_monomial_basis(self):
        """ Generates a monomial basis.
        
        This function will be useful for constructing polynomials with known coefficients
        or to construct generalized Vandermonde matrices to solve for the coefficients.
        
        Returns:
            monomials (vector): A monomial basis.
            indices_dict (dict): A dictionary that associates multi-indices with the corresponding monomial basis position
            
        """
        
        m, d = self.m, self.d
        monomials = []
        indices_dict = dict()
        counter = 0
        for n in range(d+1):
            for powers in constrained_partitions(n, m, 0, n):
                prod = 1
                idx = []
                for i, j in enumerate(powers):
                    prod = prod * self.variables[i]**j
                    idx.append(j)
                monomials.append(prod)
                idx_vec = vector(idx)
                idx_vec.set_immutable()
                indices_dict[idx_vec] = counter
                counter += 1
        return vector(monomials), indices_dict


    def construct_polynomial(self):
        """ Constructs a polynomial from a monomial basis and a list of coefficients
        
        Returns: A polynomial with coefficients from coefs

        """
        
        return self.coefs * self.monomials


    def encode_sequence(self, message_seq):
        """ Encodes a message as a polynomial given some interpolating set

        The encoding is done by mapping each character from the message to 1 or
        multiple elements from GF(q) and then interpolating the data points consisting
        of the interpolating set and the message.

        """

        q =  self.q
        
        if self.vandermonde is None:
            matr = self.comp_vandermonde()
            self.vandermonde = matr
        else:
            matr = self.vandermonde

        ev = vector(GF(q), message_seq)
        self.coefs = matr.inverse() * ev
        self.polynomial = self.monomials * self.coefs
        logging.debug('constructed an encoding polynomial')


    def compute_codeword(self, hyperplanes_condition=True):
        """ Creates the codeword by evaluating the polynomial at all points.

        Returns:
            hyperplanes (list): A list of hyperplanes (list) which partition the evaluated codeword
                (over all enumerated points) from the polynomial.

        """

        m, q = self.m, self.q

        if hyperplanes_condition:
            hyperplanes = [list() for _ in range(q)]
        else:
            hyperplanes = []

        # create enumeration where GF(q) enumerates like an integer
        sorted_GF = []
        for i in range(q):
            sorted_GF.append(to_GF(q,i))

        fast_poly = fast_callable(self.polynomial, vars=self.variables)
        for p in enum_tuple_ring(m, sorted_GF):
            eval_pt = fast_poly(*p)
            # store the evaluation in the correct hyperplane
            if hyperplanes_condition:
                hyperplanes[int(from_GF(16,p[-1]))].append(eval_pt)
            else:
                hyperplanes.append(eval_pt)

        logging.debug('partitioned codeword into {} hyperplanes'.format(str(q)))
        return hyperplanes