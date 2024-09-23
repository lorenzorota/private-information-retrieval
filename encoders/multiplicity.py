""" Custom multivariate polynomial interpolation library

This is not yet a proper library, so a proper description will be created later

Author: Lorenzo Rota.
"""

from sage.all import *
from helpers.math import *
from helpers.misc import *
from encoders.reed_muller import ReedMullerEncoder
import logging


class MultiplicityEncoder(ReedMullerEncoder):
    """A container that handles message encoding for Multiplicity codes
    through polynomial interpolation for individual Reed-Muller polynomials,
    and then constructing a outer polynomial that can be locally decoded for all
    of its partial derivatives to obtain parts of the message.
    
    The multivariate polynomial that encodes part of a message is determined by the 
    parameters m, d, q, and allows for k = (m+d) choose d elements in q to be 
    encoded. Typically a value of q = 256 would be chosen such that each
    ascii character can be mapped to a field element in the polynomial ring.
    All such polynomials can be obtained in the encoding process in a way that when
    they are used to construct the outer polynomial (of degree d_out), its derivatives
    restricted to each interpolating set yield the individual polynomial restricted to
    its interpolating set. This in turn allows for encoding messages of size k*sigma
    where sigma = (m+s-1) choose m
    
    Example: Pick m=3,d=17,q=256,s=2 to encode a message of 3420 bytes.

    Attributes:
        m (int): The number of indeterminates.
        d (int): The maximum degree of each monomial in the Reed-Muller polynomials.
        d_out (int): The maximum degree of the Multiplicity polynomial
        q (int): The finite field order.
        s (int): The multiplicity of the code. 
        k (int): The total number of monomials for each Reed-Muller polynomial
        k_out (int): The total number of monomials added up from each Reed-Muller polynomial
        sigma (int): The total number of partial derivatives upto the multiplicity.
        poly_ring: a d-degree m-variate polynomial ring with elements in GF(q).
        variables (tuple): The tuple of variables.
        monomials (vector): A monomial basis for each Reed-Muller polynomial.
        indices_dict (dict): A dictionary of multi-indices (vector) that represent
            each variable exponent in the monomial basis of any Reed-Muller polynomial.
        interpolating_sets (list): The collection of sets of k points (m-tuples)
            that map to the message associated to a specific Reed-Muller polynomial.
        vandermondes (list): The collection of vandermonde matrices (matrix)
            associated with an interpolating set of a Reed-Muller polynomial.
        coefs (list): The list of coefficients (vector) that determine every unique
            interpolating Reed-Muller polynomial.
        msg (str): The message to be encoded or which has been retrieved.
        polynomial (sage polynomial): The polynomial that stores the message.

    """

    def __init__(self, m, d, q, s):
        super().__init__(m, d, q)
        self.s = s
        self.sigma = binomial((m+s-1),m)
        #self.interpolating_sets = []
        #self.vandermondes = []
        self.k_out = self.k * self.sigma
        logging.debug('initializing a MultiplicityEncoder')
        logging.debug('parameters:')
        logging.debug('-----------------------------------')
        logging.debug('+ m: {}'.format(str(m)))
        logging.debug('+ d: {}'.format(str(d)))
        logging.debug('+ q: {}'.format(str(q)))
        logging.debug('+ s: {}'.format(str(s)))
        logging.debug('+ k_out: {} symbols or {} bytes'.format(str(self.k_out), str(self.k_out // int(log(256) / log(q)))))
        logging.debug('-----------------------------------')

    def encode_sequence(self, message_seq):
        """ Encodes a message as a polynomial given some interpolating set

        The encoding is done by mapping each character from the message to 1 or
        multiple elements from GF(q) and then interpolating the data points consisting
        of the interpolating set and the message.

        """

        m, d, q, s, k = self.m, self.d, self.q, self.s, self.k
        sigma, k_out = self.sigma, self.k_out

        # split up message into sigma many parts by incrementally modding out.
        # this way we ensure that no contigous subsequence of the message can be
        # determined over a single interpolating set.


        message_segments = [list() for _ in range(sigma)]
        for i in range(k_out):
            message_segments[i % sigma].append(message_seq[i])
        
        for i in range(sigma):
            message_segments[i] = vector(message_segments[i])
    
        polynomials = dict()

        def V_k(variables, k):
            """ Sub-function, which is only relevant in the construction
            of the outer polynomial.

            Explanation:
                The outer polynomial is a sum of the terms V_k * A_k for
                as many terms as there are partial derivatives upto some order.
                V_k is a product of (X^q - X), for all indeterminates X
                such that the whole term vanishes when when replaced with a point.

            Returns:
                prod (polynomial): The product of of vanishing terms
            """

            prod = 1
            for j in range(m):
                prod *= (variables[j]**q - variables[j])**k[j]
            return prod
            
        partials = []
        for i in constrained_bounded_partitions(s, m):
            partial = vector(i)
            partial.set_immutable()
            partials.append(partial)

        # iterate over all partial derivatives
        for pos, i in enumerate(partials):
            # base case:
            # initialize first polynomial to be zero s.th outer polynomial is the first inner polynomial
            if pos == 0:
                polynomials[i] = 0

            terms = []
            # iterate over all remaining unknown inner polynomials
            for k in all_multi_index_leq(partials, i, m):
                k = vector(k)
                k.set_immutable()
                if k != i:
                    # Leibniz' product rule for Hasse derivatives
                    for j in all_multi_index_bdd(partials,k,i,m):
                        Vk_term = V_k(self.variables, k)
                        Vk_indices_dict = comp_indices_dict(Vk_term, m)
                        
                        # compute derivative of V_k
                        Vk_derivative = hasse_derivative(Vk_term.monomials(), Vk_indices_dict, Vk_term.coefficients(), j, self.variables, m)
                        
                        Ak = polynomials[k]
                        Ak_indices_dict = comp_indices_dict(Ak, m)
                        
                        # compute derivative of A_k
                        Ak_derivative = hasse_derivative(Ak.monomials(), Ak_indices_dict, Ak.coefficients(), i-j, self.variables, m)
                        
                        # compute derivative term
                        terms.append(Vk_derivative * Ak_derivative)

            deriv_terms = sum(terms)

            if type(deriv_terms) is not sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular:
                eval_deriv_terms = deriv_terms * vector([1 for _ in self.interpolating_set])
            else:
                fast_deriv_terms = fast_callable(deriv_terms, vars=self.variables)
                eval_deriv_terms = vector([fast_deriv_terms(*p) for p in self.interpolating_set])
            polynomial_coefs = self.vandermonde.inverse()*(eval_deriv_terms + (-1)**sum(i) * message_segments[pos])
            polynomials[i] = polynomial_coefs * self.monomials
            logging.debug('constructed an inner encoding polynomial for partial derivative: {}'.format(str(i)))

        # construct the outer polynomial
        Q = 0
        for k in partials:
            Q += polynomials[k]*V_k(self.variables, k)
        
        self.coefs = Q.coefficients()
        self.polynomial = Q
        logging.debug('constructed an outer encoding polynomial')


    def compute_codeword(self, hyperplanes_condition=True):
        """ Creates the codeword by evaluating the outer polynomial and its partial derivatives
        at all points.

        Returns:
            hyperplanes (list): A list of hyperplanes (list) which partition the evaluated codeword
                (over all enumerated points) from the outer polynomial and all of its partial derivatives.

        """

        m, q, s = self.m, self.q, self.s

        if hyperplanes_condition:
            hyperplanes = [list() for _ in range(q)]
        else:
            hyperplanes = []
            
        poly_indices_dict = comp_indices_dict(self.polynomial, m)

        # create enumeration where GF(q) enumerates like an integer
        sorted_GF = []
        for i in range(q):
            sorted_GF.append(to_GF(q,i))

        partials = []
        for i in constrained_bounded_partitions(s, m):
            partial = vector(i)
            partial.set_immutable()
            partials.append(partial)

        for idx, i in enumerate(partials):
            deriv = hasse_derivative(self.polynomial.monomials(), poly_indices_dict, self.polynomial.coefficients(), i, self.variables, m)
            fast_deriv = fast_callable(deriv, vars=self.variables)
            for p in enum_tuple_ring(m, sorted_GF):
                eval_pt = fast_deriv(*p)
                # store the evaluation in the correct hyperplane
                if hyperplanes_condition:
                    hyperplanes[int(from_GF(16,p[-1]))].append(eval_pt)
                else:
                    hyperplanes.append(eval_pt)

        logging.debug('partitioned codeword into {} hyperplanes'.format(str(q)))
        return hyperplanes