""" Custom multivariate polynomial interpolation library

This is not yet a proper library, so a proper description will be created later

Author: Lorenzo Rota.
"""

from sage.all import *
from helpers.math import *
from helpers.misc import *
from decoders.local_reed_muller import ReedMullerLocalDecoder
import random
import logging


class MultiplicityLocalDecoder(ReedMullerLocalDecoder):
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
        query_complexity (int): The number of queries used in the local reconstruction.
        length (int): Length of the codeword.
        alphabet (list): The list of symbols from GF(q).
        interpolating_set (list): The collection of sets of k points (m-tuples)
            that map to the message associated to a specific Reed-Muller polynomial.
        msg (str): The message to be encoded or which has been retrieved.
        polynomial (sage polynomial): The polynomial that stores the message.

    """

    def __init__(self, m, d, q, s):
        super().__init__(m, d, q)
        self.s = s
        self.sigma = binomial((m+s-1),m)
        self.k_out = self.k * self.sigma
        logging.debug('initializing a MultiplicityLocalDecoder')
        logging.debug('parameters:')
        logging.debug('-----------------------------------')
        logging.debug('+ m: {}'.format(str(m)))
        logging.debug('+ d: {}'.format(str(d)))
        logging.debug('+ q: {}'.format(str(q)))
        logging.debug('+ s: {}'.format(str(s)))
        logging.debug('+ k_out: {} symbols or {} bytes'.format(str(self.k_out), str(self.k_out // int(log(256) / log(q)))))
        logging.debug('-----------------------------------')

    def build_query_set(self, position, previous_dirs):
        """Builds a query set, which in the case of Multiplicity codes is a line to which
        the polynomial is restricted.
        
        Args:
            previous_dirs (list): List of other random directions (to check for collinearity)
        
        Returns:
            query_set: The query_set
            rand_dir: The random direction that determined the query_set
        
        """
        m, d, q = self.m, self.d, self.q

        position_vec = vector(GF(q), position)
        logging.debug('construct query set for point: {}'.format(str(position_vec)))

        for attempt in range(100):
            not_collinear = True
            rand_dir = vector(GF(q), GF_tuple_from_decimal(q, m, random.randint(1, self.length)))
            # check that rand_dir is not colinear with the other lines
            for other_dir in previous_dirs:
                try:
                    rand_dir / other_dir
                    not_collinear = False
                except ArithmeticError:
                    not_collinear = True
            # check that rand_dir is not embedded in any hyperplane and not
            if rand_dir[-1] != to_GF(q, 0) and not_collinear:
                break
        if attempt == 99:
            raise Exception('Unable to find appropriate query set')

        query_set = [position_vec + (t * rand_dir) for t in GF(q)][1:]
        logging.debug('constructed query set: {}'.format(str(query_set)))
        return query_set, rand_dir


    def build_query_sets(self, position):
        """Construct a query set for each partial derivative, and evaluate each
        query set (a list of matrices for each line).

        Args:
            position: The position (m-tuple) at which we want to query.

        Returns:
            evaled_query_sets (list): The list of all evaluated query sets.
            rand_dirs (list): The list of random directions for each query set.

        """

        query_sets = []
        rand_dirs = []

        # create a query set for each partial derivative
        logging.debug('construct query sets for all partial derivatives')
        logging.debug('------------------------------------------------')
        for i in range(self.sigma):
            query_set, rand_dir = self.build_query_set(position, rand_dirs)
            query_sets.append(query_set)
            rand_dirs.append(rand_dir)

        return query_sets, rand_dirs


    def univariate_evaluated_hasse(self, evaluated_queries, rand_dir, b, e):
        """ Compute the derivative of the polynomial restricted to a random line.
        
        Args:
            evaluated_queries: The evaluated query set with points in GF(q)^sigma
            rand_dir (vector): The random direction of the line
            b (int): The position in the query set of the element to be evaluated
            e (int): The derivative order that the restricted polynomial needs to be evaluated at
        
        Returns:
            solution: Element in GF(q) that is the evaluation of the b-th point in the query set
            
        """
        
        m, s = self.m, self.s
        
        # create multi_index lookup table, which will be needed for retrieving specific derivatives
        multi_index_lookup = dict()
        for pos, v in enumerate(constrained_bounded_partitions(s, m)):
            v = vector(v)
            v.set_immutable()
            multi_index_lookup[v] = pos
            
        solution = 0
        for v in constrained_partitions(e, m, 0, e):
            v = vector(v)
            v.set_immutable()
            v_pos = multi_index_lookup[v]
            # print the v-th coordinate of the b-th query from the 0th line
            # Compute U_i^v
            u_power = 1
            for j in range(m):
                u_power *= rand_dir[j]**v[j]
            solution += evaluated_queries[v_pos][b] * u_power
        return solution


    def fetch_restricted_poly_derivative_table(self, evaluated_queries, rand_dir):
        """ Computes the restricted polynomial derivative table using the univariate derivative eval function
        
        Args:
            evaluated_queries: The evaluated query set with points in GF(q)^sigma
            rand_dir (vector): The random direction of the line
            
        Returns:
            deriv_table (matrix): Columns represent [x f(x) f'(x) ...]
            
        """
        
        q, s = self.q, self.s

        deriv_table = []
        query_arg = list(GF(q))[1:]
        for b, num in enumerate(query_arg):
            evals_deriv = [num]
            for e in range(s):
                evals_deriv.append(self.univariate_evaluated_hasse(evaluated_queries, rand_dir, b, e))
            deriv_table.append(evals_deriv)
        return matrix(deriv_table)


    def multiplicity_berlekamp_welch(self, deriv_table):
        """ The modified Berlekamo-Welch algorithm for Multiplicity Codes.
        [A description will be given...]

        Args:
            deriv_table (matrix): Columns represent [x f(x) f'(x) ...]

        Returns:
            Q (polynomial): The error-corrected polynomial that suports

        """

        def gen_equation(deriv_table, point_idx, deriv_order):
            """ Compute the polynomials D^(i)N(p) - D^(i)(E*l)(p) with unkown coefficients and zero on the RHS

            Args:
                deriv_table (matrix): Columns represent [x f(x) f'(x) ...]
                point_idx (int): Row index that specifies a particular point x
                deriv_order (int): The derivative order up to which we differentiate

            Returns:
                LHS (polynomial): Left-hand side of the equation (since we want to find the kernel)
                
            """
            
            # separate variables
            a = a_and_b[0:deg_E+1]
            b = a_and_b[deg_E+1:deg_E+deg_N+2]

            # fetch point from table (0th column)
            point = deriv_table[point_idx][0]
            point_vec = vector(GF(q), [point])

            # set initial term to 0 (but keep it as a polynomial object)
            term = a[0]*0

            # extract derivatives from deriv_table
            derivatives = deriv_table[:,1:]

            # Leibniz product rule
            for j in range(deriv_order + 1):
                deriv_order_vec = vector([j])
                E_deriv = hasse_derivative(E_monomials, E_indices_dict, a, deriv_order_vec, point_vec, 1)
                term += E_deriv * derivatives[point_idx][deriv_order - (j)]
            N_deriv = hasse_derivative(N_monomials, N_indices_dict, b, vector([deriv_order]), point_vec, 1)

            LHS = N_deriv - (term)
            return LHS

        d, q, s = self.d, self.q, deriv_table.ncols()-1

        # initalize variable for univariate polynomial
        x = GF(q)['x'].objgens()[1][0]

        # compute d_out and then deg_E and deg_N accordingly
        d_out = (s-1)*(q-1) + d
        deg_E = floor((s*(q-1) - d_out) / 2)
        deg_N = ceil((s*(q-1) + d_out) / 2)
        E_monomials = vector([x**i for i in range(deg_E + 1)])
        N_monomials = vector([x**i for i in range(deg_N + 1)])

        # create index dictionaries (messy but hey, that's how I implemented the Hasse derivative)
        E_indices_dict = dict()
        for i in range(deg_E + 1):
            vec = vector([i])
            vec.set_immutable()
            E_indices_dict[vec] = i

        N_indices_dict = dict()
        for i in range(deg_N + 1):
            vec = vector([i])
            vec.set_immutable()
            N_indices_dict[vec] = i

        a_and_b = vector(GF(q)[','.join('a%s'%i for i in range(deg_E + 1))+','+(','.join('b%s'%i for i in range(deg_N + 1)))].objgens()[1])

        # compute the polynomials D^(i)N(x) - D^(i)(E*l)(x) with zeros on the RHS
        eqs = []
        for i in range(q-1):
            for j in range(s):
                eq = gen_equation(deriv_table, i, j)
                eqs.append(eq)

        # create the evaluated monomials matrix
        matr = []
        for eq in eqs:
            row = []
            for var in a_and_b:
                elem = (eq.coefficient(var))
                row.append(GF(q)(elem))
            matr.append(row)
        matr = matrix(matr)
        
        # compute the kernel and find a working solution
        kernel = matr.right_kernel()
        
        for solution in kernel.basis():
            coef_vec = vector(GF(q), solution)
            
            # from the valid coefficients, compute the error-locator polynomial
            E = coef_vec[0:deg_E+1] * E_monomials
            N = coef_vec[deg_E+1:deg_E+deg_N+2] * N_monomials
            
            if E != 0:
                Q = N/E
                # The polynomial needs to actually be in the quotient ring
                if Q.denominator().is_unit():
                    return GF(q)[x](Q)
        return 0


    def reconstruct(self, position, evaled_query_sets, rand_dirs):
        """Locally reconstructs a symbol given the evaluations from a query set and their random directions
        """

        m, d, q, s = self.m, self.d, self.q, self.s

        # compute the derivative tables for each restricted polynomial (dependent on rand_dir)
        all_deriv_tables = []
        for i in range(self.sigma):
            all_deriv_tables.append(self.fetch_restricted_poly_derivative_table(evaled_query_sets[i], rand_dirs[i]))
            logging.debug('constructed univariate polynomial derivative table')

        # apparently we have to consider polynomial interpolations with varying derivative orders
        # this is because if we include evaluations of derivative order >2, then polynomial evaluated at 0 suddenly different
        restricted_poly = []
        for i in range(s):
            row = []
            for j in range(self.sigma):
                result = self.multiplicity_berlekamp_welch(all_deriv_tables[j][:,:2+i])
                if result == None:
                    raise ValueError('Berlekamp-Welch reconstruction failed.')
                row.append(result)
            restricted_poly.append(row)
            logging.debug('recovered univariate polynomial with modified berlekamp-welch')

        # compute the coefficients of the i-th derivative of all restricted polynomials
        # in total there are (sigma * s) coefficients
        deriv_coeffs = []
        for i in range(s):
            for j in range(self.sigma):
                poly = restricted_poly[i][j]
                poly_indices_dict = comp_indices_dict(poly, 1)
                deriv_order = vector([i])
                deriv = hasse_derivative(poly.monomials(), poly_indices_dict, poly.coefficients(), deriv_order, vector([0]), 1)
                deriv_coeffs.append(deriv)
        deriv_coeffs = vector(deriv_coeffs)

        # list of vector of multi-indices, separated by their weight (upto s)
        derivatives_of_weight = []
        for e in range(s):
            derivatives_of_weight.append(list(constrained_partitions(e,m,0,m)))

        # matrix of directions exponentiated by multi-indices
        direction_matr = [[None for _ in range(self.sigma)] for _ in range(self.sigma * s)]
        for i in range(self.sigma*s):
            
            e = i // self.sigma
            for j, v in enumerate(constrained_bounded_partitions(s,m)):
                if v in derivatives_of_weight[e]:
                    v = vector(v)
                    v.set_immutable()
                    u_power = 1
                    for k in range(m):
                        u_power *= rand_dirs[i % self.sigma][k] ** v[k]
                else:
                    u_power = 0
                direction_matr[i][j] = u_power
        direction_matr = matrix(direction_matr)

        # reconstruction of overdetermined system (should have exactly 1 solution due to intersection of lines)
        recon_symbols = direction_matr.solve_right(deriv_coeffs)

        logging.debug('reconstructed the symbols: {}'.format(str(recon_symbols)))
        return recon_symbols