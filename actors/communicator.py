""" Custom multivariate polynomial interpolation library

This is not yet a proper library, so a proper description will be created later

Author: Lorenzo Rota.
"""

from sage.all import *
from helpers.math import *
from helpers.misc import *
import logging

DATABASE_PATH = './database'


class Communicator:
    """
    """

    def __init__(self, mode, m, d, q, s):
        self.mode = mode
        self.m = m
        self.d = d
        self.q = q
        self.s = s
        self.k = binomial(m+d,d)
        self.sigma = binomial((m+s-1),m)
        self.interpolating_set = None


    def resolve_query(self, point):
        """Resolve the server and index associated with some query point

        Args:
            point (m-tuple): The point of evaluation.

        Returns:
            server_no: The server that stores the query evaluation.
            indices: The indices at which the evaluation over all partial derivatives is stored.

        """

        m, d, q = self.m, self.d, self.q
        indices = []
        
        data_per_chunk = q**(m-1)
        # determine server and local index
        server_no = from_GF(q, point[-1])
        for i in range(self.sigma):
            indices.append((i * data_per_chunk) + (GF_tuple_to_decimal(q, m, point) // q))

        return server_no, indices


    def carry_out_query(self, server_no, indices, block_number):
            """Carries out the query by reading from the respective database and returning the evaluation

            Args:
                server_no (int): The server that stores part of the database.
                indices (int): The index at which the evaluation is stored.
                block_number (int): The location that stores a particular encoding of the db.
            
            Returns:
                eval (GF(q)): The result of the query. 

            """

            # --- SERVER SHOULD DO THIS ----
            # at the moment this is still strongly coupled to the system
            # in a truly distributed system, this should not have knowledge of the encoding format

            q = self.q
            div_count = int(log(256) / log(q)) # represents the number of segments a byte can be divided into given the size q
            
            # store result for each partial derivative
            symbols = []
            for i in range(self.sigma):
                with open(DATABASE_PATH+'/server{}/{}.bin'.format(server_no, block_number), 'rb') as file:
                    file.seek(indices[i] // div_count, 0)
                    read_byte = int.from_bytes(file.read(1), "big")
                    read_byte = GF_tuple_from_decimal(q, div_count, read_byte)
                # if index corresponds to leading symbol
                if indices[i] % 2 == 0:
                    symbols.append(read_byte[0])
                else:
                    symbols.append(read_byte[1])

                logging.debug('server {} - block {} \t returned: {}'.format(str(server_no), str(block_number), str(symbols[-1])))
            if self.mode == 'ReedMuller':
                return symbols[0]
            else:
                return symbols


    def eval_query_set(self, query_set):
        """Evaluates the query sets
        """

        m, d, q = self.m, self.d, self.q

        evaluations = []
        for point in query_set:
            server_no, indices = self.resolve_query(point)
            logging.debug('resolved query information {{server: {}, indices: {}}}'.format(str(server_no), str(indices)))
            symbol = self.carry_out_query(server_no, indices, 0)
            evaluations.append(symbol)
        logging.debug('done querying')
        logging.debug('------------------------------------------------')

        return evaluations


    def reconstruct_multiple(self, start=None, end=None):
        """Reconstructs a contiguous subsequence of the encoded message

        Args:
            start (int): The starting point of the subsequence.
            end (int): The ending point of the subsequence.

        Returns:
            reconstruction (list): List of reconstructed points of the subsequence.

        """
        
        q = self.q
        sigma = self.sigma
        div_count = int(log(256) / log(q)) # represents the number of times 8 bits can be divided given q

        subsequence_points = []
        if start is None or end is None:
            subsequence_points = self.interpolating_set
        else:
            diff = start - end
            start_is, end_is = int(start * (div_count / sigma)), int(end * (div_count / sigma))
            subsequence_points = self.interpolating_set[start_is:end_is]

        reconstruction = []
        for point in subsequence_points:
            evaled_query_sets, rand_dirs = self.build_eval_query_matrix(point)
            symbols = self.reconstruct(point, evaled_query_sets, rand_dirs)
            for symbol in symbols:
                reconstruction.append(symbol)

        return reconstruction