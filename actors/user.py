""" Custom multivariate polynomial interpolation library

This is not yet a proper library, so a proper description will be created later

Author: Lorenzo Rota.
"""

from sage.all import *
from helpers.math import *
from helpers.misc import *
from helpers.interpolating_set import InterpolatingSet
from actors.communicator import Communicator
from decoders.local_reed_muller import ReedMullerLocalDecoder
from decoders.local_multiplicity import MultiplicityLocalDecoder
import logging


KEYS_PATH = './config/keys'


class User:
    """
    """

    def __init__(self, mode='ReedMuller', params=(3,13,16)):
        self.mode = mode
        if mode == 'ReedMuller':
            self.m = params[0]
            self.d = params[1]
            self.q = params[2]
            self.decoder = ReedMullerLocalDecoder(self.m, self.d, self.q)
            self.communicator = Communicator(mode, self.m, self.d, self.q, 1)
        elif mode == 'Multiplicity':
            self.m = params[0]
            self.d = params[1]
            self.q = params[2]
            self.s = params[3]
            self.decoder = MultiplicityLocalDecoder(self.m, self.d, self.q, self.s)
            self.communicator = Communicator(mode, self.m, self.d, self.q, self.s)
        else:
            raise ValueError('Invalid mode. The available modes are \'ReedMuller\' and \'Multiplicity\'')


    def init_decoder(self, path=KEYS_PATH):
        """
        """
        
        m, d, q = self.m, self.d, self.q
        with open(path+'/key_m{}d{}q{}.txt'.format(str(m),str(d),str(q)), 'r') as file:
            key = file.read()

        interpolator = InterpolatingSet(self.decoder)
        interpolating_set = interpolator.build_interpolating_sets_from_keys([key])
        self.decoder.interpolating_set = interpolating_set[0]


    def query_sequence(self, start=None, end=None):
        """
        
        """

        q = self.q
        div_count = int(log(256) / log(q)) # represents the number of times 8 bits can be divided given q

        if self.mode == 'ReedMuller':
            sigma = 1
        elif self.mode == 'Multiplicity':
            sigma = self.decoder.sigma

        subsequence_points = []
        if start is None or end is None:
            subsequence_points = self.decoder.interpolating_set
        else:
            diff = start - end
            start_is, end_is = int(start * (div_count / sigma)), int(end * (div_count / sigma))
            subsequence_points = self.decoder.interpolating_set[start_is:end_is]

        # locally reconstruct symbols for each query that is made (slightly different for Multiplicity codes)
        reconstruction = []
        for point in subsequence_points:
            if self.mode == 'ReedMuller':
                query_set = self.decoder.build_query_set(point)
                evaluations = self.communicator.eval_query_set(query_set)
                symbol = self.decoder.reconstruct(point, evaluations)
                reconstruction.append(symbol)
            elif self.mode == 'Multiplicity':
                query_sets, rand_dirs = self.decoder.build_query_sets(point)
                all_evaluations = []
                # evaluate all individual query sets (one for each partial derivative)
                for query_set in query_sets:
                    all_evaluations.append(matrix(self.communicator.eval_query_set(query_set)).transpose())
                symbols = self.decoder.reconstruct(point, all_evaluations, rand_dirs)
                for symbol in symbols:
                    reconstruction.append(symbol)
        logging.debug('------------------------------------------------')
        logging.debug('locally reconstructed all symbols')
        logging.debug('------------------------------------------------')
        return reconstruction


    def translate_message(self, msg):
        """ Retrieves the message from the queried sequence """

        q, k = self.q, self.decoder.k
        div_count = int(log(256) / log(q)) # represents the number of times 8 bits can be divided given q

        # if k % div_count:
        #     raise ValueError('Message length should be {} characters long'.format(k // div_count))
        if 8 % div_count:
            raise ValueError('Discrete logarithm of field size must divide 8')

        msg_in_GF = msg
        msg = ''
        for i in range(len(msg_in_GF) // div_count):
            # convert each subsequence of div_count many elements to a byte (in decimal representation)
            # then convert it to a char for the 'message' representation
            sub_seq = msg_in_GF[i*div_count:(i+1)*div_count]
            byte = GF_tuple_to_decimal(q, div_count, sub_seq)
            msg += chr(byte)
        
        logging.debug('translated sequence of symbols into message')
        return msg