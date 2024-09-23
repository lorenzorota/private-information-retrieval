""" Custom multivariate polynomial interpolation library

This is not yet a proper library, so a proper description will be created later

Author: Lorenzo Rota.
"""

from sage.all import *
from helpers.math import *
from helpers.misc import *
from helpers.interpolating_set import InterpolatingSet
from encoders.reed_muller import ReedMullerEncoder
from encoders.multiplicity import MultiplicityEncoder
import pickle
import logging


KEYS_PATH = './config/keys'
VANDERMONDES_PATH = './config/vandermondes'
DATABASE_PATH = './database'
MESSAGE_LOCATION = './messages'


class DBMS:
    """
    """

    def __init__(self, mode, params=(3,13,16)):
        self.mode = mode
        if mode == 'ReedMuller':
            self.m = params[0]
            self.d = params[1]
            self.q = params[2]
            self.encoder = ReedMullerEncoder(self.m, self.d, self.q)
        elif mode == 'Multiplicity':
            self.m = params[0]
            self.d = params[1]
            self.q = params[2]
            self.s = params[3]
            self.encoder = MultiplicityEncoder(self.m, self.d, self.q, self.s)
        else:
            raise ValueError('Invalid mode. The available modes are \'ReedMuller\' and \'Multiplicity\'')

        self.msg = ''


    def init_encoder(self, key_path=KEYS_PATH, vm_path=VANDERMONDES_PATH):
        """
        """

        m, d, q = self.m, self.d, self.q
        with open(key_path+'/key_m{}d{}q{}.txt'.format(str(m),str(d),str(q)), 'r') as file:
            key = file.read()

        with open(vm_path+'/vandermonde_m{}d{}q{}.dat'.format(str(m),str(d),str(q)), 'rb') as file:
            vandermonde = pickle.load(file)

        interpolator = InterpolatingSet(self.encoder)
        interpolating_set = interpolator.build_interpolating_sets_from_keys([key])
        self.encoder.interpolating_set = interpolating_set[0]
        self.encoder.vandermonde = vandermonde
        logging.debug('initialized the encoder')


    def translate_message(self, msg):
        """ Constructs a sequence of elements in GF(q) from the message string

        """

        q = self.q
        
        if self.mode == 'ReedMuller':
            k = self.encoder.k
        elif self.mode == 'Multiplicity':
            k = self.encoder.k_out

        div_count = int(log(256) / log(q)) # represents the number of segments a byte can be divided into given the size q

        if k % div_count or k // div_count != len(msg):
            raise ValueError('Message length should be {} characters long'.format(k // div_count))
        elif 8 % div_count:
            raise ValueError('Discrete logarithm of field size must divide 8')

        message_seq = []

        # convert byte to a sequence of field elements and append to message sequence
        for i in range(0, k // div_count):
            byte = ord(msg[i])
            sub_seq = GF_tuple_from_decimal(q, div_count, byte)
            message_seq += sub_seq
        
        if self.encoder.vandermonde is None:
            matr = self.comp_vandermonde()
            self.encoder.vandermonde = matr
        else:
            matr = self.encoder.vandermonde

        logging.debug('translated message into a sequence of {} symbols'.format(k))
        return message_seq


    def build_distributed_db(self, storages, block_number):
        """Distribute the partitioned codeword among all servers and store them in binary files.

        Args:
            storages (list): A list of codeword partitions (list) that will be distributed across all q servers
            block_number (int): The number that will index the current codeword block
        """

        q = self.q
        div_count = int(log(256) / log(q)) # represents the number of segments a byte can be divided into given the size q

        # write the data to the individual files
        for server_no, storage in enumerate(storages):
            with open('database/server{}/{}.bin'.format(server_no, block_number), 'wb') as file:
                for i in range(len(storage) // div_count):
                    # convert each subsequence of div_count many elements to a byte (in decimal representation)
                    # then convert it to a bytes object
                    sub_seq = storage[i*div_count:(i+1)*div_count]
                    byte = GF_tuple_to_decimal(q, div_count, sub_seq)
                    file.write(bytes([byte]))
                logging.debug('updated server {} - block {}'.format(str(server_no), str(block_number)))


    def write_message(self, msg):
        """
        """
        sequence = self.translate_message(msg)
        self.encoder.encode_sequence(sequence)
        hyperplanes = self.encoder.compute_codeword()
        self.build_distributed_db(hyperplanes, 0)
        print('success')