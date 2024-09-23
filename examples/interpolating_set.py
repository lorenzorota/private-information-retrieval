import pickle
import os
from sage.all import *
from encoders.reed_muller import ReedMullerEncoder
from encoders.multiplicity import MultiplicityEncoder
from helpers.interpolating_set import InterpolatingSet


KEYS_PATH = './config/keys'
VANDERMONDES_PATH = './config/vandermondes'


if __name__ == 'examples.interpolating_set':
    # specify parameters
    # eg, m=3,d=17,d=256 can encode messages of 1140 bytes
    m = 2
    d = 1
    q = 16
    s = 1


    # specify whether we work with ReedMuller or Multiplicity codes
    code_type = 'ReedMuller'

    if code_type == 'ReedMuller':
        encoder = ReedMullerEncoder(m, d, q)
        interpolator = InterpolatingSet(encoder)

    elif code_type == 'Multiplicity':
        encoder = MultiplicityEncoder(m, d, q, s)
        interpolator = InterpolatingSet(encoder)

    interpolator.find_set()
    keys = interpolator.build_keys()
    vandermondes = interpolator.vandermondes


    # write the keys into their appropriate folder(s)
    if len(keys) > 1:
        for count, key in enumerate(keys):
            new_path = KEYS_PATH+'/keys_m{}d{}q{}'.format(str(m),str(d),str(q))
            if not os.path.exists(new_path):
                os.mkdir(new_path)
            with open(new_path+'/{}.txt'.format(str(count)), 'w') as file:
                file.write(key)
    else:
        with open(KEYS_PATH+'/key_m{}d{}q{}.txt'.format(str(m),str(d),str(q)), 'w') as file:
            file.write(keys[0])

    # write the vandermonde matrix dumps into their appropriate folder(s)
    if len(vandermondes) > 1:
        for count, vandermonde in enumerate(vandermondes):
            new_path = VANDERMONDES_PATH+'/vandermondes_m{}d{}q{}'.format(str(m),str(d),str(q))
            if not os.path.exists(new_path):
                os.mkdir(new_path)
            with open(new_path+'/{}.dat'.format(str(count)), 'wb') as file:
                pickle.dump(vandermonde, file)
    else:
        with open(VANDERMONDES_PATH+'/vandermonde_m{}d{}q{}.dat'.format(str(m),str(d),str(q)), 'wb') as file:
            pickle.dump(vandermondes[0], file)