import pickle
from actors.user import User


KEYS_PATH = './config/keys'


if __name__ == 'examples.pir_multiplicity':
    # set parameters
    m = 3
    d = 11
    q = 16
    s = 2

    user = User(mode='Multiplicity', params=(m,d,q,s))
    user.init_decoder()
    sequence = user.query_sequence(0,20)
    msg = user.translate_message(sequence)
    print(msg)