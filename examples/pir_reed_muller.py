import pickle
from actors.user import User


KEYS_PATH = './config/keys'


if __name__ == 'examples.pir_reed_muller':
    # set parameters
    m = 3
    d = 13
    q = 16

    user = User(mode='ReedMuller', params=(m,d,q))
    user.init_decoder()
    sequence = user.query_sequence(0,20)
    msg = user.translate_message(sequence)
    print(msg)