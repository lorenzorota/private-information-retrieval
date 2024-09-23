import pickle
from actors.dbms import DBMS


MESSAGE_LOCATION = './messages'


if __name__ == 'examples.encode_multiplicity':
    # specify parameters
    # eg, m=3,d=17,q=256 can encode messages of 1140 bytes  
    m = 3
    d = 11
    q = 16
    s = 2
    
    with open(MESSAGE_LOCATION+'/msg-m{}d{}q{}s{}.txt'.format(m,d,q,s), 'r') as file:
        msg = file.read()

    dbms = DBMS(mode='Multiplicity', params=(m,d,q,s))
    dbms.init_encoder()
    dbms.write_message(msg)