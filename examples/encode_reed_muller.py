import pickle
from actors.dbms import DBMS


MESSAGE_LOCATION = './messages'


if __name__ == 'examples.encode_reed_muller':
    # specify parameters
    # eg, m=3,d=17,q=256 can encode messages of 1140 bytes  
    m = 3
    d = 13
    q = 16
    
    with open(MESSAGE_LOCATION+'/msg-m{}d{}q{}s{}.txt'.format(m,d,q,1), 'r') as file:
        msg = file.read()

    dbms = DBMS(mode='ReedMuller', params=(m,d,q))
    dbms.init_encoder()
    dbms.write_message(msg)