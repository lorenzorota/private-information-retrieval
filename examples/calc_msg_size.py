from sage.all import *


if __name__ == 'examples.calc_msg_size':
    # specify parameters
    # eg, m=3,d=17,q=256 can encode messages of 1140 bytes
    m = 2
    d = 4
    q = 16
    s = 3
    
    k = binomial(m+d, d)
    sigma = binomial(m+s-1, m)
    k_out = sigma*k
    # represents the number of segments a byte can be divided into given the size q
    div_count = int(log(256) / log(q))

    print('Message length: {} bytes'.format(k_out // div_count))
    print('Sequence length: {} symbols'.format(k_out))
    print('Sigma: {}'.format(sigma))
