import os
import csv
from actors.user import User
from actors.dbms import DBMS
from sage.all import *


# PIR parameters
MESSAGE_LOCATION = './messages'
KEYS_PATH = './config/keys'

# Benchmarking config
MSG_QUERY_LEN = 20
TRIALS = 5


def multiplicity_encode(m,d,q,s):
    with open(MESSAGE_LOCATION+'/msg-m{}d{}q{}s{}.txt'.format(m,d,q,s), 'r') as file:
        msg = file.read()

    dbms = DBMS(mode='Multiplicity', params=(m,d,q,s))
    dbms.init_encoder()
    dbms.write_message(msg)


def multiplicity_pir(m,d,q,s):
    user = User(mode='Multiplicity', params=(m,d,q,s))
    user.init_decoder()
    sequence = user.query_sequence(0,MSG_QUERY_LEN)
    msg = user.translate_message(sequence)


def reed_muller_encode(m,d,q):
    with open(MESSAGE_LOCATION+'/msg-m{}d{}q{}s{}.txt'.format(m,d,q,1), 'r') as file:
        msg = file.read()

    dbms = DBMS(mode='ReedMuller', params=(m,d,q))
    dbms.init_encoder()
    dbms.write_message(msg)


def reed_muller_pir(m,d,q):
    user = User(mode='ReedMuller', params=(m,d,q))
    user.init_decoder()
    sequence = user.query_sequence(0,MSG_QUERY_LEN)
    msg = user.translate_message(sequence)


if __name__ == 'examples.benchmark':

    multiplicity_encoder_file = open('multiplicity-encoder-benchmark.csv', 'w')
    multiplicity_encoder_csv = csv.writer(multiplicity_encoder_file)
    multiplicity_encoder_csv.writerow(['m', 'd', 'q', 's', 'k_out', 'msg_size', 'min_time', 'max_time', 'mean_time'])
    multiplicity_encoder_file.flush()

    multiplicity_decoder_file = open('multiplicity-decoder-benchmark.csv', 'w')
    multiplicity_decoder_csv = csv.writer(multiplicity_decoder_file)
    multiplicity_decoder_csv.writerow(['m', 'd', 'q', 's', 'k_out', 'msg_size', 'min_time', 'max_time', 'mean_time'])
    multiplicity_decoder_file.flush()

    # --- First we test all the Multiplicity codes ----
    print('')
    print('------------------------------------------------------------------------')
    print('Benchmarking the multiplicity code based PIR scheme')
    print('- msg query length: {} characters'.format(MSG_QUERY_LEN))
    print('- number of trials: {} times'.format(TRIALS))
    print('------------------------------------------------------------------------')
    print('')
    q = 16
    for m in range(2,4):
        for d in range(1,15):
            for s in range(1,5):

                k = binomial(m+d, d)
                sigma = binomial(m+s-1, m)
                k_out = sigma*k
                # represents the number of segments a byte can be divided into given the size q
                div_count = int(log(256) / log(q))

                msg_len = k_out // div_count

                with open(MESSAGE_LOCATION+'/msg-m{}d{}q{}s{}.txt'.format(m,d,q,s), 'r') as file:
                    msg = file.read()

                # only performing test on messages larger than query size
                # and messages that are exactly divisible by div_count
                if k_out % div_count or k_out // div_count != len(msg) or msg_len < MSG_QUERY_LEN:
                    pass
                else:
                    print('------------------------------------------------------------------------')
                    print('m:{} d:{} q:{} s:{}'.format(m,d,q,s))
                    try:
                        time_obj = timeit('multiplicity_encode(m,d,q,s)', globals=globals(), number=TRIALS)
                        min_time = min(time_obj.series)
                        max_time = max(time_obj.series)
                        mean_time = mean(time_obj.series)
                        multiplicity_encoder_csv.writerow([str(m), str(d), str(q), str(s), str(k_out), str(msg_len), str(min_time), str(max_time), str(mean_time)])
                        multiplicity_encoder_file.flush()
                        print('Encoding time (best): {}'.format(min_time))
                    except ValueError as err:
                        print('Could not encode: {}'.format(err))
                    except Exception as err:
                        print('Unexpected error')
                        print('Could not encode: {}'.format(err))

                    try:
                        time_obj = timeit('multiplicity_pir(m,d,q,s)', globals=globals(), number=TRIALS)
                        min_time = min(time_obj.series)
                        max_time = max(time_obj.series)
                        mean_time = mean(time_obj.series)
                        multiplicity_decoder_csv.writerow([str(m), str(d), str(q), str(s), str(k_out), str(msg_len), str(min_time), str(max_time), str(mean_time)])
                        multiplicity_decoder_file.flush()
                        print('Decoding time (best): {}'.format(min_time))
                    except ValueError as err:
                        print('Could not decode: {}'.format(err))
                    except Exception:
                        print('Unexpected error')

    multiplicity_encoder_file.close()
    multiplicity_decoder_file.close()

    reed_muller_encoder_file = open('reed-muller-encoder-benchmark.csv', 'w')
    reed_muller_encoder_csv = csv.writer(reed_muller_encoder_file)
    reed_muller_encoder_csv.writerow(['m', 'd', 'q', 'k_out', 'msg_size', 'min_time', 'max_time', 'mean_time'])
    reed_muller_encoder_file.flush()

    reed_muller_decoder_file = open('reed-muller-decoder-benchmark.csv', 'w')
    reed_muller_decoder_csv = csv.writer(reed_muller_decoder_file)
    reed_muller_decoder_csv.writerow(['m', 'd', 'q', 'k_out', 'msg_size', 'min_time', 'max_time', 'mean_time'])
    reed_muller_decoder_file.flush()

    # --- Now we test all the Reed-Muller codes ----
    print('')
    print('------------------------------------------------------------------------')
    print('Benchmarking the Reed-Muller code based PIR scheme')
    print('- msg query length: {} characters'.format(MSG_QUERY_LEN))
    print('- number of trials: {} times'.format(TRIALS))
    print('------------------------------------------------------------------------')
    print('')
    q = 16
    for m in range(2,4):
        for d in range(1,15):

                k = binomial(m+d, d)
                k_out = k
                # represents the number of segments a byte can be divided into given the size q
                div_count = int(log(256) / log(q))

                msg_len = k_out // div_count

                with open(MESSAGE_LOCATION+'/msg-m{}d{}q{}s{}.txt'.format(m,d,q,1), 'r') as file:
                    msg = file.read()

                # only performing test on messages larger than query size
                # and messages that are exactly divisible by div_count
                if k_out % div_count or k_out // div_count != len(msg) or msg_len < MSG_QUERY_LEN:
                    pass
                else:
                    print('------------------------------------------------------------------------')
                    print('m:{} d:{} q:{}'.format(m,d,q))
                    try:
                        time_obj = timeit('reed_muller_encode(m,d,q)', globals=globals(), number=TRIALS)
                        min_time = min(time_obj.series)
                        max_time = max(time_obj.series)
                        mean_time = mean(time_obj.series)
                        reed_muller_encoder_csv.writerow([str(m), str(d), str(q), str(k_out), str(msg_len), str(min_time), str(max_time), str(mean_time)])
                        reed_muller_encoder_file.flush()
                        print('Encoding time (best): {}'.format(min_time))
                    except ValueError as err:
                        print('Could not encode: {}'.format(err))
                    except Exception as err:
                        print('Unexpected error: {}'.format(err))

                    try:
                        time_obj = timeit('reed_muller_pir(m,d,q)', globals=globals(), number=TRIALS)
                        min_time = min(time_obj.series)
                        max_time = max(time_obj.series)
                        mean_time = mean(time_obj.series)
                        reed_muller_decoder_csv.writerow([str(m), str(d), str(q), str(k_out), str(msg_len), str(min_time), str(max_time), str(mean_time)])
                        reed_muller_decoder_file.flush()
                        print('Decoding time (best): {}'.format(min_time))
                    except ValueError as err:
                        print('Could not decode: {}'.format(err))
                    except Exception as err:
                        print('Unexpected error: {}'.format(err))

    reed_muller_encoder_file.close()
    reed_muller_decoder_file.close()    