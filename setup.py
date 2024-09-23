import sys, getopt, os
import logging

error_msg = """Usage: ./setup.py [-d <server-count>] [-hV]

-h, -help,          --help                  Display help.

-d, -database,      --database              Set up basic database folder and specify number of servers.

-V, -verbose,       --verbose               Run script in verbose mode.
"""

def main(argv):
    verbose = False
    try:
        opts, args = getopt.getopt(argv,"d:hV",["database=", "help", "verbose"])
    except getopt.GetoptError:
        print(error_msg)
        sys.exit(2)

    if opts == []:
        print(error_msg)
        sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(error_msg)
        elif opt in ("-d", "--database"):
            server_count = int(arg)
        elif opt in ("-V, --verbose"):
            verbose = True

    if verbose:
        logging.basicConfig(level=logging.DEBUG)

    if not os.path.exists('database'):
        os.mkdir('database')
        logging.debug('created ./database folder')
    for i in range(server_count):
        if not os.path.exists('database/server{}'.format(str(i))):
            os.mkdir('database/server{}'.format(str(i)))
            logging.debug('created ./database/server{} folder'.format(str(i)))


if __name__ == "__main__":
   main(sys.argv[1:])