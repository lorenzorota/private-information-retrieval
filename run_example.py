import sys, getopt
import logging

def main(argv):
    verbose = False
    try:
        opts, args = getopt.getopt(argv,"r:hV",["run=", "help", "verbose"])
    except getopt.GetoptError:
        print('run_example.py --run <script> [-V]')
        sys.exit(2)

    if opts == []:
        sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print('run_example.py --run <script>')
            file = arg
        elif opt in ("-r", "--run"):
            file = arg
        elif opt in ("-V", "--verbose"):
            verbose = True
        
    if verbose:
        logging.basicConfig(level=logging.DEBUG)

    __import__('examples.{}'.format(file))

if __name__ == "__main__":
   main(sys.argv[1:])