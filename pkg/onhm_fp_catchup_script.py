from fponhm import FpoNHM
import sys, getopt
from pathlib import Path
import datetime
import argparse

def main():
    numdays = None
    idir = None
    odir = None
    print(sys.argv[0:])
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:i:o:", ["help"])

    except getopt.GetoptError as err:
        print(err)
        print("onhm_fp_script.py -d <numdays>")
        sys.exit(2)
    if opts == []:
        print("missing option - Usage: onhm_fp_script.py -d <numdays>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("onhm_fp_script.py -d <numdays> -i <input_dir>, -o <output_dir>")
        elif opt in ("-d"):
            numdays = int(arg)
        elif opt in ('-i'):
            idir = Path(arg)
        elif opt in ('-o'):
            odir = Path(arg)
    print("numdays = ", numdays)

    numdays = (datetime.date.today() - datetime.date(2015,1,1)).days

    print('starting Script')
    #numdays = 2
    fp = FpoNHM(numdays)
    print('instantiated')
    ready = fp.initialize(idir, odir)
    if ready:
        print('initalized\n')
        print('running')
        fp.run_weights()
        print('finished running')
        fp.finalize()
        print('finalized')
    else:
        print('Gridmet not updated continue with numdays -1')
        fp.setNumdays(numdays-1)
        print('initalized\n')
        print('running')
        fp.run_weights()
        print('finished running')
        fp.finalize()
        print('finalized')

if __name__ == "__main__":
    main()