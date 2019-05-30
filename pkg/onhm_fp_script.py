from fponhm import FpoNHM
import sys, getopt

def main():
    numdays = None
    print(sys.argv[0:])
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:", ["help"])

    except getopt.GetoptError as err:
        print(err)
        print('test.py -d <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("onhm_fp_script.py -d <numdays>")
        elif opt in ("-d"):
            numdays = int(arg)
    print("numdays = ", numdays)



    print('starting Script')
    #numdays = 2
    fp = FpoNHM(numdays)
    print('instantiated')
    ready = fp.initialize(r'../Data', r'../Output')
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