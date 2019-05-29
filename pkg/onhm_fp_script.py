from fponhm import FpoNHM
print('starting Script')
numdays = 2
fp = FpoNHM(2)
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