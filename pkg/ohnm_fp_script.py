from pkg import fponhm

fp = fponhm.FpoNHM(2)
fp.initialize(r'../Data', r'../output')
fp.run_weights()
fp.finalize()
