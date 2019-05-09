from pkg.fponhm import FpoNHM
fp = FpoNHM(2)
fp.initialize(r'../Data', r'../output')
fp.run_weights()
fp.finalize()
