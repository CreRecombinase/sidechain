import h5py
import numpy as np
import original
import refactor
import io
import cProfile
import pstats

def Benchmark_original(fn, seed=123, N=100, scans=None):
    print("Testing " + str(fn))
    original.spectrumTest(fn, N=N, seed=seed, scans=scans)
    original.XICTest(fn, seed=seed, N=N, scans=scans)
    original.XICFromSpectrumTest(fn, N=N, seed=seed)
    return True


def Benchmark_refactor(fn, seed=123, N=100, scans=None):
    with h5py.File(fn, "r") as h5:
        print("Testing " + str(h5))
        refactor.spectrumTest(h5=h5, N=N, seed=seed, scans=scans)
        refactor.XICTest(h5, seed=seed, N=N, scans=scans)
        refactor.XICFromSpectrumTest(h5, N=N, seed=seed)
    return True



original_file = "../data/blacktea_hcd_7500_top3_1.h5"
refactor_file = "../new_data/blacktea_hcd_7500_top3_1.h5"

pr_orig = cProfile.Profile()
pr_orig.enable()
Benchmark_original(original_file)
pr_orig.disable()
pr_orig.dump_stats('../profiles/original.stat')



pr_refactor = cProfile.Profile()
pr_refactor.enable()
Benchmark_refactor(refactor_file)
pr_refactor.disable()
pr_refactor.dump_stats('../profiles/refactor.stat')
