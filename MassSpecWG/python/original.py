import pandas as pd
import h5py
from pathlib import Path
import time
import random


def getSpectrum(fn, scan):
    try:
        scan = pd.read_hdf(fn, "spectra/S" + str(scan))
        return scan
    except:
        return None


def getXIC(fn, mz, ppm=100):
    delta = ppm * mz * 0.000001
    mzlo = mz - delta
    mzhi = mz + delta
    table = pd.read_hdf(fn, "trace_table")
    idlo = int(table.query("MZLow<=" + str(mzlo)).tail(1).ID)
    idhi = int(table.query("MZHigh>=" + str(mzhi)).head(1).ID)

    df = pd.DataFrame()
    with h5py.File(fn, "r") as f:
        for i in range(idlo, idhi + 2):
            try:
                trace = pd.read_hdf(fn, f"trace/T{i}").query(
                    f"MZ <= {mzhi} and MZ>= {mzlo}"
                )
                df = df.append(trace)
            except:
                continue
    return df


def XICTest(fn, mzmin=400, mzmax=1000, ppm=100, N=100, seed=123, scans=None):
    random.seed(seed)
    mzs = getMZ(fn, seed=seed, scans=scans)
    npoints = 0
    start = time.time()
    for i in range(N):
        mz = mzs.iloc[i]
        df = getXIC(fn, mz, ppm)
        npoints += len(df)
    end = time.time()
    note(fn, "XIC test time per XIC", (end - start) / N)
    note(fn, "Average points per XIC", npoints / N)
    return (end - start) / N


def XICFromSpectrumTest(fn, mzmin=400, mzmax=1000, ppm=100, N=100, seed=123):
    random.seed(seed)
    mzs = getMZ(fn, seed=seed)
    npoints = 0
    start = time.time()
    for i in range(N):
        mz = mzs.iloc[i]
        df = getXICFromSpectrum(fn, mz, ppm)
        npoints += len(df)
    end = time.time()
    note(fn, "XIC from Spectrum test time per XIC", (end - start) / N)
    note(fn, "Average points per XIC from Spectrum", npoints / N)
    return (end - start) / N


def spectrumTest(fn, N=100, seed=123, scans=None):
    table = pd.read_hdf(fn, "spectra_table")
    maxscan = max(table.ScanNumber)
    random.seed(seed)
    start = time.time()
    if scans is None:
        scans = [int(random.random() * maxscan + 1) for i in range(N)]
    for scan in scans:
        try:
            getSpectrum(fn, scan)
        except:
            print("error in getSpectrum")
            continue
    end = time.time()
    note(fn, "Spectra test time per spectra", (end - start) / N)
    return (end - start) / N


def getXICFromSpectrum(fn, mz, ppm=100):
    delta = ppm * mz * 0.000001
    mzlo = mz - delta
    mzhi = mz + delta
    table = pd.read_hdf(fn, "spectra_table")
    fullscans = table.query("MSOrder==1").ScanNumber
    df = pd.DataFrame()
    for i in fullscans:
        x = getSpectrum(fn, i)
        for j in range(len(x)):
            xmz = x.MZ[j]
            if (xmz < mzhi) and (xmz > mzlo):
                df = df.append(x.iloc[j])
    return df


def getMZ(fn, seed=123, scans=None):
    table = pd.read_hdf(fn, "spectra_table")
    maxscan = max(table.ScanNumber)
    random.seed(seed)
    start = time.time()
    df = pd.DataFrame()

    if scans is None:
        N = 25
        scans = [int(random.random() * maxscan + 1) for i in range(N)]
    else:
        N = len(scans)

    for scan in scans:
        try:
            df = df.append(getSpectrum(fn, scan))
        except:
            print("error in getSpectrum")
            continue
    end = time.time()
    note(fn, "Spectra test time per spectra", (end - start) / N)
    return df.sort_values(by=["Intensity"], ascending=False).head(100)["MZ"]


annotations = pd.DataFrame(columns=["Source", "Comment", "Value"])


def note(src, com, val):
    dic = {"Source": str(Path(src).stem), "Comment": com, "Value": val}
    annotations.loc[len(annotations)] = dic


def clearNote():
    annotations = pd.DataFrame(columns=["Source", "Comment", "Value"])


# def Benchmark(files, seed=123, N=100):
#     for fn in files:
#         #        outfile = getH5FileName(fn)
#         print("Testing " + fn)
#         spectrumTest(fn, N=N, seed=seed)
#         XICTest(fn, seed=seed, N=N)
#         XICFromSpectrumTest(fn, N=N, seed=seed)
#     return annotations
