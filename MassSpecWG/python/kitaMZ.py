import pandas as pd
import numpy as np
import h5pyd as h5py
import random


def get_names(h5):
    return list(h5.keys())


def get_names2(h5):
    h5dt = h5.dtype
    cols = list(h5dt.names)
    return cols


def read_hdf(h5, index=None, names=None):
    if isinstance(h5, h5py._hl.dataset.Dataset):
        return read_hdf2(h5, index, names)
    cols = get_names(h5) if names is None else names
    if index is None:
        return pd.DataFrame({x: h5[x][:] for x in cols})
    else:
        if index.stop - index.start == 0:
            return None
        return pd.DataFrame({x: h5[x][index] for x in cols})


def read_hdf2(h5, index=None, names=None):
    cols = get_names2(h5) if names is None else names
    if index is None:
        tret = h5[:]
        return pd.DataFrame({x: tret[x] for x in cols})
    else:
        tret = h5[:]
        return pd.DataFrame({x: tret[x] for x in cols})


def gen_range_list(df):
    ranges = [
        slice(int(x), int(x + y)) for x, y in zip(list(df["offset"]), list(df["N"]))
    ]
    ranges = [x for x in ranges if x.stop - x.start > 0]
    return ranges


def spec_pred(target, query, delta):
    result = np.zeros(target.size)
    qbi = 0
    q_b = query[qbi]
    qei = len(query)
    min_ql = q_b - delta * (q_b)
    max_ql = q_b + delta * (q_b)
    for i in range(len(target)):
        t_t = target[i]
        while t_t > max_ql:
            if qbi == (qei - 1):
                return result
            else:
                qbi = qbi + 1
                min_ql = query[qbi] - delta * query[qbi]
                max_ql = query[qbi] + delta * query[qbi]
        if (t_t <= max_ql) and (t_t > min_ql):
            result[i] = 1
    return result


def spec_pred_range(target_low, target_high, query, delta):
    result = np.zeros(target_low.size)
    qbi = 0
    q_b = query[qbi]
    qei = len(query)
    min_ql = q_b - delta * (q_b)
    max_ql = q_b + delta * (q_b)
    for i in range(len(target_low)):
        t_tl = target_low[i]
        t_th = target_high[i]
        while t_tl > max_ql:
            if qbi == (qei - 1):
                return result
            else:
                qbi = qbi + 1
                min_ql = query[qbi] - delta * query[qbi]
                max_ql = query[qbi] + delta * query[qbi]
        if (t_tl <= max_ql) and (min_ql <= t_th):
            result[i] = 1
    return result


def filter_spec(target, query, delta):
    return target[spec_pred(target, query, delta) == 1]


def filter_df(df, mz, delta):
    mzv = df["MZ"].values
    pred = spec_pred(mzv, mz, delta) == 1
    return df[pred]


def filter_df_range(df, mz, delta):
    mzl = df["MZLow"].values
    mzh = df["MZHigh"].values
    return df[spec_pred_range(mzl, mzh, mz, delta) == 1]


def getXIC(h5, mz, ppm=100):
    delta = ppm * 0.000001
    tbl_cols = ["MZLow", "MZHigh", "offset", "N", "ID"]
    table = read_hdf(h5["trace_table"], None, tbl_cols)
    iddf = filter_df_range(table, mz, delta)
    rl = gen_range_list(iddf)
    return pd.concat(
        [filter_df(read_hdf(h5["trace"], x, ["MZ"]), mz, delta) for x in rl]
    )


def getMZ(h5, seed=123, N=25, scans=None):
    table = read_hdf(h5["spectra_table"], None, ["ScanNumber"])
    maxscan = max(table.ScanNumber)
    scans = (
        [int(random.uniform(0, maxscan) + 1) for i in range(N)]
        if scans is None
        else scans
    )
    spec_df = (
        getSpectrum(h5, scans, cols=["Intensity", "MZ"])
        .sort_values(by=["Intensity"], ascending=False)
        .head(100)["MZ"]
        .values
    )
    spec_df.sort()
    return spec_df


def getXICFromSpectrum(h5, mz, ppm=100):
    delta = ppm * 0.000001
    table = read_hdf(h5["spectra_table"])
    fullscans = table.query("MSOrder==1")
    rl = gen_range_list(fullscans)
    return pd.concat(
        [filter_df(read_hdf(h5["spectra"], x, ["MZ"]), mz, delta) for x in rl]
    )


def XICFromSpectrumTest(
    h5, mzmin=400, mzmax=1000, ppm=100, N=100, seed=123, scans=None
):
    random.seed(seed)
    mzs = getMZ(h5, seed=seed, scans=scans)
    return len(getXICFromSpectrum(h5, mzs, ppm))


def getSpectrum(h5, scans, cols=None):
    cols = get_names(h5["spectra"]) if cols is None else cols
    spec_df = read_hdf(
        h5["spectra_table"], None, names=["ScanNumber", "offset", "N"]
    ).iloc[[x - 1 for x in scans]]
    assert list(spec_df["ScanNumber"]) == scans
    ranges = gen_range_list(spec_df)
    return pd.concat([read_hdf(h5["spectra"], x, cols) for x in ranges])


def XICTest(h5, mzmin=400, mzmax=1000, ppm=100, N=100, seed=123, scans=None):
    random.seed(seed)
    mzs = getMZ(h5, seed=seed, scans=scans)
    df = getXIC(h5, mzs, ppm)
    npoints = len(df)
    return npoints


def spectrumTest(h5, seed=123, N=100, scans=None):
    ret = getMZ(h5, seed, N, scans)
    return len(ret)
