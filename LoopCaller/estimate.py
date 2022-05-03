from collections import Counter
import numpy as np
import pandas as pd

def estimateFragSize(ds, top=500):
    """
    Estimation PETs fragment size.
    @param ds: list, PET distance
    """
    ds = pd.Series(Counter(ds))
    ds.sort_values(inplace=True, ascending=False)
    ds = ds[:top]
    frags = int(np.median(ds.index))
    return frags

def estimateInterAndSelfCutFrag(di, ds, log=1):
    """
    Estimation of distance cutoff for inter-ligation and self-ligation pets.
    @param di: list,distance for inter-ligation cluster pets
    @param ds: list,distance for self-ligation cluster pets
    """
    di = np.abs(np.array(di))
    ds = np.abs(np.array(ds))
    di = di[~np.isnan(di)]
    ds = ds[~np.isnan(ds)]
    di = di[di > 0]
    ds = ds[ds > 0]
    if log:
        di = np.log2(di)
        ds = np.log2(ds)
    # self-ligation and inter-ligation distance cutoff
    cut1 = np.median(ds) + 3 * ds.std()
    cut2 = (ds.mean() * ds.std() + di.mean() * di.std()) / (ds.std() + di.std())
    cut = min([cut1, cut2])
    rcut = int(2 ** cut)
    # fragment size
    frags = np.median(ds)
    rfrags = int(2 ** frags)
    return rcut, rfrags
