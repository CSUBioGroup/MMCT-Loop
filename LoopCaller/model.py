"""
Statistical significance is tested for every chromosome using the local permutated background.
"""

import gc
import numpy as np
import pandas as pd
from scipy.stats import poisson
from ioput import parseJdFile, parseIv
from utils import mesFlush

def getCorLink(cs):
    """
    @param cs: [1,2,3,4], a list for the coordinates x or y
    @return ts_keys: [ 1,2,3,4,5]
            ts: { 1:0,2:1,3:2,4:3,5:4 }
    """
    ts = {}
    for i, c in enumerate(cs):
        ts.setdefault(c, []).append(i)
    ts_keys = np.sort(cs)
    return ts_keys, ts

def getGenCov(f, cut=0):
    """
    Build the genomic model for random access. Could use a lot of memory.
    @param f:.jd file 
    @param cut: distance cutoff for self-ligation PETs.
    """
    key, mat = parseJdFile(f, cut)
    j = mat.shape[0]
    if j < 2:
        return None, 0
    xs_keys, xs = getCorLink(mat[:, 1])
    ys_keys, ys = getCorLink(mat[:, 2])
    return [[xs_keys, xs], [ys_keys, ys]], j

def getCounts(iv, model):
    ps = []
    ts_keys, ts = model
    l_idx = np.searchsorted(ts_keys, iv[0], side="left")
    r_idx = np.searchsorted(ts_keys, iv[1], side="right")
    for i in range(l_idx, r_idx):
        ps.extend(ts[ts_keys[i]])
    return set(ps)

def getPETsforRegions(iva, ivb, model):
    raSource = getCounts(iva, model[0])
    raTarget = getCounts(iva, model[1])
    rbSource = getCounts(ivb, model[0])
    rbTarget = getCounts(ivb, model[1])
    ra = len(raSource.union(raTarget))
    rb = len(rbSource.union(rbTarget))
    rab = len(raSource.intersection(rbTarget))
    return ra, rb, rab

def getNearbyPairRegions(iva, ivb, win=5):
    """
    @param iva: [start,end] 
    Get the nearby regions for interacting two locus,win as how many nearby, 6 is enough for interacting more than 100 regions to estimate FDR and others. The mean distance of all the permutated regions is the same to that between iva and ivb.
    """
    ivas, ivbs = [], []
    ca = sum(iva) / 2
    cb = sum(ivb) / 2
    sa = (iva[1] - iva[0]) / 2
    sb = (ivb[1] - ivb[0]) / 2
    step = (sa + sb) / 2
    for i in range(0 - win, win + 1):
        if i == 0:
            continue
        niva = [iva[0], iva[1]]
        niva[0] = max([0, ca + i * step - sa])
        niva[1] = max([0, ca + i * step + sa])
        nivb = [ivb[0], ivb[1]]
        nivb[0] = max([0, cb + i * step - sb])
        nivb[1] = max([0, cb + i * step + sb])
        ivas.append(niva)
        ivbs.append(nivb)
    return ivas, ivbs

def getMultiplePsFdr1(iva, ivb, model, win=5):
    """
    LoopCallerFDR
    for the interval a and b, searching its nearby windows to estimate FDR and p-values.
    return ra, rb, rab, es,es_ra, es_rb, fdr, pop
    """
    ra, rb, rab = getPETsforRegions(iva, ivb, model)
    ivas, ivbs = getNearbyPairRegions(iva, ivb, win=win)
    # nras is a list for storing points ids for permutated regions
    nras, nrbs = [], []
    for na in ivas:
        nraSource = getCounts(na, model[0])
        nraTarget = getCounts(na, model[1])
        nra = nraSource.union(nraTarget)
        nras.append(nra)
    for nb in ivbs:
        nrbSource = getCounts(nb, model[0])
        nrbTarget = getCounts(nb, model[1])
        nrb = nrbSource.union(nrbTarget)
        nrbs.append(nrb)
    # caculating the permutated background
    rabs = []
    for nra in nras:
        for nrb in nrbs:
            nrab = float(len(nra.intersection(nrb)))
            if nrab > 0:
                # collect the value for poisson test
                rabs.append(nrab)
            else:
                rabs.append(0.0)
    if len(rabs) == 0:
        return ra, rb, rab, np.inf, 0.0, 0.0, 1e-300
    rabs = np.array(rabs)
    # local fdr
    fdr = len(rabs[rabs > rab]) / float(len(rabs))
    mrabs = float(np.mean(rabs))
    # enrichment score
    if mrabs > 0:
        es = rab / np.mean(rabs[rabs > 0])
    else:
        es = np.inf
    # simple possion test
    lam = mrabs
    pop = max([1e-300, poisson.sf(rab - 1.0, lam)])
    return ra, rb, rab, es, fdr, pop

def getBonPvalues(ps):
    """
    Return the Bonferroni corrected p-values.
    """
    ps = np.array(ps)
    ps = ps * len(ps)
    ps[ps > 1.0] = 1.0
    return ps

def checkOneEndOverlap(xa, xb, ya, yb):
    """
    check the overlap of a region for the same chromosome
    """
    if (ya <= xa <= yb) or (ya <= xb <= yb) or (ya <= xa <= xb <= yb):
        return True
    if (xa <= ya <= xb) or (xa <= yb <= xb) or (xa <= ya <= yb <= xb):
        return True
    return False

def checkOverlap(ivai, ivbi, ivaj, ivbj):
    """
    check the overlap of two anchors,ra=[chr,left_start,left_end,right_start,right_end]
    """
    if ivai[0] != ivaj[0] or ivbi[0] != ivbj[0]:
        return
    if checkOneEndOverlap(ivai[1], ivai[2], ivaj[1], ivaj[2]) \
            and checkOneEndOverlap(ivbi[1], ivbi[2], ivbj[1], ivbj[2]):
        return True
    return False

def removeDup(ds, ppcut=1e-5):
    """
    Remove overlapped called loops, keep the more significant one for multiple eps result. 
    @param:ds, from getIntSig
    """
    uniqueds = {}
    reds = {}
    rekeys = set()
    keys = list(ds.keys())
    for i in range(len(keys) - 1):
        keyi = keys[i]
        if keyi in rekeys:
            continue
        ivai = parseIv(ds[keyi]["iva"])
        ivbi = parseIv(ds[keyi]["ivb"])
        # 1 means unique loops
        flag = 1
        # collect overlapped loops
        for j in range(i + 1, len(keys)):
            keyj = keys[j]
            if keyj in rekeys:
                continue
            ivaj = parseIv(ds[keyj]["iva"])
            ivbj = parseIv(ds[keyj]["ivb"])
            flagj = checkOverlap(ivai, ivbi, ivaj, ivbj)
            # there is overlapped loops,collect them
            if flagj:
                if keyi not in reds:
                    reds[keyi] = [keyi]
                    rekeys.add(keyi)
                reds[keyi].append(keyj)
                rekeys.add(keyj)
                flag = 0
        # collect unique loops
        if flag:
            uniqueds[keyi] = ds[keyi]
    # for overlapped loops, choose the more significant ones
    for key in reds.keys():
        ts = {}
        for t in reds[key]:
            if ds[t]["poisson_p-value"] > ppcut:
                continue
            # first select the significant loops, then select the loops with smaller anchors and higher density
            ts[t] = float(ds[t]["rab"]) / ds[t]["ra"] / ds[t]["rb"]

        if len(ts) == 0:
            continue
        ts = pd.Series(ts)
        ts.sort_values(inplace=True, ascending=False)
        uniqueds[ts.index[0]] = ds[ts.index[0]]
    return uniqueds

def getIntSig(f, records, minPts, discut):
    """
    @param:discut, distance cutoff determined for self-ligation pets.
    """
    print("Starting estimate significance for %s candidate interactions in %s" % (len(records), f))
    model, N = getGenCov(f, discut)
    print("Genomic coverage model built from %s" % f)
    if N == 0:
        print("No cis-PETs parsed as requiring distance cutoff >%s from %s" %
              (discut, f))
        return None
    ds = {}
    i = 0
    for r in records:
        chrom = r[0]
        key = "%s-%s-%s" % (r[0], r[3], i)
        iva = [max(0, r[1]), r[2]]
        ivb = [max(0, r[4]), r[5]]
        # filter loops
        distance = abs(sum(ivb) / 2.0 - sum(iva) / 2.0)
        if distance < discut:
            continue
        ra, rb, rab = getPETsforRegions(iva, ivb, model)
        # filter clusters contain many self-ligation PETs within distance cutoff
        if rab < int(minPts):
            continue
        i += 1
        mesFlush("%s interaction p-values estimated for %s" % (i, f))
        # 只保留poisson_p-value
        ra, rb, rab, es, fdr, pop = getMultiplePsFdr1(iva, ivb, model)
        # this part should be furthur modified, as for most ideable data, there are no noise, so the es should be inf, however, not possible

        ds[key] = {
            "distance": distance,
            "ra": ra,
            "rb": rb,
            "rab": rab,
            "ES": es,
            "FDR": fdr,
            "poisson_p-value": pop,
            "iva": "%s:%s-%s" % (chrom, iva[0], iva[1]),
            "ivb": "%s:%s-%s" % (chrom, ivb[0], ivb[1])
        }

    print()
    # memory usage
    del model
    gc.collect()
    if len(ds.keys()) == 0:
        return None
    ds = removeDup(ds)
    if len(ds.keys()) == 0:
        return None
    ds = pd.DataFrame(ds).T
    return ds

def markIntSig(ds, escut=3.0, ppcut=1e-5,):
    # filter data according to cutoffs
    # larger enrichment score
    a = ds["ES"]
    a = a[a >= escut]
    # smaller poisson
    b = ds.loc[a.index, "poisson_p-value"]
    b = b[b <= ppcut]
    rs = b.index
    ns = pd.Series(data=np.zeros(ds.shape[0]), index=ds.index)
    ns[rs] = 1.0
    ds["significant"] = ns
    return ds