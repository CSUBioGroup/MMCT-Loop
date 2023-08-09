import sys, os
import time
from utils import getLogger, getHelp
from ioput import parseBedpeFiles, txt2jd, parseJdFile
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
from HNSW3 import HNSW3 as HNSW
from estimate import estimateInterAndSelfCutFrag
from plot import plotInterSelfCutFrag
from model import getIntSig, markIntSig
import shutil

global log

def checkOverlap(peak_value, lx, ly):
    for i in range(len(peak_value)):
        px = peak_value[i][0]
        py = peak_value[i][1]
        if (px >= lx and px <= ly) or (py >= lx and py <= ly):
            return True
    return False

def singleHNSW(f, M, ef, k, bs, cut=0):
    dataI, readI, dataS, readS, dis, dss = [], [], [], [], [], []
    key, mat = parseJdFile(f, cut=0)
    if cut > 0:
        d = mat[:, 2] - mat[:, 1]
        p = np.where(d >= cut)[0]
        mat = mat[p, :]
        dss.extend(list(d[d < cut]))
    if len(mat) == 0:
        return key, f, dataI, dataS, list(dis), list(dss)
    # data for interaction records, read for readId
    report = "Clustering %s and %s, pre-set distance cutoff as > %s\n" % (key[0], key[1], cut)
    sys.stderr.write(report)

    print("Raw PETs count: ", len(mat))
    db = HNSW(mat, M, ef, k, bs)
    labels = pd.Series(db.labels)
    # print("****after cluster:", len(db.labels))

    mat = np.array(mat)
    mat = pd.DataFrame(mat[:, 1:].astype("float"),
                       index=mat[:, 0],
                       columns=["X", "Y"])
    nlabels = set(labels.values)
    print("Clusters count: ", len(nlabels))

    # collect cluster
    for label in nlabels:
        los = list(labels[labels == label].index)
        sub = mat.loc[los, :]
        # BEDPE format,+1 to escape the error that exact the same start and end
        # changed to remove such interactions
        if int(np.min(sub["X"])) == int(np.max(sub["X"])) \
                or int(np.min(sub["Y"])) == int(np.max(sub["Y"])):
            continue
        r = [
            key[0],
            int(np.min(sub["X"])),
            int(np.max(sub["X"])),
            key[1],
            int(np.min(sub["Y"])),
            int(np.max(sub["Y"])),
        ]

        print("i:", label, end="\r")

        if r[2] < r[4]:
            dataI.append(r)
            readI.extend(los)
        else:
            dataS.append(r)
            readS.extend(los)

    report = "Clustering %s and %s finished. Estimated %s self-ligation reads and %s inter-ligation reads\n" % (
        key[0], key[1], len(readS), len(readI))
    sys.stderr.write(report)
    if len(dataI) > 0:
        dis = mat.loc[readI, "Y"] - mat.loc[readI, "X"]
    if len(dataS) > 0:
        dss.extend(list(mat.loc[readS, "Y"] - mat.loc[readS, "X"]))
    return key, f, dataI, dataS, list(dis), list(dss)

def runHNSW(fs, M, ef, k, bs, cut=0, cpu=1):
    ds = Parallel(n_jobs=cpu)(delayed(singleHNSW)(f, M, ef, k, bs, cut) for f in fs)
    dataI, dataS, dis, dss = {}, [], [], []
    for d in ds:
        if len(d[2]) == 0:
            continue
        dataI[d[0]] = {"f": d[1], "records": d[2]}
        dataS.extend(d[3])
        dis.extend(d[4])
        dss.extend(d[5])
    return dataI, dataS, dis, dss

def combineCluster(dataI, dataI_2):
    """
    Combining multiple clustering result.
    """
    for key in dataI_2.keys():
        if key not in dataI:
            dataI[key] = {
                "f": dataI_2[key]["f"],
                "records": dataI_2[key]["records"]
            }
        else:
            ds = set()
            for r in dataI[key]["records"]:
                t = [r[1], r[2], r[4], r[5]]
                ds.add(tuple(t))
            for r in dataI_2[key]["records"]:
                t = tuple([r[1], r[2], r[4], r[5]])
                if t not in ds:
                    dataI[key]["records"].append(r)
    return dataI


def filterCluster(data, cut):
    """
    Filter inter-ligation clusters by distances
    """
    for key in data:
        nr = []
        for r in data[key]["records"]:
            d = (r[4] + r[5]) / 2 - (r[1] + r[2]) / 2
            if d >= cut:
                nr.append(r)
        data[key]["records"] = nr
    return data

def callpvalue(dataI, minPts, cut, cpu, fout):
    """
    Calling p-values of interactions for all chromosomes.
    """
    log.info("Starting estimate significance for interactions using distance cutoff as %s" % cut)
    ds = Parallel(n_jobs=cpu)(delayed(getIntSig)(dataI[key]["f"], dataI[key]["records"], minPts, cut)
                              for key in dataI.keys())
    ds = [d for d in ds if d is not None]
    if len(ds) == 0:
        log.error("Something wrong, no loops found, sorry, bye.")
        return 1
    ds = pd.concat(ds)
    try:
        # ds.sort_values(by="poisson_p-value")
        ds = markIntSig(ds)
        ds.to_csv(fout + ".loop", sep="\t", index_label="loopId")
    except:
        log.warning(
            "Something wrong happend to significance estimation, only output called loops"
        )
        ds.to_csv(fout + "_raw.loop", sep="\t", index_label="loopId")
    return 0

def pipe(fs,
         fout,
         M,
         ef,
         k,
         bs,
         minPts,
         cpu=1,
         chroms="",
         tmp=False,
         cut=0,
         plot=0
         ):
    if chroms == "":
        chroms = []
    else:
        chroms = set(chroms.split(","))
    if os.path.isdir(fout):
        mes = "Working directory %s exists, return." % fout
        log.error(mes)
        return
    os.mkdir(fout)
    MMCT_Loop_cluster_start_time = time.time()
    log.info("MMCT_Loop start clustering.")
    cfs = parseBedpeFiles(fs, fout, chroms, cut, log)
    cfs = Parallel(n_jobs=cpu)(delayed(txt2jd)(f) for f in cfs)
    dataI = {}
    dataI_2, dataS_2, dis_2, dss_2 = runHNSW(cfs, M, ef, k, bs, cut, cpu)
    if len(dataI_2) == 0:
        log.error("no inter-ligation PETs found")
        MMCT_Loop_cluster_end_time = time.time()
        log.info("MMCT_Loop cluster finished.")
        log.info("MMCT_Loop clustering used %s real cpu time." % (int(MMCT_Loop_cluster_end_time - MMCT_Loop_cluster_start_time)))
        return
    cuts = [cut, ]

    if len(dis_2) == 0 or len(dss_2) == 0:
        dataI = combineCluster(dataI, dataI_2)
    else:
        cut_2, frags = estimateInterAndSelfCutFrag(np.array(dis_2), np.array(dss_2))
        if plot:
            plotInterSelfCutFrag(dis_2, dss_2, cut_2, frags, prefix=fout + "disCutoff")
        log.info("Estimated inter-ligation and self-ligation distance cutoff")
        # experimental
        cuts.append(cut_2)
        cut = cut_2
        dataI = combineCluster(dataI, dataI_2)

    cuts = [c for c in cuts if c > 0]
    cut = np.min(cuts)

    dataI = filterCluster(dataI, cut)
    MMCT_Loop_cluster_end_time = time.time()
    log.info("MMCT_Loop cluster finished.")
    log.info("MMCT_Loop clustering used %s real time." % (int(MMCT_Loop_cluster_end_time - MMCT_Loop_cluster_start_time)))

    # estimate the significance
    MMCT_Loop_estimate_start_time = time.time()
    log.info("MMCT_Loop start estimating loop significance .")
    e = callpvalue(dataI, minPts, 0, cpu, fout)
    MMCT_Loop_estimate_end_time = time.time()
    log.info("MMCT_Loop estimate loop significance finished.")
    log.info("MMCT_Loop estimating loop significance used %s real cpu time." % (int(MMCT_Loop_estimate_end_time - MMCT_Loop_estimate_start_time)))

    if e and (tmp == False):
        shutil.rmtree(fout)
        return
    # remove temple files
    if tmp == False:
        shutil.rmtree(fout)

def main():
    stTime = time.time()
    global log
    log_file = os.path.join(os.getcwd(), "MMCT_Loop.log")
    log = getLogger(log_file)
    log.info("MMCT_Loop start at: " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(int(float(stTime)))) + ".")
    op = getHelp()
    report = "Command line:python MMCT-Loop.py -f {} -o {} -M {} -ef {} -k {} " \
             "-bs {} -minPts {} -cpu {} -chr {} -s {} -cut {} -plot {}".format(
        op.in_valid_bedpe, op.out_prefix, op.hnsw_M, op.hnsw_ef, op.hnsw_k, op.base_step, op.minPts, op.cpu, op.chromes, op.save, op.distance_cut, op.plot_dis)
    log.info(report)
    pipe(op.in_valid_bedpe.split(","), op.out_prefix, int(op.hnsw_M), int(op.hnsw_ef), int(op.hnsw_k),
         int(op.base_step), int(op.minPts), int(op.cpu), op.chromes, op.save, op.distance_cut, op.plot_dis)
    endTime = time.time()
    log.info("MMCT_Loop all steps finished.")
    log.info("MMCT_Loop total used %s real cpu time." % (int(endTime - stTime)))
    log.info("MMCT_Loop finished at: " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(int(float(endTime)))) + ".")
    
if __name__ == '__main__':
    main()