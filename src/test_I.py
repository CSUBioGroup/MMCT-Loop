import os
import time
from utils import getLogger, getEstimatingArg
from ioput import parseBedpeFiles, txt2jd, parseJdFile
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
from plot import plotInterSelfCutFrag
from model_test import getIntSig, markIntSig
import shutil

global log


# checkOverlap( [chr1_peak_STi,chr1_peak_ENi], minX, maxX)

def checkOverlap(peak_value, lx, ly):
    for i in range(len(peak_value)):
        px = peak_value[i][0]
        py = peak_value[i][1]
        if (px >= lx and px <= ly) or (py >= lx and py <= ly):
            return True
    return False



def getCandidateLoop(fs, loop_filename=""):

    dataI, dataS, dis, dss = {}, {}, [], []
    dis_t = {}
    dss_t = {}

    for f in fs:
        chr = f.split('.')[0].split('-')[1]
        key = tuple([chr, chr])
        dataI[key] = {"f" : "", "records":[]}
        dataS[key] = {"f" : "", "records":[]}
        dis_t[key] = []
        dss_t[key] = []

    with open(loop_filename, 'r') as fin:
        while True:
            line = fin.readline()
            if line == "":
                break
            line = line.split('\n')[0].split('\t')
            if len(line) < 6:
                log.error("Loop file format error.")
                return
            r = [
                line[0],
                int(line[1]),
                int(line[2]),
                line[3],
                int(line[4]),
                int(line[5]),
            ]
            if r[0] != r[3]:
                log.error("Error: not a inter-ligation loop")
                return
            key = tuple([r[0], r[3]])
            if r[2] < r[4]:
                dataI[key]["records"].append(r)
                dis_t[key].append((r[4] + r[5] - r[1] - r[2]) / 2)
            else:
                dataS[key]["records"].append(r)
                dss_t[key].append((r[4] + r[5] - r[1] - r[2]) / 2)

    for f in fs:
        chr = f.split('.')[0].split('-')[1]
        key = tuple([chr, chr])
        dataI[key]["f"] = f
        dataS[key]["f"] = f
        # dataS.extend(dataS_t[key]["records"])
        if len(dis_t[key]) > 0:
            dis.extend(list(dis_t[key]))
        if len(dss_t[key]) > 0:
            dss.extend(list(dss_t[key]))
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


def callpvalue(dataI, cut, cpu, fout):
    """
    Calling p-values of interactions for all chromosomes.
    """
    log.info("Starting estimate significance for interactions using distance cutoff as %s" % cut)
    ds = Parallel(n_jobs=cpu)(delayed(getIntSig)(dataI[key]["f"], dataI[key]["records"], cut)
                              for key in dataI.keys())
    ds = [d for d in ds if d is not None]
    if len(ds) == 0:
        log.error("Something wrong, no loops found, sorry, bye.")
        return 1
    ds = pd.concat(ds)
    try:
        # ds.sort_values(by="poisson_p-value")
        ds = markIntSig(ds)
        ds.to_csv(fout + "_intra_ligation.loop", sep="\t", index_label="loopId")
    except:
        log.warning(
            "Something wrong happend to significance estimation, only output called loops"
        )
        ds.to_csv(fout + "_raw.loop", sep="\t", index_label="loopId")
    return 0

def estimateInterAndSelfCutFrag_2(di, ds, log=1):
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
    if len(ds) > 0 and len(di) > 0:
        cut1 = np.median(ds) + 3 * ds.std()
        cut2 = (ds.mean() * ds.std() + di.mean() * di.std()) / (ds.std() + di.std())
    else:
        if len(ds) > 0:
            cut1 = np.median(ds) + 3 * ds.std()
            cut2 = ds.mean() * ds.std() / ds.std()
        if len(di) > 0:
            cut1 = np.median(di) + 3 * di.std()
            cut2 = di.mean() * di.std() / di.std()
    cut = min([cut1, cut2])
    rcut = int(2 ** cut)
    # fragment size
    if len(di) > 0:
        frags = np.median(di)
    else:
        frags = 0
    rfrags = int(2 ** frags)
    return rcut, rfrags

def pipe(fs,
         fout,
         cpu=1,
         chroms="",
         cut=0,
         loop_filename=""):
    if chroms == "":
        chroms = []
    else:
        chroms = set(chroms.split(","))

    if os.path.isdir(fout):
        mes = "Working directory %s exists, remove old." % fout
        log.warning(mes)
        shutil.rmtree(fout)
    os.mkdir(fout)

    cfs = parseBedpeFiles(fs, fout, chroms, cut, log)
    cfs = Parallel(n_jobs=cpu)(delayed(txt2jd)(f) for f in cfs)
    dataI = {}
    dataS = {}
    dataI_2, dataS_2, dis_2, dss_2 = getCandidateLoop(cfs, loop_filename)
    if len(dataI_2) == 0:
        log.error("No inter-ligation PETs found.")
        return
    cuts = [cut, ]
    if len(dis_2) == 0 or len(dss_2) == 0:
        dataI = combineCluster(dataI, dataI_2)
    else:
        cut_2, frags = estimateInterAndSelfCutFrag_2(np.array(dis_2), np.array(dss_2))
        log.info("Estimated inter-ligation and self-ligation distance cutoff")
        # experimental
        cuts.append(cut_2)
        cut = cut_2
        dataI = combineCluster(dataI, dataI_2)

    cuts = [c for c in cuts if c > 0]
    cuts.append(0)
    cut = np.min(cuts)

    # dataI = filterCluster(dataI, cut)
    # estimate the significance
    e = callpvalue(dataI, 0, cpu, fout)
    shutil.rmtree(fout)

def startCalPvalue():
    global log
    log_file = os.path.join(os.getcwd(), "LoopCaller.log")
    log = getLogger(log_file)
    op = getEstimatingArg()
    st_time = time.time()
    log.info("LoopCaller start estimating inter-ligation loops from " + op.test_file + " at: " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(int(float(st_time)))))
    report = "Command line: test_I.py -f {} -o {} -cpu {} -chr {} -cut {} -tf {}".format(
        op.in_valid_bedpe, op.out_prefix, op.cpu, op.chromes, op.distance_cut, op.test_file)
    log.info(report)
    pipe(op.in_valid_bedpe.split(","), op.out_prefix, int(op.cpu), op.chromes, int(op.distance_cut), op.test_file)
    en_time = time.time()
    log.info("LoopCaller estimating inter-ligation loops finished at: " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(int(float(en_time)))))
    log.info("LoopCaller estimating inter-ligation loops total used %s real cpu time." % (int(en_time - st_time)))

startCalPvalue()