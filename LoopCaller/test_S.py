import sys, os
from time import time
from utils import getLogger, getHelp
from ioput import parseBedpeFiles, txt2jd, parseJdFile
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
from HNSW3 import HNSW3 as HNSW
from estimate import estimateInterAndSelfCutFrag
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


def singleLoadLoop(f, M, ef, k, le, r, minPts, fp, cut=0):
    dataI, readI, dataS, readS, dis, dss = [], [], [], [], [], []
    key, mat = parseJdFile(f, cut=0)
    if cut > 0:
        d = mat[:, 2] - mat[:, 1]
        p = np.where(d >= cut)[0]
        mat = mat[p, :]
        dss.extend(list(d[d < cut]))
    if len(mat) == 0:
        return key, f, dataI, dataS, list(dis), list(dss)

    db = HNSW(mat, M, ef, k, le, r, minPts)
    labels = pd.Series(db.labels)

    mat = np.array(mat)
    mat = pd.DataFrame(mat[:, 1:].astype("float"),
                       index=mat[:, 0],
                       columns=["X", "Y"])
    nlabels = set(labels.values)

    # collect cluster
    flag_peak = False
    if len(fp) != 0:
        flag_peak = True
        peak_context = open(fp)
        peaks = {}
        for line in peak_context:
            line = line.split("\n")[0].split("\t")
            peak = line[0]
            start = int(line[1])
            end = int(line[2])
            peaks.setdefault(peak, []).append([start, end])

    t = time()
    c0, c1, c2 = 0, 0, 0
    for label in nlabels:
        los = list(labels[labels == label].index)
        sub = mat.loc[los, :]
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

        if flag_peak:
            try:
                flag1 = checkOverlap(peaks[r[0]], r[1], r[2])
                flag2 = checkOverlap(peaks[r[3]], r[4], r[5])
                if flag1 is True and flag2 is True:
                    c2 += 1
                elif flag1 is False and flag2 is False:
                    c0 += 1
                    continue
                else:
                    c1 += 1
                    continue
            except:
                print("chrom no peak data!", key[0])
                # logger.warning("chrom no peak data!", key[0])

        if r[2] < r[4]:
            dataI.append(r)
            readI.extend(los)
        else:
            dataS.append(r)
            readS.extend(los)

    if len(dataI) > 0:
        dis = mat.loc[readI, "Y"] - mat.loc[readI, "X"]
    if len(dataS) > 0:
        dss.extend(list(mat.loc[readS, "Y"] - mat.loc[readS, "X"]))
    return key, f, dataI, dataS, list(dis), list(dss)


def getCandidateLoop(fs, M, ef, k, le, r, minPts, fp, cut=0, cpu=1, loop_filename=""):
    # ds = Parallel(n_jobs=cpu)(delayed(singleLoadLoop)(f, M, ef, k, le, r, minPts, fp, cut) for f in fs)
    dataI, dataS, dis, dss = {}, {}, [], []
    dis_t = {}
    dss_t = {}

    for f in fs:
        print(f)
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
                print('Loop file format error ', line)
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
                print('Error: not a inter-ligation loop ', line)
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
        ds.to_csv(fout + "_S.loop", sep="\t", index_label="loopId")
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
    if len(ds) > 0:
        frags = np.median(ds)
    else:
        frags = 0
    rfrags = int(2 ** frags)
    return rcut, rfrags

def pipe(fs,
         fout,
         M,
         ef,
         k,
         le,
         r,
         minPts,
         fpeak="",
         cpu=1,
         chroms="",
         tmp=True,
         cut=0,
         plot=0,
         loop_filename=""):
    if chroms == "":
        chroms = []
    else:
        chroms = set(chroms.split(","))
    if os.path.isdir(fout):
        mes = "Working directory %s exists, return." % fout
        log.error(mes)
        return
    os.mkdir(fout)


    cfs = parseBedpeFiles(fs, fout, chroms, cut, log)
    cfs = Parallel(n_jobs=cpu)(delayed(txt2jd)(f) for f in cfs)
    dataI = {}
    dataS = {}
    dataI_2, dataS_2, dis_2, dss_2 = getCandidateLoop(cfs, M, ef, k, le, r, minPts, fpeak, cut, cpu, loop_filename)
    if len(dataI_2) == 0:
        log.info("ERROR: no inter-ligation PETs")
        return
    cuts = [cut, ]
    if len(dis_2) == 0 or len(dss_2) == 0:
        dataI = combineCluster(dataI, dataI_2)
    else:
        cut_2, frags = estimateInterAndSelfCutFrag_2(np.array(dis_2), np.array(dss_2))
        if plot:
            plotInterSelfCutFrag(dis_2, dss_2, cut_2, frags, prefix=fout + "disCutoff")
        log.info("Estimated inter-ligation and self-ligation distance cutoff")
        # experimental
        cuts.append(cut_2)
        cut = cut_2
        dataS = combineCluster(dataS, dataS_2)

    cuts = [c for c in cuts if c > 0]
    cuts.append(0)
    cut = np.min(cuts)

    # dataI = filterCluster(dataI, cut)

    # 5.estimate the significance
    # TODO
    e = callpvalue(dataS, minPts, 0, cpu, fout)
    if e and (tmp == False):
        shutil.rmtree(fout)
        return
    # 7.remove temple files
    if tmp == False:
        shutil.rmtree(fout)

def startCalPvalue(loop_filename):
    stTime = time()
    print("Start at:", stTime, ".")
    global log
    fileName = os.path.join(os.getcwd(), "LoopCaller_test_self-ligation.log")
    log = getLogger(fileName)
    op = getHelp()
    report = "Command line: test_S -f {} -o {} -M {} -ef {} -k {} -le {} " \
             "-r {} -minPts {} -p {} -cpu {} -c {} -s {} -cut {} -plot {}".format(
        op.fnIn, op.fnOut, op.M, op.ef, op.k, op.le,
        op.resolution, op.minPts, op.fnPeak, op.cpu, op.chroms, op.tmp, op.cut, op.plot)
    log.info(report)
    pipe(op.fnIn.split(","), op.fnOut, op.M, op.ef, op.k, op.le,
         op.resolution, op.minPts, op.fnPeak, op.cpu, op.chroms, op.tmp, op.cut, op.plot, loop_filename)
    endTime = time()
    log.info("LoopCaller test self-ligation finished.\nUsed %s CPU time.\n" % (endTime-stTime))


