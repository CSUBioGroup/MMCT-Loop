import gzip
import os
import numpy as np
from utils import mesFlush
import joblib

class PET(object):
    __slots__ = [
        "chromL", "chromR",
        "startL", "startR",
        "endL", "endR",
        "strandL", "strandR",
        "cL", "cR",
        "distance"
    ]

    def __init__(self, date):
        """
        date = line.split( "\n" )[ 0 ].split( "\t" ) from BEDPE file
        """
        self.chromL = date[0]
        self.startL = int(date[1])
        self.endL = int(date[2])
        self.strandL = date[8]
        self.chromR = date[3]
        self.startR = int(date[4])
        self.endR = int(date[5])
        self.strandR = date[9]
        if self.chromL == self.chromR:
            # assert left end < right end
            if self.startL + self.endL > self.startR + self.endR:
                self.startL, self.startR = self.startR, self.startL
                self.endL, self.endR = self.endR, self.endL
                self.strandL, self.strandR = self.strandR, self.strandL
            self.cL = int((self.startL + self.endL)/2)
            self.cR = int((self.startR + self.endR)/2)
            self.distance = self.cR - self.cL
        else:
            self.cL, self.cR, self.distance = None, None, None

def parseBedpeFiles(fs, fout, cs, cut, log):
    """
    Get the cis-PETs, organized by chromosomes. Input could be mixed PETs in bedpe.gz or bedpe. Also change read id to numbers to minize memory usage.
    @param fs: bedpe files of replicates, could be .bedpe or .bedpe.gz
    @param fout: output prefix, the name for directory
    @param cs: chroms that wanted, list like ["chr1","chr2"]
    """
    # chroms data
    chroms = {}
    # cis files
    cfs = []
    # distance between PETs mapped to different strands
    i, j, = 0, 0
    for f in fs:
        r = "Parsing PETs from %s, requiring initial distance cutoff > %s" % ( f, cut)
        log.info(r)
        if f.endswith(".gz"):
            of = gzip.open(f, "rb")
        else:
            of = open(f)
        for line in of:
            i += 1
            if i % 100000 == 0:
                mesFlush("%s PETs processed from %s" % (i, f))
            line = line.split("\n")[0].split("\t")
            if "*" in line and "-1" in line:
                continue
            if len(line) < 6:
                continue
            try:
                pet = PET(line)
            except:
                continue
            # cis reads
            if pet.chromL != pet.chromR:
                continue
            # filtering unwanted PETs in chroms
            if len(cs) > 0 and (not (pet.chromL in cs and pet.chromR in cs)):
                continue
            # filtering too close PETs
            if cut > 0 and pet.distance < cut:
                continue
            if pet.chromL not in chroms:
                cf = os.path.join(fout, "%s-%s" % (pet.chromL, pet.chromR) + ".txt")
                chroms[pet.chromL] = {"f": open(cf, "w"), "c": 0}
                cfs.append(cf)
            nline = [chroms[pet.chromL]["c"], pet.cL, pet.cR]
            chroms[pet.chromL]["f"].write("\t".join(map(str, nline)) + "\n")
            chroms[pet.chromL]["c"] += 1
            j += 1
    del (chroms)
    r = "Totally %s PETs from %s, in which %s cis PETs" % (i, ",".join(fs), j)
    log.info(r)
    return cfs

def txt2jd(f):
    data = []
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        data.append([int(line[0]), int(line[1]), int(line[2])])
    data = np.array(data)
    joblib.dump(data, f.replace(".txt", ".jd"))
    os.remove(f)
    return f.replace(".txt", ".jd")

def parseJdFile(file, cut=0):
    key = os.path.split(file)[1].replace(".jd", "")
    key = tuple(key.split("-"))
    mat = joblib.load(file)
    if cut > 0:
        d = mat[:, 2] - mat[:, 1]
        p = np.where(d >= cut)[0]
        mat = mat[p, :]
    return key, mat

def parseIv(iv):
    iv = [
        iv.split(":")[0],
        int(iv.split(":")[1].split("-")[0]),
        int(iv.split(":")[1].split("-")[1])
    ]
    return iv