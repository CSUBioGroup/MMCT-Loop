import numpy as np


from configs import *
import pylab

def plotFragSize(ds, frag, log=1, prefix="test"):
    """
    Plot the estimated fragment size.
    """
    # print ("@plotFragSize")
    ds = np.abs(np.array(ds))
    ds = ds[~np.isnan(ds)]
    ds = ds[ds > 0]
    if log:
        ds = np.log2(ds)
    fig, ax = pylab.subplots()
    sns.kdeplot(ds,
                ax=ax,
                shade=True,
                label="distance between PETs",
                color=colors[0])
    ax.axvline(np.log2(frag),
               label="fragment size:%s bp" % (frag),
               color=colors[1])
    ax.set_xlabel("Distance between different strand PETs (log2(bp))")
    ax.set_ylabel("Density")
    ax.legend(loc="best")
    pylab.savefig("%s.pdf" % prefix)


def plotInterSelfCutFrag(di, ds, cut, frag, log=1, prefix="test"):
    """
    Plot the distance cutoff of self-ligation and inter-ligation reads.
    """
    # print ("@plotIntSelCutFrag")
    di = np.abs(np.array(di))
    ds = np.abs(np.array(ds))
    di = di[~np.isnan(di)]
    ds = ds[~np.isnan(ds)]
    di = di[di > 0]
    ds = ds[ds > 0]
    if log:
        di = np.log2(di)
        ds = np.log2(ds)
    fig, ax = pylab.subplots()
    sns.kdeplot(di,
                ax=ax,
                shade=True,
                label="inter-ligation PETs:%s" % len(di),
                color=colors[0])
    sns.kdeplot(ds,
                ax=ax,
                shade=True,
                label="self-ligation PETs:%s" % len(ds),
                color=colors[1])
    ax.axvline(np.log2(cut),
               label="distance cutoff:%.2f kb" % (cut / 1000.0),
               color=colors[2])
    # ax.axvline(
    #    np.log2(frag), label="fragment size:%s bp" % (frag), color=colors[3])
    leg = ax.legend(loc="best", shadow=True, fancybox=True)
    ax.set_xlabel("Distance between PETs (log2(bp))")
    ax.set_ylabel("Density")
    # ax.set_title("ratio of inter/self PETs:%.3f"%(len(di)/float(len(ds))))
    pylab.savefig("%s.pdf" % prefix)


def plotFingerPrint(data, prefix="test"):
    """
    Plot the finger print for quality control. 
    data: pd.Dataframe, raw are raw of summaried bins, columns are different datasets name
    """
    # print ("@plotFingerPrint")
    fig, ax = pylab.subplots()
    x = data.index
    for c in data.columns:
        ax.plot(x, data[c], label=c)
    ax.legend()
    ax.set_xlabel("bins of contact matrix rank from low to high")
    ax.set_ylabel("PETs ratio")
    pylab.savefig("%s_fingerprint.pdf" % prefix)
