import argparse
import sys, os, logging

EPILOG = "Any bug is welcome reported to 204712125@csu.edu.cn"

def getLogger(file_name=None):
    # set up logging, both write log info to console and log file
    logging.basicConfig(
        format='%(asctime)s %(name)-6s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=file_name,
        filemode='a')
    logger = logging.getLogger()
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger

def getHelp():
    description = """
        Intra-chromosomal loops calling for ChIA-PET,HiChIP and high-resolution Hi-C data.
        """
    parser = argparse.ArgumentParser(description=description, epilog=EPILOG)
    parser.add_argument(
        "-f",
        dest="fnIn",
        required = True,
        type=str,
        help=
        "Input file is mapped PETs, tab-delimited BEDPE format in gzip. Replicates could be input as -f A.bedpe.gz,B.bedpe.gz,C.bedpe.gz. Loops will be called with pooled data."
    )
    parser.add_argument(
        "-o",
        dest="fnOut",
        required=True,
        type=str,
        help="outputPrefix."
    )
    parser.add_argument(
        "-p",
        dest="fnPeak",
        required=False,
        type=str,
        help="Peak file."
    )
    parser.add_argument(
        "-M",
        dest="M",
        default=0,
        required=True,
        help=
        "M"
    )
    parser.add_argument(
        "-ef",
        dest="ef",
        default=0,
        required=True,
        help=
        "ef"
    )
    parser.add_argument(
        "-k",
        dest="k",
        default=0,
        required=True,
        help=
        "k"
    )
    parser.add_argument(
        "-le",
        dest="le",
        default=0,
        required=True,
        help=
        "leiden model: 0-ModularityVertexPartition; 1-CPMVertexPartition; 2-RBConfigurationVertexPartition"
    )
    parser.add_argument(
        "-r",
        dest="resolution",
        default=0,
        required=False,
        help=
        "resolution"
    )
    parser.add_argument(
        "-minPts",
        dest="minPts",
        default=0,
        help=
        "Points required in a cluster"
    )
    parser.add_argument(
        "-cpu",
        dest="cpu",
        required=False,
        default=-1,
        type=int,
        help=
        "CPU number used to run the job, default is 1,set -1 to use all cpus available. Too many CPU could cause memory error."
    )
    parser.add_argument(
        "-c",
        dest="chroms",
        required=False,
        default="",
        type=str,
        help=
        "Whether to process limited chroms, specify it as chr1,chr2,chr3, default is not. Use this to filter reads in like chr22_KI270876v1_alt"
    )
    parser.add_argument(
        "-s",
        dest="tmp",
        required=False,
        default=True,
        action="store_true",
        help=
        "Whether or not to save temp files for each chromosomes during processing. Set this flag for following calling differentially enriched loops or converting PETs to washU track or hic file load into juicebox. Default is not."
    )
    parser.add_argument(
        "-cut",
        dest="cut",
        required=False,
        default=0,
        type=int,
        help=
        "Initial distance cutoff to filter PETs, default is 0, only used for debuging."
    )
    parser.add_argument(
        "-plot",
        dest="plot",
        required=False,
        default=False,
        action="store_true",
        help=
        "Whether to plot estimated inter-ligation and self-ligation PETs distance distrbution, default is not."
    )

    op = parser.parse_args()
    return op

def mesFlush(mes):
    sys.stdout.write("\r%s" % mes)
    sys.stdout.flush()

def callSys(cmds, log = None):
    for cmd in cmds:
        try:
            log.info(cmd)
        except:
            print(cmd)
        try:
            os.system(cmd)
        except:
            try:
                log.error(cmd)
            except:
                print("ERROR !!! found in ", cmd)