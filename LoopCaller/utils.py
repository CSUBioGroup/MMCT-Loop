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
        Description: Intra-chromosomal loops calling for ChIA-PET,HiChIP and high-resolution Hi-C data.
        """
    parser = argparse.ArgumentParser(description=description, epilog=EPILOG)
    parser.add_argument(
        "-f",
        dest="in_valid_bedpe",
        required = True,
        help=
        "Input file is mapped PETs, tab-delimited BEDPE format in gzip. Replicates could be input as -f A.bedpe.gz,B.bedpe.gz,C.bedpe.gz. Loops will be called with pooled data."
    )
    parser.add_argument(
        "-o",
        dest="out_prefix",
        required=True,
        help="The prefix of the output file name."
    )
    parser.add_argument(
        "-M",
        dest="hnsw_M",
        default=32,
        required=True,
        help=
        "The maximum number of adjacency edges will be established for each point in building HNSW network stage."
    )
    parser.add_argument(
        "-ef",
        dest="hnsw_ef",
        default=500,
        required=True,
        help=
        "The maximum length of the array of cached nearest neighbor relationships in HNSW query stage."
    )
    parser.add_argument(
        "-k",
        dest="hnsw_k",
        default=4,
        required=True,
        help=
        "Query k approximate nearest neighbors for each point in HNSW query stage."
    )
    parser.add_argument(
        "-bs",
        dest="base_step",
        default=3,
        required=False,
        help=
        "The maximum steps from each start point to construct base graph before clustering."
    )
    parser.add_argument(
        "-minPts",
        dest="minPts",
        default=4,
        required=False,
        help=
        "The minimum count of PET points per cluster required after clustering. Default is 4."
    )
    parser.add_argument(
        "-cpu",
        dest="cpu",
        required=False,
        default=-1,
        help=
        "The number of CPU used to run the job, default is 1, set -1 to use all cpus available. Too many CPUs may cause memory error."
    )
    parser.add_argument(
        "-chr",
        dest="chromes",
        required=False,
        default="",
        help=
        "Whether to process limited chromes, specify it as chr1,chr2,chr3, default is not. Use this to filter reads in like chr22_KI270876v1_alt"
    )
    parser.add_argument(
        "-s",
        dest="save",
        required=False,
        default=False,
        help=
        "Whether or not to save temp files for each chromosomes during processing. Set this flag for following calling differentially enriched loops or converting PETs to WashU Track or hic file load into Juice box. Default is not."
    )
    parser.add_argument(
        "-cut",
        dest="distance_cut",
        required=False,
        default=0,
        help=
        "Initial distance cutoff to filter PETs, default is 0, only used for debug."
    )
    parser.add_argument(
        "-plot",
        dest="plot_dis",
        required=False,
        default=False,
        help=
        "Whether to plot estimated inter-ligation and self-ligation PETs distance distribution, default is not."
    )

    op = parser.parse_args()
    return op

def getEstimatingArg():
    description = """
            Description: LoopCaller estimating inter-ligation loops.
            """
    parser = argparse.ArgumentParser(description=description, epilog=EPILOG)
    parser.add_argument(
        "-f",
        dest="in_valid_bedpe",
        required=True,
        help=
        "Input file is mapped PETs, tab-delimited BEDPE format in gzip. Replicates could be input as -f A.bedpe.gz,B.bedpe.gz,C.bedpe.gz. Loops will be called with pooled data."
    )
    parser.add_argument(
        "-o",
        dest="out_prefix",
        required=True,
        help="The prefix of the output file name."
    )
    parser.add_argument(
        "-minPts",
        dest="minPts",
        default=4,
        required=False,
        help=
        "The minimum count of PET points per cluster required after clustering. Default is 4."
    )
    parser.add_argument(
        "-cpu",
        dest="cpu",
        required=False,
        default=-1,
        help=
        "The number of CPU used to run the job, default is 1, set -1 to use all cpus available. Too many CPUs may cause memory error."
    )
    parser.add_argument(
        "-chr",
        dest="chromes",
        required=False,
        default="",
        help=
        "Whether to process limited chromes, specify it as chr1,chr2,chr3, default is not. Use this to filter reads in like chr22_KI270876v1_alt"
    )
    parser.add_argument(
        "-cut",
        dest="distance_cut",
        required=False,
        default=0,
        help=
        "Initial distance cutoff to filter PETs, default is 0, only used for debug."
    )
    parser.add_argument(
        "-tf",
        dest="test_file",
        required=False,
        default=0,
        help=
        "The file name of loops in 'bedpe' format to estimate."
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