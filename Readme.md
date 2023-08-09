## MMCT-Loop

---
### Install
MMCT-Loop can only run on UNIX or Mac operating systems.

If you are familiar with conda, you can run these commands to start MMCT-Loop by creating a new conda environment:
```text
conda create -n MMCT-Loop python=3.7
conda activate MMCT-Loop
conda install matplotlib=3.5.3
conda install seaborn=0.12.0
conda install joblib=1.1.1
conda install -c conda-forge hnswlib=0.6.2
conda install scipy=1.7.3
```
Or any other Python(>=3.7) environment that meet the following python package version requirements is ok.
```text
matplotlib >= 3.5.3
seaborn >= 0.12.0
joblib >= 1.1.1
hnswlib >= 0.6.2
scipy >= 1.7.3
```

---
### Usage

#### MMCT-Loop pipeline
All the running steps of MMCT-Loop are written in `run_MMCT-Loop.sh`.

The mandatory parameters of this file are as follows, just configured and run `bash run_MMCT-Loop.sh`.
Initially, we have configured the `run_MMCT-Loop.sh` file to provide loop predictions on the example data.

```text
-out_prefix         The prefix of the output file name.

-in_valid_bedpe     Processed valid pairs in '.bedpe' or '.bedpe.gz' format, multiple files use ',' split.

-cpu                The number of CPU used to run the job, default is 1, set -1 to use all cpus available.
                    (set as the number of chomes is recommanded)

-M                  The maximum number of adjacency edges will be established for each point in building HNSW network
                    stage. Default is 32.

-k                  Query k approximate nearest neighbors for each point in HNSW query stage. Default is 4.

-ef                 The maximum length of the array of cached nearest neighbor relationships in HNSW query stage.
                    Default is 500.
                    
-bs                 The maximum steps from each start point to construct base graph before clustering. Default is 5.

-minPts             The minimum count of PET points per cluster required after clustering. Default is 4.

-dis_arr            A list of required anchor distances between the two candidate loops to be compressed.
                    Default is '(-100000 0 100000 1000000)'. (Negative distance indicate the required length of overlap)

-compress_pets      The file name of PETs used to compress candidate loops.
```

#### Evaluation of candidate loops
This function support users to calculate loops ES, FDR, Poisson_p-value for each candidate loop in `bedpe` format.

Configure follow parameters and just run `bash run_estimate.sh`.
```
-in_valid_bedpe          Processed valid pairs in '.bedpe' format.

-candidate_loops_bedpe   Candidate loops file(in 'bedpe' format) that need to be estimated.
```

---
### Input and Output
Both input valid pairs (PETs) and output loops are of '.bedpe' format, which contains at least six of the following columns:
```
chr1    start1  end1    chr2    start2  end2
```
while the output of test model with '.loop' format in MMCT-Loop evaluation part contains following nine columns:
```
chr1	start1	end1	chr2	start2	end2	ES	FDR	Poisson_p-value
```

### Extension function
| Filename                            | Function                                                                                                                             |
|-------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------|
| ./scripts/bedpe_merge.cpp           | Merge two loop files in '.bedpe' format according to the distance or overlap between the two loops.                                  |
| ./scripts/loop2csv.py               | Filter out loops whose significance is 1 (with high reliability) in the '.loop' format file and convert it to a table in csv format. |
| ./scripts/merge_I_S.cpp             | Filter the effective loops in inter- and intra-ligation loops, extract the main features and merge them.                             |
| ./src/plot.py                       | Plot the estimated fragment size and the finger print for quality control.                                                           |
| ./src/model.py                      | Calculation of ES, FDR and Poisson_p-value.                                                                                          |
| ./src/test_I.py<br/>./src/test_S.py | Calculate the eigenvalues of the respective loops for clusters of inter- and self-ligation loop types.                               |
