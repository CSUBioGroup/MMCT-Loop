## LoopCaller

---
### Install
LoopCaller can only run on UNIX or Mac operating systems.

If you are familiar with conda, you can run these commands to start LoopCaller by creating a new conda environment:
```
conda create -n LoopCaller python=3.7
conda activate LoopCaller
conda install matplotlib=3.5.3
conda install seaborn=0.12.0
conda install joblib=1.1.1
conda install hnswlib=0.6.2
conda install scipy=1.7.3
```
Or any other environment that meet the following python package version requirements is ok.
```
matplotlib >= 3.5.3
seaborn >= 0.12.0
joblib >= 1.1.1
hnswlib >= 0.6.2
scipy >= 1.7.3
```

---
### Usage

#### LoopCaller pipeline
All the running steps of LoopCaller are written in `run_LoopCaller.sh`.

The mandatory parameters of the `run_LoopCaller.sh` are as follows:
```
-out_prefix         The prefix of the output file name.
-in_valid_bedpe     Processed valid pairs in '.bedpe' or '.bedpe.gz' format, multiple files use ',' split.
-cpu                The number of CPU used to run the job, default is 1, set -1 to use all cpus available. (set as the number of chomes is recommanded)
-M                  The maximum number of adjacency edges will be established for each point in building HNSW network stage. Default is 32.
-k                  Query k approximate nearest neighbors for each point in HNSW query stage. Default is 4.
-ef                 The maximum length of the array of cached nearest neighbor relationships in HNSW query stage. Default is 500.
-bs                 The maximum steps from each start point to construct base graph before clustering.
-minPts             The minimum count of PET points per cluster required after clustering. Default is 4.
-target             The candidate loops genarated by peak-based tools.(in '.bedpe' format)
-mergeDis           A list of anchor distances between two candidate loops, used to determine whether to merge.
                    Default is '(-10000 0 10000 100000)'. (Negative distance indicate the required length of overlap)
```

#### Evaluation of candidate loops
Given valid pairs and candidate loops file in '.bedpe' format, modify `run_test.sh` and run the following script to calculate ES,FDR,P-value for each loop.
```
bash run_test.sh ${filename_valid_pairs} ${filename_loops}
```
---
### Input and Output
Both input valid pairs (PETs) and output loops are of '.bedpe' format, which contains at least six of the following columns:
```
chr1    start1  end1    chr2    start2  end2
```
while the output of test model with '.loop' format in LoopCaller contains following nine columns:
```
chr1	start1	end1	chr2	start2	end2	ES	FDR	Poisson_p-value
```
