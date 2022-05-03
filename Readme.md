## LoopCaller


### Install
If you are familar with conda, LoopCaller could be installed very easily with following after clone this project.
```
git clone https://github.com/CSUBioGroup/LoopCaller
```

Before you install a new environment, please make sure your `~/.condarc` file contains the following channels or their mirror
links:
```
channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
  - conda-forge
  - bioconda
  - defaults
```
Then, you can cd in LoopCaller, import `LoopCaller_env.yaml` and create a new environment.
```
cd LoopCaller/
conda env create -n LoopCaller -f LoopCaller_env.yaml
```
---
### Usage

#### LoopCaller pipeline

All the running steps of LoopCaller are written in `run_LoopCaller.sh`.
The following parameters are required:
```
-prefix         The name of outdir and the prefix of output files.

-validBedpe     Processed valid pairs in '.bedpe' format.

-M              The number of adjacencies established by each PET point at most in HNSW. Default is 32.

-k              Query for k closest points for each PET point in HNSW. Default is 4.

-minPts         The minimum count of PET points per cluster required after clustering. Default is 4.

-le             The clustering model of Leiden algorithm. Use 'ModularityVertexPartition' model by 0, 
                or 'CPMVertexPartition' model by 1. Default is 0.

-cpu            The number of cpu cores used. Default is 16.

-target         The candidate loops genarated by peak-based tools.(in '.bedpe' format)

-mergeDis       A list of anchor distances between two candidate loops, used to determine whether to merge.
                Default is '(-10000 0 10000 100000)'. (Negative distance indicate the required length of overlap) 
```

#### Evaluation of candidate loops
Given valid pairs and candidate loops file in '.bedpe' format, modify `run_test.sh` and run the following script to calculate ES,FDR,P-value for each loop.
```
bash run_test.sh ${valid_pairs} ${loops}
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
