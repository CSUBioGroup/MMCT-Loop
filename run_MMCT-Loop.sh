#!/bin/bash

#------ The mandatory parameters of MMCT-Loop ------
in_dir=example_workspace

in_valid_bedpe=${in_dir}/example_valid.bedpe # Processed valid pairs in '.bedpe' format.
out_prefix=example 		# The prefix of the output file name.
cpu=4 			# The number of CPU used to run the job, default is 1, set -1 to use all cpus available. Too many CPUs may cause memory error.

# HNSW parameters
M=32 			# The maximum number of adjacency edges will be established for each point in building HNSW network stage. Default is 32.
k=4 			# Query k approximate nearest neighbors for each point in HNSW query stage. Default is 4.
ef=500 			# The maximum length of the array of cached nearest neighbor relationships in HNSW query stage. Default is 500.

# Cluster parameters
bs=3 			# The maximum steps from each start point to construct base graph before clustering. Default is 3.
minPts=4 			# The minimum count of PET points per cluster required after clustering. Default is 4.


#---------------- Start MMCT-Loop ----------------
run_steps=(1 2) 	# Default used steps 1 and 2. You can also add step 3 and 4 to merge loops and estimate significance.

results_dir=${out_prefix}_MMCT-Loop_results
if [ ! -d ${results_dir} ];then
  mkdir ${results_dir}
fi

## Step 1 : Call loops
if [[ "${run_steps[*]}" =~ "1" ]];then
  python ./src/MMCT-Loop.py -f ${in_valid_bedpe} -o ${out_prefix} -cpu ${cpu} -M ${M} -k ${k} -ef ${ef} -bs ${bs} -minPts ${minPts}
  mv ${out_prefix}.loop ${results_dir}/${out_prefix}.loop
fi

## Step 2 : Loops format
if [[ "${run_steps[*]}" =~ "2" ]];then
  loop_format_bedpe=${out_prefix}.bedpe
  if [ -f ${results_dir}/${out_prefix}.loop ];then
    python ./scripts/loopFormat.py -b ${results_dir}/${out_prefix}.loop > ${results_dir}/${loop_format_bedpe}
    python ./scripts/loopFormat.py -l ${results_dir}/${out_prefix}.loop > ${results_dir}/${out_prefix}.loop2
    mv ${results_dir}/${out_prefix}.loop2 ${results_dir}/${out_prefix}.loop
  fi
fi

## Step 3 : Loops merge
dis_arr=(-100000 0 100000 1000000)  # A list of required anchor distances between the two loops to be merged.
                                    # Default is '(-100000 0 100000 1000000)'. (Negative distance indicate the required length of overlap)
if [[ "${run_steps[*]}" =~ "3" ]];then
  merge_pets=${in_dir}/example_merge_pets.bedpe      # The PETs used to merge loops

  if [ ! -e "./scripts/bedpe_merge" ];then
    cd scripts
    g++ -std=c++11 bedpe_merge.cpp -o bedpe_merge
    cd ..
  fi

  loop_format_bedpe=./${results_dir}/${out_prefix}.bedpe
  if [ -f ${merge_pets} ] && [ -f ${loop_format_bedpe} ];then
    if [ ! -d ${results_dir} ];then
      mkdir ${results_dir}
    fi
    cat ${merge_pets} > pets_set.bedpe
    cat ${loop_format_bedpe} >> pets_set.bedpe
    for dis in ${dis_arr[*]}
    do
      ./scripts/bedpe_merge pets_set.bedpe -d ${dis} > ./${results_dir}/${out_prefix}_dis_${dis}.bedpe
    done
    rm pets_set.bedpe
  fi
fi

## Step 4 : Estimate significance of merged loops
if [[ "${run_steps[*]}" =~ "4" ]];then
  if [ ! -e "./scripts/merge_I_S" ]; then
   cd scripts
   g++ -std=c++11 merge_I_S.cpp -o merge_I_S
   cd ..
  fi

  for dis in ${dis_arr[*]}
  do
    python ./src/test_I.py -f ./${in_valid_bedpe} -o ${out_prefix} -cpu ${cpu} -tf ./${results_dir}/${out_prefix}_dis_${dis}.bedpe
    python ./src/test_S.py -f ./${in_valid_bedpe} -o ${out_prefix} -cpu ${cpu} -tf ./${results_dir}/${out_prefix}_dis_${dis}.bedpe
    ./scripts/merge_I_S ${out_prefix}_intra_ligation.loop ${out_prefix}_self_ligation.loop > ${results_dir}/${out_prefix}_dis_${dis}.loop
    rm ./*_ligation.loop
  done
  cd ..
fi
