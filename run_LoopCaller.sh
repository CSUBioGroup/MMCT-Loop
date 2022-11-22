#!/bin/bash

#------ The mandatory parameters of LoopCaller ------

in_valid_bedpe=example_valid.bedpe 	# Processed valid pairs in '.bedpe' format.
out_prefix=example			            # The prefix of the output file name.
cpu=4 					                    # The number of CPU used to run the job, default is 1, set -1 to use all cpus available. Too many CPUs may cause memory error.

# HNSW parameters
M=32					# The maximum number of adjacency edges will be established for each point in building HNSW network stage. Default is 32.
k=4           # Query k approximate nearest neighbors for each point in HNSW query stage. Default is 4.
ef=500			  # The maximum length of the array of cached nearest neighbor relationships in HNSW query stage. Default is 500.

# Cluster parameters
bs=3					# The maximum steps from each start point to construct base graph before clustering. Default is 3.
minPts=4			# The minimum count of PET points per cluster required after clustering. Default is 4.


#---------------- Start LoopCaller ----------------
run_steps=(1 2 3 4)

results_dir=${out_prefix}_LoopCaller_results
if [ ! -d ${results_dir} ];then
  mkdir ${results_dir}
fi

## Step 1 : Call loops
if [[ "${run_steps[*]}" =~ "1" ]];then
  python ./src/LoopCaller.py -f ${in_valid_bedpe} -o ${out_prefix} -cpu ${cpu} -M ${M} -k ${k} -ef ${ef} -bs ${bs} -minPts ${minPts}
  mv ${out_prefix}.loop ${results_dir}/${out_prefix}.loop
fi

## Step 2 : Loops format
if [[ "${run_steps[*]}" =~ "2" ]];then
  loop_format_bedpe=${out_prefix}.bedpe
  if [ -f ${results_dir}/${out_prefix}.loop ];then
    python ./scripts/loopFormat.py ${results_dir}/${out_prefix}.loop > ${results_dir}/${loop_format_bedpe}
  fi
fi

## Step 3 : Compress loops
dis_arr=(-100000 0 100000 1000000)            	# A list of required anchor distances between the two loops to be compressed.
                                             	# Default is '(-100000 0 100000 1000000)'. (Negative distance indicate the required length of overlap)
if [[ "${run_steps[*]}" =~ "3" ]];then
  compress_pets=example_compress_pets.bedpe      # The PETs used to compress loops

  if [ ! -e "./scripts/bedpe_compress" ];then
    cd scripts
    g++ -std=c++11 bedpe_compress.cpp -o bedpe_compress
    cd ..
  fi

  loop_format_bedpe=./${results_dir}/${out_prefix}.bedpe
  if [ -f ${compress_pets} ] && [ -f ${loop_format_bedpe} ];then
    if [ ! -d ${results_dir} ];then
      mkdir ${results_dir}
    fi
    cat ${compress_pets} > pets_set.bedpe
    cat ${loop_format_bedpe} >> pets_set.bedpe
    for dis in ${dis_arr[*]}
    do
      ./scripts/bedpe_compress pets_set.bedpe -d ${dis} > ./${results_dir}/${out_prefix}_dis_${dis}.bedpe
    done
    rm pets_set.bedpe
  fi
fi

## Step 4 : Estimate significance of compressed loops
if [[ "${run_steps[*]}" =~ "4" ]];then
  if [ ! -e "./scripts/compress_I_S" ]; then
   cd scripts
   g++ -std=c++11 compress_I_S.cpp -o compress_I_S
   cd ..
  fi

  for dis in ${dis_arr[*]}
  do
    python ./src/test_I.py -f ./${in_valid_bedpe} -o ${out_prefix} -minPts ${minPts} -cpu ${cpu} -tf ./${results_dir}/${out_prefix}_dis_${dis}.bedpe
    python ./src/test_S.py -f ./${in_valid_bedpe} -o ${out_prefix} -minPts ${minPts} -cpu ${cpu} -tf ./${results_dir}/${out_prefix}_dis_${dis}.bedpe
    ./scripts/compress_I_S ${out_prefix}_intra_ligation.loop ${out_prefix}_self_ligation.loop > ${results_dir}/${out_prefix}_dis_${dis}.loop
    rm ./*_ligation.loop
  done
  cd ..
fi
