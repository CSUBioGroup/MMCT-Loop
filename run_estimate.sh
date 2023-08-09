#!/bin/bash

in_dir=example_workspace

in_valid_bedpe=${in_dir}/example_valid.bedpe    #Processed valid pairs in '.bedpe' format.
candidate_loops_bedpe=${in_dir}/example.bedpe   #Candidate loops file(in 'bedpe' format) that need to be estimated.

out_prefix=example

cpu=16

python ./src/test_I.py -f ./${in_valid_bedpe} -o ${out_prefix} -cpu ${cpu} -tf ${candidate_loops_bedpe}
python ./src/test_S.py -f ./${in_valid_bedpe} -o ${out_prefix} -cpu ${cpu} -tf ${candidate_loops_bedpe}

if [ ! -e "./scripts/merge_I_S" ]; then
   cd scripts
   g++ -std=c++11 merge_I_S.cpp -o merge_I_S
   cd ..
fi

./scripts/merge_I_S ${out_prefix}_intra_ligation.loop ${out_prefix}_self_ligation.loop > ${out_prefix}.loop
rm ./*_ligation.loop
