#!/bin/bash


work_dir=example_workspace
prefix=example
validBedpe=example_valid.bedpe
target=example_target.bedpe
k=4
minPts=4
cpu=16
mergeDis=(-10000 0 10000 100000)


### PETs cluster

echo ${mergeDis[@]}

awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ./${work_dir}/${target} > ${prefix}_temp.bedpe

python ./LoopCaller/main.py -f ./${work_dir}/${validBedpe} -o ${prefix} -M 32 -ef 500 -k ${k} -le 1 -r 0.05 -minPts ${minPts} -cpu ${cpu}
rm -rf ${prefix}/

cluster_results=${prefix}_cluster.loop

python ./scripts/loopFormat.py ${cluster_results} >> ${prefix}_temp.bedpe

rm ${cluster_results}

if [ ! -e "./scripts/bedpemerge" ]; then
 cd scripts
 g++ bedpemerge.cpp -o bedpemerge
 cd ..
fi

### merge loops

if [ -d "output_LoopCaller" ]; then
 rm -rf output_LoopCaller
fi

mkdir output_LoopCaller

for dis in ${mergeDis[@]}
do
 ./scripts/bedpemerge ${prefix}_temp.bedpe -d ${dis} > ./output_LoopCaller/${prefix}_LoopCaller_d${dis}.bedpe
done
rm ${prefix}_temp.bedpe

### test loops

if [ ! -e "./scripts/merge_I_S" ]; then
 cd scripts
 g++ merge_I_S.cpp -o merge_I_S
 cd ..
fi

cd output_LoopCaller/

echo "cd output_LoopCaller/"

for dis in ${mergeDis[@]}
do
 echo "startCalPvalue('${prefix}_LoopCaller_d${dis}.bedpe')" >> ../LoopCaller/test_I.py
 python ../LoopCaller/test_I.py -f ../${work_dir}/${validBedpe} -o ${prefix} -M 32 -ef 500 -k ${k} -le 1 -r 0.05 -minPts ${minPts} -cpu ${cpu}
 mv ${prefix}_I.loop ${prefix}_LoopCaller_d${dis}_estimated_I.bedpe
 head -n -1 ../LoopCaller/test_I.py > test_I_tmp.py
 mv test_I_tmp.py ../LoopCaller/test_I.py
 rm -rf ${prefix}/

 echo "startCalPvalue('${prefix}_LoopCaller_d${dis}.bedpe')" >> ../LoopCaller/test_S.py
 python ../LoopCaller/test_S.py -f ../${work_dir}/${validBedpe} -o ${prefix} -M 32 -ef 500 -k ${k} -le 1 -r 0.05 -minPts ${minPts} -cpu ${cpu}
 mv ${prefix}_S.loop ${prefix}_LoopCaller_d${dis}_estimated_S.bedpe
 head -n -1 ../LoopCaller/test_S.py > test_S_tmp.py
 mv test_S_tmp.py ../LoopCaller/test_S.py
 rm -rf ${prefix}/

 ../scripts/merge_I_S ${prefix}_LoopCaller_d${dis}_estimated_I.bedpe ${prefix}_LoopCaller_d${dis}_estimated_S.bedpe > ${prefix}_LoopCaller_d${dis}_estimated.bedpe
 rm *_estimated_*.bedpe

done

cd ..


