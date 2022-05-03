#!/bin/bash
valid_pairs=$1
loops=$2
prefix=prefix_example

echo "startCalPvalue('${loops}')" >> ../LoopCaller/test_I.py
python ../LoopCaller/test_I.py -f ${valid_pairs} -o ${prefix} -M 32 -ef 500 -k 4 -le 1 -r 0.05 -minPts 4 -cpu 16
mv ${prefix}_I.loop ${prefix}_LoopCaller_estimated_I.bedpe
head -n -1 ../LoopCaller/test_I.py > test_I_tmp.py
mv test_I_tmp.py ../LoopCaller/test_I.py
rm -rf ${prefix}/

echo "startCalPvalue('${prefix}_LoopCaller_d${dis}.bedpe')" >> ../LoopCaller/test_S.py
python ../LoopCaller/test_S.py -f ${valid_pairs} -o ${prefix} -M 32 -ef 500 -k 4 -le 1 -r 0.05 -minPts 4 -cpu 16
mv ${prefix}_S.loop ${prefix}_LoopCaller_estimated_S.bedpe
head -n -1 ../LoopCaller/test_S.py > test_S_tmp.py
mv test_S_tmp.py ../LoopCaller/test_S.py
rm -rf ${prefix}/

./scripts/merge_I_S ${prefix}_LoopCaller_estimated_I.bedpe ${prefix}_LoopCaller_estimated_S.bedpe > ${prefix}_tested.loop
rm *_estimated_*.bedpe

