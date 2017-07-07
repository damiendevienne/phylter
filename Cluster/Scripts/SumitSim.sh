#!/bin/sh 

for i in `seq 10001 11700`; #11700
do
        qsub -q q1day -N script.$i -o /panhome/comte/output -e /panhome/comte/output <<EOF 
/panhome/comte/scriptR-5.sh $i
EOF
done 

