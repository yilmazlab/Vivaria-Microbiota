#!/bin/bash

for i in $(seq 0 149)
do
        qsub -V -cwd -N MP_${i} -q sc03.q submit2cluster/submit2cluster_${i}.sh
done