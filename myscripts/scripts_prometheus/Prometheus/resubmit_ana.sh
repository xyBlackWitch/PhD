#!/bin/bash

while read i; do
	qsub -t $i job_ana_Xi1820.sge XiMinus1820_3half_branching_run3 1000 4.6
done < /hera/panda/jpuetz/data/XiMinus1820_branching/spin_3half/failed_id_new3.txt
