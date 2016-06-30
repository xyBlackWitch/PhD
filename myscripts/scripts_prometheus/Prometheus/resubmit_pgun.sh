#!/bin/bash

while read i; do
	qsub -t $i job_prod_pgun.sge XiMinus 2000 3312 3 1
done < /hera/panda/jpuetz/data/failed_id.txt
