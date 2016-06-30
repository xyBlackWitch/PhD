#!/bin/bash

while read i; do
	qsub -t $i job_prod.sge DPM_run23 1000 DPM 4.6 
done < /hera/panda/jpuetz/data/failed_id_run23.txt
