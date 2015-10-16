#!/bin/sh
outdir="/private/puetz/mysimulations/analysis/pbarp_Xiplus_Ximinus/run1"
logdir=$outdir"/logs"

nevts=$1
prefix=$2


echo ana
root -l -q ../analysis/AnalysisTaskXi.C"($nevts,\"$outdir/$prefix\")" >> $logdir/ana_$prefix.log 2>&1

