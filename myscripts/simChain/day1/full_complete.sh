#! /bin/bash


nevts=10
prefix=""
mom=4.6
path=""
dec=""

if test "$1"!=""; then
	prefix=$1
fi

if test "$2"!=""; then
	nevts=$2
fi

if test "$3"!=""; then
	dec=$3
fi


if test "$4"!=""; then
	mom=$4
fi

if test "$5" != ""; then
	prefix=$5"/"$prefix
fi

date

echo sim
root -l -q ~/panda/myscripts/simChain/day1/sim_complete.C"($nevts,\"$dec\",$mom,\"$prefix\")" >> $prefix"_"sim.log 2>&1 

echo digi
root -l -q ~/panda/myscripts/simChain/day1/digi_complete.C"(\"$prefix\",$nevts)" >> $prefix"_"digi.log 2>&1 

rename .root _complete.root $prefix*
mv $prefix"_"par"_"complete.root $prefix"_"simparams.root

echo reco
root -l -q ~/panda/myscripts/simChain/SimMacros/recoideal_complete.C"(\"$prefix\")" >> $prefix"_"reco.log 2>&1 

rename _complete.root .root $prefix*
mv $prefix"_"simparams.root $prefix"_"par.root

echo pid
root -l -q ~/panda/myscripts/simChain/day1/pid_complete.C"(\"$prefix\",$nevts)" >> $prefix"_"pid.log 2>&1
 
rename .root _complete.root $prefix*
mv $prefix"_"par"_"complete.root $prefix"_"simparams.root

echo ana
root  -l -q /home/ikp1/puetz/panda/myscripts/analysis/AnalysisTaskRunAntiXi1820.C"($mom,$nevts,\"$prefix\")" >> $prefix"_"ana.log 2>&1

echo END

date