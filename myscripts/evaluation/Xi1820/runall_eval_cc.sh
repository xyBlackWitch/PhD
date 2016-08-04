#!/bin/bash

prefix=""
path=""
save=true

if test "$1" != ""; then
	prefix=$1
fi

if test "$2" != ""; then
	path=$2
fi

if test "$3" != ""; then
	save=false
fi



#echo evaluation of final states
#root -l -q eval_final_states_antiparticle.C"(\"$prefix\",\"$path\",$save)"

#echo evaluation of Lambda0 and AntiLambda0
#root -l -q eval_lambda0.C"(\"$prefix\",\"$path\",$save)"

#echo evaluation of Xi or anti-Xi
#root -l -q eval_Xi.C"(\"$prefix\",\"$path\",$save)"

#echo evaluation of Xi1820 or anti-Xi1820
#root -l  -q eval_Xi1820.C"(\"$prefix\",\"$path\",$save)"

echo evaluation of XiSys
root -l  eval_XiSys.C"(\"$prefix\",\"$path\",$save)"
