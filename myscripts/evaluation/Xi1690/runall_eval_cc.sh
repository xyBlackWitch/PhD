#!/bin/bash

prefix=""
save=true

if test "$1" != ""; then
	prefix=$1
fi

if test "$2" != ""; then
	save=false
fi



echo evaluation of final states
root -l -q eval_final_states_antiparticle.C"(\"$prefix\",$save)"

echo evaluation of Lambda0 and AntiLambda0
root -l -q eval_lambda0.C"(\"$prefix\",$save)"

echo evaluation of Xi or anti-Xi
root -l -q eval_Xi.C"(\"$prefix\",$save)"

echo evaluation of Xi1690 or anti-Xi1690
root -l -q eval_Xi1690.C"(\"$prefix\",$save)"

echo evaluation of XiSys
root -l -q eval_XiSys.C"(\"$prefix\",$save)"
