#!/bin/bash

path=""
save=1

if test "$1" != ""; then
	path=$1
fi


echo evaluation of final states
root -l -q eval_final_states.C"(\"$path\", $save)"

echo evaluation of Lambda0 and AntiLambda0
root -l -q eval_lambda0.C"(\"$path\", $save)"

echo evaluation of Xi or anti-Xi
root -l -q eval_Xi.C"(\"$path\", $save)"

echo evaluation of Xi(1820) or anti-Xi(1820)
root -l -q eval_Xi1820.C"(\"$path\", $save)"