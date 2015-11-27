#!/bin/bash

path=""


if test "$1" != ""; then
	path=$1
fi


echo evaluation of final states
root -l -q eval_final_states_antiparticle.C"(\"$path\")"

echo evaluation of Lambda0 and AntiLambda0
root -l -q eval_lambda0.C"(\"$path\")"

echo evaluation of Xi or anti-Xi
root -l -q eval_Xi.C"(\"$path\")"

echo evaluation of Xi1820 or anti-Xi1820
root -l -q eval_Xi1820.C"(\"$path\")"

echo evaluation of XiSys
root -l -q eval_XiSys.C"(\"$path\")"
