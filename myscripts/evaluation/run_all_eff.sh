old=""
new=""

if test "$1"!=""; then
	old=$1
fi

if test "$2"!=""; then
	new=$2
fi

echo Lambda0
run -r -q ~/panda/myscripts/evaluation/efficiency_after_cuts.C"(\"$1\",\"$2\")"
echo AntiLambda0
run -r -q ~/panda/myscripts/evaluation/efficiency_after_cuts.C"(\"$1\",\"$2\", \"ntpAntiLambda0\", \"antiLambda0\")"
echo XiPlus
run -r -q ~/panda/myscripts/evaluation/efficiency_after_cuts.C"(\"$1\",\"$2\", \"ntpXiPlus\", \"xiplus\", false)"
echo Xi1820
run -r -q ~/panda/myscripts/evaluation/efficiency_after_cuts_Xi1820.C"(\"$1\",\"$2\")"
echo XiSys
run -r -q ~/panda/myscripts/evaluation/efficiency_after_cuts_Full_Tree.C"(\"$1\",\"$2\")"