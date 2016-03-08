prefix=""
path=""
nevts=100
mom=6.2


if test "$1"!=""; then
	prefix=$1
fi

outname=$prefix

if test "$2"!=""; then
	nevts=$2
fi

if test "$3"!=""; then
	mom=$3
fi

if test "$4" != ""; then
	path=$4
	outname=$path"/"$prefix
fi

echo sim
root  -l -q sim_complete.C"($nevts,$mom,\"$outname\")" >> sim$prefix.log 2>&1 
echo digi
root  -l -q digi_complete.C"(\"$outname\")" >> digi$prefix.log 2>&1 
echo reco
root  -l -q recoideal_GenFit2_complete.C"(\"$outname\")" >> reco$prefix.log 2>&1
echo pid
root  -l -q pid_complete_GenFit2.C"(\"$outname\")" >> pid$prefix.log 2>&1 
echo END

