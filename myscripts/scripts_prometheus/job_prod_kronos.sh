#! /bin/bash

# start pandaroot software
. "/lustre/nyx/panda/jpuetz/pandaroot/trunk/build/config.sh"

#SBATCH --get-user-env


# Task name
#SBATCH -J XiPlus1690_run1 

# Run time limit
#SBATCH --time=03:00:00

# Working directory on shared storage
#SBATCH -D /lustre/nyx/panda/jpuetz/data

# Script directory (If different from working directory!)
scripts="/lustre/nyx/panda/jpuetz/myscripts"

# Default parameters
prefix=9999
nEvts=1000
dec="XiMinus_1820_lambda0_K.dec"
mom=12.3485
res="pbarpSystem"
save=""

sig=1


tmpdir="/tmp/"$USER"_"$SLURM_JOB_ID"/"
mkdir $tmpdir

_target="/lustre/nyx/panda/jpuetz/data/"

# Parameters set by user
if test "$1" != ""; then
  prefix=$tmpdir$1
fi

if test "$2" != ""; then
  nEvts=$2
fi

if test "$3" != ""; then
  dec=$3
fi

if test "$4" != ""; then
  mom=$4
fi

if test "$5" != ""; then
  res=$5
fi

if test "$6" != ""; then
  save=$6
fi

# Test which kind of simulation should be start
if test "$dec" == "DPM"; then
  sig=0
else
  dec=$scripts"/"$dec
fi

# If task arrays are used
if test "$SLURM_ARRAY_TASK_ID" == ""; then
	outprefix=$prefix
else
	outprefix=$prefix"_"$SLURM_ARRAY_TASK_ID
fi


# Run simulation

root -l -q -b -w $scripts"/"prod_sim.C\(\"$outprefix\",$nEvts,\"$dec\",$mom,\"$res\"\) &> $outprefix"_sim.log"
root -l -b -q -w $scripts"/"prod_dig.C\(\"$outprefix\"\) &> $outprefix"_dig.log"
root -l -b -q -w $scripts"/"prod_ideal_rec.C\(\"$outprefix\"\) &> $outprefix"_rec.log"
root -l -b -q -w $scripts"/"prod_pid.C\(\"$outprefix\"\) &> $outprefix"_pid.log"
root -l -b -q -w $scripts"/"prod_ana_AntiXi1690.C\(\"$outprefix\",$nEvts,$mom\) &> $outprefix"_ana.log"

NUMEV=`grep 'Generated Events' $outprefix"_sim.log"`
echo $NUMEV >> $outprefix"_pid.log"
   
cp  $outprefix"_par.root" $_target
cp  $outprefix"_pid.log" $_target
cp  $outprefix"_rec.root" $_target
cp  $outprefix"_pid.root" $_target
cp  $outprefix"_sim.log" $_target
cp  $outprefix"_output_ana.root" $_target
cp  $outprefix"_ana.log" $_target


# Standard and error output in different files
#SBATCH -o %j_%N.out.log
#SBATCH -e %j_%N.err.log

# tidy up
rm  $outprefix"_par.root"
rm  $outprefix"_sim.log"
rm  $outprefix"_pid.log"
rm  $outprefix"_sim.root"
rm  $outprefix"_pid.root"
rm  $outprefix"_dig.root"
rm  $outprefix"_dig.log"
rm  $outprefix"_ana.root"
rm  $outprefix"_output_ana.root"
rm  $outprefix"_rec.root"
rm  $outprefix"_rec.log"
rm  $outprefix"_ana.log"


# Execute application code
hostname; uptime; sleep 30; uname -a
