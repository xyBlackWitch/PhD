#! /bin/bash

# start pandaroot software
. "/lustre/nyx/panda/jpuetz/pandaroot/trunk/build/config.sh"

# Task name
#SBATCH -J Analysis_XiMinus1820

# Run time limit
#SBATCH --time=0:20:00

#Working directory on shared storage
#SBATCH -D /lustre/nyx/panda/jpuetz/data/XiMinus1820_Lambda0_K/spin_3half/run1

# Script directory (If different from working directory!)
scripts="/lustre/nyx/panda/jpuetz/myscripts"

# Default parameters
prefix=9999
nEvts=1000
mom=12.3485


# Parameters set by user
if test "$1" != ""; then
  prefix=$1
fi

if test "$2" != ""; then
  nEvts=$2
fi

if test "$3" != ""; then
  mom=$3
fi


# If task arrays are used
if test "$SLURM_ARRAY_TASK_ID" == ""; then
	outprefix=$prefix
else
	outprefix=$prefix"_"$SLURM_ARRAY_TASK_ID
fi


# Run simulation

root -l -b -q -w $scripts"/"prod_ana_Xi1820.C\(\"$outprefix\",$nEvts,$mom\) &> $outprefix"_ana.log"


# Standard and error output in different files
#SBATCH -o %j_%N.out.log
#SBATCH -e %j_%N.err.log


# Execute application code
hostname; uptime; sleep 30; uname -a
