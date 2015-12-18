#!/bin/bash
#
# Author: bjpjones@fnal.gov from echurch@fnal.gov from dbox@fnal.gov
#
# A script to run the optical library building job
#
#
# To run this job:
#
# jobsub -X509_USER_PROXY /scratch/[user]/grid/[user].uboone.proxy -N [NoOfJobs] -q OpticalLibraryBuild_Grid.sh `whoami` `pwd`
#
# You will get outputs in the area specified by the "outstage" variable 
# which is specified below.
#
# The form of the output is one file for each few voxels. These then need 
# stitching together, which is done after all jobs are done, with a
# dedicated stitching script.
#

umask 0002
verbose=T

# Copy arguments into meaningful names.
cluster=${CLUSTER}
process=${PROCESS}
user=$1
submitdir=$2


# Directory in which to put the output files.
outstage=/uboone/data/users/bjpjones/OpticalProduction


# Library building parameters

# In each grid job, do this many voxels:
NVoxelsPerJob=240

# In each voxel, run this many photons:
NPhotonsPerVoxel=30000

NTopVoxel=2250000

# This works out which voxels this job should focus on: 
FirstVoxel=`echo "($NVoxelsPerJob * $PROCESS ) % $NTopVoxel" | bc`
LastVoxel=`echo "(($NVoxelsPerJob * $PROCESS ) + $NVoxelsPerJob - 1 ) % $NTopVoxel" | bc`


# And then tell the user about it:
echo "This job will run from voxel $FirstVoxel to $LastVoxel, generating $NPhotonsPerVoxel in each"


# Most of what comes next is from echurch's grid submission script.
 

echo "Input arguments:"
echo "Cluster:    " $cluster
echo "Process:    " $process
echo "User:       " $user
echo "SubmitDir:  " $submitdir
echo " "

# Do not change this section.
ORIGDIR=`pwd`
TMP=`mktemp -d ${_CONDOR_SCRATCH_DIR:-/var/tmp}/working_dir.XXXXXXXXXX`
TMP=${TMP:-${_CONDOR_SCRATCH_DIR:-/var/tmp}/working_dir.$$}

{ [[ -n "$TMP" ]] && mkdir -p "$TMP"; } || \
  { echo "ERROR: unable to create temporary directory!" 1>&2; exit 1; }
trap "[[ -n \"$TMP\" ]] && { cd ; rm -rf \"$TMP\"; }" 0

echo "Temp directory : $TMP"
cd $TMP


cp -r $ORIGDIR/* .
# End of the section you should not change.


# Construct the run-time configuration file for this job;
# Make the random number seeds a function of the process number.
generatorSeed=$(( $process *23 + 31))
g4Seed=$(( $process *41 + 37))


# Copy fcl file and configure for this process 
cp /uboone/app/users/bjpjones/GridReadyNew/prodsingle_buildopticallibrary.fcl thisjob.fcl

echo "physics.producers.generator.FirstVoxel: $FirstVoxel" >> thisjob.fcl
echo "physics.producers.generator.LastVoxel: $LastVoxel" >> thisjob.fcl
echo "physics.producers.generator.N: $NPhotonsPerVoxel">> thisjob.fcl

# No need to set random seeds - generated from machine state
#echo "physics.producers.generator.RandomSeed: $generatorSeed">> thisjob.fcl
#echo "physics.producers.largeant.RandomSeed: $g4Seed">> thisjob.fcl


unset LD_LIBRARY_PATH

# Establish environment and run the job.
export GROUP=uboone
export EXPERIMENT=uboone
export EXTRA_PATH=lar
export HOME=/uboone/app/users/bjpjones/GridReadyNew

## This sets all the needed FW and SRT and LD_LIBRARY_PATH envt variables. 
## Then we cd back to our TMP area. EC, 23-Nov-2010.

cd $TMP

# Setup larsoft within the uboone frozen release
source /grid/fermiapp/lbne/lar/code/larsoft/releases/S2013.03.11/setup/setup_larsoft_fnal.sh -r S2013.03.11
srt_setup SRT_BASE_RELEASE=S2013.03.11





# Run the job
echo "Starting job"
lar -c thisjob.fcl -n $NVoxelsPerJob >& thisjob.log
echo "Job completed"


# Make sure the user's output staging area exists.
test -e $outstage || mkdir $outstage
if [ ! -d $outstage ];then
   echo "File exists but is not a directory."
   outstage=/grid/data/lbne/outstage/nobody
   echo "Changing outstage directory to: " $outstage 
   exit
fi

# Make a directory in the outstage area to hold all files from this job.
echo "Trying to make directory :  ${outstage}/${cluster}_${process} "
mkdir ${outstage}/${cluster}_${process}


# Copy all files from the working directory to the output staging area.
/grid/fermiapp/minos/scripts/cpn *.root ${outstage}/${cluster}_${process}
/grid/fermiapp/minos/scripts/cpn *.log ${outstage}/${cluster}_${process}
/grid/fermiapp/minos/scripts/cpn *.fcl ${outstage}/${cluster}_${process}
# Make sure EXPERIMENTana (under which name the jobs run on the grid)
# writes these as group rw, so you can rm 'em, etc.
#chmod -R g+rw $outstage

exit 0
