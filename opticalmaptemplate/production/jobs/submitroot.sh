#!/bin/bash

### Job configuration starts here #############################################

# Export all current environment variables to the job
#SBATCH --get-user-env

##Not sure the single-threading is optimal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
##SBATCH --partition=nodes

##Sometimes certain nodes perform worse - may be necessary to exclude them
##SBATCH --exclude=lfl10
##SBATCH --mem=30G

# Request x minutes of runtime - the job will be killed if it exceeds this
##SBATCH --time=x:00

# Redirect output (optional)
##SBATCH --output=$SDC/-/jobs/slurm%A.out

### Commands to run the program start here ####################################

pwd
source /home/cbarton/.bashrc

root -b -x -q "$SDC/opticalmap/4-innermap-Dec2023/scripts/buildopticalmap.C($counter)"
#root -b -x -q "$SDC/opticalmap/4-innermap-Dec2023/scripts/mergeopticalmap.C($counter)"
