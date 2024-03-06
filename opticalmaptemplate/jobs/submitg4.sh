#!/bin/bash

### Job configuration starts here #############################################

# Export all current environment variables to the job
#SBATCH --get-user-env

# One task per node (single-threading)
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
##SBATCH --partition=nodes

# Request x minutes of runtime - the job will be killed if it exceeds this
##SBATCH --time=x:00

# Redirect output (optional)
##SBATCH --output=

### Commands to run the program start here ####################################

pwd
#Optical map simulations
source /home/cbarton/.bashrc

singularity exec /lfs/l1/legend/users/morella/legend1k_simulation/remage-base.sif /bin/bash $SDC/opticalmap/4-innermap-Dec2023/jobs/rg4.sh
