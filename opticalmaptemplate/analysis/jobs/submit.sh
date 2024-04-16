#!/bin/bash

### Job configuration starts here #############################################

# Export all current environment variables to the job
#SBATCH --get-user-env

# One task per node (single-threading)
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=20G

# Request x minutes of runtime - the job will be killed if it exceeds this
##SBATCH --time=x:00

### Commands to run the program start here ####################################

singularity exec /lfs/l1/legend/software/singularity/legendexp_legend-software_latest.sif /bin/bash run_py.sh
