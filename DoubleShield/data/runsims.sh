#!/bin/bash
 
# Example job submission script

# This is a comment.
# Lines beginning with the # symbol are comments and are not interpreted by 
# the Job Scheduler.

# Lines beginning with #SBATCH are special commands to configure the job.
		
### Job Configuration Starts Here #############################################

# Export all current environment variables to the job (Don't change this)
#SBATCH --get-user-env 
			
# The default is one task per node
#SBATCH --ntasks=1 --cpus-per-task 24
#SBATCH --nodes=1
#SBATCH --partition=nodes
#SBATCH --hint=multithread

#request 10 minutes of runtime - the job will be killed if it exceeds this
##SBATCH --time=10:00
# Change email@example.com to your real email address
##SBATCH --mail-user=CJ.Barton@coyotes.usd.edu
##SBATCH --mail-type=ALL



### Commands to run your program start here ####################################
source /home/xu_lab/shared/CJMJ/Geant4-10.7/bin/geant4.sh

pwd

/home/xu_lab/shared/CJMJ/WLGD/warwick-legend -m ../DoubleShieldStudy.mac -o DoubleShield$counter-
#$WORKINGDIR/CoBaLEP $WORKINGDIR/mac/weighted.mac $(openssl rand 3 | od -DAn)
