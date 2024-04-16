# Optical simulations analysis

## jobs

This is the directory which contains information necessary for distributed job submission to the MPIK lfs2 computing cluster. Jobs should always be submitted from this directory. All of the log files from job submission will write to this directory as well.

execute submith.sh (e.g. source submit.sh) to submit a job to the slurm queue system. submith.sh then calls run_py.sh to execute a python script. In this case pe_calulator.py, but it can be used to execute any python script.

## scripts

Available analysis scripts:

pe_calculator.py

pe_calculator.py can be run interactively from the shell for small files, otherwise submit jobs. 
It computes trigger efficiency (assuming Ar41 captures only, no background) for different photoelectron thresholds on a single light guide and different majority conditions (majority defined as number of light guide above PE_threhsold)

User-defined parameters:
- pe_calculator.py
  - the user must change the input directory to point to the output of "applyopticalmap.C" script
  - PRELIMINARY: this script was used during Vancouver general meeting in 2023 to generate the tagging efficiency vs PE threshold (at different multiplicities) plot. For more refined analysis, change this script accordingly.

## output

The previous sections have already given a good idea of the output structure. The directories are as follows:

- output/

Ensure that all scripts point to the appropriate subdirectory in THIS directory (whichever directory the README file is in).