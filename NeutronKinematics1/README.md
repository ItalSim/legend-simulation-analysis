# Neutron Kinematics Study (CJ's version, July-August 2023)

## Code Overview

This study represents a multi-staged analysis campaign to understand and characterize the muon-induced neutrons in a single reentrance tube LEGEND-1000-like geometry, including a solid neutron moderating shield. In particular, this study focuses on how the neutrons interact along the boundaries of the shield and inside of the shield. More information can be found in the included presentation and writeup.

The `code/` directory contains:

 - elossinshieldplots.C , used primarily to generate plots relating to neutron behavior at the shield border. Also outputs text information used for statistical analysis.
 - steplevelplots.C , used for gathering information on a (G4)step by step basis.
 - pathplots.C , used for plots relating to neutron mean free and total path lengths in materials of interest.
 - scatteringangleplots.C , used for plots relating to neutron scattering angle in elastic scatters.
 - transverseplots.C , used to check information about neutron penetration into the shield.

The `data/` directory contains:

 - The batch submission scripts, `jobsubmit.C` and `runsims.sh`
 - `WLGDRunAction.cc` and `WLGDSteppingAction.cc` files, for the modified data collection and output
 - The GIT hash of the stable release of `WLGD/ItalSim` modified for this study

The `results/` directory contains:

 - A PDF of a presentation of some of the results
 - A PDF of a writeup containing all of the results
 - A PDF of only the images/plots from the writeup

## Dependencies

 - ROOT v.6 or later
 - The job submission scripts rely on slurm, v.22
 - Geant4 10.7 with multithreading enabled, in order to compile the ItalSim version of WLGD
 - A C++ compiler with C++14 functionality
 - A bash environment (not sure what version)


## Instructions for execution

 - Fork and compile the relevant commit of WLGD/ItalSim
 - Swap the `WLGDRunAction.cc` and `WLGDSteppingAction.cc` files in this study with those in the `src/` directory of WLGD/ItalSim
 - Run WLGD/ItalSim, using the provided macro file and, if possible, the provided slurm job submission scripts
 - Use any of the analysis scripts in `code/`

Note on using analysis scripts in `code/`: all of the plotting commands are commented out at the bottom of each script. Uncomment whichever plots you'd like to make in the script, then execute it.