# Double Shield simulation study

This study was performed in order to quickly gauge the effects of splitting a 10 cm thick solid neutron shield into two 5 cm thick shields with some liquid argon in between. The geometry chosen is a current (as of May 2023) version of the LEGEND-1000 "single reentrance tube" geometry, with some modifications. The outer diameters and heights of the shields are 4 meters and 3 meters, respectively. This means there is a 50 cm gap between the outer edge of each shield on all sides. Since the shields are 5 cm thick, this means there is a 45 cm gap between the inner edge of the big shield and the outer edge of the small shield. This space, as well as the usual volumes, are filled with liquid argon.

This idea comes from Bernhard and was shared to the ItalSim Slack channel on the 16th of May by Natalia Di Marco.

Potential variables: material, height/radius of larger shield, and size difference between larger shield
and smaller shield (which determines the height and radius of the smaller shield).

## Code Overview

There are no analysis scripts for this simple study. Information was extracted from the output files using the ROOT CLI and exported to an external plotting program which cannot be included here. More information is available near the bottom of this readme.

The `data/` directory includes:
 - The macro file used in the Geant4 input, `DoubleShield.mac`
 - The batch submission scripts, `jobsubmit.C` and `runsims.sh`
 - The modified `WLGDDetectorConstruction.cc` file (has a double shield)
 - The GIT hash of the version of WLGD/ItalSim used in this study

The results/ directory contains only the PDF of the report which was created for the conclusion of this study.

Additional info:

The output file was also modified to only include neutron capture information, to save on data storage space (this may not be necessary on a sufficiently large cluster).
10 sets of 10 million muons were simulated. Each set was comprised of 10 batches of 1 million muons per job, for a total of 100 million muons.
Each job was executed with multithreading enabled, and thread count set to 24, the maximum available on the standard processing nodes used for simulation.


## Dependencies

 - ROOT v.6 or later
 - The job submission scripts rely on slurm, v.22
 - Geant 4 10.7 with multithreading enabled, in order to compile the ItalSim version of WLGD
 - A C++ compiler with C++14 functionality
 - A bash environment (not sure what version)


## Instructions for execution

Ensure all dependencies are available
Execute the recommended commit of WLGD, with the `DoubleShield.mac` macro in the data/ directory (make sure to replace the `src/WLGDDetectorConstruction.cc` in your fork of WLGD with the one included in data/)
If using slurm, you may use `jobsubmit.C`, which calls `runsims.sh`, to automatically submit batch jobs
To extract the number of Ge77 produced via neutron capture in the detectors, combine (if necessary) and load the ROOT data file(s) into the ROOT CLI and use these commands:

> Score->Draw("nCapture_EventID","nCapture_A==77&&nCapture_Z==32");
> cout << htemp->Integral() << endl;



This study was meant to be brief and simple, so does not have full functionality. Any user who is interested can feel free to expand on and modify the included functionality to suit their needs.
