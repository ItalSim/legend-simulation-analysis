# Optical map production

This is a subroutine of the ItalSim branch of the warwick-legend simulation module based on Geant4. Given an ordered set of coordinate pairs (X1,Y1,Z1) and (X2,Y2,Z2), construct an improved version of the LEGEND-200 concept of an "optical map".

An optical map is a tool used to simplify simulations analysis. Optical simulations in Geant4 are costly both in terms of memory and in terms of processing resources required. For a fixed geometry, it is possible to populate the entire relevant volume with photons one time in a dedicated campaign and reference the results, instead of simulating optics each time they are necessary for a new simulation campaign. These results are referenced by use of a 3D map of the volume, broken into small cubes called 'voxels'. Each voxel contains coordinates for the piece of volume which it encloses, and the average detection probability of a photon which originates in that piece. A good map will have sufficiently small voxels, with sufficiently high statistics in the probability determination.

This optical map, however, takes the process further. Instead of saving the voxel's location and integral probability, it saves the voxel's location, integral probability, and the X2,Y2,Z2 values of EVERY photon in that voxel which hits the target location. This way, the exact position of every hit is saved, and can be sampled randomly in future applications.

Below is a description of this optical map production template. It is assumed that the user has an appropriate version of WLGD/ItalSim installed.


## jobs

This is the directory which contains information necessary for distributed job submission to the MPIK lfs2 computing cluster. The macro used to generate the raw files for map generation should also be stored here, for documentation purposes. Jobs should always be submitted from this directory. All of the log files from job submission will write to this directory as well.

For Geant4 jobs (i.e. producing the raw files), 

jobsubmitg4.C calls
submitg4.sh which calls
rg4.sh

User-defined parameters:
- jobsubmitg4
  - int i: informs the script how many jobs should be submitted (default 100)
- submitg4
  - source command: should point to the user's .bashrc file
  - singularity exec: only the third argument should be changed, to point to rg4.sh in this directory
- rg4.sh
  - First source command: should point to the user's .bashrc file
  - Second source command: should point to the user's ItalSim directory, specifically, to the env.sh file inside
  - Third command: first part points to the WLGD executable, second part points to a macro for optical map generation (one version included in the 'jobs' directory for reference), third part should point to this directory's 'raw' output subdirectory (/output/raw).


For ROOT jobs (i.e. producing the maps),

jobsubmitroot.C calls
submitroot.sh

User-defined parameters:
- jobsubmitroot.C
  - int i: informs the script how many jobs should be submitted (default 100)
- submitroot
  - source command: should point to the user's .bashrc file
  - root command: should point to this directory's 'scripts' subdirectory (/scripts)
  - The user should decide whether it is appropriate to run the builder script or the merger script (typically in that order)


## scripts

The scripts should be applied in this order:

buildopticalmap.C
mergeopticalmap.C
opticalmappostproc.C

The first two are applied from the 'jobs' directory (batch executed), but the third is typically run interactively from the shell.

User-defined parameters:
- buildopticalmap.C
  - The user must change the input/output directories in the script to point to the 'raw' and the 'built' directories in /output
  - The user is responsible for defining the voxel size and the limits in the x, y, and z directions to enclose the sampling volume
- mergeopticalmap.C
  - The user must change the input/output directories in the script to point to the 'built' and the 'merge' directories in /output
  - The user is responsible for defining the voxel size and the limits in the x, y, and z directions to enclose the sampling volume (again)
  - The user must ensure the volume of the REAL sampling volume (not the enclosing prism, see below) is accurately calculated in the script by adjusting the values where necessary.
  - The user must change the numbertoprocess variable according to the instructions in the script
- opticalmappostproc.C
  - The user must change the input directory to point to the 'merged' directory, and the output should point wherever the user would like the map to be stored
  - The user should add the number of voxels in the x, y, and z direction (labeled xbins,ybins,zbins in the script)
  
Special note about the sampling volume limits: the optical map is always a rectangular prism. The sampling volume shape from the optical simulations must be fully enclosed by this prism. Examples for the sampling limits of the 'inner' optical map and 'outer' optical map along the sides of the shield have been included, but these may change as the LEGEND design changes.

Also, again, sorry about how complicated the merge script is.

## output

The previous sections have already given a good idea of the output structure. The directories are as follows:

- output/
  - raw
  - built
  - merged

Ensure that all scripts point to the appropriate subdirectory in THIS directory (whichever directory the README file is in).

# 
A final note: this optical map script has been developed for the ItalSim version of warwick-legend, and it is intended to work with output from this module. The parameters may be adjusted for processing of other ROOT-based coordinate data, but this does not come with a guarantee of functionality.