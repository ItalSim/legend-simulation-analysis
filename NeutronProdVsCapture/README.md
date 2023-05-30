# Neutron Capture Analysis

This code performs an analysis of neutron capture events in a simulation dataset. It reads ROOT files containing information about neutrons and captures on different volumes. The main purpose is to compare the neutron capture events with and without moderators and visualize the production volumes of the captured neutrons.

## Code Overview

The code is written in Python 3 and consists of the following files:

- `code/neutronProdvsCapture.py`: The main script that orchestrates the execution of the code.
- `data/`: information on how to reproduce the simulated data.
- `results/`: example of plot made using this script
- `README.md`: This file, providing documentation and instructions.

## Dependencies

The code relies on the following Python 3 packages:

- `uproot`: For reading ROOT files.
- `numpy`: For numerical operations.
- `matplotlib`: For data visualization.

Please make sure to install these dependencies before running the code.

## Usage

To use the code, follow these steps:

1. Clone the repository and navigate to the code directory.
2. Set up the required dependencies (see "Dependencies" section).
3. Modify the code as needed for your specific use case.
4. Execute the `neutron_capture_analysis.py` script using the command `python3 neutron_capture_analysis.py [-f <flag_plot>] [-p <path1>] [-d <path2>] [-n <number_files>] [-v]`.

The available options for the `main.py` script are as follows:

- `-f`, `--flag-plot`: Choose between "noMod", "Mod", or "both" to specify the type of plot to generate.
- `-p`, `--path1`: Path to the data for the specified flag-plot option.
- `-d`, `--path2`: Path to additional data (required only for the "both" flag-plot option).
- `-n`, `--number`: Number of files to analyze.
- `-v`, `--verbose`: Enable verbose output.

Ensure that you provide the correct paths to the data files and specify the desired flag-plot option to generate the appropriate plots.

## Functionality

The code provides the following functions:

- `prodPosVsCap`: Finds neutron IDs that have been captured on Ge77 detectors and extracts their positional information.
- `countVolumes`: Counts the occurrences of different volumes in the neutron capture data.
- `loadAndCalcProd`: Loads data from ROOT files, performs calculations, and returns the histogram of captured neutrons.
- `plotProdvsCap`: Plots the production volume of neutrons captured on Ge detectors.

## Examples

Here are a few example commands to run the code:

- To generate a plot without moderators:
`python main.py -f noMod -p /path/to/data -n 200`

- To generate both plots (with and without moderators):
`python main.py -f both -p /path/to/noMod/data -d /path/to/Mod/data -n 200`

Note: Make sure to replace /path/to/data with the actual paths to your data files.

Feel free to explore and adapt the code to suit your specific needs.
