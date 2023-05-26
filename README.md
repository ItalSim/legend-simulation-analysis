# Legend Simulation Analysis

This repository contains analysis code and documentation for the legend simulation.

## Description

The legend simulation analysis is based on the Warwick-LEGEND fast simulation. This project aims to analyze simulation data and extract meaningful insights in a reproducible manner.

## Getting Started

To get started with the legend simulation analysis, you can follow two different paths:

### 1. Contributing to the Repository:

1. Fork this repository.

2. Clone the repository to your local machine.

3. Create a branch to develop your code.

4. Create an analysis folder with a meaningful name (avoid spaces, start each word with a capital letter, e.g., NeutronProductionVsCapture), following the directory structure:

#### Directory Structure
- `data`: Here you must add:
        1. the macro files used to create the simulation you are going to analyse;
        2. a .txt file with the 'git hash' of the Warwick Simulation that you ran the simulation.
- `code`: Contains the analysis code and scripts;
- `requeriments`: a list of all the packages (with version) needed to run your code;
- `results`: Stores the .png generated visualizations (500kb maximum for each file, as a way to keep a history).

5. Create a `README.md` file in your analysis folder with the following information:
- Name, date, and contact information.
- A brief explanation of what triggered this analysis and the expected outcomes.
- A link to any relevant presentation, if available.
- Instructions on how to run the code.

6. Include your analysis in the Git Wiki page list.

7. After testing the code, create a Pull Request to the main repository.

### 2. Reproduce the Analysis in this Repository:

1. Clone this repository.

2. Navigate to the folder containing the analysis you want to reproduce.

3. Use the `data` folder to find the relevant information needed to reproduce the data required for the analysis.

4. Install the dependencies specified in the `requirements.txt` file.

5. Run the code using the instructions provided in the `README.md` file.

## Contributing

Contributions are welcome! If you want to contribute to the legend simulation analysis, please follow the guidelines in the `CONTRIBUTING.md` file.
