# RNA Folding Scoring Project

## Project Description
This project focuses on the development and evaluation of an objective function for the RNA folding problem. Given a ribonucleotide chain, the goal is to determine the native fold by estimating Gibbs free energy using interatomic distance distributions from experimentally determined 3D structures.

The project consists of two main Python scripts:
1. **Training Script (`training.py`)**: Computes interatomic distances from PDB files and generates scoring values based on residue interactions.
2. **Visualization Script (`visualization.py`)**: Plots the scoring profiles as a function of interatomic distances and visualizes the scoring results.

## Data and Outputs
- **Input Data**: The scripts process RNA PDB files stored in a specified folder.
- **Outputs**:
  - **Plots** are saved in `results/Plot/`.
  - **Scoring values** are stored in `results/Output/`.

## Code Execution
The scripts are executed from the terminal as follows:

### Training and Scoring
```sh
python src/training.py <path_to_pdb_folder>
```
- Computes distance distributions.
- Generates scoring files.

### Visualization
```sh
python src/visualization.py <path_to_pdb_folder>
```
- Plots the scoring profiles.

## Functionality Overview
- **Training (`src/training.py`)**:
  - Computes interatomic distances using only C3’ atoms.
  - Considers 10 distance distributions for 10 base pair combinations (AA, AU, AC, AG, UU, UC, UG, CC, CG, GG).
  - Includes only intrachain distances and residues separated by at least 3 positions.
  - Calculates observed and reference frequencies to generate scoring values.

- **Visualization (`src/visualization.py`)**:
  - Reads the computed scoring values.
  - Generates plots showing the scoring profiles for different residue interactions.

## Folder Structure
```
structural_bioinformatics/
│── pdb/               # Directory containing input PDB files
│── results/
│   ├── Plot/          # Contains generated plots
│   ├── Output/        # Contains computed scoring values
│── src/
│   ├── training.py    # Training and scoring script
│   ├── visualization.py # Visualization script
│── README.md          # Project documentation
```

## Dependencies
Ensure that Python is installed along with the following libraries:
- `numpy`
- `matplotlib`
- `pandas`

Install dependencies using:
```sh
pip install numpy matplotlib pandas
```

## Usage Notes
- The project requires a dataset of RNA PDB files.
- The output scoring values and plots can be used to analyze RNA folding energies.
- This implementation follows good coding practices, including modularization and clear documentation.

## Contact
For any inquiries, please contact: milica.urosevic307@gmail.com

