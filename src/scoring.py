import os
import argparse
import numpy as np
import math

# Class for extracting C3' atoms from PDB files
class FindC3:
    def __init__(self, pdb_file_path):
        self.pdb_file_path = pdb_file_path
        self.c3_atoms = []

    def process(self):
        """Extracts only C3' atoms from a PDB file."""
        with open(self.pdb_file_path, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith('ATOM') and "C3'" in line:
                    parts = line.split()
                    self.c3_atoms.append(parts[:11])  # Extract relevant columns
        return self.c3_atoms

# Function to calculate distance
def compute_distance(coord1, coord2):
    coord1 = np.array(list(map(float, coord1)))
    coord2 = np.array(list(map(float, coord2)))
    return np.linalg.norm(coord1 - coord2)

# Function to compute Gibbs free energy
def compute_gibbs_free_energy(c3_list, scoring_data):
    """Calculates Gibbs free energy based on scoring data."""
    gibbs_free_energy = 0.0
    for i, atom1 in enumerate(c3_list):
        for j, atom2 in enumerate(c3_list[i+1:], start=i+1):
            try:
                res1 = int(''.join(filter(str.isdigit, atom1[5])))  # Extract only numeric part
                res2 = int(''.join(filter(str.isdigit, atom2[5])))
                if res1 + 4 <= res2:
                    base_pair = ''.join(sorted((atom1[4], atom2[4])))
                    distance = compute_distance(atom1[8:11], atom2[8:11])
                    if 0.5 <= distance <= 19.5 and base_pair in scoring_data:
                        gibbs_free_energy += linear_interpolation(base_pair, distance, scoring_data)
            except ValueError:
                continue  # Skip non-numeric residue indices
    return gibbs_free_energy

# Function to perform linear interpolation
def linear_interpolation(base_pair, distance, scoring_data):
    """Performs linear interpolation to calculate energy values."""
    energy_values = scoring_data[base_pair]
    lower_idx = math.floor(distance) - 1
    upper_idx = math.floor(distance)
    energy_before = float(energy_values[lower_idx])
    energy_after = float(energy_values[upper_idx])
    return energy_before + (distance - (lower_idx + 0.5)) * (energy_after - energy_before)

# Function to load scoring data
def load_scoring_data(scoring_dir):
    """Loads precomputed scoring values from files."""
    scoring_data = {}
    for bp in ['AA', 'AU', 'AC', 'AG', 'UU', 'CU', 'GU', 'CC', 'CG', 'GG']:
        with open(f"{scoring_dir}/{bp}", 'r') as f:
            scoring_data[bp] = [float(line.strip()) for line in f.readlines()]
    return scoring_data

# Main function
def main(pdb_file, scoring_dir):
    """Computes Gibbs free energy for a given PDB file."""
    c3_atoms = FindC3(pdb_file).process()
    scoring_data = load_scoring_data(scoring_dir)
    gibbs_energy = compute_gibbs_free_energy(c3_atoms, scoring_data)
    print(f"Gibbs Free Energy: {gibbs_energy:.2f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_file", help="Path to the PDB file.")
    parser.add_argument("scoring_dir", help="Path to the directory containing scoring files.")
    args = parser.parse_args()
    main(args.pdb_file, args.scoring_dir)
