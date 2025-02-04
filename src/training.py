import os
import argparse
import numpy as np
import math
from Bio.PDB import PDBParser

# Class for extracting C3' atoms from PDB files
class PDBCleaner:
    def __init__(self, pdb_file_path):
        self.pdb_file_path = pdb_file_path
        self.c3_atoms = []

    def process(self):
        """Extracts only C3' atoms from a PDB file."""
        with open(self.pdb_file_path, 'r') as pdb_file:
            for line in pdb_file:
                parts = line.split()
                if line.startswith('ATOM') and "C3'" in line:
                    self.c3_atoms.append(parts[:11])  # Extract relevant columns
        return self.c3_atoms

# Function to calculate Euclidean distance
def compute_distance(coord1, coord2):
    """Computes Euclidean distance between two points in 3D space."""
    try:
        coord1 = np.array(list(map(float, coord1)))
        coord2 = np.array(list(map(float, coord2)))
        return np.linalg.norm(coord1 - coord2)
    except ValueError:
        return None

# Function to compute frequency distribution
def update_frequencies(c3_list, frequency_data):
    """Calculates distances between atoms and updates frequency distribution."""
    for i, atom1 in enumerate(c3_list):
        for j, atom2 in enumerate(c3_list[i+1:], start=i+1):
            try:
                if atom1[5].isdigit() and atom2[5].isdigit() and (int(atom1[5]) + 4) <= int(atom2[5]):
                    base_pair = ''.join(sorted((atom1[4], atom2[4])))
                    distance = compute_distance(atom1[8:11], atom2[8:11])
                    if distance is not None and 0 <= distance < 20:
                        bin_index = math.floor(distance)
                        if base_pair in frequency_data:
                            frequency_data[base_pair][bin_index] += 1
                            frequency_data['XX'][bin_index] += 1
            except IndexError:
                continue
    return frequency_data

# Function to normalize frequencies
def normalize_frequencies(frequency_data):
    """Converts raw frequency counts into probabilities."""
    for base_pair, freq_values in frequency_data.items():
        total_count = np.sum(freq_values)
        if total_count > 0:
            frequency_data[base_pair] = freq_values / total_count
    return frequency_data

# Function to compute pseudo-energy values
def calculate_pseudo_energies(frequency_data):
    """Derives pseudo-energy values from normalized frequencies."""
    energy_values = {}
    reference_distribution = frequency_data['XX']
    for base_pair, freq_values in frequency_data.items():
        if base_pair == 'XX':
            continue
        energy_scores = [-np.log(obs_freq / ref_freq) if obs_freq > 0 and ref_freq > 0 else 10 
                         for obs_freq, ref_freq in zip(freq_values, reference_distribution)]
        energy_values[base_pair] = [min(score, 10) for score in energy_scores]
    return energy_values

# Function to save results
def save_energy_results(energy_values):
    """Saves pseudo-energy values into separate text files."""
    os.makedirs("Energy", exist_ok=True)
    for base_pair, scores in energy_values.items():
        with open(f"Energy/{base_pair}", "w") as file:
            file.write("\n".join(map(str, scores)))

# Main function
def main(pdb_directory):
    """Processes all PDB files in a directory and computes pseudo-energies."""
    frequency_matrix = {bp: np.zeros(20) for bp in ['AA', 'AU', 'AC', 'AG', 'UU', 'CU', 'GU', 'CC', 'CG', 'GG', 'XX']}
    for pdb_filename in os.listdir(pdb_directory):
        if pdb_filename.endswith(".pdb"):
            pdb_path = os.path.join(pdb_directory, pdb_filename)
            c3_atoms = PDBCleaner(pdb_path).process()
            frequency_matrix = update_frequencies(c3_atoms, frequency_matrix)
    
    frequency_matrix = normalize_frequencies(frequency_matrix)
    pseudo_energies = calculate_pseudo_energies(frequency_matrix)
    save_energy_results(pseudo_energies)
    print("Results have been saved in the 'Energy' directory.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_directory", help="Path to the directory containing PDB files.")
    args = parser.parse_args()
    main(args.pdb_directory)
