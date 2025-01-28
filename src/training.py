import os
import argparse
import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import math

# Constants
distance_bins = np.arange(21)  # 0 to 20 Å
base_pairs = ["AA", "AU", "AC", "AG", "UU", "UC", "UG", "CC", "CG", "GG", "XX"]

# Classes
class AtomProcessor:
    """
    Processes atoms from a PDB file and extracts C3' atoms.
    """
    def __init__(self, pdb_file_path):
        self.pdb_file_path = pdb_file_path
        self.atoms = []
        self.c3_atoms = []

    def extract_atoms(self):
        """Extracts all atoms from the PDB file."""
        with open(self.pdb_file_path, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith('ATOM'):
                    cols = line.strip()
                    atom = [cols[:6], cols[6:11], cols[12:16], cols[16:17], cols[17:20], cols[21:22], cols[22:26],
                            cols[26:27], cols[30:38], cols[38:46], cols[46:54], cols[54:60], cols[60:66], cols[76:78],
                            cols[78:80]]
                    self.atoms.append([i.strip() for i in atom])

    def extract_c3_atoms(self):
        """Filters and retains only C3' atoms."""
        self.c3_atoms = [atom for atom in self.atoms if atom[2] == "C3'"]

    def process(self):
        """Runs the atom extraction and filtering process."""
        self.extract_atoms()
        self.extract_c3_atoms()
        return self.c3_atoms

# Functions
def create_frequency_matrix():
    """Creates a dictionary to store frequencies for each base pair."""
    return {bp: np.zeros(20) for bp in base_pairs}

def calculate_distance(coord1, coord2):
    """Calculates the Euclidean distance between two points."""
    return np.linalg.norm(np.array(coord1) - np.array(coord2))

def read_structure(filename):
    """Reads a PDB file and returns a BioPython structure object."""
    parser = PDBParser(QUIET=True)
    return parser.get_structure("RNA", filename)

def process_c3_atoms_from_file(pdb_file_path):
    """Uses AtomProcessor to extract C3' atoms and their coordinates."""
    processor = AtomProcessor(pdb_file_path)
    c3_list = processor.process()
    c3_coords = np.array([list(map(float, atom[8:11])) for atom in c3_list])
    residues = c3_list
    return c3_coords, residues

def compute_distance_frequencies(c3_coords, residues, frequency_matrix):
    """Computes distance-based frequencies for each base pair."""
    for i, res1 in enumerate(residues):
        for j, res2 in enumerate(residues):
            if i < j and int(res2[6]) >= int(res1[6]) + 4:
                pair = ''.join(sorted((res1[4], res2[4])))
                dist = calculate_distance(c3_coords[i], c3_coords[j])
                if 0 <= dist < 20:
                    bin_index = int(dist)
                    if pair in frequency_matrix:
                        frequency_matrix[pair][bin_index] += 1
                        frequency_matrix["XX"][bin_index] += 1
    return frequency_matrix

def normalize_frequencies(frequency_matrix):
    """Normalizes frequency values to probabilities."""
    for pair, freq_array in frequency_matrix.items():
        total = np.sum(freq_array)
        if total > 0:
            frequency_matrix[pair] = freq_array / total
    return frequency_matrix

def compute_pseudo_energies(frequency_matrix):
    """Calculates pseudo-energy scores from normalized frequencies."""
    energies = {}
    reference_freq = frequency_matrix["XX"]
    for pair, freq_array in frequency_matrix.items():
        if pair == "XX":
            continue
        scores = []
        for obs_freq, ref_freq in zip(freq_array, reference_freq):
            if obs_freq > 0 and ref_freq > 0:
                score = -np.log(obs_freq / ref_freq)
                scores.append(min(score, 10))
            else:
                scores.append(10)
        energies[pair] = scores
    return energies

def save_energy_scores(energies, output_folder):
    """Saves pseudo-energy scores to text files."""
    os.makedirs(output_folder, exist_ok=True)
    for pair, scores in energies.items():
        output_path = os.path.join(output_folder, f"{pair}_scores.txt")
        with open(output_path, "w") as file:
            file.write("\n".join(map(str, scores)))

def plot_energy_profiles(energies, output_folder):
    """Generates and saves plots of energy profiles."""
    os.makedirs(output_folder, exist_ok=True)
    for pair, scores in energies.items():
        plt.figure()
        plt.plot(range(len(scores)), scores, label=pair)
        plt.xlabel("Distance bin (Å)")
        plt.ylabel("Pseudo-energy")
        plt.title(f"Energy Profile: {pair}")
        plt.legend()
        plt.savefig(os.path.join(output_folder, f"{pair}_profile.png"))
        plt.close()

# Main Execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_folder", help="Path to the folder containing PDB files.", required=True)
    parser.add_argument("--output_folder", help="Path to the folder for saving output scores.", required=True)
    args = parser.parse_args()

    # Initialize frequency matrix
    frequency_matrix = create_frequency_matrix()

    # Process each PDB file in the folder
    for pdb_file in os.listdir(args.pdb_folder):
        if pdb_file.endswith(".pdb"):
            pdb_path = os.path.join(args.pdb_folder, pdb_file)
            c3_coords, residues = process_c3_atoms_from_file(pdb_path)
            frequency_matrix = compute_distance_frequencies(c3_coords, residues, frequency_matrix)

    # Normalize frequencies and compute energies
    frequency_matrix = normalize_frequencies(frequency_matrix)
    pseudo_energies = compute_pseudo_energies(frequency_matrix)

    # Save results and generate plots
    save_energy_scores(pseudo_energies, args.output_folder)
    plot_energy_profiles(pseudo_energies, args.output_folder)

    print(f"Results saved in {args.output_folder}")
