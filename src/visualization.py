import os
import sys
from Bio.PDB import PDBParser, Structure
from Bio.PDB.Superimposer import Superimposer
from scipy.spatial.transform import Rotation
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt

def calculate_distance(coord1, coord2):
    """Calculate the Euclidean distance between two points."""
    return np.linalg.norm(np.array(coord1) - np.array(coord2))

def read_structure(filename: str) -> Structure:
    """Reads a PDB structure file and parses it into a Structure object."""
    parser = PDBParser(QUIET=True)
    return parser.get_structure("RNA", filename)

def extract_c3_atoms(structure: Structure):
    """Extracts coordinates of C3' atoms from a given structure."""
    return np.array([
        atom.get_coord() for atom in structure.get_atoms() if atom.get_name() == "C3'"
    ])

def compute_distance_distributions(c3_coords):
    """Compute distance distributions for given C3' atom coordinates."""
    distributions = {pair: np.zeros(20) for pair in ["AA", "AU", "AC", "AG", "UU", "UC", "UG", "CC", "CG", "GG"]}
    total_counts = {pair: 0 for pair in distributions}

    for i in range(len(c3_coords)):
        for j in range(i + 4, len(c3_coords)):  # Minimum sequence separation of 4
            distance = calculate_distance(c3_coords[i], c3_coords[j])
            if distance <= 20:
                bin_index = int(distance)  # Each bin has a width of 1 Å
                pair = "XX"  # Placeholder for scoring logic
                distributions[pair][bin_index] += 1
                total_counts[pair] += 1

    return distributions, total_counts

def compute_scores(distributions, total_counts):
    """Compute log-ratio scores based on observed and reference frequencies."""
    scores = {}
    for pair, distribution in distributions.items():
        pair_scores = []
        for count in distribution:
            observed_frequency = count / total_counts[pair] if total_counts[pair] > 0 else 0
            reference_frequency = count / total_counts["XX"] if total_counts["XX"] > 0 else 0
            if observed_frequency > 0 and reference_frequency > 0:
                score = -np.log(observed_frequency / reference_frequency)
                score = min(score, 10)  # Cap the score at 10
            else:
                score = 10
            pair_scores.append(score)
        scores[pair] = pair_scores
    return scores

def plot_interaction_profiles(scores, output_folder):
    """Plot interaction profiles for each pair of bases."""
    os.makedirs(output_folder, exist_ok=True)
    for pair, pair_scores in scores.items():
        plt.figure(figsize=(8, 6))
        plt.plot(range(len(pair_scores)), pair_scores, marker="o", linestyle="-", label=f"Pair: {pair}")
        plt.title(f"Interaction Profile for {pair}")
        plt.xlabel("Distance (Å)")
        plt.ylabel("Score")
        plt.ylim(0, 10)
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(output_folder, f"{pair}_profile.png"))
        plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RNA Folding Objective Function Training Script")
    parser.add_argument("--pdb_folder", type=str, required=True, help="Path to the folder containing PDB files")
    parser.add_argument("--output_folder", type=str, required=True, help="Path to the folder to save output scores and plots")
    args = parser.parse_args()

    pdb_folder = args.pdb_folder
    output_folder = args.output_folder
    os.makedirs(output_folder, exist_ok=True)

    all_distributions = {}
    all_total_counts = {}

    for pdb_file in os.listdir(pdb_folder):
        if pdb_file.endswith(".pdb"):
            pdb_path = os.path.join(pdb_folder, pdb_file)
            structure = read_structure(pdb_path)
            c3_coords = extract_c3_atoms(structure)
            distributions, total_counts = compute_distance_distributions(c3_coords)

            for pair in distributions:
                if pair not in all_distributions:
                    all_distributions[pair] = np.zeros(20)
                    all_total_counts[pair] = 0
                all_distributions[pair] += distributions[pair]
                all_total_counts[pair] += total_counts[pair]

    scores = compute_scores(all_distributions, all_total_counts)

    for pair, pair_scores in scores.items():
        output_file = os.path.join(output_folder, f"{pair}_scores.txt")
        with open(output_file, "w") as f:
            for score in pair_scores:
                f.write(f"{score}\n")

    plot_interaction_profiles(scores, output_folder)

    print(f"Results and plots have been saved in the folder: {output_folder}")
