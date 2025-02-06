import os
import argparse
import numpy as np
import math
import matplotlib.pyplot as plt

# Class for extracting C3' atoms from PDB files
class FindC3:
    def __init__(self, pdb_file_path):
        self.pdb_file_path = pdb_file_path
        self.c3_atoms = []

    def process(self):
        """Extracts only C3' atoms from a PDB file using fixed-width parsing."""
        with open(self.pdb_file_path, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith('ATOM') and "C3'" in line:
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    chain_id = line[21].strip()
                    res_num = line[22:26].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    self.c3_atoms.append([atom_name, res_name, chain_id, res_num, x, y, z])
        return self.c3_atoms

# Function to calculate Euclidean distance
def distance_calculator(coord1, coord2):
    coord1 = np.array(list(map(float, coord1)))
    coord2 = np.array(list(map(float, coord2)))
    return np.linalg.norm(coord1 - coord2)

# Function to classify base pairs based on distance and sequence constraints
def bp_attribution(c3_list, frequency):
    for i, res1 in enumerate(c3_list):
        for j, res2 in enumerate(c3_list[i+1:], start=i+1):
            if res1[2] == res2[2] and (int(res1[3]) + 4) <= int(res2[3]):
                bp = ''.join(sorted((res1[1] + res2[1])))
                dist = distance_calculator(res1[4:7], res2[4:7])
                if 0 <= dist <= 20:
                    int_dist = math.floor(dist)
                    for bp_tot in frequency:
                        if bp == bp_tot[0]:
                            bp_tot[1][int_dist] += 1
                            frequency[-1][1][int_dist] += 1
    return frequency

# Function to normalize observed frequencies
def normalize_frequencies(frequency):
    for bp in frequency:
        total = sum(bp[1])
        if total != 0:
            bp[1] = [val / total for val in bp[1]]
    return frequency

# Function to compute pseudo-energy values
def pseudo_energy(frequency):
    score = []
    for bp in frequency:
        pe = []
        for fobs, ref in zip(bp[1], frequency[-1][1]):
            if fobs > 0 and ref > 0:
                pseudo_score = -np.log(fobs / ref)
                pe.append(min(pseudo_score, 10))
            else:
                pe.append(10)
        score.append([bp[0], pe])
    score.pop()
    return score

# Function to generate and save interaction profile plots
def interaction_profile_plot(base_pair, scores):
    x_values = np.arange(20)
    y_values = [np.nan if score == 10 else score for score in scores]
    
    fig, ax = plt.subplots()
    ax.plot(x_values, y_values, linewidth=2.0)
    plt.title(f'Interaction Profile: {base_pair}')
    ax.set_xlabel('Distance (Å)')
    ax.set_ylabel('Pairwise Score')
    ax.set_xlim(0, 20)
    
    os.makedirs('Plot', exist_ok=True)
    plt.savefig(f'Plot/{base_pair}.png')
    plt.close()

# Main function
def main(pdb_directory):
    if not os.path.exists(pdb_directory):
        print(f"Error: The directory '{pdb_directory}' does not exist.")
        return
    
    frequency = [[bp, [0] * 20] for bp in ['AA', 'AU', 'AC', 'AG', 'UU', 'CU', 'GU', 'CC', 'CG', 'GG', 'XX']]
    
    pdb_files = [f for f in os.listdir(pdb_directory) if f.endswith(".pdb")]
    if not pdb_files:
        print(f"Error: No .pdb files found in '{pdb_directory}'.")
        return
    
    for pdb_filename in pdb_files:
        pdb_path = os.path.join(pdb_directory, pdb_filename)
        c3_atoms = FindC3(pdb_path).process()
        frequency = bp_attribution(c3_atoms, frequency)
    
    frequency = normalize_frequencies(frequency)
    pseudo_energies = pseudo_energy(frequency)
    
    for base_pair, scores in pseudo_energies:
        interaction_profile_plot(base_pair, scores)
    
    print("Plots have been saved in the 'Plot' directory.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_directory", help="Path to the directory containing PDB files.")
    args = parser.parse_args()
    main(args.pdb_directory)
