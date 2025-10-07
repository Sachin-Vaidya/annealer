import sys
import os
import json
import re
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba, to_hex
from scipy.optimize import curve_fit
from collections import defaultdict
from PIL import Image

def are_close(a, b, tol=1e-6):
    """Check if two energy values are within tolerance."""
    return abs(a - b) < tol

def group_energies(data, tol=1e-6):
    """Group energies within tolerance and count occurrences."""
    energy_counts = defaultdict(int)
    energy_to_bitstrings = defaultdict(list)
    
    for entry in data:
        if isinstance(entry, (list, tuple)) and len(entry) == 2:
            try:
                energy, bitstring = entry
                if isinstance(energy, (int, float)) and isinstance(bitstring, list):
                    matched = False
                    for existing_energy in list(energy_counts.keys()):
                        if are_close(energy, existing_energy, tol):
                            energy_counts[existing_energy] += 1
                            energy_to_bitstrings[existing_energy].append(bitstring)
                            matched = True
                            break
                    if not matched:
                        energy_counts[energy] += 1
                        energy_to_bitstrings[energy].append(bitstring)
            except (ValueError, TypeError):
                print(f"Skipping invalid entry: {entry}")
        else:
            print(f"Skipping non-conforming entry: {entry}")
    
    return energy_counts, energy_to_bitstrings

def is_valid_solution(vert, nT, nV, flag_cluster_major):
    """Check if a bitstring represents a valid solution."""
    int_vert_assignment = [-1] * nT
    for l in range(nT):
        sum_val = 0
        selected_p = -1
        for p in range(nV):
            val = int(vert[p * nT + l]) if flag_cluster_major else int(vert[p + l * nV])
            if val == 1:
                selected_p = p
                sum_val += val
        if sum_val != 1:
            return False, int_vert_assignment
        else:
            int_vert_assignment[l] = selected_p
    return True, int_vert_assignment

def create_sorted_truth_clusters(truth, event_nT):
    """Create sorted truth clusters from truth assignments."""
    int_truth = [int(val) for val in truth]
    num_tracks_assigned_to_truth_vertex = {}
    for i_internal in range(event_nT):
        vertex = int_truth[i_internal]
        if vertex not in num_tracks_assigned_to_truth_vertex:
            num_tracks_assigned_to_truth_vertex[vertex] = set()
        num_tracks_assigned_to_truth_vertex[vertex].add(i_internal)
    
    valueSetsInMap_truth = set()
    for vertex, tracks in num_tracks_assigned_to_truth_vertex.items():
        valueSetsInMap_truth.add(frozenset(tracks))
    
    sorted_truth_clusters = []
    for cluster in valueSetsInMap_truth:
        tmp = list(cluster)
        tmp.sort()
        sorted_truth_clusters.append(tmp)
    sorted_truth_clusters.sort()
    
    return sorted_truth_clusters

def is_correct_solution(int_vert_assignment, nV, sorted_truth_clusters):
    """Check if the vertex assignment matches the truth clusters."""
    clusters = [[] for _ in range(nV)]
    for t in range(len(int_vert_assignment)):
        if int_vert_assignment[t] >= 0:  # Ensure valid assignment
            clusters[int_vert_assignment[t]].append(t)
    
    filtered_clusters = [c for c in clusters if c]
    for c in filtered_clusters:
        c.sort()
    filtered_clusters.sort()
    
    return filtered_clusters == sorted_truth_clusters

def main(file_path=None, json_data=None, nV=None, nT=None, sorted_truth_clusters=None, flag_cluster_major=False):
    """Process JSON data and compute valid/correct solutions."""
    if file_path:
        with open(file_path, 'r') as f:
            data = json.load(f)
    elif json_data:
        data = json.loads(json_data)
    else:
        raise ValueError("Either file_path or json_data must be provided")
    
    energy_counts, energy_to_bitstrings = group_energies(data)
    
    count_tot = 0
    valid_solutions = defaultdict(int)
    correct_solutions = defaultdict(int)
    
    if nV is not None and nT is not None and sorted_truth_clusters is not None:
        for energy, bitstrings in energy_to_bitstrings.items():
            for bitstring in bitstrings:
                is_valid, int_vert_assignment = is_valid_solution(bitstring, nT, nV, flag_cluster_major)
                if is_valid:
                    valid_solutions[energy] += 1
                    if is_correct_solution(int_vert_assignment, nV, sorted_truth_clusters):
                        correct_solutions[energy] += 1
    
    total_valid = sum(valid_solutions.values())
    total_correct = sum(correct_solutions.values())
    
    for energy in sorted(energy_counts.keys()):
        count = energy_counts[energy]
        count_tot += count
        valid_count = valid_solutions.get(energy, 0)
        correct_count = correct_solutions.get(energy, 0)
        #print(f"Energy: {energy:.6f}, Count: {count}, Valid: {valid_count}, Correct: {correct_count}")
    
    #print(f"Total Count: {count_tot}")
    #print(f"Total Valid Solutions: {total_valid}")
    #print(f"Total Correct Solutions: {total_correct}")
    
    return count_tot, total_valid, total_correct

if __name__ == "__main__":
    # For testing in Jupyter or environments without command-line arguments
    arg = sys.argv[2]
    match = re.match(r"(\d+)Vertices_(\d+)Tracks", arg)
    if match:
        nV = int(match.group(1))
        nT = int(match.group(2))
        #print(f"Vertices: {nV}, Tracks: {nT}")
    else:
        #print("Invalid format. Expected format: '<num>Vertices_<num>Tracks'")
        sys.exit(1)
    filename_prefix = sys.argv[1]
    filename_extension = sys.argv[3]
        
    num_files = 100
    
    base_dir = filename_prefix.split('/')[0]
    output_file = f"{base_dir}/QA_valid_correctvalid_efficiency.txt"
    with open(output_file, 'w') as f_out:
        pass
    
    for i_file in range(num_files):
        file_dunn = f"{filename_prefix}{i_file+1}/dunnIndex.json"
        with open(file_dunn) as f_dunn:
            dunn = json.load(f_dunn)
        
        file_path = f"{filename_prefix}{i_file+1}{filename_extension}"
        sorted_truth_clusters = [int(i_int / (nT // nV)) for i_int in range(nT)]
        sorted_truth_clusters = create_sorted_truth_clusters(sorted_truth_clusters, nT)
        
        try:
            #print(f"\nProcessing file #{i_file+1}: {file_path}")
            count_tot, total_valid, total_correct = main(file_path=file_path, nV=nV, nT=nT, sorted_truth_clusters=sorted_truth_clusters)
            print(f"Dunn Index: {dunn}, Total samples: {count_tot}, Number of Valid solutions: {total_valid}, Number of Correct|Valid solutions: {total_correct}")
        except FileNotFoundError:
            print(f"Error: File '{file_path}' not found. Please check the file path.")
            
        with open(output_file, 'a') as f_out:
            valid_ratio = total_valid / count_tot if count_tot > 0 else 0
            correct_ratio = total_correct / total_valid if total_valid > 0 else 0
            f_out.write(f"{dunn}\t{valid_ratio}\t{correct_ratio}\n")