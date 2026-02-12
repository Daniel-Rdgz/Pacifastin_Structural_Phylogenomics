#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Structural Landscape Analysis via Multidimensional Scaling (MDS).

This script performs MDS on the pairwise RMSD matrix of the Pacifastin domains
to visualize the conformational space (Figure 4). It also calculates the
Structural Dispersion Index (SDI) for the Compact vs. Extended lineages.

Usage:
    python analisis_estructural_MDS_CE.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import MDS

# --- Configuration ---
INPUT_RMSD_MATRIX = "data/RMSD_Matrix_All_vs_All.csv" # Pre-computed RMSD matrix
METADATA_FILE = "data/Pacifastin_Metadata.csv"        # Must contain 'Sequence_ID', 'Lineage'
OUTPUT_DIR = "results/structural_landscape/"

def calculate_dispersion(coords_df, lineage_name):
    """Calculates the mean Euclidean distance to the centroid (Dispersion Index)."""
    subset = coords_df[coords_df['Lineage'] == lineage_name]
    if subset.empty:
        return 0.0
    centroid = subset[['Dim1', 'Dim2']].mean()
    distances = np.sqrt(((subset[['Dim1', 'Dim2']] - centroid) ** 2).sum(axis=1))
    return distances.mean()

def run_mds_analysis():
    print("Loading data...")
    # 1. Load RMSD Matrix
    try:
        rmsd_df = pd.read_csv(INPUT_RMSD_MATRIX, index_col=0)
        metadata = pd.read_csv(METADATA_FILE)
    except FileNotFoundError:
        print(f"Error: Input files not found. Please check {INPUT_RMSD_MATRIX} and {METADATA_FILE}")
        return

    # 2. Run MDS (Metric Multidimensional Scaling)
    print("Running MDS projection...")
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42, n_init=4)
    coords = mds.fit_transform(rmsd_df)
    
    # 3. Prepare Data for Plotting
    coords_df = pd.DataFrame(coords, columns=['Dim1', 'Dim2'], index=rmsd_df.index)
    coords_df = coords_df.merge(metadata, left_index=True, right_on='Sequence_ID')
    
    # Calculate Dispersion Indices
    di_extended = calculate_dispersion(coords_df, "Extended-loop")
    di_compact = calculate_dispersion(coords_df, "Compact-loop")
    print(f"Dispersion Index (Extended): {di_extended:.2f}")
    print(f"Dispersion Index (Compact): {di_compact:.2f}")

    # 4. Plotting
    plt.figure(figsize=(10, 8))
    sns.set_style("white")
    
    # KDE Density plot for background
    sns.kdeplot(
        data=coords_df, x='Dim1', y='Dim2', hue='Lineage',
        fill=True, alpha=0.1, levels=5, thresh=0.2, palette=['blue', 'red']
    )
    
    # Scatter plot
    sns.scatterplot(
        data=coords_df, x='Dim1', y='Dim2', hue='Lineage', style='Lineage',
        palette={'Extended-loop': 'blue', 'Compact-loop': 'red'},
        s=60, alpha=0.8, edgecolor='k'
    )
    
    plt.title("Evolutionary Structural Landscape (Backbone RMSD)", fontsize=14)
    plt.xlabel("Conformational Space PC1", fontsize=12)
    plt.ylabel("Conformational Space PC2", fontsize=12)
    plt.legend(title='Architecture')
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/Figure4_MDS_Landscape.pdf")
    print(f"Plot saved to {OUTPUT_DIR}")

if __name__ == "__main__":
    import os
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    run_mds_analysis()