#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of Evolutionary Metrics (Phyletic Spread & Depth)
for the Pacifastin Family.

This script calculates the Si (Spread) and Di (Depth) metrics as defined
in the manuscript and generates the stratified scatter plot (Figure 3).

Usage:
    python analisis_Si_Di_estratificado.py --input data/taxonomic_data.csv
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# --- Configuration ---
# Update these paths relative to your repository structure
INPUT_FILE = "data/Pacifastin_Taxonomy_Validated.csv"  # Expected columns: 'Sequence_ID', 'Lineage', 'Phylum', 'Species'
OUTPUT_DIR = "results/evolutionary_metrics/"

def calculate_metrics(df):
    """Calculates Si and Di metrics stratified by Lineage and Taxonomic Group."""
    metrics = []
    
    # Stratify by Supergroup (Arthropoda vs Non-Arthropoda)
    df['Supergroup'] = df['Phylum'].apply(lambda x: 'Arthropoda' if x == 'Arthropoda' else 'Non-Arthropoda')
    
    groups = df.groupby(['Lineage', 'Supergroup'])
    
    for (lineage, supergroup), group_data in groups:
        n_sequences = len(group_data)
        n_species = group_data['Species'].nunique()
        n_phyla = group_data['Phylum'].nunique()
        
        # Calculate Metrics
        Si = n_phyla  # Phyletic Spread (Number of unique Phyla)
        Di = n_sequences / n_species if n_species > 0 else 0 # Phyletic Depth (Avg copies per species)
        
        metrics.append({
            'Lineage': lineage,
            'Supergroup': supergroup,
            'Si': Si,
            'Di': Di,
            'N_Seq': n_sequences,
            'N_Species': n_species
        })
        
    return pd.DataFrame(metrics)

def plot_landscape(metrics_df):
    """Generates the evolutionary landscape plot (Si vs Di)."""
    plt.figure(figsize=(10, 8))
    sns.set_style("whitegrid")
    
    # Custom markers and colors
    markers = {"Extended-loop": "s", "Compact-loop": "o"} # s=square, o=circle
    palette = {"Arthropoda": "#d62728", "Non-Arthropoda": "#1f77b4"}
    
    sns.scatterplot(
        data=metrics_df, 
        x='Si', y='Di',
        hue='Supergroup', style='Lineage',
        size='N_Seq', sizes=(100, 1000),
        markers=markers, palette=palette,
        alpha=0.8, edgecolor="black"
    )
    
    plt.title("Stratified Evolutionary Landscape of Pacifastin Lineages", fontsize=14)
    plt.xlabel("Phyletic Spread ($S_i$) [Number of Phyla]", fontsize=12)
    plt.ylabel("Phyletic Depth ($D_i$) [Avg. Copies per Species]", fontsize=12)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    
    # Save
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    plt.savefig(f"{OUTPUT_DIR}/Figure3_EvolutionaryLandscape.png", dpi=300)
    print(f"Plot saved to {OUTPUT_DIR}")

if __name__ == "__main__":
    print("Loading data...")
    # Dummy data generation for testing (Replace with pd.read_csv(INPUT_FILE))
    # df = pd.read_csv(INPUT_FILE)
    print("Please ensure your input CSV has columns: Sequence_ID, Lineage (Extended-loop/Compact-loop), Phylum, Species")