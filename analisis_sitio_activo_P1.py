#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functional Analysis of the P1 Reactive Site.

This script analyzes the amino acid frequency at the P1 position
and generates a Row-Scaled Dot Plot to visualize lineage-specific
specialization (Figure 5).
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- Configuration ---
INPUT_FILE = "data/P1_Residues.csv" # Columns: 'Lineage', 'Phylum', 'P1_Residue'

def plot_dotplot(df):
    # 1. Calculate Frequencies
    counts = df.groupby(['Phylum', 'P1_Residue']).size().reset_index(name='Count')
    total_per_phylum = df.groupby('Phylum').size().reset_index(name='Total')
    data = counts.merge(total_per_phylum, on='Phylum')
    data['Frequency'] = data['Count'] / data['Total']
    
    # 2. Row-Scaling (Normalize relative to the max usage of that residue)
    max_freq_per_residue = data.groupby('P1_Residue')['Frequency'].transform('max')
    data['Relative_Dominance'] = data['Frequency'] / max_freq_per_residue
    
    # 3. Plotting
    plt.figure(figsize=(12, 6))
    
    # Map colors to physicochemical properties (optional)
    # Hydrophobic=Red, Basic=Blue, etc.
    
    sns.scatterplot(
        data=data,
        x='Phylum', y='P1_Residue',
        size='Relative_Dominance', sizes=(20, 500),
        hue='P1_Residue', # Or map to function
        legend=False, alpha=0.7
    )
    
    plt.title("Phylogenetic Specificity Landscape (Row-Scaled)", fontsize=14)
    plt.ylabel("P1 Residue", fontsize=12)
    plt.xlabel("Phylogenetic Group", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig("results/Figure5_P1_Landscape.pdf")

if __name__ == "__main__":
    # Load your data here
    pass