#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Taxonomic Distribution and Detection Method Analysis.

This script analyzes the distribution of Pacifastin homologs across phyla
and visualizes the detection bias (Sequence-based vs. Structure-based)
as shown in Figure 1C.

Usage:
    python analisis_taxonomico.py
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- Configuration ---
INPUT_FILE = "data/Pacifastin_Full_Dataset.csv" # Columns: 'Phylum', 'Method' (HMM/Foldseek/Common)
OUTPUT_FILE = "results/Figure1C_Taxonomy_Distribution.png"

def plot_taxonomic_bias(df):
    # 1. Aggegate data
    # Count homologs per Phylum and Detection Method
    counts = df.groupby(['Phylum', 'Method']).size().reset_index(name='Count')
    
    # Pivot for stacked bar chart format
    pivot_df = counts.pivot(index='Phylum', columns='Method', values='Count').fillna(0)
    
    # Sort by total abundance
    pivot_df['Total'] = pivot_df.sum(axis=1)
    pivot_df = pivot_df.sort_values('Total', ascending=False).drop(columns='Total')
    
    # 2. Plotting (Horizontal Bar Chart)
    ax = pivot_df.plot(
        kind='barh', 
        stacked=True, 
        figsize=(10, 8),
        color={'Common': 'gray', 'Sequence-only': '#d62728', 'Structure-only': '#1f77b4'},
        edgecolor='black'
    )
    
    ax.set_xscale('log') # Log scale as noted in the paper
    plt.title("Stratification of Detection Methods by Taxonomic Group", fontsize=14)
    plt.xlabel("Number of Proteins (Log Scale)", fontsize=12)
    plt.ylabel("Taxonomic Group", fontsize=12)
    plt.legend(title='Detection Method')
    plt.grid(axis='x', linestyle='--', alpha=0.5, which='both')
    
    plt.tight_layout()
    plt.savefig(OUTPUT_FILE, dpi=300)
    print(f"Chart saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    # Example loading
    # df = pd.read_csv(INPUT_FILE)
    # plot_taxonomic_bias(df)
    pass