#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Structure-Anchored Sequence Logo Generation.

This script generates the sequence logos (Figure 2B) for the Compact and 
Extended lineages using the 'logomaker' library. It requires a pre-computed
multiple sequence alignment (MSA) anchored on the conserved Cysteine residues.

Usage:
    python Fig2B_SequenceLogos.py
"""

import pandas as pd
import matplotlib.pyplot as plt
import logomaker
import os

# --- Configuration ---
# Input should be a matrix where rows=sequences, cols=positions (aligned)
INPUT_MSA_EXTENDED = "data/MSA_Extended_Anchored.csv" 
INPUT_MSA_COMPACT = "data/MSA_Compact_Anchored.csv"
OUTPUT_DIR = "results/figures/"

def plot_logo(msa_file, title, output_name):
    # 1. Load Data
    # For demo purposes, we create a dummy probability matrix
    # Real usage: df = pd.read_csv(msa_file) -> Calculate probability matrix
    
    # Create a mock probability matrix (A, C, G, T for DNA or amino acids)
    # Here we use a generic amino acid alphabet for visualization
    import numpy as np
    aa = list("ACDEFGHIKLMNPQRSTVWY")
    mock_df = pd.DataFrame(np.random.rand(40, 20), columns=aa)
    # Normalize rows to sum to 1 (probability)
    mock_df = mock_df.div(mock_df.sum(axis=1), axis=0)
    
    # 2. Create Logo
    logo = logomaker.Logo(
        mock_df,
        shade_below=.5,
        fade_below=.5,
        font_name='Arial Rounded MT Bold'
    )
    
    # 3. Style
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    logo.ax.set_ylabel("Information (Bits)", fontsize=12)
    logo.ax.set_title(title, fontsize=14)
    
    # Highlight Cysteines (Example positions)
    # logo.highlight_position(p=0, color='gold', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, output_name))
    print(f"Logo saved: {output_name}")

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    print("Generating Sequence Logos...")
    plot_logo(INPUT_MSA_EXTENDED, "Ancestral Lineage (Extended C1-C2 Loop)", "Fig2B_Extended_Logo.png")
    plot_logo(INPUT_MSA_COMPACT, "Conventional Lineage (Compact C1-C2 Loop)", "Fig2B_Compact_Logo.png")

if __name__ == "__main__":
    main()