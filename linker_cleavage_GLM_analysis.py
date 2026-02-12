#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generalized Linear Model (GLM) Analysis of Linker Cleavage Probability.

This script implements a logistic regression model to predict the presence
of dibasic cleavage sites (KR, RR, RK, KK) in inter-domain linkers as a
function of linker length and taxonomic lineage. It generates the Odds Ratios
(OR) presented in Figure 7.

Usage:
    python linker_cleavage_GLM_analysis.py
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import seaborn as sns

# --- Configuration ---
INPUT_FILE = "data/Linker_Sequences_Validated.csv" # Columns: 'Sequence_ID', 'Linker_Seq', 'Lineage', 'Phylum'
OUTPUT_DIR = "results/glm_analysis/"

def detect_cleavage_sites(seq):
    """Detects canonical dibasic motifs (KR, RR, RK, KK)."""
    motifs = ['KR', 'RR', 'RK', 'KK']
    for motif in motifs:
        if motif in seq:
            return 1
    return 0

def run_glm_analysis():
    print("Loading data...")
    # df = pd.read_csv(INPUT_FILE)  # Uncomment when using real data
    
    # --- Mock Data Generation (Remove this block in production) ---
    np.random.seed(42)
    n = 200
    df = pd.DataFrame({
        'Linker_Seq': ['A'*np.random.randint(5, 50) for _ in range(n)],
        'Lineage': np.random.choice(['Extended-loop', 'Compact-loop'], n),
        'Phylum': np.random.choice(['Arthropoda', 'Mollusca', 'Cnidaria'], n)
    })
    # -------------------------------------------------------------

    # 1. Feature Engineering
    df['Length'] = df['Linker_Seq'].apply(len)
    df['Has_Cleavage'] = df['Linker_Seq'].apply(detect_cleavage_sites)
    
    # 2. GLM: Logistic Regression
    # Formula: Probability(Cleavage) ~ Length + Lineage
    print("Fitting GLM (Binomial Family)...")
    model = smf.glm(
        formula="Has_Cleavage ~ Length + Lineage", 
        data=df, 
        family=sm.families.Binomial()
    ).fit()
    
    print(model.summary())
    
    # Extract Odds Ratios
    params = model.params
    conf = model.conf_int()
    conf['OR'] = params
    conf.columns = ['Lower CI', 'Upper CI', 'OR']
    print("\nOdds Ratios:")
    print(np.exp(conf))

    # 3. Visualization (Sigmoid Curves)
    plt.figure(figsize=(8, 6))
    sns.regplot(
        x="Length", y="Has_Cleavage", 
        data=df[df['Lineage']=='Compact-loop'], 
        logistic=True, ci=None, label='Compact-loop', color='red'
    )
    sns.regplot(
        x="Length", y="Has_Cleavage", 
        data=df[df['Lineage']=='Extended-loop'], 
        logistic=True, ci=None, label='Extended-loop', color='blue'
    )
    
    plt.title("Probability of Processing Site Occurrence vs Length", fontsize=14)
    plt.xlabel("Linker Length (Amino Acids)", fontsize=12)
    plt.ylabel("Probability (P)", fontsize=12)
    plt.legend()
    plt.tight_layout()
    plt.savefig("results/Figure7A_GLM_Curves.pdf")
    print("Plot saved.")

if __name__ == "__main__":
    import os
    if not os.path.exists("results/"):
        os.makedirs("results/")
    run_glm_analysis()