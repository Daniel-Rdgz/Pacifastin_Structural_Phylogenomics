#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
High-Throughput Structural Prediction via ESMFold API.

This script automates the retrieval of predicted 3D structures for Pacifastin
candidates using the ESMFold API. It includes error handling, batch processing,
and automatic retries for robust large-scale modeling.

Note: 
- Requires an active internet connection.
- Rate limits may apply depending on the API provider.
"""

import requests
import time
import os
from Bio import SeqIO

# --- Configuration ---
INPUT_FASTA = "data/candidates_for_modeling.fasta"
OUTPUT_DIR = "data/pdb_models/"
API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MAX_RETRIES = 3
SLEEP_TIME = 2  # Seconds between requests to avoid rate limiting

def predict_structure(sequence, header):
    """Sends a sequence to ESMFold API and returns the PDB content."""
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.post(API_URL, data=sequence, timeout=30)
            if response.status_code == 200:
                return response.text
            else:
                print(f"  [Error] API returned {response.status_code} for {header}. Retrying...")
                time.sleep(SLEEP_TIME * (attempt + 1)) # Exponential backoff
        except Exception as e:
            print(f"  [Exception] {e}. Retrying...")
            time.sleep(SLEEP_TIME)
    return None

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    print(f"Reading sequences from {INPUT_FASTA}...")
    sequences = list(SeqIO.parse(INPUT_FASTA, "fasta"))
    total = len(sequences)
    
    for i, record in enumerate(sequences):
        seq_id = record.id
        output_file = os.path.join(OUTPUT_DIR, f"{seq_id}.pdb")
        
        # Skip if already exists (Resume capability)
        if os.path.exists(output_file):
            print(f"[{i+1}/{total}] Skipping {seq_id} (Already exists)")
            continue
            
        print(f"[{i+1}/{total}] Modeling {seq_id}...")
        pdb_content = predict_structure(str(record.seq), seq_id)
        
        if pdb_content:
            with open(output_file, "w") as f:
                f.write(pdb_content)
        else:
            with open("failed_sequences.log", "a") as log:
                log.write(f"{seq_id}\n")
            print(f"  [FAILED] Could not model {seq_id} after {MAX_RETRIES} attempts.")
            
        time.sleep(SLEEP_TIME)

if __name__ == "__main__":
    main()