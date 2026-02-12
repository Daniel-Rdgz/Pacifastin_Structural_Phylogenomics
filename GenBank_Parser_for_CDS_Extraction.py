#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GenBank Parser for CDS Extraction and Metadata Retrieval.

This script processes local GenBank files (.gb/.gbk) to extract:
1. Nucleotide Coding Sequences (CDS)
2. Protein Translations
3. Taxonomic Metadata (Organism, Taxonomy lineage)

Used to build the initial sequence database before structural filtering.
"""

from Bio import SeqIO
import pandas as pd
import os

# --- Configuration ---
INPUT_DIR = "data/genbank_files/"
OUTPUT_FASTA = "data/extracted_cds.fasta"
OUTPUT_METADATA = "data/genbank_metadata.csv"

def process_genbank_files():
    records_data = []
    
    with open(OUTPUT_FASTA, "w") as fasta_out:
        for filename in os.listdir(INPUT_DIR):
            if not filename.endswith((".gb", ".gbk")):
                continue
                
            filepath = os.path.join(INPUT_DIR, filename)
            print(f"Processing {filename}...")
            
            try:
                # Parse GenBank file
                for gb_record in SeqIO.parse(filepath, "genbank"):
                    organism = gb_record.annotations.get("organism", "Unknown")
                    taxonomy = "; ".join(gb_record.annotations.get("taxonomy", []))
                    
                    for feature in gb_record.features:
                        if feature.type == "CDS":
                            # Extract identifiers
                            locus_tag = feature.qualifiers.get("locus_tag", ["N/A"])[0]
                            protein_id = feature.qualifiers.get("protein_id", ["N/A"])[0]
                            
                            # Extract sequences
                            try:
                                nucleotide_seq = feature.location.extract(gb_record).seq
                                translation = feature.qualifiers.get("translation", [""])[0]
                                
                                # Write to FASTA
                                if translation:
                                    fasta_out.write(f">{protein_id} | {organism}\n{translation}\n")
                                
                                # Store metadata
                                records_data.append({
                                    "File": filename,
                                    "Protein_ID": protein_id,
                                    "Locus_Tag": locus_tag,
                                    "Organism": organism,
                                    "Taxonomy": taxonomy
                                })
                                
                            except Exception as e:
                                print(f"  Error extracting feature in {gb_record.id}: {e}")
                                
            except Exception as e:
                print(f"Error parsing file {filename}: {e}")

    # Save Metadata CSV
    df = pd.DataFrame(records_data)
    df.to_csv(OUTPUT_METADATA, index=False)
    print(f"Done. Extracted {len(df)} CDS records.")

if __name__ == "__main__":
    if not os.path.exists(INPUT_DIR):
        print(f"Please create {INPUT_DIR} and place .gb files there.")
    else:
        process_genbank_files()