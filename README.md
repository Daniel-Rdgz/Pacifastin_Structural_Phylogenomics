# Structure-guided discovery of divergent homologs unveils deep ancestry and arthropod specialization of the pacifastin family

This repository contains the custom Python scripts and data analysis pipelines used in the manuscript submitted to **Molecular Biology and Evolution (MBE)**.

## Overview
By combining sequence-based methods (HMMER) with structural alignment (Foldseek/AlphaFold), we identified a "cryptic" proteome of Pacifastin inhibitors undetectable by sequence alone. This codebase reproduces the key analytical steps of the study:

* **Iterative Mining:** Integration of HMM profiles and structural hits.
* **Architectural Classification:** Algorithms to distinguish between "Compact-loop" (Conventional) and "Extended-loop" (Ancestral) topologies.
* **Evolutionary Metrics:** Calculation of Phyletic Spread ($S_i$) and Depth ($D_i$).
* **GLM Analysis:** Logistic regression models to predict proteolytic processing sites in inter-domain linkers.

## Contents
* `src/`: Python scripts for structural mining and classification.
* `analysis/`: R/Python notebooks for the GLM and metric calculations.
* `data/`: (Optional) Sample datasets or small input files.

**Note:** The full dataset of 394 predicted structural models (PDBs) is available at Zenodo: 10.5281/zenodo.18616635.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
