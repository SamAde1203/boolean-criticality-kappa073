# Phase Transitions in Biological Organisation: The κ=0.73 Threshold

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-3100/)
[![R 4.3](https://img.shields.io/badge/R-4.3-blue.svg)](https://www.r-project.org/)

## Overview

This repository contains all code, data, and supplementary materials for the manuscript:

> **Adeyemi, S.** (2026). Phase Transitions in Biological Organisation: The κ=0.73 Threshold.
> *[Journal Name]*. DOI: [to be assigned upon publication]

### Abstract

Biological systems across scales — from gene regulatory networks to neuronal assemblies —
exhibit critical dynamics that position them at phase transition boundaries between order
and chaos. This work identifies **κ = 0.73 ± 0.02** as a universal threshold for effective
connectivity demarcating critical behaviour in Boolean network models of biological
organisation. Through analytical derivation, large-scale computational simulation
(N=10,000 networks, n=1,000 nodes), and empirical analysis of 48 biological regulatory
networks, we demonstrate that systems near κ ≈ 0.73 exhibit maximal dynamic range,
optimal information transmission, and robust adaptability.

The critical value emerges from first principles as:

κ_c^bio = β × K_c^hom

where β = ⟨k⟩²/⟨k²⟩ ≈ 0.36 for biological power-law networks (γ=2.3),
and K_c^hom = 1/[2p(1−p)] ≈ 2.02 for bias p=0.55.

---

## Repository Structure

boolean-criticality-kappa073/
├── code/python/ # Python simulation scripts
├── code/R/ # R analysis and visualisation scripts
├── data/ # Empirical network data and simulation outputs
├── figures/ # Generated manuscript figures
├── supplementary/ # Supplementary notes and tables
└── manuscript/ # Manuscript-related files

text

---

## Quick Start

### 1. Clone the Repository
```bash
git clone https://github.com/SamAdeyemi/boolean-criticality-kappa073.git
cd boolean-criticality-kappa073
2. Install Python Dependencies
bash
pip install -r requirements.txt
3. Install R Dependencies
r
install.packages(c("ggplot2", "igraph", "poweRlaw", "boot",
                   "dplyr", "tidyr", "gridExtra", "scales",
                   "viridis", "BoolNet"))
4. Run Full Analysis Pipeline
Step 1 — Derive β from first principles (R):

bash
Rscript code/R/00_beta_derivation.R
Step 2 — Generate Boolean networks and run simulations (Python):

bash
python code/python/01_network_generation.py
python code/python/02_dynamical_analysis.py
python code/python/03_avalanche_statistics.py
python code/python/04_functional_metrics.py
python code/python/05_evolutionary_simulation.py
Step 3 — Reproduce all manuscript figures (R):

bash
Rscript code/R/01_phase_diagram.R
Rscript code/R/02_avalanche_analysis.R
Rscript code/R/03_empirical_networks.R
Rscript code/R/04_functional_performance.R
Rscript code/R/05_robustness_analysis.R
Rscript code/R/06_supplementary_analyses.R
Key Results
Finding	Value
Critical threshold κ_c	0.73 ± 0.02
Biological networks mean κ	0.74 ± 0.09 (n=48)
Networks within Δκ < 0.15	83% (40/48)
Scaling factor β (derived)	⟨k⟩²/⟨k²⟩ = 0.358 ≈ 0.36
Peak information capacity	0.89 ± 0.08 bits/node
Peak dynamic range	4.2 ± 0.5 log units
Peak memory retention	152 ± 23 time steps
Biological vs random (Cohen's d)	1.61
Data Sources
Empirical biological Boolean networks sourced from:

BioModels — https://www.ebi.ac.uk/biomodels/

RegulonDB — http://regulondb.ccg.unam.mx/

KEGG Pathway — https://www.genome.jp/kegg/pathway.html

WormAtlas — https://www.wormatlas.org

Citation
If you use this code or data, please cite:

text
@article{adeyemi2026kappa,
  author    = {Adeyemi, Sam},
  title     = {Phase Transitions in Biological Organisation:
               The κ=0.73 Threshold},
  journal   = {[Journal Name]},
  year      = {2026},
  doi       = {10.5281/zenodo.XXXXXX}
}
Author
Sam Adeyemi
Independent Researcher, United Kingdom
samo.adeyemi@yahoo.co.uk
ORCID: [your ORCID here]

License
This project is licensed under the MIT License — see LICENSE for details.

Acknowledgements
Data contributions from BioModels, RegulonDB, KEGG Pathway, and WormAtlas.
To God be all the glory.