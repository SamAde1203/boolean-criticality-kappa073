# Data Directory

## Structure

data/
├── empirical_networks/ # 48 empirical Boolean network datasets
│ ├── README_empirical.md
│ ├── gene_regulatory/ # 28 GRN networks
│ ├── signalling/ # 12 signalling pathway networks
│ ├── neuronal/ # 8 C. elegans neural networks
│ └── empirical_kappa_values.csv # Summary κ_eff for all 48 networks
│
└── simulation_outputs/ # Outputs from Python simulations
├── beta_derivation_primary.csv
├── beta_sensitivity_gamma.csv
├── beta_sensitivity_kmax.csv
├── lyapunov_data.csv
├── kappa_c_bootstrap.csv
├── avalanche_data.csv
├── functional_metrics.csv
└── evolutionary_simulation.csv



## Empirical Network Sources

| Database | URL | Networks Used |
|---|---|---|
| BioModels | https://www.ebi.ac.uk/biomodels/ | 19 |
| RegulonDB | http://regulondb.ccg.unam.mx/ | 19 |
| KEGG Pathway | https://www.genome.jp/kegg/pathway.html | 10 |
| WormAtlas | https://www.wormatlas.org | 8 |

**Total screened:** 74 networks
**Excluded:** 26 (8 too small, 9 incomplete functions, 5 conflicting, 4 duplicates)
**Included:** 48 networks (see Supplementary Figure S2 for PRISMA diagram)

No network was excluded on the basis of its κ value.