This repository contains the code for the simulation study in the manuscript:

**Nonparametric spectral density estimation using interactive mechanisms under local differential privacy** 


## Execution Order
1. Create folders: `DATA/`, `FIGURES/`, `RESULTS/`
2. Run `R/simulate_data.R` to generate data (4 available processes: EX1–EX4)
3. Run `R/estimation_plots_acf.R` for covariance estimation plots
4. Run `R/estimation_plots_sdf.R` for spectral density plots

## Requirements
- R (≥ 4.0)
- Packages: `MASS`, `dtt`, `extraDistr`

## Citation
If you use this code, please cite:

C. Butucea, K. Klockmann, and T. Krivobokova. (2026). *Nonparametric Spectral Density Estimation using Interactive Mechanisms under Local Differential Privacy*. arXiv preprint arXiv:2504.00919

available at https://arxiv.org/abs/2504.00919

This repository is linked to the ArXiv preprint (v2, 2026). 
Once the paper is published, the citation will be updated to include the journal and DOI.
