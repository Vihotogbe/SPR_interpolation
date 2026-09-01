# Spatial Pattern Regression (SPR) for Meteorological Fields Interpolation

This repository contains the R code and reproducibility materials for the article:

> **Houssou, V. and Carreau, J.** (2026). *Spatial pattern regression for meteorological fields interpolation*. Hydrology and Earth System Sciences (HESS), https://doi.org/10.5194/egusphere-2026-1702, Copernicus Publications.

---

## Overview

**Spatial Pattern Regression (SPR)** is a data-driven interpolation method that reconstructs high-resolution gridded meteorological fields from sparse station observations. SPR combines:

- **Spatial patterns** extracted from Regional Climate Model (RCM) simulations via Principal Component Analysis (PCA/SVD)
- **Daily regression** of sparse station observations onto these patterns using Ordinary Least Squares (OLS)

SPR is designed for hydrological impact studies, particularly in data-sparse regions, and ensures spatial consistency between interpolated historical fields and future RCM projections.

---

## Method

SPR operates in two steps:

1. **Pattern extraction** — Spatial patterns (EOFs) are extracted from the ClimEx RCM dataset over an auxiliary period (1980–2009) using SVD/PCA. These patterns form a fixed orthonormal basis that does not vary in time.

2. **Daily regression** — For each day of the interpolation period, the amplitudes of the spatial patterns are estimated by OLS regression on available station observations. The reconstructed field is then obtained as a linear combination of the patterns.

The method is formulated as:

```
Z_d = Z̄_grid + V_grid * β_k + ε_d
```

where `Z̄_grid` is the temporal column mean of the RCM data matrix (a single constant per grid cell, distinct from monthly climatologies), `V_grid` contains the spatial patterns, and `β_k` are the daily regression coefficients estimated from observations.

---

## Data

### RCM simulations

- **Source**: [ClimEx project](https://climex-data.srv.lrz.de/Public/) — Canadian RCM (CRCM5) driven by CanESM2
- **Ensemble member**: kdj
- **Domain**: North American domain, 280×280 grid cells, ~11 km resolution
- **Auxiliary period**: 1980–2009 (n = 10,950 daily fields, standard 365-day calendar)
- **Interpolation period**: 2000–2009

### Station observations

- **Source**: Environment and Climate Change Canada (ECCC) — [Bulk data API](https://climate.weather.gc.ca/climate_data/bulk_data_e.html)
- **Domain**: Southern Quebec, ~70,000 km², latitudes 45.0–47.0°N, longitudes 75.0–70.0°W
- **Variables and periods**:
  - Precipitation: 15 stations, October 8, 2000 – October 7, 2002
  - Minimum temperature: 44 stations, December 7, 2000 – December 6, 2002
  - Maximum temperature: 44 stations, December 12, 2016 – December 11, 2018
- **Selection criterion**: Strictly continuous records with no missing values and no gap filling

---

## Experiments

### Synthetic experiments

Virtual station networks are defined as random subsets of RCM grid cells, allowing controlled evaluation under known ground truth.

**Factorial design**:

| Factor | Values |
|---|---|
| Region location | South (around Montreal), North (northern Quebec) |
| Region size | Large (S1), Medium (S2), Small (S3) |
| Network density | 10%, 30%, 50%, 70%, 90% |
| Variables | Precipitation, Tmin, Tmax |

This yields **30 experiments per variable** (2 locations × 3 sizes × 5 densities).

**Grid cell counts**:

| Size | Region South | Region North |
|---|---|---|
| Large (S1) | 3,025 | 2,970 |
| Medium (S2) | 1,332 | 1,332 |
| Small (S3) | 342 | 342 |

### Stress-test experiment

Extremely sparse network: 3 virtual stations (~0.1% density) over the large northern region (2,970 grid cells). Hyperparameters are fixed to values obtained at 10% density.

### Real data experiments

SPR is applied to real ECCC station observations. Station density is set to 10% and 30%. To assess robustness to station sampling, station selection is **repeated 100 times** with independent random seeds. Performance metrics are averaged across repetitions.

---

## Baseline methods

SPR is compared against three standard spatial interpolation methods:

| Method | Description |
|---|---|
| **IDW** | Inverse Distance Weighting — deterministic, distance-based weights |
| **OK** | Ordinary Kriging — geostatistical, variogram-based spatial correlation |
| **KED** | Kriging with External Drift — OK with RCM monthly climatologies as auxiliary drift |

---

## Hyperparameter selection

Hyperparameters are selected by maximizing average RMSE on a held-out validation set over the full interpolation period. A single globally optimal configuration is retained per experimental setup.

| Method | Hyperparameters |
|---|---|
| IDW | Distance power (1–5) |
| OK | Variogram model (Gaussian, Spherical, Exponential) |
| KED | Variogram model (Gaussian, Spherical, Exponential) |
| SPR | Number of retained patterns k (10%–90% in steps of 5%) |

**SPR-specific notes**:
- The selected k is fixed once per configuration (region × size × density × variable) and applied uniformly to all days
- For precipitation, optimal k is typically 40%–50% of available patterns under sparse conditions
- For temperature, optimal k is typically 25%–35%, reflecting smoother spatial variability

---

## Performance metrics

| Metric | Description |
|---|---|
| **RMSE** | Root Mean Squared Error — pointwise accuracy, captures bias |
| **SSIM** | Structural Similarity Index — spatial structure, captures variance and correlation |

Together RMSE and SSIM cover the three aspects of the Kling-Gupta Efficiency (bias, variability, correlation).

For precipitation, a **softplus transformation** `f(x) = log(exp(x) - 1)` is applied prior to interpolation to enforce non-negativity. Its inverse `f⁻¹(x) = log(1 + eˣ) ≈ x` for practical precipitation values, introducing negligible back-transformation bias.

---

## Repository structure

```
.
├── README.md
├── data/
│   ├── coord.csv                          # RCM grid coordinates
│   └── montecarlo_spr_ked_results.csv     # Real data experiment results (100 repetitions)
├── R/
│   ├── 00_data_processing.R               # Data loading, CRS setup, region definition
│   ├── 01_synthetic_experiments.R         # Synthetic experiment loop
│   ├── 02_real_data_experiments.R         # Real data experiments with 100 repetitions
│   ├── 03_hyperparameter_selection.R      # Hyperparameter grid search
│   ├── 04_figures_synthetic.R             # RMSE/SSIM scatter plots (Fig. 3–6)
│   ├── 05_figures_stress_test.R           # Stress-test field maps (Fig. 7–8)
│   ├── 06_figures_real_data.R             # Boxplots from 100 repetitions (Fig. 9)
│   └── 07_figures_map.R                   # Study region map (Fig. 1)
└── figures/
    └── *.png                              # All manuscript figures
```

---

## Requirements

```r
# R version 4.4.3
# Main packages
library(bigstatsr)    # Large matrix SVD
library(gstat)        # Kriging (OK and KED)
library(sf)           # Spatial data handling
library(rnaturalearth) # Base maps
library(ggplot2)      # Visualization
library(tidyr)        # Data reshaping
library(dplyr)        # Data manipulation
```

---

## Code availability

Analysis scripts are available from the corresponding author upon reasonable request.

---

## Citation

If you use this code, please cite the associated manuscript:

```bibtex
@article{houssou2026spr,
  title     = {Spatial pattern regression for meteorological fields interpolation},
  author    = {Houssou, Vihotogb{\'e} and Carreau, Julie},
  journal   = {Hydrology and Earth System Sciences},
  publisher = {Copernicus Publications},
  year      = {2026},
  doi       = {https://doi.org/10.5194/egusphere-2026-1702}
}
```
And the archived code itself:

```bibtex

```

---

## Funding

This work was supported by the Natural Sciences and Engineering Research Council of Canada (NSERC), the Fonds de Recherche du Québec — Nature et Technologies (FRQNT), and IVADO.

---

## Contact

**Vihotogbé Houssou** — vihotogbe-2.houssou@polymtl.ca  
PhD candidate, Department of Mathematics and Industrial Engineering  
Polytechnique Montréal | GERAD | IVADO
