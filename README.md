# rEGA: Reliability Exploratory Graph Analysis

[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Open Science: OSF Ready](https://img.shields.io/badge/Open%20Science-SHA256%20Audited-orange.svg)](https://osf.io/)

**`rEGA`** (*Reliability Exploratory Graph Analysis*) is an R package for network psychometrics that bridges Exploratory Graph Analysis (`EGA`; Golino & Epskamp, 2017; Golino et al., 2020) with exact measurement error decomposition and score reliability estimation directly from Gaussian Graphical Models (GGMs) and Markov Random Fields.

---

## Methodological Overview

Standard EGA pipelines estimate precision matrices **K = Sigma^-1**, extract partial correlations **Omega**, and discard the conditional precision diagonal **D_K = diag(K)**, implicitly treating all nodes as equally free of conditional error. **`rEGA`** resolves this limitation through six sequential innovations:

1. **Support-Constrained MLE Refitting (`K_refit`)**: Eliminates L1 (LASSO) shrinkage bias on active edges identified by `EBICglasso`.
2. **Neumann-Markov Resolvent Decomposition**: Partitions the reproduced covariance matrix via the resolvent operator **Q = (I - Omega)^-1**.
3. **Nodal Variance Conservation Law**: Decomposes each item's unit variance into three interpretable parameters summing to 100% of item variance:
   > **h2** *(Within-Community Signal)* + **beta2** *(Bridging Variance)* + **psi2** *(Markov Uniqueness)* = **1 (100%)**
4. **Purity-Weighted Community Detection (`A^rEGA`)**: Weights edges by topological purity (**h_G^2 = 1 - psi2**) prior to Walktrap, Leiden, or Louvain partitioning.
5. **Network Reliability Coefficient Gamma & Graph SEM**: Provides closed-form dimensional and total reliability coefficients (**Gamma_c, Gamma_t**) alongside standardized and raw-score Standard Errors of Measurement (**SEM_c**).
6. **Holistic Interpretation (No Rigid Cutoffs)**: Because `h2 + beta2 + psi2 = 1` forms a compositional conservation law, `rEGA` and `boot.rEGA` evaluate nodal variance components, structural reliability, and item replicability holistically without imposing arbitrary fixed thresholds.

---

## Installation

Install the package directly from GitHub:

```r
# If remotes is not installed: install.packages("remotes")
remotes::install_github("Huerta-Huacchahuilca/rEGA")
```

---

## Comprehensive User Guide

### 1. Model Estimation (`rEGA`)

`rEGA()` estimates the regularized network from raw item responses (or a correlation matrix with `n`), supports flexible item subsetting via `items`, and displays a side-by-side visualization combining the **EGA-style Network** and the **Nodal Variance Triad (with exact percentages summing to 100% per item)**:

```r
library(rEGA)

# Simulated 2-factor structure (6 items, N = 300)
set.seed(2026)
Sigma <- matrix(0.30, 6, 6)
Sigma[1:3, 1:3] <- 0.70
Sigma[4:6, 4:6] <- 0.70
diag(Sigma) <- 1.0
data_items <- as.data.frame(MASS::mvrnorm(n = 300, mu = rep(0, 6), Sigma = Sigma))
colnames(data_items) <- paste0("Item_", 1:6)

# Fit the rEGA model (displays EGA-style network + 100% Nodal Variance Triad)
fit <- rEGA(
  data       = data_items,
  items      = 1:6,
  cor.method = "cor_auto",
  algorithm  = "walktrap",
  plot.rEGA  = TRUE,
  seed       = 2026
)

print(fit)
```

### 2. Visualizing the Network and Standardized Nodal Percentages (`plot.rEGA`)

```r
# Combined figure (EGA-style Network + Nodal Variance Triad in %)
plot(fit, type = "both", show.values = TRUE)

# Only EGA-style Support-Refitted Network
plot(fit, type = "network")

# Only Nodal Variance Triad showing exact % (h2% + beta2% + psi2% = 100% per item)
plot(fit, type = "variance", show.values = TRUE)
```

### 3. Parallel Bootstrap Inference and Unified Stability Plot (`boot.rEGA`)

`boot.rEGA()` performs parametric or non-parametric resampling and generates a **unified side-by-side figure** combining the empirical **rEGA Network** and the **Bootstrap Item Stability Plot** (similar to `bootEGA`):

```r
boot_fit <- boot.rEGA(
  rEGA.obj       = fit,
  iter           = 500,
  type           = "parametric",
  ncores         = 2,
  seed           = 2026,
  plot.boot.rEGA = TRUE
)

print(boot_fit)

# Plot options for boot.rEGA:
plot(boot_fit, type = "both")      # Unified rEGA Network + Item Stability
plot(boot_fit, type = "all")       # 3-Panel: Network + Nodal Triad (%) + Stability
plot(boot_fit, type = "stability") # Item Stability panel only
```

### 4. Attenuation-Aware Power Planning & Open Science Audit (`rEGA.power` & `rEGA.audit`)

```r
# Sample size planning accounting for network measurement error (Henry, 2026)
rEGA.power(rEGA.obj = fit, p.target = c(6, 10, 14), r.true = 0.30)

# Export multi-sheet Excel audit workbook with SHA-256 cryptographic data hash
rEGA.audit(rEGA.obj = fit, boot.obj = boot_fit, file = "rEGA_Audit_Report.xlsx")
```

---

## Author & Maintainer

* **César Huerta-Huacchahuilca**  
  *Creator, Lead Developer & Maintainer*  
  *Email:* [andrehuerta1411@gmail.com](mailto:andrehuerta1411@gmail.com)  
  *GitHub Profile:* [@Huerta-Huacchahuilca](https://github.com/Huerta-Huacchahuilca)  
  *ORCID:* [0009-0000-8906-2678](https://orcid.org/0009-0000-8906-2678)  
  *Issue Tracker:* [rEGA Issues](https://github.com/Huerta-Huacchahuilca/rEGA/issues)

---

## Citation & Academic Attribution

If you use **`rEGA`** in academic studies, thesis projects, or methodological research, please cite the package as:

### APA 7th Edition

> Huerta-Huacchahuilca, C. (2026). *rEGA: Reliability Exploratory Graph Analysis* (R package version 0.1.0) [Computer software]. GitHub. https://github.com/Huerta-Huacchahuilca/rEGA

### BibTeX Entry

```bibtex
@manual{huerta2026rega,
  title  = {rEGA: Reliability Exploratory Graph Analysis},
  author = {Huerta-Huacchahuilca, C{\'e}sar},
  year   = {2026},
  note   = {R package version 0.1.0},
  url    = {https://github.com/Huerta-Huacchahuilca/rEGA}
}
```

---

## License

This package is free and open-source software licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.
