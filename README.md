# rEGA: Reliability Exploratory Graph Analysis <img src="man/figures/logo.png" align="right" height="120" alt="" />

[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Open Science: OSF Ready](https://img.shields.io/badge/Open%20Science-SHA256%20Audited-orange.svg)](https://osf.io/)

**`rEGA`** (*Reliability Exploratory Graph Analysis*) is an R package for network psychometrics that bridges Exploratory Graph Analysis (`EGA`; Golino & Epskamp, 2017; Golino et al., 2020) with exact measurement error decomposition and score reliability estimation directly from Gaussian Graphical Models (GGMs) and Markov Random Fields.

---

## Methodological Overview

Standard EGA pipelines estimate precision matrices $\mathbf{K} = \mathbf{\Sigma}^{-1}$, extract partial correlations $\mathbf{\Omega}$, and discard the conditional precision diagonal $\mathbf{D}_K = \text{diag}(\mathbf{K})$, implicitly treating all nodes as equally free of conditional error. **`rEGA`** resolves this limitation through six sequential innovations:

1. **Support-Constrained MLE Refitting ($\hat{\mathbf{K}}_{\text{refit}}$)**: Eliminates $\ell_1$ (LASSO) shrinkage bias on active edges $E_\lambda$ identified by `EBICglasso`.
2. **Neumann-Markov Resolvent Decomposition**: Partitions the reproduced covariance matrix via the resolvent operator $\mathbf{Q} = (\mathbf{I} - \hat{\mathbf{\Omega}})^{-1} \equiv \mathbf{I} + \hat{\mathbf{\Omega}} + \hat{\mathbf{\Omega}}(\mathbf{I} - \hat{\mathbf{\Omega}})^{-1}\hat{\mathbf{\Omega}}$.
3. **Nodal Variance Conservation Law**: Decomposes each item's standardized unit variance into three interpretable parameters summing to 100% of item variance:
   $$h_i^2 \text{ (Within-Community Signal)} + \beta_i^2 \text{ (Bridging Variance)} + \psi_i^2 \text{ (Markov Uniqueness)} \equiv 1 \text{ (100\%)}$$
4. **Purity-Weighted Community Detection ($A_{ij}^{\text{rEGA}}$)**: Weights active edges by topological purity ($h_{\mathcal{G},i}^2 = 1 - \psi_i^2$) prior to Walktrap, Leiden, or Louvain partitioning.
5. **Network Reliability Coefficient $\Gamma$ ($\Gamma_c, \Gamma_t$) & Graph $SEM_c$**: Provides closed-form dimensional and total reliability coefficients alongside standardized and raw-score Standard Errors of Measurement ($SEM_c^{std}, SEM_c^{raw}$) and Global Graph Noise ($\Psi_{\mathcal{G}}$).
6. **Compositional Evaluation Without Rigid Cutoffs**: Because $h_i^2 + \beta_i^2 + \psi_i^2 = 1$ forms an exact conservation system, `rEGA` and `boot.rEGA` evaluate nodal variance components, structural reliability, and item replicability jointly without imposing arbitrary fixed thresholds.

---

## Installation

Install the package directly from GitHub:

```r
# If remotes is not installed: install.packages("remotes")
remotes::install_github("Huerta-Huacchahuilca/rEGA")
```

Alternatively, install using `pak`:

```r
# install.packages("pak")
pak::pak("Huerta-Huacchahuilca/rEGA")
```

---

## Comprehensive User Guide

### 1. Model Estimation (`rEGA`)

`rEGA()` estimates the support-refitted network from raw item responses (or a correlation matrix with sample size `n`), supports flexible item selection via `items` (column names, indices, or ranges), sets `seed = 11` by default for full reproducibility, and displays a side-by-side figure combining the **EGA-style Network** and the **Nodal Variance Triad Decomposition (summing to 100% per item)**:

```r
library(rEGA)

# Simulated 2-dimension structure (6 items, N = 300)
set.seed(11)
Sigma <- matrix(0.28, 6, 6)
Sigma[1:3, 1:3] <- 0.68
Sigma[4:6, 4:6] <- 0.68
diag(Sigma) <- 1.0
data_items <- as.data.frame(MASS::mvrnorm(n = 300, mu = rep(0, 6), Sigma = Sigma))
colnames(data_items) <- paste0("Item_", 1:6)

# Fit the rEGA model (default seed = 11)
fit <- rEGA(
  data       = data_items,
  items      = 1:6,
  cor.method = "cor_auto",
  algorithm  = "walktrap",
  plot.rEGA  = TRUE,
  seed       = 11
)

# Multi-row console summary of reliability, nodal triad, and network loadings
print(fit)
summary(fit)
```

### 2. Visualizing the Network and Standardized Nodal Percentages (`plot.rEGA`)

`plot.rEGA()` supports **combined (`type = "both"`)** or **separate (`type = "separate"`)** panel displays, EGA-compatible palettes (**`palette = "polychrome"`** and **`palette = "colorblind"`**), network layouts (**`layout = "fr"`** for Fruchterman-Reingold and **`layout = "circle"`** for circular layout), and automatic percentage printing inside bar segments (`show.values = TRUE`, `label.min.pct = 10` so narrow segments never overflow):

```r
# 1. Combined side-by-side figure (EGA-style Network + 100% Nodal Variance Triad)
plot(fit, type = "both", palette = "polychrome", layout = "fr", show.values = TRUE, label.min.pct = 10)

# 2. Separate individual figures (exported sequentially to the RStudio Plots pane)
plot(fit, type = "separate", palette = "polychrome", layout = "fr")

# 3. Colorblind-safe palette (Okabe-Ito) and circular network layout
plot(fit, type = "both", palette = "colorblind", layout = "circle")

# 4. Single diagnostic panels
plot(fit, type = "network", palette = "polychrome", layout = "fr")
plot(fit, type = "variance", palette = "polychrome", show.values = TRUE, label.min.pct = 10)
```

### 3. Parallel Bootstrap Inference and Stability Analysis (`boot.rEGA`)

`boot.rEGA()` performs parametric or non-parametric case resampling with a real-time console progress bar and **Estimated Time of Arrival (ETA)**, using `seed = 11` by default, and generates either a unified side-by-side figure (`type = "both"`) or separate figures (`type = "separate"`):

```r
boot_fit <- boot.rEGA(
  rEGA.obj       = fit,
  iter           = 500,
  type           = "parametric",
  ncores         = 2,
  seed           = 11,
  progress       = TRUE,
  plot.boot.rEGA = TRUE
)

print(boot_fit)
summary(boot_fit)

# Visualization options for boot.rEGA:
plot(boot_fit, type = "both", palette = "polychrome")      # Unified rEGA Network + Item Stability
plot(boot_fit, type = "separate", palette = "colorblind")  # Separate Network and Stability figures
plot(boot_fit, type = "all", palette = "polychrome")       # 3-Panel: Network + Nodal Triad (%) + Stability
plot(boot_fit, type = "stability", palette = "colorblind") # Item Stability panel only
```

### 4. Attenuation-Aware Power Planning & Open Science Audit (`rEGA.power` & `rEGA.audit`)

```r
# Sample size planning accounting for network measurement error (Henry, 2026)
power_plan <- rEGA.power(rEGA.obj = fit, p.target = c(6, 9, 12), r.true = 0.30)
print(power_plan)

# Export 5-sheet Excel audit workbook with random seeds and SHA-256 cryptographic data hash
rEGA.audit(rEGA.obj = fit, boot.obj = boot_fit, file = "rEGA_Audit_Report.xlsx")
```

---

## Author & Maintainer

- **César Huerta-Huacchahuilca**  
  *Creator, Lead Developer & Maintainer*  
  *Email:* [andrehuerta1411@gmail.com](mailto:andrehuerta1411@gmail.com)  
  *GitHub Profile:* [@Huerta-Huacchahuilca](https://github.com/Huerta-Huacchahuilca)  
  *ORCID:* [0009-0000-8906-2678](https://orcid.org/0009-0000-8906-2678)  
  *Issue Tracker:* [rEGA Issues](https://github.com/Huerta-Huacchahuilca/rEGA/issues)

---

## Citation & Academic Attribution

If you use **`rEGA`** in academic studies, thesis projects, or methodological research, please cite the software package (solely developed by C. Huerta-Huacchahuilca) as:

### APA 7th Edition

> Huerta-Huacchahuilca, C. (2026). *rEGA: Reliability Exploratory Graph Analysis* (R package version 1.0.0) [Computer software]. GitHub. [https://github.com/Huerta-Huacchahuilca/rEGA](https://github.com/Huerta-Huacchahuilca/rEGA)

### BibTeX Entry

```bibtex
@manual{huerta2026rega,
  title  = {rEGA: Reliability Exploratory Graph Analysis},
  author = {Huerta-Huacchahuilca, C{\'e}sar},
  year   = {2026},
  note   = {R package version 1.0.0},
  url    = {https://github.com/Huerta-Huacchahuilca/rEGA}
}
```

### Accompanying Methodological Paper & OSF Repository

To cite the theoretical framework, Monte Carlo simulation, and empirical validation studies:

> Huerta-Huacchahuilca, C., & Ramos-Diaz, J. (2026). *Reliability Exploratory Graph Analysis (rEGA): Nodal Variance Decomposition and Coefficient Gamma for Psychometric Networks*. Open Science Framework (OSF). (https://doi.org/10.17605/OSF.IO/7QPRA)
