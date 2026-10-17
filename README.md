# rEGA: Reliability Exploratory Graph Analysis <img src="man/figures/logo.png" align="right" height="120" alt="" />

[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Open Science: OSF Ready](https://img.shields.io/badge/Open%20Science-SHA256%20Audited-orange.svg)](https://osf.io/)

**`rEGA`** (*Reliability Exploratory Graph Analysis*) is an R package for network psychometrics that bridges Exploratory Graph Analysis (`EGA`; Golino & Epskamp, 2017; Golino et al., 2020) with exact measurement error decomposition and score reliability estimation directly from Gaussian Graphical Models (GGMs) and Markov Random Fields.

---

## Why `rEGA`?

Standard EGA pipelines estimate precision matrices $\mathbf{K} = \mathbf{\Sigma}^{-1}$, extract partial correlations $\mathbf{\Omega}$, and discard the conditional precision diagonal $\mathbf{D}_K = \text{diag}(\mathbf{K})$, implicitly treating all nodes as equally free of conditional error. **`rEGA`** resolves this limitation through five sequential innovations:

1. **Support-Constrained MLE Refitting ($\hat{\mathbf{K}}_{\text{refit}}$)**: Eliminates $\ell_1$ (LASSO) shrinkage bias on active edges $E_\lambda$ identified by `EBICglasso`.
2. **Neumann-Markov Resolvent Decomposition**: Partitions the reproduced covariance matrix via the resolvent operator $\mathbf{Q} = (\mathbf{I} - \hat{\mathbf{\Omega}})^{-1} \equiv \mathbf{I} + \hat{\mathbf{\Omega}} + \hat{\mathbf{\Omega}}(\mathbf{I} - \hat{\mathbf{\Omega}})^{-1}\hat{\mathbf{\Omega}}$.
3. **Nodal Variance Conservation Law**: Decomposes each item's unit variance into three interpretable parameters:
   $$h_i^2 \text{ (Within-Community Signal)} + \beta_i^2 \text{ (Bridging Variance)} + \psi_i^2 \text{ (Markov Uniqueness)} \equiv 1$$
4. **Purity-Weighted Community Detection ($A_{ij}^{\text{rEGA}}$)**: Weights edges by topological purity ($h_{\mathcal{G},i}^2 = 1 - \psi_i^2$) prior to Walktrap, Leiden, or Louvain partitioning.
5. **Network Reliability Coefficient $\Gamma$ ($\Gamma_c, \Gamma_t$) & Graph $\text{SEM}_c$**: Provides closed-form dimensional and total reliability coefficients alongside standardized and raw-score Standard Errors of Measurement.

---

## Installation

You can install the development version of **`rEGA`** directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("yourusername/rEGA")