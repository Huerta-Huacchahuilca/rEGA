#' @title Reliability Exploratory Graph Analysis (rEGA)
#'
#' @description
#' Estimates Reliability Exploratory Graph Analysis (\code{rEGA}) from raw item responses
#' or a correlation matrix by decomposing the support-refitted Markov precision matrix
#' via the Neumann-Markov resolvent operator. Partitions each item's standardized unit
#' variance into Within-Community Network Communality (\code{h2}), Inter-Community
#' Bridging Variance (\code{beta2}), and Markov Nodal Uniqueness (\code{psi2}),
#' satisfying the exact conservation law \code{h2 + beta2 + psi2 = 1} (100 percent of nodal variance).
#' Computes closed-form structural reliability coefficients (Coefficient Gamma:
#' \code{Gamma_c} per dimension and \code{Gamma_t} total), Global Graph Noise
#' (\code{Psi_G}), and Graph Standard Error of Measurement (\code{SEM_c}).
#'
#' @param data A numeric \code{data.frame} or \code{matrix} containing raw item responses
#'   (n x p) or a symmetric correlation matrix (p x p).
#' @param items Optional character vector of column names or numeric vector of column indices
#'   specifying the subset of items to analyze from \code{data}. Default is \code{NULL} (uses all columns).
#' @param n Integer. Sample size required only when \code{data} is a correlation matrix.
#'   Automatically detected when \code{data} contains raw item responses. Default is \code{NULL}.
#' @param cor.method Character. Correlation estimation method:
#'   \code{"cor_auto"} (default), \code{"polychoric"}, \code{"pearson"}, or \code{"spearman"}.
#' @param algorithm Character. Community detection algorithm applied to the purity-weighted
#'   adjacency matrix: \code{"walktrap"} (default), \code{"leiden"}, or \code{"louvain"}.
#' @param plot.rEGA Logical. Should the diagnostic \code{rEGA} visualization (EGA-style network
#'   and Nodal Variance Triad in percentages) be displayed upon estimation? Default is \code{TRUE}.
#' @param seed Integer. Master random seed for reproducibility in unidimensionality checks,
#'   network layout, and community partitioning. Default is \code{2026L}.
#' @param ... Additional arguments for plot customization (\code{plot.type = "both"},
#'   \code{show.values = TRUE}) or advanced algorithmic control (see section below).
#'
#' @section Advanced Control Arguments:
#' Advanced algorithmic hyperparameters (\code{gamma}, \code{uni.method}, \code{lambda.min.ratio},
#' \code{nlambda}, \code{walktrap.steps}, \code{leiden.objective}, \code{resolution.parameter},
#' \code{weights.type}, and \code{refit.tol}) can be passed via \code{...}.
#' These parameters control the core Neumann-Markov resolving algorithms and EBICglasso
#' penalization. Modifying them deviates from the optimal standardized architecture and
#' should be theoretically justified in the manuscript to prevent methodological flexibility.
#'
#' @section Holistic Interpretation (No Rigid Cutoffs):
#' \code{rEGA} intentionally avoids rigid or arbitrary cutoff thresholds. Because the Nodal
#' Variance Triad forms an exact compositional conservation system (\code{h2 + beta2 + psi2 = 1},
#' or 100 percent of item variance), parameters must be interpreted holistically in conjunction with
#' structural reliability (\code{Gamma_c}, \code{Gamma_t}), Global Graph Noise (\code{Psi_G}),
#' and bootstrap structural stability in \code{\link{boot.rEGA}}:
#' \itemize{
#'   \item \strong{Essential Unidimensionality vs. Lexical Doublets}: When partitioning splits a scale
#'         into minor clusters where items display bridging variance exceeding within-community signal
#'         (\code{beta2 > h2}) alongside high total reliability (\code{Gamma_t}) and lower sub-dimension
#'         reliability (\code{Gamma_c}), the split reflects local semantic redundancy rather than substantive
#'         multidimensionality.
#'   \item \strong{Transdiagnostic Bridging Nodes}: In multidimensional instruments, an item with prominent
#'         bridging variance (\code{beta2}) and moderate within-community signal (\code{h2}) but controlled
#'         Markov uniqueness (\code{psi2}) operates as a substantive bridge connecting communities.
#'   \item \strong{Noisy or Defective Items}: An item is psychometrically problematic when convergent
#'         evidence shows dominant Markov uniqueness (\code{psi2} exceeding total network signal
#'         \code{h2 + beta2}), elevated Global Graph Noise (\code{Psi_G}), and unstable dimension
#'         replication in \code{\link{boot.rEGA}}.
#' }
#'
#' @return An object of class \code{"rEGA"} containing:
#' \itemize{
#'   \item \code{network.refit}: Unbiased partial correlation matrix (\code{Omega_refit}) after support refitting.
#'   \item \code{network.glasso}: Initial regularized network (\code{Omega_glasso}).
#'   \item \code{adjacency.rEGA}: Purity-weighted adjacency matrix (\code{A_rEGA}).
#'   \item \code{K.refit}: Unbiased support-refitted precision matrix (\code{K_refit}).
#'   \item \code{D.K}: Diagonal matrix of conditional nodal precisions (\code{D_K}).
#'   \item \code{resolvent}: Global Neumann-Markov resolvent matrix (\code{Q = (I - Omega_refit)^(-1)}).
#'   \item \code{R.obs}: Observed positive-definite correlation matrix.
#'   \item \code{wc}: Named integer vector of community memberships.
#'   \item \code{n.dim}: Number of estimated dimensions.
#'   \item \code{nodal.triad}: Data frame of item-level parameters (\code{h2}, \code{beta2},
#'         \code{psi2}, percentages \code{h2_pct}, \code{beta2_pct}, \code{psi2_pct},
#'         \code{h2_G}, \code{K_ii}, and \code{Conservation_Sum}).
#'   \item \code{net.loadings}: Matrix of standardized network loadings across communities.
#'   \item \code{reliability}: List containing dimension-level summary table (\code{$dimensions}:
#'         \code{Gamma_c}, \code{SEM_std}, \code{SEM_raw}) and global metrics
#'         (\code{$Gamma.t}, \code{$Psi.G}, \code{$SEM.total.std}, \code{$SEM.total.raw}).
#'   \item \code{audit.metadata}: List of algebraic verification checks, active edge counts,
#'         \code{gamma.used}, and SHA-256 cryptographic data hash.
#' }
#'
#' @importFrom stats cor cov2cor sd median quantile rnorm complete.cases setNames
#' @importFrom utils globalVariables packageVersion sessionInfo
#' @importFrom Matrix nearPD
#' @importFrom MASS mvrnorm
#' @importFrom glasso glasso
#' @importFrom qgraph EBICglasso cor_auto
#' @importFrom igraph graph_from_adjacency_matrix cluster_walktrap cluster_leiden cluster_louvain membership E layout_with_fr
#' @importFrom psych polychoric cor.smooth
#' @export
#'
#' @examples
#' set.seed(2026)
#' Sigma <- matrix(0.25, 6, 6)
#' Sigma[1:3, 1:3] <- 0.65
#' Sigma[4:6, 4:6] <- 0.65
#' diag(Sigma) <- 1
#' dat <- as.data.frame(MASS::mvrnorm(n = 300, mu = rep(0, 6), Sigma = Sigma))
#' colnames(dat) <- paste0("Item_", 1:6)
#' fit <- rEGA(data = dat, items = 1:6, cor.method = "pearson", plot.rEGA = FALSE)
#' print(fit)
rEGA <- function(
  data,
  items = NULL,
  n = NULL,
  cor.method = c("cor_auto", "polychoric", "pearson", "spearman"),
  algorithm = c("walktrap", "leiden", "louvain"),
  plot.rEGA = TRUE,
  seed = 2026L,
  ... # <-- Elipsis para parametros avanzados ocultos
) {

  cor.method <- match.arg(cor.method)
  algorithm <- match.arg(algorithm)

  # ---------------------------------------------------------------------------
  # EXTRACCION DE ARGUMENTOS AVANZADOS (Control Estricto / Anti p-hacking)
  # ---------------------------------------------------------------------------
  args <- list(...)

  gamma <- if (!is.null(args$gamma)) args$gamma else 0.50
  uni.method <- if (!is.null(args$uni.method)) args$uni.method else "expand"
  lambda.min.ratio <- if (!is.null(args$lambda.min.ratio)) args$lambda.min.ratio else 0.10
  nlambda <- if (!is.null(args$nlambda)) as.integer(args$nlambda) else 100L
  walktrap.steps <- if (!is.null(args$walktrap.steps)) as.integer(args$walktrap.steps) else 4L
  leiden.objective <- if (!is.null(args$leiden.objective)) args$leiden.objective else "modularity"
  resolution.parameter <- if (!is.null(args$resolution.parameter)) args$resolution.parameter else 1.0
  weights.type <- if (!is.null(args$weights.type)) args$weights.type else "unit"
  refit.tol <- if (!is.null(args$refit.tol)) args$refit.tol else 1e-6

  plot.type <- if (!is.null(args$plot.type)) args$plot.type else "both"
  show.values <- if (!is.null(args$show.values)) args$show.values else TRUE

  t0 <- Sys.time()
  raw_data <- NULL
  item_sds <- NULL
  R_obs <- NULL
  Omega_glasso <- NULL
  gamma_used <- gamma

  # ---------------------------------------------------------------------------
  # STAGE 1: DATA EXTRACTION & ITEM SELECTION (RAW DATA STRICT)
  # ---------------------------------------------------------------------------
  if (missing(data) || is.null(data)) {
    stop("Error: Debe proporcionar su base de datos (data.frame o matrix) en 'data'.")
  }

  # 1.1 Filtrado agil de items (Soporta rangos, nombres o paste0)
  if (!is.null(items)) {
    if (is.data.frame(data) || is.matrix(data)) {
      data_input <- data[, items, drop = FALSE]
    } else {
      stop("Error: El argumento 'items' solo funciona si ingresa un data.frame o matrix en 'data'.")
    }
  } else {
    data_input <- data
  }

  # 1.2 Deteccion: Es Matriz de Correlaciones o Datos Crudos?
  mat_input <- as.matrix(data_input)

  # Condicion: Matriz cuadrada, simetrica y con 1s en la diagonal
  if (nrow(mat_input) == ncol(mat_input) &&
      isSymmetric(mat_input, tol = 1e-5) &&
      all(abs(diag(mat_input) - 1) < 1e-4)) {
    R_obs <- mat_input
    if (is.null(n)) {
      stop("Error: Al ingresar una matriz de correlaciones, debe especificar el tamano de muestra 'n'.")
    }
    item_sds <- rep(1.0, ncol(R_obs))
    pathway_used <- "Correlation_Matrix"

  } else {
    # Proceso de Datos Crudos (Auto-deteccion de 'n')
    raw_data <- as.data.frame(data_input)
    raw_data <- raw_data[stats::complete.cases(raw_data), , drop = FALSE]
    n <- nrow(raw_data)

    if (n < 10) {
      stop("Error: El tamano de muestra 'n' es demasiado pequeno despues de remover valores perdidos (NA).")
    }

    item_sds <- apply(raw_data, 2, stats::sd, na.rm = TRUE)
    R_obs <- .compute_rEGA_cor(raw_data, cor.method = cor.method)
    pathway_used <- "Raw_Data"
  }

  p <- ncol(R_obs)
  item_names <- colnames(R_obs)
  if (is.null(item_names)) item_names <- paste0("V", seq_len(p))
  dimnames(R_obs) <- list(item_names, item_names)
  names(item_sds) <- item_names

  # Asegurar matriz definida positiva para continuar con el operador resolvente
  R_obs <- .ensure_pd_cor(R_obs)

  # ---------------------------------------------------------------------------
  # STAGE 2: UNIDIMENSIONALITY CHECK & GLASSO
  # ---------------------------------------------------------------------------
  is_unidimensional_override <- FALSE
  if (uni.method == "expand") {
    is_unidimensional_override <- .check_unidimensional_golino(
      raw_data = raw_data, R_obs = R_obs, n = n,
      cor.method = cor.method, gamma = gamma,
      lambda.min.ratio = lambda.min.ratio, nlambda = nlambda,
      algorithm = algorithm, walktrap.steps = walktrap.steps,
      leiden.objective = leiden.objective,
      resolution.parameter = resolution.parameter,
      seed = seed
    )
  }

  ebic_res <- .estimate_ebic_glasso_iterative(
    R_obs = R_obs, n = n, gamma = gamma,
    lambda.min.ratio = lambda.min.ratio, nlambda = nlambda
  )
  Omega_glasso <- ebic_res$Omega_glasso
  gamma_used <- ebic_res$gamma_used
  dimnames(Omega_glasso) <- list(item_names, item_names)

  # ---------------------------------------------------------------------------
  # STAGE 3: SUPPORT EXTRACTION & CONSTRAINED MLE REFITTING (K_refit)
  # ---------------------------------------------------------------------------
  refit_obj <- .refit_precision_support(R_obs = R_obs, Omega_glasso = Omega_glasso, tol = refit.tol)
  K_refit <- refit_obj$K_refit
  Omega_refit <- refit_obj$Omega_refit
  D_K <- refit_obj$D_K
  active_edges <- refit_obj$active_edges

  # ---------------------------------------------------------------------------
  # STAGE 4: NEUMANN-MARKOV RESOLVENT OPERATOR & GLOBAL NODAL PURITY
  # ---------------------------------------------------------------------------
  I_p <- diag(1, p)
  dimnames(I_p) <- list(item_names, item_names)
  Q_res <- tryCatch(
    solve(I_p - Omega_refit),
    error = function(e) as.matrix(solve(Matrix::nearPD(I_p - Omega_refit, corr = FALSE)$mat))
  )
  dimnames(Q_res) <- list(item_names, item_names)
  q_diag <- pmax(diag(Q_res), 1.0)
  psi2_i <- 1.0 / q_diag
  h2_G_i <- pmax(0.0, 1.0 - psi2_i)

  # ---------------------------------------------------------------------------
  # STAGE 5: PURITY-WEIGHTED ADJACENCY MATRIX & COMMUNITY DETECTION
  # ---------------------------------------------------------------------------
  purity_outer <- sqrt(outer(h2_G_i, h2_G_i, "*"))
  A_rEGA <- abs(Omega_refit) * purity_outer
  diag(A_rEGA) <- 0.0
  A_rEGA <- (A_rEGA + t(A_rEGA)) / 2.0
  dimnames(A_rEGA) <- list(item_names, item_names)

  if (isTRUE(is_unidimensional_override)) {
    wc <- rep(1L, p)
  } else {
    wc <- .detect_communities_rEGA(
      A_mat = A_rEGA, algorithm = algorithm, walktrap.steps = walktrap.steps,
      leiden.objective = leiden.objective, resolution.parameter = resolution.parameter, seed = seed
    )
  }
  names(wc) <- item_names
  n_dim <- length(unique(wc[!is.na(wc)]))

  # ---------------------------------------------------------------------------
  # STAGE 6: NODAL VARIANCE CONSERVATION TRIAD (h_i^2 + beta_i^2 + psi_i^2 = 1)
  # ---------------------------------------------------------------------------
  h2_i <- numeric(p)
  beta2_i <- numeric(p)
  names(h2_i) <- item_names
  names(beta2_i) <- item_names
  unique_comms <- sort(unique(wc))
  for (c_id in unique_comms) {
    idx_c <- which(wc == c_id)
    if (length(idx_c) == 1L) {
      q_intra_c <- 1.0
    } else {
      I_c <- diag(1, length(idx_c))
      Omega_cc <- Omega_refit[idx_c, idx_c, drop = FALSE]
      Q_cc <- tryCatch(solve(I_c - Omega_cc), error = function(e) diag(1, length(idx_c)))
      q_intra_c <- pmax(diag(Q_cc), 1.0)
    }
    inv_q_intra <- 1.0 / q_intra_c
    h2_i[idx_c] <- pmax(0.0, 1.0 - inv_q_intra)
    beta2_i[idx_c] <- pmax(0.0, inv_q_intra - psi2_i[idx_c])
  }

  triad_sum <- h2_i + beta2_i + psi2_i
  # Standardized percentages summing exactly to 100.0% per item
  h2_pct <- round((h2_i / triad_sum) * 100, 1)
  beta2_pct <- round((beta2_i / triad_sum) * 100, 1)
  psi2_pct <- round(100.0 - h2_pct - beta2_pct, 1)

  k_diag <- diag(K_refit)
  nodal_triad <- data.frame(
    Item = item_names,
    Community = wc,
    h2 = round(h2_i, 6),
    beta2 = round(beta2_i, 6),
    psi2 = round(psi2_i, 6),
    h2_pct = h2_pct,
    beta2_pct = beta2_pct,
    psi2_pct = psi2_pct,
    h2_G = round(h2_G_i, 6),
    K_ii = round(k_diag, 6),
    Conservation_Sum = round(triad_sum, 6),
    stringsAsFactors = FALSE
  )
  rownames(nodal_triad) <- NULL

  # ---------------------------------------------------------------------------
  # STAGE 7: STANDARDIZED NETWORK LOADINGS
  # ---------------------------------------------------------------------------
  net_loadings <- .compute_network_loadings(Omega_refit, wc, item_names)

  # ---------------------------------------------------------------------------
  # STAGE 8: STRUCTURAL RELIABILITY (GAMMA_c, GAMMA_t, PSI_G, AND GRAPH SEM)
  # ---------------------------------------------------------------------------
  Sigma_reprod <- tryCatch(solve(K_refit), error = function(e) R_obs)
  dimnames(Sigma_reprod) <- list(item_names, item_names)
  D_K_inv <- diag(1.0 / k_diag, p, p)
  D_K_inv_sqrt <- diag(1.0 / sqrt(k_diag), p, p)
  dimnames(D_K_inv) <- list(item_names, item_names)
  dimnames(D_K_inv_sqrt) <- list(item_names, item_names)

  Psi_G <- sum(diag(D_K_inv)) / sum(diag(Sigma_reprod))
  w_vec <- if (weights.type == "purity") sqrt(pmax(h2_G_i, 1e-6)) else rep(1.0, p)
  names(w_vec) <- item_names

  num_total <- as.numeric(t(w_vec) %*% D_K_inv %*% w_vec)
  den_total <- as.numeric(t(w_vec) %*% Sigma_reprod %*% w_vec)
  Gamma_t <- pmax(0.0, pmin(1.0, 1.0 - (num_total / den_total)))
  SEM_total_std <- sqrt(num_total)

  if (!is.null(raw_data)) {
    total_raw_sd <- stats::sd(as.matrix(raw_data) %*% w_vec, na.rm = TRUE)
  } else {
    D_sd <- diag(item_sds, p, p)
    cov_raw <- D_sd %*% R_obs %*% D_sd
    total_raw_sd <- sqrt(as.numeric(t(w_vec) %*% cov_raw %*% w_vec))
  }
  SEM_total_raw <- total_raw_sd * sqrt(1.0 - Gamma_t)

  dim_rows <- lapply(unique_comms, function(c_id) {
    idx_c <- which(wc == c_id)
    p_c <- length(idx_c)
    w_c <- w_vec[idx_c]
    k_c <- k_diag[idx_c]
    if (p_c == 1L) {
      Gamma_c <- h2_G_i[idx_c]
      SEM_c_std <- sqrt(1.0 / k_c)
    } else {
      I_c <- diag(1, p_c)
      Omega_cc <- Omega_refit[idx_c, idx_c, drop = FALSE]
      Q_cc <- tryCatch(solve(I_c - Omega_cc), error = function(e) I_c)
      D_Kc_inv <- diag(1.0 / k_c, p_c, p_c)
      D_Kc_inv_sqrt <- diag(1.0 / sqrt(k_c), p_c, p_c)
      Sigma_cc_intra <- D_Kc_inv_sqrt %*% Q_cc %*% D_Kc_inv_sqrt
      err_var_c <- as.numeric(t(w_c) %*% D_Kc_inv %*% w_c)
      tot_var_c <- as.numeric(t(w_c) %*% Sigma_cc_intra %*% w_c)
      Gamma_c <- pmax(0.0, pmin(1.0, 1.0 - (err_var_c / tot_var_c)))
      SEM_c_std <- sqrt(err_var_c)
    }
    if (!is.null(raw_data)) {
      sd_raw_c <- stats::sd(as.matrix(raw_data[, idx_c, drop = FALSE]) %*% w_c, na.rm = TRUE)
    } else {
      D_sd_c <- diag(item_sds[idx_c], p_c, p_c)
      cov_raw_c <- D_sd_c %*% R_obs[idx_c, idx_c, drop = FALSE] %*% D_sd_c
      sd_raw_c <- sqrt(as.numeric(t(w_c) %*% cov_raw_c %*% w_c))
    }
    SEM_c_raw <- sd_raw_c * sqrt(1.0 - Gamma_c)
    data.frame(
      Dimension = c_id, n_Items = p_c, Gamma_c = round(Gamma_c, 6),
      Mean_h2 = round(mean(h2_i[idx_c]), 6), Mean_beta2 = round(mean(beta2_i[idx_c]), 6),
      Mean_psi2 = round(mean(psi2_i[idx_c]), 6), SEM_std = round(SEM_c_std, 6),
      SD_Raw = round(sd_raw_c, 6), SEM_raw = round(SEM_c_raw, 6), stringsAsFactors = FALSE
    )
  })
  dim_reliability_df <- do.call(rbind, dim_rows)

  Level0_Error <- D_K_inv
  Level1_Direct <- D_K_inv_sqrt %*% Omega_refit %*% D_K_inv_sqrt
  Level2_Reverb <- D_K_inv_sqrt %*% (Omega_refit %*% Q_res %*% Omega_refit) %*% D_K_inv_sqrt
  dimnames(Level1_Direct) <- list(item_names, item_names)
  dimnames(Level2_Reverb) <- list(item_names, item_names)
  alg_max_err <- max(abs(Sigma_reprod - (Level0_Error + Level1_Direct + Level2_Reverb)))

  out <- list(
    network.refit = Omega_refit, network.glasso = Omega_glasso, adjacency.rEGA = A_rEGA,
    K.refit = K_refit, D.K = D_K, resolvent = Q_res, R.obs = R_obs, Sigma.reprod = Sigma_reprod,
    resolvent.levels = list(
      Level0_Conditional_Error = Level0_Error,
      Level1_Direct_Coupling = Level1_Direct,
      Level2_Systemic_Reverberation = Level2_Reverb
    ),
    wc = wc, n.dim = n_dim, nodal.triad = nodal_triad, net.loadings = net_loadings,
    reliability = list(
      dimensions = dim_reliability_df, Gamma.t = Gamma_t, Psi.G = Psi_G,
      SEM.total.std = SEM_total_std, SEM.total.raw = SEM_total_raw
    ),
    audit.metadata = list(
      pathway = pathway_used, n = n, p = p, cor.method = cor.method, algorithm = algorithm,
      gamma.used = gamma_used, active.edges = active_edges, seed = seed,
      min.eigen.Krefit = min(eigen(K_refit, symmetric = TRUE, only.values = TRUE)$values),
      max.algebraic.error = alg_max_err, max.conservation.error = max(abs(triad_sum - 1.0)),
      data.sha256 = digest::digest(if (!is.null(raw_data)) raw_data else R_obs, algo = "sha256"),
      runtime.sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ),
    raw.data = raw_data, call = match.call()
  )
  class(out) <- "rEGA"
  if (isTRUE(plot.rEGA)) {
    p_res <- plot(out, type = plot.type, show.values = show.values)
    if (inherits(p_res, "ggplot")) print(p_res)
  }
  return(out)
}

.compute_rEGA_cor <- function(df, cor.method = "cor_auto") {
  df_num <- as.data.frame(lapply(df, as.numeric))
  colnames(df_num) <- colnames(df)

  if (cor.method == "cor_auto") {
    is_ord <- vapply(df_num, function(x) {
      length(unique(x[!is.na(x)])) <= 7L
    }, logical(1))

    if (!any(is_ord)) {
      R <- stats::cor(df_num, use = "pairwise.complete.obs", method = "pearson")
    } else {
      R <- tryCatch(
        suppressMessages(suppressWarnings(qgraph::cor_auto(df_num, forcePD = TRUE))),
        error = function(e) {
          tryCatch(
            suppressMessages(suppressWarnings(psych::polychoric(df_num, smooth = TRUE, correct = 0.01)$rho)),
            error = function(e2) stats::cor(df_num, use = "pairwise.complete.obs", method = "pearson")
          )
        }
      )
    }
  } else if (cor.method == "polychoric") {
    R <- tryCatch(
      suppressMessages(suppressWarnings(psych::polychoric(df_num, smooth = TRUE, correct = 0.01)$rho)),
      error = function(e) stats::cor(df_num, use = "pairwise.complete.obs", method = "pearson")
    )
  } else if (cor.method == "spearman") {
    R <- stats::cor(df_num, use = "pairwise.complete.obs", method = "spearman")
  } else {
    R <- stats::cor(df_num, use = "pairwise.complete.obs", method = "pearson")
  }

  R <- as.matrix(R)
  colnames(R) <- rownames(R) <- colnames(df)
  return(R)
}

.ensure_pd_cor <- function(R) {
  R <- (R + t(R)) / 2.0
  diag(R) <- 1.0
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < 1e-5) R <- as.matrix(Matrix::nearPD(R, corr = TRUE, keepDiag = TRUE)$mat)
  return(R)
}

.estimate_ebic_glasso_iterative <- function(R_obs, n, gamma = 0.50, lambda.min.ratio = 0.10, nlambda = 100L) {
  candidate_gammas <- sort(unique(c(gamma, 0.25, 0.00)), decreasing = TRUE)
  candidate_gammas <- candidate_gammas[candidate_gammas <= gamma]
  Omega_out <- NULL
  g_selected <- gamma
  for (g_val in candidate_gammas) {
    net_try <- suppressMessages(suppressWarnings(
      qgraph::EBICglasso(S = R_obs, n = n, gamma = g_val, lambda.min.ratio = lambda.min.ratio,
                         nlambda = nlambda, returnAllResults = FALSE, checkPD = FALSE)
    ))
    net_try <- as.matrix(net_try)
    deg <- rowSums(abs(net_try) > 1e-7)
    Omega_out <- net_try
    g_selected <- g_val
    if (all(deg > 0)) break
  }
  return(list(Omega_glasso = Omega_out, gamma_used = g_selected))
}

.refit_precision_support <- function(R_obs, Omega_glasso, tol = 1e-6) {
  p <- ncol(R_obs)
  item_names <- colnames(R_obs)
  E_active <- (abs(Omega_glasso) > 1e-6)
  diag(E_active) <- FALSE
  active_edges <- sum(E_active[upper.tri(E_active)])
  if (active_edges == 0L) {
    K_refit <- diag(1, p)
  } else {
    zero_mat <- which(!E_active & row(E_active) != col(E_active), arr.ind = TRUE)
    if (nrow(zero_mat) == 0L) {
      K_refit <- tryCatch(solve(R_obs), error = function(e) glasso::glasso(s = R_obs, rho = 1e-4, penalize.diagonal = FALSE)$wi)
    } else {
      gl_refit <- suppressWarnings(glasso::glasso(s = R_obs, rho = 1e-8, zero = zero_mat, penalize.diagonal = FALSE, thr = tol))
      K_refit <- as.matrix(gl_refit$wi)
    }
  }
  K_refit <- (K_refit + t(K_refit)) / 2.0
  ev_k <- eigen(K_refit, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev_k) <= 1e-6) K_refit <- as.matrix(Matrix::nearPD(K_refit, corr = FALSE)$mat)
  dimnames(K_refit) <- list(item_names, item_names)
  k_diag <- pmax(diag(K_refit), 1e-6)
  D_K <- diag(k_diag, p, p)
  dimnames(D_K) <- list(item_names, item_names)
  D_K_inv_sqrt <- diag(1.0 / sqrt(k_diag), p, p)
  Omega_refit <- -D_K_inv_sqrt %*% K_refit %*% D_K_inv_sqrt
  diag(Omega_refit) <- 0.0
  Omega_refit[!E_active] <- 0.0
  Omega_refit <- (Omega_refit + t(Omega_refit)) / 2.0
  dimnames(Omega_refit) <- list(item_names, item_names)
  return(list(K_refit = K_refit, Omega_refit = Omega_refit, D_K = D_K, active_edges = active_edges))
}

.detect_communities_rEGA <- function(A_mat, algorithm = "walktrap", walktrap.steps = 4L,
                                     leiden.objective = "modularity", resolution.parameter = 1.0, seed = 2026L) {
  p <- ncol(A_mat)
  if (sum(A_mat) <= 1e-8) return(seq_len(p))
  g <- igraph::graph_from_adjacency_matrix(A_mat, mode = "max", weighted = TRUE, diag = FALSE)
  if (!is.null(seed)) set.seed(seed)
  if (algorithm == "leiden") {
    leiden_args <- list(graph = g, objective_function = leiden.objective, weights = igraph::E(g)$weight)
    if ("resolution" %in% names(formals(igraph::cluster_leiden))) {
      leiden_args$resolution <- resolution.parameter
    } else {
      leiden_args$resolution_parameter <- resolution.parameter
    }
    comm_obj <- do.call(igraph::cluster_leiden, leiden_args)
  } else if (algorithm == "louvain") {
    comm_obj <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  } else {
    comm_obj <- igraph::cluster_walktrap(g, weights = igraph::E(g)$weight, steps = walktrap.steps)
  }
  return(as.integer(igraph::membership(comm_obj)))
}

.check_unidimensional_golino <- function(raw_data, R_obs, n, cor.method, gamma, lambda.min.ratio,
                                         nlambda, algorithm, walktrap.steps, leiden.objective,
                                         resolution.parameter, seed) {
  set.seed(seed)
  p <- ncol(R_obs)
  p_sim <- 4L
  Sigma_sim <- matrix(0.49, p_sim, p_sim)
  diag(Sigma_sim) <- 1.0
  if (!is.null(raw_data)) {
    sim_dat <- MASS::mvrnorm(n = nrow(raw_data), mu = rep(0, p_sim), Sigma = Sigma_sim)
    colnames(sim_dat) <- paste0("SIM_UNI_", seq_len(p_sim))
    R_sim <- stats::cor(sim_dat)
    R_cross <- stats::cor(sim_dat, as.data.frame(lapply(raw_data, as.numeric)), use = "pairwise.complete.obs")
    R_exp <- rbind(cbind(R_sim, R_cross), cbind(t(R_cross), R_obs))
    colnames(R_exp) <- rownames(R_exp) <- c(colnames(sim_dat), colnames(R_obs))
  } else {
    R_exp <- matrix(0, p + p_sim, p + p_sim)
    R_exp[seq_len(p_sim), seq_len(p_sim)] <- Sigma_sim
    R_exp[(p_sim + 1):(p + p_sim), (p_sim + 1):(p + p_sim)] <- R_obs
    colnames(R_exp) <- rownames(R_exp) <- c(paste0("SIM_UNI_", seq_len(p_sim)), colnames(R_obs))
  }
  R_exp <- .ensure_pd_cor(R_exp)
  ebic_exp <- .estimate_ebic_glasso_iterative(R_exp, n = n, gamma = gamma, lambda.min.ratio = lambda.min.ratio, nlambda = nlambda)
  wc_exp <- .detect_communities_rEGA(abs(ebic_exp$Omega_glasso), algorithm, walktrap.steps, leiden.objective, resolution.parameter, seed)
  return(length(unique(wc_exp)) <= 2L)
}

.compute_network_loadings <- function(Omega_refit, wc, item_names) {
  unique_comms <- sort(unique(wc))
  K_dim <- length(unique_comms)
  p <- ncol(Omega_refit)
  A_abs <- abs(Omega_refit)
  diag(A_abs) <- 0.0
  load_mat <- matrix(0.0, nrow = p, ncol = K_dim)
  rownames(load_mat) <- item_names
  colnames(load_mat) <- paste0("Dim_", unique_comms)
  for (idx_k in seq_along(unique_comms)) {
    nodes_c <- which(wc == unique_comms[idx_k])
    load_mat[, idx_k] <- rowSums(A_abs[, nodes_c, drop = FALSE])
  }
  col_totals <- sqrt(colSums(load_mat))
  for (idx_k in seq_len(K_dim)) {
    if (col_totals[idx_k] > 1e-8) load_mat[, idx_k] <- load_mat[, idx_k] / col_totals[idx_k]
  }
  return(round(load_mat, 4))
}
