#' @title Parallel Bootstrap Inference and Structural Stability for rEGA
#'
#' @description
#' Performs parametric or non-parametric bootstrap resampling for \code{rEGA} models
#' with multi-core parallelization, real-time progress bar with Estimated Time of Arrival (ETA),
#' and reproducible random number streams. Computes structural consistency, item replicability,
#' 95 percent percentile confidence intervals for nodal variance components and Network Gamma
#' (\code{Gamma_c}, \code{Gamma_t}), and diagnostic agreement matrices.
#'
#' @param rEGA.obj An object of class \code{"rEGA"} fitted by \code{\link{rEGA}}.
#' @param iter Integer. Number of bootstrap replications. Default is \code{500L}.
#' @param type Character. Resampling scheme: \code{"parametric"} (default, simulates multivariate
#'   normal data from \code{R.obs}) or \code{"resampling"} (non-parametric case resampling with
#'   replacement from \code{raw.data}).
#' @param ncores Integer. Number of parallel CPU cores to use. Default is \code{1L}.
#' @param seed Integer. Master random seed for reproducible bootstrap streams. Default is \code{11L}.
#' @param progress Logical. Display real-time progress bar and estimated remaining time (ETA)? Default is \code{TRUE}.
#' @param plot.boot.rEGA Logical. Automatically display the \code{boot.rEGA} visualization upon completion?
#'   Default is \code{TRUE}.
#' @param ... Additional plotting arguments passed to \code{plot.boot.rEGA} (e.g.,
#'   \code{plot.type = "both"} or \code{"separate"}, \code{palette = "polychrome"} or \code{"colorblind"},
#'   \code{layout = "fr"} or \code{"circle"}).
#'
#' @return An object of class \code{"boot.rEGA"} containing:
#' \itemize{
#'   \item \code{rEGA.obj}: The original fitted \code{"rEGA"} model object.
#'   \item \code{seed}: The master random seed used for the bootstrap replications.
#'   \item \code{dimension.frequency}: Frequency table of estimated dimensions across replications.
#'   \item \code{structural.consistency}: Proportion of replications where each empirical community
#'         is recovered intact.
#'   \item \code{item.stability}: Data frame of item-level stability and replication proportions
#'         across all communities.
#'   \item \code{agreement.matrix}: Diagnostic dimension-by-dimension partition agreement matrix.
#'   \item \code{global.ci}: Summary table with 95 percent confidence intervals for
#'         \code{Gamma_t}, \code{Psi_G}, and \code{SEM_std}.
#'   \item \code{nodal.ci}: Item-level 95 percent confidence intervals for \code{h2},
#'         \code{beta2}, and \code{psi2}.
#'   \item \code{boot.log}: Complete replication-by-replication audit trail.
#' }
#'
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom progressr handlers progressor with_progress
#' @export
#'
#' @examples
#' set.seed(11)
#' Sigma <- matrix(0.25, 6, 6)
#' Sigma[1:3, 1:3] <- 0.65
#' Sigma[4:6, 4:6] <- 0.65
#' diag(Sigma) <- 1
#' dat <- as.data.frame(MASS::mvrnorm(n = 200, mu = rep(0, 6), Sigma = Sigma))
#' colnames(dat) <- paste0("Item", 1:6)
#' fit <- rEGA(data = dat, cor.method = "pearson", plot.rEGA = FALSE)
#' boot_fit <- boot.rEGA(fit, iter = 10, type = "parametric", ncores = 1,
#'                       progress = FALSE, plot.boot.rEGA = FALSE)
#' print(boot_fit)
boot.rEGA <- function(rEGA.obj,
                      iter = 500L,
                      type = c("parametric", "resampling"),
                      ncores = 1L,
                      seed = 11L,
                      progress = TRUE,
                      plot.boot.rEGA = TRUE,
                      ...) {

  if (!inherits(rEGA.obj, "rEGA")) {
    stop("'rEGA.obj' must be an object of class 'rEGA'.")
  }
  type <- match.arg(type)
  if (type == "resampling" && is.null(rEGA.obj$raw.data)) {
    warning("Raw data not stored in 'rEGA.obj'. Switching to 'parametric' bootstrap.")
    type <- "parametric"
  }

  seed <- if (is.null(seed)) 11L else as.integer(seed)
  iter <- as.integer(iter)
  ncores <- max(1L, as.integer(ncores))

  dots <- list(...)
  plot.type <- if (!is.null(dots$plot.type)) dots$plot.type else "both"
  layout <- if (!is.null(dots$layout)) dots$layout else "fr"
  palette <- if (!is.null(dots$palette)) dots$palette else "polychrome"
  show.values <- if (!is.null(dots$show.values)) dots$show.values else TRUE
  label.min.pct <- if (!is.null(dots$label.min.pct)) dots$label.min.pct else 10.0

  t0 <- Sys.time()
  n <- rEGA.obj$audit.metadata$n
  p <- rEGA.obj$audit.metadata$p
  R_emp <- rEGA.obj$R.obs
  raw_emp <- rEGA.obj$raw.data
  wc_emp <- rEGA.obj$wc
  item_names <- names(wc_emp)
  cor_method <- rEGA.obj$audit.metadata$cor.method
  algo <- rEGA.obj$audit.metadata$algorithm
  gamma_init <- rEGA.obj$audit.metadata$gamma.used

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  if (ncores > 1L) {
    future::plan(future::multisession, workers = ncores)
  } else {
    future::plan(future::sequential)
  }

  set.seed(seed)
  rep_seeds <- sample.int(1e8, size = iter, replace = FALSE)

  run_one_rep <- function(b) {
    options(mc.cores = 1L)
    s_b <- rep_seeds[b]
    set.seed(s_b)

    if (type == "resampling") {
      idx_b <- sample.int(n, size = n, replace = TRUE)
      dat_b <- raw_emp[idx_b, , drop = FALSE]
      fit_b <- tryCatch(
        rEGA(
          data = dat_b,
          cor.method = cor_method,
          algorithm = algo,
          gamma = gamma_init,
          plot.rEGA = FALSE,
          seed = s_b
        ),
        error = function(e) NULL
      )
    } else {
      sim_b <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = R_emp)
      colnames(sim_b) <- item_names
      R_b <- stats::cor(sim_b)
      fit_b <- tryCatch(
        rEGA(
          data = R_b,
          n = n,
          cor.method = "pearson",
          algorithm = algo,
          gamma = gamma_init,
          plot.rEGA = FALSE,
          seed = s_b
        ),
        error = function(e) NULL
      )
    }

    if (is.null(fit_b)) {
      return(list(converged = FALSE, seed = s_b))
    }

    wc_aligned <- .homogenize_communities(target = wc_emp, estimate = fit_b$wc)

    list(
      converged = TRUE,
      seed = s_b,
      n.dim = fit_b$n.dim,
      wc.aligned = wc_aligned,
      Gamma.t = fit_b$reliability$Gamma.t,
      Psi.G = fit_b$reliability$Psi.G,
      SEM.std = fit_b$reliability$SEM.total.std,
      h2 = fit_b$nodal.triad$h2,
      beta2 = fit_b$nodal.triad$beta2,
      psi2 = fit_b$nodal.triad$psi2
    )
  }

  # Real-time progress bar with Elapsed Time and Estimated Time Remaining (ETA)
  if (isTRUE(progress)) {
    cat(sprintf(
      "\nRunning boot.rEGA (%d replications | type = '%s' | cores = %d | seed = %d)...\n",
      iter, type, ncores, seed
    ))
    boot_list <- vector("list", iter)
    n_batches <- if (ncores == 1L) iter else min(iter, max(10L, ceiling(iter / ncores)))
    batch_indices <- split(seq_len(iter), cut(seq_len(iter), breaks = n_batches, labels = FALSE))
    done_count <- 0L

    for (b_idx in seq_along(batch_indices)) {
      ids <- batch_indices[[b_idx]]
      if (ncores > 1L) {
        res_batch <- future.apply::future_lapply(ids, run_one_rep, future.seed = TRUE)
      } else {
        res_batch <- lapply(ids, run_one_rep)
      }
      boot_list[ids] <- res_batch
      done_count <- done_count + length(ids)

      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      rate <- elapsed / max(1L, done_count)
      eta <- max(0.0, rate * (iter - done_count))
      pct_done <- round((done_count / iter) * 100)
      bar_width <- 30L
      filled <- round((done_count / iter) * bar_width)
      bar_str <- paste0(
        "[",
        paste(rep("=", filled), collapse = ""),
        paste(rep(" ", bar_width - filled), collapse = ""),
        "]"
      )
      cat(sprintf(
        "\r  %s %3d%% (%d/%d) | Elapsed: %5.1fs | ETA: %5.1fs   ",
        bar_str, pct_done, done_count, iter, elapsed, eta
      ))
      utils::flush.console()
    }
    cat("\n")
  } else {
    boot_list <- future.apply::future_lapply(
      seq_len(iter),
      run_one_rep,
      future.seed = TRUE
    )
  }

  valid_reps <- Filter(function(x) isTRUE(x$converged), boot_list)
  B_valid <- length(valid_reps)
  if (B_valid == 0L) {
    stop("No bootstrap replications converged successfully.")
  }

  n_dims_vec <- vapply(valid_reps, `[[`, integer(1), "n.dim")
  dim_tab <- table(n_dims_vec)
  dim_freq <- data.frame(
    n.Dimensions = as.integer(names(dim_tab)),
    Frequency = as.integer(dim_tab),
    Proportion = round(as.numeric(dim_tab) / B_valid, 4),
    stringsAsFactors = FALSE
  )

  wc_mat <- do.call(rbind, lapply(valid_reps, `[[`, "wc.aligned"))
  colnames(wc_mat) <- item_names

  emp_comms <- sort(unique(wc_emp))
  all_comms <- sort(unique(c(emp_comms, as.integer(wc_mat))))

  item_rep_mat <- matrix(0.0, nrow = p, ncol = length(all_comms))
  rownames(item_rep_mat) <- item_names
  colnames(item_rep_mat) <- paste0("Dim_", all_comms)

  for (j in seq_len(p)) {
    for (k_idx in seq_along(all_comms)) {
      c_val <- all_comms[k_idx]
      item_rep_mat[j, k_idx] <- mean(wc_mat[, j] == c_val, na.rm = TRUE)
    }
  }

  item_stab_vec <- numeric(p)
  for (j in seq_len(p)) {
    c_assigned <- wc_emp[j]
    col_key <- paste0("Dim_", c_assigned)
    item_stab_vec[j] <- if (col_key %in% colnames(item_rep_mat)) item_rep_mat[j, col_key] else 0.0
  }

  item_stability_df <- data.frame(
    Item = item_names,
    Empirical_Dim = wc_emp,
    Item_Stability = round(item_stab_vec, 4),
    round(item_rep_mat, 4),
    stringsAsFactors = FALSE
  )
  rownames(item_stability_df) <- NULL

  struct_cons <- vapply(emp_comms, function(c_id) {
    idx_c <- which(wc_emp == c_id)
    mean(apply(wc_mat[, idx_c, drop = FALSE], 1, function(row_c) {
      length(unique(row_c)) == 1L && row_c[1] == c_id
    }))
  }, numeric(1))

  struct_cons_df <- data.frame(
    Dimension = emp_comms,
    Structural_Consistency = round(struct_cons, 4),
    Mean_Item_Stability = round(vapply(emp_comms, function(c_id) {
      mean(item_stab_vec[wc_emp == c_id])
    }, numeric(1)), 4),
    stringsAsFactors = FALSE
  )

  agree_mat <- matrix(0.0, nrow = length(emp_comms), ncol = length(all_comms))
  rownames(agree_mat) <- paste0("Empirical_Dim_", emp_comms)
  colnames(agree_mat) <- paste0("Boot_Dim_", all_comms)

  for (r_i in seq_along(emp_comms)) {
    c_id <- emp_comms[r_i]
    idx_c <- which(wc_emp == c_id)
    agree_mat[r_i, ] <- colMeans(item_rep_mat[idx_c, , drop = FALSE])
  }

  off_diag_drift <- numeric(length(emp_comms))
  for (r_i in seq_along(emp_comms)) {
    c_id <- emp_comms[r_i]
    col_match <- which(all_comms == c_id)
    off_diag_drift[r_i] <- 1.0 - agree_mat[r_i, col_match]
  }
  agreement_df <- data.frame(
    Empirical_Dimension = emp_comms,
    round(agree_mat, 4),
    Off_Diagonal_Dispersion = round(off_diag_drift, 4),
    stringsAsFactors = FALSE
  )
  rownames(agreement_df) <- NULL

  gamma_t_vec <- vapply(valid_reps, `[[`, numeric(1), "Gamma.t")
  psi_g_vec <- vapply(valid_reps, `[[`, numeric(1), "Psi.G")
  sem_std_vec <- vapply(valid_reps, `[[`, numeric(1), "SEM.std")

  global_ci <- data.frame(
    Parameter = c("Gamma_t (Total Reliability)", "Psi_G (Global Noise)", "SEM_std (Standardized SEM)"),
    Empirical = round(c(rEGA.obj$reliability$Gamma.t,
                        rEGA.obj$reliability$Psi.G,
                        rEGA.obj$reliability$SEM.total.std), 4),
    Boot_Mean = round(c(mean(gamma_t_vec), mean(psi_g_vec), mean(sem_std_vec)), 4),
    Boot_SE = round(c(stats::sd(gamma_t_vec), stats::sd(psi_g_vec), stats::sd(sem_std_vec)), 4),
    CI95_Lower = round(c(stats::quantile(gamma_t_vec, 0.025),
                         stats::quantile(psi_g_vec, 0.025),
                         stats::quantile(sem_std_vec, 0.025)), 4),
    CI95_Upper = round(c(stats::quantile(gamma_t_vec, 0.975),
                         stats::quantile(psi_g_vec, 0.975),
                         stats::quantile(sem_std_vec, 0.975)), 4),
    stringsAsFactors = FALSE
  )

  h2_mat <- do.call(rbind, lapply(valid_reps, `[[`, "h2"))
  beta2_mat <- do.call(rbind, lapply(valid_reps, `[[`, "beta2"))
  psi2_mat <- do.call(rbind, lapply(valid_reps, `[[`, "psi2"))

  nodal_ci <- data.frame(
    Item           = item_names,
    h2_Emp         = as.numeric(rEGA.obj$nodal$h2),
    h2_Mean        = apply(h2_mat, 2, mean, na.rm = TRUE),
    h2_SE          = apply(h2_mat, 2, sd, na.rm = TRUE),
    h2_Low         = apply(h2_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    h2_Up          = apply(h2_mat, 2, quantile, probs = 0.975, na.rm = TRUE),
    beta2_Emp      = as.numeric(rEGA.obj$nodal$beta2),
    beta2_Mean     = apply(beta2_mat, 2, mean, na.rm = TRUE),
    beta2_SE       = apply(beta2_mat, 2, sd, na.rm = TRUE),
    beta2_Low      = apply(beta2_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    beta2_Up       = apply(beta2_mat, 2, quantile, probs = 0.975, na.rm = TRUE),
    psi2_Emp       = as.numeric(rEGA.obj$nodal$psi2),
    psi2_Mean      = apply(psi2_mat, 2, mean, na.rm = TRUE),
    psi2_SE        = apply(psi2_mat, 2, sd, na.rm = TRUE),
    psi2_Low       = apply(psi2_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    psi2_Up        = apply(psi2_mat, 2, quantile, probs = 0.975, na.rm = TRUE),
    Item_Stability = item_stab_vec,
    row.names      = NULL,
    stringsAsFactors = FALSE
  )
  rownames(nodal_ci) <- NULL

  boot_log <- data.frame(
    Replication = seq_len(B_valid),
    Seed = vapply(valid_reps, `[[`, integer(1), "seed"),
    n_Dim = n_dims_vec,
    Gamma_t = round(gamma_t_vec, 5),
    Psi_G = round(psi_g_vec, 5),
    SEM_std = round(sem_std_vec, 5),
    stringsAsFactors = FALSE
  )

  out <- list(
    rEGA.obj = rEGA.obj,
    iter = iter,
    converged.iter = B_valid,
    type = type,
    ncores = ncores,
    seed = seed,
    dimension.frequency = dim_freq,
    structural.consistency = struct_cons_df,
    item.stability = item_stability_df,
    agreement.matrix = agreement_df,
    global.ci = global_ci,
    nodal.ci = nodal_ci,
    boot.log = boot_log,
    runtime.sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
  class(out) <- "boot.rEGA"

  if (isTRUE(plot.boot.rEGA)) {
    p_boot <- plot(
      out,
      type = plot.type,
      layout = layout,
      palette = palette,
      show.values = show.values,
      label.min.pct = label.min.pct
    )
    if (inherits(p_boot, "ggplot")) print(p_boot)
  }

  return(out)
}

.homogenize_communities <- function(target, estimate) {
  tab <- table(Estimate = estimate, Target = target)
  est_labels <- as.integer(rownames(tab))
  targ_labels <- as.integer(colnames(tab))
  aligned <- estimate

  used_targets <- integer(0)
  order_est <- order(rowSums(tab), decreasing = TRUE)

  for (r_idx in order_est) {
    e_val <- est_labels[r_idx]
    best_t_idx <- which.max(tab[r_idx, ])
    best_t_val <- targ_labels[best_t_idx]
    if (!(best_t_val %in% used_targets) && tab[r_idx, best_t_idx] > 0) {
      aligned[estimate == e_val] <- best_t_val
      used_targets <- c(used_targets, best_t_val)
    } else {
      new_id <- max(c(targ_labels, used_targets)) + 1L
      aligned[estimate == e_val] <- new_id
      used_targets <- c(used_targets, new_id)
    }
  }
  return(as.integer(aligned))
}
