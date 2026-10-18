#' @title Attenuation-Aware Power and Sample Size Planning for rEGA
#'
#' @description
#' Implements the 3-step psychometric power and sample size workflow of Henry (2026)
#' adapted to \code{rEGA} network reliability (\code{Gamma_c} or \code{Gamma_t}):
#' \enumerate{
#'   \item Projects expected Network Reliability (\code{Gamma_p}) when modifying
#'         the number of items from baseline \code{p.baseline} to \code{p.target}.
#'   \item Computes the attenuated observable correlation
#'         (\code{r_obs = r_true * sqrt(Gamma_x * Gamma_y)}) and required sample
#'         size \code{N} for target statistical power (\code{1 - beta}).
#' }
#'
#' @param rEGA.obj Optional object of class \code{"rEGA"} to automatically extract
#'   baseline reliability \code{Gamma.t} and item count \code{p}.
#' @param gamma.baseline Numeric. Baseline network reliability if \code{rEGA.obj = NULL}.
#'   Default is \code{0.80}.
#' @param p.baseline Integer. Baseline number of items if \code{rEGA.obj = NULL}.
#'   Default is \code{7L}.
#' @param p.target Integer vector. Candidate target item counts to evaluate.
#'   Default is \code{c(5L, 7L, 10L, 15L)}.
#' @param r.true Numeric. Expected true (unattenuated) effect size correlation
#'   with an external criterion. Default is \code{0.30}.
#' @param gamma.criterion Numeric. Reliability of the external criterion measure.
#'   Default is \code{0.85}.
#' @param power Numeric. Target statistical power. Default is \code{0.80}.
#' @param alpha Numeric. Type I error rate. Default is \code{0.05}.
#'
#' @return A \code{data.frame} containing projected reliability,
#'   attenuated observable correlations, required sample sizes \code{N}, and sample
#'   size inflation due to measurement error.
#'
#' @importFrom pwr pwr.r.test
#' @export
#'
#' @examples
#' power_plan <- rEGA.power(gamma.baseline = 0.85, p.baseline = 7, r.true = 0.30)
#' print(power_plan)
rEGA.power <- function(rEGA.obj = NULL,
                       gamma.baseline = 0.80,
                       p.baseline = 7L,
                       p.target = c(5L, 7L, 10L, 15L),
                       r.true = 0.30,
                       gamma.criterion = 0.85,
                       power = 0.80,
                       alpha = 0.05) {

  if (!is.null(rEGA.obj) && inherits(rEGA.obj, "rEGA")) {
    gamma.baseline <- rEGA.obj$reliability$Gamma.t
    p.baseline <- rEGA.obj$audit.metadata$p
  }

  n_unattenuated <- ceiling(pwr::pwr.r.test(
    r = r.true,
    sig.level = alpha,
    power = power,
    alternative = "two.sided"
  )$n)

  rows <- lapply(p.target, function(p_new) {
    k_factor <- p_new / p.baseline
    gamma_proj <- (k_factor * gamma.baseline) / (1.0 + (k_factor - 1.0) * gamma.baseline)
    r_obs <- r.true * sqrt(gamma_proj * gamma.criterion)

    n_req <- ceiling(pwr::pwr.r.test(
      r = r_obs,
      sig.level = alpha,
      power = power,
      alternative = "two.sided"
    )$n)

    data.frame(
      Items_Target = p_new,
      Projected_Gamma = round(gamma_proj, 4),
      True_r = r.true,
      Criterion_Reliability = gamma.criterion,
      Observable_Attenuated_r = round(r_obs, 4),
      Attenuation_Pct = round((1.0 - (r_obs / r.true)) * 100, 2),
      Required_N = n_req,
      Unattenuated_Ideal_N = n_unattenuated,
      Extra_N_For_Measurement_Error = n_req - n_unattenuated,
      stringsAsFactors = FALSE
    )
  })

  res_df <- do.call(rbind, rows)
  return(res_df)
}

#' @title Export Open Science Cryptographic Audit Workbook for rEGA
#'
#' @description
#' Generates a multi-sheet Excel audit workbook (\code{.xlsx}) containing the Open Science
#' computational manifest (R version, package versions, SHA-256 cryptographic data hash),
#' stage-by-stage matrices (\code{R.obs}, \code{K.refit}, \code{network.refit},
#' \code{resolvent}, \code{adjacency.rEGA}), the Nodal Variance Conservation Triad,
#' Network Loadings, Structural Reliability (\code{Gamma_c}, \code{SEM_c}), and
#' bootstrap inference tables if provided.
#'
#' @param rEGA.obj An object of class \code{"rEGA"}.
#' @param boot.obj Optional object of class \code{"boot.rEGA"}.
#' @param file Character. Output file path for the Excel workbook.
#'   Default is \code{"rEGA_OpenScience_Audit.xlsx"}.
#'
#' @return Invisibly returns the normalized path to the generated Excel file.
#'
#' @importFrom openxlsx createWorkbook addWorksheet writeData createStyle addStyle setColWidths saveWorkbook
#' @export
rEGA.audit <- function(rEGA.obj,
                       boot.obj = NULL,
                       file = "rEGA_OpenScience_Audit.xlsx") {

  if (!inherits(rEGA.obj, "rEGA")) {
    stop("'rEGA.obj' must be an object of class 'rEGA'.")
  }

  wb <- openxlsx::createWorkbook()
  hdr_style <- openxlsx::createStyle(
    fontColour = "#FFFFFF", fgFill = "#1B365D",
    textDecoration = "bold", border = "TopBottomLeftRight"
  )

  write_sheet_table <- function(sheet, df, start_row = 1, title = NULL) {
    if (!is.null(title)) {
      openxlsx::writeData(wb, sheet = sheet, x = title, startRow = start_row, startCol = 1)
      start_row <- start_row + 1L
    }
    openxlsx::writeData(wb, sheet = sheet, x = df, startRow = start_row,
                        startCol = 1, headerStyle = hdr_style, borders = "rows")
    openxlsx::setColWidths(wb, sheet = sheet, cols = seq_len(ncol(df)), widths = "auto")
    return(start_row + nrow(df) + 3L)
  }

  mat_to_df <- function(M) {
    df <- as.data.frame(round(M, 6))
    cbind(Item = rownames(M), df, stringsAsFactors = FALSE)
  }

  openxlsx::addWorksheet(wb, "01_OpenScience_Manifest")
  manifest_df <- data.frame(
    Audit_Attribute = c(
      "Execution Timestamp",
      "R Version & Platform",
      "rEGA Package Version",
      "Estimation Pathway",
      "Sample Size (n)",
      "Number of Items (p)",
      "Correlation Method",
      "Community Algorithm",
      "EBIC Gamma Used",
      "Active Edges (E_lambda)",
      "Min Eigenvalue of K_refit",
      "Max Resolvent Algebraic Error",
      "Max Nodal Triad Conservation Error",
      "Data Cryptographic Hash (SHA-256)",
      "Computation Time (seconds)"
    ),
    Recorded_Value = c(
      as.character(Sys.time()),
      paste0(R.version$version.string, " ert{} ", R.version$platform),
      as.character(utils::packageVersion("rEGA")),
      rEGA.obj$audit.metadata$pathway,
      rEGA.obj$audit.metadata$n,
      rEGA.obj$audit.metadata$p,
      rEGA.obj$audit.metadata$cor.method,
      rEGA.obj$audit.metadata$algorithm,
      rEGA.obj$audit.metadata$gamma.used,
      rEGA.obj$audit.metadata$active.edges,
      round(rEGA.obj$audit.metadata$min.eigen.Krefit, 8),
      format(rEGA.obj$audit.metadata$max.algebraic.error, scientific = TRUE),
      format(rEGA.obj$audit.metadata$max.conservation.error, scientific = TRUE),
      rEGA.obj$audit.metadata$data.sha256,
      round(rEGA.obj$audit.metadata$runtime.sec, 4)
    ),
    stringsAsFactors = FALSE
  )
  write_sheet_table("01_OpenScience_Manifest", manifest_df, 1,
                    "rEGA OPEN SCIENCE & REPRODUCIBILITY MANIFEST")

  openxlsx::addWorksheet(wb, "02_Nodal_Triad_Loadings")
  r_next <- write_sheet_table("02_Nodal_Triad_Loadings", rEGA.obj$nodal.triad, 1,
                              "1. Nodal Variance Conservation Triad (h2 + beta2 + psi2 = 1 | 100%)")
  write_sheet_table("02_Nodal_Triad_Loadings", mat_to_df(rEGA.obj$net.loadings), r_next,
                    "2. Standardized Network Loadings")

  openxlsx::addWorksheet(wb, "03_Reliability_Gamma_SEM")
  r_next <- write_sheet_table("03_Reliability_Gamma_SEM", rEGA.obj$reliability$dimensions, 1,
                              "1. Dimensional Reliability (Coefficient Gamma_c & Graph SEM_c)")
  glob_df <- data.frame(
    Total_Gamma_t = round(rEGA.obj$reliability$Gamma.t, 6),
    Global_Noise_Psi_G = round(rEGA.obj$reliability$Psi.G, 6),
    Total_SEM_Std = round(rEGA.obj$reliability$SEM.total.std, 6),
    Total_SEM_Raw = round(rEGA.obj$reliability$SEM.total.raw, 6)
  )
  write_sheet_table("03_Reliability_Gamma_SEM", glob_df, r_next,
                    "2. Global Network Reliability Summary")

  openxlsx::addWorksheet(wb, "04_Resolvent_Matrices")
  r_next <- write_sheet_table("04_Resolvent_Matrices", mat_to_df(rEGA.obj$R.obs), 1,
                              "1. Observed Correlation Matrix (R_obs)")
  r_next <- write_sheet_table("04_Resolvent_Matrices", mat_to_df(rEGA.obj$K.refit), r_next,
                              "2. Support-Refitted Precision Matrix (K_refit)")
  r_next <- write_sheet_table("04_Resolvent_Matrices", mat_to_df(rEGA.obj$network.refit), r_next,
                              "3. Unbiased Partial Correlation Matrix (Omega_refit)")
  r_next <- write_sheet_table("04_Resolvent_Matrices", mat_to_df(rEGA.obj$resolvent), r_next,
                              "4. Neumann-Markov Resolvent Operator Q = (I - Omega)^(-1)")
  write_sheet_table("04_Resolvent_Matrices", mat_to_df(rEGA.obj$adjacency.rEGA), r_next,
                    "5. Purity-Weighted Adjacency Matrix (A^rEGA)")

  if (!is.null(boot.obj) && inherits(boot.obj, "boot.rEGA")) {
    openxlsx::addWorksheet(wb, "05_Bootstrap_Stability")
    r_next <- write_sheet_table("05_Bootstrap_Stability", boot.obj$global.ci, 1,
                                "1. Global Reliability & SEM 95% Bootstrap Confidence Intervals")
    r_next <- write_sheet_table("05_Bootstrap_Stability", boot.obj$structural.consistency, r_next,
                                "2. Structural Consistency by Dimension")
    r_next <- write_sheet_table("05_Bootstrap_Stability", boot.obj$agreement.matrix, r_next,
                                "3. Diagnostic Agreement Matrix (Holistic evaluation, no rigid cutoffs)")
    write_sheet_table("05_Bootstrap_Stability", boot.obj$nodal.ci, r_next,
                      "4. Nodal Triad 95% Bootstrap Confidence Intervals & Item Stability")
  }

  openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
  invisible(normalizePath(file, mustWork = FALSE))
}
