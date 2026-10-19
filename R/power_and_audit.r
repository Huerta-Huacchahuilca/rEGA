#' @title Planificacion de Potencia y Tamano Muestral con Atenuacion para rEGA
#'
#' @description
#' Calcula la fiabilidad proyectada de red (\code{Gamma_p}) segun la formula de
#' Spearman-Brown y estima el tamano muestral requerido (\code{N}) considerando
#' la atenuacion por error de medicion (Henry, 2026).
#'
#' @param rEGA.obj Objeto opcional de clase \code{"rEGA"}.
#' @param gamma.baseline Numerico. Fiabilidad base si \code{rEGA.obj = NULL}. Por defecto \code{0.80}.
#' @param p.baseline Entero. Numero de items base. Por defecto \code{7L}.
#' @param p.target Vector entero de cantidades de items a evaluar. Por defecto \code{c(5L, 7L, 10L, 15L)}.
#' @param r.true Numerico. Correlacion verdadera esperada. Por defecto \code{0.30}.
#' @param gamma.criterion Numerico. Fiabilidad del criterio externo. Por defecto \code{0.85}.
#' @param power Numerico. Potencia estadistica deseada. Por defecto \code{0.80}.
#' @param alpha Numerico. Nivel de significancia. Por defecto \code{0.05}.
#'
#' @return Un \code{data.frame} con la fiabilidad proyectada, correlacion observable atenuada y \code{N} requerido.
#'
#' @importFrom pwr pwr.r.test
#' @export
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

#' @title Exportacion de Libro Excel de Auditoria Criptografica para rEGA
#'
#' @description
#' Genera un archivo Excel (\code{.xlsx}) de 5 hojas con el manifiesto de reproducibilidad,
#' hash SHA-256, matrices por etapa, triada nodal, fiabilidad estructural e intervalos bootstrap.
#'
#' @param rEGA.obj Objeto de clase \code{"rEGA"}.
#' @param boot.obj Objeto opcional de clase \code{"boot.rEGA"}.
#' @param file Ruta del archivo Excel de salida. Por defecto \code{"rEGA_OpenScience_Audit.xlsx"}.
#'
#' @return Devuelve de forma invisible la ruta normalizada del archivo Excel.
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
      "Random Seed Used",
      "Min Eigenvalue of K_refit",
      "Max Resolvent Algebraic Error",
      "Max Nodal Triad Conservation Error",
      "Data Cryptographic Hash (SHA-256)",
      "Computation Time (seconds)"
    ),
    Recorded_Value = c(
      as.character(Sys.time()),
      paste0(R.version$version.string, " - ", R.version$platform),
      as.character(utils::packageVersion("rEGA")),
      rEGA.obj$audit.metadata$pathway,
      rEGA.obj$audit.metadata$n,
      rEGA.obj$audit.metadata$p,
      rEGA.obj$audit.metadata$cor.method,
      rEGA.obj$audit.metadata$algorithm,
      rEGA.obj$audit.metadata$gamma.used,
      rEGA.obj$audit.metadata$active.edges,
      if (!is.null(rEGA.obj$audit.metadata$seed)) rEGA.obj$audit.metadata$seed else 11L,
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
                              "1. Nodal Variance Conservation Triad (h2 + beta2 + psi2 = 1)")
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
                                "3. Partition Agreement Matrix (Empirical vs. Bootstrap)")
    write_sheet_table("05_Bootstrap_Stability", boot.obj$nodal.ci, r_next,
                      "4. Nodal Triad 95% Bootstrap Confidence Intervals & Item Stability")
  }

  openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
  invisible(normalizePath(file, mustWork = FALSE))
}
