# Declare NSE variables used in ggplot2::aes() to satisfy R CMD check
utils::globalVariables(c(
  "Item", "Proportion", "Component", "Item_Stability", "Dim_Label",
  "Pct_Label", "x", "y", "xend", "yend", "Weight",
  "Sign", "Node", "Stability_Label"
))

# Palettes matching EGAnet: "polychrome" (default EGA) and "colorblind" (Okabe-Ito)
.rEGA_community_palette <- function(n_comms, palette = c("polychrome", "colorblind")) {
  palette <- match.arg(palette)
  if (palette == "colorblind") {
    cols <- c(
      "#0072B2", "#E69F00", "#009E73", "#CC79A7",
      "#56B4E9", "#D55E00", "#F0E442", "#999999"
    )
  } else {
    # Classic EGAnet / ggplot2 polychrome palette
    cols <- c(
      "#F8766D", "#00BA38", "#619CFF", "#F564E3",
      "#B79F00", "#00BFC4", "#00B0F6", "#E76BF3"
    )
  }
  if (n_comms <= length(cols)) {
    return(cols[seq_len(n_comms)])
  }
  grDevices::hcl.colors(n_comms, palette = "Dynamic")
}

.rEGA_triad_palette <- function(palette = c("polychrome", "colorblind")) {
  palette <- match.arg(palette)
  if (palette == "colorblind") {
    c(
      "Within-Community Signal (h2)" = "#0072B2",
      "Bridging Variance (beta2)" = "#E69F00",
      "Markov Uniqueness (psi2)" = "#D55E00"
    )
  } else {
    c(
      "Within-Community Signal (h2)" = "#2B5C8F",
      "Bridging Variance (beta2)" = "#E68A00",
      "Markov Uniqueness (psi2)" = "#C0392B"
    )
  }
}

.build_network_ggplot <- function(rEGA.obj,
                                  layout = c("fr", "circle"),
                                  palette = c("polychrome", "colorblind")) {
  layout <- match.arg(layout)
  palette <- match.arg(palette)

  Omega <- rEGA.obj$network.refit
  wc <- rEGA.obj$wc
  p <- ncol(Omega)
  item_names <- colnames(Omega)
  seed_val <- if (!is.null(rEGA.obj$audit.metadata$seed)) rEGA.obj$audit.metadata$seed else 11L

  A_abs <- abs(Omega)
  diag(A_abs) <- 0
  g <- igraph::graph_from_adjacency_matrix(A_abs, mode = "undirected", weighted = TRUE, diag = FALSE)
  set.seed(seed_val)
  if (layout == "circle" || sum(A_abs) <= 1e-6) {
    coords <- igraph::layout_in_circle(g)
  } else {
    coords <- igraph::layout_with_fr(g, weights = igraph::E(g)$weight)
  }

  for (col_i in 1:2) {
    rng <- range(coords[, col_i])
    if (diff(rng) > 1e-6) {
      coords[, col_i] <- (coords[, col_i] - rng[1]) / diff(rng) * 2 - 1
    } else {
      coords[, col_i] <- 0
    }
  }

  unique_dims <- sort(unique(wc))
  nodes_df <- data.frame(
    Node = item_names,
    x = coords[, 1],
    y = coords[, 2],
    Dim_Label = factor(paste0("Dimension ", wc),
                       levels = paste0("Dimension ", unique_dims)),
    stringsAsFactors = FALSE
  )

  upper_idx <- which(upper.tri(Omega) & abs(Omega) > 1e-5, arr.ind = TRUE)
  if (nrow(upper_idx) > 0) {
    edges_df <- data.frame(
      x = coords[upper_idx[, 1], 1],
      y = coords[upper_idx[, 1], 2],
      xend = coords[upper_idx[, 2], 1],
      yend = coords[upper_idx[, 2], 2],
      Weight = abs(Omega[upper_idx]),
      Sign = ifelse(Omega[upper_idx] >= 0, "Positive", "Negative"),
      stringsAsFactors = FALSE
    )
  } else {
    edges_df <- data.frame(
      x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0),
      Weight = numeric(0), Sign = character(0)
    )
  }

  pal <- stats::setNames(
    .rEGA_community_palette(length(unique_dims), palette = palette),
    paste0("Dimension ", unique_dims)
  )
  edge_cols <- if (palette == "colorblind") {
    c("Positive" = "#009E73", "Negative" = "#D55E00")
  } else {
    c("Positive" = "#1B7837", "Negative" = "#C51B7D")
  }

  p_net <- ggplot2::ggplot()
  if (nrow(edges_df) > 0) {
    p_net <- p_net +
      ggplot2::geom_segment(
        data = edges_df,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                     linewidth = Weight, alpha = Weight, linetype = Sign, color = Sign),
        show.legend = FALSE
      ) +
      ggplot2::scale_linewidth_continuous(range = c(0.4, 3.2)) +
      ggplot2::scale_alpha_continuous(range = c(0.30, 0.88)) +
      ggplot2::scale_color_manual(values = edge_cols) +
      ggplot2::scale_linetype_manual(values = c("Positive" = "solid", "Negative" = "dashed"))
  }

  p_net <- p_net +
    ggplot2::geom_point(
      data = nodes_df,
      ggplot2::aes(x = x, y = y, fill = Dim_Label),
      shape = 21, size = 11.5, color = "grey15", stroke = 0.85
    ) +
    ggplot2::geom_text(
      data = nodes_df,
      ggplot2::aes(x = x, y = y, label = Node),
      color = "grey10", fontface = "bold", size = 3.3
    ) +
    ggplot2::scale_fill_manual(values = pal, name = "Dimension") +
    ggplot2::coord_equal(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
    ggplot2::labs(
      title = sprintf("rEGA Network (%d %s)", rEGA.obj$n.dim,
                      ifelse(rEGA.obj$n.dim == 1, "Dimension", "Dimensions")),
      subtitle = sprintf("Gamma_t = %.3f | Psi_G = %.3f",
                         rEGA.obj$reliability$Gamma.t, rEGA.obj$reliability$Psi.G),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey30")
    )

  return(p_net)
}

.build_variance_ggplot <- function(rEGA.obj,
                                   palette = c("polychrome", "colorblind"),
                                   show.values = TRUE,
                                   label.min.pct = 10.0) {
  palette <- match.arg(palette)
  triad <- rEGA.obj$nodal.triad

  tot_i <- triad$h2 + triad$beta2 + triad$psi2
  pct_h2 <- round((triad$h2 / tot_i) * 100, 1)
  pct_b2 <- round((triad$beta2 / tot_i) * 100, 1)
  pct_p2 <- round(100.0 - pct_h2 - pct_b2, 1)

  df_long <- rbind(
    data.frame(
      Item = triad$Item, Community = triad$Community,
      Component = "Within-Community Signal (h2)",
      Proportion = triad$h2 / tot_i,
      Pct_Label = ifelse(pct_h2 >= label.min.pct, sprintf("%.1f%%", pct_h2), ""),
      stringsAsFactors = FALSE
    ),
    data.frame(
      Item = triad$Item, Community = triad$Community,
      Component = "Bridging Variance (beta2)",
      Proportion = triad$beta2 / tot_i,
      Pct_Label = ifelse(pct_b2 >= label.min.pct, sprintf("%.1f%%", pct_b2), ""),
      stringsAsFactors = FALSE
    ),
    data.frame(
      Item = triad$Item, Community = triad$Community,
      Component = "Markov Uniqueness (psi2)",
      Proportion = triad$psi2 / tot_i,
      Pct_Label = ifelse(pct_p2 >= label.min.pct, sprintf("%.1f%%", pct_p2), ""),
      stringsAsFactors = FALSE
    )
  )

  # Clean Y-axis showing only item names
  df_long$Item <- factor(df_long$Item, levels = rev(triad$Item))
  df_long$Component <- factor(df_long$Component, levels = c(
    "Markov Uniqueness (psi2)",
    "Bridging Variance (beta2)",
    "Within-Community Signal (h2)"
  ))

  p_var <- ggplot2::ggplot(df_long, ggplot2::aes(x = Item, y = Proportion, fill = Component)) +
    ggplot2::geom_col(width = 0.70, color = "white", linewidth = 0.25) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      labels = function(v) paste0(round(v * 100), "%"),
      expand = c(0.01, 0.01)
    ) +
    ggplot2::scale_fill_manual(values = .rEGA_triad_palette(palette)) +
    ggplot2::labs(
      title = "Nodal Variance Triad Decomposition",
      subtitle = sprintf("h2 + beta2 + psi2 = 100%% | Gamma_t = %.3f | Psi_G = %.3f",
                         rEGA.obj$reliability$Gamma.t, rEGA.obj$reliability$Psi.G),
      x = "Item",
      y = "Standardized Nodal Variance (%)",
      fill = "Component"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey30"),
      axis.text.y = ggplot2::element_text(face = "bold", size = 10)
    )

  if (isTRUE(show.values)) {
    p_var <- p_var +
      ggplot2::geom_text(
        ggplot2::aes(label = Pct_Label),
        position = ggplot2::position_stack(vjust = 0.5),
        color = "white", fontface = "bold", size = 3.2
      )
  }

  return(p_var)
}

.build_stability_ggplot <- function(boot.obj,
                                    palette = c("polychrome", "colorblind"),
                                    show.values = TRUE) {
  palette <- match.arg(palette)
  df <- boot.obj$item.stability
  df$Item <- factor(df$Item, levels = rev(df$Item))
  unique_dims <- sort(unique(df$Empirical_Dim))
  df$Dim_Label <- factor(paste0("Dimension ", df$Empirical_Dim),
                         levels = paste0("Dimension ", unique_dims))
  df$Stability_Label <- sprintf("%.1f%%", df$Item_Stability * 100)

  pal <- stats::setNames(
    .rEGA_community_palette(length(unique_dims), palette = palette),
    paste0("Dimension ", unique_dims)
  )

  p_stab <- ggplot2::ggplot(df, ggplot2::aes(x = Item, y = Item_Stability, color = Dim_Label)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = Item, xend = Item, y = 0, yend = Item_Stability),
      linewidth = 1.2, alpha = 0.85
    ) +
    ggplot2::geom_point(size = 4.2) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.12),
      breaks = seq(0, 1.0, by = 0.20),
      labels = function(v) paste0(round(v * 100), "%")
    ) +
    ggplot2::scale_color_manual(values = pal, name = "Dimension") +
    ggplot2::labs(
      title = "Bootstrap Item Stability",
      subtitle = sprintf("Replication Proportion (%d Bootstrap Replications)", boot.obj$converged.iter),
      x = "Item",
      y = "Replication Proportion (%)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey30"),
      axis.text.y = ggplot2::element_text(face = "bold", size = 10)
    )

  if (isTRUE(show.values)) {
    p_stab <- p_stab +
      ggplot2::geom_text(
        ggplot2::aes(label = Stability_Label),
        hjust = -0.25, size = 3.2, fontface = "bold", show.legend = FALSE
      )
  }

  return(p_stab)
}

.draw_side_by_side <- function(plot_list) {
  n_plots <- length(plot_list)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, n_plots)))
  for (i in seq_len(n_plots)) {
    print(plot_list[[i]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
  }
  invisible(plot_list)
}

#' @title Plot Method for rEGA Objects
#' @description
#' Visualizes the EGA-style Purity-Weighted Partial Correlation Network and the Nodal Variance
#' Conservation Triad either combined side-by-side (\code{type = "both"}), as separate plots
#' (\code{type = "separate"}), or individually (\code{"network"} or \code{"variance"}).
#' @param x An object of class \code{"rEGA"}.
#' @param type Character. Display mode: \code{"both"} (default, combined side-by-side),
#'   \code{"separate"} (prints two separate figures sequentially), \code{"network"} (network only),
#'   or \code{"variance"} (Nodal Variance Triad only).
#' @param layout Character. Network layout algorithm: \code{"fr"} (default, Fruchterman-Reingold)
#'   or \code{"circle"} (circular layout).
#' @param palette Character. Color palette matching EGAnet: \code{"polychrome"} (default) or
#'   \code{"colorblind"} (Okabe-Ito colorblind-safe palette).
#' @param show.values Logical. If \code{TRUE} (default), displays standardized percentages inside
#'   each bar rectangle when the component is at least \code{label.min.pct}.
#' @param label.min.pct Numeric. Minimum percentage required to print the label inside a bar
#'   segment so text never overflows narrow rectangles. Default is \code{10}.
#' @param ... Additional arguments (unused).
#' @importFrom ggplot2 ggplot aes geom_col geom_point geom_segment geom_text coord_flip coord_equal scale_fill_manual scale_color_manual scale_linewidth_continuous scale_alpha_continuous scale_linetype_manual scale_y_continuous labs theme_minimal theme element_text element_blank position_stack
#' @importFrom grid grid.newpage pushViewport viewport grid.layout
#' @importFrom grDevices hcl.colors
#' @export
plot.rEGA <- function(x,
                      type = c("both", "separate", "variance", "network"),
                      layout = c("fr", "circle"),
                      palette = c("polychrome", "colorblind"),
                      show.values = TRUE,
                      label.min.pct = 10.0,
                      ...) {
  type <- match.arg(type)
  layout <- match.arg(layout)
  palette <- match.arg(palette)

  if (type == "network") {
    return(.build_network_ggplot(x, layout = layout, palette = palette))
  }
  if (type == "variance") {
    return(.build_variance_ggplot(x, palette = palette, show.values = show.values, label.min.pct = label.min.pct))
  }

  p_net <- .build_network_ggplot(x, layout = layout, palette = palette)
  p_var <- .build_variance_ggplot(x, palette = palette, show.values = show.values, label.min.pct = label.min.pct)

  if (type == "separate") {
    print(p_net)
    print(p_var)
    return(invisible(list(Network = p_net, Variance_Triad = p_var)))
  }

  return(.draw_side_by_side(list(Network = p_net, Variance_Triad = p_var)))
}

#' @title Plot Method for boot.rEGA Objects
#' @description
#' Visualizes the empirical \code{rEGA} network and the bootstrap item stability plot either
#' combined side-by-side (\code{type = "both"}), as separate figures (\code{type = "separate"}),
#' or individually (\code{"stability"}, \code{"network"}, \code{"variance"}, or \code{"all"}).
#' @param x An object of class \code{"boot.rEGA"}.
#' @param type Character. Display mode: \code{"both"} (default, combined side-by-side Network +
#'   Item Stability), \code{"separate"} (prints Network and Stability as separate figures),
#'   \code{"stability"} (Item Stability only), \code{"network"} (Network only),
#'   \code{"variance"} (Nodal Variance Triad only), or \code{"all"} (3-panel figure).
#' @param layout Character. Network layout algorithm: \code{"fr"} (default, Fruchterman-Reingold)
#'   or \code{"circle"} (circular layout).
#' @param palette Character. Color palette matching EGAnet: \code{"polychrome"} (default) or
#'   \code{"colorblind"} (Okabe-Ito colorblind-safe palette).
#' @param show.values Logical. Display percentages on items? Default is \code{TRUE}.
#' @param label.min.pct Numeric. Minimum percentage required to print labels inside Nodal Triad bars.
#'   Default is \code{10}.
#' @param ... Additional arguments (unused).
#' @export
plot.boot.rEGA <- function(x,
                           type = c("both", "separate", "stability", "network", "variance", "all"),
                           layout = c("fr", "circle"),
                           palette = c("polychrome", "colorblind"),
                           show.values = TRUE,
                           label.min.pct = 10.0,
                           ...) {
  type <- match.arg(type)
  layout <- match.arg(layout)
  palette <- match.arg(palette)

  p_stab <- .build_stability_ggplot(x, palette = palette, show.values = show.values)

  if (type == "stability" || is.null(x$rEGA.obj)) {
    return(p_stab)
  }
  if (type == "network") {
    return(.build_network_ggplot(x$rEGA.obj, layout = layout, palette = palette))
  }
  if (type == "variance") {
    return(.build_variance_ggplot(x$rEGA.obj, palette = palette, show.values = show.values, label.min.pct = label.min.pct))
  }

  p_net <- .build_network_ggplot(x$rEGA.obj, layout = layout, palette = palette)

  if (type == "separate") {
    print(p_net)
    print(p_stab)
    return(invisible(list(Network = p_net, Stability = p_stab)))
  }

  if (type == "all") {
    p_var <- .build_variance_ggplot(x$rEGA.obj, palette = palette, show.values = show.values, label.min.pct = label.min.pct)
    return(.draw_side_by_side(list(Network = p_net, Variance_Triad = p_var, Stability = p_stab)))
  }

  return(.draw_side_by_side(list(Network = p_net, Stability = p_stab)))
}

#' @title Print Method for rEGA Objects
#' @param x An object of class \code{"rEGA"}.
#' @param ... Additional arguments (unused).
#' @export
print.rEGA <- function(x, ...) {
  cat("\n====================================================================\n")
  cat("  Reliability Exploratory Graph Analysis (rEGA)\n")
  cat("====================================================================\n")
  cat(sprintf("  Data Pathway                : %s\n", x$audit.metadata$pathway))
  cat(sprintf("  Sample Size (n)             : %d\n", x$audit.metadata$n))
  cat(sprintf("  Number of Items (p)         : %d\n", x$audit.metadata$p))
  cat(sprintf("  Correlation Method          : %s\n", x$audit.metadata$cor.method))
  cat(sprintf("  Community Algorithm         : %s\n", x$audit.metadata$algorithm))
  cat(sprintf("  EBIC Tuning Gamma           : %.2f\n", x$audit.metadata$gamma.used))
  cat(sprintf("  Active Network Edges        : %d\n", x$audit.metadata$active.edges))
  cat(sprintf("  Random Seed Used            : %d\n", x$audit.metadata$seed))
  cat(sprintf("  Estimated Dimensions        : %d\n", x$n.dim))
  cat("--------------------------------------------------------------------\n")
  cat("  Global Network Reliability & Measurement Error:\n")
  cat(sprintf("  * Total Network Gamma (Gamma_t) : %.4f\n", x$reliability$Gamma.t))
  cat(sprintf("  * Global Graph Noise (Psi_G)    : %.4f\n", x$reliability$Psi.G))
  cat(sprintf("  * Standardized Graph SEM        : %.4f\n", x$reliability$SEM.total.std))
  cat(sprintf("  * Raw-Score Graph SEM           : %.4f\n", x$reliability$SEM.total.raw))
  cat("--------------------------------------------------------------------\n")
  cat("  1. Dimensional Structural Reliability & Mean Nodal Triad:\n\n")
  dim_tab <- data.frame(
    Dim = x$reliability$dimensions$Dimension,
    Items = x$reliability$dimensions$n_Items,
    Gamma_c = sprintf("%.4f", x$reliability$dimensions$Gamma_c),
    Mean_h2 = sprintf("%.4f", x$reliability$dimensions$Mean_h2),
    Mean_beta2 = sprintf("%.4f", x$reliability$dimensions$Mean_beta2),
    Mean_psi2 = sprintf("%.4f", x$reliability$dimensions$Mean_psi2),
    SEM_std = sprintf("%.4f", x$reliability$dimensions$SEM_std),
    SEM_raw = sprintf("%.4f", x$reliability$dimensions$SEM_raw),
    stringsAsFactors = FALSE
  )
  print(dim_tab, row.names = FALSE)
  cat("--------------------------------------------------------------------\n")
  cat("  2. Item-Level Nodal Variance Conservation Triad (Proportions & %):\n\n")
  triad_tab <- data.frame(
    Item = x$nodal.triad$Item,
    Dim = x$nodal.triad$Community,
    h2 = sprintf("%.4f", x$nodal.triad$h2),
    beta2 = sprintf("%.4f", x$nodal.triad$beta2),
    psi2 = sprintf("%.4f", x$nodal.triad$psi2),
    h2_Pct = sprintf("%5.1f%%", x$nodal.triad$h2_pct),
    beta2_Pct = sprintf("%5.1f%%", x$nodal.triad$beta2_pct),
    psi2_Pct = sprintf("%5.1f%%", x$nodal.triad$psi2_pct),
    stringsAsFactors = FALSE
  )
  print(triad_tab, row.names = FALSE)
  cat("====================================================================\n")
  invisible(x)
}

#' @title Summary Method for rEGA Objects
#' @param object An object of class \code{"rEGA"}.
#' @param ... Additional arguments (unused).
#' @export
summary.rEGA <- function(object, ...) {
  print.rEGA(object, ...)
  cat("  3. Standardized Network Loadings:\n\n")
  print(object$net.loadings)
  cat("--------------------------------------------------------------------\n")
  cat(sprintf("  Max Resolvent Algebraic Error : %.2e\n", object$audit.metadata$max.algebraic.error))
  cat(sprintf("  Max Triad Conservation Error  : %.2e\n", object$audit.metadata$max.conservation.error))
  cat(sprintf("  Data SHA-256 Hash (first 16)  : %s\n", substr(object$audit.metadata$data.sha256, 1, 16)))
  cat("====================================================================\n")
  invisible(object)
}

#' @title Print Method for boot.rEGA Objects
#' @param x An object of class \code{"boot.rEGA"}.
#' @param ... Additional arguments (unused).
#' @export
print.boot.rEGA <- function(x, ...) {
  cat("\n====================================================================\n")
  cat("  Bootstrap Inference & Stability Analysis for rEGA (boot.rEGA)\n")
  cat("====================================================================\n")
  cat(sprintf("  Resampling Scheme           : %s\n", x$type))
  cat(sprintf("  Converged Replications      : %d / %d\n", x$converged.iter, x$iter))
  cat(sprintf("  Parallel CPU Cores          : %d\n", x$ncores))
  cat(sprintf("  Random Seed Used            : %d\n", x$seed))
  cat(sprintf("  Computation Time            : %.2f sec\n", x$runtime.sec))
  cat("--------------------------------------------------------------------\n")
  cat("  1. Frequency of Estimated Dimensions:\n\n")
  print(x$dimension.frequency, row.names = FALSE)
  cat("--------------------------------------------------------------------\n")
  cat("  2. Structural Consistency & Mean Item Stability by Dimension:\n\n")
  print(x$structural.consistency, row.names = FALSE)
  cat("--------------------------------------------------------------------\n")
  cat("  3. Global Network Parameters (95% Bootstrap Confidence Intervals):\n\n")
  print(x$global.ci, row.names = FALSE)
  cat("--------------------------------------------------------------------\n")
  cat("  4. Partition Agreement Matrix (Empirical vs. Bootstrap):\n\n")
  print(x$agreement.matrix, row.names = FALSE)
  cat("====================================================================\n")
  invisible(x)
}

#' @title Summary Method for boot.rEGA Objects
#' @param object An object of class \code{"boot.rEGA"}.
#' @param ... Additional arguments (unused).
#' @export
summary.boot.rEGA <- function(object, ...) {
  print.boot.rEGA(object, ...)
  cat("  5. Item Replication Stability Across Dimensions:\n\n")
  print(object$item.stability, row.names = FALSE)
  cat("--------------------------------------------------------------------\n")
  cat("  6. Nodal Variance Triad 95% Bootstrap Confidence Intervals:\n\n")
  print(object$nodal.ci, row.names = FALSE)
  cat("====================================================================\n")
  invisible(object)
}
