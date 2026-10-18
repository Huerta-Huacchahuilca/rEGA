# Declare NSE variables used in ggplot2::aes() to satisfy R CMD check
utils::globalVariables(c(
  "Item", "Proportion", "Component", "Item_Stability", "Dim_Label",
  "Pct_Label", "Triad_Summary", "x", "y", "xend", "yend", "Weight",
  "Sign", "Node", "Stability_Label"
))

.rEGA_palette <- function(n_comms) {
  base_cols <- c(
    "#0072B2", "#E69F00", "#009E73", "#CC79A7",
    "#56B4E9", "#D55E00", "#F0E442", "#332288"
  )
  if (n_comms <= length(base_cols)) {
    return(base_cols[seq_len(n_comms)])
  }
  grDevices::hcl.colors(n_comms, palette = "Dynamic")
}

.build_network_ggplot <- function(rEGA.obj) {
  Omega <- rEGA.obj$network.refit
  wc <- rEGA.obj$wc
  p <- ncol(Omega)
  item_names <- colnames(Omega)
  seed_val <- if (!is.null(rEGA.obj$audit.metadata$seed)) rEGA.obj$audit.metadata$seed else 2026L

  A_abs <- abs(Omega)
  diag(A_abs) <- 0
  g <- igraph::graph_from_adjacency_matrix(A_abs, mode = "undirected", weighted = TRUE, diag = FALSE)
  set.seed(seed_val)
  if (sum(A_abs) > 1e-6) {
    coords <- igraph::layout_with_fr(g, weights = igraph::E(g)$weight)
  } else {
    angles <- seq(0, 2 * pi, length.out = p + 1)[seq_len(p)]
    coords <- cbind(cos(angles), sin(angles))
  }

  for (col_i in 1:2) {
    rng <- range(coords[, col_i])
    if (diff(rng) > 1e-6) {
      coords[, col_i] <- (coords[, col_i] - rng[1]) / diff(rng) * 2 - 1
    } else {
      coords[, col_i] <- 0
    }
  }

  nodes_df <- data.frame(
    Node = item_names,
    x = coords[, 1],
    y = coords[, 2],
    Dim_Label = factor(paste0("Dimension ", wc),
                       levels = paste0("Dimension ", sort(unique(wc)))),
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

  n_comms <- length(unique(wc))
  pal <- stats::setNames(.rEGA_palette(n_comms), paste0("Dimension ", sort(unique(wc))))

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
      ggplot2::scale_alpha_continuous(range = c(0.25, 0.85)) +
      ggplot2::scale_color_manual(values = c("Positive" = "#1B5E20", "Negative" = "#B71C1C")) +
      ggplot2::scale_linetype_manual(values = c("Positive" = "solid", "Negative" = "dashed"))
  }

  p_net <- p_net +
    ggplot2::geom_point(
      data = nodes_df,
      ggplot2::aes(x = x, y = y, fill = Dim_Label),
      shape = 21, size = 11.5, color = "grey15", stroke = 0.9
    ) +
    ggplot2::geom_text(
      data = nodes_df,
      ggplot2::aes(x = x, y = y, label = Node),
      color = "white", fontface = "bold", size = 3.3
    ) +
    ggplot2::scale_fill_manual(values = pal, name = "Community") +
    ggplot2::coord_equal(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
    ggplot2::labs(
      title = sprintf("rEGA Network (%d %s)", rEGA.obj$n.dim,
                      ifelse(rEGA.obj$n.dim == 1, "Dimension", "Dimensions")),
      subtitle = sprintf("Support-Refitted GGM | Gamma_t = %.3f | Psi_G = %.3f",
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

.build_variance_ggplot <- function(rEGA.obj, show.values = TRUE) {
  triad <- rEGA.obj$nodal.triad
  p <- nrow(triad)

  # Ensure exact 100.0% conservation per item for percentage labels
  tot_i <- triad$h2 + triad$beta2 + triad$psi2
  pct_h2 <- round((triad$h2 / tot_i) * 100, 1)
  pct_b2 <- round((triad$beta2 / tot_i) * 100, 1)
  pct_p2 <- round(100.0 - pct_h2 - pct_b2, 1)

  df_long <- rbind(
    data.frame(
      Item = triad$Item, Community = triad$Community,
      Component = "Within-Community Signal (h2)",
      Proportion = triad$h2 / tot_i,
      Pct = pct_h2,
      Pct_Label = ifelse(pct_h2 >= 5.5, sprintf("%.1f%%", pct_h2), ""),
      stringsAsFactors = FALSE
    ),
    data.frame(
      Item = triad$Item, Community = triad$Community,
      Component = "Bridging Variance (beta2)",
      Proportion = triad$beta2 / tot_i,
      Pct = pct_b2,
      Pct_Label = ifelse(pct_b2 >= 5.5, sprintf("%.1f%%", pct_b2), ""),
      stringsAsFactors = FALSE
    ),
    data.frame(
      Item = triad$Item, Community = triad$Community,
      Component = "Markov Uniqueness (psi2)",
      Proportion = triad$psi2 / tot_i,
      Pct = pct_p2,
      Pct_Label = ifelse(pct_p2 >= 5.5, sprintf("%.1f%%", pct_p2), ""),
      stringsAsFactors = FALSE
    )
  )

  # If show.values = TRUE, append exact 3-component percentage readout to item label
  if (isTRUE(show.values)) {
    item_labels <- sprintf(
      "%s  [h2: %4.1f%% | b2: %4.1f%% | psi2: %4.1f%%]",
      triad$Item, pct_h2, pct_b2, pct_p2
    )
    names(item_labels) <- triad$Item
    df_long$Item <- factor(df_long$Item, levels = rev(triad$Item), labels = rev(item_labels))
  } else {
    df_long$Item <- factor(df_long$Item, levels = rev(triad$Item))
  }

  df_long$Component <- factor(df_long$Component, levels = c(
    "Markov Uniqueness (psi2)",
    "Bridging Variance (beta2)",
    "Within-Community Signal (h2)"
  ))

  p_var <- ggplot2::ggplot(df_long, ggplot2::aes(x = Item, y = Proportion, fill = Component)) +
    ggplot2::geom_col(width = 0.72, color = "white", linewidth = 0.25) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      labels = function(v) paste0(round(v * 100), "%"),
      expand = c(0.01, 0.01)
    ) +
    ggplot2::scale_fill_manual(values = c(
      "Within-Community Signal (h2)" = "#1B5E20",
      "Bridging Variance (beta2)" = "#F57F17",
      "Markov Uniqueness (psi2)" = "#B71C1C"
    )) +
    ggplot2::labs(
      title = sprintf("Nodal Variance Triad (Gamma_t = %.3f, Psi_G = %.3f)",
                      rEGA.obj$reliability$Gamma.t, rEGA.obj$reliability$Psi.G),
      subtitle = "Exact Conservation per Item: h2% + beta2% + psi2% = 100.0% (No Rigid Cutoffs)",
      x = "Item [Exact Triad Percentages: h2 | beta2 | psi2]",
      y = "Percentage of Standardized Nodal Variance (100%)",
      fill = "Variance Component"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(family = "mono", face = "bold", size = 9.5)
    )

  if (isTRUE(show.values)) {
    p_var <- p_var +
      ggplot2::geom_text(
        ggplot2::aes(label = Pct_Label),
        position = ggplot2::position_stack(vjust = 0.5),
        color = "white", fontface = "bold", size = 3.1
      )
  }

  return(p_var)
}

.build_stability_ggplot <- function(boot.obj, show.values = TRUE) {
  df <- boot.obj$item.stability
  df$Item <- factor(df$Item, levels = rev(df$Item))
  unique_dims <- sort(unique(df$Empirical_Dim))
  df$Dim_Label <- factor(paste0("Dimension ", df$Empirical_Dim),
                         levels = paste0("Dimension ", unique_dims))
  df$Stability_Label <- sprintf("%.1f%%", df$Item_Stability * 100)

  pal <- stats::setNames(.rEGA_palette(length(unique_dims)),
                         paste0("Dimension ", unique_dims))

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
    ggplot2::scale_color_manual(values = pal, name = "Empirical Community") +
    ggplot2::labs(
      title = sprintf("Item Replicability (%d Bootstrap Replications)", boot.obj$converged.iter),
      subtitle = "Replication proportion in assigned community (Holistic evaluation, no rigid cutoffs)",
      x = "Psychometric Node (Item)",
      y = "Bootstrap Replication Stability (%)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold")
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
#' Visualizes the EGA-style Purity-Weighted Partial Correlation Network, the Nodal Variance
#' Conservation Triad (with standardized item percentages summing to 100 percent), or both side-by-side.
#' @param x An object of class \code{"rEGA"}.
#' @param type Character. Type of plot: \code{"both"} (default, side-by-side EGA network and
#'   Nodal Variance Triad), \code{"network"} (EGA-style network), or \code{"variance"}
#'   (Nodal Variance Triad bar chart).
#' @param show.values Logical. If \code{TRUE} (default), displays the exact standardized
#'   percentages for each of the 3 nodal variance components (\code{h2 pct + beta2 pct + psi2 pct = 100 percent})
#'   on every item.
#' @param ... Additional arguments (unused).
#' @importFrom ggplot2 ggplot aes geom_col geom_point geom_segment geom_text coord_flip coord_equal scale_fill_manual scale_color_manual scale_linewidth_continuous scale_alpha_continuous scale_linetype_manual scale_y_continuous labs theme_minimal theme element_text element_blank position_stack
#' @importFrom grid grid.newpage pushViewport viewport grid.layout
#' @importFrom grDevices hcl.colors
#' @export
plot.rEGA <- function(x, type = c("both", "variance", "network"), show.values = TRUE, ...) {
  type <- match.arg(type)

  if (type == "network") {
    return(.build_network_ggplot(x))
  }
  if (type == "variance") {
    return(.build_variance_ggplot(x, show.values = show.values))
  }

  p_net <- .build_network_ggplot(x)
  p_var <- .build_variance_ggplot(x, show.values = show.values)
  return(.draw_side_by_side(list(Network = p_net, Variance_Triad = p_var)))
}

#' @title Plot Method for boot.rEGA Objects
#' @description
#' Generates a unified side-by-side figure combining the empirical \code{rEGA} network
#' (EGA-style) and the bootstrap item stability plot (similar to \code{bootEGA}), or
#' individual diagnostic panels, without rigid cutoff lines.
#' @param x An object of class \code{"boot.rEGA"}.
#' @param type Character. Visualization layout: \code{"both"} (default, side-by-side
#'   \code{rEGA} Network + Item Stability), \code{"all"} (3-panel: Network + Nodal Triad +
#'   Item Stability), \code{"stability"} (Item Stability only), \code{"network"} (Network only),
#'   or \code{"variance"} (Nodal Variance Triad only).
#' @param show.values Logical. Display exact percentages on items? Default is \code{TRUE}.
#' @param ... Additional arguments (unused).
#' @export
plot.boot.rEGA <- function(x, type = c("both", "stability", "network", "variance", "all"),
                           show.values = TRUE, ...) {
  type <- match.arg(type)
  p_stab <- .build_stability_ggplot(x, show.values = show.values)

  if (type == "stability" || is.null(x$rEGA.obj)) {
    return(p_stab)
  }
  if (type == "network") {
    return(.build_network_ggplot(x$rEGA.obj))
  }
  if (type == "variance") {
    return(.build_variance_ggplot(x$rEGA.obj, show.values = show.values))
  }
  if (type == "all") {
    p_net <- .build_network_ggplot(x$rEGA.obj)
    p_var <- .build_variance_ggplot(x$rEGA.obj, show.values = show.values)
    return(.draw_side_by_side(list(Network = p_net, Variance_Triad = p_var, Stability = p_stab)))
  }

  # Default "both": Unified rEGA Network + Item Stability (bootEGA style)
  p_net <- .build_network_ggplot(x$rEGA.obj)
  return(.draw_side_by_side(list(Network = p_net, Stability = p_stab)))
}

#' @title Print Method for rEGA Objects
#' @param x An object of class \code{"rEGA"}.
#' @param ... Additional arguments (unused).
#' @export
print.rEGA <- function(x, ...) {
  cat("======================================================================\n")
  cat(" Reliability Exploratory Graph Analysis (rEGA)\n")
  cat("======================================================================\n")
  cat(sprintf(" Pathway: %s | Sample Size (n): %d | Items (p): %d\n",
              x$audit.metadata$pathway, x$audit.metadata$n, x$audit.metadata$p))
  cat(sprintf(" Correlation: %s | Community Algorithm: %s | Active Edges: %d\n",
              x$audit.metadata$cor.method, x$audit.metadata$algorithm,
              x$audit.metadata$active.edges))
  cat(sprintf(" Estimated Dimensions: %d | Total Gamma_t: %.4f | Global Noise Psi_G: %.4f\n",
              x$n.dim, x$reliability$Gamma.t, x$reliability$Psi.G))
  cat(sprintf(" Standardized SEM: %.4f | Raw-Score SEM: %.4f\n",
              x$reliability$SEM.total.std, x$reliability$SEM.total.raw))
  cat("----------------------------------------------------------------------\n")
  cat("1. Dimensional Structural Reliability (Coefficient Gamma_c & Graph SEM):\n")
  print(x$reliability$dimensions, row.names = FALSE)
  cat("----------------------------------------------------------------------\n")
  cat("2. Nodal Variance Conservation Triad (h2 + beta2 + psi2 = 1 | 100%):\n")
  print(x$nodal.triad[, c("Item", "Community", "h2", "beta2", "psi2",
                          "h2_pct", "beta2_pct", "psi2_pct")], row.names = FALSE)
  cat("======================================================================\n")
  invisible(x)
}

#' @title Summary Method for rEGA Objects
#' @param object An object of class \code{"rEGA"}.
#' @param ... Additional arguments (unused).
#' @export
summary.rEGA <- function(object, ...) {
  print.rEGA(object, ...)
  cat("3. Standardized Network Loadings:\n")
  print(object$net.loadings)
  cat("======================================================================\n")
  invisible(object)
}

#' @title Print Method for boot.rEGA Objects
#' @param x An object of class \code{"boot.rEGA"}.
#' @param ... Additional arguments (unused).
#' @export
print.boot.rEGA <- function(x, ...) {
  cat("======================================================================\n")
  cat(" Parallel Bootstrap Stability & Inference for rEGA (boot.rEGA)\n")
  cat("======================================================================\n")
  cat(sprintf(" Resampling Type: %s | Converged Replications: %d / %d (%.2f sec)\n",
              x$type, x$converged.iter, x$iter, x$runtime.sec))
  cat("----------------------------------------------------------------------\n")
  cat("1. Dimension Frequency Distribution:\n")
  print(x$dimension.frequency, row.names = FALSE)
  cat("----------------------------------------------------------------------\n")
  cat("2. Structural Consistency & Mean Item Stability per Dimension:\n")
  print(x$structural.consistency, row.names = FALSE)
  cat("----------------------------------------------------------------------\n")
  cat("3. Global Reliability & Graph SEM 95% Bootstrap Confidence Intervals:\n")
  print(x$global.ci, row.names = FALSE)
  cat("----------------------------------------------------------------------\n")
  cat("4. Diagnostic Agreement Matrix (Holistic evaluation, no rigid cutoffs):\n")
  print(x$agreement.matrix, row.names = FALSE)
  cat("======================================================================\n")
  invisible(x)
}

#' @title Summary Method for boot.rEGA Objects
#' @param object An object of class \code{"boot.rEGA"}.
#' @param ... Additional arguments (unused).
#' @export
summary.boot.rEGA <- function(object, ...) {
  print.boot.rEGA(object, ...)
  cat("5. Nodal Variance Triad 95% Bootstrap CIs & Item Stability:\n")
  print(object$nodal.ci, row.names = FALSE)
  cat("======================================================================\n")
  invisible(object)
}
