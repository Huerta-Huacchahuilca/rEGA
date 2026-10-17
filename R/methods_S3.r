# Declare NSE variables used in ggplot2::aes() to satisfy R CMD check
utils::globalVariables(c("Item", "Proportion", "Component", "Item_Stability", "Dim_Label"))

#' @title Plot Method for rEGA Objects
#' @param x An object of class \code{"rEGA"}.
#' @param type Character. Type of plot: \code{"variance"} (Nodal Variance Triad) or \code{"network"} (Purity-Weighted Graph).
#' @param ... Additional arguments passed to \code{qgraph} (if type = "network").
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_fill_manual labs theme_minimal theme element_text
#' @export
plot.rEGA <- function(x, type = c("variance", "network"), ...) {
  type <- match.arg(type)

  if (type == "variance") {
    triad <- x$nodal.triad
    df_long <- rbind(
      data.frame(Item = triad$Item, Community = triad$Community,
                 Component = "Within-Community Signal (h2)", Proportion = triad$h2),
      data.frame(Item = triad$Item, Community = triad$Community,
                 Component = "Bridging Variance (beta2)", Proportion = triad$beta2),
      data.frame(Item = triad$Item, Community = triad$Community,
                 Component = "Markov Uniqueness (psi2)", Proportion = triad$psi2)
    )
    df_long$Item <- factor(df_long$Item, levels = rev(triad$Item))
    df_long$Component <- factor(df_long$Component, levels = c(
      "Markov Uniqueness (psi2)", "Bridging Variance (beta2)", "Within-Community Signal (h2)"
    ))

    p_out <- ggplot2::ggplot(df_long, ggplot2::aes(x = Item, y = Proportion, fill = Component)) +
      ggplot2::geom_col(width = 0.72, color = "white", linewidth = 0.2) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c(
        "Within-Community Signal (h2)" = "#1B5E20",
        "Bridging Variance (beta2)" = "#F57F17",
        "Markov Uniqueness (psi2)" = "#B71C1C"
      )) +
      ggplot2::labs(
        title = sprintf("rEGA Nodal Variance Decomposition (Gamma_t = %.3f, Psi_G = %.3f)",
                        x$reliability$Gamma.t, x$reliability$Psi.G),
        subtitle = "Exact Conservation Law: h2 (Within-Community) + beta2 (Bridging) + psi2 (Markov Uniqueness) = 1",
        x = "Psychometric Node (Item)", y = "Proportion of Standardized Nodal Variance", fill = "Variance Component"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold"))

    return(p_out)
  }

  if (type == "network") {
    if (!requireNamespace("qgraph", quietly = TRUE)) stop("The 'qgraph' package is required to plot the network.")

    cat("Plotting rEGA Purity-Weighted Partial Correlation Network...\n")
    p_out <- qgraph::qgraph(
      x$network.refit,
      groups = as.factor(x$wc),
      layout = "spring",
      theme = "colorblind",
      vsize = 6,
      node.width = 1.5,
      label.scale = FALSE,
      label.cex = 1.1,
      title = "rEGA Support-Constrained Network",
      DoNotPlot = FALSE,
      ...
    )
    return(invisible(p_out))
  }
}

#' @title Plot Method for boot.rEGA Objects
#' @param x An object of class \code{"boot.rEGA"}.
#' @param ... Additional arguments (unused).
#' @importFrom ggplot2 ggplot aes geom_point geom_segment geom_hline coord_flip labs theme_minimal ylim
#' @export
plot.boot.rEGA <- function(x, ...) {
  df <- x$item.stability
  df$Item <- factor(df$Item, levels = rev(df$Item))
  df$Dim_Label <- paste0("Dim ", df$Empirical_Dim)

  p_stab <- ggplot2::ggplot(df, ggplot2::aes(x = Item, y = Item_Stability, color = Dim_Label)) +
    ggplot2::geom_segment(ggplot2::aes(x = Item, xend = Item, y = 0, yend = Item_Stability), linewidth = 1.1) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::geom_hline(yintercept = 0.70, linetype = "dashed", color = "firebrick") +
    ggplot2::coord_flip() +
    ggplot2::ylim(0, 1.05) +
    ggplot2::labs(
      title = "boot.rEGA Item Replicability Across Bootstrap Replications",
      subtitle = "Dashed line indicates recommended structural stability threshold (0.70)",
      x = "Item", y = "Replication Proportion in Assigned Dimension", color = "Empirical Dimension"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  return(p_stab)
}
