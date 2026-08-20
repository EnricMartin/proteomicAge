#' Convert clock outputs to long format
#'
#' @keywords internal
.clock_outputs_long <- function(clock_outputs) {
  if (is.data.frame(clock_outputs)) {
    required <- c("id", "chronological_age", "proteomic_age", "age_acceleration")
    if (!all(required %in% names(clock_outputs))) {
      stop("clock_outputs must contain id, chronological_age, proteomic_age, and age_acceleration")
    }

    rows <- list(
      global = data.frame(
        id = clock_outputs$id,
        chronological_age = clock_outputs$chronological_age,
        clock = "global",
        proteomic_age = clock_outputs$proteomic_age,
        age_acceleration = clock_outputs$age_acceleration,
        n_proteins_matched = NA_integer_,
        n_proteins_missing = NA_integer_,
        stringsAsFactors = FALSE
      )
    )

    component_age_cols <- grep("_age$", names(clock_outputs), value = TRUE)
    component_age_cols <- setdiff(component_age_cols, c("proteomic_age", "chronological_age"))
    for (age_col in component_age_cols) {
      clock <- sub("_age$", "", age_col)
      gap_col <- paste0(clock, "_age_acceleration")
      if (!gap_col %in% names(clock_outputs)) next
      rows[[clock]] <- data.frame(
        id = clock_outputs$id,
        chronological_age = clock_outputs$chronological_age,
        clock = clock,
        proteomic_age = clock_outputs[[age_col]],
        age_acceleration = clock_outputs[[gap_col]],
        n_proteins_matched = NA_integer_,
        n_proteins_missing = NA_integer_,
        stringsAsFactors = FALSE
      )
    }
    return(do.call(rbind, rows))
  }

  if (!is.list(clock_outputs) || is.null(names(clock_outputs))) {
    stop("clock_outputs must be a named list of clock result data.frames or a Global Age result")
  }

  rows <- lapply(names(clock_outputs), function(clock) {
    result <- clock_outputs[[clock]]
    required <- c("id", "chronological_age", "proteomic_age", "age_acceleration")
    if (!is.data.frame(result) || !all(required %in% names(result))) {
      stop("Each clock result must contain id, chronological_age, proteomic_age, and age_acceleration")
    }
    matched <- if ("n_proteins_matched" %in% names(result)) result$n_proteins_matched else NA_integer_
    missing <- if ("n_proteins_missing" %in% names(result)) result$n_proteins_missing else NA_integer_
    data.frame(
      id = result$id,
      chronological_age = result$chronological_age,
      clock = clock,
      proteomic_age = result$proteomic_age,
      age_acceleration = result$age_acceleration,
      n_proteins_matched = matched,
      n_proteins_missing = missing,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

.clock_wide <- function(clock_outputs, value) {
  long <- .clock_outputs_long(clock_outputs)
  ids <- unique(long$id)
  wide <- data.frame(id = ids, stringsAsFactors = FALSE)
  for (clock in unique(long$clock)) {
    sub <- long[long$clock == clock, c("id", value)]
    names(sub)[2] <- clock
    wide <- merge(wide, sub, by = "id", all.x = TRUE, sort = FALSE)
  }
  wide
}

.cor_p_value <- function(x, y, method) {
  ok <- stats::complete.cases(x, y)
  if (sum(ok) < 3 || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
    return(NA_real_)
  }
  stats::cor.test(x[ok], y[ok], method = method)$p.value
}

#' Summarize QC metrics for proteomic age clock outputs
#'
#' @param clock_outputs A named list of outputs from `compute_*_age()` functions,
#'   or the output of `compute_global_age()`.
#' @return A data.frame with one row per clock.
#' @export
summarize_clock_qc <- function(clock_outputs) {
  long <- .clock_outputs_long(clock_outputs)
  rows <- lapply(split(long, long$clock), function(x) {
    ok <- stats::complete.cases(x$chronological_age, x$proteomic_age)
    has_variation <- sum(ok) >= 3 &&
      stats::sd(x$chronological_age[ok]) > 0 &&
      stats::sd(x$proteomic_age[ok]) > 0
    agecor <- if (has_variation) stats::cor(x$chronological_age[ok], x$proteomic_age[ok]) else NA_real_
    agecorp <- if (has_variation) stats::cor.test(x$chronological_age[ok], x$proteomic_age[ok])$p.value else NA_real_
    data.frame(
      clock = x$clock[1],
      n = sum(ok),
      n_proteins_matched = suppressWarnings(max(x$n_proteins_matched, na.rm = TRUE)),
      n_proteins_missing = suppressWarnings(max(x$n_proteins_missing, na.rm = TRUE)),
      agecor = agecor,
      agecorp = agecorp,
      mean_chronological_age = mean(x$chronological_age, na.rm = TRUE),
      mean_omic_age = mean(x$proteomic_age, na.rm = TRUE),
      mean_age_acceleration = mean(x$age_acceleration, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$n_proteins_matched[is.infinite(out$n_proteins_matched)] <- NA
  out$n_proteins_missing[is.infinite(out$n_proteins_missing)] <- NA
  rownames(out) <- NULL
  out
}

#' Compute a correlation matrix between proteomic clocks
#'
#' @param clock_outputs A named list of outputs from `compute_*_age()` functions,
#'   or the output of `compute_global_age()`.
#' @param value `"proteomic_age"` or `"age_acceleration"`.
#' @param method Correlation method passed to `cor.test()`.
#' @param p_adjust_method Method passed to `p.adjust()`.
#' @return A list containing correlation, p-value, adjusted p-value, and
#'   significance matrices.
#' @export
clock_correlation_matrix <- function(clock_outputs,
                                     value = c("proteomic_age", "age_acceleration"),
                                     method = c("pearson", "spearman"),
                                     p_adjust_method = "none") {
  value <- match.arg(value)
  method <- match.arg(method)
  wide <- .clock_wide(clock_outputs, value)
  mat <- as.matrix(wide[, setdiff(names(wide), "id"), drop = FALSE])
  cor_mat <- stats::cor(mat, use = "pairwise.complete.obs", method = method)
  p_mat <- matrix(NA_real_, nrow = ncol(mat), ncol = ncol(mat), dimnames = dimnames(cor_mat))
  for (i in seq_len(ncol(mat))) {
    for (j in seq_len(ncol(mat))) {
      p_mat[i, j] <- if (i == j) 0 else .cor_p_value(mat[, i], mat[, j], method)
    }
  }
  padj <- matrix(stats::p.adjust(as.vector(p_mat), method = p_adjust_method),
                 nrow = nrow(p_mat), dimnames = dimnames(p_mat))
  stars <- matrix("", nrow = nrow(padj), ncol = ncol(padj), dimnames = dimnames(padj))
  stars[padj < 0.05] <- "*"
  stars[padj < 0.01] <- "**"
  stars[padj < 0.001] <- "***"
  list(correlation = cor_mat, p_value = p_mat, p_adjusted = padj,
       significance = stars, value = value, method = method)
}

#' Plot a clock correlation matrix
#'
#' @param clock_outputs A named list of clock outputs or a Global Age result.
#' @param value `"proteomic_age"` or `"age_acceleration"`.
#' @param method Correlation method.
#' @param p_adjust_method Method passed to `p.adjust()`.
#' @return Invisibly returns the result from `clock_correlation_matrix()`.
#' @export
plot_clock_correlation_matrix <- function(clock_outputs,
                                          value = c("proteomic_age", "age_acceleration"),
                                          method = c("pearson", "spearman"),
                                          p_adjust_method = "none") {
  result <- clock_correlation_matrix(clock_outputs, value, method, p_adjust_method)
  cor_mat <- result$correlation
  stars <- result$significance
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(7, 7, 3, 1))
  graphics::image(seq_len(nrow(cor_mat)), seq_len(ncol(cor_mat)), cor_mat[nrow(cor_mat):1, ],
                  axes = FALSE, zlim = c(-1, 1), col = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
                  xlab = "", ylab = "", main = paste(result$method, result$value, "correlation"))
  labels <- colnames(cor_mat)
  graphics::axis(1, at = seq_along(labels), labels = labels, las = 2)
  graphics::axis(2, at = seq_along(labels), labels = rev(labels), las = 2)
  for (i in seq_len(nrow(cor_mat))) {
    for (j in seq_len(ncol(cor_mat))) {
      graphics::text(i, ncol(cor_mat) - j + 1, paste0(sprintf("%.2f", cor_mat[j, i]), stars[j, i]))
    }
  }
  invisible(result)
}

#' Plot proteomic age or age acceleration against chronological age
#'
#' @param clock_outputs A named list of clock outputs or a Global Age result.
#' @param value `"proteomic_age"` or `"age_acceleration"`.
#' @return Invisibly returns the long-format data used for plotting.
#' @export
plot_clock_scatter <- function(clock_outputs,
                               value = c("proteomic_age", "age_acceleration")) {
  value <- match.arg(value)
  long <- .clock_outputs_long(clock_outputs)
  clocks <- unique(long$clock)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = grDevices::n2mfrow(length(clocks)))
  for (clock in clocks) {
    x <- long[long$clock == clock, ]
    graphics::plot(x$chronological_age, x[[value]], pch = 19, col = "#2F6F9F",
                   xlab = "Chronological age", ylab = value, main = clock)
    if (value == "proteomic_age") graphics::abline(0, 1, col = "gray60", lty = 2)
    graphics::abline(stats::lm(x[[value]] ~ x$chronological_age), col = "#B23A48", lwd = 2)
  }
  invisible(long)
}

#' Plot clock outputs by a categorical variable
#'
#' @param clock_outputs A named list of clock outputs or a Global Age result.
#' @param group A grouping vector, or a column name in `sample_data`.
#' @param sample_data Optional data.frame containing the grouping column and IDs.
#' @param value `"proteomic_age"` or `"age_acceleration"`.
#' @return Invisibly returns the long-format data used for plotting.
#' @export
plot_clock_violin <- function(clock_outputs,
                              group,
                              sample_data = NULL,
                              value = c("proteomic_age", "age_acceleration")) {
  value <- match.arg(value)
  long <- .clock_outputs_long(clock_outputs)
  group_values <- if (length(group) == 1 && is.character(group) && is.data.frame(sample_data)) {
    sample_data[[group]][match(long$id, sample_data$id)]
  } else {
    rep(group, length.out = length(unique(long$id)))[match(long$id, unique(long$id))]
  }
  long$group <- as.factor(group_values)
  clocks <- unique(long$clock)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = grDevices::n2mfrow(length(clocks)))
  for (clock in clocks) {
    x <- long[long$clock == clock, ]
    groups <- levels(x$group)
    y_range <- range(x[[value]], na.rm = TRUE)
    graphics::plot(seq_along(groups), rep(NA_real_, length(groups)), ylim = y_range,
                   xaxt = "n", xlab = "", ylab = value, main = clock)
    graphics::axis(1, at = seq_along(groups), labels = groups)
    for (i in seq_along(groups)) {
      vals <- x[[value]][x$group == groups[i]]
      vals <- vals[is.finite(vals)]
      if (length(vals) < 2 || stats::sd(vals) == 0) {
        graphics::points(i, vals, pch = 19, col = "#2F6F9F")
        next
      }
      dens <- stats::density(vals)
      width <- 0.35 * dens$y / max(dens$y)
      graphics::polygon(c(i - width, rev(i + width)), c(dens$x, rev(dens$x)),
                        col = "#9ECAE1", border = "#2F6F9F")
      graphics::points(rep(i, length(vals)), vals, pch = 16, cex = 0.55, col = "#33333366")
    }
  }
  invisible(long)
}
