# Optional QC and visualization bundle for compute_*_age() outputs.

.clock_resolve_group <- function(predictions, group, sample_data, source_data, id_col) {
  if (is.null(group)) return(NULL)

  if (length(group) == 1 && is.character(group)) {
    if (is.data.frame(sample_data) && group %in% names(sample_data)) {
      if ("id" %in% names(sample_data)) {
        return(sample_data[[group]][match(predictions$id, sample_data$id)])
      }
      if (id_col %in% names(sample_data)) {
        return(sample_data[[group]][match(predictions$id, sample_data[[id_col]])])
      }
      return(rep(sample_data[[group]], length.out = nrow(predictions)))
    }
    if (is.data.frame(source_data) && group %in% names(source_data)) {
      return(source_data[[group]][match(predictions$id, source_data[[id_col]])])
    }
  }

  group[match(predictions$id, unique(predictions$id))]
}

.clock_bundle_input <- function(predictions, clock_name) {
  if (is.null(clock_name)) predictions else stats::setNames(list(predictions), clock_name)
}

.clock_scatter_plot <- function(predictions, value = "proteomic_age", clock_name = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  long <- .clock_outputs_long(.clock_bundle_input(predictions, clock_name))
  long$value_plot <- long[[value]]
  ggplot2::ggplot(long, ggplot2::aes(x = chronological_age, y = value_plot)) +
    ggplot2::geom_point(color = "#2F6F9F", alpha = 0.75, na.rm = TRUE) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "#B23A48", na.rm = TRUE) +
    ggplot2::facet_wrap(stats::as.formula("~ clock"), scales = "free_y") +
    ggplot2::labs(x = "Chronological age", y = value) +
    ggplot2::theme_bw()
}

.clock_group_plot <- function(predictions, group_values, value = "age_acceleration", clock_name = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || is.null(group_values)) return(NULL)
  long <- .clock_outputs_long(.clock_bundle_input(predictions, clock_name))
  long$group <- as.factor(group_values)
  long$value_plot <- long[[value]]
  ggplot2::ggplot(long, ggplot2::aes(x = group, y = value_plot)) +
    ggplot2::geom_violin(trim = FALSE, fill = "#9ECAE1", color = "#2F6F9F", na.rm = TRUE) +
    ggplot2::geom_jitter(width = 0.12, height = 0, size = 0.9, alpha = 0.35, na.rm = TRUE) +
    ggplot2::facet_wrap(stats::as.formula("~ clock"), scales = "free_y") +
    ggplot2::labs(x = NULL, y = value) +
    ggplot2::theme_bw()
}

.clock_group_comparison <- function(predictions, group_values, value = "age_acceleration", clock_name = NULL) {
  if (is.null(group_values)) return(NULL)
  long <- .clock_outputs_long(.clock_bundle_input(predictions, clock_name))
  long$group <- as.factor(group_values)
  rows <- lapply(split(long, long$clock), function(x) {
    groups <- stats::na.omit(unique(x$group))
    if (length(groups) != 2) {
      return(data.frame(
        clock = x$clock[1], test = NA_character_, p_value = NA_real_,
        group_1 = NA_character_, group_2 = NA_character_,
        mean_1 = NA_real_, mean_2 = NA_real_, difference_2_minus_1 = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    v1 <- x[[value]][x$group == groups[1]]
    v2 <- x[[value]][x$group == groups[2]]
    ok <- stats::complete.cases(v1)
    v1 <- v1[ok]
    ok <- stats::complete.cases(v2)
    v2 <- v2[ok]
    p <- if (length(v1) >= 2 && length(v2) >= 2) stats::t.test(v2, v1)$p.value else NA_real_
    data.frame(
      clock = x$clock[1],
      test = "Two-tailed Welch t-test",
      p_value = p,
      group_1 = as.character(groups[1]),
      group_2 = as.character(groups[2]),
      mean_1 = mean(v1, na.rm = TRUE),
      mean_2 = mean(v2, na.rm = TRUE),
      difference_2_minus_1 = mean(v2, na.rm = TRUE) - mean(v1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.clock_result_bundle <- function(predictions,
                                 group = NULL,
                                 sample_data = NULL,
                                 source_data = NULL,
                                 id_col = "SampleID",
                                 return_list = FALSE,
                                 clock_name = "clock") {
  group_values <- .clock_resolve_group(predictions, group, sample_data, source_data, id_col)
  if (!return_list && is.null(group_values)) return(predictions)
  qc_input <- .clock_bundle_input(predictions, clock_name)

  list(
    predictions = predictions,
    qc = summarize_clock_qc(qc_input),
    scatter_plot = .clock_scatter_plot(predictions, clock_name = clock_name),
    group_plot = .clock_group_plot(predictions, group_values, clock_name = clock_name),
    group_comparison = .clock_group_comparison(predictions, group_values, clock_name = clock_name),
    group = group_values
  )
}
