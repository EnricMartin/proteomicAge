test_that("summarize_clock_qc returns one row per clock", {
  clock_outputs <- list(
    clock_a = data.frame(
      id = paste0("S", 1:5),
      chronological_age = 50:54,
      proteomic_age = c(51, 52, 52, 55, 56),
      age_acceleration = c(-1, -0.5, 0, 0.5, 1),
      n_proteins_matched = 10,
      n_proteins_missing = 2
    ),
    clock_b = data.frame(
      id = paste0("S", 1:5),
      chronological_age = 50:54,
      proteomic_age = c(49, 50, 51, 53, 54),
      age_acceleration = c(-0.2, -0.1, 0, 0.1, 0.2),
      n_proteins_matched = 20,
      n_proteins_missing = 0
    )
  )

  summary <- summarize_clock_qc(clock_outputs)

  expect_equal(summary$clock, c("clock_a", "clock_b"))
  expect_equal(summary$n, c(5, 5))
  expect_equal(summary$n_proteins_matched, c(10, 20))
  expect_true(all(is.finite(summary$agecor)))
  expect_true(all(is.finite(summary$agecorp)))
})

test_that("clock_correlation_matrix returns correlations and significance", {
  clock_outputs <- list(
    clock_a = data.frame(
      id = paste0("S", 1:5),
      chronological_age = 50:54,
      proteomic_age = 1:5,
      age_acceleration = c(-1, -0.5, 0, 0.5, 1)
    ),
    clock_b = data.frame(
      id = paste0("S", 1:5),
      chronological_age = 50:54,
      proteomic_age = 2:6,
      age_acceleration = c(-2, -1, 0, 1, 2)
    )
  )

  mat <- clock_correlation_matrix(clock_outputs)

  expect_equal(dim(mat$correlation), c(2L, 2L))
  expect_equal(mat$correlation["clock_a", "clock_b"], 1)
  expect_equal(mat$significance["clock_a", "clock_b"], "***")
})

test_that("clock summaries can consume compute_global_age output", {
  global <- data.frame(
    id = paste0("S", 1:4),
    chronological_age = 61:64,
    proteomic_age = c(60, 62, 63, 65),
    age_acceleration = c(-1, -0.2, 0.2, 1),
    n_clocks = 2,
    clocks_used = "clock_a;clock_b",
    match_by = "seqid_dot",
    clock_a_age = c(59, 61, 63, 64),
    clock_b_age = c(61, 63, 63, 66),
    clock_a_age_acceleration = c(-1, -0.4, 0.1, 0.8),
    clock_b_age_acceleration = c(-1, 0, 0.3, 1.2),
    n_proteins_matched = 30,
    n_proteins_missing = 2,
    clock_a_n_proteins_matched = 10,
    clock_a_n_proteins_missing = 2,
    clock_b_n_proteins_matched = 20,
    clock_b_n_proteins_missing = 0
  )

  summary <- summarize_clock_qc(global)

  expect_true(all(c("global", "clock_a", "clock_b") %in% summary$clock))
  expect_equal(nrow(summary), 3)
  expect_equal(summary$n_proteins_matched, c(10, 20, 30))
  expect_equal(summary$n_proteins_missing, c(2, 0, 2))
})

test_that("plot_clock_violin returns a ggplot without fixed limits", {
  skip_if_not_installed("ggplot2")
  clock_outputs <- list(
    clock_a = data.frame(
      id = paste0("S", 1:6),
      chronological_age = 50:55,
      proteomic_age = c(40, 42, 45, 60, 75, 90),
      age_acceleration = c(-5, -3, -1, 1, 3, 5)
    )
  )

  p <- plot_clock_violin(
    clock_outputs,
    group = c("A", "A", "A", "B", "B", "B")
  )

  expect_s3_class(p, "ggplot")
  expect_null(p$coordinates$limits$y)
})

test_that("compute outputs can return bundled QC and plots", {
  skip_if_not_installed("ggplot2")
  clock_outputs <- data.frame(
    id = paste0("S", 1:6),
    chronological_age = 50:55,
    proteomic_age = c(40, 42, 45, 60, 75, 90),
    age_acceleration = c(-5, -3, -1, 1, 3, 5),
    n_proteins_matched = 3,
    n_proteins_missing = 0,
    clock = "demo"
  )

  bundle <- proteomicAge:::.clock_result_bundle(
    clock_outputs,
    group = c("A", "A", "A", "B", "B", "B"),
    return_list = TRUE
  )

  expect_named(bundle, c(
    "predictions", "qc", "scatter_plot", "group_plot",
    "group_comparison", "group"
  ))
  expect_s3_class(bundle$scatter_plot, "ggplot")
  expect_s3_class(bundle$group_plot, "ggplot")
  expect_equal(bundle$qc$n_proteins_matched, 3)
  expect_equal(bundle$group_comparison$test, "Two-tailed Welch t-test")
})
