test_that("compute_global_age combines five conventional clocks", {
  add_seq_prefix <- function(x) ifelse(grepl("^seq[.]", x), x, paste0("seq.", x))

  seqids <- unique(c(
    add_seq_prefix(tanaka2018_proteins()$seqid_dot),
    add_seq_prefix(lehallier2019_proteins()$seqid_dot),
    add_seq_prefix(sathyan2020_proteins()$seqid_dot),
    oh2023_conventional_proteins()$seqid_dot,
    wang2024_aric_proteins()$seqid_dot
  ))

  demo <- data.frame(
    SampleID = paste0("S", 1:3),
    Age = c(55, 65, 75),
    Sex_F = c(0, 1, 0),
    stringsAsFactors = FALSE
  )

  set.seed(42)
  for (sid in seqids) {
    demo[[sid]] <- rlnorm(3, meanlog = log(2000), sdlog = 0.4)
  }

  result <- suppressWarnings(compute_global_age(
    demo,
    id_col = "SampleID",
    age_col = "Age",
    sex_col = "Sex_F",
    match_by = "seqid_dot"
  ))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_equal(result$n_clocks, rep(5, 3))
  expect_true(all(is.finite(result$proteomic_age)))
  expect_true(all(is.finite(result$age_acceleration)))
  expect_true(all(c(
    "tanaka2018_age",
    "lehallier2019_age",
    "sathyan2020_age",
    "oh2023_conventional_age",
    "wang2024_aric_age"
  ) %in% names(result)))
})
