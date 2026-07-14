test_that("oh2023_conventional_proteins returns complete model metadata", {
  prots <- oh2023_conventional_proteins()

  expect_s3_class(prots, "data.frame")
  expect_equal(nrow(prots), 4778)
  expect_true(all(c("SOMAID", "Gene", "UniProt", "Weight", "seqid_sl", "seqid_dot") %in% names(prots)))
})

test_that("compute_oh2023_conventional_age uses the full ensemble", {
  prots <- oh2023_conventional_proteins()
  demo <- data.frame(
    SampleID = paste0("S", 1:3),
    Age = c(60, 70, 80),
    Sex_F = c(0, 1, 0),
    stringsAsFactors = FALSE
  )

  set.seed(42)
  for (sid in prots$seqid_dot) {
    demo[[sid]] <- rlnorm(3, meanlog = log(2000), sdlog = 0.4)
  }

  result <- compute_oh2023_conventional_age(demo, match_by = "seqid_dot")

  expect_equal(result$n_proteins_matched, rep(4778, 3))
  expect_equal(result$n_proteins_missing, rep(0, 3))
  expect_equal(result$match_by, rep("seqid_dot", 3))
  expect_true(all(is.finite(result$proteomic_age)))
})
