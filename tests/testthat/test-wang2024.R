test_that("wang2024_aric_proteins exposes complete seqid_dot mapping", {
  prots <- wang2024_aric_proteins()

  expect_s3_class(prots, "data.frame")
  expect_equal(nrow(prots), 788)
  expect_true(all(!is.na(prots$seqid_dot) & prots$seqid_dot != ""))
})

test_that("compute_wang2024_aric_age matches all seqid_dot proteins", {
  prots <- wang2024_aric_proteins()
  demo <- data.frame(
    SampleID = paste0("S", 1:3),
    Age = c(50, 60, 70),
    stringsAsFactors = FALSE
  )

  set.seed(42)
  for (sid in prots$seqid_dot) {
    demo[[sid]] <- rlnorm(3, meanlog = log(2000), sdlog = 0.4)
  }

  result <- compute_wang2024_aric_age(demo, match_by = "seqid_dot")

  expect_equal(result$n_proteins_matched, rep(788, 3))
  expect_equal(result$n_proteins_missing, rep(0, 3))
  expect_equal(result$match_by, rep("seqid_dot", 3))
})
