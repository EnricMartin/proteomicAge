test_that("compute_sathyan2020_age preserves duplicate coefficient contributions", {
  prots <- sathyan2020_proteins()
  dup_ids <- names(which(table(prots$SOMAID) > 1))
  expect_true("SL012540" %in% dup_ids)

  dup <- prots[prots$SOMAID == "SL012540", ]
  demo <- data.frame(
    SampleID = paste0("S", 1:4),
    Age = c(55, 60, 65, 70),
    stringsAsFactors = FALSE
  )
  demo[[dup$seqid_sl[[1]]]] <- c(1, 2, 3, 4)

  result <- suppressWarnings(compute_sathyan2020_age(
    demo, id_col = "SampleID", age_col = "Age", match_by = "seqid_sl"
  ))

  proteomicAge:::load_sathyan2020_coefs()
  expected <- proteomicAge:::.sathyan2020_cache$intercept +
    sum(dup$Weight) * demo[[dup$seqid_sl[[1]]]]
  expect_equal(result$proteomic_age, expected)
})
