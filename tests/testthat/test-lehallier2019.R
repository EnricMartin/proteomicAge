test_that("compute_lehallier2019_age matches seqid_dot columns exactly", {
  prots <- lehallier2019_proteins()
  demo <- data.frame(
    SampleID = paste0("S", 1:4),
    Age = c(40, 50, 60, 70),
    Sex = c(0, 1, 0, 1),
    stringsAsFactors = FALSE
  )
  for (sid in paste0("seq.", prots$seqid_dot)) {
    demo[[sid]] <- 10
  }

  result <- compute_lehallier2019_age(
    demo, id_col = "SampleID", age_col = "Age", match_by = "seqid_dot"
  )

  expect_equal(result$n_proteins_matched, rep(373, 4))
  expect_equal(result$n_proteins_missing, rep(0, 4))
})
