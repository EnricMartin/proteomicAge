test_that("kuo2024_pac_proteins exposes the 128 protein predictors", {
  prots <- kuo2024_pac_proteins()

  expect_s3_class(prots, "data.frame")
  expect_equal(nrow(prots), 128)
  expect_true(all(c("predictor", "weight") %in% names(prots)))
})

test_that("compute_kuo2024_pac_age computes PAC with complete Olink data", {
  prots <- kuo2024_pac_proteins()
  demo <- data.frame(
    SampleID = paste0("S", 1:4),
    Age = c(45, 55, 65, 70),
    stringsAsFactors = FALSE
  )

  set.seed(42)
  for (protein in prots$predictor) {
    demo[[toupper(protein)]] <- stats::rnorm(4)
  }

  result <- compute_kuo2024_pac_age(demo)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4)
  expect_equal(result$n_proteins_matched, rep(128, 4))
  expect_true(all(is.finite(result$proteomic_age)))
  expect_true(all(is.finite(result$age_acceleration)))
})

test_that("goeminne2025_organaging_proteins exposes model coefficients", {
  prots <- goeminne2025_organaging_proteins()

  expect_s3_class(prots, "data.frame")
  expect_true(all(c("organ", "model_type", "fold", "protein", "weight") %in% names(prots)))
  expect_true("Conventional" %in% prots$organ)
  expect_true("chronological" %in% prots$model_type)
  expect_true("mortality" %in% prots$model_type)
})

test_that("compute_goeminne2025_organ_age computes chronological and mortality scores", {
  prots <- goeminne2025_organaging_proteins()
  needed <- unique(prots$protein[prots$organ == "Conventional" &
                                  prots$model_type == "chronological" &
                                  prots$fold == 1])
  needed <- head(needed, 25)
  demo <- data.frame(
    SampleID = paste0("S", 1:4),
    Age = c(45, 55, 65, 70),
    stringsAsFactors = FALSE
  )

  set.seed(42)
  for (protein in needed) {
    demo[[tolower(protein)]] <- stats::rnorm(4)
  }

  age_result <- compute_goeminne2025_organ_age(demo, organs = "Conventional")
  mortality_result <- compute_goeminne2025_organ_age(
    demo,
    organs = "Conventional",
    model_type = "mortality"
  )

  expect_equal(nrow(age_result), 4)
  expect_equal(age_result$organ, rep("Conventional", 4))
  expect_true(all(is.finite(age_result$proteomic_age)))
  expect_true(all(is.finite(age_result$age_acceleration)))
  expect_equal(nrow(mortality_result), 4)
  expect_true(all(is.finite(mortality_result$mortality_score)))
})

test_that("Olink gene matching accepts R-syntactic gene column names", {
  lookup <- proteomicAge:::.olink_gene_lookup(c("SampleID", "Age", "ERVV.1"), "ERVV-1")

  expect_equal(unname(lookup[["ERVV-1"]]), "ERVV.1")
})
