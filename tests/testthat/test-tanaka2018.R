test_that("tanaka2018_proteins returns correct structure", {
  prots <- tanaka2018_proteins()
  
  expect_s3_class(prots, "data.frame")
  expect_equal(names(prots), c("SOMAID", "Gene", "Weight"))
  expect_equal(nrow(prots), 76)
  
  # Check intercept exists in cache
  load_tanaka2018_coefs()  # ensure loaded
  expect_true(.tanaka2018_cache$intercept > 80)
})

test_that("compute_tanaka2018_age rejects invalid input", {
  expect_error(
    compute_tanaka2018_age("not_a_dataframe"),
    "must be a data.frame"
  )
  expect_error(
    compute_tanaka2018_age(mtcars, id_col = "nonexistent"),
    "not found in data"
  )
})

test_that("compute_tanaka2018_age works with demo data", {
  # Create synthetic data with all 76 proteins
  prots <- tanaka2018_proteins()
  n_prots <- nrow(prots)
  
  # Build a data.frame with 10 samples
  demo <- data.frame(
    SampleID = paste0("S", 1:10),
    Age = seq(22, 85, length.out = 10),
    stringsAsFactors = FALSE
  )
  
  # Add protein columns with realistic-ish RFU values
  set.seed(42)
  for (i in seq_len(n_prots)) {
    sid <- prots$SOMAID[i]
    coef <- prots$Weight[i]
    # Generate values that produce realistic ages
    # RFU values typically in hundreds to tens of thousands
    demo[[sid]] <- runif(10, 500, 15000)
  }
  
  result <- compute_tanaka2018_age(demo, id_col = "SampleID", age_col = "Age")
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 10)
  expect_equal(names(result), c(
    "id", "chronological_age", "proteomic_age",
    "age_acceleration", "n_proteins_matched", "n_proteins_missing"
  ))
  expect_equal(result$chronological_age, demo$Age)
  expect_equal(result$n_proteins_matched, rep(76, 10))
  expect_equal(result$n_proteins_missing, rep(0, 10))
  
  # Age acceleration should be computed
  expect_equal(result$age_acceleration,
               result$proteomic_age - result$chronological_age)
})

test_that("compute_tanaka2018_age warns on missing proteins", {
  prots <- tanaka2018_proteins()
  
  demo <- data.frame(
    SampleID = "S1",
    Age = 50,
    stringsAsFactors = FALSE
  )
  
  # Only add 10 of 76 proteins
  for (i in 1:10) {
    sid <- prots$SOMAID[i]
    demo[[sid]] <- 1000
  }
  
  expect_warning(
    compute_tanaka2018_age(demo, id_col = "SampleID", age_col = "Age"),
    "not found in data"
  )
})

test_that("tanaka2018_proteins returns known top proteins", {
  prots <- tanaka2018_proteins()
  
  # GDF15 is the strongest positive predictor
  gdf15 <- prots[prots$Gene == "GDF15", ]
  expect_equal(nrow(gdf15), 1)
  expect_gt(gdf15$Weight, 4.0)
  
  # CDON is strongly negative
  cdon <- prots[prots$Gene == "CDON", ]
  expect_equal(nrow(cdon), 1)
  expect_lt(cdon$Weight, -8.0)
  
  # ALB is an important positive predictor
  alb <- prots[prots$Gene == "ALB", ]
  expect_equal(nrow(alb), 1)
  expect_gt(alb$Weight, 3.0)
})

test_that("preprocess_somascan detects SOMAID columns", {
  demo <- data.frame(
    SampleID = 1:5,
    Age = c(30, 40, 50, 60, 70),
    SL003869 = c(1000, 2000, 3000, 4000, 5000),
    SL000045 = c(500, 600, 700, 800, 900),
    NotAProtein = letters[1:5]
  )
  
  expect_message(
    result <- preprocess_somascan(demo, report_missingness = FALSE),
    "Auto-detected 2 protein columns"
  )
  
  # Should preserve non-protein columns
  expect_true("SampleID" %in% names(result))
  expect_true("Age" %in% names(result))
  expect_true("NotAProtein" %in% names(result))
  
  # Protein columns should be log2-transformed
  expect_equal(result$SL003869, log2(c(1000, 2000, 3000, 4000, 5000)))
})

test_that("preprocess_somascan handles outliers", {
  set.seed(123)
  vals <- c(rnorm(98, mean = 10, sd = 1), 50, -20)  # two extreme outliers
  demo <- data.frame(SampleID = 1:100, SL003869 = vals)
  
  result <- preprocess_somascan(demo, 
    log_transform = FALSE, 
    handle_outliers = TRUE,
    report_missingness = FALSE
  )
  
  # Outliers should be winsorized
  expect_lt(max(result$SL003869, na.rm = TRUE), 50)
  expect_gt(min(result$SL003869, na.rm = TRUE), -20)
})

test_that("validate_somascan_input detects issues", {
  prots <- tanaka2018_proteins()
  
  # Good data
  good <- data.frame(SampleID = "S1", Age = 50)
  for (i in 1:nrow(prots)) {
    good[[prots$SOMAID[i]]] <- 1000
  }
  expect_output(
    res <- validate_somascan_input(good),
    "VALID"
  )
  
  # Missing age
  bad <- data.frame(SampleID = "S1")
  expect_output(
    validate_somascan_input(bad),
    "INVALID"
  )
  
  # Implausible age
  weird <- data.frame(SampleID = "S1", Age = 999)
  expect_output(
    validate_somascan_input(weird, required_proteins = character(0)),
    "plausible range"
  )
})
