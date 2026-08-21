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

test_that("compute_kuo2024_pac_age accepts UniProt columns with protein map", {
  prots <- kuo2024_pac_proteins()
  demo <- data.frame(
    SampleID = paste0("S", 1:4),
    Age = c(45, 55, 65, 70),
    stringsAsFactors = FALSE
  )
  protein_map <- data.frame(
    gene = prots$predictor,
    uniprot = paste0("UP", seq_len(nrow(prots))),
    stringsAsFactors = FALSE
  )

  set.seed(42)
  for (uniprot in protein_map$uniprot) {
    demo[[uniprot]] <- stats::rnorm(4)
  }

  result <- compute_kuo2024_pac_age(
    demo,
    match_by = "uniprot",
    protein_map = protein_map
  )

  expect_equal(result$n_proteins_matched, rep(128, 4))
  expect_equal(result$match_by, rep("uniprot", 4))
  expect_true(all(is.finite(result$proteomic_age)))
})

test_that("compute_kuo2024_pac_age uses built-in Olink UniProt map", {
  protein_map <- olink_protein_map()
  prots <- kuo2024_pac_proteins()
  map <- protein_map[toupper(protein_map$gene) %in% toupper(prots$predictor), ]
  map <- map[!duplicated(toupper(map$gene)), ]
  demo <- data.frame(
    SampleID = paste0("S", 1:4),
    Age = c(45, 55, 65, 70),
    stringsAsFactors = FALSE
  )

  set.seed(42)
  for (uniprot in map$uniprot) {
    demo[[uniprot]] <- stats::rnorm(4)
  }

  result <- compute_kuo2024_pac_age(demo, match_by = "uniprot")

  expect_equal(result$n_proteins_matched, rep(128, 4))
  expect_equal(result$match_by, rep("uniprot", 4))
  expect_true(all(is.finite(result$proteomic_age)))
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

test_that("compute_goeminne2025_organ_age accepts UniProt columns with protein map", {
  prots <- goeminne2025_organaging_proteins()
  needed <- unique(prots$protein[prots$organ == "Conventional" &
                                  prots$model_type == "chronological" &
                                  prots$fold == 1])
  needed <- setdiff(head(needed, 25), "Intercept")
  protein_map <- data.frame(
    gene = needed,
    uniprot = paste0("UP", seq_along(needed)),
    stringsAsFactors = FALSE
  )
  demo <- data.frame(
    SampleID = paste0("S", 1:4),
    Age = c(45, 55, 65, 70),
    stringsAsFactors = FALSE
  )

  set.seed(42)
  for (uniprot in protein_map$uniprot) {
    demo[[uniprot]] <- stats::rnorm(4)
  }

  result <- compute_goeminne2025_organ_age(
    demo,
    organs = "Conventional",
    match_by = "uniprot",
    protein_map = protein_map
  )

  expect_equal(nrow(result), 4)
  expect_equal(result$match_by, rep("uniprot", 4))
  expect_true(all(result$n_proteins_matched > 0))
  expect_true(all(is.finite(result$proteomic_age)))
})

test_that("compute_goeminne2025_organ_age uses built-in Olink UniProt map", {
  protein_map <- olink_protein_map()
  prots <- goeminne2025_organaging_proteins()
  needed <- unique(prots$protein[prots$organ == "Conventional" &
                                  prots$model_type == "chronological" &
                                  prots$fold == 1])
  needed <- setdiff(head(needed, 25), "Intercept")
  map <- protein_map[toupper(protein_map$gene) %in% toupper(needed), ]
  map <- map[!duplicated(toupper(map$gene)), ]
  demo <- data.frame(
    SampleID = paste0("S", 1:4),
    Age = c(45, 55, 65, 70),
    stringsAsFactors = FALSE
  )

  set.seed(42)
  for (uniprot in map$uniprot) {
    demo[[uniprot]] <- stats::rnorm(4)
  }

  result <- compute_goeminne2025_organ_age(
    demo,
    organs = "Conventional",
    match_by = "uniprot"
  )

  expect_equal(nrow(result), 4)
  expect_equal(result$match_by, rep("uniprot", 4))
  expect_true(all(result$n_proteins_matched > 0))
  expect_true(all(is.finite(result$proteomic_age)))
})

test_that("Olink gene matching accepts R-syntactic gene column names", {
  lookup <- proteomicAge:::.olink_gene_lookup(c("SampleID", "Age", "ERVV.1"), "ERVV-1")

  expect_equal(unname(lookup[["ERVV-1"]]), "ERVV.1")
})

test_that("olink_protein_map_from_long extracts AIFI-style Olink maps", {
  long <- data.frame(
    olink.assay = c("GDF15", "GDF15", "IL6"),
    olink.uniprot_id = c("Q99988", "Q99988", "P05231"),
    stringsAsFactors = FALSE
  )

  protein_map <- olink_protein_map_from_long(long)

  expect_equal(nrow(protein_map), 2)
  expect_true(all(c("gene", "uniprot") %in% names(protein_map)))
  expect_true("GDF15" %in% protein_map$gene)
})

test_that("built-in Olink map covers PAC and Goeminne proteins", {
  protein_map <- olink_protein_map()
  pac <- toupper(kuo2024_pac_proteins()$predictor)
  goeminne <- toupper(goeminne2025_organaging_proteins()$protein)

  expect_true(all(pac %in% toupper(protein_map$gene)))
  expect_true(all(goeminne %in% toupper(protein_map$gene)))
})
