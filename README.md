# proteomicAge

An R package for computing biological age from plasma proteomic data using
published proteomic aging clocks. **All computation is local — your proteomic
data never leaves your machine.**

## Currently Implemented Clocks

| Clock | Proteins | Platform | Method | r | Reference |
|-------|----------|----------|--------|---|-----------|
| **Tanaka 2018** | 76 | SOMAscan 1.3K | Elastic Net | 0.94 | [Tanaka et al. (2018) Aging Cell](https://doi.org/10.1111/acel.12799) |
| **Lehallier 2019** | 373 | SOMAscan 2.9K | LASSO | 0.93-0.97 | [Lehallier et al. (2019) Nature Medicine](https://doi.org/10.1038/s41591-019-0673-2) |
| **Lehallier 2020** | 491 | SOMAscan 2.9K | LASSO | 0.98 | [Lehallier et al. (2020) Aging Cell](https://doi.org/10.1111/acel.13256) |

### Planned Clocks

- ProtAge / ProtAge20 (Nature Medicine 2024, Olink)
- Sathyan 2020 (162 proteins, SOMAscan)

## Installation

```r
# Install from GitHub
remotes::install_github("ClawBio/ClawBio", subdir = "proteomicAge")

# Or install locally
devtools::install("proteomicAge")
```

## Quick Start

```r
library(proteomicAge)

# List the 76 proteins in the Tanaka clock
tanaka2018_proteins()

# Read SOMAscan data (RFU values, samples x proteins)
data <- read.csv("my_somascan_data.csv")

# Validate input
validate_somascan_input(data, id_col = "SampleID", age_col = "Age")

# Compute proteomic age
results <- compute_tanaka2018_age(data,
  id_col = "SampleID",
  age_col = "Age"
)

head(results)
#   id chronological_age proteomic_age age_acceleration
# 1 S1                50          52.3              2.3
# 2 S2                65          61.1             -3.9
# ...
```

## Input Format

The input data.frame should contain:

- **ID column** — unique sample identifier
- **Age column** — chronological age in years
- **Protein columns** — named with SOMAscan SeqIds (e.g., `SL003869`, `SL000045`)

Protein values should be in RFU (relative fluorescence units). The function
will log2-transform them internally, matching the original paper's methodology.

```
  SampleID Age SL003869 SL000045 SL000254 ...
1     S001  45     1234      890     5678 ...
2     S002  62     1567      750    12345 ...
```

## Preprocessing

For more control over preprocessing:

```r
# Preprocess raw SOMAscan data
processed <- preprocess_somascan(raw_data,
  log_transform = TRUE,
  handle_outliers = TRUE
)

# Then compute age
results <- compute_tanaka2018_age(processed,
  log_transform = FALSE  # already transformed
)
```

## Output

`compute_tanaka2018_age()` returns a data.frame with:

| Column | Description |
|--------|-------------|
| `id` | Sample identifier |
| `chronological_age` | Input age |
| `proteomic_age` | Predicted biological age (PROage) |
| `age_acceleration` | PROaccel = proteomic_age - chronological_age |
| `n_proteins_matched` | Number of the 76 required proteins found |
| `n_proteins_missing` | Number of required proteins missing |

> **Positive age acceleration**: predicted older than actual age  
> **Negative age acceleration**: predicted younger than actual age

## Methodology

The Tanaka et al. (2018) clock was developed using:

- **Data**: 1,301 SOMAscan proteins measured in 240 healthy adults (22-93 years)
  from BLSA and GESTALT cohorts
- **Training**: 120 subjects, stratified by 15-year age strata
- **Validation**: 120 subjects
- **Model**: Elastic Net regression (R package `glmnet`)
  - α (elastic net mixing) = 0.5
  - λ (regularization) = 0.8767859, selected by 10-fold cross-validation
  - 76 proteins selected from 1,301 candidates
- **Preprocessing**: Each protein was log-transformed; outliers ±4 SD were removed
- **Performance**: Pearson r = 0.94 between predicted and observed age

## Citation

If you use this package in your research, please cite:

```
Tanaka T, Basisty N, Fantoni G, Candia J, Moore AZ, Biancotto A,
Schilling B, Bandinelli S, Ferrucci L. Plasma proteomic signature
of age in healthy humans. Aging Cell. 2018 Oct;17(5):e12799.
doi: 10.1111/acel.12799.

Tanaka T, Basisty N, Fantoni G, et al. Plasma proteomic biomarker
signature of age predicts health and life span. eLife. 2020;9:e61073.
doi: 10.7554/eLife.61073.
```

## License

MIT
