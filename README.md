# proteomicAge

An R package for computing biological age from plasma proteomic data using
published proteomic aging clocks. **All computation is local — your proteomic
data never leaves your machine.**

## Currently Implemented Clocks

| Clock | Proteins | Platform | Method | r | Reference |
|-------|----------|----------|--------|---|-----------|
| **Tanaka 2018** | 76 | SOMAscan 1.3K | Elastic Net | 0.94 | [Tanaka et al. (2018) Aging Cell](https://doi.org/10.1111/acel.12799) |
| **Lehallier 2019** | 373 | SOMAscan 2.9K | LASSO | 0.93-0.97 | [Lehallier et al. (2019) Nature Medicine](https://doi.org/10.1038/s41591-019-0673-2) |

## Installation

```r
remotes::install_github("ClawBio/ClawBio", subdir = "proteomicAge")
```

For full protein name mapping (recommended):
```r
BiocManager::install("SomaScan.db")
```

## Quick Start

```r
library(proteomicAge)

# Your data can use any protein identifier format
data <- read.csv("my_somascan_data.csv")

# Compute proteomic age with automatic format detection
result <- compute_tanaka2018_age(data, protein_format = "auto")
result <- compute_lehallier2019_age(data, protein_format = "auto")

head(result)
#   id chronological_age proteomic_age age_acceleration
# 1 S1                50          52.3              1.8
# 2 S2                65          61.1             -4.2
```

## Protein Name Mapping

The package accepts protein columns in **any common format** and automatically
normalizes them to the clock's expected format.

### Supported Input Formats

| `protein_format` | Example column names | Primary mapper |
|---|---|---|
| `"gene"` | `GDF15`, `CHI3L1` | SomaScan.db ALIAS |
| `"uniprot"` | `Q99988`, `P36222` | SomaScan.db UNIPROT |
| `"seqid_dot"` | `4374.45.2` | SomaScan.db PROBEID |
| `"seqid_sl"` | `SL003869` | Built-in mapping |
| `"seqid_full"` | `GDF15.4374.45.2` | Direct use |
| `"auto"` (default) | (automatic detection) | Best match |

### Architecture

```
Your data columns (gene / uniprot / seqid_dot / seqid_sl / seqid_full)
        │
        ▼
┌──────────────────────────────┐
│  SomaScan.db (PRIMARY)       │  ← Install with BiocManager::install("SomaScan.db")
│  • gene → PROBEID → gene     │    Covers all SOMAscan panels (5k, 7k, 11k)
│  • uniprot → gene + PROBEID  │
│  • seqid_dot ↔ PROBEID       │
└──────────────────────────────┘
        │  (SL format is v3-only — not in SomaScan.db)
        ▼
┌──────────────────────────────┐
│  Built-in CSV (FALLBACK)     │  ← Always available, zero dependencies
│  • SL format ↔ gene symbols  │    Covers the 76 + 373 clock proteins
│  • Multi-gene disambiguation │
└──────────────────────────────┘
        │
        ▼
   Clock's expected format (seqid_sl or seqid_full)
```

### Usage examples

```r
# Gene symbol input → Tanaka clock
data <- data.frame(SampleID = "S1", Age = 50, GDF15 = 5000, CHI3L1 = 3000)
result <- compute_tanaka2018_age(data, protein_format = "gene")

# UniProt input → Lehallier clock
data <- data.frame(SampleID = "S1", Age = 50, Q99988 = 5000, P36222 = 3000)
result <- compute_lehallier2019_age(data, protein_format = "uniprot")

# Dot-format SeqId input (seq.xxx.xx) → Lehallier clock
data <- data.frame(SampleID = "S1", Age = 50, `4374.45.2` = 5000)
result <- compute_lehallier2019_age(data, protein_format = "seqid_dot")
```

## Log Transform

Different clocks use different log transforms — the package applies the correct one automatically:

| Clock | Transform | Details |
|-------|-----------|---------|
| **Tanaka 2018** | log₂ | Each protein was log₂-transformed before elastic net fitting (Tanaka et al. 2018) |
| **Lehallier 2019** | log₁₀ | LASSO model fitted on Z-scaled log₁₀ RFU values + sex covariate (Lehallier et al. 2019) |

You can disable the built-in transform if your data is already log-scaled:

```r
compute_tanaka2018_age(data, log_transform = FALSE)
```

## Age Acceleration

Age acceleration (PROaccel) is computed as the **residuals of the linear
regression of predicted proteomic age on chronological age**:

```
PROaccel = residuals( lm(PROage ~ chronological_age) )
```

This follows the original papers (Tanaka et al. 2020 eLife; Lehallier et al.
2019 Nature Medicine). The mean age acceleration across a population is
approximately zero by construction.

> **Positive age acceleration**: predicted older than expected for one's age  
> **Negative age acceleration**: predicted younger than expected for one's age

## Preprocessing

For manual control over preprocessing:

```r
processed <- preprocess_somascan(raw_data,
  log_transform = TRUE,
  handle_outliers = TRUE
)
result <- compute_tanaka2018_age(processed, log_transform = FALSE)
```

## Output

All `compute_*_age()` functions return a data.frame with:

| Column | Description |
|--------|-------------|
| `id` | Sample identifier |
| `chronological_age` | Input age |
| `proteomic_age` | Predicted biological age (PROage) |
| `age_acceleration` | PROaccel = residuals(PROage ~ chronological_age) |
| `n_proteins_matched` | Number of required proteins found in data |
| `n_proteins_missing` | Number of required proteins missing |

## Methodology

### Tanaka et al. (2018)

- **Data**: 1,301 SOMAscan proteins measured in 240 healthy adults (22-93 years)
  from BLSA and GESTALT cohorts
- **Training**: 120 subjects, stratified by 15-year age strata
- **Model**: Elastic Net regression (`glmnet`, α = 0.5, λ = 0.8767859, 10-fold CV)
- **Preprocessing**: log₂ transformation; outliers ±4 SD removed
- **Coefficients**: 76 proteins selected from 1,301. Weights from eLife (2020) Supplementary File 1D

### Lehallier et al. (2019)

- **Data**: 2,925 SOMAscan proteins measured in 4,263 healthy adults (18-95 years)
  from INTERVAL and LonGenity cohorts
- **Training**: 2,817 subjects; validation: 1,446 subjects
- **Model**: LASSO regression (`glmnet`, α = 1.0, λ.min, 10-fold CV)
- **Preprocessing**: Z-scaled log₁₀ RFU values + sex as covariate
- **Coefficients**: 373 proteins. Weights from Supplementary Table 7

## Citation

If you use this package, please cite both the package and the original clock paper:

```
Tanaka T, Basisty N, Fantoni G, Candia J, Moore AZ, Biancotto A,
Schilling B, Bandinelli S, Ferrucci L. Plasma proteomic signature
of age in healthy humans. Aging Cell. 2018 Oct;17(5):e12799.
doi: 10.1111/acel.12799.

Lehallier B, Gate D, Schaum N, Nanasi T, Lee SE, Yousef H,
Moran Losada P, Berdnik D, Keller A, Verghese J, Sathyan S,
Franceschi C, Milman S, Barzilai N, Wyss-Coray T. Undulating
changes in human plasma proteome profiles across the lifespan.
Nature Medicine. 2019 Dec;25(12):1843-1850.
doi: 10.1038/s41591-019-0673-2.
```

## License

MIT
