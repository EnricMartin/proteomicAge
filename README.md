# proteomicAge

Compute biological age from plasma proteomic data using published proteomic
aging clocks.

## Installation

```r
remotes::install_github("EnricMartin/proteomicAge")
```

## Supported Clocks

| Clock | Proteins | Method | r | Reference |
|-------|----------|--------|---|-----------|
| **Tanaka 2018** | 76 | Elastic Net | 0.94 | Tanaka et al. Aging Cell (2018) |
| **Lehallier 2019** | 373 | LASSO | 0.93-0.97 | Lehallier et al. Nat Med (2019) |

## Input Format

The input data.frame must contain the following columns:

| Column | Description | Example |
|--------|-------------|---------|
| Sample ID | Unique sample identifier (any name) | `"P001"` |
| Age | Chronological age in years | `50` |
| Sex | (optional) 0 = male, 1 = female | `0` |
| Protein columns | Named using one of the supported conventions below | |

**Supported protein naming conventions:**

| Convention | Example | Description |
|---|---|---|
| UniProt accession | `P36222`, `Q99988` | Standard UniProtKB identifier |
| Gene symbol | `CHI3L1`, `GDF15` | Entrez/NCBI gene symbol |
| SeqId (dot format) | `seq.11104.13.3`, `11104.13.3` | SomaScan SeqId with or without `seq.` prefix |
| SeqId (SL format) | `SL003340`, `SL003869` | SomaScan legacy SL-format SeqId |

Column naming is detected automatically or can be specified explicitly via
the `match_by` parameter.

## Quick Start

```r
library(proteomicAge)

# Input with any supported protein naming
data <- data.frame(
  SampleID = "S1",
  Age = 50,
  CHI3L1 = 3000,
  GDF15  = 5000,
  SOST   = 2000
)

# Compute — package auto-detects gene symbols
result <- compute_tanaka2018_age(data, match_by = "gene")
result <- compute_lehallier2019_age(data, match_by = "gene")

head(result)
```

## Log Transform

| Clock | Transform |
|-------|-----------|
| **Tanaka 2018** | log₂ |
| **Lehallier 2019** | log₁₀ |

## Age Acceleration

PROaccel = residuals of `lm(PROage ~ chronological_age)`. Mean ≈ 0.
Positive = predicted older than expected.

## Output

| Column | Description |
|--------|-------------|
| `id` | Sample identifier |
| `chronological_age` | Input age |
| `proteomic_age` | Predicted biological age |
| `age_acceleration` | PROaccel |
| `n_proteins_matched` | Clock proteins found |
| `n_proteins_missing` | Clock proteins missing |
| `match_by` | Naming convention used |

## Methodology

**Tanaka et al. (2018):** 1,301 SOMAscan proteins, 240 healthy adults (BLSA + GESTALT).
Elastic Net (α=0.5, λ=0.8767859, 10-fold CV), log₂ transform. 76 proteins; r=0.94.

**Lehallier et al. (2019):** 2,925 SOMAscan proteins, 4,263 adults (INTERVAL + LonGenity).
LASSO (α=1.0, λ.min, 10-fold CV), log₁₀ transform + Z-scaling. 373 proteins; r=0.93-0.97.

## Citation

```
Tanaka T, et al. Plasma proteomic signature of age in healthy humans.
Aging Cell. 2018;17(5):e12799.

Lehallier B, et al. Undulating changes in human plasma proteome profiles
across the lifespan. Nat Med. 2019;25(12):1843-1850.
```

## License

MIT
