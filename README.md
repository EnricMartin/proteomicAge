# proteomicAge

Compute biological age from plasma proteomic data using published proteomic
aging clocks. **UniProt accession is the required input format.**

## Installation

```r
remotes::install_github("EnricMartin/proteomicAge")
```

## Supported Clocks

| Clock | Proteins | Method | r | Reference |
|-------|----------|--------|---|-----------|
| **Tanaka 2018** | 76 | Elastic Net | 0.94 | Tanaka et al. Aging Cell (2018) |
| **Lehallier 2019** | 373 | LASSO | 0.93-0.97 | Lehallier et al. Nat Med (2019) |

## Quick Start

```r
library(proteomicAge)

# Input: data.frame with UniProt-named protein columns
data <- data.frame(
  SampleID = "S1",
  Age = 50,
  P36222 = 3000,   # CHI3L1
  Q99988 = 5000,   # GDF15
  Q9BQB4 = 2000    # SOST
)

# Tanaka 2018 clock
result <- compute_tanaka2018_age(data)
print(result)

# Lehallier 2019 clock
result_lh <- compute_lehallier2019_age(data)
print(result_lh)
```

## Input Format

**All protein columns must be named by UniProt accession.**

| Column | Description | Example |
|--------|-------------|---------|
| `SampleID` | Sample identifier | `"P001"` |
| `Age` | Chronological age (years) | `50` |
| `Sex` | Sex (0=male, 1=female) | `0` |
| `P36222` | CHI3L1 protein | `3000` |
| `Q99988` | GDF15 protein | `5000` |

Protein values should be in RFU (relative fluorescence units).
Log transformation is applied automatically per the clock's methodology.

## Log Transform

| Clock | Transform |
|-------|-----------|
| **Tanaka 2018** | log₂ |
| **Lehallier 2019** | log₁₀ |

## Age Acceleration

PROaccel = residuals of `lm(PROage ~ chronological_age)` — matches the original papers.

## Output

| Column | Description |
|--------|-------------|
| `id` | Sample identifier |
| `chronological_age` | Input age |
| `proteomic_age` | Predicted biological age |
| `age_acceleration` | PROaccel (mean ≈ 0) |
| `n_proteins_matched` | Clock proteins found in data |
| `n_proteins_missing` | Clock proteins not found |

## Methodology

### Tanaka et al. (2018)
- 1,301 SOMAscan proteins, 240 healthy adults (BLSA + GESTALT)
- Elastic Net (α=0.5, λ=0.8767859, 10-fold CV), log₂ transform
- 76 proteins selected; r=0.94

### Lehallier et al. (2019)
- 2,925 SOMAscan proteins, 4,263 adults (INTERVAL + LonGenity)
- LASSO (α=1.0, λ.min, 10-fold CV), log₁₀ transform + Z-scaling
- 373 proteins; r=0.93-0.97

## Citation

```
Tanaka T, et al. Plasma proteomic signature of age in healthy humans.
Aging Cell. 2018;17(5):e12799. doi: 10.1111/acel.12799.

Lehallier B, et al. Undulating changes in human plasma proteome profiles
across the lifespan. Nat Med. 2019;25(12):1843-1850.
doi: 10.1038/s41591-019-0673-2.
```

## License

MIT
