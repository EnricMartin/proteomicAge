# proteomicAge

Compute biological age from plasma proteomic data using published proteomic
aging clocks.

**Contributed by:** Han Xiao (hx624@ic.ac.uk), Esther Herrera, Arias Julian,
Juan-Carlos Rivilla, Oliver Robinson.

## Installation

```r
remotes::install_github("EnricMartin/proteomicAge")
```

## Supported Clocks

| Clock | Proteins | Method | r | Reference |
|-------|----------|--------|---|-----------|
| **Tanaka 2018** | 76 | Elastic Net | 0.94 | Tanaka et al. Aging Cell (2018) |
| **Lehallier 2019** | 373 | LASSO | 0.93-0.97 | Lehallier et al. Nat Med (2019) |
| **Sathyan 2020** | 4,265 | Elastic Net | — | Sathyan et al. Aging Cell (2020) |

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

## Usage

### Step 1: Detect protein naming format

```r
library(proteomicAge)

# Your data — protein columns in any supported convention
data <- read.csv("my_somascan_data.csv")

# Auto-detect the naming convention
protein_cols <- setdiff(names(data), c("SampleID", "Age", "Sex"))
fmt <- detect_format(protein_cols)
print(fmt)  # e.g., "gene", "uniprot", "seqid_sl", "seqid_dot"
```

### Step 2: (Optional) Convert to a different format

```r
# Convert gene symbols to UniProt accessions
data <- convert_format(data, target_format = "uniprot",
                       id_col = "SampleID", age_col = "Age")
```

### Step 3: Compute proteomic age

```r
# Tanaka 2018 clock — specify how to match your column names
result_tanaka <- compute_tanaka2018_age(
  data,
  id_col = "SampleID",
  age_col = "Age",
  match_by = fmt
)

# Lehallier 2019 clock (includes sex as covariate)
result_leh <- compute_lehallier2019_age(
  data,
  id_col = "SampleID",
  age_col = "Age",
  sex_col = "Sex",
  match_by = fmt
)

head(result_tanaka)
```

### Demo with synthetic data

```r
# Generate realistic demo data with many proteins
prots <- tanaka2018_proteins()
demo <- data.frame(
  SampleID = paste0("P", 1:10),
  Age = c(32, 45, 51, 63, 71, 38, 55, 67, 42, 78),
  Sex = c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0)
)
set.seed(123)
for (i in 1:min(nrow(prots), 20)) {
  demo[[prots$Gene[i]]] <- round(rlnorm(10, meanlog = log(2000), sdlog = 0.5))
}

result <- compute_tanaka2018_age(demo, match_by = "gene")
print(result[, c("id", "chronological_age", "proteomic_age", "age_acceleration")])
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

**Sathyan et al. (2020):** 4,265 SOMAscan v4 proteins, 1,025 older adults (LonGenity).
Elastic Net (α=0.5, 10-fold CV), log transform. Model: log(SOMAmer) ~ age + gender + cohort.

## Citation

```
Tanaka T, et al. Plasma proteomic signature of age in healthy humans.
Aging Cell. 2018;17(5):e12799.

Lehallier B, et al. Undulating changes in human plasma proteome profiles
across the lifespan. Nat Med. 2019;25(12):1843-1850.

Sathyan S, et al. Plasma proteomic profile of age, health span, and
all-cause mortality in older adults. Aging Cell. 2020;19(11):e13250.
```

## License

MIT
