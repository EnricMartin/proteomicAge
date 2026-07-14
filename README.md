# proteomicAge

Compute biological age from plasma proteomic data using published proteomic
aging clocks.

**Contributed by:** Han Xiao (hx624@ic.ac.uk), Esther Herrera, Arias Julian,
Juan-Carlos Rivilla, Oliver Robinson.

## Installation

```r
remotes::install_github("EnricMartin/proteomicAge")
```

## Supported clocks

`proteomicAge` currently provides five published proteomic aging clocks.
Each clock has a `compute_*_age()` function for prediction and a matching
`*_proteins()` function for listing the proteins used by that model.

| Clock | Function | Proteins | Default `match_by` | Transform | Notes |
|-------|----------|----------|--------------------|-----------|-------|
| Tanaka 2018 | `compute_tanaka2018_age()` | 76 | `uniprot` | `log2` | Elastic net clock trained on SOMAscan 1.3K data |
| Lehallier 2019 | `compute_lehallier2019_age()` | 373 | `uniprot` | `log10` | LASSO clock; accepts an optional sex column |
| Sathyan 2020 | `compute_sathyan2020_age()` | 162 | `uniprot` | natural log | Elastic net clock from LonGenity SOMAscan v4 data |
| Oh 2023 conventional | `compute_oh2023_conventional_age()` | 4,778 | `seqid_dot` | natural log | Conventional whole-body clock from the organ aging study |
| Wang 2024 ARIC midlife | `compute_wang2024_aric_age()` | 788 | `seqid_dot` | `log2` | Elastic net clock trained in the ARIC midlife cohort |

## Input format

The input data frame should contain one row per sample, sample metadata columns,
and one column per protein.

| Column | Required | Description | Example |
|--------|----------|-------------|---------|
| Sample ID | Yes | Unique sample identifier; name is set with `id_col` | `"P001"` |
| Age | Yes | Chronological age in years; name is set with `age_col` | `50` |
| Sex | Optional | Used by Lehallier 2019 if supplied; default male value is `0` | `0` |
| Protein columns | Yes | Protein abundance values named with one supported convention | `GDF15` |

Supported protein naming conventions:

| `match_by` value | Example | Description |
|------------------|---------|-------------|
| `uniprot` | `P36222`, `Q99988` | UniProtKB accession |
| `gene` | `CHI3L1`, `GDF15` | Gene symbol |
| `seqid_dot` | `seq.11104.13.3`, `11104.13.3` | SomaScan dot-format SeqId, with or without `seq.` prefix |
| `seqid_sl` | `SL003340`, `SL003869` | SomaScan legacy SL-format SeqId |

Use `detect_format()` to infer the naming convention from protein column names,
or pass `match_by` explicitly to any clock function.

## Quick start

```r
library(proteomicAge)

dat <- read.csv("my_somascan_data.csv")

protein_cols <- setdiff(names(dat), c("SampleID", "Age", "Sex"))
fmt <- detect_format(protein_cols)
fmt
```

Run one clock:

```r
tanaka_age <- compute_tanaka2018_age(
  dat,
  id_col = "SampleID",
  age_col = "Age",
  match_by = fmt
)

head(tanaka_age)
```

Run all five clocks:

```r
tanaka_age <- compute_tanaka2018_age(
  dat,
  id_col = "SampleID",
  age_col = "Age",
  match_by = fmt
)

lehallier_age <- compute_lehallier2019_age(
  dat,
  id_col = "SampleID",
  age_col = "Age",
  sex_col = "Sex",
  male_value = 0,
  match_by = fmt
)

sathyan_age <- compute_sathyan2020_age(
  dat,
  id_col = "SampleID",
  age_col = "Age",
  sex_col = "Sex",
  match_by = fmt
)

oh_age <- compute_oh2023_conventional_age(
  dat,
  id_col = "SampleID",
  age_col = "Age",
  match_by = fmt
)

wang_age <- compute_wang2024_aric_age(
  dat,
  id_col = "SampleID",
  age_col = "Age",
  match_by = fmt
)
```

List the proteins required by a clock:

```r
tanaka2018_proteins()
lehallier2019_proteins()
sathyan2020_proteins()
oh2023_conventional_proteins()
wang2024_aric_proteins()
```

Convert protein column names when needed:

```r
dat_uniprot <- convert_format(
  dat,
  target_format = "uniprot",
  id_col = "SampleID",
  age_col = "Age"
)
```

## Output

Each `compute_*_age()` function returns a data frame with the same standard
columns.

| Column | Description |
|--------|-------------|
| `id` | Sample identifier |
| `chronological_age` | Input age |
| `proteomic_age` | Predicted biological age |
| `age_acceleration` | Proteomic age acceleration |
| `n_proteins_matched` | Number of clock proteins found in the input data |
| `n_proteins_missing` | Number of clock proteins not found in the input data |
| `match_by` | Protein naming convention used for matching |

Age acceleration is computed as the residual from:

```r
lm(proteomic_age ~ chronological_age)
```

Positive values indicate a predicted proteomic age older than expected for the
sample's chronological age.

## Demo with synthetic data

```r
prots <- tanaka2018_proteins()

demo <- data.frame(
  SampleID = paste0("P", 1:10),
  Age = c(32, 45, 51, 63, 71, 38, 55, 67, 42, 78),
  Sex = c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0)
)

set.seed(123)
for (i in seq_len(nrow(prots))) {
  demo[[prots$Gene[i]]] <- round(rlnorm(10, meanlog = log(2000), sdlog = 0.5))
}

result <- compute_tanaka2018_age(demo, match_by = "gene")
result[, c("id", "chronological_age", "proteomic_age", "age_acceleration")]
```

## Methodology

**Tanaka et al. (2018):** 1,301 SOMAscan proteins, 240 healthy adults from
BLSA and GESTALT. Elastic net model with 76 selected proteins.

**Lehallier et al. (2019):** 2,925 SOMAscan proteins, 4,263 adults from
INTERVAL and LonGenity. LASSO model with 373 selected proteins.

**Sathyan et al. (2020):** 4,265 SOMAscan v4 proteins, 1,025 older adults from
the LonGenity cohort. Elastic net model with 162 selected proteins.

**Oh et al. (2023):** SOMAscan v4 organ aging study. The conventional
proteomic age model is implemented as `compute_oh2023_conventional_age()`.

**Wang et al. (2024):** Population-based proteomic aging clock development
and replication study. The ARIC midlife model is implemented as
`compute_wang2024_aric_age()`.

## Citation

```text
Tanaka T, et al. Plasma proteomic signature of age in healthy humans.
Aging Cell. 2018;17(5):e12799.

Lehallier B, et al. Undulating changes in human plasma proteome profiles
across the lifespan. Nature Medicine. 2019;25(12):1843-1850.

Sathyan S, et al. Plasma proteomic profile of age, health span, and
all-cause mortality in older adults. Aging Cell. 2020;19(11):e13250.

Oh HS, et al. Organ aging signatures in the plasma proteome track health
and disease. Nature. 2023;624:164-172.

Wang S, et al. Development, characterization, and replication of proteomic
aging clocks: analysis of 2 population-based cohorts. PLOS Medicine.
2024;21(9):e1004464.
```

## License

MIT
