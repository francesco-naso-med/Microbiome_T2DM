# Data

**No raw or patient-level data are committed to this repository, by design.**

All inputs are pulled at runtime from the public
[`curatedMetagenomicData`](https://waldronlab.io/curatedMetagenomicData/) resource
(Pasolli et al., *Nat Methods* 2017; doi:10.1038/nmeth.4468) via `R/01_fetch_curatedMetagenomicData.R`.

## Why nothing is stored here

1. **Privacy / re-identification.** Shotgun metagenomes are individually
   identifying (Franzosa et al., *PNAS* 2015) and can carry host-genome reads.
   See `REPORT.md` §7e.
2. **Licensing.** `curatedMetagenomicData` and the underlying studies retain
   their own licenses and citation requirements.
3. **Reproducibility.** Pulling from a versioned resource at runtime beats
   committing a frozen, undocumented copy.

## Expected layout after running the pipeline

```
data/
├── README.md                      (this file)
└── processed/                     (git-ignored, created by R/01_ and R/02_)
    ├── QinJ_2012_species.csv
    ├── QinJ_2012_species_clr.csv
    ├── QinJ_2012_metadata.csv
    ├── KarlssonFH_2013_species.csv
    ├── KarlssonFH_2013_species_clr.csv
    └── KarlssonFH_2013_metadata.csv
```

## Cohorts used

| study_name | Condition | Population | Role |
|------------|-----------|-----------|------|
| `QinJ_2012` | T2DM | China | Discovery / training |
| `KarlssonFH_2013` | T2DM, IGT, NGT | Sweden (women) | External validation |

To list every available study and condition in your installed package version:

```r
library(curatedMetagenomicData)
sampleMetadata |>
  dplyr::filter(study_condition %in% c("T2D", "IGT")) |>
  dplyr::count(study_name, study_condition)
```

Confirm sample counts against your package release — the figures cited in
`REPORT.md` flagged `[VERIFY]` should be checked here.
