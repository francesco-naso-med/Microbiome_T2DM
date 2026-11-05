# Gut microbiome & type 2 diabetes — a leakage-free ML analysis across two populations

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-curatedMetagenomicData-276DC3?logo=r)](R/)
[![Python](https://img.shields.io/badge/Python-scikit--learn%20%7C%20XGBoost%20%7C%20SHAP-3776AB?logo=python)](python/)

Can the gut microbiome tell people with type 2 diabetes (T2D) from controls, is the signal really the disease rather than metformin, and does it hold in a different population? This repository answers these three questions on public shotgun-metagenomic data with a reproducible, leakage-free pipeline. Full analysis: [`REPORT.md`](REPORT.md).

## Results at a glance

| | Random Forest | XGBoost | Elastic net |
|---|---|---|---|
| Internal — `QinJ_2012`, China (170 T2D vs 174 controls, 10×10 CV) — AUROC | **0.73** ± 0.01 | **0.73** ± 0.01 | 0.61 ± 0.01 |
| External — `KarlssonFH_2013`, Sweden (53 T2D vs 43 controls) — AUROC [95% CI] | 0.49 [0.38–0.62] | 0.49 [0.38–0.61] | 0.35 [0.25–0.46] |

- **Within a cohort, the signal is real and biologically coherent:** depletion of butyrate producers (*Roseburia*, *Faecalibacterium*, *Anaerostipes*) drives the prediction.
- **It is not a metformin artifact:** treated cases are not easier to classify than untreated or drug-free ones.
- **It does not transfer across populations:** at chance in Swedish women even after excluding impaired glucose tolerance and restricting to participants on neither metformin nor statins, and significantly inverted for the linear model.

## Repository contents

| Path | Contents |
|------|----------|
| [`REPORT.md`](REPORT.md) | The analysis report: background, methods, results, discussion, limitations, references. |
| [`R/01_fetch_curatedMetagenomicData.R`](R/01_fetch_curatedMetagenomicData.R) | Downloads `QinJ_2012` and `KarlssonFH_2013`; exports species profiles and metadata, with three-state metformin and statin flags. |
| [`python/clr.py`](python/clr.py) | Compositional preprocessing as a scikit-learn transformer: prevalence filter, multiplicative zero replacement, CLR — fitted on training folds only. |
| [`python/03_train_rf_shap.py`](python/03_train_rf_shap.py) | T2D vs control; repeated stratified CV (nested with `--tune`); Random Forest, XGBoost, elastic net; medication subgroup AUROCs with bootstrap CIs; external validation. |
| [`python/04_interpret_shap.py`](python/04_interpret_shap.py) | TreeSHAP importance and direction, beeswarm, per-patient waterfalls, dependence plot. |
| [`R/02_preprocess_clr.R`](R/02_preprocess_clr.R) | Optional, exploratory only; not used by the ML pipeline. |
| [`docs/references.bib`](docs/references.bib) | BibTeX for the references in the report. |

## Quickstart

**1. Data (R ≥ 4.3, Bioconductor)**
```r
# BiocManager::install(c("curatedMetagenomicData", "mia"))
source("R/01_fetch_curatedMetagenomicData.R")   # writes data/processed/*.csv
```

**2. Models and interpretation (Python ≥ 3.10)**
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python python/03_train_rf_shap.py --cohort QinJ_2012 --external KarlssonFH_2013
python python/04_interpret_shap.py --cohort QinJ_2012
```
A full 10×10 run takes several minutes on a laptop; add `--folds 5 --repeats 2` for a quick check. On macOS, XGBoost needs OpenMP (`brew install libomp`); without it the XGBoost model is skipped. Outputs: `artifacts/metrics.json`, `artifacts/shap_ranked.csv`, `figures/`.

## Scope and limitations

- Two cohorts; treatment not recorded for most of `QinJ_2012`; fixed hyperparameters; no batch or covariate correction; species-level taxonomy only. See [`REPORT.md` §5](REPORT.md#5-limitations).
- No patient-level data are stored here: all data are downloaded at runtime from the public `curatedMetagenomicData` resource, under its own license.
- A research artifact, not clinical guidance.

## Citation

Please cite the repository URL: https://github.com/francesco-naso-med/Microbiome_T2DM

## License

[MIT](LICENSE). Underlying datasets and tools (`curatedMetagenomicData`, MetaPhlAn) retain their own licenses.
