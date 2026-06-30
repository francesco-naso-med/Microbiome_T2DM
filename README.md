# Gut Microbiome & Type 2 Diabetes — Mechanisms, Evidence, and a Reproducible ML Roadmap [W.I.P.]

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Made with R](https://img.shields.io/badge/R-curatedMetagenomicData-276DC3?logo=r)](R/)
[![Made with Python](https://img.shields.io/badge/Python-RandomForest%20%2B%20SHAP-3776AB?logo=python)](python/)
[![Status](https://img.shields.io/badge/status-research%20artifact-success)](REPORT.md)

> A self-contained, citation-disciplined review of the gut microbiome in Type 2 Diabetes Mellitus (T2DM), **paired with a runnable bioinformatics pipeline** for building an interpretable predictive model from public shotgun-metagenomic data. Designed to double as the theoretical + computational foundation for a Master's/MD thesis.

---

## Why this repository exists

Most microbiome reviews stop at biology. Most ML tutorials stop at toy data. This repo connects the two: a **Nature-Medicine-level mechanistic review** ([`REPORT.md`](REPORT.md)) whose §6 "Analysis Roadmap" maps one-to-one onto **working code skeletons** ([`/R`](R/), [`/python`](python/)) that load `curatedMetagenomicData`, transform it correctly (CLR for compositional data), train Random Forest / ElasticNet / XGBoost under honest nested cross-validation, and explain predictions with SHAP.

The central scientific thesis it advances: **much of the published "T2DM microbiome signature" is confounded by metformin** (Forslund et al., *Nature* 2015) — and the open, fundable frontier is whether the **baseline microbiome predicts response to GLP-1 receptor agonists** (semaglutide, tirzepatide), for which no public shotgun cohort yet exists.

## What's inside

| Path | Contents |
|------|----------|
| [`REPORT.md`](REPORT.md) | The full review: executive summary, 5 mechanistic axes (SCFA/GLP-1, LPS/TLR4, bile acids/FXR-TGR5, BCAA/mTOR, barrier), a 22-row taxa evidence table (Oxford CEBM grading), a 25-study landmark timeline (2006–2025), a deep GLP-1-RA × microbiome section, the ML roadmap, clinical translation, 7 research gaps, 50+ Vancouver references, and a thesis-integration note. |
| [`R/01_fetch_curatedMetagenomicData.R`](R/01_fetch_curatedMetagenomicData.R) | Pull `QinJ_2012` (discovery) + `KarlssonFH_2013` (external validation) as `TreeSummarizedExperiment`; export species + pathway matrices and metadata. |
| [`R/02_preprocess_clr.R`](R/02_preprocess_clr.R) | Prevalence filtering, zero handling, **CLR transform**, train/test export — the compositionally-correct preprocessing. |
| [`python/03_train_rf_shap.py`](python/03_train_rf_shap.py) | Nested, stratified CV; RF vs ElasticNet vs XGBoost; **metformin-stratified** evaluation; AUROC/AUPRC/calibration; cross-cohort transfer test. |
| [`python/04_interpret_shap.py`](python/04_interpret_shap.py) | TreeSHAP global importance, beeswarm, per-patient waterfall, dependence plots. |
| [`docs/references.bib`](docs/references.bib) | BibTeX for every reference in the report. |
| [`data/README.md`](data/README.md) | Data provenance; why no raw data is committed. |
| [`CITATION.cff`](CITATION.cff) | Makes the repo formally citable. |

## The analysis, in one diagram

```
curatedMetagenomicData            preprocessing (R)                   modeling (Python)
┌───────────────────┐   01   ┌──────────────────────┐   02   ┌────────────────────────────┐
│ QinJ_2012 (disc.) │ ─────▶ │ prevalence filter     │ ─────▶ │ nested CV: RF/ENet/XGB      │
│ KarlssonFH_2013   │        │ zero impute → CLR     │        │ metformin-stratified        │
│ (external valid.) │        │ train / test split    │        │ AUROC · AUPRC · calibration │
└───────────────────┘        └──────────────────────┘        │ cross-cohort transfer       │
                                                              └─────────────┬──────────────┘
                                                                            │ 04
                                                                  ┌─────────▼──────────┐
                                                                  │ SHAP: global +     │
                                                                  │ per-patient + dep. │
                                                                  └────────────────────┘
```

## Quickstart

> The code is a **scaffold**: function signatures, the correct pipeline order, and methodological guardrails are in place; some steps are marked `TODO` so a thesis student fills them in deliberately (and learns the method). It is intentionally not a one-click black box.

**1. R — fetch & preprocess** (needs R ≥ 4.3, Bioconductor):
```r
# install.packages("BiocManager")
# BiocManager::install(c("curatedMetagenomicData", "mia", "zCompositions"))
source("R/01_fetch_curatedMetagenomicData.R")   # writes data/processed/*.csv
source("R/02_preprocess_clr.R")                 # writes data/processed/*_clr.csv
```

**2. Python — model & explain** (needs Python ≥ 3.10):
```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python python/03_train_rf_shap.py --cohort QinJ_2012 --external KarlssonFH_2013
python python/04_interpret_shap.py --model artifacts/rf_model.joblib
```

## Scope, honesty, and limitations (read before citing)

- This is a **research/educational artifact**, not clinical guidance.
- The review uses strict citation discipline: claims are tied to a numbered reference or flagged `[INFERENCE]`; figures/DOIs the author could not fully confirm are flagged `[VERIFY]`. **Confirm every `[VERIFY]` against the primary source before reuse.**
- No patient-level or raw sequencing data are stored here; all data are pulled at runtime from the public `curatedMetagenomicData` resource under its own license.
- GLP-1 RA response prediction is presented as **motivation and future work**, not a delivered result — the public data to settle it does not yet exist (see `REPORT.md` §5c, §8 Gap 1).

## How to use this as a thesis foundation

See the **Thesis Integration Note** at the end of [`REPORT.md`](REPORT.md): exact datasets, five concrete comparisons to run, and a month-by-month 6-month plan. Short version: reproduce the metformin-vs-disease disentangling inside a modern RF+SHAP pipeline, quantify the cross-population generalization gap, and end with a fundable prospective GLP-1 RA cohort proposal.

## Citation

If this scaffold or review helps your work, please cite it via [`CITATION.cff`](CITATION.cff).

## License

[MIT](LICENSE) for the code and prose. Underlying datasets and tools (`curatedMetagenomicData`, MetaPhlAn, HUMAnN) retain their own licenses.
