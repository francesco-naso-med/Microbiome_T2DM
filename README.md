# [W.I.P.] Gut Microbiome & Type 2 Diabetes — Mechanisms, Evidence, and a Reproducible ML Roadmap 

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Made with R](https://img.shields.io/badge/R-curatedMetagenomicData-276DC3?logo=r)](R/)
[![Made with Python](https://img.shields.io/badge/Python-RandomForest%20%2B%20SHAP-3776AB?logo=python)](python/)
[![Status](https://img.shields.io/badge/status-research%20artifact-success)](REPORT.md)

> A self-contained, citation-disciplined review of the gut microbiome in Type 2 Diabetes Mellitus (T2DM), **paired with a runnable bioinformatics pipeline** for building an interpretable predictive model from public shotgun-metagenomic data. Designed to double as the theoretical + computational foundation for a Master's/MD thesis.

---

## Why this repository exists

Most microbiome reviews stop at biology. Most ML tutorials stop at toy data. This repo connects the two: a **A mechanistic review** ([`REPORT.md`](REPORT.md)) whose §6 "Analysis Roadmap" maps one-to-one onto **working code skeletons** ([`/R`](R/), [`/python`](python/)) that load `curatedMetagenomicData`, transform it correctly (CLR for compositional data), train Random Forest / ElasticNet / XGBoost under leakage-free repeated cross-validation, test transfer to an independent cohort, and explain predictions with SHAP.

The central scientific thesis it advances: **much of the published "T2DM microbiome signature" is confounded by metformin** (Forslund et al., *Nature* 2015) — and the open, fundable frontier is whether the **baseline microbiome predicts response to GLP-1 receptor agonists** (semaglutide, tirzepatide), for which no public shotgun cohort yet exists.

## Results at a glance (v1.1)

| | Random Forest | XGBoost | ElasticNet |
|---|---|---|---|
| Internal, `QinJ_2012` (170 T2DM vs 174 controls, 10×10 CV) — AUROC | **0.73** ± 0.01 | **0.73** ± 0.01 | 0.61 ± 0.01 |
| External, `KarlssonFH_2013` (53 T2DM vs 43 controls) — AUROC [95% CI] | 0.49 [0.38–0.62] | 0.49 [0.38–0.61] | 0.35 [0.25–0.46] |

- **Within a cohort, the signal is real and biologically coherent:** depletion of butyrate producers (*Roseburia*, *Faecalibacterium*, *Anaerostipes*) and enrichment of opportunistic Clostridia drive the prediction.
- **It is not a metformin artifact in Qin:** metformin-treated cases are *not* easier to classify than untreated or drug-free ones (subgroup AUROCs with bootstrap CIs, `REPORT.md` §6h).
- **It does not transfer across populations:** at chance in Swedish women even after excluding IGT and restricting to drug-free participants (no metformin, no statins), and **significantly inverted** for the linear model. Population-specific signal, not labels or drugs, is the most plausible explanation.

Full tables, subgroup analyses, SHAP directions and limitations: [`REPORT.md` §6h](REPORT.md#6h-results-v11-november-2026).

## What's inside

| Path | Contents |
|------|----------|
| [`REPORT.md`](REPORT.md) | The full review: executive summary, 5 mechanistic axes (SCFA/GLP-1, LPS/TLR4, bile acids/FXR-TGR5, BCAA/mTOR, barrier), a 22-row taxa evidence table (Oxford CEBM grading), a 25-study landmark timeline (2006–2025), a deep GLP-1-RA × microbiome section, the ML roadmap, clinical translation, 7 research gaps, 50+ Vancouver references, and a thesis-integration note. |
| [`R/01_fetch_curatedMetagenomicData.R`](R/01_fetch_curatedMetagenomicData.R) | Pull `QinJ_2012` (discovery) + `KarlssonFH_2013` (external validation) as `TreeSummarizedExperiment`; export species + pathway matrices and metadata. |
| [`R/02_preprocess_clr.R`](R/02_preprocess_clr.R) | *Optional / exploratory.* CLR on the whole cohort for descriptive analyses. **Not used by the ML pipeline**: in `python/03_*` filtering, zero replacement and CLR are re-fitted inside each CV training fold (`python/clr.py`) to avoid leakage. |
| [`python/clr.py`](python/clr.py) | Leakage-free compositional preprocessing as a scikit-learn transformer: prevalence filter, multiplicative zero replacement and CLR, fitted on training folds only. |
| [`python/03_train_rf_shap.py`](python/03_train_rf_shap.py) | T2DM vs control (IGT excluded); repeated stratified CV (nested with `--tune`); RF vs ElasticNet vs XGBoost; **metformin/statin subgroup scoring** with bootstrap 95% CI; AUROC/AUPRC/Brier; cross-cohort transfer test. |
| [`python/04_interpret_shap.py`](python/04_interpret_shap.py) | TreeSHAP global importance with **direction** (rank correlation between abundance and SHAP), beeswarm, per-patient waterfall, dependence plots. |
| [`docs/references.bib`](docs/references.bib) | BibTeX for every reference in the report. |
| [`data/README.md`](data/README.md) | Data provenance; why no raw data is committed. |
| [`CITATION.cff`](CITATION.cff) | Makes the repo formally citable. |

## The analysis, in one diagram

```
curatedMetagenomicData            export (R)                          modeling (Python)
┌───────────────────┐   01   ┌──────────────────────┐   02   ┌────────────────────────────┐
│ QinJ_2012 (disc.) │ ─────▶ │ export raw species +  │ ─────▶ │ repeated CV: RF/ENet/XGB    │
│ KarlssonFH_2013   │        │ metadata (02 = EDA    │        │ CLR re-fitted in each fold  │
│ (external valid.) │        │ only, not in the ML)  │        │ drug subgroups + boot. CI   │
└───────────────────┘        └──────────────────────┘        │ cross-cohort transfer       │
                                                              └─────────────┬──────────────┘
                                                                            │ 04
                                                                  ┌─────────▼──────────┐
                                                                  │ SHAP: global +     │
                                                                  │ per-patient + dep. │
                                                                  └────────────────────┘
```

## Quickstart

> The pipeline runs end to end and reproduces the results in `REPORT.md` §6h. A full 10×10 run takes several minutes on a laptop; add `--folds 5 --repeats 2` for a quick check.

**1. R — fetch & preprocess** (needs R ≥ 4.3, Bioconductor):
```r
# install.packages("BiocManager")
# BiocManager::install(c("curatedMetagenomicData", "mia", "zCompositions"))
source("R/01_fetch_curatedMetagenomicData.R")   # writes data/processed/*.csv
# source("R/02_preprocess_clr.R")               # optional, exploratory only
```

**2. Python — model & explain** (needs Python ≥ 3.10):
```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python python/03_train_rf_shap.py --cohort QinJ_2012 --external KarlssonFH_2013
python python/04_interpret_shap.py --cohort QinJ_2012
```
On macOS, XGBoost needs OpenMP: `brew install libomp` (without it the XGBoost model is skipped, not the run).

## Scope, honesty, and limitations (read before citing)

- This is a **research/educational artifact**, not clinical guidance.
- The review uses strict citation discipline: claims are tied to a numbered reference or flagged `[INFERENCE]`; figures/DOIs the author could not fully confirm are flagged `[VERIFY]`. **Confirm every `[VERIFY]` against the primary source before reuse.**
- No patient-level or raw sequencing data are stored here; all data are pulled at runtime from the public `curatedMetagenomicData` resource under its own license.
- Results come from **two cohorts only**, with sparse medication metadata in `QinJ_2012`, untuned models and no batch correction (`REPORT.md` §6h, Limitations).
- GLP-1 RA response prediction is presented as **motivation and future work**, not a delivered result — the public data to settle it does not yet exist (see `REPORT.md` §5c, §8 Gap 1).

## How to use this as a thesis foundation

See the **Thesis Integration Note** at the end of [`REPORT.md`](REPORT.md): exact datasets, five concrete comparisons to run, and a month-by-month 6-month plan. Short version: reproduce the metformin-vs-disease disentangling inside a modern RF+SHAP pipeline, quantify the cross-population generalization gap, and end with a fundable prospective GLP-1 RA cohort proposal.

## Citation

If this scaffold or review helps your work, please cite it via [`CITATION.cff`](CITATION.cff).

## License

[MIT](LICENSE) for the code and prose. Underlying datasets and tools (`curatedMetagenomicData`, MetaPhlAn, HUMAnN) retain their own licenses.
