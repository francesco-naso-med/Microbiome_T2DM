# Gut microbiome classification of type 2 diabetes: a within-cohort signal that does not cross populations

**Version 2.0 — November 2026** · Francesco Naso · Code: [`R/`](R/), [`python/`](python/) · Reproduce: see [README](README.md#quickstart)

---

## Summary

**Background.** Shotgun metagenomic studies report a gut microbiome signature of type 2 diabetes (T2D), but the two landmark cohorts disagree on its markers [1,2], and part of the signal attributed to T2D has been shown to reflect metformin treatment [3]. Machine-learning studies in this field are also prone to over-optimistic evaluation and poor cross-study transfer [6,7].

**Aim.** To estimate, with a leakage-free pipeline, (i) how well species-level profiles separate T2D from controls within a cohort, (ii) whether this depends on metformin, and (iii) whether the signature transfers to an independent population.

**Methods.** MetaPhlAn species profiles of `QinJ_2012` (China; 170 T2D, 174 controls) and `KarlssonFH_2013` (Sweden, women aged ~70; 53 T2D, 43 controls; impaired glucose tolerance excluded) from `curatedMetagenomicData` [8]. Prevalence filtering, zero replacement and centred log-ratio (CLR) transform were refitted inside every training fold. Random Forest, XGBoost and elastic-net logistic regression were evaluated by 10×10 repeated stratified cross-validation on Qin, then trained on all of Qin and tested on Karlsson. Medication effects were assessed by scoring T2D subgroups (metformin, statins) against controls with bootstrap 95% confidence intervals. Models were interpreted with TreeSHAP.

**Results.** Within Qin, tree ensembles reached AUROC 0.73 (Random Forest 0.732 ± 0.008; XGBoost 0.734 ± 0.012). Metformin-treated cases were not classified better than untreated ones (0.64 [0.48–0.79] vs 0.74 [0.66–0.81]). On Karlsson, discrimination was at chance (AUROC 0.49 [0.38–0.62]), also when both groups were restricted to participants on neither metformin nor statins; the elastic-net signature was significantly inverted (0.35 [0.25–0.46]). The most influential species were butyrate producers (*Roseburia intestinalis*, *Faecalibacterium prausnitzii*, *Anaerostipes hadrus*), lower abundance predicting T2D.

**Conclusions.** The Qin T2D signature is reproducible within its cohort, biologically coherent and not explained by metformin, but it does not transfer to a Swedish cohort, and the failure is not explained by outcome definition or medication. Species-level microbiome classifiers of T2D appear population-specific.

---

## 1. Background

Qin *et al.* described T2D in a Chinese cohort as a moderate dysbiosis with depletion of butyrate-producing bacteria and enrichment of opportunistic pathogens [1]. A year later, Karlsson *et al.* profiled 145 Swedish women with normal, impaired or diabetic glucose control; applying their model to the Chinese cohort, they found that the discriminant markers for T2D differed between the two populations and argued that predictive tools should be specific to age and geography [2].

Forslund *et al.* then re-analysed these cohorts with a Danish one and showed that antidiabetic treatment confounds the picture: metformin has its own strong microbiome signature, while, after controlling for it, a shared T2D signature with depletion of butyrate producers remained [3]. Metformin's microbiome effects were confirmed in a randomised trial in treatment-naive patients [4]. Statins are a second, less often considered confounder [5].

Methodologically, Pasolli *et al.* established Random Forest as a strong default for metagenomic classification [6], and Wirbel *et al.* showed that models naively transferred across studies lose accuracy, while within-study estimates are often inflated by evaluation errors [7].

This project asks three questions on public data: how well T2D can be classified within a cohort when the evaluation is leakage-free; whether the result depends on metformin; and whether it transfers to a different population.

## 2. Methods

### 2.1 Data

Species-level relative abundances (MetaPhlAn 3, bioBakery 3 [9]) and sample metadata were obtained from `curatedMetagenomicData` (release 2021-10-14) [8] with [`R/01_fetch_curatedMetagenomicData.R`](R/01_fetch_curatedMetagenomicData.R). Both cohorts contain one sample per subject.

### 2.2 Outcome and exclusions

The outcome was T2D versus control as recorded in `study_condition`. Excluded were 19 Qin samples with missing condition and the 49 Karlsson participants with impaired glucose tolerance (IGT), which are neither cases nor normoglycaemic controls.

### 2.3 Medication

Metformin and statin use were derived from the `treatment` field as three-state variables (yes / no / not recorded); a missing record was never treated as "untreated". In Qin, treatment is recorded for 145 of 344 samples; in Karlsson for all but one.

### 2.4 Compositional preprocessing

Microbiome profiles are compositional, so raw proportions induce spurious correlations [10]. Each training fold was processed by a single scikit-learn transformer ([`python/clr.py`](python/clr.py)): closure to proportions; retention of species with abundance > 10⁻⁴ in ≥ 10% of training samples; re-closure; multiplicative replacement of zeros [11] with, for each species, 0.65 × its smallest positive training value; CLR. All parameters were learned on the training fold and applied unchanged to test samples and to the external cohort, where species absent from the table were treated as undetected. On the full Qin cohort, 155 species were retained.

### 2.5 Models

Random Forest [13] (1000 trees, √p features per split, minimum leaf 3); XGBoost [14] (600 trees, depth 4, learning rate 0.03, row/column subsampling 0.8/0.5); elastic-net logistic regression [15] (α-mix 0.5, C = 1, standardised inputs). Class imbalance was handled by balanced class weights or `scale_pos_weight`. Hyperparameters were fixed; nested tuning is available (`--tune`) but was not used for the results below. Implementation: scikit-learn [16], XGBoost.

### 2.6 Validation

- **Internal:** 10×10 repeated stratified cross-validation on Qin; metrics averaged across repeats (mean ± SD).
- **Medication subgroups:** out-of-fold probabilities, averaged across repeats, were scored separately for T2D subgroups (metformin yes / no / not recorded; free of both metformin and statins) against controls. Models were not retrained within subgroups, which would be unstable at n = 15. Rationale: if a model detected the drug rather than the disease, treated cases would be the easiest to identify [3].
- **External:** pipeline trained on all of Qin and applied unchanged to Karlsson, overall and by subgroup.
- **Metrics:** AUROC with stratified percentile bootstrap 95% CI (2000 resamples); AUPRC (chance level = prevalence); Brier score. Subgroups with fewer than 10 cases or controls were not scored.

### 2.7 Interpretation

TreeSHAP [12] was computed for the Random Forest trained on all of Qin. Importance is the mean |SHAP| per species; direction is the Spearman correlation between a species' CLR value and its SHAP value (the mean signed SHAP is not a measure of direction, as it averages to approximately zero).

## 3. Results

### 3.1 Cohorts

| | QinJ_2012 | KarlssonFH_2013 |
|---|---|---|
| Population | China, both sexes | Sweden, women aged ~70 |
| T2D / controls analysed | 170 / 174 | 53 / 43 |
| Excluded | 19 (condition missing) | 49 (IGT) |
| T2D on metformin / not / not recorded | 15 / 56 / 99 | 20 / 33 / 0 |
| Statin users (T2D / controls) | not recorded | 26 / 10 |

### 3.2 Internal validation (Qin, 10×10 CV)

| Model | AUROC | AUPRC (chance 0.49) | Brier |
|---|---|---|---|
| Random Forest | 0.732 ± 0.008 | 0.721 ± 0.011 | 0.212 |
| XGBoost | 0.734 ± 0.012 | 0.718 ± 0.016 | 0.227 |
| Elastic net | 0.607 ± 0.012 | 0.571 ± 0.017 | 0.342 |

The two tree ensembles are equivalent. The untuned elastic net discriminates worse and is poorly calibrated.

### 3.3 Medication subgroups (Qin, out-of-fold)

| T2D subgroup vs controls | n | Random Forest | XGBoost |
|---|---|---|---|
| All T2D | 170 vs 174 | 0.736 [0.680–0.788] | 0.741 [0.689–0.793] |
| On metformin | 15 vs 174 | 0.636 [0.479–0.794] | 0.623 [0.448–0.793] |
| Not on metformin | 56 vs 174 | 0.741 [0.662–0.814] | 0.762 [0.687–0.831] |
| Treatment not recorded | 99 vs 174 | 0.748 [0.689–0.809] | 0.748 [0.689–0.806] |
| Drug-free, both groups* | 56 vs 74 | 0.777 [0.691–0.854] | 0.796 [0.714–0.872] |

\*Controls restricted to the 74 with a recorded, drug-free treatment status; not directly comparable with the rows above.

Metformin-treated cases were classified no better, and nominally worse, than untreated cases. With 15 treated cases the interval is wide.

### 3.4 External validation (trained on Qin, tested on Karlsson)

| T2D subgroup vs controls | n | Random Forest | XGBoost | Elastic net |
|---|---|---|---|---|
| All T2D | 53 vs 43 | 0.494 [0.376–0.617] | 0.494 [0.378–0.614] | 0.349 [0.245–0.461] |
| On metformin | 20 vs 43 | 0.540 [0.392–0.693] | 0.535 [0.384–0.688] | 0.329 [0.192–0.484] |
| Not on metformin | 33 vs 43 | 0.467 [0.335–0.600] | 0.469 [0.339–0.606] | 0.361 [0.233–0.495] |
| No metformin, no statins (both groups) | 15 vs 32 | 0.463 [0.275–0.637] | 0.467 [0.285–0.654] | 0.333 [0.175–0.510] |

Three of the 155 model species were absent from the Karlsson profiles. Tree ensembles performed at chance in every subgroup. The elastic-net AUROC was below 0.5 with a confidence interval excluding 0.5: its signature is inverted in this cohort.

### 3.5 Model interpretation (SHAP, Qin)

| Higher abundance → lower predicted risk | Higher abundance → higher predicted risk |
|---|---|
| *Roseburia intestinalis* (ρ = −0.93) | *Erysipelatoclostridium ramosum* (ρ = 0.82) |
| *Anaerostipes hadrus* (−0.80) | *Intestinimonas butyriciproducens* (0.92) |
| *Haemophilus parainfluenzae* (−0.88) | *Gordonibacter pamelaeae* (0.84) |
| *Faecalibacterium prausnitzii* (−0.86) | *Clostridium citroniae* (0.84) |
| *Agathobaculum butyriciproducens* (−0.92) | *Lawsonibacter asaccharolyticus* (0.81) |
| *Roseburia inulinivorans* (−0.79) | *Coprobacillus cateniformis* (0.72) |
| *Ruminococcus callidus* (−0.82) | *Bacteroides dorei* (0.77) |

Species are listed in order of mean |SHAP| within each column. *Lactobacillus mucosae* ranked third overall with a non-monotonic effect.

## 4. Discussion

**A real within-cohort signal.** An AUROC of 0.73 under leakage-free evaluation is in the range reported by earlier machine-learning analyses of this cohort [6]. The model relies mainly on the depletion of butyrate producers, which is the most reproducible T2D-associated change [1,3]. Two species do not fit this picture and should not be over-interpreted: *Intestinimonas butyriciproducens*, a butyrate producer associated here with higher risk, and *Haemophilus parainfluenzae*, associated with lower risk.

**Metformin does not explain the Qin result.** If the classifier had learned a drug signature, treated cases would be the easiest to classify; the opposite trend was observed. This differs from Forslund *et al.*, who found metformin-treated T2D more separable from controls than untreated T2D [3]. Our treated group is small (n = 15) and the comparison uses species profiles and a different classifier, so the two results are not directly in conflict. The proportion of treated cases among Qin patients with a recorded treatment (15/71, 21%) matches the 21% reported for the Chinese cohort by Forslund *et al.* [3].

**The signature does not transfer.** Transfer to Karlsson failed even after removing IGT from the controls and restricting both groups to participants on neither metformin nor statins. Outcome definition and medication therefore do not explain the failure. The significantly inverted elastic-net signature indicates that some species are associated with T2D in opposite directions in the two cohorts. This agrees with Karlsson *et al.*'s original observation [2] and with the general loss of accuracy of naively transferred microbiome models [7]. Remaining candidate explanations are population differences (ancestry, sex, age, diet) and technical differences between studies (sample handling, DNA extraction, sequencing), which these data cannot separate.

**Implication.** A species-level T2D classifier trained on one population should not be expected to work on another. Claims of a "microbiome test" for T2D require validation across several independent populations.

## 5. Limitations

- Two cohorts only; external validation in one direction (Qin → Karlsson).
- Treatment not recorded for 58% of Qin samples; small treated subgroups produce wide intervals.
- Fixed hyperparameters; no correction for batch, BMI, age or sex.
- Species-level taxonomy only (no functional pathways or strains).
- SHAP describes the model, not causal biology.

## 6. Next steps

1. **Leave-one-cohort-out validation** with the other T2D cohorts in `curatedMetagenomicData`.
2. **Identify the inverted species** by compositional differential abundance within each cohort (e.g. ANCOM-BC [18]).
3. **Batch harmonisation** (MMUPHin [17]), noting that it may remove genuine population differences.
4. **Nested tuning**, mainly to improve elastic-net calibration.
5. **Functional profiles** (HUMAnN pathways), which may transfer better than taxonomy.

A longer-term question is whether the baseline microbiome predicts response to glucose-lowering drugs such as GLP-1 receptor agonists. No public cohort with pre- and post-treatment metagenomes and response phenotypes is currently available in `curatedMetagenomicData`, so this is outside the scope of the present analysis.

## Version history

- **2.0 (November 2026).** Rewritten as a focused analysis report. The previous broad literature review (mechanisms, GLP-1, clinical translation) was removed: it contained unverified references and was not supported by the analysis. All references verified.
- **Code corrections behind these results.** Compositional preprocessing moved inside cross-validation folds (previously fitted on the whole cohort); cross-validation made genuinely repeated; external CLR computed on the discovery species set (previously, separately computed CLR matrices were intersected); IGT and unlabelled samples excluded (previously counted as controls); missing treatment no longer treated as untreated; SHAP direction computed as a rank correlation.

## References

1. Qin J, Li Y, Cai Z, et al. A metagenome-wide association study of gut microbiota in type 2 diabetes. *Nature*. 2012;490(7418):55-60. doi:10.1038/nature11450
2. Karlsson FH, Tremaroli V, Nookaew I, et al. Gut metagenome in European women with normal, impaired and diabetic glucose control. *Nature*. 2013;498(7452):99-103. doi:10.1038/nature12198
3. Forslund K, Hildebrand F, Nielsen T, et al. Disentangling type 2 diabetes and metformin treatment signatures in the human gut microbiota. *Nature*. 2015;528(7581):262-266. doi:10.1038/nature15766
4. Wu H, Esteve E, Tremaroli V, et al. Metformin alters the gut microbiome of individuals with treatment-naive type 2 diabetes, contributing to the therapeutic effects of the drug. *Nat Med*. 2017;23(7):850-858. doi:10.1038/nm.4345
5. Vieira-Silva S, Falony G, Belda E, et al. Statin therapy is associated with lower prevalence of gut microbiota dysbiosis. *Nature*. 2020;581(7808):310-315. doi:10.1038/s41586-020-2269-x
6. Pasolli E, Truong DT, Malik F, Waldron L, Segata N. Machine learning meta-analysis of large metagenomic datasets: tools and biological insights. *PLoS Comput Biol*. 2016;12(7):e1004977. doi:10.1371/journal.pcbi.1004977
7. Wirbel J, Zych K, Essex M, et al. Microbiome meta-analysis and cross-disease comparison enabled by the SIAMCAT machine learning toolbox. *Genome Biol*. 2021;22(1):93. doi:10.1186/s13059-021-02306-1
8. Pasolli E, Schiffer L, Manghi P, et al. Accessible, curated metagenomic data through ExperimentHub. *Nat Methods*. 2017;14(11):1023-1024. doi:10.1038/nmeth.4468
9. Beghini F, McIver LJ, Blanco-Míguez A, et al. Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. *eLife*. 2021;10:e65088. doi:10.7554/eLife.65088
10. Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. Microbiome datasets are compositional: and this is not optional. *Front Microbiol*. 2017;8:2224. doi:10.3389/fmicb.2017.02224
11. Martín-Fernández JA, Barceló-Vidal C, Pawlowsky-Glahn V. Dealing with zeros and missing values in compositional data sets using nonparametric imputation. *Math Geol*. 2003;35(3):253-278. doi:10.1023/A:1023866030544
12. Lundberg SM, Erion G, Chen H, et al. From local explanations to global understanding with explainable AI for trees. *Nat Mach Intell*. 2020;2(1):56-67. doi:10.1038/s42256-019-0138-9
13. Breiman L. Random forests. *Mach Learn*. 2001;45(1):5-32. doi:10.1023/A:1010933404324
14. Chen T, Guestrin C. XGBoost: a scalable tree boosting system. In: *Proc 22nd ACM SIGKDD Int Conf Knowl Discov Data Min*. 2016:785-794. doi:10.1145/2939672.2939785
15. Zou H, Hastie T. Regularization and variable selection via the elastic net. *J R Stat Soc Series B Stat Methodol*. 2005;67(2):301-320. doi:10.1111/j.1467-9868.2005.00503.x
16. Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn: machine learning in Python. *J Mach Learn Res*. 2011;12:2825-2830.
17. Ma S, Shungin D, Mallick H, et al. Population structure discovery in meta-analyzed microbial communities and inflammatory bowel disease using MMUPHin. *Genome Biol*. 2022;23(1):208. doi:10.1186/s13059-022-02753-4
18. Lin H, Peddada SD. Analysis of compositions of microbiomes with bias correction. *Nat Commun*. 2020;11(1):3514. doi:10.1038/s41467-020-17041-7
