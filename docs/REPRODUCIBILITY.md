# What reproduces from this repository

`python run_all.py` runs fifteen stages from `data/` alone. All fifteen complete, and
`python code/checks/verify_results.py` confirms 26 reported values against the manuscript.

Every resampling routine draws from its own generator seeded with the fixed seed in
`.phi-allowlist`, so repeated runs are identical. An earlier version shared one module-level
generator across routines, which made a confidence interval depend on the order in which
panels were drawn; each routine now constructs its own.

## Reported values checked

| value | reported | regenerated |
|---|---|---|
| Figure 2A nuclear density, Spearman rho | +0.743 | +0.7430 |
| Figure 2A exact p, biopsies permuted | 0.001 | 0.0007 |
| Figure 2A n | 33 slices, 13 biopsies | 33, 13 |
| Figure 2B Ki-67, rho / exact p | +0.651 / 0.015 | +0.6507 / 0.0153 |
| Figure 2C cleaved caspase-3, rho / exact p | −0.548 / 0.020 | −0.5484 / 0.0203 |
| Figure 3A exact p, outcome permuted across patients | 0.20 | 0.2000 |
| Figure 3B rho, 95% CI, exact p | −0.88, −1.00 to −0.51, 0.06 | −0.8787, −1.0000 to −0.5143, 0.0556 |
| Figure 3B n | 9 courses, 6 patients | 9, 6 |
| 2×2 at the 30% cutoff | TP 4, FP 1, FN 1, TN 4 | 4, 1, 1, 4 |
| 2×2 sensitivity, specificity, PPV, NPV, accuracy | 80% each | 80% each |
| 2×2 Fisher exact p | 0.206 | 0.206 |
| Threshold, DMSO CV route | 32% inhibition | 31.5% |
| Threshold, half-normal route (sigma) | 34% (17.0 pp) | 34.1% (17.05 pp) |
| Negative absorbance readings | 2.1% | 2.1% |

The adopted cutoff is >30% inhibition, slightly more conservative than either derivation.
The classification is unchanged for any cutoff between 27.7% and 39.3% inhibition, because no
observed viability falls between 60.7% and 72.3% of control.

## Stages

| stage | script | output |
|---|---|---|
| `qc_gates` | `analysis/qc_gates.py` | DMSO replicate CV per experiment, negative readings, replication depth |
| `threshold_figure` | `figures/make_threshold_figure.py` | `figS_threshold_derivation` |
| `clustered_p` | `analysis/clustered_p.py` | the patient-level tests behind Figure 3A |
| `recompute` | `analysis/recompute.py` | all three analyses under patient-level clustering |
| `sensitivity` | `analysis/sensitivity.py` | `sensitivity_grid.csv` |
| `cellularity` | `analysis/adjust_for_cellularity.py` | tumor-content adjustment |
| `fig2_mts_vs_histology` | `figures/make_fig2_v2.py` | `fig2_mts_vs_histology_v2` |
| `fig3_concordance` | `figures/make_fig3_v5.py` | `fig3_concordance_response_row` |
| `fig3_confusion` | `figures/make_confusion_response.py` | `fig3_confusion_response` |
| `figS_tumor_content` | `figures/make_figS_cellularity.py` | `figS_tumour_content` |
| `figS_replicate_agreement` | `figures/make_fig_reproducibility.py` | `figS_replicate_agreement` |
| `fig_flow_diagram` | `figures/make_fig_flow.py` | `fig_flow_diagram` |
| `sap_forest` | `figures/make_fig_forest.py` | `fig_rho_forest` |
| `sap_exvivo_vs_recist` | `figures/make_fig_exvivo_vs_recist.py` | `fig_exvivo_vs_recist` |
| `sap_exvivo_vs_ttnt` | `figures/make_fig_exvivo_vs_ttnt.py` | `fig_exvivo_vs_ttnt` |

`make_fig3_v5.py` is run with `FIG3_OUTCOME=response FIG3_LAYOUT=row FIG3_STATS=computed
FIG3_TICKS=short`, which is the variant in the manuscript. `run_all.py` sets these.

# What does not reproduce here

`code/restricted/` holds the stages that need inputs this repository does not contain. They
are included as the provenance of the released tables and figures, not as a runnable path.

| script | needs |
|---|---|
| `deidentify.py` | the identified source workbooks; writes the re-identification crosswalk |
| `deid_chart_abstraction.py` | the identified chart abstraction workbook and the crosswalk |
| `extract_prism.py` | the GraphPad Prism files and the crosswalk |
| `extract_mts_wells.py` | the two identified MTS workbooks |
| `build_supp_data.py` | the Prism files for the supplementary series |
| `build_release_dataset.py` | the restricted de-identified directory |
| `master_sheet_data.py` | the crosswalk and the withheld clinical tables |
| `build_regimens_merged.py` | the withheld therapy log |
| `build_course_csvs.py` | the withheld therapy, imaging and marker tables |
| `recompute_with_chart.py` | the withheld therapy log |
| `adjudicate_and_recompute.py` | the output of `recompute_with_chart.py`; produces `pairs_final.csv` |
| `derive_threshold_local.py` | the identified master workbook |
| `validate_table1.py` | the serum marker columns withheld from `table1_source_deid.csv` |
| `make_fig1_v4.py` | the Figure 1 micrographs |
| `make_suppfig1_v4.py` | the supplementary-series micrographs |
| `make_fig4_cases.py` | the Figure 4 micrographs and the withheld therapy and marker tables |
| `make_fig4_panels_big.py` | the enlarged Figure 4 panel tiles |
| `build_fig4_v6.py` | the Figure 4 page PDF, the enlarged panel PDFs and `qpdf` |
| `make_fig_erk_cases_wide.py` | the immunohistochemistry and ERK micrographs, and the withheld tables |
| `make_fig_sts.py` | the identified MTS workbook |

`lib/course_panel.py` draws the clinical-course strip and reads the withheld therapy,
imaging and marker tables. It is imported only by the restricted figure scripts.

Two specific gaps follow from the release policy:

- **Figure 4 panels D and H**, the clinical-course plots, need the serial tumor marker and
  scan-level RECIST series. Both are withheld. The submitted figure is in `figures/`.
- **The serum marker rows of Table 1** need the withheld CA19-9, CEA and CA125 columns.

Micrographs are not included in this repository. The quantitative panels of every figure are
rebuilt from `data/`; the photomicrograph panels are not. The figure files as submitted are
in `figures/`, as PDF throughout and as PNG and SVG where the file size allows.

## Analysis conventions

These are stated in `docs/statistical_analysis_plan.md` and are repeated here because they
determine whether a reimplementation will agree.

- The unit of analysis is the drug–patient pair. Both the number of pairs and the number of
  patients are reported for every statistic.
- Inference is by patient-level cluster bootstrap with 4,000 resamples, together with
  enumeration over one pair per patient. Where a permutation p-value is reported, whole
  clusters of outcomes are permuted between clusters of equal size and every such permutation
  is enumerated; the cluster is the patient for Figure 3 and the biopsy for Figure 2.
  Enumerate within each size class and take the product, rather than over all orderings.
- Drugs are matched at class level. FOLFOX and FOLFIRI are not matched to FOLFIRINOX.
- The primary cutoff of >30% inhibition is data-derived, not pre-specified, and is described
  that way in the report.
- The switch of the primary clinical outcome to binary clinical benefit is a post-hoc
  redefinition of the frozen plan and is reported as one.
- Rows flagged tissue-exhausted are drawn in the figures but held out of every statistic.
