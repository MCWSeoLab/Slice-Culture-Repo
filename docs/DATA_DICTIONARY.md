# Data dictionary

All files are UTF-8 CSV with a header row, in `data/`.

Conventions used throughout:

- `study_id` — patient identifier, `PT-001` to `PT-049`, stable across every file.
- `sample` / `sample_label` / `biopsy` — laboratory specimen label, e.g. `Liver Bx 9`.
- `day`, `study_day`, `day_start`, `day_stop`, `day_progression` — integer days relative to
  that patient's own biopsy. Negative values precede the biopsy. No calendar date appears in
  any file.
- `pct_of_control` — treated absorbance as a percentage of the DMSO control on the same plate
  and timepoint. Percent inhibition is `100 - pct_of_control`.
- `[DATE]`, `[NAME]`, `[ID]` — markers left by the de-identification scrubber recording that
  an identifier was removed at that position.

## Clinical

| file | rows | contents |
|---|---|---|
| `table1_source_deid.csv` | 19 | Cohort characteristics per specimen: `sample`, `group`, `age_years`, `sex`, `stage`, `grade`, `approach`, `prior_systemic`. Serum tumor marker columns are withheld. |
| `table2_source_deid.csv` | 19 | Per-biopsy characteristics: `sample`, `study_day`, `age_years`, `sex`, `histologic_type`, `genomic`, `met_site`. |
| `genomics_deid.csv` | 53 | Tidy form of the `genomic` column above: one row per specimen per reported gene, with `gene` and `variant`. Variants are recorded verbatim as the study worksheet reports them; `wild-type` is recorded where the worksheet states it, and an empty `variant` means the worksheet named the gene without a protein-level change. |
| `per_sheet_annotations_deid.csv` | 94 | The laboratory's own per-worksheet response call (`lab_resp2`) against the recorded clinical response (`clinic_resp`), keyed by `study_id`, `sheet`, `row`, `which_matched`. Used only to document what the un-normalized laboratory metric yields; it is not the primary endpoint. |
| `pairs_final.csv` | 13 | The adjudicated analysis set: one row per matched drug–patient pair, with ex vivo `mts_pct_of_control`, the matched `regimen`, its `day`/`day_stop`/`day_progression`, best `response`, the `flag` marking tissue-exhausted rows held out of every statistic, and the derived `sens`, `benefit` and `rank`. |
| `pairs_recist.csv` | 11 | The same pairs carrying the RECIST ordinal (`recist_rank`) and the drug-class match. |
| `pairs_pfs.csv` | 11 | The same pairs extended with the time-to-next-therapy surrogate: `start`, `end`, `ttnt`, `event`, `time`, `started_after_biopsy`. |

## Ex vivo assay

| file | rows | contents |
|---|---|---|
| `mts_wells_deid.csv` | 311 | Well-level MTS readings: `sample`, `timepoint`, `drug`, `concentration`, `well`, `normalised` (blank-corrected absorbance). The primary assay source. |
| `mts_treated_pct_deid.csv` | 224 | Treated wells only, with the matched `control_mean`, `pct_of_control` and the `qc_fail` flag. |
| `mts_dmso_cv_wells_deid.csv` | 26 | DMSO replicate coefficient of variation per experiment, with the pass/fail verdict at the 20% gate. |
| `mts_normalized_final.csv` | 28 | One row per patient × drug class × timepoint, aggregated from the wells. |
| `prism_raw_long_deid.csv` | 2191 | Every value extracted from the GraphPad Prism files, one row per file × sheet × row label × column × replicate. |
| `prism_mts_normalized_deid.csv` | 830 | The MTS subset of the above, normalized to the control on the same sheet and timepoint. |
| `prism_index_deid.csv` | 46 | Inventory of the Prism files and sheets found, with value counts. |
| `qc_dmso_cv.csv` | 14 | QC verdicts as computed by `qc_gates.py` from the Prism route. |
| `qc_by_patient.csv` | 4 | Patient-level QC verdict used to gate the primary analysis. |

## Histology and imaging quantification

| file | rows | contents |
|---|---|---|
| `masterlist_linked_deid.csv` | 34 | One row per slice: `sample_label`, `tissue_type`, `biopsy_site`, `treatment`, `kras`, `dose`, `timepoint_hr`, and the four paired endpoints `mts_pct_of_control`, `nuclear_density_pct`, `caspase3_pct`, `ki67_pct`. Source for Figure 2. |
| `ihc_tidy.csv` | 652 | Immunohistochemistry scoring, one row per field: `marker`, `biopsy`, `treatment`, `replicate`, `value`. |
| `erk_tidy.csv` | 32 | ERK immunohistochemistry by `allele`, `treatment`, `marker` (total and phosphorylated) and `replicate`. |
| `tumor_content_deid.csv` | 10 | EpCAM-positive fraction per biopsy: `epcam_pct_mean`, `epcam_pct_sd`, `n_slices`. |
| `tumour_content_deid.csv` | 10 | Identical content under the spelling written by `make_figS_cellularity.py`, which `adjust_for_cellularity.py` reads. Both are kept so neither script needs editing. |
| `fig1_timecourse_deid.csv` | 27 | Slice surface area and MTS over seven days, one row per `panel` × `biopsy` × `day`. |

## Supplementary series

| file | rows | contents |
|---|---|---|
| `suppfig1_pdac_mts_deid.csv` | 8 | Extended-culture MTS, one specimen, by day. |
| `suppfig1_pdac_nuclear_density_deid.csv` | 30 | Nuclear density for the same series, six fields per day. |
| `suppfig1_pdac_ihc_deid.csv` | 80 | Immunohistochemistry for the same series, by `marker` and `group`. |
| `suppfig2_peritoneal_mts_deid.csv` | 40 | Peritoneal specimen MTS by timepoint and agent, with `pct_of_dmso`. |
| `suppfig2_peritoneal_nuclear_density_deid.csv` | 30 | Nuclear density for the peritoneal specimen. |
| `suppfig2_peritoneal_ihc_deid.csv` | 80 | Immunohistochemistry for the peritoneal specimen. |

## Derived analysis outputs

These are written by scripts in this repository and are included so a reader can compare a
fresh run against them.

| file | rows | contents |
|---|---|---|
| `sensitivity_grid.csv` | 48 | Every combination of cutoff, exclusion rule and response definition, with the resulting 2×2, accuracy, Fisher p and Spearman rho. Written by `analysis/sensitivity.py`. |
| `forest_rho.csv` | 6 | Spearman rho and cluster-bootstrap interval for each analysis set. Written by `figures/make_fig_forest.py`. |
| `table3_assaylog_source_deid.csv` | 19 | Assay log per specimen: `assay` and `details` by `study_day`. |
