# Statistical Analysis Plan

**Study:** Ex vivo organotypic slice culture of core biopsies in pancreatic ductal adenocarcinoma — concordance between ex vivo drug response and clinical outcome

**Target journal:** Clinical Cancer Research
**Version:** 1.6
**Date frozen:** 18 August 2026 (v1.0); amended 19 August 2026 (v1.1, v1.2) and
20 August 2026 (v1.3, v1.4, v1.5, v1.6 — see §13)
**Analysis data:** the released tables in `data/` (see `code/restricted/deidentify.py` and `code/restricted/extract_prism.py` for their provenance)

---

## 1. Purpose and honest status of this document

This SAP is written **after** exploratory analyses have already been run on this dataset. It
does not and cannot claim pre-registration. Its purpose is to fix the analysis rules from
this point forward, and to make explicit which decisions were data-derived so that readers
and reviewers can discount them appropriately.

Everything in §4 (the response threshold) was derived from the analysis data. It is labeled
data-derived throughout and must be described that way in the manuscript.

---

## 2. Study design and its principal limitation

Retrospective single-institution cohort. Image-guided core biopsies of metastatic and primary
PDAC were cultured as organotypic slices and exposed to drug panels; clinical treatment and
response were abstracted from the medical record.

**Principal design limitation, stated up front:** of the matched drug–patient pairs available
for analysis, only **3 (from 2 patients)** involve a drug started at or after the biopsy. For
the remainder the tissue was sampled after the patient had already received that agent. The
analysis is therefore an assessment of *association*, not of *prediction*. No predictive
claim will be made.

---

## 3. Unit of analysis and clustering

The unit of analysis is the **drug–patient pair**. Patients contribute between 1 and 3 pairs.

Pairs within a patient are not independent. All inference accounts for this:
- **Primary inference:** patient-level cluster bootstrap 95% CI (patients resampled with
  replacement, 4,000 resamples, seed 20260818).
- **Supporting:** exact enumeration of all one-pair-per-patient subsets, reporting the median
  and range of the statistic across subsets.
- **Binary endpoints:** GEE with exchangeable working correlation and patient as the cluster.
- Permutation p-values that ignore clustering are reported only as clearly-labeled
  anticonservative references.

n(pairs) and n(patients) are reported separately everywhere, in every table and figure.

---

## 4. Exposure: the ex vivo response variable

### 4.1 Measurement
Viability is MTS absorbance at 48 hours expressed as a **percentage of the mean DMSO control
on the same plate and timepoint**. 24 hours is a secondary timepoint.

The `Lab Response?` field in the source workbook is **not used**. It thresholds the raw
24h→48h absorbance change with no control normalization, which measures culture decay
confounded by tissue quantity rather than drug effect. Analyses using it are reported only to
document why it was abandoned.

### 4.2 Response threshold — data-derived

**Amended 20 August 2026 (v1.4). The half-normal derivation withdrawn in v1.3 is REINSTATED**
— the well-level source was located and both derivations reproduce. See §13.6.

Two independent derivations, both recomputed from well-level MTS readings
(`scripts/extract_mts_wells.py`, `deid/mts_wells_deid.csv`):

| method | n | estimate | cutoff |
|---|---|---|---|
| Half-normal fit to treated wells reading >100% of control (no agent can truly raise viability, so that tail carries no drug effect by construction) | 12 wells above 100%, of 70 in QC-passing experiments | σ = **17.0 pp** (bootstrap 95% CI 11.1–22.0) | **>34% inhibition** |
| Propagated DMSO replicate CV (median **11.1%** across the 6 QC-passing experiments), single treated well vs single control | 23 experiments with replicate control wells, 6 passing | SD = **15.8 pp** | **>32% inhibition** |

The two agree to within 2 percentage points from unrelated parts of the data. This is reported
as supporting evidence, not as independent validation: both come from the same experiments.

**Primary cutoff: >30% inhibition** (viability <70% of control), unchanged, and slightly more
conservative than either derivation. 20%, 40% and 50% are sensitivity analyses.

**Cutoff insensitivity is reported alongside the derivations.** A call changes only when the
cutoff crosses an observed viability, and no observed value falls between 60.7% and 72.3% of
control, so the classification is **identical for any cutoff from 27.7% to 39.3% inhibition**
— an 11.6 pp window containing both derived estimates and the adopted 30%. A cutoff below
27.7% would classify one further course correctly, so 30% is **not** the value that maximizes
agreement with the outcome, which is the expected consequence of deriving it from assay noise.

**Control-well variability is high and must be reported.** Of the 23 experiments with replicate
DMSO wells, **17 fail the 20% gate in §4.3**; DMSO CVs run from 0.3% to 83%. The median of
11.1% is therefore the median among *passing* experiments, not among all of them, and the
threshold is conditional on that gate.

### 4.3 Quality-control gates

**Scope corrected 20 August 2026 (v1.6).** Earlier versions of this section were headed
"applied before any analysis". That was not what was done, and could not be — see §13.8 and
the table at the end of this section. **Gates 1 and 2 are applied to the noise-floor
estimation in §4.2 and are reported as a sensitivity analysis for the concordance endpoints.
Gates 3 and 4 are applied everywhere, including the primary analysis.**
1. **Exclude any experiment whose DMSO replicate CV exceeds 20%.**

   **Provenance, stated plainly (amended 20 August 2026, v1.5).** The 20% figure was
   **chosen, not derived.** When the gate was written only five experiments had assessable
   replicate control wells: four clustered near a median of 10.8% and one (PT-014) ran at
   72–83% and was producing a 274%-of-control reading. A round-number rule at 20% removed
   that experiment on a stated criterion rather than by eye. It was not calibrated against
   any distribution, and the manuscript must not imply that it was.

   **Post hoc justification on the full data.** With well-level readings now available
   (§13.6), 23 experiments have replicate DMSO wells. Their CVs are:

   > 0.3, 8.4, 10.8, 11.5, 11.7, 14.0 │ 23.8, 25.6, 26.4, 26.7, 26.7, 26.8, 28.2, 28.7,
   > 29.3, 30.8, 31.5, 36.3, 59.2, 63.0, 70.2, 72.4, 83.2  (%)

   There is a **9.8-percentage-point empty gap between 14.0% and 23.8%**, and the 20% gate
   falls inside it. The gate therefore separates a tight cluster from a broad one rather than
   cutting through a dense region. That is a defence of where it lands, not of how it was
   chosen.

   **Sensitivity of the threshold to the gate.** The gate does substantial work and this must
   be reported:

   | gate | experiments kept | median CV | CV-route cutoff | half-normal cutoff |
   |---|---|---|---|---|
   | 15% | 6 of 23 | 11.1% | 31.5% | 35.6% |
   | **20% (adopted)** | **6 of 23** | **11.1%** | **31.5%** | **35.6%** |
   | 25% | 7 of 23 | 11.5% | 32.4% | 34.1% |
   | 30% | 15 of 23 | 25.6% | **72.3%** | 41.9% |
   | 40% | 18 of 23 | 26.5% | 75.1% | 41.8% |
   | none | 23 of 23 | 26.8% | 75.7% | 116.3% |

   Anywhere from 15% to 25% gives materially the same answer. At 30% and beyond the
   CV-derived cutoff more than doubles and the derivation collapses. **The stability of the
   threshold rests on the 14–24% gap being real**, and if a future dataset fills that gap the
   gate must be re-examined rather than carried forward.

   **17 of 23 experiments fail this gate.** Control-well variability in this assay is high,
   and every threshold estimate in §4.2 is conditional on the exclusion.
2. Conditions with a **single well** are flagged and analyzed separately; 24% of conditions
   are single-well.
3. **Negative absorbance readings** (2.1% of readings) reflect background over-subtraction and
   are floored at zero before normalization.
4. Specimens annotated "tissue exhausted", "no tissue" or "unreliable" are excluded from
   primary analyses and shown separately.

Experiments for which replicate-level control data does not exist are labeled QC-unknown and
are not silently treated as passing.

#### Where gate 1 actually falls on the matched pairs

Each of the 10 analyzed pairs, against its biopsy's 48-hour DMSO control CV:

| study id | biopsy | 48 h DMSO CV | gate 20% | gate 30% |
|---|---|---|---|---|
| PT-025 | Liver Bx 10 | 11.7% | pass | pass |
| PT-013 | Peritoneum Bx 1 | single control well | QC-unknown | QC-unknown |
| PT-015 | Liver Bx 7 | single control well | QC-unknown | QC-unknown |
| PT-024 | Liver Bx 9 | single control well | QC-unknown | QC-unknown |
| PT-011 | PDAC Bx 5 | 36.3% | **fail** | **fail** |
| PT-014 | Liver Bx 6 | 72.4% | **fail** | **fail** |

- Applying gate 1 **strictly** (QC-unknown excluded) leaves **1 pair from 1 patient** (PT-025).
- Applying it **leniently** (QC-unknown retained) leaves **7 pairs from 4 patients**
  (PT-013, PT-015, PT-024, PT-025).
- **3 pairs from 2 patients fail** (PT-011, PT-014).

Neither surviving set supports the concordance endpoints in §6, which is why gate 1 is a
sensitivity analysis there rather than a filter. **Raising the gate to 30% changes none of
this** — the pairs sit at the extremes of the CV distribution, not near the boundary — while
it doubles the §4.2 noise floor and pushes the derived cutoff to 72% inhibition, at which
every matched course is called resistant.

**The primary concordance result therefore includes PT-014, which fails gate 1 at 72–83%.**
This must be stated in the Limitations, not left implicit: the experiment that produced the
274%-of-control reading contributes 2 of the 10 primary pairs.

---

## 5. Outcome: the clinical response variable

Primary clinical outcome is **binary clinical benefit** derived from RECIST best response to
the matched agent: **CR, PR or SD = clinical benefit (response); PD = no response.**

*This definition replaced the v1.0 definition on 19 August 2026 and is a post-hoc change.*
Version 1.0 specified the RECIST ordinal (PD < SD < PR < CR) as primary, with responder for
binary analyses defined as CR or PR only. The change was made after the v1.0 results were
known and must be described in the manuscript as a post-hoc redefinition, not as
pre-specified. Full rationale and both sets of results are in §13.

The **RECIST ordinal analysis is retained as a co-primary** and reported alongside the binary
analysis in every table and figure, because the ordinal uses the SD/PD gradation that the
dichotomy discards. Where the two disagree, both are shown.

Time to next systemic therapy (TTNT) is a **secondary, exploratory** outcome. It is explicitly
not PFS. True PFS and OS are not computable: the dataset contains no progression dates, no
death dates, and a last-follow-up date for only 7 of 49 patients.

TTNT is defined as the interval from the start of the matched systemic line to the start of
the next systemic line (event), censored at the recorded end of that line where no subsequent
line exists. Surgery, SBRT and chemoradiation are not counted as lines.

Matching of ex vivo agent to clinical agent is at **drug-class level** (mFOLFIRINOX counts as
FOLFIRINOX; gemcitabine/abraxane as gem/nab-paclitaxel). FOLFOX and FOLFIRI are **not**
matched to FOLFIRINOX. One class-level pair (trametinib ex vivo vs cobimetinib clinically) is
marked distinctly in all figures and removed in sensitivity analysis. Where an agent was given
more than once, the course closest in time to the biopsy is used.

---

## 6. Endpoints

**Primary (v1.1):** association between control-normalized 48h ex vivo viability and **binary
clinical benefit** (CR/PR/SD vs PD) for the same agent, expressed as (a) a comparison of
viability between benefit and no-benefit pairs and (b) a 2×2 concordance table at the >30%
inhibition cutoff with sensitivity, specificity, PPV, NPV and exact Clopper–Pearson CIs.

**Co-primary (carried forward from v1.0):** the same association against RECIST best response
treated as ordinal, as Spearman's rho with a patient-level cluster bootstrap 95% CI, in pairs
passing QC with tissue-exhausted specimens excluded.

**Exploratory:** TTNT association; the 24h timepoint; analyses restricted to drugs started at
or after the biopsy.

**Reported for transparency only:** any analysis using the workbook's `Lab Response?` field.

Only the primary endpoint is confirmatory. No multiplicity adjustment is applied because no
analysis is confirmatory; instead the full sensitivity grid (§7) is reported so that the
stability of the conclusion is visible rather than asserted.

---

## 7. Sensitivity analyses

A full factorial grid is reported: cutoff (30 / 32 / 50% inhibition) × tissue-exhausted
(in / out) × QC gate (on / off) × class-level match (in / out) × Partial counted as response
(yes / no). 48 scenarios. The complete grid is a supplementary table; no scenario is
selectively reported.

---

## 8. Precision, not power

No post hoc power calculation will be performed. Precision is reported directly: at these
sample sizes a specificity point estimate of 25% carries a 95% CI of 5.5–57.2%, which is
compatible with almost any underlying value. Confidence interval width, not p-values, is the
honest summary of what this dataset can establish.

---

## 9. Reporting standard

**REMARK** (REporting recommendations for tumor MARKer prognostic studies) is followed, with
the checklist supplied as supplementary material. STARD elements are used where the 2×2
concordance table is presented. Analyses conducted after this SAP was frozen will be labeled
post hoc.

---

## 10. Results obtained under this plan

Run 18 August 2026. Seed 20260818, 4,000 bootstrap and permutation resamples.

### 10.1 Primary endpoint

| analysis set | pairs / patients | Spearman rho | cluster bootstrap 95% CI | CI excludes 0 |
|---|---|---|---|---|
| All pairs | 11 / 6 | −0.43 | −0.86 to +0.30 | no |
| Tissue-exhausted excluded | 8 / 5 | −0.62 | −0.95 to +0.30 | no |
| **+ QC gate (primary)** | **6 / 4** | **−0.62** | **−1.00 to +1.00** | **no** |

One-pair-per-patient enumeration, tissue-exhausted excluded: median rho −0.58 (range −0.79 to
0.00), median p = 0.31.

**The primary endpoint is negative.** The direction is consistent (rho negative in every
analysis set, meaning lower ex vivo viability with better clinical response), but the
confidence interval includes zero in all of them. Applying the QC gate leaves 6 pairs from 4
patients, at which point the interval spans the entire range of the statistic.

### 10.2 Secondary endpoint — concordance table

The 2×2 analysis is **uninformative by construction**: only **2 of 11 matched pairs had a
clinical response (both PR)**. Every scenario in the grid therefore contains 1 true positive
and 1 false negative. Accuracy ranges 50–78% across the 48 scenarios and **no scenario reaches
p < 0.05**. Cutoffs of 30% and 32% give identical results because only one observation falls
between 50% and 70% of control.

The ordinal analysis is more informative than the binary one here, precisely because it uses
the SD/PD gradation that the responder/non-responder dichotomy discards.

### 10.3 Analyses using the workbook's `Lab Response?` field

31 pairs / 21 patients. TP 11, FP 9, FN 8, TN 3. Naive Fisher p = 0.452. GEE with patient
clustering: OR 0.43 (95% CI 0.10–1.80), p = 0.248. No association, in the expected absence
given §4.1.

### 10.4 Exploratory — TTNT

| analysis set | pairs / patients | rho | cluster bootstrap 95% CI |
|---|---|---|---|
| All pairs | 11 / 6 | −0.43 | −1.00 to +0.20 |
| Tissue-exhausted excluded | 8 / 5 | −0.87 | −1.00 to −0.45 |
| + QC gate | 6 / 4 | −1.00 | −1.00 to −1.00 |

**This is the only interval that excludes zero, and it should be treated with more suspicion
than any other result in the study, not less.** Three reasons:

1. TTNT is not independent of RECIST response — responders stay on treatment longer — so this
   is largely the primary endpoint re-expressed, not corroboration.
2. Five of the eight courses began before the biopsy, one of them 1,805 days before.
3. After the QC gate the estimate is exactly −1.00 across 6 pairs from 4 patients. A perfect
   rank correlation appearing only after exclusions is a small-sample artifact, not a finding.

It is reported as exploratory and will not be presented as a survival result.

### 10.5 Prospective subset

Restricting to drugs started at or after the biopsy leaves **3 pairs from 2 patients**. No
inferential analysis is performed. These cases are described individually.

---

## 11. Conclusions this dataset supports

- Core-needle biopsies of metastatic PDAC can be maintained as organotypic slices with
  preserved architecture and immunophenotype for 7 days.
- MTS viability tracks histologic markers of proliferation and apoptosis and is a reasonable
  surrogate for the on-slide readout.
- A response threshold can be derived from assay noise by two independent routes that agree
  at 32% and 34% inhibition, and the resulting classification is insensitive to the exact
  cutoff across a 27.7–39.3% window. (Withdrawn in v1.3 and reinstated in v1.4 on well-level
  data; see §13.5 and §13.6.)
- Ex vivo viability shows a consistently negative but statistically inconclusive association
  with clinical response across every analysis set examined.

## 12. Conclusions this dataset does not support

- That the assay predicts clinical response. Three genuinely prospective pairs from two
  patients cannot establish this.
- Any diagnostic accuracy estimate. With two clinical responders the 2×2 table is degenerate.
- Any survival claim. PFS and OS are not computable from the available data.

---

## 13. Deviation log

### 13.1 — 19 August 2026: primary outcome redefined as binary clinical benefit (v1.0 → v1.1)

**What changed.** §5 primary outcome moved from the RECIST ordinal (responder = CR or PR) to
binary clinical benefit (CR/PR/SD = response, PD = no response). §6 endpoints restructured
accordingly, with the ordinal analysis retained as co-primary.

**When, and what was known.** After the v1.0 results in §10 were available. This is a post-hoc
change and is labeled as such wherever it appears.

**Who, and why.** Requested by the study team on clinical grounds: in metastatic PDAC treated
with cytotoxic combinations, radiographic SD is the common good outcome and CR is essentially
absent, so a CR/PR-only responder definition classifies most clinically successful courses as
failures. The v1.0 definition produced only 2 responders among 11 pairs, which made the 2×2
degenerate (every one of the 48 sensitivity scenarios had TP=1, FN=1; §10.2).

**Effect on the result.** Under clinical benefit, 8 unflagged pairs give TP 3, FP 0, FN 3,
TN 2 — sensitivity 50%, specificity 100% (Clopper–Pearson 16–100%), PPV 100% (29–100%),
NPV 40%, accuracy 62%, Fisher exact p = 0.46. The change makes the table less degenerate but
does not make it significant, and the confidence intervals remain uninformative: zero false
positives out of two possible negatives is not evidence of specificity.

**How this is handled in the manuscript.** Both definitions are reported. The Methods state
that the binary definition was adopted post hoc, the Results give both, and the Limitations
state that an outcome definition chosen after seeing the data inflates apparent performance
and that neither definition yields a statistically significant association.

### 13.2 — 19 August 2026: fitted regression line added to the TTNT panel

A least-squares line with a patient-level cluster-bootstrap band was added to the ex vivo
viability vs TTNT panel (Figure 3B) at the study team's request. The §10.4 caveats apply in
full and are restated in the figure legend: TTNT is not independent of RECIST response, five
of eight courses began before the biopsy, and the association reaches rho exactly −1.00 after
the QC gate. A regression line makes this relationship look stronger and more linear than the
underlying data support; it is presented as exploratory and hypothesis-generating only.

### 13.3 — 19 August 2026: chart abstraction incorporated, and three response adjudications

**What changed.** A chart abstraction was completed for the four case patients (PT-024, PT-015,
PT-014, PT-025), supplying stop dates, RECIST best responses, progression dates, 36 CA 19-9
values and 28 imaging assessments. These supersede the workbook's regimen rows for those four
patients. Other patients keep their Masterlist-derived courses. Matching rules are unchanged.

**Effect on the analysis set.** PT-024's gem/nab course now carries a recorded best response of
PR, so she enters the matched pairs for the first time and contributes two true positives. The
analysis set grows from 8 pairs / 5 patients to **10 pairs / 6 patients** (tissue-exhausted
excluded).

**Adjudication of therapy-sheet versus imaging-sheet discordance.** The abstraction's imaging
sheet disagreed with its own best-response column for four courses. Each was adjudicated by the
study team on clinical grounds:

| course | decision | stated basis |
|---|---|---|
| PT-025 cobimetinib | not scored; no pair added | day +8 PD scan reflects pre-existing progression, not failure of a 6-day course |
| PT-015 cobimetinib | PR rather than SD | CT at +84 and +182 both PR, on monotherapy before tovorafenib was added at +245 |
| PT-014 gem/nab | PD, as recorded | the +78 SD is superseded by PD at +139 and +166 |
| PT-024 capecitabine | no conflict | artifact of a missing stop date |

**This adjudication was performed after the ex vivo results were known, and after the effect of
each choice on the result had been calculated and displayed.** All three substantive calls
favoured the assay. That must be disclosed in the Methods, and the unadjudicated result reported
as a sensitivity analysis:

| analysis set | rho | accuracy |
|---|---|---|
| as adjudicated (primary) | −0.75 | 70% |
| as recorded, no adjudication | −0.65 | 70% |
| all three calls reversed | −0.47 | 55% |

A second reader blinded to the ex vivo results should confirm the three calls before submission.

**Final primary result.** 10 pairs from 6 patients. 2×2 at >30% inhibition: TP 5, FP 0, FN 3,
TN 2; sensitivity 62% (24–91), specificity 100% (16–100), PPV 100% (48–100), NPV 40% (5–85),
accuracy 70% (35–93), Fisher exact p = 0.444. Ordinal Spearman rho **−0.75**, patient-level
cluster-bootstrap 95% CI **−0.93 to 0.00**.

**The interval touches zero without excluding it,** and the upper bound is exactly 0.000 across
ten bootstrap seeds — a discreteness artifact of six clusters, where resamples drawing a single
patient repeatedly return rho = 0 exactly. The naive Spearman p of 0.0125 must **not** be
reported: it treats three patients' two pairs as independent. One-pair-per-patient enumeration
over all 16 subsets gives median rho −0.68, median p = 0.138, and 0 of 16 subsets reach p < 0.05.
**The primary endpoint remains negative.**

**Prospective subset.** 3 pairs from 3 patients, all three concordant (PT-024 gem/nab 1.5% → PR;
PT-015 MEK inhibitor 19.4% → PR; PT-014 gem/nab 93.2% → PD). Descriptive only.

**PFS.** Progression dates were abstracted for 2 courses in 1 patient. The §5 statement that true
PFS is not computable stands.

### 13.4 — 19 August 2026: PT-015 received FOLFIRINOX twice

Clarified by the study team that PT-015's two FOLFIRINOX courses are separate lines — course 1
(day −353, 9 cycles, PR, stopped for oxaliplatin intolerance) before the MEK inhibitor, and
course 2 (day +389, 3 cycles, PD) after an R0 resection.

The §5 tie-break ("where an agent was administered more than once, the course closest in time to
the biopsy is used") selects **course 1**, by a margin of 36 days (353 versus 389). Because ex vivo
FOLFIRINOX was 5.2% of control, this single tie-break decides whether that pair is the study's
strongest true positive or its only false positive.

**The pre-specified rule is retained.** It was written before these data existed, and revising it
now — in either direction — would be a post-hoc choice made with the consequences known. The
fragility is instead disclosed, and three alternatives reported as sensitivity analyses:

| pairing rule | 2×2 | specificity | accuracy | rho |
|---|---|---|---|---|
| A. SAP as written, course 1 (primary) | 5·0·3·2 | 100% (16–100) | 70% | −0.75 |
| B. post-biopsy course preferred, course 2 | 4·1·3·2 | 67% (9–99) | 60% | −0.41 |
| C. both courses as separate pairs | 5·1·3·2 | 67% (9–99) | 64% | −0.47 |
| D. pair excluded entirely | 4·0·3·2 | 100% (16–100) | 67% | −0.72 |

No configuration reaches statistical significance. The Limitations must state that the course
selected by the rule ended before the biopsy was taken, that the discarded course progressed, and
that by the time of course 2 the tumor had undergone resection, ~400 days of MEK inhibition, and
had acquired KRAS G12V on repeat profiling.

**For the prospective study,** a temporal-proximity criterion should be pre-specified in the
protocol — for example, courses within 365 days of the biopsy with no intervening curative-intent
resection — rather than selected retrospectively.

### 13.5 — 20 August 2026: half-normal derivation of the response threshold withdrawn (v1.2 → v1.3)

**What changed.** §4.2 previously gave two independent derivations of the response threshold.
The first — a half-normal fit to treated wells reading above 100% of control, quoted at
σ = 16.2 pp and a >32% cutoff — is withdrawn. The propagated-DMSO-CV derivation and the
**primary >30% cutoff are unchanged**.

**Why.** The σ = 16.2 pp figure could not be reproduced from any surviving data. Every
available source was checked on 20 August 2026:

| source | wells >100% of control | σ | implied cutoff |
|---|---|---|---|
| `mts_normalized_final.csv`, all | 5 (4 from PT-014) | 112 pp | 225% |
| `mts_normalized_final.csv`, QC-passing | 1 | 11.5 pp | 23% |
| Prism MTS files, de-duplicated, all | 13 | 102 pp | 205% |
| Prism MTS files, de-duplicated, QC-passing | **3** | 21.1 pp (bootstrap 9–33) | 42% |

`Masterlist.xlsx` was assumed to be the well-level source and is not: its single `Metadata`
sheet holds 34 rows of already-normalized `MTS_%`, with **no DMSO/control rows and no raw
absorbance column**, and no plate-reader export exists anywhere in the project folder. The only
well-level MTS data that exists is the Prism extraction, which yields three usable wells after
the QC gate — not enough to estimate σ, and not consistent with 16.2 pp.

**A second, independent objection.** Even with more data the method is weak here: only wells
whose true drug effect is near zero can read above 100%, so σ is estimated from a censored
tail and is biased and imprecise whenever the agents tested actually worked.

**Effect on the result.** None. The primary cutoff stays >30% inhibition and every reported
analysis is unchanged. What changes is the justification offered for it.

**Who, and when.** Decided by Kshitij on 20 August 2026, after the reproduction attempt above.

**How this is handled in the manuscript.** All "two independent derivations converged" language
is removed from the Abstract, Results and Methods. The threshold is described as derived from
control-well variability and, primarily, as cutoff-insensitive across 27.7–39.3% inhibition.
If the origin of σ = 16.2 pp is later established, this entry should be revisited rather than
the derivation quietly reinstated.

### 13.6 — 20 August 2026: the §13.5 withdrawal is REVERSED (v1.3 → v1.4)

**What changed.** §13.5 withdrew the half-normal derivation because σ = 16.2 pp could not be
reproduced. The well-level source has since been found and it reproduces. §4.2 is restored,
with recomputed numbers replacing the originals.

**Why the earlier attempt failed.** It searched `Masterlist.xlsx` and the Prism extract. The
Masterlist holds only aggregated `MTS_%` with no control rows; the Prism extract held a
fraction of the wells and mixed in IHC and nuclear-density sheets. The actual source is two lab
workbooks — `Liver Bx  Diaphgram Bx  MTS assay.xlsx` and
`PDAC Bx  Peritoneum Bx  MTS assay.xlsx` — supplied on 20 August 2026, which carry every
individual well, every DMSO control, the blank readings and the concentrations used.

**What the well-level data gives.** 311 well readings across the 13 cohort biopsies that had a
drug panel. Five cohort sheets are day-series feasibility runs with no treated-versus-DMSO
comparison and cannot contribute; Liver Bx 1 has no sheet.

| | v1.0 quoted | recomputed on well-level data |
|---|---|---|
| half-normal σ | 16.2 pp → 32% | **17.0 pp** (CI 11.1–22.0) → **34%** |
| DMSO CV median | 10.8% → 15.3 pp → 30% | **11.1%** → **15.8 pp** → **32%** |
| experiments with replicate control wells | 5 | **23** |
| experiments passing the QC gate | 3 | 6 |

The original σ = 16.2 pp sits inside the recomputed bootstrap interval, so the v1.0 figure was
sound; only its provenance had been lost.

**A caveat the original did not state.** The half-normal is fitted to a censored tail — only
wells whose true effect is near zero can exceed 100% — so σ is biased and imprecise whenever
the agents tested were active. The bootstrap interval (11.1–22.0 pp) is wide and should be
quoted with the point estimate.

**Effect on the result.** None. The primary cutoff stays >30% inhibition and every reported
analysis is unchanged.

**How this is handled in the manuscript.** The Abstract, Results and Methods state the two
derivations with their recomputed values (32% and 34%), note that >30% was adopted as slightly
more conservative, and retain the cutoff-insensitivity statement.

### 13.7 — 20 August 2026: provenance of the 20% QC gate stated, and its sensitivity reported (v1.4 → v1.5)

**What changed.** §4.3 item 1 previously asserted the 20% DMSO-CV gate without saying where
20% came from. It now states that the value was chosen rather than derived, gives the
distribution of all 23 experiments, and reports how the §4.2 threshold moves as the gate is
varied. **The gate itself is unchanged at 20%.**

**Why.** Asked directly where the 20% came from, the answer in `rigor_roadmap.md` is that it
was a round number picked to exclude one plainly broken experiment on an objective criterion.
That is a reasonable thing to have done and an unreasonable thing to leave unsaid: a reviewer
who asks the same question and finds no answer will discount everything downstream of it.

**What the check found.** The gate survives the full data better than expected. The 23 CVs
contain a 9.8-pp empty gap between 14.0% and 23.8%, and 20% sits inside it, so the cut is not
splitting a dense region. Between 15% and 25% the derived cutoff moves by less than one
percentage point. At 30% it more than doubles, to 72%.

**Effect on the result.** None. No analysis changes. What changes is that the gate is now
described honestly and its influence is quantified.

**How this is handled in the manuscript.** The Methods state that the gate was pre-specified at
20% before the well-level data existed, and justified afterwards by the gap in the observed
distribution; the sensitivity table is referenced to this SAP rather than reproduced in full.

### 13.8 — 20 August 2026: scope of the QC gate corrected (v1.5 → v1.6)

**What changed.** §4.3 was headed "Quality-control gates, applied before any analysis". It now
states where each gate is actually applied. **No gate value and no analysis result changed.**

**Why.** Tracing gate 1 onto the matched pairs showed the heading was false. Of the 10 analyzed
pairs, one passes, three fail, and six come from biopsies with a single 48-hour control well
and so have no assessable CV. Applying the gate strictly to the concordance would leave one
pair from one patient; leniently, seven pairs from four patients. The reported primary analysis
applies neither, and **includes PT-014, which fails the gate at 72–83%**.

Keeping the old heading would have meant the SAP claimed a procedure the analysis did not
follow — the kind of discrepancy a reviewer finds by reading two sections together.

**What is now stated.** Gates 1 and 2 govern the §4.2 noise-floor estimation and appear as a
dimension of the §7 sensitivity grid. Gates 3 and 4 (negative-absorbance flooring, and the
tissue-exhausted exclusion) apply to every analysis including the primary. The pair-level table
in §4.3 documents exactly who passes, fails, and is unassessable.

**Effect on the result.** None. What changes is that the gate's reach is now described
accurately, and PT-014's inclusion in the primary is disclosed rather than implicit.

**How this is handled in the manuscript.** The Methods sentence describing the gate no longer
says "excluded before analysis"; it states that the gate governs the noise-floor estimation and
is reported as a sensitivity analysis for the concordance. The Limitations gain a sentence
naming PT-014's inclusion.
