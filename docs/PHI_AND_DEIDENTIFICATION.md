# De-identification and release policy

## Statement

The data in `data/` contain no protected health information. No name, medical record
number, tissue bank accession, initials, calendar date, geographic subdivision, contact
detail or device identifier is present. Ages are reported in years and the cohort contains
no patient over 89.

## Procedure

De-identification runs on the study machine against the identified source workbooks. Those
workbooks and the re-identification crosswalk are not part of this repository and are not
distributable. `code/restricted/deidentify.py`, `deid_chart_abstraction.py`,
`extract_prism.py` and `extract_mts_wells.py` are the scripts that produce the
de-identified tables.

Three properties of the procedure matter for anyone reading the data:

**Dates become intervals, not shifted dates.** Every date is converted to an integer number
of days from that patient's own biopsy. Every interval needed for the analysis is preserved
and no calendar date is emitted, including dates embedded in free text, which become
`[D+123]` or `[D-45]` tokens.

**Identifiers are assigned, not derived.** Patients are `PT-001` through `PT-049`, assigned
from the normalized record number and stable across all source workbooks so that tables
join. The mapping exists only in the crosswalk on the study machine.

**Scrubbing is applied by field class.** Free-text fields are swept for names, record
numbers and accession numbers. Short categorical fields are scrubbed for dates only, and a
stop list of clinical abbreviations with a three-character minimum protects them: an early
version with a two-character minimum turned the literal value `No` into a redaction marker
wherever a patient's initials were `NO`, which silently destroyed outcome data. A final
global pass removes every crosswalk token from every cell, catching an identifier that
appears in another patient's row.

Sheet titles in the source workbooks sometimes carry a collection date. No such title is
written into this repository: `extract_mts_wells.py` resolves a worksheet by the stable
prefix of its title instead.

## Release policy for clinical variables

The de-identified tables held on the study machine carry more clinical metadata than this
repository releases. The released clinical variables are:

- demographics — age in years, sex, stage, grade, surgical approach, prior systemic therapy
- treatment intervals and progression — start, stop and progression day per regimen, and the
  time-to-next-therapy surrogate
- best overall response per regimen
- genomics — the genes and variants recorded in the study worksheet

Withheld, with the reason:

| withheld | reason |
|---|---|
| serial tumor marker results (CA19-9, CEA, CA125) | longitudinal clinical metadata beyond the released variables |
| scan-level RECIST target measurements | longitudinal clinical metadata beyond the released variables |
| line-by-line therapy log with free-text reasons for stopping | free text, and detail beyond the released variables |
| follow-up, vital status and cause of death | not required by the analysis released here |
| serial circulating tumor DNA results | not part of this report |
| per-patient notes and the clinical correlation worksheet | free text |
| extended sequencing report text | the released genomics are the driver findings recorded in the study worksheet; the fuller report text is not released |
| program-wide cohort and accrual tables | outside this report |

`code/restricted/build_release_dataset.py` applies this policy. It reads the restricted
de-identified directory and writes `data/`, printing every file it copies, every column it
removes and every file it withholds.

Two consequences for reproducibility are stated plainly in `docs/REPRODUCIBILITY.md`: the
clinical-course panels of Figure 4 and the serum marker rows of Table 1 cannot be rebuilt
from this repository.

## Screening

`code/checks/check_no_phi.py` screens every text file in the repository.

    python code/checks/check_no_phi.py --repo .

The pattern pass needs no restricted input. It flags calendar dates in any common format,
digit runs long enough to be a record number, government-identifier shapes and ages over 89.
Digit runs inside decimal fractions are measurements and are not flagged. Literals that are
part of the analysis are listed with their justification in `.phi-allowlist`; the file
currently holds one entry, the fixed random seed.

The token pass requires `--crosswalk` and is run on the study machine only. It loads every
identifier in the crosswalk and searches the repository for each one at token boundaries, so
a four-digit record number is not reported because those digits fall inside a measurement.
It prints per-file match counts and nothing else; no value from the crosswalk is printed,
stored or copied. Alphabetic fragments of four characters or fewer are reported for manual
review rather than as failures, because a short name fragment collides with ordinary English
and with identifiers in source code.

Vector figure files are screened as text. Their coordinate attributes are machine-generated
geometry and are removed before screening, so the record-number pattern is not applied to
them; every other pattern is.

Both passes were run against this repository before release. The pattern pass is clean. The
token pass reports zero failures and five items for review, each confirmed to be an ordinary
English or source-code word that coincides with a short name fragment: four are a variable
name or an everyday verb in source comments, and the fifth is that verb in a figure caption.

The figure files were screened separately. Their text layers and document metadata contain no
date, accession or record number, and optical character recognition of the rendered pages
found none burned into the images.
