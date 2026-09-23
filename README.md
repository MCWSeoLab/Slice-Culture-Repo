# Organotypic slice culture of pancreatic cancer core biopsies — analysis code and released data

Code and de-identified data for the analysis reported in the manuscript on organotypic
slice culture of image-guided core biopsies in pancreatic ductal adenocarcinoma
(Seo Laboratory, Medical College of Wisconsin).

Every statistic, table and quantitative figure panel in the report is regenerated from
the tables in `data/` by `run_all.py`. `code/checks/verify_results.py` then compares the
regenerated values against the values printed in the manuscript; all 26 currently agree.

## Layout

    data/                de-identified released data (CSV)
    figures/             figure files as submitted
    code/lib/            shared plotting and panel helpers
    code/analysis/       statistics; runs from data/
    code/figures/        figure builds; run from data/
    code/restricted/     stages requiring inputs not released here (see below)
    code/checks/         PHI screen and results verification
    docs/                data dictionary, de-identification statement, reproducibility notes,
                         statistical analysis plan
    run_all.py           runs every stage that can run from data/
    results/             created by run_all.py; not tracked

## Requirements

Python 3.10 or later and the packages in `requirements.txt`:

    python -m pip install -r requirements.txt

Arial is used for all figure type. Where Arial is absent, matplotlib substitutes a
metrically identical face and the vector output is still tagged Arial, so the PDF and SVG
render with Arial on a machine that has it.

## Reproducing the analysis

    python run_all.py
    python code/checks/verify_results.py

`run_all.py` copies `data/` into `results/work/`, runs each stage there, and collects
console output into `results/logs/`, figures into `results/figures/` and derived tables
into `results/tables/`. The released tables are never modified. `--list` prints the stages
and `--only <substring>` runs a subset.

Every resampling routine uses a fixed seed, so repeated runs give identical numbers.

## What runs here, and what does not

`code/analysis/` and `code/figures/` run from `data/` alone. `code/restricted/` holds the
stages that cannot: they read either the identified source workbooks, the clinical tables
withheld under the release policy, or the histology image files, none of which are part of
this repository. Those scripts are included because they are the provenance of the released
tables and of the figure panels built from images. `docs/REPRODUCIBILITY.md` lists each one
and states what it needs.

## Data

`data/` contains no protected health information. Dates are carried only as integers
relative to each patient's own biopsy, so intervals are preserved and no calendar date is
present. Patients are identified as `PT-001` through `PT-049` and specimens by laboratory
label. The clinical variables released are demographics, treatment intervals and
progression, best overall response, and the genomic findings recorded in the study
worksheet; other clinical metadata is withheld. `docs/PHI_AND_DEIDENTIFICATION.md` states
the procedure and the release policy in full, and `code/checks/check_no_phi.py` is the
screen that enforces it.

## Citation and license

See `CITATION.cff`. Code is released under the MIT License (`LICENSE`). The data files are
released under CC BY 4.0.
