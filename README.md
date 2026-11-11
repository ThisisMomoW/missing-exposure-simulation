# Causal inference with missing exposures

This study compares estimation of a population average risk difference when a
binary exposure is partly observed. It includes MCAR, MAR and MNAR mechanisms,
treatment–covariate interactions, exclusion of the outcome from missing-exposure
models, and nonlinear covariates with GLM and random-forest analyses.

The methods are complete-data and complete-case AIPW, single and multiple
imputation followed by AIPW, and MR (`drcmd`). Estimation uses 10 outer folds.
MICE uses 10 imputations and 10 iterations; all RF models use 500 trees. The study
uses 500 fixed attempts per analysis. Equations, model specifications, inference
and failure handling are in [docs/STUDY.md](docs/STUDY.md).

MICE attempts all 10 completed-data fits and retains their results, warnings,
errors and counts. Rubin pooling uses the actual number of finite estimates
paired with finite, nonnegative variances. At least 2 valid pairs are required;
complete and partial successes are distinguished, and fewer than 2 is a failed
replication. Two is a computational minimum, not a guarantee of adequate
imputation precision. Override with `--min-valid-imputations N` (2–10).

Use R 4.4.3 and restore the project library:

```sh
Rscript scripts/setup.R
```

The lockfile records CRAN dependencies and the exact `drcmd` GitHub commit.
Compilation requires an R development toolchain; PDF figures require Cairo support.
Installed libraries, run checkpoints and scratch outputs are excluded from version
control. The bundled figure PDFs and numerical CSVs are included in
`outputs/figures/bundled/`.

The included replication records reproduce the numerical summaries in Figures
2–5; Figure 1 is a diagram. Generate all figures and their numerical CSVs without
fitting models:

```sh
Rscript scripts/aggregate.R --run bundled
Rscript scripts/figures.R --run bundled
```

A small execution check through the main runner is:

```sh
Rscript scripts/run.R --run check --cell nonlinear_1_500 --method complete_data --specification glm --reps 1:3
```

This checks execution, not statistical performance. Run all methods in one cell
or the entire grid:

```sh
Rscript scripts/run.R --run nonlinear --cell nonlinear_1_500
Rscript scripts/run.R --run study --workers 4
Rscript scripts/aggregate.R --run study
Rscript scripts/figures.R --run study
```

Run the same command again to resume missing records. Keep run names distinct
when changing selections or settings. `--reps` accepts ranges such as `1:3` or
`1:25,51:75`; omitting it retains 500 attempts. `--method`, `--specification` and
`--predictors` select analyses listed in `config/design.R`. A Slurm example is in
`slurm/study.sh`. Focused checks are described in [docs/VALIDATION.md](docs/VALIDATION.md).

- `R/`: data generation, estimators, pooling, summaries and plots.
- `config/`: study cells, analysis selections, numerical settings and seeds.
- `scripts/`: setup, execution, aggregation, figures and checks.
- `data/bundled/`: replication records, selection metadata and checksums.
- `outputs/`: run checkpoints, summary tables, figures and diagnostics.

See [NOTICE](NOTICE) for methodological and software attribution. A repository
license remains to be selected by the authors.
