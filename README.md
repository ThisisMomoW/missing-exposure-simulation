# Missing-exposure simulation

Simulation study of estimators for the marginal causal risk difference, E[Y(1) − Y(0)], when a binary exposure is partly missing. A denotes exposure, Y the binary outcome, X three observed covariates, and R = 1 an observed exposure. The study examines bias, precision, interval coverage and fitting failures under missing completely at random (MCAR) and missing at random (MAR).

## Scenarios and estimators

| Scenario argument | Variation |
|---|---|
| `missing_share` | Missing exposure proportions 0.2, 0.4, 0.5, 0.6; MAR |
| `gamma_y` | Outcome coefficient in the observation model: 0, 0.5, 0.8 |
| `outcome_interaction` | Treatment–covariate interaction: 0, 0.5, 1; MAR |
| `missing_mechanism` | MCAR versus MAR |
| `outcome_interaction_mcar` | Treatment–covariate interaction: 0, 0.5, 1; MCAR |
| `nonlinear_binary` | Transformed covariates and latent covariate interactions |

The design uses sample sizes of 500, 1,000 and 2,000, with **500 Monte Carlo replications per cell**. Methods are full-data AIPW, complete-case AIPW, single-imputation AIPW, MICE AIPW, package `drcmd`, and TR-AIPW. Families 1–5 use main-effects and interaction models; `outcome_interaction` also includes random forests. `nonlinear_binary` uses GLM and RF specifications, with GLM only for TR-AIPW.

AIPW and `drcmd` use 10 outer cross-fitting folds. MICE uses **10 imputations and 10 iterations**, with **500 trees** for RF imputation; SuperLearner forests use **500 trees**. TR-AIPW estimates uncertainty from **100 bootstrap resamples**, requiring at least 80 valid estimates. See [analysis settings and limitations](report/ANALYSIS_NOTES.md) for nuisance models, inference and summary definitions.

## Installation

Use R 4.4.3 and run from the project root:

```sh
Rscript scripts/setup.R
```

This restores the project library from `renv.lock`, including AIPW, SuperLearner, drcmd, mice, ranger and rootSolve. Installation requires network access and may require compilers. HTML reports need Pandoc; PDF reports also need LaTeX. See [environment details](dependencies/ENVIRONMENT.md).

## Run simulations

A short execution check uses two replications with the full model settings. It does not assess statistical performance.

```sh
Rscript scripts/run.R scenario missing_share --cell=1 --rep=1,2 --run=check
```

Run one complete scenario family, then generate its tables and figures:

```sh
Rscript scripts/run.R scenario missing_share --run=missing_share
Rscript scripts/aggregate.R missing_share
Rscript scripts/figures.R missing_share
```

Results are written to `outputs/<run>/`, with summaries in `tables/` and plots in `figures/`. Seeds are determined by the cell, replication and method. Repeat the same command to resume completed checkpoints; choose a new `--run=NAME` when changing selections, code or settings. The RF-MICE tree count is recorded in run settings, result rows and summaries. Incomplete runs, duplicate results and incompatible configurations are rejected during aggregation.

To run every scenario, use `Rscript scripts/run.R all --run=study`, then aggregate and plot `study` as above. A complete six-family run can also be rendered with `Rscript scripts/render.R html study`.

## Use the included dataset

`report/inputs/mice_rf_10/` contains 342,000 replication records with **10-tree RF-MICE imputation**, including failed fits. Load this dataset explicitly to generate its tables, figures and report. A complete 500-tree RF-MICE result set is not included.

```sh
Rscript scripts/aggregate.R bundled_rf10
Rscript scripts/figures.R bundled_rf10
Rscript scripts/render.R html bundled_rf10
```

Outputs are under `outputs/bundled_rf10/`, and the report identifies its input tree count. Use `Rscript scripts/render.R pdf bundled_rf10` for PDF. `tables/failure_audit.csv` gives attempted, valid and failed counts and coverage denominators.

## Project files

`R/` contains statistical functions; `config/` defines scenarios and settings; `scripts/` provides command-line entry points; `report/` contains the report source, inputs and analysis notes. `dependencies/` documents the environment and references. Generated outputs and installed libraries are excluded from Git.

Run `Rscript scripts/validate.R` to check the included 10-tree dataset and summary calculations. See [validation scope](validation/VALIDATION.md), [software citations](dependencies/software_citations.txt) and [attribution and licensing](NOTICE). A project-wide license has not been selected.
