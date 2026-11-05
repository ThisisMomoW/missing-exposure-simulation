# Missing-exposure simulation

Research code for estimating the marginal causal risk difference, E[Y(1) − Y(0)], when a binary exposure is partly missing. A is exposure, Y is the binary outcome, X contains three observed covariates, and R = 1 denotes an observed exposure.

## Analyses and methods

The six simulation families vary missing exposure share (`missing_share`), outcome dependence in observation (`gamma_y`), treatment–covariate interaction under MAR or MCAR (`outcome_interaction`, `outcome_interaction_mcar`), observation mechanism (`missing_mechanism`), and nonlinear covariates (`nonlinear_binary`). Each method/design cell has **500 replications** at n = 500, 1,000 or 2,000.

The manuscript compares complete-data AIPW, complete-case AIPW, single-imputation AIPW, MICE AIPW, and **MR**, implemented by `drcmd::drcmd()` with code identifier `drcmd`. TR-AIPW (`tr_aipw`) is an additional analysis. The nonlinear manuscript cells use η = 1 (cells 7–9); η = 0 and 0.5 (cells 1–6) are available through `--cell` but are not included in the distributed dataset.

In Scenario 6, GLM single imputation uses logistic probabilities followed by Bernoulli draws, and GLM MICE uses `logreg`. The RF branches use RF probabilities for single imputation and MICE `rf` with **500 trees**. Outcome regression follows the GLM/RF branch; treatment propensity uses main-effects logistic regression throughout. MICE uses **10 imputations, 10 iterations and Rubin t intervals**. AIPW and `drcmd` use **10 outer folds**. Scenario 3 RF-MICE also uses 500 trees. [Analysis notes](report/ANALYSIS_NOTES.md) specify internal CV, seeds, bootstrap inference and summary denominators.

## Install and reproduce the figures

Use R 4.4.3 and run these commands from the repository root:

```sh
Rscript scripts/setup.R
Rscript scripts/aggregate.R bundled
Rscript scripts/manuscript_figures.R bundled
```

Setup restores the project library from `renv.lock`. See [system requirements and package sources](dependencies/ENVIRONMENT.md). The included `report/inputs/simulation/` records support all three manuscript comparisons without fitting models. Outputs in `outputs/bundled/manuscript_figures/` include PDF/PNG figures and the numerical summaries for each figure, with attempted and valid counts. The nonlinear figure uses η = 1, MAR, 40% missing exposure and all three sample sizes.

For the broader simulation report and its figures:

```sh
Rscript scripts/figures.R bundled
Rscript scripts/render.R html bundled
```

HTML requires Pandoc; replace `html` with `pdf` if LaTeX is installed.

## Run simulations

Run one scenario family, then summarize it:

```sh
Rscript scripts/run.R scenario nonlinear_binary --run=nonlinear
Rscript scripts/aggregate.R nonlinear
Rscript scripts/figures.R nonlinear
```

A two-replication execution check uses the same model settings; it does not assess statistical performance:

```sh
Rscript scripts/run.R scenario nonlinear_binary --cell=7 --rep=1,2 --method=mice_aipw --spec=glm --run=check
```

`Rscript scripts/run.R all --run=study` runs all six families, using the three manuscript cells for Scenario 6. Seeds depend on cell, method and replication ID. Repeat a command to resume completed checkpoints. Use a new run name when changing the selected grid, code or settings. Aggregation rejects missing or duplicate attempts and incompatible configurations; fitting failures remain recorded attempts.

## Files and validation

`R/` contains statistical functions, `config/` defines designs and settings, `scripts/` contains entry points, and `report/` contains the report and included data. `outputs/` holds generated results and is excluded from Git. The included dataset contains 309,000 attempts across 618 method/design cells; [analysis notes](report/ANALYSIS_NOTES.md) describe its scope and limitations.

Run `Rscript scripts/validate.R` for data, pooling and manuscript-summary checks. [Validation](validation/VALIDATION.md) describes focused estimator checks. See [software citations](dependencies/software_citations.txt) and [attribution](NOTICE). A project-wide license has not been selected.
