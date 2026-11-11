# Validation

Run the focused numerical checks and regenerate the bundled figures:

```sh
Rscript scripts/check_design.R
Rscript scripts/check_estimators.R
Rscript scripts/aggregate.R --run bundled
Rscript scripts/figures.R --run bundled
Rscript scripts/check_bundle.R
```

`check_design.R` checks all 19 design cells and 176 analysis specifications,
calibration, observed-exposure masking, figure selections and summary denominators.
Independent quadrature checks the population risk-difference integration for the
Gaussian experiments and the zero-interaction nonlinear experiment.

`check_estimators.R` checks Rubin pooling against direct calculations, including
zero between-imputation variance, complete and partial success, too few valid
components, all failures, warning/error retention and continuation after an error.
Tests cover the default minimum of 2 and stricter overrides, including exactly
two valid pairs, one valid pair, and unchanged pooling when all 10 fits succeed.

`check_bundle.R` checks unique replication keys, all 500 attempts per analysis,
normal intervals, Rubin intervals where components are recorded, RF-MICE tree
metadata and independently recalculated summary values. Missing component
information remains missing. Figure exports contain the numerical values and
actual point coordinates used by ggplot2.

For a longer execution check at the configured sample sizes and fitting settings:

```sh
Rscript scripts/check_execution.R
```

This fits one replication through the main-effects, interaction, outcome-exclusion
and nonlinear RF branches, and both MICE imputation models. It checks actual
500-tree RF-MICE fits and paired completed-data folds. It uses the default pooling
threshold and checks that all 10 fits succeed. Results go to `outputs/checks/`;
they are not Monte Carlo performance results. This command can take several minutes.

The project does not run the complete Monte Carlo study during setup or checking.
Numerical figure reproduction uses the included replication records. Refitting is
sensitive to software versions, including versions not recorded for every
included analysis. Vector PDFs match the manuscript's selections, method ordering,
colors, symbols, panel structure and reference lines; font rendering and panel
spacing may differ across graphics systems.

Checks completed with R 4.4.3 and the locked packages:

- The 58,000 bundled attempts have unique keys, 500 attempts per analysis and no
  invalid estimates or intervals. All 390 plotted summaries agree with the
  manuscript inputs to within 3.6e-13.
- Single-replication checks across 24 selected analyses reproduced their bundled
  estimates, SEs and confidence limits exactly. All completed-data fits succeeded.
  The RF-MICE check inspected 100 fitted forests, each with 500 trees.
- Two GLM MICE replications completed through the runner and aggregation with all
  component diagnostics retained. Four GLM analysis records were identical with
  one and two workers. Resuming a completed run reused its checkpoints.
- Pooling and failure checks passed, including a controlled initial-imputation
  error. Setup, all figure exports, and execution from another working directory
  were checked. No full Monte Carlo study was run as part of these checks.
