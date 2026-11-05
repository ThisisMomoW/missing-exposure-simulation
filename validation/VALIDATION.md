# Validation

Run the data checks without fitting models:

```sh
Rscript scripts/validate.R
```

The checks verify checksums, design settings, unique replication keys and 500
attempts in each of 618 method/design cells. The included 309,000 records contain
298 failed or invalid fits: 190 in `missing_share`, 36 in `gamma_y`, 32 in
`outcome_interaction`, 15 in `missing_mechanism`, 22 in
`outcome_interaction_mcar`, and 3 in `nonlinear_binary`. Failures remain in the
data and attempted counts. All 298 failures occur in the additional TR-AIPW analyses. Known-value fixtures verify summary denominators.
Duplicate, incomplete, altered-design and incompatible-imputation records must
be rejected.

All 3,000 included Scenario 6 MICE intervals are recalculated from the saved
within- and between-imputation variances and Rubin degrees of freedom. The
resulting SEs and endpoints agree within 1e-14. Logistic rows have no RF tree
count; RF-MICE rows identify 500 trees, ten imputations and ten iterations.

`manuscript_reference.csv` records the 105 plotted method/design cells.
Validation recomputes bias, RMSE, coverage and valid counts directly from the
replication records. The figure script exports matching summaries for all three
comparisons. With the pinned Linux environment, all three PNG exports are
byte-identical to the manuscript figure files. PDF exports and HTML report
rendering also succeed. Image rendering on another platform may use different
fonts; the numerical checks do not depend on image rendering.

## Estimator checks

```sh
Rscript scripts/check_estimators.R
```

This fits one replication at n = 500 for Scenario 6 single imputation and MICE
under both GLM and RF specifications, and Scenario 3 RF-MICE. It uses the full
settings: ten imputations, ten iterations, ten outer AIPW folds and 500 RF-MICE
trees. The checks inspect both the requested tree count and the fitted ranger
forest and terminal-node matrix, and count 100 imputation forests per RF-MICE
fit. Logistic MICE must use 100 logistic imputation steps and no RF donor
forests. Estimates, SEs and interval endpoints are compared with the matching
included replication. The saved imputation vectors must reproduce the pooled
result exactly. A two-replication run through `scripts/run.R` also passes
aggregation, preserves its imputation vectors and resumes without refitting.

Additional fixtures check AIPW extraction with reordered columns and nonfinite
values, and Rubin pooling with zero between-imputation variance, invalid inputs
and fewer than two usable imputations. The shared and nonlinear pooling rules
must agree.

The five focused fits agree exactly with their included replication; see
[the recorded estimator checks](estimator_checks.csv).

Detailed check results and session information are written to
`outputs/validation/`. These checks establish execution and numerical agreement
for selected fits; they do not establish statistical performance across a new
Monte Carlo study. Included pooling components are available for Scenario 6;
other included MICE rows retain pooled results only. See
[analysis notes](../report/ANALYSIS_NOTES.md) and
[the environment record](../dependencies/ENVIRONMENT.md).
