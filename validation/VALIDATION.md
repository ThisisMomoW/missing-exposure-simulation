# Validation scope

## RF-MICE execution checks

RF-MICE uses 500 trees with the ranger backend, 10 imputations and 10 iterations. Instrumented calls verified the argument at `mice::mice()`, `mice.impute.rf()`, `.ranger.donors()` and `ranger::ranger()`. The predictor set remained Y, X1, X2 and X3.

| Path | Check |
|---|---|
| Nonlinear binary, GLM outcome | One complete n=500 replication through `scripts/run.R`; all 100 imputation forests received 500 trees; fitting, pooling, aggregation and plotting succeeded without fitting warnings |
| Outcome interaction, RF outcome | Wrapper-to-backend probe using n=500 data; the first actual imputation forest received 500 trees |
| Nonlinear binary, RF outcome | Wrapper-to-backend probe using n=500 data; the first actual imputation forest received 500 trees |

The two RF-outcome probes stopped after the first imputation forest; their downstream RF AIPW fits were not executed. These checks establish parameter forwarding and limited execution, not Monte Carlo performance. No full simulation was run.

A manifest with a mismatched RF-MICE tree count is rejected before fitting. A cached result row with an incompatible count is rejected by both the runner and aggregation. Repeating an identical completed run resumes without refitting. A partial run cannot be rendered as a complete six-family report.

## Included dataset

`Rscript scripts/validate.R` checks the 342,000 records in `report/inputs/mice_rf_10/`, a dataset with 10-tree RF-MICE imputation. It verifies expected keys, design settings, truth consistency and failure accounting across six families and 684 method/design combinations. Deliberately duplicated, incomplete and altered grids must be rejected. A known-value fixture checks bias, RMSE, sample SD, mean SE and coverage/failure denominators. Results are written to `outputs/validation/`, with the input tree count identified.

This dataset contains 302 failed or invalid fits, including seven in the nonlinear binary family. They remain in attempted counts. See [analysis definitions and limitations](../report/ANALYSIS_NOTES.md) for validity criteria, conditional coverage and missing error details.

Aggregation and HTML rendering succeed when this dataset is explicitly selected as `bundled_rf10`. Reports print its recorded RF-MICE tree count. A complete 500-tree RF-MICE result set is not included. See [environment details](../dependencies/ENVIRONMENT.md) for installation and platform limitations.
