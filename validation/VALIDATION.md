# Validation scope

## Scenario 6 Rubin t intervals

Scenario 6 uses the same pooling function body and arguments as Scenarios 1–5. Eight cases check positive, zero and near-zero between-imputation variance; zero within-imputation or total variance; invalid Q/U pairs; and fewer than two usable imputations. Outputs agree with the shared implementation.

Both MICE-GLM and MICE-RF dispatch paths were exercised with fixed per-imputation estimates and variances. These controlled checks replace imputation and AIPW fitting with test inputs: pooled estimates and SEs are identical under normal and t intervals, while endpoints follow Rubin's degrees-of-freedom formula. They do not constitute full RF-outcome fitting runs.

One live n=500 MICE-GLM replication completed with 10 imputations, 10 iterations and 500 RF trees, without fitting warnings. Its estimate (0.1851836276) and SE (0.0409829478) are unchanged. With Rubin df = 145.5434, the 95% interval is [0.1041850333, 0.2661822218], compared with the normal interval [0.1048570500, 0.2655102052]. Its per-imputation `mice_q` and `mice_u` are saved for further calculations. The shared estimator file and all numerical configuration settings are unchanged.

Normal-interval checkpoints cannot resume as t-based runs. Rows with mismatched interval-method labels are rejected during reuse or aggregation. Study-scale MICE t-based coverage still requires the [Scenario 6 runs listed in the analysis notes](../report/ANALYSIS_NOTES.md#scenario-6-mice-runs-needed-for-t-based-coverage); no full simulation was launched.

## RF-MICE execution checks

RF-MICE uses 500 trees with the ranger backend, 10 imputations and 10 iterations. Instrumented calls verified the argument at `mice::mice()`, `mice.impute.rf()`, `.ranger.donors()` and `ranger::ranger()`. The predictor set remained Y, X1, X2 and X3.

| Path | Check |
|---|---|
| Nonlinear binary, GLM outcome | One complete n=500 replication through `scripts/run.R`; all 100 imputation forests received 500 trees; fitting, pooling, aggregation and plotting succeeded without fitting warnings |
| Outcome interaction, RF outcome | Wrapper-to-backend probe using n=500 data; the first actual imputation forest received 500 trees |
| Nonlinear binary, RF outcome | Wrapper-to-backend probe using n=500 data; the first actual imputation forest received 500 trees |

A manifest with a mismatched RF-MICE tree count is rejected before fitting. A cached result row with an incompatible count is rejected by both the runner and aggregation. Repeating an identical completed run resumes without refitting. A partial run cannot be rendered as a complete six-family report.

## Included dataset

`Rscript scripts/validate.R` checks the 342,000 records in `report/inputs/mice_rf_10/`, a dataset with 10-tree RF-MICE imputation and normal-based Scenario 6 MICE intervals. It verifies expected keys, design settings, truth consistency and failure accounting across six families and 684 method/design combinations. Deliberately duplicated, incomplete and altered grids must be rejected. A known-value fixture checks bias, RMSE, sample SD, mean SE and coverage/failure denominators. Results are written to `outputs/validation/`, with the input tree count identified.

This dataset contains 302 failed or invalid fits, including seven in the nonlinear binary family. They remain in attempted counts. See [analysis definitions and limitations](../report/ANALYSIS_NOTES.md) for validity criteria, conditional coverage and missing error details.

Aggregation and HTML rendering succeed when this dataset is explicitly selected as `bundled_rf10`. Reports print its recorded RF-MICE tree count and Scenario 6 MICE interval method. Its pooled results do not contain the imputation-level variance information needed to recalculate t intervals, so its coverage remains normal-based. A complete 500-tree RF-MICE result set is not included. See [environment details](../dependencies/ENVIRONMENT.md) for installation and platform limitations.
