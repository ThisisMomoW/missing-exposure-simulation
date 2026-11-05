# Analysis settings and interpretation

## Design and estimators

A is binary exposure, Y is binary outcome, and X contains three observed covariates. R = 1 means A is observed. The target is the population marginal risk difference E[Y(1)] − E[Y(0)]. Truth is integrated from intervention outcome probabilities, not equated with a logistic coefficient. The [report](report.Rmd) specifies the generating equations and each scenario grid; executable definitions are in `../config/`.

Families 1–5 use main-effects and treatment–covariate interaction specifications. The interaction outcome model includes A × (X1 + X2 + X3); propensity models use main effects. Single and multiple logistic imputation can include Y-by-covariate interactions. Scenario 3 adds RF specifications. Scenario 6 generates latent Z and gives estimators nonlinear transforms X; its GLM/RF specifications change the outcome model, while propensity models remain GLMs. In Scenario 6 the GLM branch uses logistic exposure imputation, and the RF branch uses RF exposure imputation. Single imputation draws each missing exposure from a Bernoulli distribution; MICE uses `logreg` or `rf`, respectively. Both use Y, X1, X2 and X3 as predictors.

Full-data AIPW uses the simulated complete exposure. Complete-case AIPW uses observed exposures and averages over complete cases; under MAR this selected population can differ from the target population. Single-imputation AIPW draws missing exposures once and treats the completed values as observed when estimating uncertainty. MICE pools estimates and variances using Rubin's rules. `drcmd` calls `drcmd::drcmd()`; RF specifications use RF outcome and pseudo-outcome learners, with GLM treatment and missingness models. TR-AIPW uses Yuliang Shi's estimating equations adapted to the risk-difference scale. Its internal lowercase r indicates missing exposure, whereas the DGP's R indicates observed exposure.

## Settings

| Component | Families 1–5 | Scenario 6 |
|---|---|---|
| Monte Carlo design | 500 replications; n = 500, 1,000, 2,000 | Same |
| Truth integration | 1,000,000 draws; seed 9001 | 2,000,000 draws in blocks of 100,000; seed 761239 |
| Missingness calibration | 200,000 draws; seed 991; `uniroot` defaults, interval [−20,20]; MCAR closed form | 200,000 draws; seed `123 + cell*1000000 + 777`; tolerance 1e-10, interval [−30,30] |
| AIPW | 10 outer folds; 9 internal SuperLearner CV folds | Same |
| drcmd | 10 outer folds; 5 internal CV folds | Same |
| Probability bounds | AIPW and drcmd: 0.025/0.975; imputation probabilities: epsilon 1e-6 | Same |
| SuperLearner RF | 500 trees; mtry = floor(sqrt(p)); one thread; binomial probability forest; default node sizes | Same |
| RF single imputation | Scenario 3: fit on complete cases; 10 internal CV folds; no outer imputation cross-fit | RF specification only; GLM uses `glm(A ~ X1 + X2 + X3 + Y, family=binomial())` |
| MICE | 10 imputations, 10 iterations; logistic imputation for main/interaction models, RF for Scenario 3 RF | 10 imputations, 10 iterations; `logreg` for GLM and `rf` for RF |
| MICE RF | 500 trees; `rfPackage = "ranger"` | Same |
| MICE intervals | Rubin large-sample t degrees of freedom | Same, for both GLM and RF |
| TR-AIPW | 100 bootstrap resamples; at least 80 valid; 0.99 inverse-weight quantile shrinkage | Same; GLM main-effects specification only |

RF-MICE wrappers pass `ntree = 500` explicitly to `mice::mice()`, which forwards it to `mice.impute.rf()` and the ranger donor forest. The count comes from `execution_settings$mice_rf_trees`. SuperLearner forests use a separate 500-tree setting. The ranger backend and predictor matrix are specified in the estimator functions; other ranger tuning parameters use their package defaults.

AIPW uses package influence-curve standard errors and Wald intervals. `drcmd` uses the returned ATE standard error and a normal interval. TR-AIPW's SE is the sample SD of valid bootstrap estimates; its confidence interval uses empirical 2.5% and 97.5% quantiles with R's `type = 1`. Bootstrap intervention risks must lie in [0,1] and the risk difference in [−1,1]. Weight shrinkage applies to the point estimate and every bootstrap resample. These intervals are percentile intervals, so their coverage is not determined by a Wald SE ratio.

For MICE, let Q and U be the imputation estimates and their estimated variances. Pooling uses Qbar = mean(Q), Ubar = mean(U), B = var(Q), and T = Ubar + (1 + 1/m) B, where m counts usable imputations. With r = (1 + 1/m) B / Ubar, Rubin's degrees of freedom are df = (m − 1)(1 + 1/r)^2. The 95% interval is Qbar ± qt(0.975, df) sqrt(T). If B is at most machine precision, df is infinite and the critical value is qnorm(0.975). No additional degrees-of-freedom adjustment is applied.

Pooling excludes nonfinite Q/U pairs and negative U. With fewer than two usable imputations, SE and interval endpoints are NA; one usable estimate is returned as the point estimate, and zero usable estimates give NA. These rules are the same in all scenarios. Simulation MICE functions save the per-imputation vectors as `mice_q` and `mice_u` with each successful MICE result, allowing intervals to be recalculated without fitting the models again.

## Seeds and execution

The RNG is Mersenne–Twister with Inversion normal draws and Rejection sampling. Data seeds are `123 + cell*100000 + replication`. For families 1–5, GLM fit seeds add `method_id*10000000`, with method IDs 1–6 for main effects and 7–12 for interactions. MICE and TR receive internal seed `123 + replication`. Scenario 3 RF uses method IDs 1–5 and adds 90000000; its internal seed equals the fit seed. Scenario 6 uses IDs 1–10 for five methods × GLM/RF and the same 90000000 offset. Its TR fit seed is `500000000 + cell*1000000 + replication*1000`. Each TR bootstrap uses the internal seed plus its resample index.

The seed scheme reuses streams across families, and TR bootstrap streams can overlap across replications in families 1–5. Results across these groups should not be treated as independent simulation experiments. Matching software versions and RNG settings is necessary for reproducible numerical results.

Each named run records its selected cells, replications, methods, settings and package versions, including `mice_rf_trees = 500` and `s6_mice_ci_method = "rubin_t"`. RF-MICE result rows and aggregated tables also record that count; it is NA for methods without RF-MICE imputation. Results include seed identifiers, warnings and errors. Repeating the same command skips completed replication files; an unfinished fit restarts from its assigned seed. A different selection or code/configuration fingerprint requires a new run name. MICE rows also carry `mice_ci_method`; aggregation and cache reuse reject missing or incompatible interval methods. Execution is sequential by default. The TR bootstrap uses `SLURM_CPUS_PER_TASK` when available and assigns seeds by resample index.

Aggregation checks the expected replication keys and rejects duplicate, missing or unexpected records and incompatible configurations. A subset run is complete only for its selected grid. A fitting failure is a recorded attempt, not a missing file.

## Manuscript comparisons and included data

| Manuscript comparison | Code selection |
|---|---|
| Treatment–covariate interaction under MCAR/MAR | `outcome_interaction_mcar` and `outcome_interaction`; n = 2,000; missing share 0.4; gamma_int = 0, 0.5, 1; main/interaction GLMs |
| Outcome dependence under MAR | `gamma_y`; n = 2,000; missing share 0.4; gamma_int = 0; gamma_y = 0, 0.5, 0.8; main-effects GLM |
| Nonlinear covariates | `nonlinear_binary`; gamma_int = η = 1; n = 500, 1,000, 2,000; MAR; missing share 0.4; gamma_y = 0.5; GLM/RF |

MR is the manuscript name for method `drcmd`, which calls the package function `drcmd::drcmd()`. Its plot legend retains the identifier `drcmd`. TR-AIPW (`tr_aipw`) is a separate estimator and is not plotted in these three comparisons. In Scenario 6, η is a latent covariate–covariate interaction coefficient; it is not the treatment–covariate coefficient used in Scenario 3.

`report/inputs/simulation/` contains one RDS data frame per family, selected with run name `bundled`. `design.csv` specifies the included cells; `metadata.csv` gives shared inference settings; `checksums.csv` identifies the data files. Package versions recorded for the RF-MICE and Scenario 6 MICE fits appear in `package_versions.csv`. The environment lockfile specifies the tested runtime; package-version records are not available for every other replication.

| File | Attempts |
|---|---:|
| `missing_share.rds` | 72,000 |
| `gamma_y.rds` | 54,000 |
| `outcome_interaction.rds` | 76,500 |
| `missing_mechanism.rds` | 36,000 |
| `outcome_interaction_mcar.rds` | 54,000 |
| `nonlinear_binary.rds` | 16,500 |

The included nonlinear cells are 7–9 (η = 1), with 500 attempts for each of eleven method/specification combinations. The other six nonlinear cells are optional analyses supported by the runner. They have no distributed results. Scenario 3 RF specifications and TR-AIPW are also additional analyses.

Scenario 6 MICE records retain `mice_m_valid`, `mice_ubar`, `mice_b`, `mice_t` and `mice_df`. Along with the pooled estimate, these suffice to recompute their Rubin t intervals without fitting models. The other included MICE records contain pooled results but not component variances. New MICE fits in every family retain `mice_q` and `mice_u`; aggregation preserves these list columns in RDS outputs. Summary CSVs contain operating characteristics, not the imputation vectors.

The dataset retains 298 failed or invalid TR-AIPW fits. All other included methods have 500 jointly valid estimates and intervals per cell. Some failed fits have no recorded diagnostic message. `tables/failure_audit.csv` reports their counts and denominators. A successful pooled estimate need not mean all imputation-level fits were valid. The runner records warnings and errors where available.

## Summary definitions

Bias is mean(estimate − truth); RMSE is sqrt(mean((estimate − truth)^2)); Monte Carlo SD is the sample SD with denominator m−1; mean SE averages estimated SEs. SE ratio is mean SE divided by Monte Carlo SD. Coverage uses inclusive interval endpoints. Families 1–5 use the full-data main-effects Monte Carlo variance as the efficiency reference, including for RF comparisons.

For families 1–5, estimate-based summaries omit NA estimates, mean SE omits NA SEs, and coverage omits unavailable interval comparisons. Their denominators can differ; `n_finite` counts finite estimates rather than joint validity. Scenario 6 summaries require status `ok`, finite estimate/SE/endpoints, nonnegative SE and ordered endpoints. Its summaries are conditional on passing those checks.

`tables/failure_audit.csv` uses joint validity for every family. It reports attempted, valid and failed counts; failure rate = failed/attempted; coverage among valid fits; and `coverage_all`, which counts failures as unsuccessful coverage. Failed fits remain in attempted counts. A Monte Carlo SD requires at least two usable estimates. Very short runs check execution, not statistical performance.

The manuscript figures use finite point estimates with successful status for bias and RMSE, and finite ordered intervals for coverage. Point and interval counts are exported separately. This differs from the report's joint-validity rule for Scenario 6 when a point estimate is usable but its interval is unavailable. No finite outliers are removed.
