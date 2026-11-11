# Study specification

The target is the population average risk difference, E[Y(1) − Y(0)]. A and Y
are binary; X contains three observed covariates; R = 1 means A is observed.
Covariates and outcomes are fully observed.

## Gaussian-covariate experiments

Independent X1, X2, X3 have standard normal distributions. The linear predictors
for A and Y are

```
logit P(A = 1 | X) = -0.2 + 0.9 X1 + X2 + 0.8 X3
logit P(Y = 1 | A, X) = 0.7 + A + 0.5 X1 + 0.9 X2 + 0.7 X3
                        + gamma_Y,AX A (X1 + X2 + X3)
logit P(R = 1 | X, Y, A) = gamma_R,0
                           + delta_R (0.6 X1 + 0.7 X2 + 0.8 X3 + gamma_R,Y Y)
                           + gamma_R,A A
```

Every cell has n = 2,000 and a target missing proportion of 0.4. MCAR sets
delta_R = gamma_R,A = 0. MAR sets delta_R = 1 and gamma_R,A = 0. MNAR sets
delta_R = 1, gamma_R,A = 2 and gamma_R,Y = 0.5.

| Figure | Design |
|---|---|
| 1 | Directed graphs for MCAR, MAR given X, MAR given X,Y and MNAR |
| 2 | MCAR; MAR given X; MAR given X,Y with gamma_R,Y = 0.5; MNAR. gamma_Y,AX = 0. |
| 3 | MCAR and MAR given X,Y with gamma_R,Y = 0.5; gamma_Y,AX = 0, 0.5, 1; main-effects and interaction specifications. |
| 4 | MAR with gamma_R,Y = 0 or 2; gamma_Y,AX = 0; exclude/include Y in imputation and both MR observation/pseudo-outcome regressions. |

The outcome interaction specification uses one coefficient for A(X1+X2+X3).
Exposure imputation adds YX1, YX2 and YX3 with separate coefficients. Exposure
propensity models use main-effects logistic regressions of A on X throughout.
Complete-data and complete-case AIPW do not change between the two predictor
columns in Figure 4.

## Nonlinear experiment

Independent Z1, Z2, Z3 are standard normal. Estimators observe

```
X1 = exp(Z1 / 2)
X2 = Z2 / (1 + exp(Z1)) + 10
X3 = (Z1 Z3 / 25 + 0.6)^3
logit P(A = 1 | Z) = 0.3 - 0.8 Z1 + 0.8 Z2 - 0.25 Z3
logit P(Y = 1 | A, Z) = 1 + 2 A + 0.5 Z1 + Z2 + 0.5 Z3 + gamma_Y,Z1Z2 Z1 Z2
logit P(R = 1 | Z, Y) = gamma_NL_R,0 + 0.5 Z1 - 0.5 Z2 + 0.25 Z3 + 0.5 Y
```

The grid crosses n = 500, 1,000, 2,000 with gamma_Y,Z1Z2 = 0, 0.5, 1.
Figure 5 uses gamma_Y,Z1Z2 = 1. GLM/RF specifies both the exposure imputation
and outcome models; exposure propensity and observation models remain logistic.
MR pseudo-outcome models are Gaussian identity-link GLMs or regression forests.

## Calibration and summaries

Observation intercepts target 40% missingness using 200,000 generated units per
scenario. The population ATE integrates conditional intervention risk differences
using 1,000,000 independent covariate draws for the Gaussian DGP and 2,000,000
for the nonlinear DGP. Each interaction value has its own target ATE.

For S valid estimates, relative bias is 100 mean(estimate − truth)/truth, RMSE
is sqrt(mean((estimate − truth)^2)), and coverage is the fraction of intervals
containing truth, inclusive of endpoints. Estimates and intervals are retained
at replication level, together with attempted and valid counts.

The design equations are specified in Sections 4.1 and Appendix B; calibration,
cross-fitting and performance definitions are in Appendix A. Figure captions
specify the displayed subsets. The design table retains the full nonlinear grid.

## Estimation

| Analysis | Fitting population | Outcome model | Exposure propensity |
|---|---|---|---|
| Complete-data AIPW | All units with A_true | Logistic GLM or probability forest | Main-effects logistic, A on X |
| Complete-case AIPW | R = 1 | Same specification | Same specification |
| Single imputation + AIPW | All units after one Bernoulli exposure draw | Same specification | Same specification |
| Multiple imputation + AIPW | All units in each completed dataset | Same specification | Same specification |
| MR (`drcmd`) | All units, with observation weights in designated regressions | Logistic GLM or probability forest | Main-effects logistic, A on X |

AIPW outcome regressions use A and X. These four AIPW analyses have no additional
observation weights; their influence functions retain the usual exposure
propensity weighting. Complete-case estimates and influence-function SEs average
over the complete cases.

MR uses main-effects logistic regressions for P(R = 1 | X,Y) and P(A = 1 | X).
Outcome and propensity fits use weights R/pi. The exposure-specific pseudo-outcome
regressions use X,Y, Gaussian identity-link losses and weights R, so only observed
exposures contribute to those regressions. Their RF specification uses regression
forests. Excluding Y removes it from the observation and pseudo-outcome predictors;
it remains the outcome response. MR interaction fits add A(X1+X2+X3) to the outcome
model only. The `drcmd` options are k = 10, cv_folds = 5, cutoff = 0.025,
eem_ind = FALSE, tml = FALSE, Rprobs = NA, parallel = FALSE, att = atc = FALSE.

AIPW and MR use 10 outer folds. AIPW's internal SuperLearner cross-validation uses
9 folds; MR requests 5, with the package's binary-response fold adjustment when
classes are sparse. Single-imputation RF uses 10 internal folds. Imputation takes
place before the completed-data AIPW cross-fitting.

Both packages calculate influence-function SEs as the sample SD of the risk-
difference influence values divided by the square root of the analysis sample
size. Their intervals use estimate +/- 1.96 SE. Estimated exposure propensity
and observation probabilities are bounded to [0.025, 0.975]. Single-imputation
probabilities are bounded to [1e-6, 1-1e-6] before Bernoulli draws.

## Imputation and forests

Single imputation fits A on X,Y among observed exposures and uses Bernoulli draws
for missing A. Multiple imputation uses MICE with 10 imputations and 10 iterations,
`method = "logreg"` for GLM and `method = "rf"` for RF. Only A is imputed. The
Y-exclusion experiment removes Y from the imputation predictors. Interaction
imputation includes the three separate YX terms.

All forests use ranger and 500 trees. SuperLearner probability forests use
min.node.size = 1; regression forests use 5. They use floor(sqrt(p)) candidate
variables, replacement sampling with sample.fraction = 1, and one thread. RF-MICE
explicitly passes ntree = 500 and rfPackage = "ranger". It fits classification
forests with min.node.size = 1 and floor(sqrt(p)) candidate variables and samples
observed donors from terminal nodes. This is distinct from a Bernoulli draw from
RF probabilities. The pinned mice/ranger implementations determine donor sampling
and the default RF-MICE thread count.

For usable completed-data estimates Q and variances U, let m denote the number
of usable pairs. Rubin pooling is

```
Qbar = mean(Q); Ubar = mean(U); B = var(Q)
T = Ubar + (1 + 1/m) B
nu = (m - 1) [1 + Ubar / ((1 + 1/m) B)]^2
CI = Qbar +/- qt(0.975, nu) sqrt(T)
```

When B <= machine epsilon, nu is infinite and the normal quantile is used.
Finite Q and finite, nonnegative U are required. Nonpositive pooled variance is
invalid. There is no finite-sample degrees-of-freedom adjustment.

MICE defaults to `min_valid_imputations = 2`. The optional
`--min-valid-imputations` argument accepts integers from 2 to 10. Two is a
computational minimum, not a guarantee of adequate imputation precision.
Every completed-data analysis is attempted independently. Errors do not prevent
subsequent analyses. Full-length Q/U vectors, component warnings and errors, and
requested/valid/failed counts are saved. Eligible finite pairs determine every m
in the pooling formula. Successful pools are classified as complete or partial;
partial successes remain eligible for summaries. Below-threshold pools and failed
initial imputation produce failed records without usable intervals. All 500
replication IDs are attempted once; failures are never replenished.

## Random numbers and execution

The RNG is Mersenne-Twister with Inversion normal generation and Rejection
sampling. For replication s in seed cell c, the data seed is 123 + 100000c + s.
Fitting adds 10000000 times the method index. Method indices follow the order in
`study_methods`; interaction fits add 6 and nonlinear RF fits add 5. Nonlinear
fitting also adds 90000000. MICE seeds are 123 + s in Gaussian experiments and the
fit seed in nonlinear experiments. Figure 2's MAR given X,Y cell uses seed cell 6;
Figure 3's zero-interaction MAR cell uses seed cell 3.

Truth seeds are 9001 (Gaussian) and 761239 (nonlinear). Calibration seeds are 991
and 123 + 1000000c + 777, respectively. MCAR uses logit(0.6) directly. Nonlinear
truth integration uses blocks of 100000 covariate draws. The include/exclude-Y
MICE comparison shares completed-data AIPW fold RNG states. Select both predictor
variants when running that comparison.

Parallel workers receive complete replication tasks and reset seeds before each
method. Worker count and job order do not determine random streams. A run name
fixes the analysis selection, replication IDs, software and configuration.
Repeating exactly the same command resumes missing records. Recorded failures are
retained. A changed configuration requires a different run name. A `.claim`
directory marks an active replication; remove it only after confirming its worker
has stopped. Aggregation rejects duplicate or unexpected keys and incompatible
configurations. It rejects absent attempts unless `--incomplete allow` is explicit.

## Included results

`data/bundled/replications.rds` contains 58000 unique attempts across 116 analyses,
500 per analysis, supporting Figures 2–5. Shared reference analyses are stored
once. All included estimates and intervals are valid. The remaining nonlinear
interaction settings are implemented in the grid and can be simulated.

`manifest.rds` identifies the analysis selection and recorded software versions;
`checksums.csv` checks file integrity. Saved component results and pooling
statistics are retained when available. NA denotes an unrecorded diagnostic.
3,000 records contain completed-data Q/U vectors and 6,000 contain sufficient
pooling statistics to recalculate intervals. Component counts are not inferred
from a successful pooled interval. The
`renv.lock` file defines the executable environment. Some distributed records lack
run-level package versions, so exact numerical agreement after refitting is not
established for every analysis.

Summary CSV files report attempted, valid and failed counts, bias, relative bias,
RMSE, Monte Carlo SD, mean estimated SE, coverage and Monte Carlo coverage SE.
Performance measures condition on valid results; `coverage_all` also reports
covered intervals divided by all attempts. MICE complete, partial, pooling-failed
and imputation-failed counts are separate. Figure CSVs contain every plotted
summary, with deterministic horizontal offsets recorded in coordinate CSVs.
