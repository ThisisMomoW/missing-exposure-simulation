# Data generation, exposure missingness, intercept calibration and risk-difference
# truth for the Gaussian-covariate scenarios. R = 1 means exposure is observed.

# Small numerical bound used only to avoid exact probabilities of 0 or 1.
# This is not the simulation-level trimming discussed in the analysis.
prob_eps <- 1e-6

# Restrict estimated probabilities to (eps, 1 - eps) for numerical stability.
bound_probability <- function(p, eps = prob_eps) {
  pmin(pmax(p, eps), 1 - eps)
}

# Generate one simulated dataset.
#
# Outcome model:
# logit P(Y = 1 | A, X) = 0.7 + A + 0.5 X1 + 0.9 X2 + 0.7 X3
#                        + gamma_int * A * (X1 + X2 + X3).
#
# Missingness model:
# MCAR: logit P(R = 1) = R_intercept.
# MAR:  logit P(R = 1 | X, Y) = R_intercept
#       + 0.6 X1 + 0.7 X2 + 0.8 X3 + gamma_y Y.
#
# R = 1 means that A is observed.
simulate_dataset <- function(
    n,
    gamma_int = 0,
    missing_mechanism = "MAR",
    gamma_y = 0.5,
    R_intercept = -0.5
) {
  # Generate three continuous baseline covariates.
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rnorm(n)

  # Generate binary treatment/exposure A.
  pA <- plogis(-0.2 + 0.9 * X1 + 1.0 * X2 + 0.8 * X3)
  A <- rbinom(n, size = 1, prob = pA)

  # Generate binary outcome Y, optionally including treatment-covariate interaction.
  etaY <- 0.7 +
    1.0 * A +
    0.5 * X1 +
    0.9 * X2 +
    0.7 * X3 +
    gamma_int * A * (X1 + X2 + X3)

  pY <- plogis(etaY)
  Y <- rbinom(n, size = 1, prob = pY)

  # Construct the missingness linear predictor under MCAR or MAR.
  etaR <- switch(
    missing_mechanism,
    MCAR = rep(R_intercept, n),
    MAR = R_intercept +
      0.6 * X1 +
      0.7 * X2 +
      0.8 * X3 +
      gamma_y * Y,
    stop("missing_mechanism must be 'MCAR' or 'MAR'.")
  )

  # Generate the indicator that the exposure is observed.
  pR <- plogis(etaR)
  R <- rbinom(n, size = 1, prob = pR)

  # A_true is retained for full-data analyses; A with missing values is created later.
  data.frame(
    X1 = X1,
    X2 = X2,
    X3 = X3,
    A_true = A,
    Y = Y,
    R = R,
    pA = pA,
    pY = pY,
    pR = pR
  )
}

# Create the three analysis datasets from one generated dataset.
make_analysis_data <- function(dat) {
  # Full-data version: use the true exposure for everyone.
  full <- transform(dat, A = A_true)

  # Missing-exposure version: replace A by NA whenever R = 0.
  missing <- transform(
    dat,
    A = ifelse(R == 1, A_true, NA_real_)
  )

  # Complete-case version: retain only observations with observed A.
  complete_cases <- missing[missing$R == 1, , drop = FALSE]

  list(
    full = full,
    missing = missing,
    complete_cases = complete_cases
  )
}

# Calibrate the missingness intercept so that the expected missing proportion
# is approximately target_missing_share.
calibrate_R_intercept <- function(
    target_missing_share,
    gamma_int = 0,
    missing_mechanism = "MAR",
    gamma_y = 0.5,
    nval = 2e5,
    seed = 991L
) {
  # Under MCAR, the intercept has a closed-form solution.
  if (missing_mechanism == "MCAR") {
    return(qlogis(1 - target_missing_share))
  }

  # Under MAR, use a large calibration sample and numerical root-finding.
  set.seed(seed)

  X1 <- rnorm(nval)
  X2 <- rnorm(nval)
  X3 <- rnorm(nval)

  pA <- plogis(-0.2 + 0.9 * X1 + 1.0 * X2 + 0.8 * X3)
  A <- rbinom(nval, size = 1, prob = pA)

  pY <- plogis(
    0.7 +
      A +
      0.5 * X1 +
      0.9 * X2 +
      0.7 * X3 +
      gamma_int * A * (X1 + X2 + X3)
  )
  Y <- rbinom(nval, size = 1, prob = pY)

  target_observed_share <- 1 - target_missing_share

  calibration_equation <- function(intercept) {
    mean(
      plogis(
        intercept +
          0.6 * X1 +
          0.7 * X2 +
          0.8 * X3 +
          gamma_y * Y
      )
    ) - target_observed_share
  }

  uniroot(calibration_equation, interval = c(-20, 20))$root
}

# Calculate the true marginal risk difference by Monte Carlo integration.
# This must be recalculated when gamma_int changes because the logistic model
# is nonlinear and the treatment effect becomes heterogeneous.
true_rd <- function(gamma_int = 0, nval = 1e6, seed = 9001L) {
  set.seed(seed)

  X1 <- rnorm(nval)
  X2 <- rnorm(nval)
  X3 <- rnorm(nval)

  # Counterfactual outcome probability under A = 0.
  eta0 <- 0.7 +
    0.5 * X1 +
    0.9 * X2 +
    0.7 * X3

  # Counterfactual outcome probability under A = 1.
  eta1 <- 1.7 +
    0.5 * X1 +
    0.9 * X2 +
    0.7 * X3 +
    gamma_int * (X1 + X2 + X3)

  mean(plogis(eta1) - plogis(eta0))
}
