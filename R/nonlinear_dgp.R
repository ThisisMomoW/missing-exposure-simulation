# Binary nonlinear DGP, missingness calibration and risk-difference truth.

expit <- function(x) plogis(x)

make_covariates <- function(n) {
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  z3 <- rnorm(n)

  data.frame(
    Z1 = z1,
    Z2 = z2,
    Z3 = z3,
    X1 = exp(z1 / 2),
    X2 = z2 / (1 + exp(z1)) + 10,
    X3 = (z1 * z3 / 25 + 0.6)^3
  )
}

outcome_probability <- function(d, a, gamma_int) {
  expit(outcome_intercept + outcome_beta_a * a +
          0.5 * d$Z1 + d$Z2 + 0.5 * d$Z3 +
          gamma_int * d$Z1 * d$Z2)
}

calculate_truth <- function(gamma_int) {
  set.seed(truth_seed)
  sum_delta <- sum_delta2 <- sum_p0 <- sum_p1 <- 0
  remaining <- truth_n
  while (remaining > 0L) {
    b <- min(100000L, remaining)
    d <- make_covariates(b)
    p0 <- outcome_probability(d, 0, gamma_int)
    p1 <- outcome_probability(d, 1, gamma_int)
    delta <- p1 - p0
    sum_delta <- sum_delta + sum(delta)
    sum_delta2 <- sum_delta2 + sum(delta^2)
    sum_p0 <- sum_p0 + sum(p0)
    sum_p1 <- sum_p1 + sum(p1)
    remaining <- remaining - b
  }
  rd <- sum_delta / truth_n
  list(true_rd = rd, risk0 = sum_p0 / truth_n, risk1 = sum_p1 / truth_n,
       mc_se = sqrt(max(0, (sum_delta2 - truth_n * rd^2) /
                             (truth_n - 1)) / truth_n),
       n = truth_n, seed = truth_seed)
}

generate_complete_data <- function(n, gamma_int) {
  d <- make_covariates(n)

  # Treatment is generated parametrically from latent Z. The fitted treatment
  # propensity model sees only nonlinear X and is always a main-effects GLM.
  p_a <- expit(0.3 - 0.8 * d$Z1 + 0.8 * d$Z2 - 0.25 * d$Z3)
  a <- rbinom(n, 1, p_a)

  # Binary outcome with latent-Z predictors and a covariate interaction.
  # A constant log-odds coefficient does not imply a constant risk difference.
  p_y <- outcome_probability(d, a, gamma_int)
  y <- rbinom(n, 1, p_y)

  transform(d, A = a, Y = y, pA_true = p_a, pY_true = p_y)
}

missingness_linear_predictor <- function(d, gamma_y) {
  # MAR for the exposure: R depends on always-observed covariates and outcome.
  # R = 1 means that A is observed.
  0.5 * d$Z1 - 0.5 * d$Z2 + 0.25 * d$Z3 + gamma_y * d$Y
}

calibrate_r_intercept <- function(target_missing_share, gamma_int, gamma_y,
                                  seed, n_cal = 200000L) {
  set.seed(seed)
  d <- generate_complete_data(n_cal, gamma_int)
  eta <- missingness_linear_predictor(d, gamma_y)
  target_observed <- 1 - target_missing_share

  uniroot(
    function(alpha) mean(expit(alpha + eta)) - target_observed,
    interval = c(-30, 30),
    tol = 1e-10
  )$root
}

generate_one_dataset <- function(n, gamma_int, gamma_y, r_intercept) {
  d <- generate_complete_data(n, gamma_int)
  p_r <- expit(r_intercept + missingness_linear_predictor(d, gamma_y))
  r <- rbinom(n, 1, p_r)
  a_obs <- d$A
  a_obs[r == 0L] <- NA_integer_

  transform(d, R = r, A_obs = a_obs, pR_true = p_r)
}
