covariates <- function(n, dgp) {
  z <- matrix(stats::rnorm(3L * n), ncol = 3L)
  if (dgp == "gaussian") {
    x <- z
  } else if (dgp == "nonlinear") {
    x <- cbind(exp(z[, 1] / 2), z[, 2] / (1 + exp(z[, 1])) + 10,
               (z[, 1] * z[, 3] / 25 + 0.6)^3)
  } else stop("Unknown data-generating process: ", dgp)
  colnames(z) <- paste0("Z", 1:3)
  colnames(x) <- paste0("X", 1:3)
  list(latent = z, measured = as.data.frame(x))
}

exposure_probability <- function(z, dgp) {
  if (dgp == "gaussian") stats::plogis(-0.2 + .9*z[,1] + z[,2] + .8*z[,3])
  else if (dgp == "nonlinear") stats::plogis(.3 - .8*z[,1] + .8*z[,2] - .25*z[,3])
  else stop("Unknown data-generating process: ", dgp)
}

outcome_probability <- function(a, z, dgp, interaction) {
  if (dgp == "gaussian") {
    lp <- .7 + a + .5*z[,1] + .9*z[,2] + .7*z[,3] + interaction*a*(z[,1]+z[,2]+z[,3])
  } else if (dgp == "nonlinear") {
    lp <- 1 + 2*a + .5*z[,1] + z[,2] + .5*z[,3] + interaction*z[,1]*z[,2]
  } else stop("Unknown data-generating process: ", dgp)
  as.numeric(stats::plogis(lp))
}

observation_predictor <- function(a, y, z, design) {
  if (design$dgp == "gaussian") {
    as.numeric(design$delta_r * (.6*z[,1] + .7*z[,2] + .8*z[,3] + design$gamma_ry * y) +
      design$gamma_ra * a)
  } else {
    as.numeric(.5*z[,1] - .5*z[,2] + .25*z[,3] + .5*y)
  }
}

generate_complete <- function(n, design) {
  cv <- covariates(n, design$dgp)
  a <- stats::rbinom(n, 1L, exposure_probability(cv$latent, design$dgp))
  y <- stats::rbinom(n, 1L, outcome_probability(a, cv$latent, design$dgp, design$interaction))
  list(data = data.frame(cv$measured, A = a, Y = y), latent = cv$latent)
}

# R = 1 denotes an observed exposure. Only measured X reaches the estimators.
generate_data <- function(design, observation_intercept) {
  full <- generate_complete(design$n, design)
  p <- stats::plogis(observation_intercept +
    observation_predictor(full$data$A, full$data$Y, full$latent, design))
  r <- stats::rbinom(design$n, 1L, p)
  full$data$A_true <- full$data$A
  full$data$R <- r
  full$data$A[r == 0L] <- NA_integer_
  full$data
}

# Calibration averages over generated exposures and outcomes, including A under MNAR.
calibrate_observation <- function(design, draws = 200000L) {
  if (design$dgp == "gaussian" && design$delta_r == 0 && design$gamma_ra == 0)
    return(list(intercept = stats::qlogis(1-design$missing_share),
      expected_missing = design$missing_share, draws = 0L, residual = 0))
  full <- generate_complete(draws, design)
  lp <- observation_predictor(full$data$A, full$data$Y, full$latent, design)
  target <- 1 - design$missing_share
  f <- function(intercept) mean(stats::plogis(intercept + lp)) - target
  intercept <- if (design$dgp == "gaussian") stats::uniroot(f, c(-20, 20))$root
    else stats::uniroot(f, c(-30, 30), tol = 1e-10)$root
  list(intercept = intercept, expected_missing = 1 - mean(stats::plogis(intercept + lp)),
       draws = draws, residual = f(intercept))
}

# The estimand integrates intervention risks, rather than a logistic coefficient.
population_ate <- function(design, draws = if (design$dgp == "gaussian") 1000000L else 2000000L) {
  if (design$dgp == "gaussian") {
    z <- covariates(draws, design$dgp)$latent
    rd <- outcome_probability(1, z, design$dgp, design$interaction) -
          outcome_probability(0, z, design$dgp, design$interaction)
    return(list(ate = mean(rd), integration_se = stats::sd(rd) / sqrt(draws), draws = draws))
  }
  total <- squares <- 0
  remaining <- draws
  while (remaining > 0L) {
    n <- min(100000L, remaining)
    z <- covariates(n, design$dgp)$latent
    rd <- outcome_probability(1, z, design$dgp, design$interaction) -
          outcome_probability(0, z, design$dgp, design$interaction)
    total <- total + sum(rd); squares <- squares + sum(rd^2)
    remaining <- remaining - n
  }
  ate <- total/draws
  list(ate = ate, integration_se = sqrt(max(0,(squares-draws*ate^2)/(draws-1))/draws), draws = draws)
}
