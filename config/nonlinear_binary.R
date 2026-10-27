base_seed <- 123L
nsim <- 500L

mice_m <- 10L
mice_maxit <- 10L
trim_cutoff <- 0.025

scenario_name <- "nonlinear_binary"
dgp_version <- "ks_binary_v1"
# Log-odds coefficient, NOT the marginal risk difference.
outcome_intercept <- 1
outcome_beta_a <- 2
truth_n <- 2000000L
truth_seed <- 761239L

chunk_size <- 25L
chunks_per_cell <- nsim %/% chunk_size

if (chunks_per_cell * chunk_size != nsim) {
  stop("nsim must be divisible by chunk_size.")
}

scenario_grid <- function() {
  g <- expand.grid(
    target_missing_share = 0.4,
    n = c(500L, 1000L, 2000L),
    gamma_int = c(0, 0.5, 1),
    missing_mechanism = "MAR",
    gamma_y = 0.5,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  row.names(g) <- NULL
  g
}

method_grid <- function() {
  expand.grid(
    method = c(
      "full_aipw",
      "cc_aipw",
      "single_imp_aipw",
      "mice_aipw",
      "drcmd"
    ),
    nuisance_spec = c("glm", "rf"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}
