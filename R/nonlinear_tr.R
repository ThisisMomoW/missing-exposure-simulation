# Result schema for nonlinear binary TR-AIPW fits, including failed attempts.

empty_tr <- function() data.frame(
  estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
  tau1 = NA_real_, tau0 = NA_real_, n_boot_valid = NA_integer_,
  n_boot_invalid = NA_integer_, boot_median = NA_real_,
  boot_q025 = NA_real_, boot_q975 = NA_real_, boot_min = NA_real_, boot_max = NA_real_,
  q01_treated_ps = NA_real_, q99_treated_inv_weight = NA_real_,
  q99_control_inv_weight = NA_real_, q99_missingness_weight = NA_real_,
  max_missingness_weight = NA_real_, status = "failed", error = NA_character_)
