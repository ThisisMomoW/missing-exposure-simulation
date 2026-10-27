base_seed <- 123L
nsim <- 500L
mice_m <- 10L
mice_maxit <- 10L
trim_cutoff <- 0.025
trme_boot_B <- 100L
trme_shrink_rate <- 0.99
scenario_name <- "missing_mechanism"
scenario_grid <- function() {
  g <- expand.grid(target_missing_share=0.4, n=c(500L,1000L,2000L), gamma_int=0, missing_mechanism=c("MCAR","MAR"), gamma_y=0.5, KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  row.names(g) <- NULL
  g
}
group_vars <- c("missing_mechanism","n")
x_var <- "missing_mechanism"
x_label <- "Missingness mechanism"
