base_seed <- 123L
nsim <- 500L
mice_m <- 10L
mice_maxit <- 10L
trim_cutoff <- 0.025
trme_boot_B <- 100L
trme_shrink_rate <- 0.99
scenario_name <- "missing_share"
scenario_grid <- function() {
  g <- expand.grid(target_missing_share=c(0.2,0.4,0.5,0.6), n=c(500L,1000L,2000L), gamma_int=0, missing_mechanism="MAR", gamma_y=0.5, KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  row.names(g) <- NULL
  g
}
group_vars <- c("target_missing_share","n")
x_var <- "target_missing_share"
x_label <- "Target missing share"
