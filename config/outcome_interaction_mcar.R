base_seed <- 123L
nsim <- 500L
mice_m <- 10L
mice_maxit <- 10L
trim_cutoff <- 0.025
trme_boot_B <- 100L
trme_shrink_rate <- 0.99
scenario_name <- "outcome_interaction_mcar"
scenario_grid <- function() {
 g<-expand.grid(target_missing_share=0.4,n=c(500L,1000L,2000L),gamma_int=c(0,0.5,1),
                missing_mechanism="MCAR",gamma_y=0,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
 row.names(g)<-NULL; g
}
group_vars<-c("gamma_int","n")
x_var<-"gamma_int"
x_label<-"Outcome interaction coefficient (MCAR)"
