arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(arg), "../R/project.R"))
source(file.path(root, "R/plots.R"))
out <- file.path(root, "outputs/checks")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
design <- study_design(); analyses <- study_analyses()
stopifnot(nrow(design) == 19L, nrow(analyses) == 176L, !anyDuplicated(analyses$analysis),
          all(design$replications == 500L), all(design$missing_share == .4))
stopifnot(identical(vapply(2:5, function(i) nrow(figure_design(i)), integer(1)), c(20L,60L,20L,30L)))
selection <- figure_design(4)
refs <- selection[selection$method %in% c("complete_data", "complete_case"), ]
stopifnot(all(table(refs$analysis) == 2L))

z <- matrix(c(1, 2, -1), nrow = 1)
stopifnot(abs(exposure_probability(z, "gaussian") - plogis(1.9)) < 1e-14,
          abs(exposure_probability(z, "nonlinear") - plogis(1.35)) < 1e-14,
          abs(outcome_probability(0, z, "gaussian", .5) - plogis(2.3)) < 1e-14,
          abs(outcome_probability(1, z, "gaussian", .5) - plogis(4.3)) < 1e-14,
          abs(outcome_probability(0, z, "nonlinear", .5) - plogis(4)) < 1e-14,
          abs(outcome_probability(1, z, "nonlinear", .5) - plogis(6)) < 1e-14)
expected_observation <- c(mcar_0=0, mar_xy2_0=3.2, mnar_0=3.7, nonlinear_1_500=-.25)
for(id in names(expected_observation)) {
  cell <- design[design$cell==id,]
  stopifnot(abs(observation_predictor(1,1,z,cell)-expected_observation[[id]])<1e-14)
}
set.seed(481)
cv <- covariates(20, "nonlinear")
recovered <- cbind(2*log(cv$measured$X1),
  (cv$measured$X2-10)*(1+cv$measured$X1^2),
  25*(sign(cv$measured$X3)*abs(cv$measured$X3)^(1/3)-.6)/(2*log(cv$measured$X1)))
stopifnot(max(abs(recovered-cv$latent)) < 1e-10)

# Quadrature supplies an independent target for normal linear predictors.
logistic_normal <- function(mu, sd) integrate(function(t) plogis(mu+sd*t)*dnorm(t),
                                            -Inf, Inf, rel.tol = 1e-10)$value
truth_checks <- lapply(c("mcar_0", "mcar_05", "mcar_1", "nonlinear_0_500"), function(id) {
  cell <- design[design$cell == id, ]
  set.seed(217)
  value <- population_ate(cell)
  exact <- if (cell$dgp == "gaussian") {
    logistic_normal(1.7, sqrt(sum((c(.5,.9,.7)+cell$interaction)^2))) -
      logistic_normal(.7, sqrt(1.55))
  } else logistic_normal(3, sqrt(1.5)) - logistic_normal(1, sqrt(1.5))
  stopifnot(abs(value$ate-exact) < 6*value$integration_se)
  data.frame(cell = id, estimate = value$ate, quadrature = exact,
             integration_se = value$integration_se, draws = value$draws)
})
write.csv(do.call(rbind, truth_checks), file.path(out, "truth.csv"), row.names = FALSE)

calibration_checks <- lapply(seq_len(nrow(design)), function(i) {
  cell <- design[i, ]
  set.seed(715+i)
  calibrated <- calibrate_observation(cell)
  stopifnot(abs(calibrated$expected_missing-.4) < 1e-5)
  set.seed(917+i)
  d <- generate_data(cell, calibrated$intercept)
  stopifnot(nrow(d) == cell$n, identical(is.na(d$A),d$R==0),
    all(d$A[d$R==1] == d$A_true[d$R==1]), !anyNA(d[c("X1","X2","X3","Y")]))
  if (cell$delta_r == 0) stopifnot(abs(calibrated$intercept-qlogis(.6)) < 1e-8)
  data.frame(cell = cell$cell, intercept = calibrated$intercept,
    expected_missing = calibrated$expected_missing, realized_missing = mean(d$R==0),
    calibration_draws = calibrated$draws)
})
write.csv(do.call(rbind, calibration_checks), file.path(out, "calibration.csv"), row.names = FALSE)

fixture <- data.frame(analysis = "example", replication = 1:3, estimate = c(.1,.3,NA),
  se = c(.1,.2,NA), ci_lower = c(0,.2,NA), ci_upper = c(.2,.4,NA), truth = .2,
  valid = c(TRUE,TRUE,FALSE))
s <- operating_characteristics(fixture)
stopifnot(s$attempted==3L, s$valid==2L, s$failed==1L, abs(s$bias)<1e-14,
  abs(s$rmse-.1)<1e-14, s$coverage==1, s$coverage_all==2/3)
stopifnot(inherits(try(operating_characteristics(rbind(fixture,fixture[1,])), silent=TRUE), "try-error"))
mechanism_diagram(file.path(out, "figure_1.pdf"))
writeLines(capture.output(sessionInfo()), file.path(out, "session.txt"))
message("PASS: study grid, figure selections, DGP equations, independent truth checks, calibration and summary denominators.")
