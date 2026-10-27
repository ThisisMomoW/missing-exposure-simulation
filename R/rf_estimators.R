# Random-forest AIPW estimators for the outcome-interaction scenario.

run_aipw_rf_result <- function(df, trim = DEFAULT_TRIM) {
  invisible(capture.output(
    fit <- AIPW::AIPW$new(
      Y = df$Y,
      A = df$A,
      W = df[, c("X1", "X2", "X3")],
      Q.SL.library = "SL.ranger",
      g.SL.library = "SL.glm"
    )$fit()
  ))

  summary_text <- capture.output(
    sm <- fit$summary(g.bound = trim)
  )

  # Read risk-difference estimates and standard errors from the package result.
  candidates <- list(sm, fit$result, fit$results)
  for (obj in candidates) {
    if (is.data.frame(obj) || is.matrix(obj)) {
      rn <- rownames(obj)
      hit <- grep("Risk Difference", rn, fixed = TRUE)
      if (length(hit)) {
        z <- as.numeric(obj[hit[1], ])
        z <- z[is.finite(z)]
        if (length(z) >= 2) {
          return(data.frame(
            estimate = z[1],
            se = z[2],
            ci_lower = if (length(z) >= 3) z[3] else z[1] - 1.96 * z[2],
            ci_upper = if (length(z) >= 4) z[4] else z[1] + 1.96 * z[2]
          ))
        }
      }
    }
  }

  line <- summary_text[
    grepl("Risk Difference", summary_text, fixed = TRUE)
  ][1]

  if (is.na(line)) {
    stop("Risk Difference row not found in AIPW summary.")
  }

  nums <- as.numeric(
    regmatches(
      line,
      gregexpr(
        "[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?",
        line,
        perl = TRUE
      )
    )[[1]]
  )

  out <- empty_result()
  out$estimate <- nums[1]

  if (length(nums) >= 2) {
    out$se <- nums[2]
  }

  if (length(nums) >= 4) {
    out$ci_lower <- nums[3]
    out$ci_upper <- nums[4]
  } else if (is.finite(out$se)) {
    out$ci_lower <- out$estimate - 1.96 * out$se
    out$ci_upper <- out$estimate + 1.96 * out$se
  }

  out
}
# Single exposure imputation using a complete-case random forest.

run_single_imp_rf_result <- function(df_missing, trim=DEFAULT_TRIM) {
  cc <- df_missing[!is.na(df_missing$A),,drop=FALSE]
  Xcc <- cc[,c("X1","X2","X3","Y"),drop=FALSE]
  Xall <- df_missing[,c("X1","X2","X3","Y"),drop=FALSE]

  imp_fit <- SuperLearner::SuperLearner(
    Y=cc$A, X=Xcc, family=binomial(),
    SL.library="SL.ranger"
  )

  ph <- as.numeric(predict(imp_fit, newdata=Xall, onlySL=TRUE)$pred)
  ph <- bound_probability(ph)

  d <- df_missing
  mis <- is.na(d$A)
  d$A[mis] <- rbinom(sum(mis), 1, ph[mis])

  run_aipw_rf_result(d, trim=trim)
}
# MICE random-forest exposure imputation followed by AIPW and Rubin pooling.

run_mice_full_rf_result <- function(df_missing,m=10L,maxit=10L,seed=1L,trim=DEFAULT_TRIM,ntree=500L) {
  d <- df_missing[,c("A","Y","X1","X2","X3")]
  d$A <- factor(d$A,levels=c(0,1))

  meth <- mice::make.method(d)
  meth[] <- ""
  meth["A"] <- "rf"

  pm <- mice::make.predictorMatrix(d)
  pm[,] <- 0
  pm["A",c("Y","X1","X2","X3")] <- 1

  imp <- mice::mice(
    d,m=m,maxit=maxit,
    method=meth,predictorMatrix=pm,
    seed=seed,printFlag=FALSE,
    rfPackage="ranger",ntree=ntree
  )

  rr <- lapply(seq_len(m), function(k){
    z <- mice::complete(imp,k)
    z$A <- as.integer(as.character(z$A))
    run_aipw_rf_result(z,trim=trim)
  })

  q <- vapply(rr,`[[`,numeric(1),"estimate")
  se <- vapply(rr,`[[`,numeric(1),"se")

  if(any(!is.finite(se))) return(data.frame(
    estimate=mean(q,na.rm=TRUE),
    se=NA_real_,ci_lower=NA_real_,ci_upper=NA_real_
  ))

  pool_rubin(q,se^2)
}
# Package drcmd with random-forest outcome and pseudo-outcome learners.

run_drcmd <- function(df_missing, trim=DEFAULT_TRIM) {
  fit <- drcmd::drcmd(
    Y=df_missing$Y,
    A=df_missing$A,
    X=df_missing[,c("X1","X2","X3")],
    m_learners="SL.ranger",
    g_learners="SL.glm",
    r_learners="SL.glm",
    po_learners="SL.ranger",
    cutoff=trim,
    k=10L
  )
  eo <- fit$results$estimates; so <- fit$results$ses
  est <- if(is.data.frame(eo)) eo$psi_hat_ate[1] else unlist(eo)["psi_hat_ate"]
  se <- if(is.data.frame(so)) so$psi_hat_ate[1] else unlist(so)["psi_hat_ate"]
  if(!is.finite(est) || !is.finite(se)) stop("non-finite drcmd result")
  data.frame(
    estimate=as.numeric(est), se=as.numeric(se),
    ci_lower=as.numeric(est)-1.96*as.numeric(se),
    ci_upper=as.numeric(est)+1.96*as.numeric(se)
  )
}
