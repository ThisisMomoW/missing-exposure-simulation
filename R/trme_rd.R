# Adapted from Yuliang Shi's trme implementation for marginal risk differences.
#' @title Triple Robust Estimation with Missing Binary Exposure
#' @author Yuliang Shi
#' @description Estimate E[Y(1)] - E[Y(0)] using treatment, outcome,
#'   missingness and exposure-imputation models. Covariates and outcome must
#'   be observed. Under MAR, exposure missingness is independent of exposure
#'   conditional on the observed covariates and outcome.
#' @param covs Names of the observed confounders.
#' @param Y Name of the binary outcome column.
#' @param A Name of the binary exposure column, with NA for missing exposure.
#' @param data Data frame containing the analysis variables.
#' @param imp_model Fit exposure imputation directly if TRUE; use Bayes' rule otherwise.
#' @param shrink_rate Quantile for winsorizing inverse weights; the runners use 0.99.
#' @param ci_alpha Confidence level, 0.95 for 95% intervals.
#' @param method Estimating equations: TR-AIPW or TR-WEE. The study uses TR-AIPW.
#' @param bootstrap Whether to estimate uncertainty by resampling observations.
#' @param B Number of bootstrap resamples. The simulation runners pass 100.
#' @param seed Bootstrap seed; resample j uses seed + j.
#' @param nuisance_spec Main-effects or treatment-by-covariate interaction outcome model.
#' @param missing_mechanism MAR or MCAR observation model.
#' @return A list with the risk difference Estimate, intervention risks tau1/tau0,
#'   fitted propensity scores and missingness weights. With bootstrap inference,
#'   also returns valid resample estimates, their SD, percentile interval and
#'   valid/invalid resample counts. The summary includes a normal interval and
#'   Wald test of a zero risk difference; the study uses the percentile interval.
#' @details Requires at least max(20, ceiling(0.8 * B)) valid bootstrap estimates.
#'   Intervention risks must be in [0,1] and risk differences in [-1,1].
#'   Consistency requires two correct model groups among the treatment,
#'   outcome and missingness/imputation groups under the stated assumptions.
#' @references Yuliang Shi, Yeying Zhu, Joel Dubin. Causal Inference on Missing
#'   Exposure via Robust Estimator. https://github.com/yuliang-shi/trme
#'   Zhang, Z., Liu, W., Zhang, B., Tang, L., and Zhang, J. (2016). Causal inference
#'   with missing exposure information: Methods and applications to an obstetric
#'   study. Statistical Methods in Medical Research 25, 2053–2066.
#'   Williamson, E. J., Forbes, A., and Wolfe, R. (2012). Doubly robust estimators
#'   of causal exposure effects with missing data in the outcome, exposure or
#'   a confounder. Statistics in Medicine 31, 4382–4400.
#' @importFrom rootSolve multiroot
#' @importFrom parallel mclapply detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @export

trme_rd=function(covs,Y,A,data,imp_model=TRUE,shrink_rate=0.99,ci_alpha=0.95,
              method=c("TR-AIPW","TR-WEE"),bootstrap=TRUE,B=500,seed=1L,
              nuisance_spec=c("main","interaction"),
              missing_mechanism=c("MAR","MCAR"))
{

  ##match the method TR AIPW or TR WEE
  method=match.arg(method)
  nuisance_spec=match.arg(nuisance_spec)
  missing_mechanism=match.arg(missing_mechanism)

  cov_main <- paste(covs, collapse = "+")

  # DGP-informed comparison: treatment propensity and missingness models retain
  # their true main-effects forms. Only the outcome regression changes.
  ps_rhs <- cov_main
  imp_rhs <- paste0("Y+", cov_main)
  if (nuisance_spec == "interaction") {
    interaction_sum <- paste(covs, collapse = "+")
    or_rhs <- paste0("A+", cov_main, "+A:I(", interaction_sum, ")")
  } else {
    or_rhs <- paste0("A+", cov_main)
  }
  miss_rhs <- if (missing_mechanism == "MCAR") "1" else paste0("Y+", cov_main)
  f_miss <- as.formula(paste("r~", miss_rhs))
  f_ps   <- as.formula(paste("A~", ps_rhs))
  f_imp  <- as.formula(paste("A~", imp_rhs))
  f_or   <- as.formula(paste("Y~", or_rhs))
  f_ps_mm <- as.formula(paste("~", ps_rhs))
  f_or_mm <- as.formula(paste("~", or_rhs))


  ###start the main functions
  ##rename variables in data
  names(data)[names(data)==Y]="Y"
  names(data)[names(data)==A]="A"


  ###check input
  if((nrow(data)<1)|(ncol(data)<1))
    stop("Input dataset has either zero row or column.")

  if(length(levels(factor(data$Y)))>2)
    stop("not a binary outcome.")

  if(length(levels(factor(data$A)))>2)
    stop("not a binary treatment.")

  if((method!="TR-WEE")&(method!="TR-AIPW"))
    stop("You specify a wrong method.")


  ##test missing values
  data_covomit=na.omit(data[,covs])
  data_outomit=na.omit(data[,"Y"])

  if(nrow(data_covomit)<nrow(data))
    stop("covariates have missing values. Impute missing covariates first.")

  if(length(data_outomit)<nrow(data))
    stop("The outcome has missing values. Impute missing outcome first.")


  #####trme core function starts####
  trme_core=function(covs=covs,Y=Y,A=A,data=data,imp_model=imp_model,
                     shrink_rate=shrink_rate,ci_alpha=ci_alpha,method=method,
                     nuisance_spec=nuisance_spec,missing_mechanism=missing_mechanism)
  {
    ##convert data without factor
    df_type=lapply(data[,c("A","Y",covs)], FUN=class)
    fac_index=c(which(df_type=="character"),which(df_type=="factor"))

    if(length(fac_index)==1)
    {
      data[,c("A","Y",covs)][,fac_index]=as.numeric(as.factor(data[,c("A","Y",covs)][,fac_index]))-1

    }else{

      data[,c("A","Y",covs)][,fac_index]=lapply(data[,c("A","Y",covs)][,fac_index],
                                                FUN=function(x)
                                                {x=as.factor(x)
                                                y=as.numeric(x)-1
                                                return(y)})
    }


    ##create miss ind
    data$r=ifelse(is.na(data$A)==T,1,0)

    ##create baseline
    data$x0=1

    ##fit miss model
    glm_mar=glm(f_miss,data=data,family = binomial)
    data$fit_mar=pmin(pmax(glm_mar$fitted.values, 1e-6), 1-1e-6)

    ##inverse weights of missingessness
    data$weight_miss=1/(1-data$fit_mar)

    #shrinkage
    data$weight_miss[which(data$weight_miss>=quantile(data$weight_miss,shrink_rate,na.rm = TRUE))]=quantile(data$weight_miss,shrink_rate,na.rm = TRUE)


    ##fit naive ps model without adjusting for missingness
    glm_ps=glm(f_ps,data=data,family = binomial(link="logit"))
    glm_ps
    alpha_ps=coefficients(glm_ps)

    ##fit naive or model without adjusting for missingness
    glm_or=glm(f_or,data=data,family = binomial)
    glm_or_beta=coefficients(glm_or)


    ##fit imp model
    ##imp model A|X,Y by MAR, not affected by missing data
    if(imp_model==T){
      ##else: imp model is right.
      glm_ps_xy=glm(f_imp,data=data,family = binomial(link="logit")) #for A|X,Y

      ###Fit on complete data. Then, predict A|X,Y for all subjects
      data$pred_ps_xy=pmin(pmax(predict(glm_ps_xy,newdata=data,type="response"),
                                    1e-6), 1-1e-6)

    }else{

      ##note: for users, as long as they use imp_model=F, we will apply Bayes rule.
      ##if imp model wrong, but ps and or models are correct.
      ##Use Bayes to predict A|X,Y
      data$pred_ps_xy=rep(0,nrow(data))

      ##or model: subset on trt & observed
      ##trt model: subset on observed alpha_ps

      ##Bayes transfer A|X,Y. predict for all subjects
      link_glm_or_trt=exp(X_all%*%tr_estt)/(1+exp(X_all%*%tr_estt))
      link_glm_or_con=exp(X_all%*%beta_con)/(1+exp(X_all%*% beta_con))
      link_glm_ps=exp(X_all%*%alpha_ps)/(1+exp(X_all%*%alpha_ps))

      ##predict for all subjects in two different y=1 and y=0 model
      pred_ps_xy_y1=(link_glm_or_trt*link_glm_ps)/(link_glm_or_trt*link_glm_ps+link_glm_or_con*(1-link_glm_ps))
      pred_ps_xy_y0=((1-link_glm_or_trt)*link_glm_ps)/(1-(link_glm_or_trt*link_glm_ps+link_glm_or_con*(1-link_glm_ps)))

      ##match the fitted value with correct position with observed y=1 or y=0
      data$pred_ps_xy[data$Y==1]=pred_ps_xy_y1[data$Y==1]
      data$pred_ps_xy[data$Y==0]=pred_ps_xy_y0[data$Y==0]


      ##for extreme case, when sample is very small, and glm does not converge,
      ##so pred_ps_xy is na or nan. bayes rule can't help. we use the wrong imp model
      if(sum(is.nan(data$pred_ps_xy))>=1|sum(is.na(data$pred_ps_xy))>=1)
      {
        ##warning
        warning("Bayes method does not converged. Apply wrong imputation model.")

        ##fit wrong imp models without Y
        glm_ps_xy=glm(A~x1+x2+x3,data=data,family = binomial(link="logit")) #for A|X,Y

        ###Fit on complete data. Then, predict A|X,Y for all subjects
        data$pred_ps_xy=pmin(pmax(predict(glm_ps_xy,newdata=data,type="response"),
                                    1e-6), 1-1e-6)

      }
    }


    ##remove all NA
    data_naomit=na.omit(data)


    ##Design matrices matching the selected parametric nuisance specification.
    X_all=model.matrix(f_ps_mm, data=data)
    X_obs=model.matrix(f_ps_mm, data=data_naomit)

    XA_obs=model.matrix(f_or_mm, data=data_naomit)
    data_con=data; data_con$A=0
    data_trt=data; data_trt$A=1
    XA_con_all=model.matrix(f_or_mm, data=data_con)
    XA_trt_all=model.matrix(f_or_mm, data=data_trt)
    # Only counterfactual versions are needed for all subjects.
    XA_all=XA_con_all

    ##set diff of F all-Fobs
    F_diff=rep(0,ncol(X_all))


    #######Fit alpha EE ######

    ##score functions for alpha ee
    score_alpha_ee=function(alpha)
    {

      #dimnesion problem must use observed values to input the data
      link=X_obs%*%alpha
      link_all=X_all%*%alpha #for all subjects

      #five beta_ee_miss score functions
      score_alphaj=function(x_obs,x_all)
      {
        #only observed subjects
        u=(data_naomit$A*x_obs-(1-data_naomit$A)*x_obs*exp(link))/(1+exp(link))
        F_obs=sum((1-data_naomit$r)*data_naomit$weight_miss*u)


        ##all subjects
        u_fit=(data$pred_ps_xy*x_all-(1-data$pred_ps_xy)*x_all*exp(link_all))/(1+exp(link_all))
        F_all=sum((glm_mar$fitted.values-data$r)*data$weight_miss*u_fit)

        return(F_obs-F_all)

      }

      ##loop all columns for all covs and observed covs
      ##return F diff=F_obs-F_all to multiroot
      ##Then try the possible values s.t F diff is small enough to 0 and not change much
      for (j in 1:ncol(X_all)) {

        F_diff[j]=score_alphaj(x_obs=X_obs[,j],x_all=X_all[,j])

      }

      return(F_diff)
    }


    #solve alpha
    alpha_ee_sol=multiroot(f = score_alpha_ee, start =coefficients(glm_ps))
    alpha_ee_converge=alpha_ee_sol$estim.precis
    alpha_ee=alpha_ee_sol$root
    if (!is.finite(alpha_ee_converge) || alpha_ee_converge > 1e-6 ||
        any(!is.finite(alpha_ee))) {
      stop("Treatment propensity estimating equations did not converge.")
    }




    #######Fit beta EE ######

    ##score functions for beta miss
    score_beta_miss=function(beta)
    {
      # Weighted score for each outcome-design column, including A:I(sum X)
      # when the interaction specification is selected.
      mu_obs <- plogis(drop(XA_obs %*% beta))
      mu_con <- plogis(drop(XA_con_all %*% beta))
      mu_trt <- plogis(drop(XA_trt_all %*% beta))

      obs_multiplier <- (1-data_naomit$r) * data_naomit$weight_miss *
        (data_naomit$Y - mu_obs)
      F_obs <- colSums(XA_obs * obs_multiplier)

      expected_score <-
        XA_con_all * ((1-data$pred_ps_xy) * (data$Y-mu_con)) +
        XA_trt_all * (data$pred_ps_xy * (data$Y-mu_trt))

      all_multiplier <- (glm_mar$fitted.values-data$r) * data$weight_miss
      F_all <- colSums(expected_score * all_multiplier)

      score <- F_obs - F_all
      if (length(score) != length(beta) || any(!is.finite(score))) {
        stop("Outcome estimating equation returned an invalid score vector.")
      }
      score
    }

    #beta
    beta_ee_miss_sol=multiroot(f = score_beta_miss, start =coefficients(glm_or))
    beta_ee_miss_coverage=beta_ee_miss_sol$estim.precis
    beta_ee_miss=beta_ee_miss_sol$root
    if (!is.finite(beta_ee_miss_coverage) || beta_ee_miss_coverage > 1e-6 ||
        any(!is.finite(beta_ee_miss))) {
      stop("Outcome estimating equations did not converge.")
    }




    ###fitted values for PS
    ###A correctly specified propensity model supplies fitted exposure probabilities.
    ##if PS model is wrong, it is biased. So we require Miss+OR or OR+imp is correct.
    linear_dr_alpha=X_obs%*%alpha_ee
    fit_ps_dr=exp(linear_dr_alpha)/(1+exp(linear_dr_alpha))

    link_dr_alpha_all=X_all%*%alpha_ee
    fit_ps_dr_all=exp(link_dr_alpha_all)/(1+exp(link_dr_alpha_all))


    ####inverse fit ps for all subjects in control and treatment
    inv_ps_dr=1/fit_ps_dr
    inv_ps_dr_con=1/(1-fit_ps_dr)

    inv_ps_dr_all=1/fit_ps_dr_all
    inv_ps_dr_all_con=1/(1-fit_ps_dr_all)

    #shrink
    inv_ps_dr[which(inv_ps_dr>=quantile(inv_ps_dr,shrink_rate,na.rm = TRUE))]=quantile(inv_ps_dr,shrink_rate,na.rm = TRUE)
    inv_ps_dr_con[which(inv_ps_dr_con>=quantile(inv_ps_dr_con,shrink_rate,na.rm = TRUE))]=quantile(inv_ps_dr_con,shrink_rate,na.rm = TRUE)
    inv_ps_dr_all[which(inv_ps_dr_all>=quantile(inv_ps_dr_all,shrink_rate,na.rm = TRUE))]=quantile(inv_ps_dr_all,shrink_rate,na.rm = TRUE)
    inv_ps_dr_all_con[which(inv_ps_dr_all_con>=quantile(inv_ps_dr_all_con,shrink_rate,na.rm = TRUE))]=quantile(inv_ps_dr_all_con,shrink_rate,na.rm = TRUE)


    #####TR AIPW Method######


    ##reset A=NA as -100
    data$A[is.na(data$A)] = -100

    if(method=="TR-AIPW"){


      ###predicted response for trt and control groups
      data$m1=exp(XA_trt_all%*%beta_ee_miss)/(1+exp(XA_trt_all%*%beta_ee_miss))
      data$m0=exp(XA_con_all%*%beta_ee_miss)/(1+exp(XA_con_all%*%beta_ee_miss))

      if((is.na(alpha_ee_converge)==F)&(is.na(beta_ee_miss_coverage)==F))
      {
        ##for trt group
        dr_tau1=data$A*data$Y*inv_ps_dr_all-(data$A-fit_ps_dr_all)*inv_ps_dr_all*data$m1
        dr_aug_tau1=data$pred_ps_xy*data$Y*inv_ps_dr_all-(data$pred_ps_xy-fit_ps_dr_all)*inv_ps_dr_all*data$m1

        ##for control group
        dr_tau0=(1-data$A)*data$Y*inv_ps_dr_all_con-(fit_ps_dr_all-data$A)*inv_ps_dr_all_con*data$m0
        dr_aug_tau0=(1-data$pred_ps_xy)*data$Y*inv_ps_dr_all_con-(fit_ps_dr_all-data$pred_ps_xy)*inv_ps_dr_all_con*data$m0

        ###TR AIPW for trt group
        tr_est_tau1=mean((1-data$r)*data$weight_miss*dr_tau1)-
          mean((data$fit_mar-data$r)*data$weight_miss*dr_aug_tau1)

        ###TR AIPW for control group
        tr_est_tau0=mean((1-data$r)*data$weight_miss*dr_tau0)-
          mean((data$fit_mar-data$r)*data$weight_miss*dr_aug_tau0)


        ##TR AIPW
        tr_est=tr_est_tau1-tr_est_tau0
        tr_est

      }else{

        ##if alpha ee or beta ee is not converged
        tr_est=NA
        warning("Algorithm does not converged! Either enlarge sample size or specify the correct model.")
      }

    }



    #####TR WEE Method######
    if(method=="TR-WEE"){

      ##predicted response for trt or control

      exbit_trt_fit=exp(XA_trt_all%*%beta_ee_miss)/(1+exp(XA_trt_all%*%beta_ee_miss))
      exbit_con_fit=exp(XA_con_all%*%beta_ee_miss)/(1+exp(XA_con_all%*%beta_ee_miss))

      ##score functions for beta miss
      score_beta_wee1=function(beta)
      {
        #dimenesion problem must use observed values to input the data
        link_beta=X_all%*%beta
        exbit_beta=exp(link_beta)/(1+exp(link_beta))

        #five beta_ee_miss score functions
        beta_wee1=function(x_all)
        {
          ##trt group difference term
          F_trt=data$A*(data$Y-exbit_beta)-data$pred_ps_xy*(data$Y-exbit_trt_fit)

          ##dr term
          F_dr=data$pred_ps_xy*inv_ps_dr_all*(data$Y-exbit_trt_fit)

          ##final equation for all subjects
          F_all=sum(x_all*((1-data$r)*data$weight_miss*inv_ps_dr_all*F_trt+F_dr))

          return(F_all)
        }

        ##loop all columns for all covs and observed covs
        ##return F diff=F_obs-F_all to multiroot
        ##Then try the possible values s.t F diff is small enough to 0 and not changes much
        for (j in 1:ncol(X_all)) {
          F_diff[j]=beta_wee1(x_all=X_all[,j])
        }

        ##combine results with baseline, expousure, and covariates
        return(F_diff)
      }

      #beta
      beta_wee1_sol=multiroot(f = score_beta_wee1,
                              start =c(beta_ee_miss[1]+beta_ee_miss[2],beta_ee_miss[c(-1,-2)]))


      ###score function for beta wee0
      score_beta_wee0=function(beta){

        link_beta=X_all%*%beta
        exbit_beta=exp(link_beta)/(1+exp(link_beta))

        beta_wee0=function(x_all){

          ##control group difference term
          F_con=(1-data$A)*(data$Y-exbit_beta)-(1-data$pred_ps_xy)*(data$Y-exbit_con_fit)

          ##dr term
          F_dr=(1-data$pred_ps_xy)*inv_ps_dr_all_con*(data$Y-exbit_con_fit)

          ##final equation for all subjects
          F_all=sum(x_all*((1-data$r)*data$weight_miss*inv_ps_dr_all_con*F_con+F_dr))

          return(F_all)
        }

        for (j in 1:ncol(X_all)) {
          F_diff[j]=beta_wee0(x_all=X_all[,j])
        }

        ##combine results with baseline, expousure, and covariates
        return(F_diff)
      }

      #beta
      beta_wee0_sol=multiroot(f = score_beta_wee0,
                              start =c(beta_ee_miss[1],beta_ee_miss[c(-1,-2)]))

      ##estimated beta WEE
      beta_tr_wee1=beta_wee1_sol$root
      beta_tr_wee1_converge=beta_wee1_sol$estim.precis

      beta_tr_wee0=beta_wee0_sol$root
      beta_tr_wee0_converge=beta_wee0_sol$estim.precis

      ##pred response use WEE
      link_wee1=exp(X_all%*%beta_tr_wee1)/(1+exp(X_all%*%beta_tr_wee1))
      link_wee0=exp(X_all%*%beta_tr_wee0)/(1+exp(X_all%*%beta_tr_wee0))

      #####final TR WEE estimate####
      if((is.na(beta_tr_wee1_converge)==T)|(is.na(beta_tr_wee0_converge)==T))
      {
        ##if estimate tau is too large, remove TR WEE
        tr_est=NA
        warning("Algorithm does not converged! Either enlarge sample size or specify the correct model.")

      }else{

        ##estimate TR WEE only when alpha ee, beta ee, beta wee are converged
        ##and tr_tau_wee1,0 are not larger than 1
        tr_tau_wee1=mean((1-data$r)*data$weight_miss*(link_wee1-exbit_trt_fit)+exbit_trt_fit)

        tr_tau_wee0=mean((1-data$r)*data$weight_miss*(link_wee0-exbit_con_fit)+exbit_con_fit)

        tr_est=tr_tau_wee1-tr_tau_wee0
      }

    }

    data$A[which(data$A==-100)] =NA


    ##Point estimate on the marginal risk-difference scale.
    tau1 <- if(method=="TR-AIPW") tr_est_tau1 else tr_tau_wee1
    tau0 <- if(method=="TR-AIPW") tr_est_tau0 else tr_tau_wee0
    out=list("Estimate"=tr_est,"tau1"=tau1,"tau0"=tau0,
             "fit_ps"=as.vector(fit_ps_dr_all),"miss_weights"=data$weight_miss,data=data)
    return(out)

  }


  ##point estimate
  tr_est=trme_core(covs=covs,Y=Y,A=A,data=data,imp_model=imp_model,
                   shrink_rate=shrink_rate,ci_alpha=ci_alpha,method=method,
                   nuisance_spec=nuisance_spec,missing_mechanism=missing_mechanism)
  point_est=as.numeric(tr_est$Estimate)
  point_tau1=as.numeric(tr_est$tau1)
  point_tau0=as.numeric(tr_est$tau0)
  rd_from_components=point_tau1-point_tau0
  if (!is.finite(point_est) || !is.finite(point_tau1) || !is.finite(point_tau0) ||
      point_tau1 < 0 || point_tau1 > 1 ||
      point_tau0 < 0 || point_tau0 > 1 ||
      point_est < -1 || point_est > 1 ||
      abs(point_est-rd_from_components) > 1e-8) {
    stop("TR-AIPW point estimate failed binary-outcome validity checks.")
  }

  ###########Bootstrap Function ##############

  if(bootstrap==T){

    ##variance
    boot_fun=function(j) {

      set.seed(as.integer(seed) + j)

      ##size of data
      n=nrow(data)

      ##step1: bootstrap sample with replacement
      data_boot=data[sample.int(n,size=n,replace=TRUE),,drop=FALSE]

      ##A failed root solve in one bootstrap sample should not invalidate the
      ##entire Monte Carlo replicate. Return NA and assess the finite fraction.
      boot_fit=tryCatch(
        trme_core(covs=covs,Y=Y,A=A,data=data_boot,imp_model=imp_model,
                  shrink_rate=shrink_rate,ci_alpha=ci_alpha,method=method,
                  nuisance_spec=nuisance_spec,missing_mechanism=missing_mechanism),
        error=function(e) NULL
      )

      ## Numerical validity check for a binary outcome.  A finite root is not
      ## enough: intervention risks must lie in [0,1], and hence the RD in
      ## [-1,1]. Invalid roots are excluded from bootstrap inference.
      if (is.null(boot_fit)) return(NA_real_)

      tau1_boot <- suppressWarnings(as.numeric(boot_fit$tau1))
      tau0_boot <- suppressWarnings(as.numeric(boot_fit$tau0))
      rd_boot   <- suppressWarnings(as.numeric(boot_fit$Estimate))

      valid_boot <- length(tau1_boot)==1L && length(tau0_boot)==1L &&
                    length(rd_boot)==1L &&
                    is.finite(tau1_boot) && is.finite(tau0_boot) &&
                    is.finite(rd_boot) &&
                    tau1_boot >= 0 && tau1_boot <= 1 &&
                    tau0_boot >= 0 && tau0_boot <= 1 &&
                    rd_boot >= -1 && rd_boot <= 1 &&
                    abs(rd_boot - (tau1_boot - tau0_boot)) <= 1e-8

      if (!valid_boot) return(NA_real_)
      rd_boot
    }

    ##detect system
    sys=(Sys.info()[['sysname']])

    ##Respect the CPU allocation supplied by Slurm; never oversubscribe a node.
    slurm_cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
    if (!is.finite(slurm_cores) || slurm_cores < 1L) slurm_cores <- 1L
    n_cores <- min(slurm_cores, max(1L, detectCores()))

    if(sys=="Windows"){

      ##run foreach for windows.run parallel
      cl = makeCluster(n_cores)
      registerDoParallel(cl)

      boot_mat=foreach(i=1:B,.packages="rootSolve",.combine=cbind) %dopar% boot_fun(i)
      boot_vec=as.vector(boot_mat)

      stopCluster(cl)

    }else{

      ##run mclapply for mac and linux. run parallel
      boot_ls=mclapply(1:B, FUN=boot_fun,mc.cores=n_cores)

      ##cbind for each list, to get matrix
      boot_vec=as.vector(do.call(cbind,boot_ls))

    }

    ##Remove failed or theoretically invalid bootstrap roots.
    boot_vec=boot_vec[is.finite(boot_vec)]
    min_boot <- max(20L, ceiling(0.8 * B))
    if (length(boot_vec) < min_boot) {
      stop("Too few valid TR-AIPW bootstrap estimates: ",
           length(boot_vec), " of ", B, ".")
    }

    ##Bootstrap SD estimates uncertainty; percentile quantiles define the interval.
    boot_se=sd(boot_vec)

    alpha_tail <- (1-ci_alpha)/2
    ci_low_per=as.numeric(quantile(boot_vec,na.rm=TRUE,probs=alpha_tail,type=1))
    ci_up_per=as.numeric(quantile(boot_vec,na.rm=TRUE,probs=1-alpha_tail,type=1))
    ci_per=paste0("(",format(ci_low_per,drop0Trailing = F),",",format(ci_up_per,drop0Trailing = F),")")


    ###way2: CI BSE.
    ci_low_bse=round(point_est-qnorm(0.5+0.5*ci_alpha,0,1)*boot_se,3)
    ci_up_bse=round(point_est+qnorm(0.5+0.5*ci_alpha,0,1)*boot_se,3)
    ci_bse=paste0("(",format(ci_low_bse,drop0Trailing = F),",",format(ci_up_bse,drop0Trailing = F),")")

    ##test and pvalue using BSE
    ##H0: risk difference = 0 versus a two-sided alternative.
    # Display p-values below 0.01 as a threshold.
    z_obs=point_est/boot_se
    pvalue_bse=2*(1-pnorm(q=abs(z_obs),mean=0,sd=1))

    if(pvalue_bse<0.01){

      pvalue_bse="<0.01"
    }else{

      pvalue_bse=round(pvalue_bse,3)
    }

    ###summary df
    df_sum=data.frame(round(point_est,3),round(boot_se,3),ci_bse,ci_per,pvalue_bse)
    rownames(df_sum)=paste0(method,": ",A)
    colnames(df_sum)=c("Estimate","BSE",paste0(100*ci_alpha,"% CI Normal"),
                       paste0(100*ci_alpha,"% CI Percent"),"p.value")

    ##final output
    final=list("results"=df_sum,"method"=method,"estimand"="risk_difference",
               "Estimate"=point_est,"tau1"=tr_est$tau1,"tau0"=tr_est$tau0,
               "fit_ps"=tr_est$fit_ps,"miss_weights"=tr_est$miss_weights,
               "boot_est"=boot_vec,
               "boot_se"=boot_se,
               "ci_percentile"=c(lower=ci_low_per,upper=ci_up_per),
               "n_boot_valid"=length(boot_vec),
               "n_boot_invalid"=B-length(boot_vec),
               data=data)
  }

  if(bootstrap==F)
  {

    ##if not use bootstrap, only return point est
    df_sum=data.frame("Estimate"=point_est)
    df_sum=round(df_sum,3)
    rownames(df_sum)=paste0(method,": ",A)

    ##final output
    final=list("results"=df_sum,"method"=method,"estimand"="risk_difference",
               "Estimate"=point_est,"tau1"=tr_est$tau1,"tau0"=tr_est$tau0,
               "fit_ps"=tr_est$fit_ps,"miss_weights"=tr_est$miss_weights,
               data=data)
  }

  ##structure as trme class
  structure(final,class=c("trme_rd","trme"))

}


