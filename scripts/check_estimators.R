# Focused execution checks at the configured sample size and model settings.
arg0 <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
source(file.path(dirname(arg0), "../R/project.R"))
source(file.path(root,"R/workflow.R"))
if(length(commandArgs(TRUE))) stop("Usage: Rscript scripts/check_estimators.R")
out <- file.path(root,"outputs/validation");dir.create(out,recursive=TRUE,showWarnings=FALSE)
load_implementation("missing_share")
shared_pool <- pool_rubin
# Named columns preserve field identity even when a field is missing or columns move.
m <- matrix(c(.2,.1,0,.4),nrow=1,dimnames=list("Risk Difference",c("Estimate","SE","95% LCL","95% UCL")))
mock <- list(result=m[,c(4,2,1,3),drop=FALSE],summary=function(...) invisible(NULL))
stopifnot(identical(extract_aipw_result(mock,.025),data.frame(estimate=.2,se=.1,ci_lower=0,ci_upper=.4)))
for(column in colnames(m)) {
  bad <- mock;bad$result[,column]<-NA_real_
  stopifnot(inherits(try(extract_aipw_result(bad,.025),silent=TRUE),"try-error"))
}
load_implementation("nonlinear_binary")
stopifnot(identical(body(shared_pool),body(pool_rubin)))
for (q in list(c(.1,.2,.3),c(.2,.2,.2),c(NA,.2,.3),c(NA,NA,.3),c(NA,NA,NA))) {
  u <- c(.01,.02,.03)
  stopifnot(identical(shared_pool(q,u),pool_rubin(q,u)))
}
q<-c(.1,.2,.3);u<-c(.01,.02,.03);z<-pool_rubin(q,u)
stopifnot(z$estimate==mean(q),z$se==sqrt(mean(u)+(1+1/3)*var(q)),
  z$ci_lower<z$estimate-1.96*z$se,z$ci_upper>z$estimate+1.96*z$se)
.audit <- new.env();.audit$rf<-0L;.audit$logreg<-0L
trace(".ranger.donors",where=asNamespace("mice"),print=FALSE,
  tracer=quote(stopifnot(ntree==500L)),exit=quote({
    stopifnot(fit$num.trees==500L,ncol(nodes)==500L)
    .GlobalEnv$.audit$rf <- .GlobalEnv$.audit$rf+1L
  }))
trace("mice.impute.logreg",where=asNamespace("mice"),print=FALSE,
  tracer=quote(.GlobalEnv$.audit$logreg <- .GlobalEnv$.audit$logreg+1L))
cases <- data.frame(family=c(rep("nonlinear_binary",4),"outcome_interaction"),
  sid=c(rep(7L,4),1L),method=c("single_imp_aipw","single_imp_aipw",rep("mice_aipw",3)),
  spec=c("glm","rf","glm","rf","rf"))
checks<-list();results<-list()
for(i in seq_len(nrow(cases))) {
  c<-cases[i,];cfg<-load_implementation(c$family);sc<-cfg$scenario_grid()[c$sid,,drop=FALSE]
  constants<-cell_constants(c$family,c$sid,sc)
  dat<-generate_replication(c$family,c$sid,1L,sc,constants)
  .audit$rf<-0L;.audit$logreg<-0L
  message("Checking ",c$family," / ",c$method," / ",c$spec)
  rr<-fit_replication(c$family,c$sid,1L,c$method,c$spec,dat,sc)
  stopifnot(nrow(rr)==1L,rr$status=="ok",all(is.finite(unlist(rr[c("estimate","se","ci_lower","ci_upper")]))))
  want_rf<-if(c$method=="mice_aipw" && c$spec=="rf")100L else 0L
  want_glm<-if(c$method=="mice_aipw" && c$spec=="glm")100L else 0L
  stopifnot(.audit$rf==want_rf,.audit$logreg==want_glm)
  if(c$method=="mice_aipw") {
    stopifnot(length(rr$mice_q[[1]])==10L,length(rr$mice_u[[1]])==10L)
    rebuilt<-pool_rubin(rr$mice_q[[1]],rr$mice_u[[1]])
    stopifnot(identical(rr[names(rebuilt)],rebuilt))
  }
  ref<-readRDS(file.path(bundled_input_dir,paste0(c$family,".rds")))
  ref<-ref[ref$scenario_id==c$sid & ref$sim==1L & ref$method==c$method & ref$nuisance_spec==c$spec,]
  cols<-c("estimate","se","ci_lower","ci_upper")
  delta<-unlist(rr[cols])-unlist(ref[cols])
  stopifnot(nrow(ref)==1L,max(abs(delta))<1e-12)
  checks[[i]]<-data.frame(c,replication=1L,rf_forests=.audit$rf,logistic_imputations=.audit$logreg,
    max_absolute_difference=max(abs(delta)),status=rr$status)
  results[[i]]<-rr
  atomic_save(results,file.path(out,"estimator_results.rds"))
  write.csv(dplyr::bind_rows(checks),file.path(out,"estimator_checks.csv"),row.names=FALSE)
}
untrace(".ranger.donors",where=asNamespace("mice"))
untrace("mice.impute.logreg",where=asNamespace("mice"))
writeLines(capture.output(sessionInfo()),file.path(out,"estimator_session.txt"))
message("PASS: named extraction, pooling edge cases, imputation dispatch, runtime tree counts and five saved-data checks.")
