# Controlled component failures check pooling, diagnostics and continuation.
stopifnot(study_settings()$min_valid_imputations==2L)
q <- c(.12,.24,.18,.31,.2,.16,.22,.25,.19,.23)
u <- seq(.001,.01,length.out=10)
fixture <- lapply(1:10,function(j)list(index=j))
called <- integer()
errors <- integer(); invalid <- integer()
analyzer <- function(d,specification) {
  j <- d$index;called <<- c(called,j)
  if(j==3L)warning("component warning")
  if(j %in% errors)stop(paste("fit",j,"failed"))
  list(value=c(estimate=if(j %in% invalid)NA_real_ else q[j],se=sqrt(u[j]),ci_lower=0,ci_upper=1),diagnostics=list())
}
set.seed(517)
complete <- analyze_imputations(fixture,"glm",analyzer=analyzer)
m <- 10L;total <- mean(u)+(1+1/m)*var(q)
df <- (m-1)*(1+mean(u)/((1+1/m)*var(q)))^2
expected <- c(estimate=mean(q),se=sqrt(total),ci_lower=mean(q)-qt(.975,df)*sqrt(total),ci_upper=mean(q)+qt(.975,df)*sqrt(total))
stopifnot(identical(called,1:10),max(abs(complete$value-expected))<1e-14,
  complete$diagnostics$mice_pooling_status=="complete",complete$diagnostics$min_valid_imputations==2L,
  complete$diagnostics$mice_requested==10L,complete$diagnostics$mice_valid==10L,complete$diagnostics$mice_failed==0L)
strict <- analyze_imputations(fixture,"glm",minimum=10L,analyzer=analyzer)
stopifnot(identical(complete$value,strict$value),identical(complete$diagnostics$pooling,strict$diagnostics$pooling))

errors <- c(2L,7L);invalid <- 4L;called <- integer()
partial <- analyze_imputations(fixture,"glm",analyzer=analyzer)
v <- setdiff(1:10,c(errors,invalid));m <- length(v)
total <- mean(u[v])+(1+1/m)*var(q[v]);df <- (m-1)*(1+mean(u[v])/((1+1/m)*var(q[v])))^2
expected <- c(estimate=mean(q[v]),se=sqrt(total),ci_lower=mean(q[v])-qt(.975,df)*sqrt(total),ci_upper=mean(q[v])+qt(.975,df)*sqrt(total))
stopifnot(identical(called,1:10),max(abs(partial$value-expected))<1e-14,
  abs(partial$diagnostics$pooling$df-df)<1e-14,partial$diagnostics$pooling$m==m,
  partial$diagnostics$mice_requested==10L,partial$diagnostics$mice_valid==7L,partial$diagnostics$mice_failed==3L,
  length(partial$diagnostics$mice_q)==10L,length(partial$diagnostics$mice_u)==10L,
  all(is.na(partial$diagnostics$mice_q[c(2,4,7)])),all(is.na(partial$diagnostics$mice_u[c(2,7)])),
  partial$diagnostics$mice_pooling_status=="partial",all(partial$diagnostics$component_fits$attempted),
  identical(partial$diagnostics$component_fits$error[c(2,7)],c('fit 2 failed','fit 7 failed')),
  identical(partial$diagnostics$component_fits$warnings[[3]],'component warning'))
strict_partial <- analyze_imputations(fixture,"glm",minimum=10L,analyzer=analyzer)
stopifnot(all(is.na(strict_partial$value)),strict_partial$diagnostics$mice_valid==7L,
  strict_partial$diagnostics$mice_pooling_status=='failed')

# Exactly two valid pairs are eligible, including a zero within-imputation variance.
errors <- 1:8;invalid <- integer();called <- integer();u[9] <- 0
pair <- analyze_imputations(fixture,"rf",analyzer=analyzer)
v <- 9:10;m <- 2L
ubar <- mean(u[v]);between <- var(q[v]);total <- ubar+(1+1/m)*between
df <- (m-1)*(1+ubar/((1+1/m)*between))^2
expected <- c(estimate=mean(q[v]),se=sqrt(total),ci_lower=mean(q[v])-qt(.975,df)*sqrt(total),ci_upper=mean(q[v])+qt(.975,df)*sqrt(total))
stopifnot(identical(called,1:10),max(abs(pair$value-expected))<1e-14,
  abs(pair$diagnostics$pooling$df-df)<1e-14,pair$diagnostics$pooling$m==2L,
  pair$diagnostics$mice_valid==2L,pair$diagnostics$mice_failed==8L,pair$diagnostics$mice_pooling_status=='partial')
errors <- 1:9;called <- integer()
below <- analyze_imputations(fixture,"rf",analyzer=analyzer)
stopifnot(identical(called,1:10),all(is.na(below$value)),below$diagnostics$mice_valid==1L,
  below$diagnostics$mice_failed==9L,below$diagnostics$mice_pooling_status=='failed',
  below$diagnostics$mice_q[10]==q[10],is.finite(below$diagnostics$mice_u[10]),
  identical(below$diagnostics$component_fits$warnings[[3]],'component warning'))
errors <- 1:10;called <- integer()
failed <- analyze_imputations(fixture,"glm",analyzer=analyzer)
stopifnot(identical(called,1:10),all(is.na(failed$value)),failed$diagnostics$mice_valid==0L,
  failed$diagnostics$mice_failed==10L,all(!is.na(failed$diagnostics$component_fits$error)))

mixed <- rubin_pool(c(q[v],NA,.2,.2,.2,.2),c(u[v],.1,-1,NA,Inf,NaN))
stopifnot(max(abs(mixed$value-pair$value))<1e-14,mixed$components$m==2L)
zero <- rubin_pool(rep(.2,10),rep(.01,10))
stopifnot(is.infinite(zero$components$df),zero$value['se']==.1,
  abs(zero$value['ci_lower']-(.2-qnorm(.975)*.1))<1e-14)
for(minimum in list(1,11,NA,2.5,NULL,Inf,'2'))
  stopifnot(inherits(try(analyze_imputations(fixture,"glm",minimum,analyzer=analyzer),silent=TRUE),'try-error'))
stopifnot(inherits(try(rubin_pool(c(.1,.2),c(-1,NA)),silent=TRUE),'try-error'))

local({
  trace('mice',where=asNamespace('mice'),print=FALSE,tracer=quote({
    warning('Controlled imputation warning');stop('Controlled imputation failure')
  }))
  on.exit(untrace('mice',where=asNamespace('mice')))
  d <- data.frame(A=c(0,1,NA,1,0,NA),Y=c(1,1,0,0,1,0),X1=1:6,X2=6:1,X3=c(1,3,2,5,4,6))
  for(specification in c('main','interaction','glm','rf')) {
    z <- fit_multiple_imputation(d,specification,TRUE,123L)
    stopifnot(all(is.na(z$value)),z$diagnostics$mice_pooling_status=='imputation_failed',
      z$diagnostics$min_valid_imputations==2L,z$diagnostics$mice_requested==10L,
      z$diagnostics$mice_valid==0L,z$diagnostics$mice_failed==10L,
      length(z$diagnostics$mice_q)==10L,all(!z$diagnostics$component_fits$attempted),
      identical(z$diagnostics$imputation_error,'Controlled imputation failure'),
      identical(z$diagnostics$imputation_warnings,'Controlled imputation warning'))
  }
  invalid <- try(fit_multiple_imputation(d,'glm',TRUE,123L,minimum=1L),silent=TRUE)
  stopifnot(inherits(invalid,'try-error'),grepl('min_valid_imputations',as.character(invalid)))
})

# Successful partial pools remain eligible after checkpointing and aggregation.
pooled <- list(complete,partial,pair,below,failed)
rows <- lapply(seq_along(pooled),function(j) {
  z <- pooled[[j]];d <- z$diagnostics
  r <- data.frame(analysis='fixture',replication=j,estimate=z$value['estimate'],se=z$value['se'],
    ci_lower=z$value['ci_lower'],ci_upper=z$value['ci_upper'],truth=.2,valid=valid_estimate(z$value),
    mice_requested=d$mice_requested,mice_valid=d$mice_valid,mice_failed=d$mice_failed,
    mice_pooling_status=d$mice_pooling_status)
  r$mice_q <- I(list(d$mice_q));r$mice_u <- I(list(d$mice_u));r$mice_components <- I(list(d$component_fits));r
})
records <- do.call(rbind,rows)
local({
  path <- tempfile(fileext='.rds');on.exit(unlink(path))
  atomic_rds(records,path)
  restored <- readRDS(path);stopifnot(identical(records,restored))
  s <- operating_characteristics(restored)
  stopifnot(s$attempted==5L,s$valid==3L,s$failed==2L,s$mice_complete==1L,
    s$mice_partial==2L,s$mice_pooling_failed==2L)
})
message('PASS: default minimum 2; complete, partial, exactly-two, insufficient and failed pools; unchanged complete estimates/SEs/CIs; component diagnostics and aggregation.')
