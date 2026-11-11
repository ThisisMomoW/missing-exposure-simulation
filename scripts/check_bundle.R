arg <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
source(file.path(dirname(arg),"../R/project.R"))
source(file.path(root,"R/workflow.R"))
source(file.path(root,"R/estimators.R"))
x <- read_run("bundled");r<-x$records
stopifnot(!length(x$missing),nrow(r)==58000L,all(table(r$analysis)==500L),all(r$valid))
a <- x$manifest$analyses[match(r$analysis,x$manifest$analyses$analysis),]
mi <- a$method=="multiple_imputation";rf<-mi&a$specification=="rf"
stopifnot(all(r$mice_rf_ntree[rf]==500L),all(is.na(r$mice_rf_ntree[!rf])))
stopifnot(all(r$imputation_forests_completed[rf]==100L))
normal <- !mi
stopifnot(max(abs(r$ci_lower[normal]-(r$estimate[normal]-1.96*r$se[normal])))<1e-12,
 max(abs(r$ci_upper[normal]-(r$estimate[normal]+1.96*r$se[normal])))<1e-12)
components <- which(lengths(r$mice_q)>0)
for(i in components) {
  z <- rubin_pool(r$mice_q[[i]],r$mice_u[[i]])
  stopifnot(max(abs(z$value-unlist(r[i,c('estimate','se','ci_lower','ci_upper')])))<1e-12)
}
stats <- which(mi&is.finite(r$mice_ubar)&is.finite(r$mice_between)&is.finite(r$mice_valid))
m<-r$mice_valid[stats];u<-r$mice_ubar[stats];b<-r$mice_between[stats]
t<-u+(1+1/m)*b;df<-ifelse(b<=.Machine$double.eps,Inf,(m-1)*(1+u/((1+1/m)*b))^2)
stopifnot(max(abs(r$se[stats]-sqrt(t)))<1e-12,
 max(abs(r$ci_lower[stats]-(r$estimate[stats]-qt(.975,df)*sqrt(t))))<1e-12,
 max(abs(r$ci_upper[stats]-(r$estimate[stats]+qt(.975,df)*sqrt(t))))<1e-12)
s <- operating_characteristics(r)
reference <- utils::read.csv(file.path(root,'data/bundled/summary.csv'))
j<-match(s$analysis,reference$analysis)
fields<-c('bias','relative_bias','rmse','coverage','mc_sd','mean_se','failure_rate')
stopifnot(max(abs(as.matrix(s[fields])-as.matrix(reference[j,fields])),na.rm=TRUE)<1e-12)
for(i in 2:5) {
  d<-figure_summaries(s,i)
  stopifnot(all(d$attempted==500L),all(d$failed==0L))
  if(i==5L)stopifnot(all(d$interaction==1),setequal(d$n,c(500,1000,2000)),
    all(d$specification[d$method=='multiple_imputation'] %in% c('glm','rf')))
  path<-file.path(root,'outputs/figures/bundled',paste0('figure_',i,'_summary.csv'))
  if(file.exists(path)) {
    shown<-utils::read.csv(path);j<-match(paste(d$analysis,d$panel),paste(shown$analysis,shown$panel))
    stopifnot(!anyNA(j),max(abs(as.matrix(d[fields])-as.matrix(shown[j,fields])),na.rm=TRUE)<1e-12)
  }
}
message('PASS: 58000 unique attempts, complete figure selections, RF labels, normal and Rubin intervals, summary denominators and exported values.')
message(length(components),' records retain Q/U vectors; ',length(stats),' retain sufficient pooling statistics.')
