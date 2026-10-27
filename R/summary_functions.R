suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2)})
base_labels<-c(full_aipw="Full AIPW",cc_aipw="Complete-case AIPW",single_imp_aipw="Single-imputation AIPW",
 mice_aipw="MICE + AIPW",drcmd="drcmd",tr_aipw="TR-AIPW (Shi trme, RD)")

nuisance_labels<-c(
 main="Main effects GLM",
 interaction="Interaction GLM",
 rf="Random Forest"
)

add_estimator_label<-function(x){
 nuisance_label<-unname(nuisance_labels[x$nuisance_spec])
 nuisance_label[is.na(nuisance_label)]<-x$nuisance_spec[is.na(nuisance_label)]
 x$Method<-paste0(base_labels[x$method]," (",nuisance_label,")")
 x
}
# Missing estimates, SEs and intervals are omitted separately in these summaries.
# Attempted counts include failed fits; failure_audit reports joint validity.
summarise_oc<-function(x,true_col="true_rd",group_vars){
 x%>%group_by(across(all_of(c(group_vars,"method","nuisance_spec","estimator"))))%>%summarise(
 n_total=n(),n_finite=sum(is.finite(estimate)),mean_estimate=mean(estimate,na.rm=TRUE),
 bias=mean(estimate-.data[[true_col]],na.rm=TRUE),mc_sd=sd(estimate,na.rm=TRUE),
 mean_se=mean(se,na.rm=TRUE),se_ratio=mean(se,na.rm=TRUE)/sd(estimate,na.rm=TRUE),
 rmse=sqrt(mean((estimate-.data[[true_col]])^2,na.rm=TRUE)),
 coverage=mean(ci_lower<=.data[[true_col]]&ci_upper>=.data[[true_col]],na.rm=TRUE),
 ci_available=sum(is.finite(ci_lower)&is.finite(ci_upper)),.groups="drop")%>%
 group_by(across(all_of(group_vars)))%>%mutate(full_var=mc_sd[method=="full_aipw"&nuisance_spec=="main"][1]^2,
 relative_efficiency=full_var/(mc_sd^2),variance_inflation=(mc_sd^2)/full_var)%>%ungroup()
}
