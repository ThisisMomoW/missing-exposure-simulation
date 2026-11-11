# Run once at the configured sample sizes; this checks execution, not performance.
design <- study_design();analyses <- study_analyses()
selected <- c(
  'mcar_0__complete_data__main__with_y','mcar_0__complete_case__main__with_y',
  'mcar_0__single_imputation__main__with_y','mcar_0__mr__main__with_y',
  'mar_xy_1__complete_data__interaction__with_y','mar_xy_1__complete_case__interaction__with_y',
  'mar_xy_1__single_imputation__interaction__with_y','mar_xy_1__mr__interaction__with_y',
  'mar_xy2_0__single_imputation__main__without_y','mar_xy2_0__mr__main__without_y',
  'nonlinear_1_500__complete_data__rf__with_y','nonlinear_1_500__complete_case__rf__with_y',
  'nonlinear_1_500__single_imputation__rf__with_y','nonlinear_1_500__mr__rf__with_y')
folder <- file.path(root,'outputs/checks/execution')
dir.create(folder,recursive=TRUE,showWarnings=FALSE)
for(id in selected) {
  a <- analyses[analyses$analysis==id,];d<-design[design$cell==a$cell,]
  p<-prepare_cell(d)
  result<-fit_replication(p,a,1L)
  atomic_rds(result,file.path(folder,paste0(id,'.rds')))
  message(id,': ',result$status)
  if(!result$valid)stop(result$error)
}

mice_selected <- c('mcar_0__multiple_imputation__main__with_y',
  'mar_xy_1__multiple_imputation__interaction__with_y',
  'mar_xy2_0__multiple_imputation__main__with_y',
  'mar_xy2_0__multiple_imputation__main__without_y',
  'nonlinear_1_500__multiple_imputation__glm__with_y',
  'nonlinear_1_500__multiple_imputation__rf__with_y')
.rf_check <- new.env(parent=emptyenv());.rf_check$inside<-FALSE;.rf_check$trees<-integer()
trace('.ranger.donors',where=asNamespace('mice'),print=FALSE,
 tracer=quote({stopifnot(ntree==500L);.GlobalEnv$.rf_check$inside<-TRUE}),
 exit=quote(.GlobalEnv$.rf_check$inside<-FALSE))
trace('ranger',where=asNamespace('ranger'),print=FALSE,
 exit=quote(if(.GlobalEnv$.rf_check$inside) {
   forest<-returnValue();stopifnot(forest$num.trees==500L)
   .GlobalEnv$.rf_check$trees<-c(.GlobalEnv$.rf_check$trees,forest$num.trees)
 }))
paired <- NULL
for(id in mice_selected) {
  a<-analyses[analyses$analysis==id,];d<-design[design$cell==a$cell,];p<-prepare_cell(d)
  seeds<-replication_seeds(d,a,1L);set.seed(seeds$data);dat<-generate_data(d,p$calibration$intercept);set.seed(seeds$fit)
  result<-fit_multiple_imputation(dat,a$specification,a$include_y,seeds$imputation,
    shared_rng=if(a$include_y)NULL else paired)
  atomic_rds(result,file.path(folder,paste0(id,'.rds')))
  stopifnot(valid_estimate(result$value),result$diagnostics$mice_valid==10L)
  if(a$cell=='mar_xy2_0'&&a$include_y)paired<-result$diagnostics$component_rng
  message(id,': ',result$diagnostics$mice_pooling_status)
}
untrace('ranger',where=asNamespace('ranger'))
untrace('.ranger.donors',where=asNamespace('mice'))
stopifnot(length(.rf_check$trees)==100L,all(.rf_check$trees==500L))
write.csv(data.frame(forest=1:100,trees=.rf_check$trees),file.path(folder,'rf_mice_trees.csv'),row.names=FALSE)
message('PASS: configured estimator branches, paired MICE folds and 100 RF-MICE forests with 500 trees each.')
