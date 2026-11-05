# Manuscript comparisons from validated replication records; no models are fitted.
arg0 <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
source(file.path(dirname(arg0), "../R/project.R"))
source(file.path(root,"R/aggregation.R"))
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
a <- commandArgs(TRUE); run <- if(length(a)) a[1] else bundled_run
if(length(a)>1L || !grepl("^[a-zA-Z0-9_-]+$",run)) stop("Usage: manuscript_figures.R [RUN]")
result_dir <- file.path(root,"outputs",run,"Results")
settings <- readRDS(file.path(result_dir,"settings.rds"))
stopifnot(identical(settings$run,run),settings$mice_rf_trees==500L,
  settings$s6_mice_ci_method=="rubin_t",settings$s6_drcmd_folds=="10")
out <- file.path(root,"outputs",run,"manuscript_figures")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
labels <- c(full_aipw="Complete-data AIPW", cc_aipw="Complete-case AIPW",
  single_imp_aipw="Single imputation + AIPW", mice_aipw="Multiple imputation + AIPW", drcmd="drcmd")
colors <- c(full_aipw="#777777", cc_aipw="#D55E00", single_imp_aipw="#0072B2", mice_aipw="#009E73", drcmd="#AA3377")
shapes <- setNames(c(16,17,15,18,8), names(labels))
inputs <- character()
raw <- bind_rows(lapply(c("outcome_interaction","outcome_interaction_mcar","gamma_y","nonlinear_binary"),function(f) {
  path <- if(f=="nonlinear_binary") file.path(result_dir,f,"combined_with_tr_aipw/all_results.rds") else
    file.path(result_dir,paste0(f,"_combined.rds"))
  inputs <<- c(inputs,path)
  z <- readRDS(path); if(f!="nonlinear_binary") z <- z$results
  validate_rf_mice_settings(z,f,500L)
  validate_mice_intervals(z,f,"rubin_t")
  stopifnot(all(z$drcmd_k[z$method=="drcmd"]==10L))
  mutate(z,family=f)
}))
raw <- raw |> filter(method %in% names(labels))
if (!"status" %in% names(raw)) raw$status <- NA_character_
raw <- raw |> mutate(point_ok=(is.na(status)|status=="ok") & is.finite(estimate),
  interval_ok=point_ok & is.finite(ci_lower) & is.finite(ci_upper) & ci_lower<=ci_upper)
keys <- c("family","scenario_id","method","nuisance_spec","sim")
stopifnot(!anyDuplicated(raw[keys]), all(is.finite(raw$true_rd)))
oc <- raw |> group_by(family,scenario_id,n,target_missing_share,missing_mechanism,gamma_y,gamma_int,method,nuisance_spec,true_rd) |>
  summarise(attempted=n(), valid_point=sum(point_ok), valid_interval=sum(interval_ok),
    complete=setequal(sim,1:500),
    bias=mean(estimate[point_ok]-true_rd[point_ok]),
    rmse=sqrt(mean((estimate[point_ok]-true_rd[point_ok])^2)),
    coverage=mean(ci_lower[interval_ok]<=true_rd[interval_ok] & ci_upper[interval_ok]>=true_rd[interval_ok]),
    .groups="drop") |>
  mutate(point_failure_rate=1-valid_point/attempted, interval_failure_rate=1-valid_interval/attempted)
stopifnot(all(oc$attempted==500), all(oc$complete))
specs <- function(d) mutate(d, Specification=factor(nuisance_spec,levels=c("main","interaction"),
  labels=c("Main-effects GLM","DGP-informed interaction GLM")))
d1 <- oc |> filter(family %in% c("outcome_interaction","outcome_interaction_mcar"), n==2000,
  target_missing_share==.4, gamma_int %in% c(0,.5,1),
  (missing_mechanism=="MCAR" | (missing_mechanism=="MAR" & gamma_y==.5)),
  nuisance_spec %in% c("main","interaction")) |> specs() |>
  mutate(Mechanism=factor(missing_mechanism,levels=c("MCAR","MAR")),x=gamma_int)
d2 <- oc |> filter(family=="gamma_y",n==2000,target_missing_share==.4,gamma_int==0,
  gamma_y %in% c(0,.5,.8),missing_mechanism=="MAR",nuisance_spec=="main") |> mutate(x=gamma_y)
d3 <- oc |> filter(family=="nonlinear_binary", target_missing_share==.4,
  missing_mechanism=="MAR", gamma_y==.5, gamma_int==1, n %in% c(500,1000,2000),
  nuisance_spec %in% c("glm","rf")) |>
  mutate(Specification=factor(nuisance_spec,levels=c("glm","rf"),labels=c("GLM","RF")),x=n,
    mice_ci_method=ifelse(method=="mice_aipw","rubin_t",NA_character_),
    mice_rf_ntree=ifelse(method=="mice_aipw" & nuisance_spec=="rf",500L,NA_integer_))
# Check every plotted method and design cell before exporting figures.
expected <- bind_rows(
  expand_grid(figure=1L, family=c("outcome_interaction","outcome_interaction_mcar"),
    n=2000, gamma_int=c(0,.5,1), nuisance_spec=c("main","interaction"), method=names(labels)),
  expand_grid(figure=2L, family="gamma_y", n=2000, gamma_int=0,
    gamma_y=c(0,.5,.8), nuisance_spec="main", method=names(labels)),
  expand_grid(figure=3L, family="nonlinear_binary", n=c(500,1000,2000), gamma_int=1,
    gamma_y=.5, nuisance_spec=c("glm","rf"), method=names(labels))) |>
  mutate(target_missing_share=.4,
    missing_mechanism=ifelse(family=="outcome_interaction_mcar","MCAR","MAR"),
    gamma_y=ifelse(figure==1, ifelse(missing_mechanism=="MCAR",0,.5),gamma_y))
displayed <- bind_rows(mutate(d1,figure=1L),mutate(d2,figure=2L),mutate(d3,figure=3L))
audit_keys <- c("figure","family","n","gamma_int","gamma_y","nuisance_spec","method",
                "target_missing_share","missing_mechanism")
stopifnot(!anyDuplicated(displayed[audit_keys]))
audit <- expected |> left_join(displayed |> select(all_of(audit_keys),attempted,valid_point,valid_interval,complete,true_rd),by=audit_keys) |>
  mutate(availability=case_when(is.na(attempted) ~ "MISSING",
    !complete | attempted!=500 ~ "INCOMPLETE",
    valid_point==0 | valid_interval==0 ~ "NO VALID METRICS", TRUE ~ "available"))
write.csv(audit,file.path(out,"data_availability.csv"),row.names=FALSE)
missing <- audit |> filter(availability!="available")
write.csv(missing,file.path(out,"missing_combinations.csv"),row.names=FALSE)
if(nrow(missing)) {
  print(missing, n=Inf)
  stop("Required combinations are missing/incomplete; see missing_combinations.csv. No substitutions or simulations performed.")
}
stopifnot(nrow(d1)==60,nrow(d2)==15,nrow(d3)==30)
message("Availability verified: 60 Figure 1, 15 Figure 2, and 30 Figure 3 cells; no missing combinations.")
metrics <- c(bias="Empirical bias",rmse="RMSE",coverage="95% CI coverage")
# Typography and panel dimensions are set for a 6.5-inch manuscript column.
font_sizes <- c(axis_title=14.3, axis_tick=11.44, strip=13, legend=13)
base_plot <- function(d, breaks, ticklabels, xlab, xpad=.09) {
  z <- d |> mutate(method=factor(method,levels=names(labels))) |>
    pivot_longer(all_of(names(metrics)),names_to="Metric",values_to="value") |>
    mutate(Metric=factor(Metric,levels=names(metrics),labels=unname(metrics)))
  refs <- data.frame(Metric=factor(c("Empirical bias","95% CI coverage"),levels=unname(metrics)),ref=c(0,.95))
  anchors <- data.frame(Metric=factor(c("Empirical bias","RMSE","95% CI coverage","95% CI coverage"),levels=unname(metrics)),
    x=min(d$x),value=c(0,0,0,1))
  ggplot(z,aes(x,value,color=method,shape=method,group=method)) +
    geom_blank(data=anchors,aes(x,value),inherit.aes=FALSE) +
    geom_hline(data=refs,aes(yintercept=ref),linetype="dashed",color="grey40",linewidth=.4) +
    geom_line(linewidth=.65) + geom_point(size=2.1) +
    scale_color_manual(values=colors,labels=labels,drop=FALSE) +
    scale_shape_manual(values=shapes,labels=labels,drop=FALSE) +
    scale_x_continuous(breaks=breaks,labels=ticklabels,expand=expansion(mult=c(xpad,xpad))) +
    scale_y_continuous(expand=expansion(mult=c(.05,.07))) +
    labs(x=xlab,y=NULL,color=NULL,shape=NULL) +
    guides(color=guide_legend(nrow=2,byrow=TRUE),shape=guide_legend(nrow=2,byrow=TRUE)) +
    theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(color="grey93"),
      strip.background=element_rect(fill="#F2F4F6"),strip.text=element_text(face="bold",size=font_sizes[["strip"]]),
      legend.position="bottom",legend.text=element_text(size=font_sizes[["legend"]]),
      legend.title=element_text(size=font_sizes[["axis_title"]]),legend.key.width=grid::unit(.65,"cm"),
      axis.title=element_text(size=font_sizes[["axis_title"]]),
      axis.text=element_text(color="#333333",size=font_sizes[["axis_tick"]]),plot.caption.position="plot",
      plot.caption=element_text(hjust=0,size=10,lineheight=1.15,margin=margin(t=10)),
      plot.margin=margin(12,14,12,12))
}
counts <- function(d) {
  if (all(d$valid_point==500 & d$valid_interval==500))
    return("All displayed summaries use 500 valid replications per method and design cell (500 attempted).")
  paste0("500 replications were attempted per method and cell; valid estimates: ",paste(range(d$valid_point),collapse="-"),
    "; valid intervals: ",paste(range(d$valid_interval),collapse="-"),". Exact counts and failures are in displayed_results.csv.")
}
common <- paste("Bias and RMSE use valid point estimates relative to each cell's stored true ATE (marginal risk difference); coverage uses valid 95% intervals.",
  "Complete-data AIPW uses exposure before missingness, with estimated nuisance functions. MR (code identifier drcmd) uses k=10. No finite outliers are removed.")
titles <- c("Performance under MCAR and MAR with varying treatment–covariate interaction",
  "Performance under MAR with varying outcome dependence in exposure observation",
  "Performance under a nonlinear DGP using GLM and random forest nuisance estimation")
captions <- c(
  paste("Empirical bias, RMSE, and coverage of 95% confidence intervals at $n=2000$ and 40% missing exposure.",
    "Results compare main-effects GLM and GLM with the DGP-informed interaction term under MCAR and MAR, with $\\gamma_Y=0.5$ under MAR.",
    "The horizontal axis shows the true treatment–covariate interaction coefficient.",
    "At $\\gamma_{\\mathrm{int}}=0$, the fitted interaction specification retains the interaction term although its true coefficient is zero."),
  paste("Empirical bias, RMSE, and coverage of 95% confidence intervals at $n=2000$, 40% missing exposure, and no true treatment–covariate interaction.",
    "Main-effects GLM nuisance specifications are used throughout. The coefficient $\\gamma_Y$ controls outcome dependence in the exposure observation model.",
    "When $\\gamma_Y=0$, observation depends on covariates only; all settings remain MAR."),
  paste("Empirical bias, RMSE, and coverage of 95% confidence intervals under the nonlinear transformed-covariate DGP with $\\eta=1$, MAR, and 40% missing exposure.",
    "Columns compare GLM and RF specifications, and the horizontal axis shows sample size.",
    "In the RF column, random forests are used for AIPW outcome regression and DRCMD outcome and pseudo-outcome regressions; these regressions use GLMs in the GLM column.",
    "Treatment propensity models and the DRCMD exposure observation model use GLMs in both columns.",
    "For the imputation-based approaches, the GLM column uses logistic exposure imputation followed by AIPW with GLM outcome regression, whereas the RF column uses RF exposure imputation followed by AIPW with RF outcome regression.",
    "All exposure imputation models use covariates and outcome as predictors. Single imputation uses stochastic Bernoulli draws from the predicted exposure probabilities in both columns.",
    "Multiple imputation uses logistic regression in the GLM column and RF with 500 trees in the RF column, with 10 imputations and 10 iterations. Both columns use Rubin t-based 95% intervals.",
    "Thus, for imputation-based approaches the columns change both exposure imputation and AIPW outcome regression; differences should not be attributed solely to outcome regression.",
    "Here, $\\eta$ is the covariate–covariate interaction coefficient, stored as gamma\\_int in the results."))
plots <- list(
  base_plot(d1,c(0,.5,1),c("0","0.5","1"),expression(gamma[int]),xpad=.17) +
    facet_grid(Metric~Mechanism+Specification,scales="free_y",switch="y"),
  base_plot(d2,c(0,.5,.8),c("0","0.5","0.8"),expression(gamma[Y])) +
    facet_wrap(~Metric,nrow=1,scales="free_y"),
  base_plot(d3,c(500,1000,2000),c("500","1000","2000"),"Sample size") +
    facet_grid(Metric~Specification,scales="free_y",switch="y"))
stems <- c("Performance_under_MCAR_and_MAR_with_varying_treatment_covariate_interaction",
  "Performance_under_MAR_with_varying_outcome_dependence",
  "Performance_under_a_nonlinear_DGP_using_GLM_and_RF")
widths <- c(16,15,10); heights <- c(8,5.5,8)
# Target sizes are measured after inclusion at manuscript_width inches.
# Preserve panel dimensions while allocating space for axis and facet labels.
manuscript_width <- as.numeric(Sys.getenv("MANUSCRIPT_WIDTH_IN", "6.5"))
stopifnot(is.finite(manuscript_width), manuscript_width > 0)
final_fonts <- c(axis_title=10.5, axis_tick=9.5, strip=10.5, legend=9.5)
layout_audit <- list()
compact_labels <- c("Complete-data\nAIPW", "Complete-case\nAIPW",
                    "Single imputation\n+ AIPW", "Multiple imputation\n+ AIPW", "drcmd")
for(i in 1:3) {
  # Measure the pre-revision plotting area on the same vector device.
  grDevices::cairo_pdf(tempfile(fileext=".pdf"), width=widths[i], height=heights[i])
  baseline <- ggplotGrob(plots[[i]])
  panel_rows <- unique(baseline$layout$t[grepl("^panel",baseline$layout$name)])
  panel_cols <- unique(baseline$layout$l[grepl("^panel",baseline$layout$name)])
  panel_width <- widths[i] - grid::convertWidth(sum(baseline$widths),"in",valueOnly=TRUE)
  panel_height <- heights[i] - grid::convertHeight(sum(baseline$heights),"in",valueOnly=TRUE)
  p <- plots[[i]] + scale_color_manual(values=colors,labels=compact_labels,drop=FALSE) +
    scale_shape_manual(values=shapes,labels=compact_labels,drop=FALSE)
  if(i==1) p <- p + scale_x_continuous(breaks=c(0,.5,1),
    labels=c("0", "0.5", "1"),
    expand=expansion(mult=c(.17,.17))) +
    facet_grid(Metric~Mechanism+Specification,scales="free_y",switch="y",
      labeller=labeller(Specification=c("Main-effects GLM"="Main-effects\nGLM",
        "DGP-informed interaction GLM"="DGP-informed\ninteraction GLM"),
        Metric=c("Empirical bias"="Empirical\nbias", "RMSE"="RMSE", "95% CI coverage"="95% CI\ncoverage")))
  if(i==2) p <- p + scale_x_continuous(breaks=c(0,.5,.8),
    labels=c("0", "0.5", "0.8"),
    expand=expansion(mult=c(.09,.09)))
  # Solve for the export width: larger text changes the space around panels.
  export_width <- widths[i]
  for(iter in 1:12) {
    fs <- final_fonts * export_width / manuscript_width
    q <- p + theme(axis.title=element_text(size=fs[["axis_title"]]),
      axis.text=element_text(size=fs[["axis_tick"]],lineheight=1),
      strip.text=element_text(size=fs[["strip"]]),
      legend.text=element_text(size=fs[["legend"]],lineheight=1),
      legend.title=element_text(size=fs[["axis_title"]]),
      plot.margin=margin(12,if(i==2) 3*fs[["axis_tick"]] else if(i==1) 1.8*fs[["axis_tick"]] else 14,12,12)) +
      labs(title=NULL,subtitle=NULL,caption=NULL)
    g <- ggplotGrob(q)
    cols <- unique(g$layout$l[grepl("^panel",g$layout$name)])
    rows <- unique(g$layout$t[grepl("^panel",g$layout$name)])
    g$widths[cols] <- grid::unit(rep(panel_width/length(cols),length(cols)),"in")
    g$heights[rows] <- grid::unit(rep(panel_height/length(rows),length(rows)),"in")
    new_width <- grid::convertWidth(sum(g$widths),"in",valueOnly=TRUE)
    export_height <- grid::convertHeight(sum(g$heights),"in",valueOnly=TRUE)
    if(abs(new_width-export_width)<1e-6) break
    export_width <- new_width
  }
  dev.off()
  ggsave(file.path(out,paste0(stems[i],".pdf")),g,width=new_width,height=export_height,device=cairo_pdf)
  ggsave(file.path(out,paste0(stems[i],".png")),g,width=new_width,height=export_height,dpi=200,device=ragg::agg_png)
  layout_audit[[i]] <- data.frame(figure=stems[i],manuscript_width_in=manuscript_width,
    export_width_in=new_width,export_height_in=export_height,
    panel_width_in=panel_width/length(cols),panel_height_in=panel_height/length(rows),
    axis_title_pt=fs[["axis_title"]],tick_pt=fs[["axis_tick"]],
    facet_pt=fs[["strip"]],legend_pt=fs[["legend"]])
  message("Saved ",stems[i]," at ",round(new_width,3)," x ",round(export_height,3)," inches")
}
write.csv(bind_rows(layout_audit),file.path(out,"font_layout_audit.csv"),row.names=FALSE)
write.csv(bind_rows(mutate(d1,figure=1),mutate(d2,figure=2),mutate(d3,figure=3)),file.path(out,"displayed_results.csv"),row.names=FALSE)
writeLines(unlist(lapply(1:3,function(i)c(paste0("% Figure ",i),
  "\\begin{figure}[htbp]", "\\centering",paste0("\\includegraphics[width=\\linewidth]{",stems[i],".pdf}"),
  paste0("\\caption{\\textbf{",titles[i],".} ",gsub("%","\\%",captions[i],fixed=TRUE),"}"),
  paste0("\\label{fig:missing-exposure-",i,"}"),"\\end{figure}",""))),file.path(out,"overleaf_figures.tex"))
write.csv(data.frame(figure=1:3,title=titles,caption=captions),file.path(out,"captions.csv"),row.names=FALSE)
writeLines(unlist(lapply(1:3,function(i)c(titles[i],captions[i],""))),file.path(out,"captions.txt"))

for (i in 1:3) {
  z <- displayed[displayed$figure==i,]
  write.csv(z,file.path(out,paste0("figure_",i,"_summary.csv")),row.names=FALSE)
}
write.csv(data.frame(file=basename(inputs),md5=unname(tools::md5sum(inputs))),
  file.path(out,"input_checksums.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"session_info.txt"))
