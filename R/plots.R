method_labels <- c(complete_data = "Complete-data\nAIPW", complete_case = "Complete-case\nAIPW",
  single_imputation = "Single imputation\n+ AIPW", multiple_imputation = "Multiple imputation\n+ AIPW", mr = "MR")
method_colors <- c(complete_data = "#777777", complete_case = "#D55E00",
  single_imputation = "#0072B2", multiple_imputation = "#009E73", mr = "#AA3377")

figure_dimensions <- function(number) switch(as.character(number),
  `1` = c(10, 1.8), `2` = c(9, 8), `3` = c(17.6397148164784, 10.6619961847504),
  `4` = c(9, 8), `5` = c(10.1666666666667, 8.58474358974359))

# Offsets depend on method identity, including when a method has no valid result.
performance_figure <- function(data, number) {
  if (!number %in% 2:5 || anyDuplicated(data[c("analysis","panel")])) stop("Invalid figure selection")
  data$Method <- factor(data$method, levels = names(method_labels))
  if (anyNA(data$Method)) stop("Unexpected method")
  offset <- switch(as.character(number), `2`=.048, `3`=.032, `4`=.048, `5`=30)
  if (number == 4L) data$x <- match(data$panel,
    c("mar_x_0__exclude","mar_x_0__include","mar_xy2_0__exclude","mar_xy2_0__include"))
  data$display_x <- data$x + (as.integer(data$Method)-3)*offset
  if (number == 4L) data$display_x <- 1+(as.integer(data$Method)-3)*offset
  metrics <- c(relative_bias="Relative\nbias (%)", rmse="RMSE", coverage="95% CI\ncoverage")
  z <- do.call(rbind, lapply(names(metrics), function(k) {
    a <- data; a$Metric <- k; a$value <- a[[k]]; a
  }))
  z$Metric <- factor(z$Metric, levels=names(metrics), labels=unname(metrics))
  refs <- data.frame(Metric=factor(unname(metrics[c(1,3)]), levels=unname(metrics)), value=c(0,.95))
  p <- ggplot2::ggplot(z, ggplot2::aes(display_x, value, color=Method, shape=Method, group=Method)) +
    ggplot2::geom_hline(data=refs, ggplot2::aes(yintercept=value), linetype="dashed", color="grey50", linewidth=.4)
  if (number %in% c(3L,5L)) p <- p + ggplot2::geom_line(linewidth=.65, na.rm=TRUE)
  shapes <- setNames(if (number %in% c(2L,4L)) c(16,17,15,3,7) else c(16,17,15,18,8), names(method_labels))
  p <- p + ggplot2::geom_point(size=2.1, na.rm=TRUE) +
    ggplot2::scale_color_manual(values=method_colors, labels=method_labels, drop=FALSE) +
    ggplot2::scale_shape_manual(values=shapes, labels=method_labels, drop=FALSE) +
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.05,.07))) +
    ggplot2::labs(x=NULL,y=NULL,color=NULL,shape=NULL) +
    ggplot2::theme_bw(base_size=11) + ggplot2::theme(
      panel.grid.minor=ggplot2::element_blank(), panel.grid.major=ggplot2::element_line(color="grey93"),
      strip.background=ggplot2::element_rect(fill="#F2F4F6"), strip.text=ggplot2::element_text(face="bold"),
      legend.position="bottom", legend.text=ggplot2::element_text(lineheight=1),
      plot.margin=ggplot2::margin(8,10,8,8)) +
    ggplot2::guides(color=ggplot2::guide_legend(nrow=2,byrow=TRUE), shape=ggplot2::guide_legend(nrow=2,byrow=TRUE))
  if (number == 2L) p <- p + ggplot2::facet_grid(Metric ~ ., scales="free_y") +
    ggplot2::scale_x_continuous(breaks=1:4, limits=c(.4,4.6), expand=c(0,0), labels=expression(
      bold(MCAR), bold("MAR given X"), atop(bold("MAR given X,Y"), gamma[R*","*Y]==.5),
      atop(bold(MNAR),gamma[R*","*A]==2)))
  if (number == 4L) {
    z$Mechanism <- factor(z$cell, levels=c("mar_x_0","mar_xy2_0"),
      labels=c("atop(bold('MAR given X'),gamma[R*','*Y]==0)","atop(bold('MAR given X,Y'),gamma[R*','*Y]==2)"))
    z$Predictors <- factor(z$include_y, levels=c(FALSE,TRUE), labels=c("Exclude Y","Include Y"))
    # Each panel uses the same local axis so categorical columns cannot imply a trend.
    p$data <- z
    p <- p + ggplot2::facet_grid(Metric ~ Mechanism+Predictors, scales="free_y", switch="y",
      labeller=ggplot2::labeller(Mechanism=ggplot2::label_parsed)) +
      ggplot2::scale_x_continuous(breaks=NULL, limits=c(.4,1.6),expand=c(0,0))
  }
  if (number %in% c(3L,5L)) {
    anchors <- data.frame(Metric=factor(unname(metrics[c(1,2,3,3)]),levels=unname(metrics)),
                          display_x=min(data$x),value=c(0,0,0,1))
    p <- p + ggplot2::geom_blank(data=anchors, ggplot2::aes(display_x,value),inherit.aes=FALSE)
    if (number == 3L) {
      z$Mechanism <- factor(ifelse(z$delta_r==0,"MCAR","MAR"),levels=c("MCAR","MAR"))
      z$Specification <- factor(z$specification, levels=c("main","interaction"),
        labels=c("Main-effects\nGLM","DGP-informed\ninteraction GLM"))
      p$data <- z
      p <- p + ggplot2::facet_grid(Metric ~ Mechanism+Specification,scales="free_y",switch="y") +
        ggplot2::scale_x_continuous(breaks=c(0,.5,1),labels=c("0","0.5","1"),expand=ggplot2::expansion(mult=.17)) +
        ggplot2::labs(x=expression(gamma[Y*","*AX]))
      facet_size <- 28.4949233876708; tick_size <- 25.7811211602736
    } else {
      z$Specification <- factor(z$specification,levels=c("glm","rf"),labels=c("GLM","RF"))
      p$data <- z
      p <- p + ggplot2::facet_grid(Metric ~ Specification,scales="free_y",switch="y") +
        ggplot2::scale_x_continuous(breaks=c(500,1000,2000),labels=c("500","1000","2000"),expand=ggplot2::expansion(mult=.09)) +
        ggplot2::labs(x="Sample size")
      facet_size <- 16.4230769230769; tick_size <- 14.8589743589744
    }
    p <- p + ggplot2::theme(strip.text=ggplot2::element_text(size=facet_size,face="bold"),
      strip.text.y=ggplot2::element_text(size=.78*facet_size), axis.title=ggplot2::element_text(size=facet_size),
      axis.text=ggplot2::element_text(size=tick_size,color="#333333"),
      legend.text=ggplot2::element_text(size=tick_size),legend.key.width=grid::unit(.65,"cm"),
      plot.margin=ggplot2::margin(12, if(number==3L)1.8*tick_size else 14,12,12))
  }
  list(plot=p,data=data)
}

mechanism_diagram <- function(file) {
  grDevices::pdf(file,width=10,height=1.8,useDingbats=FALSE)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow=c(1,4),mar=c(0,0,0,0))
  for (i in 1:4) {
    graphics::plot.new(); graphics::plot.window(xlim=c(0,1),ylim=c(0,1))
    graphics::text(c(.14,.48,.82,.82),c(.64,.64,.64,.27),c("X","A","Y","R"),font=3,cex=1.2)
    arrow <- function(x0,y0,x1,y1,lty=1) graphics::arrows(x0,y0,x1,y1,length=.07,lwd=1,lty=lty)
    arrow(.21,.64,.41,.64); arrow(.55,.64,.75,.64)
    curve <- function(x0,y0,x1,y1,bend,lty=1) {
      t <- seq(0,1,length.out=80)
      xc <- (1-t)^2*x0+2*(1-t)*t*(x0+x1)/2+t^2*x1
      yc <- (1-t)^2*y0+2*(1-t)*t*bend+t^2*y1
      graphics::lines(xc,yc,lwd=1,lty=lty); arrow(xc[77],yc[77],xc[80],yc[80],lty)
    }
    curve(.18,.73,.77,.73,1)
    if (i>1) curve(.18,.55,.75,.27,.3,if(i==4)2 else 1)
    if (i>2) arrow(.82,.56,.82,.35,if(i==4)2 else 1)
    if (i==4) arrow(.53,.57,.76,.33)
    graphics::text(.5,.05,c("(a) MCAR","(b) MAR given X","(c) MAR given X,Y","(d) MNAR")[i],cex=1)
  }
  invisible(file)
}
