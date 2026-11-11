# gamma_Y,AX and gamma_Y,Z1Z2 are stored as interaction, with dgp identifying which.
study_design <- function() {
  main <- data.frame(
    cell = c("mcar_0", "mar_x_0", "mar_xy_0", "mnar_0",
             "mcar_05", "mcar_1", "mar_xy_05", "mar_xy_1", "mar_xy2_0"),
    dgp = "gaussian", n = 2000L,
    interaction = c(0, 0, 0, 0, 0.5, 1, 0.5, 1, 0),
    delta_r = c(0, 1, 1, 1, 0, 0, 1, 1, 1),
    gamma_ry = c(0, 0, 0.5, 0.5, 0, 0, 0.5, 0.5, 2),
    gamma_ra = c(0, 0, 0, 2, 0, 0, 0, 0, 0),
    missing_share = 0.4, replications = 500L,
    seed_cell = c(3L,3L,3L,6L,6L,9L,6L,9L,12L))
  mechanism <- main[main$cell == "mar_xy_0", ]
  mechanism$cell <- "mar_xy_0_mechanism"
  mechanism$seed_cell <- 6L
  main <- rbind(main, mechanism)
  nonlinear <- expand.grid(n = c(500L, 1000L, 2000L), interaction = c(0, 0.5, 1))
  nonlinear$cell <- paste0("nonlinear_", gsub("[.]", "", nonlinear$interaction), "_", nonlinear$n)
  nonlinear$dgp <- "nonlinear"
  nonlinear$delta_r <- 1
  nonlinear$gamma_ry <- 0.5
  nonlinear$gamma_ra <- 0
  nonlinear$missing_share <- 0.4
  nonlinear$replications <- 500L
  nonlinear$seed_cell <- seq_len(nrow(nonlinear))
  rbind(main, nonlinear[names(main)])
}

study_methods <- c("complete_data", "complete_case", "single_imputation", "multiple_imputation", "mr")

study_analyses <- function(design = study_design()) {
  rows <- lapply(seq_len(nrow(design)), function(i) {
    cell <- design[i, ]
    specs <- if (cell$dgp == "nonlinear") c("glm", "rf") else "main"
    if (cell$cell %in% c("mcar_0", "mcar_05", "mcar_1", "mar_xy_0", "mar_xy_05", "mar_xy_1"))
      specs <- c(specs, "interaction")
    a <- expand.grid(method = study_methods, specification = specs, stringsAsFactors = FALSE)
    a$include_y <- TRUE
    if (cell$cell %in% c("mar_x_0", "mar_xy2_0")) {
      b <- data.frame(method = c("single_imputation", "multiple_imputation", "mr"),
                      specification = "main", include_y = FALSE)
      a <- rbind(a, b)
    }
    a$cell <- cell$cell
    a
  })
  result <- do.call(rbind, rows)
  result$analysis <- paste(result$cell, result$method, result$specification,
                           ifelse(result$include_y, "with_y", "without_y"), sep = "__")
  rownames(result) <- NULL
  result
}

# Repeated reference columns map to the same analysis instead of fitting twice.
figure_design <- function(number, analyses = study_analyses(), design = study_design()) {
  x <- merge(analyses, design, by = "cell", sort = FALSE)
  if (number == 2L) {
    cells <- c("mcar_0", "mar_x_0", "mar_xy_0_mechanism", "mnar_0")
    x <- x[x$cell %in% cells & x$specification == "main" & x$include_y, ]
    x$panel <- "Observation mechanism"
    x$x <- match(x$cell, cells)
  } else if (number == 3L) {
    cells <- c("mcar_0", "mcar_05", "mcar_1", "mar_xy_0", "mar_xy_05", "mar_xy_1")
    x <- x[x$cell %in% cells & x$include_y, ]
    x$panel <- paste(ifelse(x$delta_r == 0, "MCAR", "MAR"), x$specification, sep = "__")
    x$x <- x$interaction
  } else if (number == 4L) {
    x <- x[x$cell %in% c("mar_x_0", "mar_xy2_0") & x$specification == "main", ]
    refs <- x[x$method %in% c("complete_data", "complete_case"), ]
    refs$include_y <- FALSE
    x <- rbind(x, refs)
    x$panel <- paste(x$cell, ifelse(x$include_y, "include", "exclude"), sep = "__")
    x$x <- 1
  } else if (number == 5L) {
    x <- x[x$dgp == "nonlinear" & x$interaction == 1, ]
    x$panel <- x$specification
    x$x <- x$n
  } else stop("Choose a performance figure from 2 to 5.")
  x$figure <- number
  x
}
