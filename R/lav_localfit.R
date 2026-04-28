#' Local fit diagnostics for fitted \code{lavaan} models
#'
#' Evaluate local fit in a fitted \code{lavaan} model from residual-based
#' diagnostics. The function extracts residual summaries, computes standardized 
#' covariance and mean residual diagnostics, identifies the largest local 
#' discrepancies, and builds a heatmap of covariance residual misfit. 
#' By default, the function computes Bentler-standardized residual summaries and
#' standardized covariance residual diagnostics suitable for screening localized
#' strain in the model. 
#' The default plot is a traffic-light heatmap of the absolute standardized
#' covariance residuals.
#'
#' @usage
#' lav_localfit(
#'   fit,
#'   group = 1L,
#'   type = "cor.bentler",
#'   thresholds = c(1.96, 2.58),
#'   top_n = 10L,
#'   triangle = "lower",
#'   include_diagonal = FALSE,
#'   plot = TRUE,
#'   plot_style = "trafficlight",
#'   show_values = FALSE,
#'   digits = 3L,
#'   good_color = "#009E73",
#'   moderate_color = "#F0E442",
#'   poor_color = "#D55E00",
#'   neg_color = "#3B4CC0",
#'   pos_color = "#B40426",
#'   na_color = "grey90",
#'   return_data = FALSE
#' )
#'
#' @param fit A fitted \code{lavaan} object.
#' @param group Group to inspect in a multigroup model. May be either a numeric
#'   group index or a character group label. Ignored for single-group models.
#' @param type Character string passed to \code{lavaan::lavResiduals()}.
#'   Common choices are \code{"cor.bentler"} (default), \code{"raw"},
#'   \code{"cor"}, and \code{"cor.bollen"}.
#' @param thresholds Numeric vector of length 2 giving the descriptive cutoffs
#'   used for standardized residual magnitudes. The first value separates
#'   \code{"minor"} from \code{"moderate"} discrepancies; the second separates
#'   \code{"moderate"} from \code{"notable"} discrepancies. Used for the
#'   traffic-light heatmap and for printed summaries. Defaults to
#'   \code{c(1.96, 2.58)}, as per common z-based guidelines.
#' @param top_n Integer. Number of largest absolute standardized covariance
#'   residuals to report in the output table.
#' @param triangle Character string. Either \code{"lower"} or \code{"full"}.
#'   If \code{"lower"}, only the lower triangle of the covariance residual
#'   matrix is used in the heatmap and top-residual table.
#' @param include_diagonal Logical. If \code{TRUE}, include diagonal elements
#'   (variance residuals) in the covariance residual diagnostics.
#' @param plot Logical. If \code{TRUE}, produce a heatmap of covariance residual
#'   misfit.
#' @param plot_style Character string. Either \code{"trafficlight"} or
#'   \code{"signed"}. The default \code{"trafficlight"} plots the absolute
#'   standardized covariance residuals in descriptive magnitude categories. The
#'   \code{"signed"} option plots signed standardized residuals on a diverging
#'   color scale.
#' @param show_values Logical. If \code{TRUE}, overlay numeric residual values on
#'   the heatmap tiles.
#' @param digits Non-negative integer giving the default number of digits used
#'   in printed output and numeric labels.
#' @param good_color Fill color for \code{"minor"} discrepancies in the
#'   traffic-light heatmap.
#' @param moderate_color Fill color for \code{"moderate"} discrepancies in the
#'   traffic-light heatmap.
#' @param poor_color Fill color for \code{"notable"} discrepancies in the
#'   traffic-light heatmap.
#' @param neg_color Fill color for negative residuals in the signed heatmap.
#' @param pos_color Fill color for positive residuals in the signed heatmap.
#' @param na_color Fill color for cells omitted from the heatmap, such as the
#'   upper triangle when \code{triangle = "lower"}.
#' @param return_data Logical. If \code{TRUE}, include the long-format heatmap
#'   data frame in the returned object.
#'
#' @return A list of class \code{"lav_localfit"} with elements:
#' \itemize{
#'   \item \code{summary}: residual summary returned by
#'     \code{lavaan::lavResiduals()} for the selected group;
#'   \item \code{cov_residuals}: covariance residual matrix for the selected
#'     group;
#'   \item \code{cov_z}: standardized covariance residual matrix;
#'   \item \code{mean_residuals}: mean residual vector, if available;
#'   \item \code{mean_z}: standardized mean residual vector, if available;
#'   \item \code{top_cov_residuals}: data frame of the largest absolute
#'     standardized covariance residuals;
#'   \item \code{top_mean_residuals}: data frame of the largest absolute
#'     standardized mean residuals;
#'   \item \code{counts}: list with counts of covariance residuals exceeding the
#'     descriptive thresholds, the maximum absolute standardized covariance
#'     residual, and the corresponding variable pair;
#'   \item \code{plot}: a \code{ggplot} heatmap object, or \code{NULL} if
#'     \code{plot = FALSE};
#'   \item \code{group}: numeric index of the selected group;
#'   \item \code{group_label}: character label of the selected group;
#'   \item \code{type}: residual type used;
#'   \item \code{thresholds}: thresholds used for descriptive classification;
#'   \item \code{digits}: default number of digits for printing;
#'   \item \code{call}: matched function call;
#'   \item \code{plot_data}: only when \code{return_data = TRUE}, the long-format
#'     data used to build the heatmap.
#' }
#'
#' @section Notes:
#' In this function, local fit is evaluated from residuals of the sample summary
#' statistics rather than from casewise prediction errors. Standardized residuals
#' are used to screen where the model reproduces specific covariances or means
#' poorly. The default thresholds are descriptive aids only. They should not be
#' treated as formal decision rules or as substitutes for theory-based model
#' evaluation. When the fitted object uses bootstrap standard errors,
#' standardized residual diagnostics are obtained from an auxiliary refit of the
#' same model without bootstrap.
#'
#' @references
#' Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
#' \emph{Journal of Statistical Software}, \emph{48}(2), 1-36.
#' \doi{10.18637/jss.v048.i02}
#'
#' Schermelleh-Engel, K., Moosbrugger, H., & Muller, H. (2003). Evaluating the
#' fit of structural equation models: Tests of significance and descriptive
#' goodness-of-fit measures. \emph{Methods of Psychological Research Online},
#' \emph{8}(2), 23-74.
#'
#' @seealso \code{\link[lavaan]{lavResiduals}} for residual extraction and
#'   \code{\link{lav_deltaR2}} for local structural effect-size diagnostics.
#'
#' @examples
#' model <- "
#' ind60 =~ x1 + x2 + x3
#' dem60 =~ y1 + y2 + y3 + y4
#' dem65 =~ y5 + y6 + y7 + y8
#' dem60 ~ ind60
#' dem65 ~ ind60 + dem60
#' y1 ~~ y5
#' y2 ~~ y6
#' "
#' fit <- lavaan::sem(
#' model = model, 
#' data = lavaan::PoliticalDemocracy,
#' std.lv = TRUE, 
#' estimator = "MLR", 
#' meanstructure = TRUE)
#'
#' lav_localfit(fit)
#'
#' @importFrom lavaan lavInspect lavResiduals
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_manual
#'   scale_fill_gradient2 coord_equal labs theme_minimal theme element_text
#'   element_blank
#' @importFrom rlang %||%
#' @importFrom utils head
#' @export
lav_localfit <- function(
    fit,
    group = 1L,
    type = "cor.bentler",
    thresholds = c(1.96, 2.58),
    top_n = 10L,
    triangle = "lower",
    include_diagonal = FALSE,
    plot = TRUE,
    plot_style = "trafficlight",
    show_values = FALSE,
    digits = 3L,
    good_color = "#009E73",
    moderate_color = "#F0E442",
    poor_color = "#D55E00",
    neg_color = "#3B4CC0",
    pos_color = "#B40426",
    na_color = "grey90",
    return_data = FALSE
) {
  
  if (!inherits(fit, "lavaan")) {
    stop("'fit' must be a fitted 'lavaan' object.", call. = FALSE)
  }
  
  fit_converged <- tryCatch(lavaan::lavInspect(fit, "converged"),
                            error = function(e) FALSE)
  if (!isTRUE(fit_converged)) {
    stop("'fit' did not converge.", call. = FALSE)
  }
  
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("'type' must be a single character string.", call. = FALSE)
  }
  
  triangle <- match.arg(triangle, choices = c("lower", "full"))
  plot_style <- match.arg(plot_style, choices = c("trafficlight", "signed"))
  
  if (!is.numeric(thresholds) || length(thresholds) != 2L || any(is.na(thresholds))) {
    stop("'thresholds' must be a numeric vector of length 2.", call. = FALSE)
  }
  thresholds <- as.numeric(thresholds)
  if (thresholds[1L] <= 0 || thresholds[2L] <= thresholds[1L]) {
    stop("'thresholds' must be positive and strictly increasing.", call. = FALSE)
  }
  
  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) {
    stop("'top_n' must be a positive integer.", call. = FALSE)
  }
  top_n <- as.integer(top_n)
  
  if (!is.logical(include_diagonal) || length(include_diagonal) != 1L || is.na(include_diagonal)) {
    stop("'include_diagonal' must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("'plot' must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(show_values) || length(show_values) != 1L || is.na(show_values)) {
    stop("'show_values' must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0) {
    stop("'digits' must be a non-negative integer.", call. = FALSE)
  }
  digits <- as.integer(digits)
  
  ngroups <- tryCatch(lavaan::lavInspect(fit, "ngroups"),
                      error = function(e) 1L)
  group_labels <- tryCatch(lavaan::lavInspect(fit, "group.label"),
                           error = function(e) NULL)
  
  if (is.null(group_labels) || length(group_labels) != ngroups) {
    group_labels <- as.character(seq_len(ngroups))
  } else {
    group_labels <- as.character(group_labels)
  }
  
  if (ngroups == 1L) {
    group_index <- 1L
    group_label <- group_labels[1L]
  } else {
    if (is.numeric(group) && length(group) == 1L && !is.na(group)) {
      group_index <- as.integer(group)
      if (group_index < 1L || group_index > ngroups) {
        stop("'group' is outside the available group range.", call. = FALSE)
      }
      group_label <- group_labels[group_index]
    } else if (is.character(group) && length(group) == 1L && !is.na(group)) {
      pos <- match(group, group_labels)
      if (is.na(pos)) {
        stop("'group' does not match any available group label.", call. = FALSE)
      }
      group_index <- pos
      group_label <- group_labels[group_index]
    } else {
      stop("'group' must be a single numeric index or a single group label.", call. = FALSE)
    }
  }
  
  fmt_num <- function(x, k = digits, apa = FALSE) {
    if (length(x) == 0L || is.na(x)) return("NA")
    out <- formatC(x, digits = k, format = "f")
    if (apa) {
      out <- sub("^(-?)0\\.", "\\1.", out)
    }
    out
  }
  
  mag_label <- function(abs_z, th = thresholds) {
    if (is.na(abs_z)) {
      return(NA_character_)
    }
    if (abs_z < th[1L]) {
      return("minor")
    }
    if (abs_z < th[2L]) {
      return("moderate")
    }
    "notable"
  }
  
  residual_type_label <- function(x) {
    if (identical(x, "cor.bentler")) return("Bentler-standardized residuals")
    if (identical(x, "cor.bollen")) return("Bollen-standardized residuals")
    if (identical(x, "cor")) return("correlation residuals")
    if (identical(x, "raw")) return("raw residuals")
    paste(x, "residuals")
  }
  
  extract_main_metric <- function(summary_obj, type) {
    if (is.null(summary_obj)) {
      return(list(name = NA_character_, value = NA_real_))
    }
    
    sm <- as.matrix(summary_obj)
    rn <- rownames(sm)
    cn <- colnames(sm)
    
    target_name <- if (identical(type, "cor.bentler")) {
      "srmr"
    } else if (identical(type, "raw")) {
      "rmr"
    } else if (identical(type, "cor.bollen")) {
      "crmr"
    } else {
      NA_character_
    }
    
    if (!is.na(target_name) && target_name %in% rn) {
      target_col <- if ("total" %in% cn) "total" else cn[length(cn)]
      return(list(name = toupper(target_name), value = as.numeric(sm[target_name, target_col])))
    }
    
    if (length(rn) > 0L && length(cn) > 0L) {
      target_col <- if ("total" %in% cn) "total" else cn[length(cn)]
      return(list(name = rn[1L], value = as.numeric(sm[1L, target_col])))
    }
    
    list(name = NA_character_, value = NA_real_)
  }
  
  recover_fit_data <- function(fit) {
    out <- tryCatch(lavaan::lavInspect(fit, "data"), error = function(e) NULL)
    if (is.null(out)) {
      out <- tryCatch(lavaan::lavInspect(fit, "data.original"), error = function(e) NULL)
    }
    if (is.null(out)) {
      stop("Could not recover data from `fit` for auxiliary residual refit.", call. = FALSE)
    }
    
    ov_names <- tryCatch(lavaan::lavNames(fit, type = "ov"), error = function(e) NULL)
    
    if (is.data.frame(out)) {
      out <- as.data.frame(out)
    } else if (is.matrix(out)) {
      out <- as.data.frame(out)
      if (!is.null(ov_names) && length(ov_names) == ncol(out)) {
        names(out) <- ov_names
      }
    } else if (is.list(out)) {
      ok_list <- all(vapply(out, function(x) is.data.frame(x) || is.matrix(x), logical(1)))
      if (!ok_list) {
        stop("Recovered grouped data could not be coerced to a data frame.", call. = FALSE)
      }
      
      group_labels <- tryCatch(lavaan::lavInspect(fit, "group.label"), error = function(e) NULL)
      group_var <- tryCatch(fit@Options$group, error = function(e) NULL)
      if (length(group_var) != 1L || is.na(group_var) || !nzchar(group_var)) {
        group_var <- NULL
      }
      
      out_list <- vector("list", length(out))
      for (g in seq_along(out)) {
        dg <- out[[g]]
        if (is.matrix(dg)) {
          dg <- as.data.frame(dg)
          if (!is.null(ov_names) && length(ov_names) == ncol(dg)) {
            names(dg) <- ov_names
          }
        } else {
          dg <- as.data.frame(dg)
        }
        
        if (!is.null(group_var) && !(group_var %in% names(dg))) {
          lab_g <- if (!is.null(group_labels) && length(group_labels) >= g) group_labels[g] else g
          num_lab <- suppressWarnings(as.numeric(lab_g))
          if (!is.na(num_lab)) {
            dg[[group_var]] <- num_lab
          } else {
            dg[[group_var]] <- as.character(lab_g)
          }
        }
        
        out_list[[g]] <- dg
      }
      
      out <- do.call(rbind, out_list)
      rownames(out) <- NULL
    } else {
      stop("Recovered data is not a supported object type.", call. = FALSE)
    }
    
    out
  }
  
  make_aux_resid_fit <- function(fit) {
    OPT <- fit@Options
    data_in <- recover_fit_data(fit)
    ordered_vars <- tryCatch(lavaan::lavNames(fit, type = "ov.ord"), error = function(e) character(0))
    
    refit_args <- list(
      model = lavaan::parTable(fit),
      data = data_in,
      group = if (!is.null(OPT$group) && nzchar(OPT$group)) OPT$group else NULL,
      estimator = OPT$estimator,
      missing = OPT$missing,
      std.lv = isTRUE(OPT$std.lv),
      meanstructure = isTRUE(OPT$meanstructure),
      parameterization = OPT$parameterization,
      fixed.x = OPT$fixed.x,
      conditional.x = OPT$conditional.x,
      orthogonal = isTRUE(OPT$orthogonal),
      ordered = if (length(ordered_vars)) ordered_vars else NULL,
      se = "standard",
      test = "standard",
      warn = FALSE
    )
    
    refit_args <- refit_args[!vapply(refit_args, is.null, logical(1L))]
    
    do.call(lavaan::lavaan, refit_args)
  }
  
  fit_for_resid <- fit
  use_aux_fit <- isTRUE(identical(fit@Options$se, "bootstrap"))
  
  if (use_aux_fit) {
    fit_for_resid <- make_aux_resid_fit(fit)
  }
  
  res_try <- tryCatch(
    {
      if (ngroups == 1L) {
        lavaan::lavResiduals(
          fit_for_resid,
          type = type,
          zstat = TRUE,
          summary = TRUE
        )
      } else {
        tmp <- lavaan::lavResiduals(
          fit_for_resid,
          type = type,
          zstat = TRUE,
          summary = TRUE,
          drop.list.single.group = FALSE
        )
        tmp[[group_index]]
      }
    },
    error = function(e) NULL
  )
  
  if (is.null(res_try)) {
    stop(
      "Could not compute standardized residual diagnostics for this fitted object. ",
      "This often occurs when lavaan cannot compute the residual asymptotic covariance matrix.",
      call. = FALSE
    )
  }
  
  res_obj <- res_try
  
  cov_res <- tryCatch(as.matrix(res_obj$cov), error = function(e) NULL)
  cov_z <- tryCatch(as.matrix(res_obj$cov.z), error = function(e) NULL)
  mean_res <- tryCatch(as.numeric(res_obj$mean), error = function(e) NULL)
  mean_z <- tryCatch(as.numeric(res_obj$mean.z), error = function(e) NULL)
  if (!is.null(mean_res) && !is.null(names(res_obj$mean))) {
    names(mean_res) <- names(res_obj$mean)
  }
  if (!is.null(mean_z) && !is.null(names(res_obj$mean.z))) {
    names(mean_z) <- names(res_obj$mean.z)
  }
  sum_obj <- tryCatch(res_obj$summary, error = function(e) NULL)
  
  if (is.null(cov_res) || is.null(cov_z)) {
    stop("Could not extract covariance residual matrices from 'fit'.", call. = FALSE)
  }
  
  var_names <- rownames(cov_z)
  if (is.null(var_names)) {
    var_names <- colnames(cov_z)
  }
  if (is.null(var_names)) {
    var_names <- paste0("V", seq_len(nrow(cov_z)))
  }
  rownames(cov_res) <- var_names
  colnames(cov_res) <- var_names
  rownames(cov_z) <- var_names
  colnames(cov_z) <- var_names
  
  # build long table for covariance residuals
  long_rows <- list()
  idx_row <- 0L
  
  nr <- nrow(cov_z)
  nc <- ncol(cov_z)
  
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      
      if (triangle == "lower" && j > i) next
      if (!include_diagonal && i == j) next
      
      idx_row <- idx_row + 1L
      z_val <- cov_z[i, j]
      raw_val <- cov_res[i, j]
      abs_z <- abs(z_val)
      
      long_rows[[idx_row]] <- data.frame(
        var1 = rownames(cov_z)[i],
        var2 = colnames(cov_z)[j],
        row_index = i,
        col_index = j,
        residual = as.numeric(raw_val),
        z_resid = as.numeric(z_val),
        abs_z = as.numeric(abs_z),
        magnitude = mag_label(abs_z),
        stringsAsFactors = FALSE
      )
    }
  }
  
  cov_long <- do.call(rbind, long_rows)
  rownames(cov_long) <- NULL
  cov_long <- cov_long[order(cov_long$abs_z, decreasing = TRUE), , drop = FALSE]
  
  if (nrow(cov_long) > 0L) {
    top_cov <- head(cov_long, top_n)
  } else {
    top_cov <- data.frame(
      var1 = character(0),
      var2 = character(0),
      row_index = integer(0),
      col_index = integer(0),
      residual = numeric(0),
      z_resid = numeric(0),
      abs_z = numeric(0),
      magnitude = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  # mean residual table
  if (!is.null(mean_res) && !is.null(mean_z) && length(mean_res) == length(mean_z)) {
    mean_df <- data.frame(
      variable = names(mean_z) %||% paste0("V", seq_along(mean_z)),
      residual = as.numeric(mean_res),
      z_resid = as.numeric(mean_z),
      abs_z = abs(as.numeric(mean_z)),
      magnitude = vapply(abs(as.numeric(mean_z)), mag_label, character(1L)),
      stringsAsFactors = FALSE
    )
    mean_df <- mean_df[order(mean_df$abs_z, decreasing = TRUE), , drop = FALSE]
    top_mean <- head(mean_df, top_n)
  } else {
    mean_df <- data.frame(
      variable = character(0),
      residual = numeric(0),
      z_resid = numeric(0),
      abs_z = numeric(0),
      magnitude = character(0),
      stringsAsFactors = FALSE
    )
    top_mean <- mean_df
  }
  
  n_cov_ge_t1 <- sum(cov_long$abs_z >= thresholds[1L], na.rm = TRUE)
  n_cov_ge_t2 <- sum(cov_long$abs_z >= thresholds[2L], na.rm = TRUE)
  
  max_abs_cov <- if (nrow(cov_long)) max(cov_long$abs_z, na.rm = TRUE) else NA_real_
  max_cov_pair <- if (nrow(cov_long)) {
    c(cov_long$var1[1L], cov_long$var2[1L])
  } else {
    c(NA_character_, NA_character_)
  }
  
  counts <- list(
    n_cov_total = nrow(cov_long),
    n_cov_abs_ge_threshold1 = n_cov_ge_t1,
    n_cov_abs_ge_threshold2 = n_cov_ge_t2,
    max_abs_cov_z = max_abs_cov,
    max_cov_pair = max_cov_pair
  )
  
  # build plot data
  plot_df <- NULL
  if (isTRUE(plot)) {
    
    plot_df <- cov_long
    # complete matrix cells for plotting when lower triangle requested
    all_cells <- expand.grid(
      row_index = seq_len(nr),
      col_index = seq_len(nc),
      stringsAsFactors = FALSE
    )
    all_cells$var1 <- rownames(cov_z)[all_cells$row_index]
    all_cells$var2 <- colnames(cov_z)[all_cells$col_index]
    
    plot_df <- merge(
      all_cells,
      plot_df,
      by = c("var1", "var2", "row_index", "col_index"),
      all.x = TRUE,
      sort = FALSE
    )
    
    # blank unwanted cells
    if (triangle == "lower") {
      plot_df$masked <- plot_df$col_index > plot_df$row_index
    } else {
      plot_df$masked <- FALSE
    }
    if (!include_diagonal) {
      plot_df$masked <- plot_df$masked | (plot_df$row_index == plot_df$col_index)
    }
    
    plot_df$fill_value <- if (plot_style == "trafficlight") {
      as.character(plot_df$magnitude)
    } else {
      plot_df$z_resid
    }
    
    plot_df$fill_value[plot_df$masked] <- NA
    plot_df$label_value <- if (plot_style == "trafficlight") {
      vapply(plot_df$abs_z, fmt_num, character(1L), k = digits, apa = FALSE)
    } else {
      vapply(plot_df$z_resid, fmt_num, character(1L), k = digits, apa = FALSE)
    }
    plot_df$label_value[plot_df$masked] <- ""
    
    plot_df$var1 <- factor(plot_df$var1, levels = rev(var_names))
    plot_df$var2 <- factor(plot_df$var2, levels = var_names)
    
    if (plot_style == "trafficlight") {
      plot_df$fill_value <- factor(
        plot_df$fill_value,
        levels = c("minor", "moderate", "notable")
      )
      
      p_out <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = !!rlang::sym("var2"), y = !!rlang::sym("var1"))
      ) +
        ggplot2::geom_tile(
          ggplot2::aes(fill = !!rlang::sym("fill_value")),
          color = "white",
          linewidth = 0.4
        ) +
        ggplot2::scale_fill_manual(
          values = c(
            "minor" = good_color,
            "moderate" = moderate_color,
            "notable" = poor_color
          ),
          na.value = na_color,
          name = "|z| magnitude"
        )
    } else {
      p_out <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = !!rlang::sym("var2"), y = !!rlang::sym("var1"))
      ) +
        ggplot2::geom_tile(
          ggplot2::aes(fill = !!rlang::sym("fill_value")),
          color = "white",
          linewidth = 0.4
        ) +
        ggplot2::scale_fill_gradient2(
          low = neg_color,
          mid = "white",
          high = pos_color,
          midpoint = 0,
          na.value = na_color,
          name = "z residual"
        )
    }
    
    if (isTRUE(show_values)) {
      p_out <- p_out +
        ggplot2::geom_text(
          ggplot2::aes(label = !!rlang::sym("label_value")),
          size = 3
        )
    }
    
    plot_title <- if (plot_style == "trafficlight") {
      "Absolute standardized covariance residuals"
    } else {
      "Signed standardized covariance residuals"
    }
    
    if (ngroups > 1L) {
      plot_subtitle <- paste0("Group ", group_label)
    } else {
      plot_subtitle <- NULL
    }
    
    p_out <- p_out +
      ggplot2::coord_equal() +
      ggplot2::labs(
        title = plot_title,
        subtitle = plot_subtitle,
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.title = ggplot2::element_blank()
      )
  } else {
    p_out <- NULL
  }
  
  summary_df <- NULL
  if (!is.null(sum_obj)) {
    sm <- as.matrix(sum_obj)
    summary_df <- data.frame(
      metric = rownames(sm),
      sm,
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    summary_df <- data.frame(
      metric = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  res <- list(
    summary = summary_df,
    cov_residuals = cov_res,
    cov_z = cov_z,
    mean_residuals = mean_res,
    mean_z = mean_z,
    top_cov_residuals = top_cov,
    top_mean_residuals = top_mean,
    counts = counts,
    plot = p_out,
    group = group_index,
    group_label = group_label,
    type = type,
    thresholds = thresholds,
    digits = digits,
    call = match.call()
  )
  
  if (isTRUE(return_data)) {
    res$plot_data <- plot_df
  }
  
  class(res) <- "lav_localfit"
  return(res)
}


# --- S3 methods: print, summary, plot ---------------------------------------

#' @param x A \code{"lav_localfit"} object.
#' @param ... Additional arguments; unused.
#' @rdname lav_localfit
#' @export
print.lav_localfit <- function(x, ...) {
  digits <- x$digits %||% 3L
  
  fmt <- function(v, k = digits) {
    if (length(v) == 0L || is.na(v)) "NA" else formatC(v, digits = k, format = "f")
  }
  
  fmt_apa <- function(v, k = digits) {
    if (length(v) == 0L || is.na(v)) return("NA")
    out <- formatC(v, digits = k, format = "f")
    sub("^(-?)0\\.", "\\1.", out)
  }
  
  if (!is.null(x$summary) && nrow(x$summary) > 0L) {
    main_metric_name <- if ("metric" %in% names(x$summary)) {
      x$summary$metric[1L]
    } else {
      NA_character_
    }
  } else {
    main_metric_name <- NA_character_
  }
  
  cat("\nLocal fit diagnostics from residual-based evaluation\n")
  cat(rep("-", 74), "\n", sep = "")
  cat("Residual type: ", x$type, "\n", sep = "")
  cat("Group: ", x$group_label, "\n", sep = "")
  
  if (!is.null(x$summary) && nrow(x$summary) > 0L) {
    sm <- x$summary
    total_col <- if ("total" %in% names(sm)) "total" else names(sm)[ncol(sm)]
    if ("metric" %in% names(sm) && nrow(sm) > 0L) {
      cat("Primary residual summary: ",
          toupper(sm$metric[1L]),
          " = ",
          fmt_apa(sm[[total_col]][1L]),
          "\n", sep = "")
    }
  }
  
  cat("Thresholds: |z| < ", fmt(x$thresholds[1L], k = 2),
      " = minor; ", fmt(x$thresholds[1L], k = 2), " to < ",
      fmt(x$thresholds[2L], k = 2), " = moderate; >= ",
      fmt(x$thresholds[2L], k = 2), " = notable\n", sep = "")
  cat("Covariance residuals screened: ", x$counts$n_cov_total, "\n", sep = "")
  cat("Count with |z| >= ", fmt(x$thresholds[1L], k = 2), ": ",
      x$counts$n_cov_abs_ge_threshold1, "\n", sep = "")
  cat("Count with |z| >= ", fmt(x$thresholds[2L], k = 2), ": ",
      x$counts$n_cov_abs_ge_threshold2, "\n", sep = "")
  
  if (!is.na(x$counts$max_abs_cov_z)) {
    cat("Largest absolute standardized covariance residual: ",
        fmt_apa(x$counts$max_abs_cov_z),
        " (", x$counts$max_cov_pair[1L], " with ", x$counts$max_cov_pair[2L], ")\n",
        sep = "")
  }
  
  cat(rep("-", 74), "\n", sep = "")
  
  if (!is.null(x$top_cov_residuals) && nrow(x$top_cov_residuals) > 0L) {
    cat("Largest absolute standardized covariance residuals\n")
    top_cov <- x$top_cov_residuals
    
    var1_txt <- as.character(top_cov$var1)
    var2_txt <- as.character(top_cov$var2)
    z_txt <- vapply(top_cov$z_resid, fmt_apa, character(1L), k = digits)
    abs_txt <- vapply(top_cov$abs_z, fmt_apa, character(1L), k = digits)
    mag_txt <- as.character(top_cov$magnitude)
    
    hdr1 <- "Var1"
    hdr2 <- "Var2"
    hdr3 <- "z"
    hdr4 <- "Magnitude"
    
    w1 <- max(nchar(c(hdr1, var1_txt), type = "width"))
    w2 <- max(nchar(c(hdr2, var2_txt), type = "width"))
    w3 <- max(nchar(c(hdr3, z_txt), type = "width"))
    w4 <- max(nchar(c(hdr4, mag_txt), type = "width"))
    
    row_fmt <- paste0(
      "%-", w1, "s  ",
      "%-", w2, "s  ",
      "%",  w3, "s  ",
      "%-", w4, "s\n"
    )
    
    hdr_line <- sprintf(row_fmt, hdr1, hdr2, hdr3, hdr4)
    
    rule <- paste(rep("-", nchar(sub("\n$", "", hdr_line), type = "width")),
                  collapse = "")
    
    cat(rule, "\n", sep = "")
    cat(hdr_line)
    cat(rule, "\n", sep = "")
    
    for (i in seq_len(nrow(top_cov))) {
      cat(sprintf(
        row_fmt,
        var1_txt[i],
        var2_txt[i],
        z_txt[i],
        mag_txt[i]
      ))
    }
    cat(rule, "\n", sep = "")
  }
  
  # solo se esiste almeno un res delle medie
  if (!is.null(x$top_mean_residuals) &&
      nrow(x$top_mean_residuals) > 0L &&
      any(abs(x$top_mean_residuals$z_resid) > 0, na.rm = TRUE)) { 
    cat("Largest absolute standardized mean residuals\n")
    top_mean <- x$top_mean_residuals
    
    var_txt <- as.character(top_mean$variable)
    z_txt <- vapply(top_mean$z_resid, fmt_apa, character(1L), k = digits)
    abs_txt <- vapply(top_mean$abs_z, fmt_apa, character(1L), k = digits)
    mag_txt <- as.character(top_mean$magnitude)
    
    hdr1 <- "Variable"
    hdr2 <- "z"
    hdr3 <- "Magnitude"
    
    w1 <- max(nchar(c(hdr1, var_txt), type = "width"))
    w2 <- max(nchar(c(hdr2, z_txt), type = "width"))
    w3 <- max(nchar(c(hdr3, mag_txt), type = "width"))
    
    row_fmt <- paste0(
      "%-", w1, "s  ",
      "%",  w2, "s  ",
      "%-", w3, "s\n"
    )
    
    hdr_line <- sprintf(row_fmt, hdr1, hdr2, hdr3)
    
    rule <- paste(rep("-", nchar(sub("\n$", "", hdr_line), type = "width")),
                  collapse = "")
    
    cat(rule, "\n", sep = "")
    cat(hdr_line)
    cat(rule, "\n", sep = "")
    
    for (i in seq_len(nrow(top_mean))) {
      cat(sprintf(
        row_fmt,
        var_txt[i],
        z_txt[i],
        mag_txt[i]
      ))
    }
    cat(rule, "\n", sep = "")
  }
  
  if (!is.null(x$plot)) {
    print(x$plot)
  }
  
  invisible(x)
}

#' @param object A \code{"lav_localfit"} object.
#' @param ... Additional arguments; unused.
#' @rdname lav_localfit
#' @export
summary.lav_localfit <- function(object, ...) {
  print.lav_localfit(object, ...)
  invisible(object)
}

#' @param x A \code{"lav_localfit"} object.
#' @param y Ignored.
#' @param ... Additional arguments; unused.
#' @rdname lav_localfit
#' @export
plot.lav_localfit <- function(x, y = NULL, ...) {
  if (is.null(x$plot)) {
    stop("No plot is stored in this object. Refit with 'plot = TRUE'.", call. = FALSE)
  }
  print(x$plot)
  invisible(x)
}