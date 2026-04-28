#' Simple slopes and interaction plots for fitted 'lavaan' models
#'
#' Computes conditional (simple) slopes of a focal predictor across values
#' of a moderator from a fitted 'lavaan' model that includes their explicit 
#' product term. Plots predicted lines with Wald confidence ribbons and prints
#' an APA-style test of the interaction for easy reporting and interpretation,
#' together with a simple slopes table.
#' 
#' The model should include a main effect for the predictor, a main effect for
#' the moderator, and one explicit product term between them.
#' 
#' The moderator must enter the fitted model as a single numeric
#' variable, which may be continuous, a binary observed moderator coded as a
#' single numeric dummy variable, an observed numeric moderator with a small
#' number of values treated as discrete probe points, or a latent moderator
#' treated as a single continuous latent variable. 
#' 
#' Standard errors use the delta method with the model covariance matrix of the
#' estimates. When moderator values are derived automatically for latent
#' moderators, probe points are based on the estimated latent mean and
#' model-implied latent standard deviation.
#' 
#' @usage
#' lav_slopes(
#'   fit,
#'   outcome,
#'   pred,
#'   modx,
#'   interaction,
#'   data = NULL,
#'   modx.values = NULL,
#'   modx.labels = NULL,
#'   pred.range = NULL,
#'   conf.level = 0.95,
#'   x.label = NULL,
#'   y.label = NULL,
#'   legend.title = NULL,
#'   colors = NULL,
#'   line.size = 0.80,
#'   alpha = 0.20,
#'   table = TRUE,
#'   digits = 3L,
#'   modx_n_unique_cutoff = 4L,
#'   return_data = FALSE
#' )
#'
#' @param fit A fitted 'lavaan' object that includes the product term (required).
#' @param outcome Character. Name of the dependent variable in \code{fit} (required).
#' @param pred Character. Name of the focal predictor whose simple slopes are probed (required).
#' @param modx Character. Name of the moderator. The moderator must appear in
#' the fitted model as a single numeric variable. General nominal moderators
#' with more than two categories are not supported.
#' @param interaction Character. Name of the product term in \code{fit} (e.g., \code{"X_Z"}) (required).
#' @param data \code{data.frame}. Optional raw data (not needed, retained for backward compatibility). 
#' The function automatically recovers observed data from \code{fit} when needed.   
#' @param modx.values Numeric vector of moderator values at which to compute
#' simple slopes. If \code{NULL} and \code{modx} is observed and numeric, the
#' function uses mean minus 1 SD, the mean, and mean plus 1 SD for moderators
#' with more than \code{modx_n_unique_cutoff} unique values, and otherwise
#' uses the observed numeric values as discrete probe points. If \code{modx}
#' is latent in a single-group model, the function uses the estimated latent
#' mean minus 1 latent SD, the latent mean, and the latent mean plus 1 latent
#' SD. 
#' @param modx.labels Character vector. Legend and table labels for \code{modx.values}.
#' By default, the labels are \code{c("-1 SD", "Mean", "+1 SD")} when values are derived as 
#' mean plus or minus 1 SD, and \code{as.character(modx.values)} otherwise.
#' @param pred.range Numeric vector of length 2. Range \code{c(min, max)} for the x-axis
#' for the focal predictor. If \code{NULL} and \code{pred} is observed, the observed 
#' range recovered from \code{fit} is used. If \code{pred} is latent in a single-group model, 
#' the function uses the estimated latent mean minus 2 latent SD and the latent mean plus 2 
#' latent SD. Otherwise, \code{c(-2, 2)} is used.
#' @param conf.level Numeric in (0,1). Confidence level for Wald confidence intervals and ribbons (default: 0.95).
#' @param x.label Character. X-axis label (default: \code{pred}).
#' @param y.label Character. Y-axis label (default: \code{outcome}).
#' @param legend.title Character. Legend title; if \code{NULL}, the legend shows only levels (default: NULL).
#' @param colors Character vector. Colors for lines and ribbons; named vector recommended with names matching \code{modx.labels} (default: Okabe-Ito palette).
#' @param line.size Numeric > 0. Line width (default: 0.80).
#' @param alpha Numeric in (0,1). Ribbon opacity (default 0.20).
#' @param table Logical. Print APA-style interaction test and simple-slopes table (default: \code{TRUE}).
#' @param digits Integer \code{>= 0}. Decimal digits in printed output (default: 3).
#' @param modx_n_unique_cutoff Integer \code{>= 1}. Threshold for treating a numeric moderator
#' as continuous and using mean ± SD (default: 4).
#' @param return_data Logical. If \code{TRUE}, include the plotting data.frame in the returned list (default: FALSE).
#'
#' @return A list of class \code{"lav_slopes"} with elements:
#' \describe{
#'   \item{\code{plot}}{\code{ggplot} object with lines and confidence ribbons.}
#'   \item{\code{slope_table}}{Data frame with moderator levels, simple slopes, SE, z, and CI.}
#'   \item{\code{plot_data}}{Only when \code{return_data = TRUE}: data used to build the plot.}
#' }
#'
#' @section Notes:
#' Estimates are unstandardized; a standardized coefficient for the interaction is
#' also reported for reference. Wald tests assume large-sample normality of the
#' parameter estimates. Multigroup fitted models are not supported. 
#' 
#' @references
#' Aiken, L. S., & West, S. G. (1991). \emph{Multiple regression: Testing and
#' interpreting interactions}. Sage.
#'
#' Preacher, K. J., Curran, P. J., & Bauer, D. J. (2006). Computational tools
#' for probing interactions in multiple linear regression, multilevel modeling,
#' and latent curve analysis. \emph{Journal of Educational and Behavioral
#' Statistics}, \emph{31}(4), 437-448.
#' \doi{10.3102/10769986031004437}
#'
#' Rogosa, D. (1980). Comparing nonparallel regression lines.
#' \emph{Psychological Bulletin}, \emph{88}(2), 307-321.
#' \doi{10.1037/0033-2909.88.2.307}
#'
#' @seealso \code{\link{lav_jn}} for Johnson-Neyman regions of significance.
#' 
#' @examples
#' set.seed(42)
#' X <- rnorm(100); Z <- rnorm(100); X_Z <- X*Z
#' Y <- 0.6*X + 0.6*Z + 0.3*X_Z + rnorm(100, sd = 0.7) 
#' dataset <- data.frame(Y, X, Z, X_Z)
#' fit <- lavaan::sem("Y ~ X + Z + X_Z", data = dataset)
#' lav_slopes(
#' fit = fit, 
#' outcome = "Y", 
#' pred = "X", 
#' modx = "Z", 
#' interaction = "X_Z")
#'
#' @importFrom lavaan parameterEstimates parTable lavInspect lavNames
#' @importFrom stats vcov qnorm sd pnorm
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon scale_colour_manual scale_fill_manual labs
#' @importFrom rlang sym %||%
#' @export
lav_slopes <- function(
    fit,
    outcome,
    pred,
    modx,
    interaction,
    data = NULL,
    modx.values = NULL,
    modx.labels = NULL,
    pred.range = NULL,
    conf.level = 0.95,
    x.label = NULL,
    y.label = NULL,
    legend.title = NULL,
    colors = NULL,
    line.size = 0.80,
    alpha = 0.20,
    table = TRUE,
    digits = 3L,
    modx_n_unique_cutoff = 4L,
    return_data = FALSE
) {
  if (!inherits(fit, "lavaan")) stop("'fit' must be a 'lavaan' object.", call. = FALSE)
  if (missing(interaction) || is.null(interaction))
    stop("'interaction' (product-term name) must be supplied.", call. = FALSE)
  
  ngroups <- tryCatch(lavaan::lavInspect(fit, "ngroups"), error = function(e) 1L)
  if (length(ngroups) != 1L || is.na(ngroups)) ngroups <- 1L
  if (ngroups > 1L) {
    stop("'lav_slopes' does not currently support multigroup fitted models. Please fit one group at a time.", call. = FALSE)
  }
  
# recupera i dati!!!
  fetch_data <- function(f) {
    out <- tryCatch(lavaan::lavInspect(f, "data"), error = function(e) NULL)
    if (is.null(out)) {
      out <- tryCatch(lavaan::lavInspect(f, "data.original"), error = function(e) NULL)
    }
    if (is.null(out)) {
      out <- tryCatch(f@Data@X, error = function(e) NULL)
    }
    if (is.null(out)) {
      return(NULL)
    }
    ov_names <- tryCatch(lavaan::lavNames(f, type = "ov"), error = function(e) NULL)
    if (is.data.frame(out)) {
      out <- as.data.frame(out)
      return(out)
    }
    if (is.matrix(out)) {
      out <- as.data.frame(out)
      if (!is.null(ov_names) && length(ov_names) == ncol(out)) {
        names(out) <- ov_names
      }
      return(out)
    }
    if (is.list(out)) {
      ok_list <- all(vapply(out, function(x) is.data.frame(x) || is.matrix(x), logical(1)))
      if (!ok_list) {
        return(NULL)
      }
      group_labels <- tryCatch(lavaan::lavInspect(f, "group.label"), error = function(e) NULL)
      group_var <- NULL
      group_call <- tryCatch(f@call$group, error = function(e) NULL)
      if (!is.null(group_call)) {
        group_var <- paste(as.character(group_call), collapse = "")
        group_var <- gsub('^"|"$', "", group_var)
      }
      if (is.null(group_var) || length(group_var) != 1L || is.na(group_var) || !nzchar(group_var)) {
        tmp_group <- tryCatch(f@Options$group, error = function(e) NULL)
        if (length(tmp_group) == 1L && !is.na(tmp_group) && nzchar(tmp_group)) {
          group_var <- tmp_group
        } else {
          group_var <- NULL
        }
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
        if (!is.null(group_var) && nzchar(group_var) && !(group_var %in% names(dg))) {
          if (!is.null(group_labels) && length(group_labels) >= g) {
            lab_g <- group_labels[g]
          } else {
            lab_g <- g
          }
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
      return(out)
    }
    NULL
  }
  dat <- fetch_data(fit)
  
  # solo per compatibilita
  if (is.null(dat) && !is.null(data)) {
    dat <- as.data.frame(data)
  }
  
  # estrazione pick a point delle latenti in automatico!!!
  lv_names <- tryCatch(lavaan::lavNames(fit, type = "lv"), error = function(e) character(0))
  # media e ds dilatente dal modello stimato:
  prendi_media_sd_latente <- function(fit, lv) {
    mean_lv <- tryCatch(lavaan::lavInspect(fit, "mean.lv"), error = function(e) NULL)
    cov_lv  <- tryCatch(lavaan::lavInspect(fit, "cov.lv"), error = function(e) NULL)
    if (is.list(mean_lv)) {
      mean_lv <- mean_lv[[1L]]
    }
    if (is.list(cov_lv)) {
      cov_lv <- cov_lv[[1L]]
    }
    if (is.matrix(mean_lv) && nrow(mean_lv) >= 1L) {
      tmp <- as.numeric(mean_lv[1L, ])
      names(tmp) <- colnames(mean_lv)
      mean_lv <- tmp
    }
    mu <- NA_real_
    var_lv <- NA_real_
    if (!is.null(mean_lv) && length(mean_lv) > 0L && lv %in% names(mean_lv)) {
      mu <- as.numeric(mean_lv[lv])
    }
    if (!is.null(cov_lv) && is.matrix(cov_lv)) {
      rn <- rownames(cov_lv)
      cn <- colnames(cov_lv)
      if (!is.null(rn) && !is.null(cn) && lv %in% rn && lv %in% cn) {
        var_lv <- as.numeric(cov_lv[lv, lv])
      }
    }
    if (!is.finite(mu)) {
      mu <- 0
    }
    if (!is.finite(var_lv) || var_lv < 0) {
      return(NULL)
    }
    list(mean = mu, sd = sqrt(var_lv))
  }
# estrai moderatori e label 
  if (is.null(modx.values)) {
    # caso 1: moderatore osservato!!!
    if (!is.null(dat) && modx %in% names(dat)) {
      z <- dat[[modx]]
      if (is.numeric(z) && length(unique(z)) > modx_n_unique_cutoff) {
        m <- mean(z, na.rm = TRUE)
        s <- sd(z, na.rm = TRUE)
        modx.values <- round(c(m - s, m, m + s), 2)
        modx.labels <- modx.labels %||% c("-1 SD", "Mean", "+1 SD")
      } else {
        modx.values <- sort(unique(z))
        modx.labels <- modx.labels %||% as.character(modx.values)
      }
      # caso 2: moderatore latente!!!!
    } else if (modx %in% lv_names) {
      lv_info <- prendi_media_sd_latente(fit, modx)
      if (is.null(lv_info)) {
        stop(
          "Could not derive latent mean and variance for moderator `", modx, "`. ",
          "Supply `modx.values` explicitly.",
          call. = FALSE
        )
      }
      modx.values <- round(
        c(
          lv_info$mean - lv_info$sd,
          lv_info$mean,
          lv_info$mean + lv_info$sd
        ),
        2
      )
      modx.labels <- modx.labels %||% c("-1 SD", "Mean", "+1 SD")
      # caso 3: moderatore non trovato!!!!!!!
    } else {
      stop(
        "`modx` was not found as either an observed or latent variable in `fit`. ",
        "Supply a valid moderator name or specify `modx.values` explicitly.",
        call. = FALSE
      )
    }
  }
  if (!is.null(modx.labels) && length(modx.labels) != length(modx.values))
    stop("'modx.labels' must match length of 'modx.values'.", call. = FALSE)
  modx.labels <- modx.labels %||% as.character(modx.values)
  
  if (!is.numeric(modx.values)) {
    stop(
      "'modx.values' must be numeric. General nominal moderators with more than two categories are not supported; binary observed moderators must be represented as a single numeric dummy-coded variable.",
      call. = FALSE
    )
  }

# range del predittore focale:
  if (is.null(pred.range)) {
    if (!is.null(dat) && pred %in% names(dat)) {
      pred.range <- range(dat[[pred]], na.rm = TRUE)
    } else if (pred %in% lv_names) {
      lv_info_pred <- prendi_media_sd_latente(fit, pred)
      if (!is.null(lv_info_pred)) {
        pred.range <- c(
          lv_info_pred$mean - 2 * lv_info_pred$sd,
          lv_info_pred$mean + 2 * lv_info_pred$sd
        )
      } else {
        pred.range <- c(-2, 2)
      }
    } else {
      pred.range <- c(-2, 2)
    }
  }
  x_seq <- seq(pred.range[1L], pred.range[2L], length.out = 1000L)
  
# stime dei parametri e vcov:
  pe <- lavaan::parameterEstimates(fit, standardized = FALSE)
  pt <- lavaan::parTable(fit)
  vc <- tryCatch(stats::vcov(fit),
                 error = function(e) lavaan::lavInspect(fit, "vcov"))
  if (!is.matrix(vc)) stop("Could not retrieve covariance matrix from 'fit'.", call. = FALSE)
  
  find_par <- function(lhs, rhs, op = "~", rhs_regex = NULL,
                       required = TRUE, default_est = 0) {
    idx_pe <- which(pe$lhs == lhs & pe$rhs == rhs & pe$op == op)
    if (length(idx_pe) == 0L && !is.null(rhs_regex))
      idx_pe <- which(pe$lhs == lhs & pe$op == op & grepl(rhs_regex, pe$rhs))
    if (length(idx_pe) == 0L) {
      if (!required) { est_val <- default_est; idx_pe <- NA_integer_ } else {
        stop(sprintf("Parameter '%s %s %s' not found.", lhs, op, rhs), call. = FALSE)
      }
    } else {
      idx_pe  <- idx_pe[1L]; est_val <- pe$est[idx_pe]
    }
    idx_pt <- which(pt$lhs == lhs & pt$rhs == rhs & pt$op == op)
    if (length(idx_pt) == 0L && !is.null(rhs_regex))
      idx_pt <- which(pt$lhs == lhs & pt$op == op & grepl(rhs_regex, pt$rhs))
    free_idx <- if (length(idx_pt) == 0L) 0L else {
      fi <- pt$free[idx_pt[1L]]; if (is.na(fi)) 0L else as.integer(fi)
    }
    list(est = est_val, free = free_idx, idx_pe = idx_pe)
  }
  
# estrai inetercetta + effetti principali + interazione
  p0 <- find_par(outcome, "",   op = "~1", required = FALSE, default_est = 0)
  p1 <- find_par(outcome, pred, op = "~")
  p2 <- find_par(outcome, modx, op = "~")
  p3 <- find_par(outcome, interaction, op = "~",
                 rhs_regex = paste0("(", pred, ").*(", modx, ")|(",
                                    modx, ").*(", pred, ")"))
  
  safe_free <- function(x) { if (is.na(x) || length(x) == 0L) 0L else as.integer(x) }
  p0$free <- safe_free(p0$free); p1$free <- safe_free(p1$free)
  p2$free <- safe_free(p2$free); p3$free <- safe_free(p3$free)
  
  vc_ <- function(a, b) {
    a <- safe_free(a); b <- safe_free(b)
    if (a <= 0L || b <= 0L) return(0)
    if (a > nrow(vc) || b > ncol(vc)) return(0)
    val <- vc[a, b]; if (is.na(val)) 0 else val
  }
  
  var_b1   <- vc_(p1$free, p1$free)
  var_b3   <- vc_(p3$free, p3$free)
  cov_b1b3 <- vc_(p1$free, p3$free)
  
  var_b0   <- vc_(p0$free, p0$free)
  var_b2   <- vc_(p2$free, p2$free)
  cov_b0b2 <- vc_(p0$free, p2$free)
  cov_b0b1 <- vc_(p0$free, p1$free)
  cov_b0b3 <- vc_(p0$free, p3$free)
  cov_b1b2 <- vc_(p1$free, p2$free)
  cov_b2b3 <- vc_(p2$free, p3$free)
  
# slope e anche predizioni
  z_crit <- stats::qnorm(1 - (1 - conf.level) / 2)
  slope_tbl <- data.frame(
    Moderator = modx.labels,
    Mod_value = modx.values,
    Slope = NA_real_, SE = NA_real_, z = NA_real_,
    CI_low = NA_real_, CI_high = NA_real_,
    stringsAsFactors = FALSE
  )
  plot_df <- NULL
  
  for (i in seq_along(modx.values)) {
    z0 <- modx.values[i]
    
    slope     <- p1$est + p3$est * z0
    var_slope <- var_b1 + 2 * z0 * cov_b1b3 + (z0 ^ 2) * var_b3
    se_slope  <- sqrt(max(var_slope, 0))
    
    slope_tbl$Slope[i] <- slope
    slope_tbl$SE[i] <- se_slope
    slope_tbl$z[i] <- if (se_slope > 0) slope / se_slope else NA_real_
    slope_tbl$CI_low[i] <- slope - z_crit * se_slope
    slope_tbl$CI_high[i] <- slope + z_crit * se_slope
    
    intercept <- p0$est + p2$est * z0
    var_int <- max(var_b0 + 2 * z0 * cov_b0b2 + (z0 ^ 2) * var_b2, 0)
    
    y_hat <- intercept + slope * x_seq
    se_yhat <- sqrt(pmax((x_seq ^ 2) * var_slope + var_int +
                           2 * x_seq * (cov_b0b1 + z0 * cov_b0b3 +
                                          z0 * cov_b1b2 + (z0 ^ 2) * cov_b2b3), 0))
    
    plot_df <- rbind(plot_df,
                     data.frame(
                       X = x_seq, Y_hat = y_hat,
                       CI_low = y_hat - z_crit * se_yhat,
                       CI_high= y_hat + z_crit * se_yhat,
                       LegendGroup = factor(modx.labels[i], levels = modx.labels),
                       stringsAsFactors = FALSE
                     ))
  }
  
# colori okabe ito
  if (is.null(colors)) {
    okabe_ito <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    colors <- okabe_ito[seq_along(modx.labels)]
  }
  if (is.null(names(colors)) || !all(modx.labels %in% names(colors)))
    names(colors) <- modx.labels
  
# plotting con ggplot 
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = !!rlang::sym("X"),
      y = !!rlang::sym("Y_hat"),
      colour = !!rlang::sym("LegendGroup"),
      fill   = !!rlang::sym("LegendGroup")
    )
  ) +
    ggplot2::geom_line(linewidth = line.size) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = !!rlang::sym("CI_low"),
        ymax = !!rlang::sym("CI_high")
      ),
      alpha = alpha, colour = NA
    ) +
    ggplot2::scale_colour_manual(values = colors, name = legend.title) +
    ggplot2::scale_fill_manual(values = colors, name = legend.title) +
    ggplot2::labs(x = x.label %||% pred, y = y.label %||% outcome)
  
  
  # output in console 
  
  # risultati per gli S3
  interaction_test <- list(
    b = p3$est,
    se = if (!is.na(p3$idx_pe)) pe$se[p3$idx_pe] else NA_real_,
    z = {
      se <- if (!is.na(p3$idx_pe)) pe$se[p3$idx_pe] else NA_real_
      if (!is.na(se) && se > 0) p3$est / se else NA_real_
    },
    p = {
      se <- if (!is.na(p3$idx_pe)) pe$se[p3$idx_pe] else NA_real_
      if (!is.na(se) && se > 0) 2 * stats::pnorm(abs(p3$est / se), lower.tail = FALSE) else NA_real_
    },
    ci = {
      se <- if (!is.na(p3$idx_pe)) pe$se[p3$idx_pe] else NA_real_
      c(p3$est - z_crit * se, p3$est + z_crit * se)
    },
    beta_std = {
      pe_std <- lavaan::parameterEstimates(fit, standardized = TRUE)
      idx_std <- which(pe_std$lhs == outcome & pe_std$rhs == interaction & pe_std$op == "~")
      if (length(idx_std)) pe_std$std.all[idx_std[1L]] else NA_real_
    }
  )
  
  res <- list(
    plot = p,
    slope_table  = slope_tbl,
    labels = list(outcome = outcome, pred = pred, modx = modx, interaction = interaction),
    conf.level = conf.level,
    digits = digits,
    interaction = interaction_test,
    print_table = table,
    call = match.call()
  )
  if (isTRUE(return_data)) res$plot_data <- plot_df
  
  class(res) <- "lav_slopes"
  return(res)
}

### metodi S3 print e summary:

#' @param x A 'lav_slopes' object.
#' @param ... Additional arguments; unused.
#' @rdname lav_slopes
#' @export
print.lav_slopes <- function(x, ...) {
  lab <- x$labels
  it <- x$interaction
  digits <- x$digits
  cl <- x$conf.level
  
  fmt <- function(v, k = digits) {
    if (is.na(v)) "NA" else formatC(v, format = "f", digits = k)
  }
  
  fmtp <- function(p, k = digits) {
    if (is.na(p)) return("NA")
    cutoff <- 10^(-k)
    cutoff_txt <- sub("^0\\.", ".", formatC(cutoff, digits = k, format = "f"))
    if (p < cutoff) return(paste0("< ", cutoff_txt))
    out <- formatC(p, digits = k, format = "f")
    sub("^0\\.", ".", out)
  }
  
  zfmt <- function(z) {
    if (is.na(z)) "NA" else formatC(z, format = "f", digits = digits)
  }
  
  cat(
    "\nInteraction effect (", lab$pred, " * ", lab$modx, " -> ", lab$outcome, "): ",
    "b = ", fmt(it$b), ", SE = ", fmt(it$se),
    ", beta = ", fmt(it$beta_std), ", z = ", zfmt(it$z),
    ", p = ", fmtp(it$p),
    ", ", sprintf("%.0f%% CI", cl * 100), " [", fmt(it$ci[1L]), ", ", fmt(it$ci[2L]), "]\n",
    sep = ""
  )
  
  if (isTRUE(x$print_table)) {
    tbl <- x$slope_table
    
    cat(
      "\nSimple slopes of", lab$pred, "predicting", lab$outcome,
      "at levels of", lab$modx, sprintf("(%.0f%% CI)", cl * 100), "\n",
      sep = " "
    )
    
    mod_txt <- as.character(tbl$Moderator)
    modval_txt <- as.character(tbl$Mod_value)
    slope_txt <- vapply(tbl$Slope, fmt, character(1L))
    se_txt <- vapply(tbl$SE, fmt, character(1L))
    z_txt <- vapply(tbl$z, zfmt, character(1L))
    ci_low_txt <- vapply(tbl$CI_low, fmt, character(1L))
    ci_high_txt <- vapply(tbl$CI_high, fmt, character(1L))
    ci_txt <- paste0("[", ci_low_txt, ", ", ci_high_txt, "]")
    
    hdr_mod <- "Moderator"
    hdr_modval <- "Mod. value"
    hdr_slope <- "Slope"
    hdr_se <- "SE"
    hdr_z <- "z"
    hdr_ci <- sprintf("%.0f%% CI", cl * 100)
    
    w_mod <- max(nchar(c(hdr_mod, mod_txt), type = "width"))
    w_modval <- max(nchar(c(hdr_modval, modval_txt), type = "width"))
    w_slope <- max(nchar(c(hdr_slope, slope_txt), type = "width"))
    w_se <- max(nchar(c(hdr_se, se_txt), type = "width"))
    w_z <- max(nchar(c(hdr_z, z_txt), type = "width"))
    w_ci <- max(nchar(c(hdr_ci, ci_txt), type = "width"))
    
    row_fmt <- paste0(
      "%-", w_mod, "s  ",
      "%-", w_modval, "s  ",
      "%",  w_slope, "s  ",
      "%",  w_se, "s  ",
      "%",  w_z, "s  ",
      "%-", w_ci, "s\n"
    )
    
    header_line <- sprintf(
      row_fmt,
      hdr_mod, hdr_modval, hdr_slope, hdr_se, hdr_z, hdr_ci
    )
    
    rule_width <- nchar(sub("\n$", "", header_line), type = "width") + 3L
    rule <- paste(rep("-", rule_width), collapse = "")
    
    cat(rule, "\n", sep = "")
    cat(header_line)
    cat(rule, "\n", sep = "")
    
    for (i in seq_len(nrow(tbl))) {
      cat(sprintf(
        row_fmt,
        mod_txt[i],
        modval_txt[i],
        slope_txt[i],
        se_txt[i],
        z_txt[i],
        ci_txt[i]
      ))
    }
    
    cat(rule, "\n", sep = "")
  }
  
  cat("\n")
  print(x$plot)
  invisible(x)
}

#' @param object A 'lav_slopes' object.
#' @param ... Additional arguments; unused.
#' @rdname lav_slopes
#' @export
summary.lav_slopes <- function(object, ...) {
  print.lav_slopes(object, ...)
  invisible(object)
}
