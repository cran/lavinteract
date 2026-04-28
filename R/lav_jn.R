#' Johnson-Neyman regions of significance for fitted 'lavaan' models
#'
#' Computes the Johnson-Neyman (JN) interval for a continuous moderator from a
#' fitted 'lavaan' model that includes an explicit product term. 
#' Plots the conditional slope with a confidence band, shaded significance
#' regions, and an observed-moderator density strip.
#'
#' The model should include a main effect for the predictor, a main effect for
#' the moderator, and their product term. Standard errors are obtained via the
#' delta method using the model-implied covariance matrix of the parameter
#' estimates. If the model was fitted with a robust estimator (e.g., MLR),
#' the robust covariance matrix returned by \code{vcov(fit)} is used
#' automatically, so no separate correction is needed.
#'
#' The JN boundary values are the real roots of the quadratic equation that
#' equates the squared conditional slope to the squared critical value times
#' its variance. When the discriminant is negative, the conditional slope is
#' either significant or non-significant across the entire moderator range and
#' no finite boundary exists.
#'
#' This function is restricted to continuous moderators. For categorical
#' moderators, use \code{\link{lav_slopes}}.
#'
#' @usage
#' lav_jn(
#'   fit,
#'   outcome,
#'   pred,
#'   modx,
#'   interaction,
#'   data = NULL,
#'   conf.level = 0.95,
#'   modx.range = NULL,
#'   n_grid = 1000L,
#'   x.label = NULL,
#'   y.label = NULL,
#'   sig.color = "#009E73",
#'   nonsig.color = "#D55E00",
#'   line.color = "#0072B2",
#'   alpha = 0.20,
#'   line.size = 0.80,
#'   rug = TRUE,
#'   plot = TRUE,
#'   digits = 3L,
#'   return_data = FALSE
#' )
#'
#' @param fit A fitted 'lavaan' object that includes the product term (required).
#' @param outcome Character. Name of the dependent variable in \code{fit} (required).
#' @param pred Character. Name of the focal predictor whose conditional slope
#'   is probed (required).
#' @param modx Character. Name of the continuous moderator (required). If
#'   \code{modx} is categorical (fewer than 4 unique values), the function
#'   stops with a message directing the user to \code{\link{lav_slopes}}.
#' @param interaction Character. Name of the product term in \code{fit}
#'   (e.g., \code{"X_Z"}) (required).
#' @param data \code{data.frame}. Optional raw data (not needed; retained for
#'   backward compatibility). The function automatically recovers observed data
#'   from \code{fit} when needed.
#' @param conf.level Numeric in (0, 1). Confidence level for the Wald band and
#'   the JN critical value (default: 0.95).
#' @param modx.range Numeric vector of length 2. Range \code{c(min, max)} for
#'   the moderator axis. If \code{NULL}, the observed range of \code{modx} in
#'   the data is used, extended by 5\% on each side for visual padding
#'   (default: NULL).
#' @param n_grid Integer. Number of moderator values at which the conditional
#'   slope is evaluated for plotting (default: 1000).
#' @param x.label Character. X-axis label (default: \code{modx}).
#' @param y.label Character. Y-axis label (default: \code{"Conditional slope of <pred>"}).
#' @param sig.color Character. Fill color for significant-slope regions
#'   (default: \code{"#009E73"}, Okabe-Ito bluish green).
#' @param nonsig.color Character. Fill color for non-significant regions
#'   (default: \code{"#D55E00"}, Okabe-Ito vermillion).
#' @param line.color Character. Color for the conditional-slope line
#'   (default: \code{"#0072B2"}, Okabe-Ito blue).
#' @param alpha Numeric in (0, 1). Opacity of the confidence ribbon and the
#'   region shading (default: 0.20).
#' @param line.size Numeric > 0. Width of the conditional-slope line
#'   (default: 0.80).
#' @param rug Logical. If \code{TRUE}, draw a rug (density strip) of observed
#'   moderator values along the x-axis (default: TRUE).
#' @param plot Logical. If \code{TRUE}, produce the JN plot (default: TRUE).
#' @param digits Integer \code{>= 0}. Decimal digits in printed output
#'   (default: 3).
#' @param return_data Logical. If \code{TRUE}, include the plotting
#'   \code{data.frame} in the returned list (default: FALSE).
#'
#' @return A list of class \code{"lav_jn"} with elements:
#' \describe{
#'   \item{\code{jn_points}}{Numeric vector of length 0, 1, or 2 giving the
#'     Johnson-Neyman boundary values of the moderator (the roots of the
#'     quadratic). \code{numeric(0)} when no real root exists.}
#'   \item{\code{signif_regions}}{A data frame with columns \code{lower},
#'     \code{upper}, and \code{significance} (\code{"significant"} or
#'     \code{"non-significant"}), describing the regions over the observed
#'     moderator range.}
#'   \item{\code{plot}}{A \code{ggplot} object (or \code{NULL} when
#'     \code{plot = FALSE}).}
#'   \item{\code{observed_support}}{Numeric length-2: the minimum and maximum of
#'     the observed moderator values.}
#'   \item{\code{interaction_test}}{List with the unstandardized interaction
#'     coefficient (\code{b}), its \code{se}, \code{z}, \code{p}, Wald
#'     confidence interval (\code{ci}), and standardized beta (\code{beta_std}).}
#'   \item{\code{labels}}{List of the user-supplied variable names.}
#'   \item{\code{conf.level}}{Confidence level used.}
#'   \item{\code{digits}}{Default digits for printing.}
#'   \item{\code{call}}{Matched call.}
#'   \item{\code{plot_data}}{Only when \code{return_data = TRUE}: data frame
#'     used to build the plot.}
#' }
#'
#' @section Notes:
#' All estimates are unstandardized; a standardized coefficient for the interaction is
#' also reported for reference. Wald tests assume large-sample normality of the
#' parameter estimates. When the discriminant of the JN quadratic is negative,
#' the conditional slope is either significant or non-significant everywhere;
#' the function reports which case applies rather than returning a boundary.
#'
#' @references
#' Johnson, P. O., & Neyman, J. (1936). Tests of certain linear hypotheses and
#' their application to some educational problems. \emph{Statistical Research
#' Memoirs}, \emph{1}, 57-93.
#'
#' Bauer, D. J., & Curran, P. J. (2005). Probing interactions in fixed and
#' multilevel regression: Inferential and graphical techniques.
#' \emph{Multivariate Behavioral Research}, \emph{40}(3), 373-400.
#' \doi{10.1207/s15327906mbr4003_5}
#'
#' Carden, S. W., Holtzman, N. S., & Strube, M. J. (2017). CAHOST: An Excel
#' workbook for facilitating the Johnson-Neyman technique for two-way
#' interactions in multiple regression. \emph{Frontiers in Psychology},
#' \emph{8}, Article 1293. \doi{10.3389/fpsyg.2017.01293}
#'
#' @seealso \code{\link{lav_slopes}} for pick-a-point simple slopes.
#'
#' @examples
#' set.seed(42)
#' X <- rnorm(100); Z <- rnorm(100); X_Z <- X*Z
#' Y <- 0.6*X + 0.6*Z + 0.3*X_Z + rnorm(100, sd = 0.7) 
#' dataset <- data.frame(Y, X, Z, X_Z)
#' fit <- lavaan::sem("Y ~ X + Z + X_Z", data = dataset)
#' lav_jn(
#' fit = fit,
#' outcome = "Y",
#' pred = "X",
#' modx = "Z",
#' interaction = "X_Z"
#' )
#'
#' @importFrom lavaan parameterEstimates parTable lavInspect lavNames
#' @importFrom stats vcov qnorm pnorm 
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_hline geom_vline
#'   geom_rect geom_rug labs scale_fill_manual scale_x_continuous
#' @importFrom rlang sym %||%
#' @export
lav_jn <- function(
    fit,
    outcome,
    pred,
    modx,
    interaction,
    data = NULL,
    conf.level = 0.95,
    modx.range = NULL,
    n_grid = 1000L,
    x.label = NULL,
    y.label = NULL,
    sig.color = "#009E73",
    nonsig.color = "#D55E00",
    line.color = "#0072B2",
    alpha = 0.20,
    line.size = 0.80,
    rug = TRUE,
    plot = TRUE,
    digits = 3L,
    return_data = FALSE
) {

  if (!inherits(fit, "lavaan"))
    stop("'fit' must be a 'lavaan' object.", call. = FALSE)
  if (missing(interaction) || is.null(interaction))
    stop("'interaction' (product-term name) must be supplied.", call. = FALSE)

  fetch_data <- function(f) {
    out <- tryCatch(lavaan::lavInspect(f, "data"), error = function(e) NULL)
    if (is.null(out)) {
      out <- tryCatch(lavaan::lavInspect(f, "data.original"),
                       error = function(e) NULL)
    }
    if (is.null(out)) {
      out <- tryCatch(f@Data@X, error = function(e) NULL)
    }
    if (is.null(out)) return(NULL)

    ov_names <- tryCatch(lavaan::lavNames(f, type = "ov"),
                          error = function(e) NULL)
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
      ok_list <- all(vapply(out,
                            function(x) is.data.frame(x) || is.matrix(x),
                            logical(1)))
      if (!ok_list) return(NULL)

      group_labels_fd <- tryCatch(lavaan::lavInspect(f, "group.label"),
                                   error = function(e) NULL)
      group_var_fd <- NULL
      group_call <- tryCatch(f@call$group, error = function(e) NULL)
      if (!is.null(group_call)) {
        group_var_fd <- paste(as.character(group_call), collapse = "")
        group_var_fd <- gsub('^"|"$', "", group_var_fd)
      }
      if (is.null(group_var_fd) || length(group_var_fd) != 1L ||
          is.na(group_var_fd) || !nzchar(group_var_fd)) {
        tmp_group <- tryCatch(f@Options$group, error = function(e) NULL)
        if (length(tmp_group) == 1L && !is.na(tmp_group) && nzchar(tmp_group)) {
          group_var_fd <- tmp_group
        } else {
          group_var_fd <- NULL
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
        if (!is.null(group_var_fd) && nzchar(group_var_fd) &&
            !(group_var_fd %in% names(dg))) {
          if (!is.null(group_labels_fd) && length(group_labels_fd) >= g) {
            lab_g <- group_labels_fd[g]
          } else {
            lab_g <- g
          }
          num_lab <- suppressWarnings(as.numeric(lab_g))
          if (!is.na(num_lab)) {
            dg[[group_var_fd]] <- num_lab
          } else {
            dg[[group_var_fd]] <- as.character(lab_g)
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
  if (is.null(dat) && !is.null(data)) {
    dat <- as.data.frame(data)
  }

  if (!is.null(dat) && modx %in% names(dat)) {
    modx_obs <- dat[[modx]]
    if (!is.numeric(modx_obs)) {
      stop(
        "The moderator '", modx, "' is not numeric. ",
        "lav_jn() is designed for continuous moderators only. ",
        "For categorical moderators, use lav_slopes().",
        call. = FALSE
      )
    }
    n_unique <- length(unique(modx_obs[!is.na(modx_obs)]))
    if (n_unique < 4L) {
      stop(
        "The moderator '", modx, "' has only ", n_unique, " unique value(s). ",
        "The Johnson-Neyman technique requires a continuous moderator. ",
        "For categorical or quasi-categorical moderators, use lav_slopes().",
        call. = FALSE
      )
    }
  } else {

        lv_names <- tryCatch(lavaan::lavNames(fit, type = "lv"),
                          error = function(e) character(0))
    if (!(modx %in% lv_names)) {
      stop(
        "'modx' ('", modx, "') was not found as an observed or latent variable ",
        "in 'fit'. Supply a valid moderator name.",
        call. = FALSE
      )
    }
    modx_obs <- NULL
  }

  if (!is.null(modx_obs)) {
    obs_min <- min(modx_obs, na.rm = TRUE)
    obs_max <- max(modx_obs, na.rm = TRUE)
    observed_support <- c(obs_min, obs_max)
  } else {

    ngroups <- tryCatch(lavaan::lavInspect(fit, "ngroups"),
                         error = function(e) 1L)
    mean_lv <- tryCatch(lavaan::lavInspect(fit, "mean.lv"),
                          error = function(e) NULL)
    cov_lv  <- tryCatch(lavaan::lavInspect(fit, "cov.lv"),
                          error = function(e) NULL)
    if (is.list(mean_lv)) mean_lv <- mean_lv[[1L]]
    if (is.list(cov_lv))  cov_lv  <- cov_lv[[1L]]
    if (is.matrix(mean_lv) && nrow(mean_lv) >= 1L) {
      tmp <- as.numeric(mean_lv[1L, ])
      names(tmp) <- colnames(mean_lv)
      mean_lv <- tmp
    }
    mu_z <- if (!is.null(mean_lv) && modx %in% names(mean_lv)) {
      as.numeric(mean_lv[modx])
    } else { 0 }
    sd_z <- if (!is.null(cov_lv) && is.matrix(cov_lv) &&
                modx %in% rownames(cov_lv)) {
      sqrt(as.numeric(cov_lv[modx, modx]))
    } else { 1 }
    observed_support <- c(mu_z - 3 * sd_z, mu_z + 3 * sd_z)
    obs_min <- observed_support[1L]
    obs_max <- observed_support[2L]
  }

  if (is.null(modx.range)) {
    pad <- 0.05 * (obs_max - obs_min)
    if (pad < 1e-8) pad <- 0.5
    modx.range <- c(obs_min - pad, obs_max + pad)
  }

  pe <- lavaan::parameterEstimates(fit, standardized = FALSE)
  pt <- lavaan::parTable(fit)
  vc <- tryCatch(stats::vcov(fit),
                 error = function(e) lavaan::lavInspect(fit, "vcov"))
  if (!is.matrix(vc))
    stop("Could not retrieve covariance matrix from 'fit'.", call. = FALSE)

  find_par <- function(lhs, rhs, op = "~", rhs_regex = NULL,
                       required = TRUE, default_est = 0) {
    idx_pe <- which(pe$lhs == lhs & pe$rhs == rhs & pe$op == op)
    if (length(idx_pe) == 0L && !is.null(rhs_regex))
      idx_pe <- which(pe$lhs == lhs & pe$op == op &
                        grepl(rhs_regex, pe$rhs))
    if (length(idx_pe) == 0L) {
      if (!required) {
        est_val <- default_est
        idx_pe  <- NA_integer_
      } else {
        stop(sprintf("Parameter '%s %s %s' not found.", lhs, op, rhs),
             call. = FALSE)
      }
    } else {
      idx_pe <- idx_pe[1L]
      est_val <- pe$est[idx_pe]
    }
    idx_pt <- which(pt$lhs == lhs & pt$rhs == rhs & pt$op == op)
    if (length(idx_pt) == 0L && !is.null(rhs_regex))
      idx_pt <- which(pt$lhs == lhs & pt$op == op &
                        grepl(rhs_regex, pt$rhs))
    free_idx <- if (length(idx_pt) == 0L) 0L else {
      fi <- pt$free[idx_pt[1L]]
      if (is.na(fi)) 0L else as.integer(fi)
    }
    list(est = est_val, free = free_idx, idx_pe = idx_pe)
  }

  p1 <- find_par(outcome, pred, op = "~")
  p2 <- find_par(outcome, modx, op = "~")  
  p3 <- find_par(outcome, interaction, op = "~",
                 rhs_regex = paste0("(", pred, ").*(", modx, ")|(",
                                    modx, ").*(", pred, ")"))

  safe_free <- function(x) {
    if (is.na(x) || length(x) == 0L) 0L else as.integer(x)
  }
  p1$free <- safe_free(p1$free)
  p3$free <- safe_free(p3$free)

  vc_ <- function(a, b) {
    a <- safe_free(a); b <- safe_free(b)
    if (a <= 0L || b <= 0L) return(0)
    if (a > nrow(vc) || b > ncol(vc)) return(0)
    val <- vc[a, b]
    if (is.na(val)) 0 else val
  }

  b1 <- p1$est
  b3 <- p3$est

  var_b1 <- vc_(p1$free, p1$free)
  var_b3 <- vc_(p3$free, p3$free)
  cov_b1b3 <- vc_(p1$free, p3$free)

  # quadratica di johnson neyman!!!
  z_crit <- stats::qnorm(1 - (1 - conf.level) / 2)

  A <- b3^2 - z_crit^2 * var_b3
  B <- 2 * b1 * b3 - 2 * z_crit^2 * cov_b1b3
  C <- b1^2 - z_crit^2 * var_b1

  discriminant <- B^2 - 4 * A * C
  jn_points <- numeric(0)

  if (abs(A) < 1e-14) {
    # degenera
    if (abs(B) > 1e-14) {
      jn_points <- -C / B
    }
  } else if (discriminant >= 0) {
    sqrt_disc <- sqrt(discriminant)
    root1 <- (-B - sqrt_disc) / (2 * A)
    root2 <- (-B + sqrt_disc) / (2 * A)
    jn_points <- sort(c(root1, root2))
  }

  slope_sig_at <- function(z_val) {
    slope <- b1 + b3 * z_val
    var_slope <- var_b1 + 2 * z_val * cov_b1b3 + z_val^2 * var_b3
    se_slope  <- sqrt(max(var_slope, 0))
    if (se_slope < 1e-15) return(TRUE)
    abs(slope / se_slope) > z_crit
  }

  jn_in_range <- jn_points[jn_points >= obs_min & jn_points <= obs_max]
  breaks <- sort(unique(c(obs_min, jn_in_range, obs_max)))

  signif_regions <- data.frame(
    lower = numeric(0),
    upper = numeric(0),
    significance = character(0),
    stringsAsFactors = FALSE
  )
  for (k in seq_len(length(breaks) - 1L)) {
    midpoint <- (breaks[k] + breaks[k + 1L]) / 2
    is_sig <- slope_sig_at(midpoint)
    signif_regions <- rbind(signif_regions, data.frame(
      lower = breaks[k],
      upper = breaks[k + 1L],
      significance = if (is_sig) "significant" else "non-significant",
      stringsAsFactors = FALSE
    ))
  }

  interaction_test <- list(
    b = p3$est,
    se = if (!is.na(p3$idx_pe)) pe$se[p3$idx_pe] else NA_real_,
    z = {
      se_int <- if (!is.na(p3$idx_pe)) pe$se[p3$idx_pe] else NA_real_
      if (!is.na(se_int) && se_int > 0) p3$est / se_int else NA_real_
    },
    p = {
      se_int <- if (!is.na(p3$idx_pe)) pe$se[p3$idx_pe] else NA_real_
      if (!is.na(se_int) && se_int > 0) {
        2 * stats::pnorm(abs(p3$est / se_int), lower.tail = FALSE)
      } else { NA_real_ }
    },
    ci = {
      se_int <- if (!is.na(p3$idx_pe)) pe$se[p3$idx_pe] else NA_real_
      c(p3$est - z_crit * se_int, p3$est + z_crit * se_int)
    },
    beta_std = {
      pe_std <- lavaan::parameterEstimates(fit, standardized = TRUE)
      idx_std <- which(pe_std$lhs == outcome &
                         pe_std$rhs == interaction &
                         pe_std$op == "~")
      if (length(idx_std)) pe_std$std.all[idx_std[1L]] else NA_real_
    }
  )

  # griglia di plot
  z_seq <- seq(modx.range[1L], modx.range[2L], length.out = n_grid)

  plot_df <- data.frame(
    modx_value = z_seq,
    slope = NA_real_,
    se = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(n_grid)) {
    z0 <- z_seq[i]
    slope <- b1 + b3 * z0
    var_slope <- var_b1 + 2 * z0 * cov_b1b3 + z0^2 * var_b3
    se_slope <- sqrt(max(var_slope, 0))
    plot_df$slope[i] <- slope
    plot_df$se[i] <- se_slope
    plot_df$ci_lower[i] <- slope - z_crit * se_slope
    plot_df$ci_upper[i] <- slope + z_crit * se_slope
  }

  p_out <- NULL
  if (isTRUE(plot)) {

    y_lo <- min(plot_df$ci_lower, na.rm = TRUE)
    y_hi <- max(plot_df$ci_upper, na.rm = TRUE)
    y_pad <- 0.05 * (y_hi - y_lo)
    if (y_pad < 1e-8) y_pad <- 0.1
    rect_ymin <- y_lo - y_pad
    rect_ymax <- y_hi + y_pad

    all_breaks <- sort(unique(c(modx.range[1L], jn_in_range, modx.range[2L])))
    region_df <- data.frame(
      xmin = numeric(0), xmax = numeric(0),
      fill = character(0), stringsAsFactors = FALSE
    )
    for (k in seq_len(length(all_breaks) - 1L)) {
      mid  <- (all_breaks[k] + all_breaks[k + 1L]) / 2
      is_s <- slope_sig_at(mid)
      region_df <- rbind(region_df, data.frame(
        xmin = all_breaks[k],
        xmax = all_breaks[k + 1L],
        fill = if (is_s) "Significant" else "Non-significant",
        stringsAsFactors = FALSE
      ))
    }
    region_df$fill <- factor(region_df$fill,
                             levels = c("Significant", "Non-significant"))

    region_colors <- c("Significant"     = sig.color,
                       "Non-significant"  = nonsig.color)

    p_out <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = !!rlang::sym("modx_value"))
    ) +
      ggplot2::geom_rect(
        data = region_df,
        ggplot2::aes(
          xmin = !!rlang::sym("xmin"),
          xmax = !!rlang::sym("xmax"),
          fill = !!rlang::sym("fill"),
          ymin = rect_ymin,
          ymax = rect_ymax
        ),
        alpha = alpha,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = region_colors,
        name = "Region"
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = !!rlang::sym("ci_lower"),
          ymax = !!rlang::sym("ci_upper")
        ),
        fill = line.color, alpha = alpha
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = !!rlang::sym("slope")),
        color = line.color,
        linewidth = line.size
      ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey40", linewidth = 0.4)

    if (length(jn_in_range) > 0L) {
      p_out <- p_out +
        ggplot2::geom_vline(
          xintercept = jn_in_range,
          linetype = "dotted",
          color = "grey30",
          linewidth = 0.5
        )
    }

    if (isTRUE(rug) && !is.null(modx_obs)) {
      rug_df <- data.frame(modx_rug = modx_obs[!is.na(modx_obs)])
      p_out <- p_out +
        ggplot2::geom_rug(
          data = rug_df,
          ggplot2::aes(x = !!rlang::sym("modx_rug")),
          sides = "b",
          alpha = 0.3,
          color = "grey30",
          inherit.aes = FALSE
        )
    }

    x_lab <- x.label %||% modx
    y_lab <- y.label %||% paste0("Conditional slope of ", pred)

    p_out <- p_out +
      ggplot2::labs(x = x_lab, y = y_lab) +
      ggplot2::scale_x_continuous(expand = c(0, 0))
  }

  res <- list(
    jn_points = jn_points,
    signif_regions = signif_regions,
    plot = p_out,
    observed_support = observed_support,
    interaction_test = interaction_test,
    labels = list(outcome = outcome,
                  pred = pred,
                  modx = modx,
                  interaction = interaction),
    conf.level = conf.level,
    digits = digits,
    call = match.call()
  )
  if (isTRUE(return_data)) res$plot_data <- plot_df

  class(res) <- "lav_jn"
  return(res)
}


### metodi S3 print e summary:

#' @param x A 'lav_jn' object.
#' @param ... Additional arguments; unused.
#' @rdname lav_jn
#' @export
print.lav_jn <- function(x, ...) {
  lab <- x$labels
  it <- x$interaction_test
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

  cat("\nInteraction effect (", lab$pred, " * ", lab$modx,
      " -> ", lab$outcome, "): ",
      "b = ", fmt(it$b), ", SE = ", fmt(it$se),
      ", beta = ", fmt(it$beta_std), ", z = ", zfmt(it$z),
      ", p = ", fmtp(it$p),
      ", ", sprintf("%.0f%% CI", cl * 100),
      " [", fmt(it$ci[1L]), ", ", fmt(it$ci[2L]), "]\n",
      sep = "")

  cat("\nJohnson-Neyman boundaries (alpha = ",
      formatC(1 - cl, format = "f", digits = 2), ")\n", sep = "")
  cat(rep("-", 72), "\n", sep = "")

  obs_sup <- x$observed_support
  jn <- x$jn_points
  jn_in <- jn[jn >= obs_sup[1L] & jn <= obs_sup[2L]]

  if (length(jn) == 0L) {
    test_sig <- x$signif_regions$significance[1L] == "significant"
    if (test_sig) {
      cat("No JN boundary found. The conditional slope of ", lab$pred,
          " is significant across the entire observed range of ", lab$modx,
          " [", fmt(obs_sup[1L]), ", ", fmt(obs_sup[2L]), "].\n", sep = "")
    } else {
      cat("No JN boundary found. The conditional slope of ", lab$pred,
          " is non-significant across the entire observed range of ", lab$modx,
          " [", fmt(obs_sup[1L]), ", ", fmt(obs_sup[2L]), "].\n", sep = "")
    }
  } else {
    if (length(jn_in) > 0L) {
      for (j in seq_along(jn_in)) {
        cat("  JN boundary ", j, ": ", lab$modx, " = ", fmt(jn_in[j]),
            "\n", sep = "")
      }
    } else {
      cat("Both roots fall outside the observed range.\n")
    }
    cat("Observed range of ", lab$modx, ": [",
        fmt(obs_sup[1L]), ", ", fmt(obs_sup[2L]), "]\n", sep = "")
  }

  sr <- x$signif_regions
  if (nrow(sr) > 0L) {
    cat("\nSignificance regions (within observed range):\n")
    for (k in seq_len(nrow(sr))) {
      cat("  ", lab$modx, " in [", fmt(sr$lower[k]),
          ", ", fmt(sr$upper[k]), "]: slope is ",
          sr$significance[k], "\n", sep = "")
    }
  }
  cat(rep("-", 72), "\n", sep = "")

  cat("\n")
  if (!is.null(x$plot)) print(x$plot)
  invisible(x)
}

#' @param object A 'lav_jn' object.
#' @param ... Additional arguments; unused.
#' @rdname lav_jn
#' @export
summary.lav_jn <- function(object, ...) {
  print.lav_jn(object, ...)
  invisible(object)
}
