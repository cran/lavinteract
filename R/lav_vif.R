#' Variance Inflation Factors for 'lavaan' Structural Predictors
#'
#' Compute VIF for each predictor that appears in structural regressions
#' with two or more predictors, refitting the necessary sub-models so that
#' latent predictors are handled at the latent level (i.e., with their original
#' measurement models). It returns also the \eqn{R^2} of each eligible endogenous variable 
#' from the original fit for context. 
#' 
#' Each auxiliary refitted model includes the measurement structure for latent 
#' predictors, rebuilt from the original model syntax, together with observed 
#' residual covariances among the involved indicators when specified in the 
#' original model. Then, regresses the focal predictor on the remaining predictors 
#' at the latent level when applicable.
#' 
#' \eqn{\mathrm{VIF}_i = 1 / (1 - R_i^2)} generalizes VIF to SEM while respecting 
#' measurement models.
#'   
#' The function reuses the estimator, missing-data handling, and several options 
#' from \code{fit}. 
#' In multigroup models, VIF is computed separately within each group. 
#' Higher-order latent predictors are supported by recursively rebuilding
#' lower-order measurement blocks in the auxiliary refitted models.
#' 
#' @usage
#' lav_vif(
#'   fit,
#'   data = NULL,
#'   quiet = TRUE
#' )
#'
#' @param fit A fitted \code{lavaan} object.
#' @param data Optional. The data frame used to fit \code{fit}. If \code{NULL},
#'   the function attempts to extract the data from \code{fit} via
#'   \code{lavInspect(fit, "data")}. For multigroup models, grouped data are
#'   reconstructed into a single data frame when possible.
#' @param quiet Logical. If \code{TRUE} suppresses lavaan refit messages.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{vif_table}: data.frame with columns
#'         \code{outcome}, \code{predictor}, \code{group}, \code{r2_predictor},
#'         \code{vif}, \code{k_predictors}.
#'   \item \code{outcome_r2}: data.frame with \eqn{R^2} per eligible endogenous
#'         \code{outcome} and \code{group} from the original \code{fit}.
#' }
#'
#' @references
#' Fox, J., & Monette, G. (1992). Generalized collinearity diagnostics.
#' \emph{Journal of the American Statistical Association}, \emph{87}(417),
#' 178-183. \doi{10.1080/01621459.1992.10475190}
#'
#' Belsley, D. A., Kuh, E., & Welsch, R. E. (1980).
#' \emph{Regression diagnostics: Identifying influential data and sources of
#' collinearity}. Wiley.
#'
#' @examples
#' set.seed(42)
#' x1 <- rnorm(100); x2 <- 0.85*x1 + rnorm(100, sd = sqrt(1 - 0.85^2)); x3 <- rnorm(100)
#' y <- 0.5*x1 + 0.3*x2 + 0.1*x3 + rnorm(100, sd = 0.7)
#' dataset <- data.frame(y, x1, x2, x3)
#' fit <- lavaan::sem("y ~ x1 + x2 + x3", data = dataset)
#' lav_vif(
#' fit = fit,
#' data = dataset)
#' 
#' @importFrom rlang %||%
#' @export
lav_vif <- function(fit, data = NULL, quiet = TRUE) {
  if (!inherits(fit, "lavaan")) stop("`fit` must be a fitted lavaan object.")
  
  # trova i dati e lascia errori se non va
  group_labels <- tryCatch(lavaan::lavInspect(fit, "group.label"), error = function(e) NULL)
  group_var <- NULL
  group_call <- tryCatch(fit@call$group, error = function(e) NULL)
  if (!is.null(group_call)) {
    group_var <- paste(as.character(group_call), collapse = "")
    group_var <- gsub('^"|"$', "", group_var)
  }
  if (is.null(group_var) || length(group_var) != 1L || is.na(group_var) || !nzchar(group_var)) {
    tmp_group <- tryCatch(fit@Options$group, error = function(e) NULL)
    if (length(tmp_group) == 1L && !is.na(tmp_group) && nzchar(tmp_group)) {
      group_var <- tmp_group
    } else {
      group_var <- NULL
    }
  }
  data_in <- data %||%
    tryCatch(lavaan::lavInspect(fit, "data"), error = function(e) NULL)
  if (is.null(data_in)) {
    stop("Could not locate data from `fit`; pass it via `data =`.")
  }
  if (is.data.frame(data_in)) {
    data_in <- as.data.frame(data_in)
  } else if (is.matrix(data_in)) {
    data_in <- as.data.frame(data_in)
    ov_names <- tryCatch(lavaan::lavNames(fit, type = "ov"), error = function(e) NULL)
    if (!is.null(ov_names) && length(ov_names) == ncol(data_in)) {
      names(data_in) <- ov_names
    }
  } else if (is.list(data_in)) {
    ok_list <- all(vapply(data_in, function(x) is.data.frame(x) || is.matrix(x), logical(1)))
    if (!ok_list) {
      stop("Recovered grouped data from `fit`, but could not coerce it to a data frame; please pass it via `data =`.")
    }
    ov_names <- tryCatch(lavaan::lavNames(fit, type = "ov"), error = function(e) NULL)
    data_list <- vector("list", length(data_in))
    for (g in seq_along(data_in)) {
      dg <- data_in[[g]]
      if (is.matrix(dg)) {
        dg <- as.data.frame(dg)
        if (!is.null(ov_names) && length(ov_names) == ncol(dg)) {
          names(dg) <- ov_names
        }
      } else {
        dg <- as.data.frame(dg)
      }
      if (!is.null(group_var) && !(group_var %in% names(dg))) {
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
      data_list[[g]] <- dg
    }
    data_in <- do.call(rbind, data_list)
    rownames(data_in) <- NULL
  } else {
    stop("Recovered data is not a data.frame, matrix, or supported grouped list; please pass it via `data =`.")
  }
  
  # estrai tutto
  OPT <- fit@Options
  estimator <- OPT$estimator
  missing_opt <- OPT$missing
  std_lv <- isTRUE(OPT$std.lv)
  parameterization <- OPT$parameterization
  aux_se <- "none"
  aux_test <- "none"

  ngroups <- tryCatch(
    as.integer(lavaan::lavInspect(fit, "ngroups")),
    error = function(e) 1L
  )
  if (length(ngroups) != 1L || is.na(ngroups)) ngroups <- 1L

  group_labels <- tryCatch(
    lavaan::lavInspect(fit, "group.label"),
    error = function(e) NULL
  )
  
  group_var <- NULL
  
  group_call <- tryCatch(fit@call$group, error = function(e) NULL)
  if (!is.null(group_call)) {
    group_var <- paste(as.character(group_call), collapse = "")
    group_var <- gsub('^"|"$', "", group_var)
  }
  
  if (is.null(group_var) || !nzchar(group_var)) {
    tmp_group <- tryCatch(OPT$group, error = function(e) NULL)
    if (length(tmp_group) == 1L && !is.na(tmp_group) && nzchar(tmp_group)) {
      group_var <- tmp_group
    } else {
      group_var <- NULL
    }
  }
  ordered_vars <- tryCatch(lavaan::lavNames(fit, type = "ov.ord"), error = function(e) character(0))
  
  # mappa strutturale
  pt <- lavaan::parTable(fit)
  is_struct <- pt$op == "~" & pt$rhs != "1"
  is_interaction <- grepl(":", pt$rhs, fixed = TRUE) | pt$op == "XWITH"
  struct <- pt[is_struct & !is_interaction, c("lhs", "rhs", "group", "free")]
  struct$rhs_free <- struct$free != 0
  struct_k <- stats::aggregate(rhs_free ~ lhs + group, data = struct, sum)
  eligible <- struct_k[struct_k$rhs_free >= 2, c("lhs", "group")]
  if (!nrow(eligible)) stop("No endogenous variables with two or more predictors.")
  key <- unique(struct[struct$rhs_free, c("lhs", "rhs", "group")])
  
  # mapping di misura nuovo con anche second'ordine
  latents_all <- tryCatch(lavaan::lavNames(fit, type = "lv"), error = function(e) character(0))
  loadings <- pt[pt$op == "=~", c("lhs", "rhs")]
  
  latent_lhs <- unique(loadings$lhs)
  indicators_by_latent <- lapply(latent_lhs, function(L) {
    unique(loadings$rhs[loadings$lhs == L])})
  names(indicators_by_latent) <- latent_lhs
  
  ov_all <- tryCatch(lavaan::lavNames(fit, type = "ov"), error = function(e) unique(c(pt$lhs, pt$rhs)))
  resid_cov <- pt[pt$op == "~~" & pt$lhs != pt$rhs, c("lhs", "rhs")]
  resid_cov <- resid_cov[resid_cov$lhs %in% ov_all & resid_cov$rhs %in% ov_all, , drop = FALSE]
  espandi_blocco_latenti <- function(lat_set) {
    lat_set <- intersect(unique(lat_set), names(indicators_by_latent))
    if (!length(lat_set)) return(character(0))
    
    out <- character(0)
    seen <- character(0)
    
    walk <- function(L) {
      if (L %in% seen) return(invisible(NULL))
      seen <<- c(seen, L)
      
      rhs <- indicators_by_latent[[L]]
      lower_latents <- intersect(rhs, names(indicators_by_latent))
      
      if (length(lower_latents)) {
        for (j in seq_along(lower_latents)) {
          walk(lower_latents[j])
        }
      }
      
      out <<- c(out, L)
      invisible(NULL)
    }
    
    for (j in seq_along(lat_set)) {
      walk(lat_set[j])
    }
    
    out
  }
  
  build_measurement_syntax <- function(lat_set) {
    lat_set <- espandi_blocco_latenti(lat_set)
    if (!length(lat_set)) return(character(0))
    
    meas_lines <- vapply(lat_set, function(L) {
      inds <- indicators_by_latent[[L]]
      paste0(L, " =~ ", paste(inds, collapse = " + "))
    }, character(1))
    
    used_inds <- unique(unlist(lapply(lat_set, function(L) {
      rhs <- indicators_by_latent[[L]]
      rhs[rhs %in% ov_all]
    }), use.names = FALSE))
    
    rc_lines <- character(0)
    if (length(used_inds) > 1L && nrow(resid_cov)) {
      rc_keep <- resid_cov[
        resid_cov$lhs %in% used_inds & resid_cov$rhs %in% used_inds,
        ,
        drop = FALSE
      ]
      if (nrow(rc_keep)) {
        rc_lines <- unique(paste0(rc_keep$lhs, " ~~ ", rc_keep$rhs))
      }
    }
    
    unique(c(meas_lines, rc_lines))
  } # stop blocco nuovo fixato
  
  # R2 ottenuto dal fit originale!!!
  r2_orig <- lavaan::inspect(fit, "r2")
  outcome_r2 <- list()
  
  if (is.list(r2_orig)) {
    for (i in seq_len(nrow(eligible))) {
      y <- eligible$lhs[i]
      g <- eligible$group[i]
      r2_g <- r2_orig[[g]]
      
      if (length(r2_g) && y %in% names(r2_g)) {
        outcome_r2[[length(outcome_r2) + 1L]] <- data.frame(
          outcome = y,
          group = g,
          r2_outcome = unname(r2_g[y])
        )
      }
    }
  } else {
    for (i in seq_len(nrow(eligible))) {
      y <- eligible$lhs[i]
      
      if (length(r2_orig) && y %in% names(r2_orig)) {
        outcome_r2[[length(outcome_r2) + 1L]] <- data.frame(
          outcome = y,
          group = 1L,
          r2_outcome = unname(r2_orig[y])
        )
      }
    }
  }
  
  outcome_r2 <- if (length(outcome_r2)) {
    do.call(rbind, outcome_r2)
  } else {
    data.frame(
      outcome = character(0),
      group = integer(0),
      r2_outcome = numeric(0)
    )
  }
  
  # VIF ottenuti da regressioni accessorie
  out_rows <- list()
  for (i in seq_len(nrow(eligible))) {
    y <- eligible$lhs[i]; g <- eligible$group[i]
    preds <- unique(key$rhs[key$lhs == y & key$group == g])
    if (length(preds) < 2L) next
    lat_preds <- intersect(preds, latents_all)
    meas_syntax <- build_measurement_syntax(lat_preds)
    
    for (px in preds) {
      if (grepl(":", px, fixed = TRUE)) next
      others <- setdiff(preds, px); if (!length(others)) next
      reg_line <- paste0(px, " ~ ", paste(others, collapse = " + "))
      aux_model <- paste(c(meas_syntax, reg_line), collapse = "\n")
      
      aux_fit <- lavaan::sem(
        model = aux_model,
        data = data_in,
        group = group_var,
        estimator = estimator,
        missing = missing_opt,
        se = aux_se,
        test = aux_test,
        std.lv = std_lv,
        parameterization = parameterization,
        ordered = if (length(ordered_vars)) ordered_vars else NULL,
        warn = !quiet
      )
      
      r2_aux <- lavaan::inspect(aux_fit, "r2")
      
      if (is.list(r2_aux)) {
        if (!is.null(r2_aux[[g]]) && px %in% names(r2_aux[[g]])) {
          r2v <- unname(r2_aux[[g]][px])
        } else {
          r2v <- NA_real_
        }
        
        r2v <- if (is.na(r2v)) NA_real_ else min(max(r2v, 0), 0.999999)
        vifv <- if (is.na(r2v)) NA_real_ else 1 / (1 - r2v)
        
        out_rows[[length(out_rows) + 1L]] <- data.frame(
          outcome = y,
          predictor = px,
          group = g,
          r2_predictor = r2v,
          vif = vifv,
          k_predictors = length(preds)
        )
      } else {
        if (px %in% names(r2_aux)) {
          r2v <- unname(r2_aux[px])
        } else {
          r2v <- NA_real_
        }
        
        r2v <- if (is.na(r2v)) NA_real_ else min(max(r2v, 0), 0.999999)
        vifv <- if (is.na(r2v)) NA_real_ else 1 / (1 - r2v)
        
        out_rows[[length(out_rows) + 1L]] <- data.frame(
          outcome = y,
          predictor = px,
          group = 1L,
          r2_predictor = r2v,
          vif = vifv,
          k_predictors = length(preds)
        )
      }
    }
  }
  
  vif_table <- if (length(out_rows)) do.call(rbind, out_rows) else
    data.frame(outcome = character(0), predictor = character(0), group = integer(0),
               r2_predictor = numeric(0), vif = numeric(0), k_predictors = integer(0))
  if (nrow(vif_table)) {
    vif_table <- vif_table[order(vif_table$outcome, vif_table$group, -vif_table$vif), ]
    rownames(vif_table) <- NULL
  }
  if (nrow(outcome_r2)) {
    keep <- unique(vif_table$outcome)
    outcome_r2 <- outcome_r2[outcome_r2$outcome %in% keep, , drop = FALSE]
    rownames(outcome_r2) <- NULL
  }
  
  res <- list(
    vif_table    = vif_table,
    outcome_r2   = outcome_r2,
    group_var    = group_var,
    group_labels = group_labels,
    call         = match.call()
  )
  class(res) <- c("lav_vif", "list")
  res
}

### metodi S3 print e summary:

#' @param x A 'lav_vif' object.
#' @param digits Integer number of digits to print.
#' @param cutoff Numeric length-2 thresholds used to flag VIF values.
#' @param ... Additional arguments; unused.
#' @rdname lav_vif
#' @export
print.lav_vif <- function(x, digits = 3, cutoff = c(5, 10), ...) {
  stopifnot(is.list(x), !is.null(x$vif_table))
  vt  <- x$vif_table
  or2 <- x$outcome_r2
  gl  <- x$group_labels
  u2  <- if (isTRUE(l10n_info()[["UTF-8"]])) "\u00B2" else "^2"
  
  if (!nrow(vt)) { cat("No eligible structural regressions (k < 2).\n"); return(invisible(x)) }
  
  fmt  <- function(z) ifelse(is.na(z), "NA", formatC(z, digits = digits, format = "f"))
  keys <- unique(vt[, c("outcome", "group")])
  
  for (ii in seq_len(nrow(keys))) {
    y <- keys$outcome[ii]; g <- keys$group[ii]
    gname <- if (length(gl) && !is.na(g) && g <= length(gl)) gl[g] else paste0("Group ", g)
    
    cat(sprintf("Dependent variable: %s  |  Group: %s\n", y, gname))
    
    r2y <- or2$r2_outcome[or2$outcome == y & or2$group == g]
    if (length(r2y)) cat(sprintf("  Model R%s = %s\n\n", u2, fmt(r2y)))
    
    sub <- vt[vt$outcome == y & vt$group == g, c("predictor", "r2_predictor", "vif")]
    sub <- sub[order(-sub$vif), , drop = FALSE]
    for (i in seq_len(nrow(sub))) {
      v <- sub$vif[i]
      flag <- if (is.na(v)) "" else if (v >= cutoff[2]) " [VIF>=10]" else if (v >= cutoff[1]) " [VIF>=5]" else ""
      tol <- if (is.na(v)) NA_real_ else 1 / v
      cat(sprintf("  %s: Predictor R%s = %s, VIF = %s, Tolerance = %s%s\n",
                  sub$predictor[i], u2, fmt(sub$r2_predictor[i]), fmt(v), fmt(tol), flag))
    }
    cat("\n")
  }
  ge <- if (isTRUE(l10n_info()[["UTF-8"]])) "\u2265" else ">="
  cat(sprintf(
    "Note. VIF = 1/(1 - Predictor R%s); Tolerance = 1 - Predictor R%s = 1/VIF. Flags mark VIF %s %s and %s %s.\n",
    u2, u2, ge, cutoff[1], ge, cutoff[2]
  ))
  invisible(x)
}

#' @param object A 'lav_vif' object.
#' @param ... Passed to `print.lav_vif()` (e.g., `digits`, `cutoff`).
#' @rdname lav_vif
#' @export
summary.lav_vif <- function(object, ...) {
  print.lav_vif(object, ...)
  invisible(object)
}
