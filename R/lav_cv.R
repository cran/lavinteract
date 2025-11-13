#' Repeated holdout (Monte Carlo) cross-validation of R^2 for structural equation models ('lavaan' objects)
#'
#' Estimate out-of-sample predictive performance for structural relations in a
#' fitted 'lavaan' model using repeated holdout (Monte Carlo cross-validation,
#' leave-group-out CV). At each repetition, the model is refitted on a random
#' training subset and evaluated on a disjoint test subset.
#'
#' For observed outcomes, R^2 is computed by comparing test-set observed values
#' with predictions obtained by applying the training-set structural coefficients
#' to the test-set predictors.
#' 
#' For latent outcomes, the outcome is not directly observed in the test set.
#' Factor scores for the outcome are first computed in the test set using the
#' measurement model learned on the training set; these scores serve as the
#' outcome values. Predictions are then formed by applying the training-set
#' structural coefficients to the test-set predictors (including factor scores
#' for any latent predictors). R^2 is computed by comparing the test-set factor
#' scores of the outcome with these predicted scores.
#' 
#' The in-sample baseline R^2 is computed on the full dataset using the same
#' metric as in cross-validation: observed outcomes use observed-versus-predicted
#' R^2; latent outcomes use score-versus-predicted-score R^2.
#' 
#' By default, repetitions continue until the running mean R^2 for each outcome
#' stabilizes within a specified tolerance over a trailing window of successful
#' splits, or until a maximum number of splits is reached.
#' 
#' The summary table reports the in-sample baseline R^2, the median cross-validated
#' R^2, its standard deviation, and the percent drop (baseline vs. median CV) with
#' heuristic threshold markers. The percent drop is suppressed when the in-sample
#' R^2 is very small. 
#'
#' @usage
#' lav_cv(
#'   fit,
#'   data = NULL,
#'   times = "auto",
#'   train_prop = 0.8,
#'   seed = 42L,
#'   quiet = TRUE,
#'   digits = 3L,
#'   plot = TRUE,
#'   tol = 0.001,
#'   window = 50L,
#'   max_times = 3000L,
#'   min_r2_for_pct = 0.05
#' )
#'
#' @param fit A fitted 'lavaan' object (required).
#' @param data The data frame used to fit the model; if NULL, it is extracted from 'fit' when available (default: NULL).
#' @param times Integer indicating the number of random splits, or "auto" for stabilization-based early stopping (default: "auto").
#' @param train_prop Numeric in (0,1). Proportion of cases in the training split for each repetition (default: 0.8).
#' @param seed Integer. Random seed for reproducibility of the splits (default: 42).
#' @param quiet Logical. Suppress 'lavaan' refit messages when TRUE (default: TRUE).
#' @param digits Integer. Number of digits to print in summaries (default: 3).
#' @param plot Logical. Show convergence plots of the running mean R^2 per outcome (default: TRUE).
#' @param tol Numeric. Tolerance for the auto-stop rule on the running mean (default: 0.001).
#' @param window Integer. Trailing window size (number of successful splits) used by the auto-stop rule (default: 50).
#' @param max_times Integer. Maximum number of splits when \code{times} = "auto" (default: 3000).
#' @param min_r2_for_pct Numeric in (0,1). Minimum in-sample R^2 required to compute percent drop; below this, \%_drop is set to NA (default: 0.05).
#'
#' @return A list with class 'lav_cv' and elements:
#' \describe{
#'   \item{\code{table}}{Data frame with columns:
#'     \code{outcome}, \code{type} ("observed" or "latent"),
#'     \code{r2_in}, \code{r2_cv_mean}, \code{r2_cv_median}, \code{r2_cv_sd},
#'     \code{drop_mean_pct}, \code{drop_med_pct}, \code{splits_used}.}
#'   \item{\code{split_matrix}}{Matrix of split-wise test-set R^2 values (rows = splits, columns = outcomes).}
#'   \item{\code{times}}{Character or integer indicating the number of splits used (e.g., \code{"auto(534)"} or \code{500}).}
#'   \item{\code{train_prop}}{Numeric. Training proportion used in each split.}
#'   \item{\code{N}}{Integer. Number of rows in the input data.}
#'   \item{\code{seed}}{Integer. Random seed used to generate the splits.}
#'   \item{\code{tol}}{Numeric. Tolerance used by the auto-stop rule.}
#'   \item{\code{window}}{Integer. Trailing window size for the auto-stop rule.}
#'   \item{\code{min_r2_for_pct}}{Numeric. Minimum in-sample R^2 required to compute percent drop.}
#'   \item{\code{call}}{\code{match.call()} of the function call.}
#'   \item{\code{digits}}{Integer. Default number of digits for printing.}
#' }
#' 
#' @seealso \code{\link[lavaan]{sem}}, \code{\link[lavaan]{lavPredict}},
#'   \code{\link[lavaan]{inspect}}
#'   
#' @references 
#' Cudeck, R., & Browne, M. W. (1983). Cross-Validation Of Covariance Structures. Multivariate Behavioral Research, 18(2), 147-167. \doi{10.1207/s15327906mbr1802_2}
#' 
#' Hastie, T., Friedman, J., & Tibshirani, R. (2001). The Elements of Statistical Learning. In Springer Series in Statistics. Springer New York. \doi{10.1007/978-0-387-21606-5} 
#'    
#' Kvalseth, T. O. (1985). Cautionary Note about R2. The American Statistician, 39(4), 279-285. \doi{10.1080/00031305.1985.10479448}   
#'    
#' Shmueli, G. (2010). To Explain or to Predict? Statistical Science, 25(3). \doi{10.1214/10-sts330}
#'  
#' Yarkoni, T., & Westfall, J. (2017). Choosing Prediction Over Explanation in Psychology: Lessons From Machine Learning. Perspectives on Psychological Science, 12(6), 1100-1122. \doi{10.1177/1745691617693393}
#'
#' @examplesIf requireNamespace("lavaan", quietly = TRUE)
#' library("lavaan")
#' model <- "
#' ind60 =~ x1 + x2 + x3
#' dem60 =~ y1 + y2 + y3 + y4
#' dem65 =~ y5 + y6 + y7 + y8
#' 
#' dem60 ~ ind60
#' dem65 ~ ind60 + dem60
#' 
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
#' result <- lav_cv(
#' fit = fit, 
#' data = lavaan::PoliticalDemocracy, 
#' times = 5)
#' print(result)
#'
#' @export
lav_cv <- function(
    fit,
    data = NULL,
    times = "auto",
    train_prop = 0.8,
    seed = 42L,
    quiet = TRUE,
    digits = 3L,
    plot = TRUE,
    tol = 0.001,
    window = 50L,
    max_times = 3000L,
    min_r2_for_pct = 0.05   
) {
  if (!inherits(fit, "lavaan")) stop("`fit` must be a fitted lavaan object.", call. = FALSE)
  if (!(is.numeric(times) || identical(times, "auto"))) stop("`times` must be integer or 'auto'.", call. = FALSE)
  if (train_prop <= 0 || train_prop >= 1) stop("`train_prop` must be in (0,1).", call. = FALSE)
  OPT <- fit@Options
  if (!is.null(OPT$group) && nzchar(OPT$group)) stop("Multigroup models are not supported yet.", call. = FALSE)
  
  data_in <- data %||%
    tryCatch(lavaan::lavInspect(fit, "data"), error = function(e) NULL) %||%
    tryCatch(lavaan::lavInspect(fit, "data.original"), error = function(e) NULL)
  if (is.null(data_in)) stop("Could not locate data; pass it via `data =`.", call. = FALSE)
  if (is.matrix(data_in)) data_in <- as.data.frame(data_in)
  
  pt <- lavaan::parTable(fit)
  is_struct <- pt$op == "~" & pt$rhs != "1"
  is_interaction <- grepl(":", pt$rhs, fixed = TRUE) | pt$op == "XWITH"
  struct <- pt[is_struct & !is_interaction, c("lhs", "rhs")]
  outcomes_all <- unique(struct$lhs)
  if (!length(outcomes_all)) stop("No structural regressions found.", call. = FALSE)
  
  latents_all <- tryCatch(lavaan::lavNames(fit, type = "lv"), error = function(e) character(0))
  outcome_is_latent <- outcomes_all %in% latents_all
  any_latents <- length(latents_all) > 0L
  
  estimator_use <- "ML" # vabbe non serve ancora ma ok
  missing_opt <- OPT$missing
  std_lv <- isTRUE(OPT$std.lv)
  parameterization <- OPT$parameterization
  ordered_vars <- tryCatch(lavaan::lavNames(fit, type = "ov.ord"), error = function(e) character(0))
  model_obj <- lavaan::parTable(fit)
  
  r2_fun <- function(y, yhat) {
    ok <- is.finite(y) & is.finite(yhat)
    y <- y[ok]; yhat <- yhat[ok]
    if (length(y) < 3L) return(NA_real_)
    ss_res <- sum((y - yhat)^2)
    ss_tot <- sum((y - mean(y))^2)
    if (ss_tot <= 0) return(NA_real_)
    1 - ss_res / ss_tot
  }
  get_struct_coefs <- function(f) {
    pe <- lavaan::parameterEstimates(f, standardized = FALSE)
    out <- lapply(outcomes_all, function(y) {
      b  <- pe[pe$op == "~" & pe$lhs == y, c("rhs", "est")]
      ic <- pe[pe$op == "~1" & pe$lhs == y, "est"]; if (!length(ic)) ic <- 0
      list(betas = stats::setNames(b$est, b$rhs), intercept = as.numeric(ic))
    })
    names(out) <- outcomes_all
    out
  }
  compute_lv_scores <- function(fit_obj, newdata) {
    if (!any_latents) return(NULL)
    fs <- tryCatch(lavaan::lavPredict(fit_obj, newdata = newdata, type = "lv", method = "regression"),
                   error = function(e) NULL)
    if (is.null(fs)) fs <- tryCatch(lavaan::lavPredict(fit_obj, newdata = newdata, type = "lv", method = "Bartlett"),
                                    error = function(e) NULL)
    fs
  }
  build_X <- function(pred_names, fs, newdata) {
    M <- NULL
    if (any_latents && !is.null(fs)) {
      X_lat <- intersect(pred_names, colnames(fs))
      if (length(X_lat)) M <- cbind(M, fs[, X_lat, drop = FALSE])
    }
    X_obs <- intersect(pred_names, colnames(newdata))
    if (length(X_obs)) M <- cbind(M, as.matrix(newdata[, X_obs, drop = FALSE]))
    M
  }
  
  # estrazione R^2 di baseline
  coefs_full <- get_struct_coefs(fit)
  fs_full <- if (any_latents) compute_lv_scores(fit, data_in) else NULL
  baseline_r2 <- stats::setNames(rep(NA_real_, length(outcomes_all)), outcomes_all) # ok call diretta da namespace
  for (y in outcomes_all) {
    cY <- coefs_full[[y]]
    preds <- names(cY$betas)
    X <- build_X(preds, fs_full, data_in)
    if (!is.null(X)) {
      yhat <- drop(cY$intercept + X %*% cY$betas[colnames(X)])
      if (!outcome_is_latent[match(y, outcomes_all)]) {
        if (y %in% colnames(data_in)) baseline_r2[y] <- r2_fun(data_in[[y]], yhat)
      } else if (!is.null(fs_full) && y %in% colnames(fs_full)) {
        baseline_r2[y] <- r2_fun(fs_full[, y], yhat)
      }
    }
  }
  
  # aiuto, funziona questo! lascia cosi
  set.seed(seed)
  N <- nrow(data_in); n_train <- max(3L, floor(train_prop * N))
  maxN <- if (is.numeric(times)) as.integer(times) else as.integer(max_times)
  res_mat <- matrix(NA_real_, nrow = maxN, ncol = length(outcomes_all))
  colnames(res_mat) <- outcomes_all
  used <- 0L
  done <- function(vmat) {
    for (j in seq_len(ncol(vmat))) {
      v <- vmat[, j]; v <- v[is.finite(v)]
      if (length(v) < window + 5L) return(FALSE)
      cm <- cumsum(v)/seq_along(v)
      if (abs(cm[length(cm)] - cm[length(cm) - window]) > tol) return(FALSE)
    }
    TRUE
  }
  target_splits <- if (is.numeric(times)) as.integer(times) else as.integer(max_times)
  
  while (used < target_splits) {
    tr_idx <- sort(sample.int(N, size = n_train, replace = FALSE))
    te_idx <- setdiff(seq_len(N), tr_idx); if (length(te_idx) < 3L) next
    dat_train <- data_in[tr_idx, , drop = FALSE]
    dat_test <- data_in[te_idx,  , drop = FALSE]
    
    fit_train <- tryCatch(
      lavaan::sem(
        model = model_obj, data = dat_train,
        estimator = estimator_use, missing = missing_opt,
        se = "none", test = "none",
        std.lv = std_lv, parameterization = parameterization,
        ordered = if (length(ordered_vars)) ordered_vars else NULL,
        warn = !quiet
      ),
      error = function(e) NULL
    )
    if (is.null(fit_train) || !isTRUE(lavaan::inspect(fit_train, "converged"))) next
    
    used <- used + 1L
    
    coefs_tr <- get_struct_coefs(fit_train)
    fs_te <- if (any_latents) compute_lv_scores(fit_train, dat_test) else NULL
    
    for (j in seq_along(outcomes_all)) {
      y <- outcomes_all[j]
      cY <- coefs_tr[[y]]
      preds <- names(cY$betas)
      X <- build_X(preds, fs_te, dat_test)
      if (is.null(X)) { res_mat[used, j] <- NA_real_; next }
      yhat <- drop(cY$intercept + X %*% cY$betas[colnames(X)])
      if (!outcome_is_latent[j]) {
        if (y %in% colnames(dat_test)) res_mat[used, j] <- r2_fun(dat_test[[y]], yhat)
      } else {
        if (!is.null(fs_te) && y %in% colnames(fs_te)) res_mat[used, j] <- r2_fun(fs_te[, y], yhat)
      }
    }
    
    if (identical(times, "auto") && used >= (window + 5L)) {
      if (done(res_mat[seq_len(used), , drop = FALSE])) break
    }
  }
  
  if (used < nrow(res_mat)) res_mat <- res_mat[seq_len(used), , drop = FALSE]
  
# estrai sommari (esplicitare stats senno fa casino non so perche)
  r2_cv_mean <- apply(res_mat, 2L, function(v) { v <- v[is.finite(v)]; if (!length(v)) NA_real_ else mean(v) })
  r2_cv_median <- apply(res_mat, 2L, function(v) { v <- v[is.finite(v)]; if (!length(v)) NA_real_ else stats::median(v) })
  r2_cv_sd <- apply(res_mat, 2L, function(v) { v <- v[is.finite(v)]; if (!length(v)) NA_real_ else stats::sd(v) })
  n_ok <- apply(res_mat, 2L, function(v) sum(is.finite(v)))
  
  drop_mean_pct   <- rep(NA_real_, length(outcomes_all))
  drop_median_pct <- rep(NA_real_, length(outcomes_all))
  for (j in seq_along(outcomes_all)) {
    if (is.finite(baseline_r2[j]) && baseline_r2[j] >= min_r2_for_pct) {
      drop_mean_pct[j]   <- 100 * (baseline_r2[j] - r2_cv_mean[j])   / baseline_r2[j]
      drop_median_pct[j] <- 100 * (baseline_r2[j] - r2_cv_median[j]) / baseline_r2[j]
    } 
  }
  
  # ok maybe it works...spero v5
  out_type <- ifelse(outcome_is_latent, "latent", "observed")
  res_table <- data.frame(
    outcome = outcomes_all,
    type = out_type,
    r2_in = as.numeric(baseline_r2),
    r2_cv_mean = as.numeric(r2_cv_mean),
    r2_cv_median = as.numeric(r2_cv_median),
    r2_cv_sd = as.numeric(r2_cv_sd),
    drop_mean_pct = as.numeric(drop_mean_pct),
    drop_med_pct = as.numeric(drop_median_pct),
    splits_used = as.integer(n_ok),
    stringsAsFactors = FALSE
  )
  
  res <- list(
    table = res_table[order(res_table$type, res_table$outcome), ],
    split_matrix = res_mat,
    times = if (identical(times, "auto")) paste0("auto(", used, ")") else as.integer(times),
    train_prop = train_prop,
    N = nrow(data_in),
    seed = seed,
    digits = digits,
    call = match.call(),
    tol = tol,
    window = window,
    min_r2_for_pct = min_r2_for_pct
  )
  class(res) <- c("lav_cv", "list")
  
  # devi indicare graphics::!
  if (isTRUE(plot)) {
    op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)
    for (y in res$table$outcome) {
      v <- res$split_matrix[, y]; v <- v[is.finite(v)]
      if (!length(v)) next
      m  <- seq_along(v)
      cm <- cumsum(v)/m
      s2 <- c(NA_real_, cumsum((v[-1] - cm[-length(cm)]) * (v[-1] - cm[-1])) / pmax(1, m[-1] - 1))
      se <- sqrt(s2 / m); se[1] <- NA_real_
      upper <- cm + 1.96 * se
      lower <- cm - 1.96 * se
      graphics::plot(m, cm, type = "l", xlab = "Splits used", ylab = "Running mean of R^2",
           main = paste0("CV convergence: ", y))
      graphics::lines(m, upper, lty = 2)
      graphics::lines(m, lower, lty = 2)
      graphics::abline(h = utils::tail(cm, 1L), lty = 3) # mancava la call da utils!!!
    }
  }
  
  res
}

### metodi S3 print e summary:

#' @param x A 'lav_cv' object.
#' @param ... Additional arguments; unused.
#' @rdname lav_cv
#' @export
print.lav_cv <- function(x, digits = x$digits %||% 3L, ...) {
  stopifnot(is.list(x), !is.null(x$table))
  u2 <- if (isTRUE(l10n_info()[["UTF-8"]])) "\u00B2" else "^2"
  
  # header con i dettagli giusti!!
  auto_txt <- if (grepl("^auto\\(", as.character(x$times))) {
    paste0("auto-stop at ", gsub("[^0-9]", "", x$times), " splits")
  } else {
    paste0(x$times, " splits")
  }
  cat(sprintf(
    "\nRepeated holdout (Monte Carlo) cross-validation of R%s\n", u2
  ))
  cat(sprintf("Splits: %s (train proportion = %.2f; N = %d; seed = %s)\n\n",
              auto_txt, x$train_prop, x$N, as.character(x$seed)))
  
  tbl <- x$table
  
  star_fun <- function(d){
    if (is.na(d)) "" else if (d >= 50) "***" else if (d >= 30) "**" else if (d >= 10) "*" else ""
  }
  stars <- vapply(tbl$drop_med_pct, star_fun, character(1))
  
  dash <- "-"  # al posto di NA meglio un trattino
  fmt  <- function(z, d = digits) {
    ifelse(is.na(z), dash, formatC(z, digits = d, format = "f"))
  }
  fmt1 <- function(z) {
    ifelse(is.na(z), dash, formatC(z, digits = 1, format = "f"))
  }
  
  for (tp in c("observed", "latent")) {
    sub <- tbl[tbl$type == tp, , drop = FALSE]
    if (!nrow(sub)) next
    
    out_names <- sub$outcome
    r2_infit <- round(sub$r2_in, digits)
    r2_med <- round(sub$r2_cv_median, digits)
    sd_cv <- round(sub$r2_cv_sd, digits)
    drop_pct <- ifelse(is.finite(sub$drop_med_pct), round(sub$drop_med_pct, 1), NA)
    stars_tp <- stars[tbl$type == tp]
    
    cat(sprintf("%s outcomes\n", tools::toTitleCase(tp)))
    cat(rep("-", 86), "\n", sep = "")
    cat(sprintf("%-16s %14s %16s %10s %9s %3s\n",
                "Outcome", paste0("R", u2, " in-sample"),
                paste0("CV R", u2), "CV SD", "% drop", ""))
    cat(rep("-", 86), "\n", sep = "")
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("%-16s %14s %16s %10s %9s %3s\n",
                  out_names[i],
                  fmt(r2_infit[i]),
                  fmt(r2_med[i]),
                  fmt(sd_cv[i]),
                  fmt1(drop_pct[i]),
                  stars_tp[i]))
    }
    cat("\n")
  }
  
  # notine a pie tabella
  cat(sprintf("Auto-stop rule: window = %d, tol = %.3f (applied when times = 'auto').\n", x$window, x$tol))
  cat("Stars indicate % drop thresholds: * >=10%, ** >=30%, *** >=50%.\n")
  cat(sprintf("%% drop compares in-sample R%s to median CV R%s; suppressed when in-sample R%s < %.2f.\n",
              u2, u2, u2, x$min_r2_for_pct %||% 0.05))
  cat(sprintf("CV SD is the standard deviation of split-wise R%s across successful splits.\n", u2))
  invisible(x)
}

#' @param object A 'lav_cv' object.
#' @param ... Additional arguments; unused.
#' @rdname lav_cv
#' @export
summary.lav_cv <- function(object, ...) { 
  print.lav_cv(object, ...); 
  invisible(object) 
  }
