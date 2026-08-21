#' Repeated holdout (Monte Carlo) cross-validation of \eqn{R^2} for structural equation models ('lavaan' objects)
#'
#' Estimate out-of-sample predictive performance for structural relations in a
#' fitted 'lavaan' model using repeated holdout (Monte Carlo cross-validation,
#' leave-group-out CV). At each repetition, the model is refitted on a random
#' training subset and evaluated on a disjoint test subset. 
#' 
#' A split is retained only if the training-set solution converges and is 
#' admissible (no negative estimated variance and no non-positive-definite 
#' parameter matrix); splits failing either check are discarded and replaced.
#'
#' For observed outcomes, \eqn{R^2} is computed by comparing test-set observed values
#' with predictions obtained by applying the training-set structural coefficients
#' to the test-set predictors.
#' 
#' For latent outcomes, the outcome is not directly observed in the test set.
#' Factor scores for the outcome are first computed in the test set using the
#' measurement model learned on the training set; these scores serve as the
#' outcome values. Predictions are then formed by applying the training-set
#' structural coefficients to the test-set predictors (including factor scores
#' for any latent predictors). \eqn{R^2} is computed by comparing the test-set factor
#' scores of the outcome with these predicted scores. The score method can be 
#' chosen between Bartlett and regression scores.
#' 
#' The in-sample baseline \eqn{R^2} is computed on the full dataset using the same
#' metric as in cross-validation: observed outcomes use observed-versus-predicted
#' \eqn{R^2}; latent outcomes use score-versus-predicted-score R^2.
#' 
#' By default, repetitions continue until the running mean \eqn{R^2} for each outcome
#' stabilizes within a specified tolerance over a trailing window of successful
#' splits, or until a maximum number of splits is reached.
#' 
#' The summary table reports the in-sample baseline \eqn{R^2}, the selected mean or median
#' cross-validated \eqn{R^2}, its standard deviation, the absolute overfitting index
#' \eqn{\Delta} (in-sample minus cross-validated \eqn{R^2}), and the proportional
#' overfitting index \eqn{\Delta\%} = 100 * (in-sample \eqn{R^2} - CV \eqn{R^2}) /
#' in-sample \eqn{R^2}. 
#' \eqn{\Delta\%} is near zero when the structural relation generalizes, and 
#' exceeds 100% exactly when the cross-validated \eqn{R^2} is negative, that is, 
#' when the in-sample \eqn{R^2} does not generalize out of sample. 
#' \eqn{\Delta\%} is suppressed when the in-sample \eqn{R^2} is very small.
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
#'   min_r2_for_pct = 0.05,
#'   scores = c("bartlett", "regression"),
#'   aggregation = c("median", "mean")
#' )
#'
#' @param fit A fitted 'lavaan' object (required).
#' @param data The data frame used to fit the model; if NULL, it is extracted from 'fit' when available (default: NULL).
#' @param times Integer indicating the number of random splits, or "auto" for stabilization-based early stopping (default: "auto").
#' @param train_prop Numeric in (0,1). Proportion of cases in the training split for each repetition (default: 0.8).
#' @param seed Integer. Random seed for reproducibility of the splits (default: 42).
#' @param quiet Logical. Suppress 'lavaan' refit messages when TRUE (default: TRUE).
#' @param digits Integer. Number of digits to print in summaries (default: 3).
#' @param plot Logical. Show convergence plots of the running mean \eqn{R^2} per outcome (default: TRUE).
#' @param tol Numeric. Tolerance for the auto-stop rule on the running mean (default: 0.001).
#' @param window Integer. Trailing window size (number of successful splits) used by the auto-stop rule (default: 50).
#' @param max_times Integer. Maximum number of splits when \code{times} = "auto" (default: 3000).
#' @param min_r2_for_pct Numeric in (0,1). Minimum in-sample \eqn{R^2} required to compute the proportional index \eqn{\Delta\%}; below this, \eqn{\Delta\%} is set to NA (default: 0.05).
#' @param scores Character. Factor-score method for latent variables: \code{"bartlett"}
#'   (default) or \code{"regression"}. The selected method is used for latent outcomes
#'   and latent predictors.
#' @param aggregation Character. Split-wise \eqn{R^2} summary used in printed output and
#'   overfitting-index calculations: \code{"median"} (default) or \code{"mean"}. Both summaries
#'   remain available in the returned \code{table}. This option does not affect the
#'   running-mean auto-stop rule.
#'
#' @return A list with class 'lav_cv' and elements:
#' \describe{
#'   \item{\code{table}}{Data frame with columns:
#'     \code{outcome}, \code{type} ("observed" or "latent"),
#'     \code{r2_in}, \code{r2_cv_mean}, \code{r2_cv_median}, \code{r2_cv_sd},
#'     \code{drop_mean_abs}, \code{drop_med_abs} (absolute index \eqn{\Delta}),
#'     \code{drop_mean_pct}, \code{drop_med_pct} (proportional index \eqn{\Delta\%}),
#'     \code{splits_used}.}
#'   \item{\code{split_matrix}}{Matrix of split-wise test-set \eqn{R^2} values (rows = splits, columns = outcomes).}
#'   \item{\code{times}}{Character or integer indicating the number of splits used (e.g., \code{"auto(534)"} or \code{500}).}
#'   \item{\code{train_prop}}{Numeric. Training proportion used in each split.}
#'   \item{\code{N}}{Integer. Number of rows in the input data.}
#'   \item{\code{seed}}{Integer. Random seed used to generate the splits.}
#'   \item{\code{splits_attempted}}{Integer. Total number of splits attempted, including those discarded for non-convergence or inadmissibility.}
#'   \item{\code{tol}}{Numeric. Tolerance used by the auto-stop rule.}
#'   \item{\code{window}}{Integer. Trailing window size for the auto-stop rule.}
#'   \item{\code{min_r2_for_pct}}{Numeric. Minimum in-sample \eqn{R^2} required to compute \eqn{\Delta\%}.}
#'   \item{\code{scores}}{Character. Factor-score method used for latent variables.}
#'   \item{\code{aggregation}}{Character. Split-wise \eqn{R^2} summary selected for printing and overfitting-index calculations.}
#'   \item{\code{call}}{\code{match.call()} of the function call.}
#'   \item{\code{digits}}{Integer. Default number of digits for printing.}
#' }
#' 
#' @seealso \code{\link[lavaan]{sem}}, \code{\link[lavaan]{lavPredict}},
#'   \code{\link[lavaan]{lavInspect}}
#'   
#' @references
#' Bollen, K. A., & Stine, R. A. (1992). Bootstrapping goodness-of-fit measures
#' in structural equation models. \emph{Sociological Methods & Research},
#' \emph{21}(2), 205-229. \doi{10.1177/0049124192021002004}
#'
#' Cudeck, R., & Browne, M. W. (1983). Cross-validation of covariance
#' structures. \emph{Multivariate Behavioral Research}, \emph{18}(2), 147-167.
#' \doi{10.1207/s15327906mbr1802_2}
#'
#' Hastie, T., Tibshirani, R., & Friedman, J. (2001). \emph{The elements of
#' statistical learning: Data mining, inference, and prediction}. Springer.
#' \doi{10.1007/978-0-387-21606-5}
#'
#' Kvålseth, T. O. (1985). Cautionary note about \eqn{R^2}.
#' \emph{The American Statistician}, \emph{39}(4), 279-285.
#' \doi{10.1080/00031305.1985.10479448}
#'
#' Shmueli, G. (2010). To explain or to predict?
#' \emph{Statistical Science}, \emph{25}(3), 289-310.
#' \doi{10.1214/10-STS330}
#'
#' Yarkoni, T., & Westfall, J. (2017). Choosing prediction over explanation in
#' psychology: Lessons from machine learning. \emph{Perspectives on
#' Psychological Science}, \emph{12}(6), 1100-1122.
#' \doi{10.1177/1745691617693393}
#'
#' @examples 
#' set.seed(42)
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
#' \donttest{
#' fit <- lavaan::sem(
#' model = model, 
#' data = lavaan::PoliticalDemocracy,
#' std.lv = TRUE, 
#' estimator = "MLR", 
#' meanstructure = TRUE)
#' result <- lav_cv(
#' fit = fit, 
#' data = lavaan::PoliticalDemocracy, 
#' times = 100)
#' print(result)
#' }
#'
#' @importFrom rlang %||%
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
    min_r2_for_pct = 0.05,
    scores = c("bartlett", "regression"),
    aggregation = c("median", "mean")
) {
  if (!inherits(fit, "lavaan")) stop("`fit` must be a fitted lavaan object.", call. = FALSE)
  if (!(is.numeric(times) || identical(times, "auto"))) stop("`times` must be integer or 'auto'.", call. = FALSE)
  if (train_prop <= 0 || train_prop >= 1) stop("`train_prop` must be in (0,1).", call. = FALSE)
  scores <- match.arg(scores)
  aggregation <- match.arg(aggregation)
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
  if (any(is_interaction)) stop("Interaction terms (colon terms such as x1:x2, or latent XWITH terms) are not supported: they cannot be reconstructed in the held-out test set, so the out-of-sample prediction would silently omit them. Refit the model without interaction terms.", call. = FALSE)
  struct <- pt[is_struct & !is_interaction, c("lhs", "rhs")]
  outcomes_all <- unique(struct$lhs)
  if (!length(outcomes_all)) stop("No structural regressions found.", call. = FALSE)
  
  latents_all <- tryCatch(lavaan::lavNames(fit, type = "lv"), error = function(e) character(0))
  outcome_is_latent <- outcomes_all %in% latents_all
  any_latents <- length(latents_all) > 0L
  
  estimator_use <- OPT$estimator
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
    score_method <- if (identical(scores, "bartlett")) "Bartlett" else "regression"
    tryCatch(
      lavaan::lavPredict(fit_obj, newdata = newdata, type = "lv", method = score_method),
      error = function(e) NULL
    )
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
  
  # in-sample R^2 di baseline, sulla stessa metrica della cross-validazione
  coefs_full <- get_struct_coefs(fit)
  fs_full <- if (any_latents) compute_lv_scores(fit, data_in) else NULL
  baseline_r2 <- stats::setNames(rep(NA_real_, length(outcomes_all)), outcomes_all)
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
  
  set.seed(seed)
  N <- nrow(data_in); n_train <- max(3L, floor(train_prop * N))
  maxN <- if (is.numeric(times)) as.integer(times) else as.integer(max_times)
  res_mat <- matrix(NA_real_, nrow = maxN, ncol = length(outcomes_all))
  colnames(res_mat) <- outcomes_all
  used <- 0L
  attempted <- 0L
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
  # cap di sicurezza per casi heywood e altri
  attempt_cap <- max(as.integer(max_times), 20L * target_splits)
  
  while (used < target_splits && attempted < attempt_cap) {
    attempted <- attempted + 1L
    tr_idx <- sort(sample.int(N, size = n_train, replace = FALSE))
    te_idx <- setdiff(seq_len(N), tr_idx); if (length(te_idx) < 3L) next
    dat_train <- data_in[tr_idx, , drop = FALSE]
    dat_test <- data_in[te_idx,  , drop = FALSE]
    
    fit_train <- tryCatch(
      lavaan::sem(
        model = model_obj, data = dat_train,
        estimator = estimator_use, missing = missing_opt,
        se = "none", 
        test = "none",
        std.lv = std_lv, parameterization = parameterization,
        ordered = if (length(ordered_vars)) ordered_vars else NULL,
        warn = !quiet
      ),
      error = function(e) NULL
    )
    
    # tiene uno split solo se la soluzione al training-set converge E SE e' ammissibile
    if (is.null(fit_train)) next
    conv_ok <- isTRUE(tryCatch(lavaan::lavInspect(fit_train, "converged"), error = function(e) FALSE))
    if (!conv_ok) next
    admissible_ok <- isTRUE(tryCatch(lavaan::lavInspect(fit_train, "post.check"), error = function(e) FALSE))
    if (!admissible_ok) next
    
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
  
  if (attempted >= attempt_cap && used < target_splits) {
    warning(sprintf(
      "Only %d usable splits obtained after %d attempts (requested %d). Many training-set solutions failed to converge or were inadmissible.",
      used, attempted, target_splits), call. = FALSE)
  }
  
  if (used < nrow(res_mat)) res_mat <- res_mat[seq_len(used), , drop = FALSE]
  
  # sommari split-wise
  r2_cv_mean <- apply(res_mat, 2L, function(v) { v <- v[is.finite(v)]; if (!length(v)) NA_real_ else mean(v) })
  r2_cv_median <- apply(res_mat, 2L, function(v) { v <- v[is.finite(v)]; if (!length(v)) NA_real_ else stats::median(v) })
  r2_cv_sd <- apply(res_mat, 2L, function(v) { v <- v[is.finite(v)]; if (!length(v)) NA_real_ else stats::sd(v) })
  n_ok <- apply(res_mat, 2L, function(v) sum(is.finite(v)))
  
  base_num <- as.numeric(baseline_r2)
  mean_num <- as.numeric(r2_cv_mean)
  median_num <- as.numeric(r2_cv_median)
  
  # indice di overfitting assoluto Delta (cioe' in-sample - CV)
  drop_mean_abs <- base_num - mean_num
  drop_med_abs  <- base_num - median_num
  
  # indice di overfitting proporzionale Delta% (cioe' 100 * (in-sample - CV) / in-sample)
  ok_pct <- is.finite(base_num) & (base_num >= min_r2_for_pct)
  drop_mean_pct   <- ifelse(ok_pct, 100 * (base_num - mean_num)   / base_num, NA_real_)
  drop_median_pct <- ifelse(ok_pct, 100 * (base_num - median_num) / base_num, NA_real_)
  
  out_type <- ifelse(outcome_is_latent, "latent", "observed")
  res_table <- data.frame(
    outcome = outcomes_all,
    type = out_type,
    r2_in = base_num,
    r2_cv_mean = mean_num,
    r2_cv_median = median_num,
    r2_cv_sd = as.numeric(r2_cv_sd),
    drop_mean_abs = as.numeric(drop_mean_abs),
    drop_med_abs = as.numeric(drop_med_abs),
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
    splits_attempted = attempted,
    tol = tol,
    window = window,
    min_r2_for_pct = min_r2_for_pct,
    scores = scores,
    aggregation = aggregation
  )
  class(res) <- c("lav_cv", "list")
  
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
      graphics::abline(h = utils::tail(cm, 1L), lty = 3)
    }
  }
  
  res
}

### S3 methods print and summary:

#' @param x A 'lav_cv' object.
#' @param ... Additional arguments; unused.
#' @rdname lav_cv
#' @export
print.lav_cv <- function(x, digits = x$digits %||% 3L, ...) {
  stopifnot(is.list(x), !is.null(x$table))
  u2  <- if (isTRUE(l10n_info()[["UTF-8"]])) "\u00B2" else "^2"
  dlt <- if (isTRUE(l10n_info()[["UTF-8"]])) "\u0394" else "Delta"
  
  auto_txt <- if (grepl("^auto\\(", as.character(x$times))) {
    paste0("auto-stop at ", gsub("[^0-9]", "", x$times), " splits")
  } else {
    paste0(x$times, " splits")
  }
  cat(sprintf("\nRepeated holdout (Monte Carlo) cross-validation of R%s\n", u2))
  cat(sprintf("Splits: %s (train proportion = %.2f; N = %d; seed = %s)\n\n",
              auto_txt, x$train_prop, x$N, as.character(x$seed)))
  
  tbl <- x$table
  aggregation <- x$aggregation %||% "median"
  if (!aggregation %in% c("median", "mean")) aggregation <- "median"
  scores <- x$scores %||% "regression"
  
  r2_cv     <- if (identical(aggregation, "mean")) tbl$r2_cv_mean    else tbl$r2_cv_median
  delta_abs <- if (identical(aggregation, "mean")) tbl$drop_mean_abs else tbl$drop_med_abs
  delta_pct <- if (identical(aggregation, "mean")) tbl$drop_mean_pct else tbl$drop_med_pct
  
  dash <- "-"
  fmt  <- function(z, d = digits) ifelse(is.na(z), dash, formatC(z, digits = d, format = "f"))
  fmt1 <- function(z) ifelse(is.na(z), dash, formatC(z, digits = 1, format = "f"))
  
  cat(sprintf("Factor-score method: %s; CV aggregation: %s.\n\n", scores, aggregation))
  
  fmtstr <- "%-16s %11s %11s %9s %14s %9s\n"
  nrule  <- 75L
  for (tp in c("observed", "latent")) {
    idx <- tbl$type == tp
    sub <- tbl[idx, , drop = FALSE]
    if (!nrow(sub)) next
    
    out_names <- sub$outcome
    r2_infit  <- round(sub$r2_in, digits)
    r2_agg    <- round(r2_cv[idx], digits)
    sd_cv     <- round(sub$r2_cv_sd, digits)
    dabs_tp   <- round(delta_abs[idx], digits)
    dpct_tp   <- ifelse(is.finite(delta_pct[idx]), round(delta_pct[idx], 1), NA)
    
    cat(sprintf("%s outcomes\n", tools::toTitleCase(tp)))
    cat(rep("-", nrule), "\n", sep = "")
    cat(sprintf(fmtstr,
                "Outcome",
                paste0("R", u2, " in"),
                paste0("R", u2, " CV"),
                "CV SD",
                paste0(dlt, " (in-CV)"),
                paste0(dlt, "%")))
    cat(rep("-", nrule), "\n", sep = "")
    for (i in seq_len(nrow(sub))) {
      cat(sprintf(fmtstr,
                  out_names[i],
                  fmt(r2_infit[i]),
                  fmt(r2_agg[i]),
                  fmt(sd_cv[i]),
                  fmt(dabs_tp[i]),
                  fmt1(dpct_tp[i])))
    }
    cat("\n")
  }
  
  # footnotes
  cat(sprintf("%s%% is the proportional overfitting index and the main diagnostic:\n", dlt))
  cat(sprintf("  %s%% = 100 * (in-sample R%s - CV R%s) / in-sample R%s.\n", dlt, u2, u2, u2))
  cat(sprintf("%s (in-CV) is the same gap in raw R%s units.\n", dlt, u2))
  cat(sprintf("%s%% near 0 indicates the structural relation generalizes; %s%% > 100%% is\n", dlt, dlt))
  cat(sprintf("equivalent to a negative CV R%s: the in-sample R%s does not generalize.\n", u2, u2))
  cat(sprintf("%s%% is suppressed when in-sample R%s < %.2f (near-zero baseline).\n",
              dlt, u2, x$min_r2_for_pct %||% 0.05))
  cat(sprintf("CV SD is the split-wise standard deviation of R%s; aggregation across splits: %s.\n", u2, aggregation))
  cat("Splits retained only if the training solution converged and was admissible.\n")
  cat(sprintf("Auto-stop rule: window = %d, tol = %.3f (applied when times = 'auto').\n", x$window, x$tol))
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
