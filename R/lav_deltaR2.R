#' Incremental effect sizes for structural predictors in fitted \code{lavaan} models
#'
#' Compute outcome-specific incremental effect sizes for structural predictors in
#' a fitted \code{lavaan} model by comparing the fitted model with nested reduced
#' models in which one predictor, or a block of predictors, is fixed to zero.
#' For each tested predictor set, the function returns the reduction in explained
#' variance (\eqn{\Delta R^2}), the corresponding part \eqn{R^2}, and Cohen's
#' \eqn{f^2 = \Delta R^2 / (1 - R^2_{full})}.
#'
#' Reduced models preserve all untargeted model parameters and differ from the
#' fitted model only in that the tested structural regression path(s) are fixed
#' to zero. Accordingly, the resulting \eqn{\Delta R^2} quantifies the unique
#' contribution of the tested predictor or predictor block to the explained
#' variance of a given endogenous variable, conditional on the rest of the model.
#'
#' If \code{data} is supplied, reduced models are refit from that data set.
#' Otherwise, the function reuses the internal data and fitting options stored in
#' \code{fit}. The function is intended for converged \code{lavaan} models with
#' free structural regressions.
#'
#' @usage
#' lav_deltaR2(
#'   fit,
#'   data = NULL,
#'   outcome = NULL,
#'   terms = NULL,
#'   block = NULL,
#'   quiet = TRUE,
#'   digits = 3L
#' )
#'
#' @param fit A fitted \code{lavaan} object.
#' @param data Optional data frame used to fit \code{fit}. If supplied, this data
#'   frame is used to refit reduced models. If \code{NULL}, reduced models are
#'   refit using the internal data and options stored in \code{fit}.
#' @param outcome Optional character vector naming endogenous variables for which
#'   incremental effect sizes should be computed. If \code{NULL}, all endogenous
#'   variables with at least one free structural predictor are used.
#' @param terms Optional character vector naming predictors to test one at a
#'   time. Ignored if \code{block} is supplied. If \code{NULL}, all eligible
#'   structural predictors are tested individually.
#' @param block Optional character vector naming predictors to remove jointly.
#'   When supplied, one reduced-model comparison is computed per selected outcome
#'   for the entire predictor block. A block is evaluated only for outcomes in
#'   which all named predictors appear as free structural regressors.
#' @param quiet Logical. If \code{TRUE}, suppress reduced-model refit messages.
#' @param digits Non-negative integer giving the default number of digits used in
#'   printing.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{delta_r2_table}: A data frame containing one row per
#'   outcome-by-test comparison, with columns:
#'   \itemize{
#'     \item \code{outcome}: endogenous variable whose \eqn{R^2} is evaluated;
#'     \item \code{tested_set}: predictor or predictor block fixed to zero in the
#'       reduced model;
#'     \item \code{test_type}: either \code{"single"} or \code{"block"};
#'     \item \code{n_terms}: number of predictors removed jointly;
#'     \item \code{group}: numeric group index;
#'     \item \code{group_label}: group label, if available;
#'     \item \code{r2_full}: \eqn{R^2} from the fitted model;
#'     \item \code{r2_reduced}: \eqn{R^2} from the reduced model;
#'     \item \code{delta_r2}: \eqn{R^2_{full} - R^2_{reduced}};
#'     \item \code{part_r2}: equal to \code{delta_r2} in this implementation;
#'     \item \code{f2}: Cohen's \eqn{f^2};
#'     \item \code{f2_magnitude}: conventional descriptive interpretation of
#'       \code{f2} based on Cohen's benchmarks (\code{"negligible"},
#'       \code{"small"}, \code{"medium"}, or \code{"large"});
#'     \item \code{converged}: logical indicator of whether the reduced model
#'       converged.
#'   }
#'   \item \code{settings}: list of user-supplied function settings;
#'   \item \code{group_var}: grouping variable name, or \code{NULL};
#'   \item \code{group_labels}: group labels, if available;
#'   \item \code{call}: matched function call;
#'   \item \code{digits}: default number of digits used for printing.
#' }
#'
#' @section Notes:
#' For a given endogenous variable, \eqn{\Delta R^2} is the reduction in
#' explained variance obtained when the tested predictor or predictor block is
#' fixed to zero, and the reduced model is refit with all untargeted parameters
#' left free.
#' In this implementation, part \eqn{R^2} is numerically identical to \eqn{\Delta R^2}. 
#' Cohen's \eqn{f^2} expresses the same quantity relative to the residual variance of the
#' full model. The \code{f2_magnitude} column applies Cohen's conventional 
#' benchmarks to \eqn{f^2} and is intended only as a descriptive aid.
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences}
#' (2nd ed.). Lawrence Erlbaum Associates.
#'
#' Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
#' \emph{Journal of Statistical Software}, \emph{48}(2), 1-36.
#' \doi{10.18637/jss.v048.i02}
#'
#' @seealso \code{\link{lav_cv}} for cross-validated \eqn{R^2}, \code{\link{lav_vif}}
#'   for variance inflation factors.
#'
#' @examples
#' library("lavaan")
#' data("PoliticalDemocracy", 
#' package = "lavaan")
#' model <- "
#' ind60 =~ x1 + x2 + x3
#' dem60 =~ y1 + y2 + y3 + y4
#' dem65 =~ y5 + y6 + y7 + y8
#' dem60 ~ ind60
#' dem65 ~ ind60 + dem60
#' y1 ~~ y5
#' y2 ~~ y4 + y6
#' y3 ~~ y7
#' y4 ~~ y8
#' y6 ~~ y8
#' "
#' fit <- lavaan::sem(model, 
#' data = PoliticalDemocracy)
#' lav_deltaR2(fit = fit)
#' lav_deltaR2(fit = fit, outcome = "dem65")
#' lav_deltaR2(fit = fit, outcome = "dem65", terms = "ind60")
#' lav_deltaR2(fit = fit, outcome = "dem65", block = c("ind60", "dem60"))
#'
#' @importFrom lavaan lavInspect lavNames parTable
#' @importFrom rlang %||%
#' @export
lav_deltaR2 <- function(
    fit,
    data = NULL,
    outcome = NULL,
    terms = NULL,
    block = NULL,
    quiet = TRUE,
    digits = 3L
) {
  
  if (!inherits(fit, "lavaan")) {
    stop("'fit' must be a fitted 'lavaan' object.", call. = FALSE)
  }
  
  fit_converged <- tryCatch(lavaan::lavInspect(fit, "converged"),
                            error = function(e) FALSE)
  if (!isTRUE(fit_converged)) {
    stop("'fit' did not converge.", call. = FALSE)
  }
  
  if (!is.null(outcome) && !is.character(outcome)) {
    stop("'outcome' must be NULL or a character vector.", call. = FALSE)
  }
  
  if (!is.null(terms) && !is.character(terms)) {
    stop("'terms' must be NULL or a character vector.", call. = FALSE)
  }
  
  if (!is.null(block) && !is.character(block)) {
    stop("'block' must be NULL or a character vector.", call. = FALSE)
  }
  
  if (!is.null(block) && !is.null(terms)) {
    stop("Supply either 'terms' or 'block', not both.", call. = FALSE)
  }
  
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("'quiet' must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0) {
    stop("'digits' must be a non-negative integer.", call. = FALSE)
  }
  
  digits <- as.integer(digits)
  
  data_in <- NULL
  if (!is.null(data)) {
    data_in <- as.data.frame(data)
  }
  
  OPT <- fit@Options
  estimator <- OPT$estimator
  missing_opt <- OPT$missing
  
  as_one_string_or_null <- function(x) {
    if (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)) {
      return(x)
    }
    NULL
  }
  
  group_var <- as_one_string_or_null(OPT$group)
  cluster_var <- as_one_string_or_null(OPT$cluster)
  sampling_weights <- as_one_string_or_null(OPT$sampling.weights)
  
  group_labels <- tryCatch(lavaan::lavInspect(fit, "group.label"),
                           error = function(e) NULL)
  
  ordered_vars <- tryCatch(lavaan::lavNames(fit, type = "ov.ord"),
                           error = function(e) character(0))
  
  # trova regr. strutturali
  pt <- lavaan::parTable(fit)
  
  is_struct <- pt$op == "~" & pt$rhs != "1" & pt$free != 0L
  struct <- pt[is_struct, c("lhs", "rhs", "group"), drop = FALSE]
  
  if (!nrow(struct)) {
    stop("No free structural regressions found.", call. = FALSE)
  }
  
  struct <- unique(struct)
  outcomes_all <- unique(struct$lhs)
  
  if (is.null(outcome)) {
    outcome_use <- outcomes_all
  } else {
    outcome_use <- unique(outcome[outcome %in% outcomes_all])
    if (!length(outcome_use)) {
      stop("None of the requested 'outcome' values are endogenous variables ",
           "with free structural predictors.", call. = FALSE)
    }
  }
  
  struct <- struct[struct$lhs %in% outcome_use, , drop = FALSE]
  
  if (!is.null(block)) {
    block <- unique(block)
    tests <- list(block)
    test_labels <- paste(block, collapse = " + ")
    test_type <- "block"
  } else {
    terms_all <- unique(struct$rhs)
    
    if (is.null(terms)) {
      term_use <- terms_all
    } else {
      term_use <- unique(terms[terms %in% terms_all])
      if (!length(term_use)) {
        stop("None of the requested 'terms' appear as free structural ",
             "predictors in the selected outcomes.", call. = FALSE)
      }
    }
    
    tests <- as.list(term_use)
    test_labels <- term_use
    test_type <- "single"
  }
  
  r2_to_df <- function(x) {
    if (is.list(x)) {
      rows <- list()
      for (g in seq_along(x)) {
        v <- x[[g]]
        if (is.null(v) || !length(v)) next
        rows[[length(rows) + 1L]] <- data.frame(
          outcome = names(v),
          group   = g,
          r2      = as.numeric(v),
          stringsAsFactors = FALSE
        )
      }
      if (length(rows)) {
        out <- do.call(rbind, rows)
      } else {
        out <- data.frame(outcome = character(0), group = integer(0),
                          r2 = numeric(0), stringsAsFactors = FALSE)
      }
    } else {
      out <- data.frame(outcome = names(x), group = 1L,
                        r2 = as.numeric(x), stringsAsFactors = FALSE)
    }
    rownames(out) <- NULL
    out
  }
  
  r2_full <- r2_to_df(lavaan::lavInspect(fit, "r2"))
  r2_full <- r2_full[r2_full$outcome %in% outcome_use, , drop = FALSE]
  
  out_rows <- list()
  idx_row  <- 0L
  
  for (y in outcome_use) {
    avail <- unique(struct$rhs[struct$lhs == y])
    
    for (j in seq_along(tests)) {
      remove_terms <- unique(tests[[j]])
      
      if (test_type == "block") {
        if (!all(remove_terms %in% avail)) next
      } else {
        if (!(remove_terms[1L] %in% avail)) next
      }
      
      idx_fix <- which(
        pt$op == "~" &
          pt$lhs == y &
          pt$rhs %in% remove_terms &
          pt$rhs != "1" &
          pt$free != 0L
      )
      
      if (!length(idx_fix)) next
      
      pt_red <- pt
      pt_red$free[idx_fix] <- 0L
      if ("est"    %in% names(pt_red)) pt_red$est[idx_fix]    <- 0
      if ("ustart" %in% names(pt_red)) pt_red$ustart[idx_fix] <- 0
      if ("start"  %in% names(pt_red)) pt_red$start[idx_fix]  <- 0
      if ("label"  %in% names(pt_red)) pt_red$label[idx_fix]  <- ""
      
      pt_red <- pt_red[pt_red$op != ":=", , drop = FALSE]
      
      refit_error <- NULL
      red_fit <- NULL
      
      if (is.null(data_in)) {
        refit_opt <- fit@Options
        refit_opt$se <- "none"
        refit_opt$test <- "none"
        
        red_fit <- tryCatch(
          lavaan::lavaan(
            model       = pt_red,
            slotData    = fit@Data,
            slotOptions = refit_opt,
            warn        = !quiet
          ),
          error = function(e) {
            refit_error <<- conditionMessage(e)
            NULL
          }
        )
      } else {
        refit_args <- list(
          model = pt_red,
          data = data_in,
          group = group_var,
          cluster = cluster_var,
          sampling.weights = sampling_weights,
          estimator = estimator,
          missing = missing_opt,
          fixed.x = OPT$fixed.x,
          conditional.x = OPT$conditional.x,
          meanstructure = OPT$meanstructure,
          ordered = if (length(ordered_vars)) ordered_vars else NULL,
          se = "none",
          test = "none",
          warn = !quiet
        )
        
        keep_arg <- !vapply(refit_args, is.null, logical(1L))
        refit_args <- refit_args[keep_arg]
        
        red_fit <- tryCatch(
          do.call(lavaan::lavaan, refit_args),
          error = function(e) {
            refit_error <<- conditionMessage(e)
            NULL
          }
        )
      }
      
      conv <- FALSE
      if (!is.null(red_fit)) {
        conv <- isTRUE(tryCatch(lavaan::lavInspect(red_fit, "converged"),
                                error = function(e) FALSE))
      }
      
      if (!conv && !quiet && !is.null(refit_error)) {
        message("lav_deltaR2: reduced model for '", y, " ~ ",
                paste(remove_terms, collapse = " + "),
                "' failed: ", refit_error)
      }
      
      if (conv) {
        r2_red <- r2_to_df(lavaan::lavInspect(red_fit, "r2"))
        r2_red <- r2_red[r2_red$outcome == y, , drop = FALSE]
      } else {
        r2_red <- data.frame(outcome = character(0), group = integer(0),
                             r2 = numeric(0), stringsAsFactors = FALSE)
      }
      
      full_sub   <- r2_full[r2_full$outcome == y, , drop = FALSE]
      all_groups <- sort(unique(full_sub$group))
      if (!length(all_groups)) all_groups <- 1L
      
      for (g in all_groups) {
        pos_full <- match(g, full_sub$group)
        pos_red  <- match(g, r2_red$group)
        
        r2f <- if (is.na(pos_full)) NA_real_ else full_sub$r2[pos_full]
        r2r <- if (is.na(pos_red))  NA_real_ else r2_red$r2[pos_red]
        if (!conv) r2r <- NA_real_
        
        delta <- if (!is.na(r2f) && !is.na(r2r)) r2f - r2r else NA_real_
        tol <- 0.5 * 10^(-digits)
        if (is.finite(delta) && abs(delta) < tol) {
          delta <- 0
        }
        
        denom <- if (!is.na(r2f)) 1 - r2f else NA_real_
        f2 <- if (!is.na(delta) && !is.na(denom) && denom > 0) {
          delta / denom
        } else { NA_real_ }
        
        f2_magnitude <- if (is.na(f2)) {
          NA_character_
        } else if (f2 < 0.02) {
          "negligible"
        } else if (f2 < 0.15) {
          "small"
        } else if (f2 < 0.35) {
          "medium"
        } else {
          "large"
        }
        
        idx_row <- idx_row + 1L
        out_rows[[idx_row]] <- data.frame(
          outcome = y,
          tested_set = if (length(test_labels) == 1L && test_type == "block") {
            test_labels
          } else { test_labels[j] },
          test_type = if (length(remove_terms) == 1L) "single" else "block",
          n_terms = length(remove_terms),
          group = g,
          group_label = if (!is.null(group_labels) && length(group_labels) >= g) {
            as.character(group_labels[g])
          } else { as.character(g) },
          r2_full = r2f,
          r2_reduced = r2r,
          delta_r2 = delta,
          part_r2 = delta,
          f2 = f2,
          f2_magnitude = f2_magnitude,
          converged = conv,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (!length(out_rows)) {
    stop("No eligible reduced-model comparisons could be formed from the ",
         "requested arguments.", call. = FALSE)
  }
  
  delta_r2_table <- do.call(rbind, out_rows)
  rownames(delta_r2_table) <- NULL
  
  delta_r2_table <- delta_r2_table[
    order(delta_r2_table$group, delta_r2_table$outcome,
          delta_r2_table$tested_set), , drop = FALSE
  ]
  
  res <- list(
    delta_r2_table = delta_r2_table,
    settings = list(
      outcome = if (is.null(outcome)) NULL else outcome,
      terms = if (is.null(block)) terms else NULL,
      block = if (is.null(block)) NULL else block,
      quiet = quiet
    ),
    group_var = group_var,
    group_labels = if (is.null(group_labels)) NULL else group_labels,
    call = match.call(),
    digits = digits
  )
  
  class(res) <- "lav_deltaR2"
  return(res)
}


# --- S3 methods: print and summary ------------------------------------------

#' @param x A 'lav_deltaR2' object.
#' @param ... Additional arguments; unused.
#' @rdname lav_deltaR2
#' @export
print.lav_deltaR2 <- function(x, ...) {
  tbl <- x$delta_r2_table
  digits <- x$digits %||% 3L
  u2 <- if (isTRUE(l10n_info()[["UTF-8"]])) "\u00B2" else "^2"
  
  if (!nrow(tbl)) {
    cat("No incremental effect sizes in delta_r2_table.\n")
    return(invisible(x))
  }
  
  fmt_num <- function(z, k = digits) {
    if (is.na(z)) "NA" else formatC(z, digits = k, format = "f")
  }
  
  outcome_txt <- as.character(tbl$outcome)
  tested_txt <- as.character(tbl$tested_set)
  group_txt <- if (!is.null(x$group_var)) {
    as.character(tbl$group_label)
  } else {
    as.character(tbl$group)
  }
  r2_full_txt <- vapply(tbl$r2_full, fmt_num, character(1L), k = digits)
  r2_red_txt <- vapply(tbl$r2_reduced, fmt_num, character(1L), k = digits)
  part_r2_txt <- vapply(tbl$part_r2, fmt_num, character(1L), k = digits)
  f2_txt <- vapply(tbl$f2, fmt_num, character(1L), k = digits)
  
  if ("f2_magnitude" %in% names(tbl)) {
    mag_txt <- ifelse(is.na(tbl$f2_magnitude), "NA", as.character(tbl$f2_magnitude))
  } else {
    mag_txt <- rep("", nrow(tbl))
  }
  
  hdr_outcome <- "Outcome"
  hdr_tested <- "Tested"
  hdr_group <- "Group"
  hdr_r2_full <- paste0("R", u2, " full")
  hdr_r2_red <- paste0("R", u2, " reduced")
  hdr_part_r2 <- paste0("part R", u2)
  hdr_f2 <- paste0("f", u2)
  hdr_mag <- "Magnitude"
  
  w_outcome <- max(nchar(c(hdr_outcome, outcome_txt), type = "width"))
  w_tested <- max(nchar(c(hdr_tested, tested_txt), type = "width"))
  w_group <- max(nchar(c(hdr_group, group_txt), type = "width"))
  w_r2_full <- max(nchar(c(hdr_r2_full, r2_full_txt), type = "width"))
  w_r2_red <- max(nchar(c(hdr_r2_red, r2_red_txt), type = "width"))
  w_part_r2 <- max(nchar(c(hdr_part_r2, part_r2_txt), type = "width"))
  w_f2 <- max(nchar(c(hdr_f2, f2_txt), type = "width"))
  w_mag <- max(nchar(c(hdr_mag, mag_txt), type = "width"))
  
  row_fmt <- paste0(
    "%-", w_outcome, "s  ",
    "%-", w_tested, "s  ",
    "%-", w_group, "s  ",
    "%",  w_r2_full, "s  ",
    "%",  w_r2_red, "s  ",
    "%",  w_part_r2, "s  ",
    "%",  w_f2, "s  ",
    "%-", w_mag, "s\n"
  )
  
  header_line <- sprintf(
    row_fmt,
    hdr_outcome,
    hdr_tested,
    hdr_group,
    hdr_r2_full,
    hdr_r2_red,
    hdr_part_r2,
    hdr_f2,
    hdr_mag
  )
  
  line_width <- nchar(sub("\n$", "", header_line), type = "width")
  rule <- paste(rep("-", line_width), collapse = "")
  
  cat("\nIncremental effect sizes from reduced-model comparisons\n")
  cat(rule, "\n", sep = "")
  cat(header_line)
  cat(rule, "\n", sep = "")
  
  for (i in seq_len(nrow(tbl))) {
    cat(sprintf(
      row_fmt,
      outcome_txt[i],
      tested_txt[i],
      group_txt[i],
      r2_full_txt[i],
      r2_red_txt[i],
      part_r2_txt[i],
      f2_txt[i],
      mag_txt[i]
    ))
  }
  
  cat(rule, "\n", sep = "")
  delta_chr <- if (isTRUE(l10n_info()[["UTF-8"]])) "\u0394" else "Delta"
  
  cat(
    paste0(
      "part R", u2, " equals ", delta_chr, "R", u2,
      ". Cohen's f", u2, " = ", delta_chr, "R", u2,
      " / (1 - R", u2, " full).\n"
    )
  )
  
  cat("\n")
  invisible(x)
}

#' @param object A 'lav_deltaR2' object.
#' @param ... Additional arguments; unused.
#' @rdname lav_deltaR2
#' @export
summary.lav_deltaR2 <- function(object, ...) {
  print.lav_deltaR2(object, ...)
  invisible(object)
}
