#' False Discovery Rate (FDR) Correction for 'lavaan' parameter p-values
#'
#' Apply a false discovery rate correction (Benjamini-Yekutieli by default) to
#' the p-values of selected parameters from a fitted \code{lavaan} object.
#' 
#' Useful when a SEM includes many structural paths (or many other parameters 
#' of substantive interest) and there is the need to control the expected 
#' proportion of false positives among the parameters declared 'statistically 
#' significant'.
#'  
#' With many simultaneous tests, using \code{p < .05} for each parameter 
#' inflates the expected number of false positives (about \code{m * .05} under 
#' all true null hypotheses, where \code{m} is the number of tested parameters).
#' Benjamini-Yekutieli (BY) controls the False Discovery Rate (FDR) under 
#' arbitrary dependence structures, which is suitable for SEMs where structural 
#' paths are inherently dependent through shared latent variables, covariance 
#' matrices, and model constraints.
#'
#' @usage
#' lav_fdr(
#'   fit,
#'   ops = c("reg", "load", "var.cov"),
#'   family = c("by_group", "selected"),
#'   method = "BY",
#'   alpha = 0.05,
#'   standardized = c("std.all", "std.lv", "std.nox", "none")
#' )
#'
#' @param fit A fitted \code{lavaan} object.
#' @param ops Character. One of \code{"reg"} (regressions), \code{"load"} (factor loadings),
#'   or \code{"var.cov"} (variances/covariances/residual variances). Default is \code{"reg"}.
#' @param family Character. If \code{"selected"}, FDR is applied across all selected
#'   parameters jointly. If \code{"by_group"}, FDR is applied separately within each group.
#'   Default is \code{"by_group"}.
#' @param method Character method passed to \code{stats::p.adjust} (default \code{"BY"}).
#' @param alpha Numeric significance threshold for adjusted p-values (default 0.05).
#' @param standardized Which standardized column to include, or \code{"none"}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{fdr_table}: data.frame with raw and FDR-adjusted p-values.
#'   \item \code{settings}: list of settings used.
#'   \item \code{group_var}: group variable name (or \code{NULL}).
#'   \item \code{group_labels}: group labels if available.
#'   \item \code{call}: matched call.
#' }
#' The returned object has class \code{"lav_fdr"}.
#'
#' @examplesIf requireNamespace("lavaan", quietly = TRUE)
#' library("lavaan")
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
#' lav_fdr(fit = fit)
#'
#' @export
lav_fdr <- function(
    fit,
    ops = c("reg", "load", "var.cov"),
    family = c("by_group", "selected"),
    method = "BY",
    alpha = 0.05,
    standardized = c("std.all", "std.lv", "std.nox", "none")
) {
  if (!inherits(fit, "lavaan")) stop("`fit` must be a fitted lavaan object.")
  
  ops <- match.arg(ops)
  family <- match.arg(family)
  standardized <- match.arg(standardized)
  
  if (!is.character(method) || length(method) != 1L) stop("`method` must be a single character value.")
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).")
  }
  
  OPT <- fit@Options
  group_var <- OPT$group
  if (is.null(group_var) || !nzchar(group_var)) group_var <- NULL
  group_labels <- tryCatch(lavaan::lavInspect(fit, "group.label"), error = function(e) NULL)
  
  # estrai stime parametri
  pe <- lavaan::parameterEstimates(fit, standardized = TRUE)
  
  op_map <- list(
    reg = "~",
    load = "=~",
    var.cov = "~~"
  )
  op_target <- op_map[[ops]]

  sel <- rep(TRUE, nrow(pe))
  sel <- sel & (pe$op == op_target)
  
  # no intercette di default
  sel <- sel & !(pe$op == "~1" | (pe$op == "~" & pe$rhs == "1"))
  # no threshold
  sel <- sel & !(pe$op == "|")
  # no parametri definiti
  sel <- sel & !(pe$op == ":=")
  # no interazioni
  is_int <- (pe$op == "XWITH") |
    (!is.na(pe$rhs) & grepl(":", pe$rhs, fixed = TRUE)) |
    (!is.na(pe$lhs) & grepl(":", pe$lhs, fixed = TRUE))
  sel <- sel & !is_int
  
  pe2 <- pe[sel, , drop = FALSE]
  if (!nrow(pe2)) stop("No parameters matched the selection criteria for FDR correction.")
  
  # cerca che esistano le colonne
  if (!("group" %in% names(pe2))) pe2$group <- 1L
  if (all(is.na(pe2$group))) pe2$group <- 1L
  
  # trova i pvalue che vanno corretti
  ok_p <- !is.na(pe2$pvalue) & is.finite(pe2$pvalue)
  p_adj <- rep(NA_real_, nrow(pe2))
  
  if (family == "by_group" && length(unique(pe2$group)) > 1L) {
    groups <- sort(unique(pe2$group))
    for (g in groups) {
      idx <- which(pe2$group == g & ok_p)
      if (length(idx)) p_adj[idx] <- stats::p.adjust(pe2$pvalue[idx], method = method)
    }
  } else {
    idx <- which(ok_p)
    if (length(idx)) p_adj[idx] <- stats::p.adjust(pe2$pvalue[idx], method = method)
  }
  
  # selezione della colonna standardizzata
  std_col <- NA_character_
  std_val <- rep(NA_real_, nrow(pe2))
  if (standardized != "none") {
    std_col <- standardized
    if (std_col %in% names(pe2)) std_val <- pe2[[std_col]]
  }
  
  fdr_table <- data.frame(
    lhs = pe2$lhs,
    op = pe2$op,
    rhs = pe2$rhs,
    group = pe2$group,
    est = pe2$est,
    se = pe2$se,
    z = pe2$z,
    p = pe2$pvalue,
    p_fdr = p_adj,
    sig_fdr = (!is.na(p_adj)) & (p_adj < alpha),
    stringsAsFactors = FALSE
  )
  
  if (standardized != "none") {
    fdr_table[[std_col]] <- std_val
  }
  
  fdr_table <- fdr_table[order(fdr_table$group, fdr_table$lhs, fdr_table$op, fdr_table$rhs), , drop = FALSE]
  rownames(fdr_table) <- NULL
  
  res <- list(
    fdr_table = fdr_table,
    settings = list(
      ops = ops,
      family = family,
      method = method,
      alpha = alpha,
      standardized = standardized
    ),
    group_var = group_var,
    group_labels = group_labels %||% NULL,
    call = match.call()
  )
  class(res) <- c("lav_fdr", "list")
  res
}

### metodi S3 print e summary:

#' @param x A 'lav_fdr' object.
#' @param ... Additional arguments; unused.
#' @rdname lav_fdr
#' @export
print.lav_fdr <- function(x, ...) {
  tbl <- x$fdr_table
  gl  <- x$group_labels
  st  <- x$settings
  
  if (!nrow(tbl)) {
    cat("No parameters in fdr_table.\n")
    return(invisible(x))
  }
  
  fmt  <- function(v, k = 3) if (is.na(v)) "NA" else formatC(v, format = "f", digits = k)
  fmtp <- function(p) ifelse(is.na(p), "NA",
                             ifelse(p < .001, "< .001", formatC(p, digits = 3, format = "f")))
  zfmt <- function(z) if (is.na(z)) "NA" else formatC(z, format = "f", digits = 4)
  
  std_name <- intersect(names(tbl), c("std.all", "std.lv", "std.nox"))
  std_name <- if (length(std_name)) std_name[1] else NULL
  
  groups <- sort(unique(tbl$group))
  for (g in groups) {
    gname <- if (length(gl) && !is.na(g) && g <= length(gl)) gl[g] else paste0("Group ", g)
    sub <- tbl[tbl$group == g, , drop = FALSE]
    par_texts <- paste(sub$lhs, sub$op, sub$rhs)
    par_width <- max(nchar(par_texts) + 2, 9)
    
    cat("\nFDR-corrected parameters\n")
    cat(sprintf("Ops: %s | Family: %s | Method: %s | alpha: %.3f | %s\n",
                st$ops, st$family, st$method, st$alpha, gname))
    
    dash_n <- (if (is.null(std_name)) 76 else 88) + 12 + (par_width - 30)
    
    cat(rep("-", dash_n), "\n", sep = "")
    if (is.null(std_name)) {
      cat(sprintf("%-*s %10s %8s %8s %10s %10s\n", par_width,
                  "Parameter", "est", "SE", "z", "p", "p_FDR"))
    } else {
      cat(sprintf("%-*s %10s %8s %8s %10s %10s %10s\n", par_width,
                  "Parameter", "est", "SE", "z", "p", "p_FDR", std_name))
    }
    cat(rep("-", dash_n), "\n", sep = "")
    
    for (i in seq_len(nrow(sub))) {
      par_txt <- paste(sub$lhs[i], sub$op[i], sub$rhs[i])
      flag <- if (isTRUE(sub$sig_fdr[i])) "*" else ""
      if (is.null(std_name)) {
        cat(sprintf("%-*s %10s %8s %8s %10s %10s %s\n", par_width,
                    par_txt,
                    fmt(sub$est[i]),
                    fmt(sub$se[i]),
                    zfmt(sub$z[i]),
                    fmtp(sub$p[i]),
                    fmtp(sub$p_fdr[i]),
                    flag))
      } else {
        cat(sprintf("%-*s %10s %8s %8s %10s %10s %10s %s\n", par_width,
                    par_txt,
                    fmt(sub$est[i]),
                    fmt(sub$se[i]),
                    zfmt(sub$z[i]),
                    fmtp(sub$p[i]),
                    fmtp(sub$p_fdr[i]),
                    fmt(sub[[std_name]][i]),
                    flag))
      }
    }
    cat(rep("-", dash_n), "\n", sep = "")
    cat(sprintf("FDR-significant (p_FDR < alpha): %d / %d\n",
                sum(sub$sig_fdr, na.rm = TRUE), nrow(sub)))
  }
  
  cat("\nNote. '*' marks parameters with FDR-adjusted p-value < alpha.\n")
  invisible(x)
}

#' @param object A 'lav_fdr' object.
#' @param ... Passed to \code{print.lav_fdr()}.
#' @rdname lav_fdr
#' @export
summary.lav_fdr <- function(object, ...) {
  print.lav_fdr(object, ...)
  invisible(object)
}