#' Post-Estimation Utilities for 'lavaan' Fitted Models
#' 
#' @description
#' \if{html}{\figure{logo.png}{options: style="float: right; margin: 0 0 8px 12px;" alt="lavinteract logo" width="110"}}
#' 
#' Post-estimation tools for structural equation models fitted with 'lavaan'.
#' Provides methods for probing observed and latent interactions, diagnosing local
#' misfit and multicollinearity, quantifying incremental effect sizes for 
#' structural predictors, assessing predictive performance by repeated holdout 
#' cross-validation, and adjusting selected parameter p-values for multiple 
#' testing. Functions operate from a fitted model object and, when needed, refit 
#' auxiliary or reduced models while preserving the original SEM specification.
#' 
#' @details 
#' The functions are:
#' \itemize{
#'   \item \code{\link{lav_slopes}}: simple slopes and interaction plots from a fitted 'lavaan' model.
#'   \item \code{\link{lav_jn}}: Johnson-Neyman regions of significance for continuous moderators in a fitted 'lavaan' model.
#'   \item \code{\link{lav_deltaR2}}: incremental effect sizes (part \eqn{R^2} and Cohen's \eqn{f^2}) for structural predictors via reduced-model comparisons.
#'   \item \code{\link{lav_localfit}}: residual-based local fit diagnostics and heatmaps for fitted 'lavaan' models.
#'   \item \code{\link{lav_vif}}: variance inflation factors for structural predictors with measurement preserved.
#'   \item \code{\link{lav_cv}}: repeated holdout (Monte Carlo) cross-validation of \eqn{R^2} for SEM outcomes.
#'   \item \code{\link{lav_fdr}}: false discovery rate correction for selected 'lavaan' parameter p-values.
#' }
#'
#' @section Note:
#' The development of this package grew from ongoing discussions and interactions (sic)
#' with colleagues, in particular Dr. Cataldo Giuliano Gemmano, whose steady 
#' feedback and support helped shape it. 
#'
#' @author
#' Giuseppe Corbelli (<giuseppe.corbelli@uninettunouniversity.net>)
#'
#' @docType package
#' @name lavinteract
#' @aliases lavinteract-package lavinteract
#' @keywords SEM lavaan moderation interactions diagnostics localfit residuals multicollinearity plotting R2 multipletesting validation
"_PACKAGE"
