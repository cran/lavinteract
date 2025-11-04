#' Post-Estimation Utilities for 'lavaan' Fitted Models
#'
#' Companion toolbox for structural equation models fitted with 'lavaan'.
#' Operates directly on a fitted object using its estimates and covariance.
#' Refits auxiliary models when needed to compute estimates, diagnostics, and plots. 
#'
#' The functions are:
#' \itemize{
#'   \item \code{\link{lav_slopes}}: simple slopes and interaction plots from a fitted 'lavaan' model.
#'   \item \code{\link{lav_vif}}: variance inflation factors for structural predictors with measurement preserved.
#' }
#'
#' @section Note:
#' The development of this package grew from ongoing discussions and interactions (sic)
#' with colleagues, in particular Dr. Cataldo Giuliano Gemmano, whose steady 
#' feedback and support helped shape it. 
#'
#' @author
#' Giuseppe Corbelli (<giuseppe.corbelli@uniroma1.it>, or <giuseppe.corbelli@uninettunouniversity.net>)
#'
#' @docType package
#' @name lavinteract
#' @aliases lavinteract-package lavinteract
#' @keywords SEM lavaan moderation interactions diagnostics multicollinearity plotting R2 
"_PACKAGE"
