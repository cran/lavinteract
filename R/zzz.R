#' @noRd
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("lavinteract: post-estimation utilities for 'lavaan' fitted models.")
  packageStartupMessage("For more information, type ?lavinteract.")
  packageStartupMessage(
    "For suggestions or to report issues, please contact Giuseppe Corbelli at ",
    "giuseppe.corbelli@uniroma1.it or giuseppe.corbelli@uninettunouniversity.net."
  )
}
