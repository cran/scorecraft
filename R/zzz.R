#' @keywords internal
#' @noRd
.onAttach <- function(libname, pkgname) {
  missing <- character()
  for (p in c("glmnet", "ranger", "lightgbm", "openxlsx")) {
    if (!nzchar(system.file(package = p))) missing <- c(missing, p)  # installed? without loading it
  }
  if (length(missing)) {
    packageStartupMessage(
      "scorecraft: optional package(s) not installed: ", paste(missing, collapse = ", "),
      ".\n  - glmnet, ranger and lightgbm are optional consensus voters (xgboost is required);",
      "\n  - openxlsx is required by scr_export().")
  }
}
