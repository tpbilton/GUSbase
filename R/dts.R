#### Some functions from the kutils package for removing trailing spaces for filenames.
#' @keywords internal
#' @export dts
dts <- function (name)
  gsub("/$", "", dms(name))

dms <- function(name)
  gsub("(/)\\1+", "/", name)
