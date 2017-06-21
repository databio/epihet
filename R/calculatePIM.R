#' Given a Bisulfite data.table (BSDT), calculates the proportion
#' of intermediate methylation sites.
#'
#' TODO: Make the alpha level a parameter, both for confidence intervals,
#' and also for IM definition
#'
#' @param cache Logical indicating whether or not to use caching via \code{\link{simpleCache}}; default is TRUE
#' @param cacheDir If using caching, this argument specifies the directory to use for storing the cache; defaults to global option for \code{RESOURCES.RACHE}, if no such option has been specified you must provide one
#' @export

calculatePIM = function(BSDT, cache = TRUE, cacheDir = getOption("RESOURCES.RCACHE")) {

  imtab = prepIM(BSDT, cache = cache, cacheDir = cacheDir)

  sum(imtab$IM == TRUE) / nrow(imtab)

}
