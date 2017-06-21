#' Given a Bisulfite data.table (BSDT), calculates the proportion
#' of intermediate methylation sites.
#'
#' TODO: Make the alpha level a parameter, both for confidence intervals,
#' and also for IM definition
#'
#' @param cache Logical indicating whether or not to use caching via \code{\link{simpleCache}}; default is TRUE
#' @export

calculatePIM = function(BSDT, cache = TRUE) {

  imtab = prepIM(BSDT, cache = cache)

  sum(imtab$IM == TRUE) / nrow(imtab)

}
