#' Helper function to get the relative proportion of flagged sites for a single sample versus all other samples in a list of bisulfite data tables.
#'
#' @param sampleName The sample (which should specify a name in BSDTlist) to return
#' the proportion of sites for.
#' @param BSDTsplit A BSDT (bisulfite data.table) that has been split with
#' splitDataTable (so, a list of BSDTs); one corresponds to each sample to test.
#' @param cache Logical indicating whether or not to use caching via \code{\link{simpleCache}}; default is TRUE
#' @param cacheDir If using caching, this argument specifies the directory to use for storing the cache; defaults to global option for \code{RESOURCES.RACHE}, if no such option has been specified you must provide one
calculateRPIM = function(sampleName, BSDTsplit, cache = TRUE, cacheDir = getOption("RESOURCES.RCACHE")) {

  message(sampleName)

  result = vector()

  sampleBaseline = prepIM(BSDTsplit[[sampleName]], cache = cache, cacheDir = cacheDir)

  for (y in names(BSDTsplit)) {

    sampleRelative = prepIM(BSDTsplit[[y]], cache = cache, cacheDir = cacheDir)
    result[y] = merge(sampleBaseline, sampleRelative)[,log(sum(IM.x/.N)/sum(IM.y/.N))]

  }
  return(result)
}

#' Get the relative proportion of flagged sites. This is a general version
#' of a method to get the RPIM (Relative Proportion of Intermediate Methylation).
#' Given a DT with bisulfite reads, and a flag column,
#' and then a huge table with this data for lots of samples, this will calculate
#' the relative proportion of CGs that are flagged, for pairwise comparisons
#' between all samples, subsetted to the CGs present in both samples.
#' @export
#' @param BSDTsplit A BSDT (bisulfite data.table) that has been split with
#' splitDataTable (so, a list of BSDTs); one corresponds to each sample to test.
#' @param cache Logical indicating whether or not to use caching via \code{\link{simpleCache}}; default is TRUE
#' @param cacheDir If using caching, this argument specifies the directory to use for storing the cache; defaults to global option for \code{RESOURCES.RACHE}, if no such option has been specified you must provide one
getRPIM = function(BSDTsplit, cache = TRUE, cacheDir = getOption("RESOURCES.RCACHE")) {

  x = sapply(names(BSDTsplit), calculateRPIM, BSDTsplit, cache, cacheDir)

  diag(x) = NA

  colMeans(x, na.rm = T)

}
