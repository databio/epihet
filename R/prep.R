#' Given bisulfite sequencing data, prepares the intermediate methylation (IM) table
#' @param bsData Bisulfite sequencing data
#' @param cacheDir If using caching, this argument specifies the directory to use for storing the cache; defaults to global option for \code{RESOURCES.RACHE}, if no such option has been specified you must provide one
#' @export
prepIM = function(bsData, cacheDir = getOption("RESOURCES.RCACHE")) {

    if (requireNamespace("simpleCache", quietly=TRUE)) {
        simpleCache::simpleCache("cachedBinomialIntervals95", {
        cachedBinomialIntervals95 = cacheBinomConfIntervals(2000, 2000, .95)
      }, cacheDir = cacheDir)

        cachedBinomialIntervals = cachedBinomialIntervals95

        # Make the memory use smaller by eliminating unnecessary columns
        cachedBinomialIntervals[, c("method","mean","shape1","shape2","sig"):=NULL]

        # Calculate the credibility interval
        CI = BScredIntervalCache(bsData, cachedBinomialIntervals)

    } else {
      message("install simplecache to speed up...")
      CI = BScredInterval(bsData)
    }

    data.table::setkey(CI, "chr", "start")

    # Prep matrix for relative PIM calculation
    # We define a site as IM (Intermediate Methylation) if its credibility
    # interval is not completely below .25, or above .75. Other sites are
    # more likely to be either 0 or 1 (or very close).

    # only keep columns if they exist in input data
    IM = CI[, IM := !(upper < .25 | lower > .75) ]
    keepCols = intersect(colnames(CI), c("chr", "start", "id", "IM"))

    IM = CI[, ..keepCols]

    IM

}
