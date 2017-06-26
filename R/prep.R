#' Given bisulfite sequencing data, prepares the intermediate methylation (IM) table
#' @param bsData Bisulfite sequencing data
#' @param cacheDir If using caching, this argument specifies the directory to use for storing the cache; defaults to global option for \code{RESOURCES.RACHE},
#' if no such option has been specified you must provide one
#' @param imLower The lower boundary for intermediate methylation (IM); if a site is entirely below this threshold (or if any part of a its binomial credibilty interval overlaps this boundary) it is not considered IM; defaults to .25
#' @param imUpper The upper boundary for intermediate methylation (IM); if a site is entirely above this threshold (or if any part of a its binomial credibilty interval overlaps this boundary) it is not considered IM; defaults to .75
#' @export
prepIM = function(bsData,
                    cacheDir = getOption("RESOURCES.RCACHE"),
                    imLower = .25,
                    imUpper = .75) {

    # allowed = c("data.table", "list", "BSseq")
    #
    # if(!inherits(bsData, allowed))
    #     stop(c("the following are allowed:\n",
    #             paste0(allowed, collapse = "\n")))

    # stopifnot(inherits(bsData, allowed))

    # if(inherits(bsData, "BSseq"))
    #     bsData = MIRA::bsseqToDataTable(bsData)

    # if(inherits(bsData, "BSseq")) {
    #     bsData = MIRA::bsseqToDataTable(bsData)
    # } else {
    #     bsData
    # }

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

    # only keep columns if they exist in input data
    IM = CI[, IM := !(upper < imLower | lower > imUpper) ]

    keepCols = intersect(colnames(CI), c("chr", "start", "id", "IM"))

    IM = CI[, ..keepCols]

    IM

}

#' Helper function to check input Bisulfite sequencing data and convert as necessary
#' @param bsData Bisulfite sequencing data
bsDataCheck = function(bsData) {

    allowed = c("data.table", "list", "BSseq")

    if(!inherits(bsData, allowed))
        stop(c("the following are allowed:\n",
               paste0(allowed, collapse = "\n")))

    if(inherits(bsData, "BSseq"))
        bsData = MIRA::bsseqToDataTable(bsData)

    bsData
}
