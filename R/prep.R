#' Prepares the intermediate methylation (IM) table
#'
#' @param bsData Bisulfite sequencing data
#' @param cacheDir If using caching, this argument specifies the directory to
#' use for storing the cache;
#' defaults to global option for \code{RESOURCES.RACHE},
#' if no such option has been specified you must provide one
#' @param imLower The lower boundary for intermediate methylation (IM);
#' if a site is entirely below this threshold
#' (or if any part of its binomial credibility interval overlaps this boundary)
#' it is not considered IM;
#' default is .25
#' @param imUpper The upper boundary for intermediate methylation (IM);
#' if a site is entirely above this threshold
#' (or if any part of its binomial credibility interval overlaps this boundary)
#' it is not considered IM;
#' default is .75
#' @param confLevel A decimal indicating the level of confidence
#' to be used while creating cached the binomial bayes credibility interval;
#' default is .95 for 95 percent confidence
#'
#'@return A \code{data.table} object with the following columns:
#'\itemize{
#'  \item{chr} {chromosome of methylation read}
#'  \item{start} {starting position for methylation read}
#'  \item{IM} {boolean indicator of itermediate methylation status}
#'  }
#'
#'@examples
#'
#'data("exampleBSDT", package = "RPIM")
#'
#'prepIM(exampleBSDT)
#'
#' @export
prepIM = function(bsData,
                    cacheDir = getOption("RESOURCES.RCACHE"),
                    imLower = .25,
                    imUpper = .75,
                    confLevel = .95) {

    binomCacheName = paste0("cachedBinomialIntervals", round(confLevel*100))

    if (requireNamespace("simpleCache", quietly=TRUE)) {
        simpleCache::simpleCache(binomCacheName, {
            res = cacheBinomConfIntervals(2000, 2000,confLevel)
        }, cacheDir = cacheDir, buildEnvir = list(confLevel = confLevel))

        cachedBinomialIntervals = eval(parse(
            text = paste0("cachedBinomialIntervals", round(confLevel*100))))

        # Make the memory use smaller by eliminating unnecessary columns
        cachedBinomialIntervals[, c("method","mean","shape1","shape2","sig"):=NULL]

        # Calculate the credibility interval
        CI = BScredIntervalCache(bsData,
                                cachedBinomialIntervals,
                                confLevel = confLevel)

    } else {

        message("install simplecache to speed up...")

        CI = BScredInterval(bsData, confLevel = confLevel)
    }

    data.table::setkey(CI, "chr", "start")

    # only keep columns if they exist in input data
    IM = CI[, IM := !(upper < imLower | lower > imUpper) ]

    keepCols = intersect(colnames(CI), c("chr", "start", "id", "IM"))

    IM = CI[, ..keepCols]

    IM

}

#' Check bisulfite sequencing data and convert if needed
#'
#' @param bsData Bisulfite sequencing data
#'
#'@return A \code{list} of \code{data.table} objects that each contain bisulfite
#' sequencing data
bsDataCheck = function(bsData) {

    allowed = c("data.table", "list", "BSseq")

    if(!inherits(bsData, allowed))
        stop(c("the following are allowed:\n",
                paste0(allowed, collapse = "\n")))

    if(inherits(bsData, "BSseq"))
        bsData = MIRA::bsseqToDataTable(bsData)

    bsData

}
