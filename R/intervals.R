#' Calculate credibility interval based on binomial bayes distribution
#'
#' @param bsData Bisulfite sequencing data
#' @param methylCol Name of column containing methylation count; defaults to "methylCount"
#' @param coverageCol Name of column containing coverage (i.e. number of reads); defaults to "coverage"
#' @param confLevel The level of confidence to be used in the confidence interval; default is 0.95
#' @export
BScredInterval = function(bsData,
                            methylCol="methylCount",
                            coverageCol="coverage",
                            confLevel=.95) {

    conf = binom::binom.bayes(bsData[,get(methylCol)],
                                bsData[,get(coverageCol)],
                                conf.level = confLevel,
                                tol=.005,
                                type="central")

    conf = data.table::data.table(conf)

    bsData = cbind(conf, bsData)

    bsData

}

#' Calculate credibility interval based on binomial bayes distribution with caching
#'
#' @param bsData Bisulfite sequencing data
#' @param methylCol Name of column containing methylation count; defaults to "methylCount"
#' @param coverageCol Name of column containing coverage (i.e. number of reads); defaults to "coverage"
#' @param cachedBinomialIntervals cachedBinomialIntervals
#' @param confLevel The level of confidence to be used in the confidence interval; default is 0.95
#' @export
BScredIntervalCache = function(bsData,
                                cachedBinomialIntervals,
                                methylCol="methylCount",
                                coverageCol="coverage",
                                confLevel=.95){

    storeKey = data.table::key(bsData)

    if(length(storeKey) != 2) {

        message("Key temporarily set to ", methylCol, " and ", coverageCol)

        data.table::setkeyv(bsData, c(methylCol, coverageCol))

    }

    keepCols = colnames(bsData)

    # Use cache where you can

    bsData = cachedBinomialIntervals[bsData,]

    data.table::setnames(bsData,c("x", "n"), c(methylCol, coverageCol))

    #And otherwise, count directly.

    if(nrow(bsData[is.na(upper),]) > 0) {

        a = BScredInterval(bsData[is.na(upper),keepCols, with=FALSE])

        bsData[is.na(upper),] = a[,colnames(bsData), with=FALSE]

    }

    bsData

}

#' Cache binomial confidence intervals
#'
#' @param maxHits maxHits
#' @param maxTotal maxTotal
#' @param confLevel The level of confidence to be used in the confidence interval; default is 0.95
cacheBinomConfIntervals = function(maxHits, maxTotal, confLevel) {

    allComb = cbind(hits=rep(0:maxHits, each=maxTotal),
                            total=rep(1:maxTotal))

    allPossibleComb = allComb[allComb[,"hits"] <= allComb[,"total"],]

    conf = binom::binom.bayes(allPossibleComb[,"hits"],
                                allPossibleComb[,"total"],
                                conf.level = confLevel,
                                tol=.005,
                                type="central")

    confdt = data.table::data.table(conf)

    data.table::setkey(confdt, "x", "n")

    confdt

}
