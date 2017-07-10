#' Get credibility interval via binomial bayes distribution
#'
#' @param bsData Bisulfite sequencing data
#' @param methylCol Name of column containing methylation count;
#' defaults to "methylCount"
#' @param coverageCol Name of column containing coverage (i.e. number of reads);
#' defaults to "coverage"
#' @param confLevel A decimal indicating the level of confidence
#' to be used while creating cached the binomial bayes credibility interval;
#' default is .95 for 95 percent confidence
#'
#' @return Bisulfite sequencing data as a \code{data.table} object with columns
#' indicating upper and lower limits of Bayesian binomial confidence interval
#' for methylation
#'
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

#' Get credibility interval via binomial bayes distribution and caching
#'
#' @param bsData Bisulfite sequencing data
#' @param methylCol Name of column containing methylation count;
#' defaults to "methylCount"
#' @param coverageCol Name of column containing coverage (i.e. number of reads);
#' defaults to "coverage"
#' @param cachedBinomialIntervals cachedBinomialIntervals
#' @param confLevel A decimal indicating the level of confidence
#' to be used while creating cached the binomial bayes credibility interval;
#' default is .95 for 95 percent confidence
#'
#' @return Bisulfite sequencing data as a \code{data.table} object with columns
#' indicating upper and lower limits of Bayesian binomial confidence interval
#' for methylation
#'
BScredIntervalCache = function(bsData,
                                cachedBinomialIntervals,
                                methylCol="methylCount",
                                coverageCol="coverage",
                                confLevel=.95){

    data.table::setkeyv(bsData, c(methylCol, coverageCol))

    keepCols = colnames(bsData)

    # Use cache where you can

    bsData = cachedBinomialIntervals[bsData,]

    data.table::setnames(bsData,c("x", "n"), c(methylCol, coverageCol))

    #And otherwise, count directly.

    if(nrow(bsData[is.na(upper),]) > 0) {

        a = BScredInterval(bsData[is.na(upper),keepCols, with=FALSE],
                            confLevel = confLevel)

        bsData[is.na(upper),] = a[,colnames(bsData), with=FALSE]

    }

    bsData

}

#' Cache binomial confidence intervals
#'
#' @param maxHits Maximum methylation count
#' @param maxTotal Maximum coverage
#' @param confLevel A decimal indicating the level of confidence
#' to be used while creating cached the binomial bayes credibility interval;
#' default is .95 for 95 percent confidence
#'
#' @return A \code{data.table} object containing columns for
#' upper and lower limits of a Bayesian binomial
#' confidence interval for maximum methylation count and coverage;
#' this serves as a cache that can replace the need to perform the
#' computationally expensive probability estimation
#'
#'@examples
#'
#'cacheBinomConfIntervals(100,100, confLevel = .95)
#'
#' @export
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
