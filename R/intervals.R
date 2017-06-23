#' Calculate credibility interval based on binomial bayes distribution
#'
#' @param BSDT Bisulfite sequencing data in a data.table format
#' @param hitCol Name of column containing methylation count; defaults to "methylCount"
#' @param readCol Name of column containing coverage; defaults to "coverage"
#' @param confLevel The level of confidence to be used in the confidence interval; default is 0.95
#' @export
BScredInterval = function(BSDT, hitCol="methylCount", readCol="coverage", confLevel=.95) {

  conf = binom::binom.bayes(BSDT[,get(hitCol)], BSDT[,get(readCol)], conf.level = confLevel, tol=.005, type="central")

  conf = data.table::data.table(conf)

  BSDT = cbind(conf, BSDT)

  BSDT

}

#' Calculate credibility interval based on binomial bayes distribution with caching
#'
#' @param BSDT Bisulfite sequencing data in a data.table format
#' @param hitCol Name of column containing methylation count; defaults to "methylCount"
#' @param cachedBinomialIntervals cachedBinomialIntervals
#' @param readCol Name of column containing coverage; defaults to "coverage"
#' @param confLevel The level of confidence to be used in the confidence interval; default is 0.95
#' @export
BScredIntervalCache = function(BSDT, cachedBinomialIntervals, hitCol="methylCount", readCol="coverage", confLevel=.95){
  storeKey = data.table::key(BSDT)

  if(length(storeKey) != 2) {

    message("Key temporarily set to ", hitCol, " and ", readCol)

    data.table::setkeyv(BSDT, c(hitCol, readCol))

  }

  keepCols = colnames(BSDT)

  # Use cache where you can

  BSDT = cachedBinomialIntervals[BSDT,]

  data.table::setnames(BSDT,c("x", "n"), c(hitCol, readCol))

  #And otherwise, count directly.

  if(nrow(BSDT[is.na(upper),]) > 0) {

    a = BScredInterval(BSDT[is.na(upper),keepCols, with=FALSE])

    BSDT[is.na(upper),] = a[,colnames(BSDT), with=FALSE]

  }

  BSDT

}

#' Cache binomial confidence intervals
#'
#' @param maxHits maxHits
#' @param maxTotal maxTotal
#' @param confLevel The level of confidence to be used in the confidence interval; default is 0.95
#' @export
cacheBinomConfIntervals = function(maxHits, maxTotal, confLevel) {

  allCombinations = cbind(hits=rep(0:maxHits, each=maxTotal), total=rep(1:maxTotal))

  allPossibleCombinations = allCombinations[allCombinations[,"hits"] <= allCombinations[,"total"],]

  conf = binom::binom.bayes(allPossibleCombinations[,"hits"], allPossibleCombinations[,"total"], conf.level = confLevel, tol=.005, type="central")

  confdt = data.table::data.table(conf)

  data.table::setkey(confdt, "x", "n")

  confdt

}
