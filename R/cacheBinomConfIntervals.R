#' @export

cacheBinomConfIntervals = function(maxHits, maxTotal, confLevel) {

  allCombinations = cbind(hits=rep(0:maxHits, each=maxTotal), total=rep(1:maxTotal))

  allPossibleCombinations = allCombinations[allCombinations[,"hits"] <= allCombinations[,"total"],]

  conf = binom::binom.bayes(allPossibleCombinations[,"hits"], allPossibleCombinations[,"total"], conf.level = confLevel, tol=.005, type="central")

  confdt = data.table(conf)

  setkey(confdt, "x", "n")

  return(confdt)
}
