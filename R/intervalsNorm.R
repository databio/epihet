#' Get the normalized score via binomial bayes distribution and normal distribution
#'
#' @param bsData Bisulfite sequencing data
#' @param methylCol Name of column containing methylation count;
#'     defaults to "methylCount"
#' @param coverageCol Name of column containing coverage (i.e. number of reads);
#'     defaults to "coverage"
#' @param confLevel A decimal indicating the level of confidence
#'     to be used while creating cached the binomial bayes credibility interval;
#'     default is .95 for 95 percent confidence
#' @param sdNorm A decimal indicating the standard deviation of the normal 
#'     distribution density function used for score weighting (the CpGs 
#'     with 50% methylation get always highest weight indicationg highest 
#'     heterogeneity, a low sdNorm gives a steeper weighting function which
#'     gives only samples close to 50% methylation level high score); 
#'     default is 0.13
#' @param samplingRate An integer deciding how finely are the generated 
#'     density distributions sampled (how many points do we get in the 
#'     interval 0-1 for each distribution); default is 100
#'
#' @return Bisulfite sequencing data as a \code{data.table} object with a column
#'     indicating score obtained by sum of Beta distribution density function (given
#'     by Byesian binomial shape parameters for given methylation) weighted by 
#'     predefined normal ditribution density function
#'
BSnormScore = function(bsData,
                          methylCol="methylCount",
                          coverageCol="coverage",
                          confLevel=.95,
                          sdNorm=.13,
                          samplingRate=100) {
  
  conf = binom::binom.bayes(bsData[,get(methylCol)],
                            bsData[,get(coverageCol)],
                            conf.level=confLevel,
                            tol=.005,
                            type="central")
  
  conf = BetaNormScore(conf, sdNorm=sdNorm, samplingRate=samplingRate)
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
#' @param sdNorm A decimal indicating the standard deviation of the normal 
#'     distribution density function used for score weighting (the CpGs 
#'     with 50% methylation get always highest weight indicationg highest 
#'     heterogeneity, a low sdNorm gives a steeper weighting function which
#'     gives only samples close to 50% methylation level high score); 
#'     default is 0.13
#' @param samplingRate An integer deciding how finely are the generated 
#'     density distributions sampled (how many points do we get in the 
#'     interval 0-1 for each distribution); default is 100
#'
#' @return Bisulfite sequencing data as a \code{data.table} object with a column
#'     indicating score obtained by sum of Beta distribution density function (given
#'     by Byesian binomial shape parameters for given methylation) weighted by 
#'     predefined normal ditribution density function
#'
BSnormScoreCache = function(bsData,
                               cachedBinomialIntervals,
                               methylCol="methylCount",
                               coverageCol="coverage",
                               confLevel=.95,
                               sdNorm=.13,
                               samplingRate=100){
  
  data.table::setkeyv(bsData, c(methylCol, coverageCol))
  
  keepCols = colnames(bsData)
  
  # Use cache where you can
  
  bsData = cachedBinomialIntervals[bsData,]
  
  data.table::setnames(bsData,c("x", "n"), c(methylCol, coverageCol))
  
  #And otherwise, count directly.
  
  if(nrow(bsData[is.na(score),]) > 0) {
    
    a = BSnormScore(bsData[is.na(score),keepCols, with=FALSE],
                       confLevel=confLevel, 
                       sdNorm=sdNorm,
                       samplingRate=samplingRate)
    
    bsData[is.na(score),] = a[,colnames(bsData), with=FALSE]
    
  }
  
  bsData
  
}

#' Cache normalized score
#'
#' @param maxHits Maximum methylation count
#' @param maxTotal Maximum coverage
#' @param confLevel A decimal indicating the level of confidence
#'     to be used while creating cached the binomial bayes credibility interval;
#'     default is .95 for 95 percent confidence
#' @param sdNorm A decimal indicating the standard deviation of the normal 
#'     distribution density function used for score weighting (the CpGs 
#'     with 50% methylation get always highest weight indicationg highest 
#'     heterogeneity, a low sdNorm gives a steeper weighting function which
#'     gives only samples close to 50% methylation level high score); 
#'     default is 0.13
#' @param samplingRate An integer deciding how finely are the generated 
#'     density distributions sampled (how many points do we get in the 
#'     interval 0-1 for each distribution); default is 100
#'
#' @return A \code{data.table} object containing a column for
#' score obtained by sum of weighted Bayesian binomial density function 
#' by normal distribution density function for maximum methylation count 
#' and coverage; this serves as a cache that can replace the need to 
#' perform the computationally expensive probability estimation
#'
#'@examples
#'
#'normCacheBinomScore(100,100, confLevel=.95, sdNorm=.2)
#'
#' @export

normCacheBinomScore = function(maxHits, 
                               maxTotal, 
                               confLevel, 
                               sdNorm,
                               samplingRate=100) {
  
  allComb = cbind(hits=rep(0:maxHits, each=maxTotal),
                  total=rep(1:maxTotal))
  
  allPossibleComb = allComb[allComb[,"hits"] <= allComb[,"total"],]
  
  conf = binom::binom.bayes(allPossibleComb[,"hits"],
                            allPossibleComb[,"total"],
                            conf.level=confLevel,
                            tol=.005,
                            type="central")
  
  confdt = data.table::data.table(conf)
  
  data.table::setkey(confdt, "x", "n")
  
  confdt = BetaNormScore(confdt, sdNorm=sdNorm, samplingRate=samplingRate)
  
  confdt
}



#' Compute the beta-norm score
#'
#' @param binomBayesDT Data table with binomial intervals
#' @param sdNorm A decimal indicating sigma of standard distribution normaliztion function
#' @param samplingRate Number of probability density function samples
#'
#' @return A \code{data.table} object containing score which is a 
#' result of a Bayesian binomial density probability function weighted
#' by normal distribution probability density function
#'
#'@examples
#'
#'conf = binom::binom.bayes(x=0:10,
#'                          n=10,
#'                          conf.level=0.95,
#'                          tol=.005,
#'                          type="central")
#'confDT = data.table::data.table(conf)
#'BetaNormScore(confDT, sdNorm=.2, samplingRate=100)
#'
#' @export


BetaNormScore = function(binomBayesDT, 
                         sdNorm,
                         samplingRate){
  
  # get normal distribution density function and normalize values to fall into 0-1 interval
  
  discretize = seq(0, 1, length=samplingRate)
  
  normDensity = dnorm(discretize, mean=.5, sd=sdNorm)
  
  normDensity = (normDensity - min(normDensity)) / (max(normDensity) - min(normDensity))
  
  # get beta density function and divide by 100: sum under the curve close to 1 (not 100)
  
  betaDist = t(apply(cbind(binomBayesDT$shape1, binomBayesDT$shape2), 
                     MARGIN=1, 
                     function(x) stats::dbeta(discretize, x[1], x[2]))) / 100
  
  # set the Inf values to 0, these are at margins and get downweighted anyways
  
  betaDist[!is.finite(betaDist)] = 0
  
  # multiply the normal density function with the beta density function and sum
  
  score = betaDist %*% normDensity
  
  binomBayesDT = cbind(binomBayesDT, score)
  
  colnames(binomBayesDT) = c(colnames(binomBayesDT)[-length(colnames(binomBayesDT))], "score")
  
  binomBayesDT
}







