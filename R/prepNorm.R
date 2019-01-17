#' Prepares the normalized intermediate methylation (normIM) table
#'
#' @param bsData Bisulfite sequencing data
#' @param cacheDir If using caching, this argument specifies the directory to
#'     use for storing the cache;
#'     defaults to global option for \code{RESOURCES.RACHE},
#'     if no such option has been specified you must provide one
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
#'@return A \code{data.table} object with the following columns:
#'\itemize{
#'  \item{chr} {chromosome of methylation read}
#'  \item{start} {starting position for methylation read}
#'  \item{score} {normalized intermediate methylation score}
#'  }
#'
prepNormIM = function(bsData,
                  cacheDir=getOption("RESOURCES.RCACHE"),
                  confLevel=.95,
                  sdNorm=.13,
                  samplingRate=100) {
  
  binomCacheName = paste0("cachedBinomialIntervals", round(confLevel*100), "NormSd", round(sdNorm*100))
  
  if (requireNamespace("simpleCache", quietly=TRUE)) {
    
    suppressMessages ({
      simpleCache::simpleCache(binomCacheName, {
        res = normCacheBinomScore(2000, 2000,confLevel, sdNorm, samplingRate)
      }, cacheDir=cacheDir, buildEnvir=list(confLevel=confLevel), loadEnvir=globalenv())
    })
    
    cachedBinomialIntervals = eval(parse(
      text = paste0("cachedBinomialIntervals", round(confLevel*100), "NormSd", round(sdNorm*100))))
    
    # Make the memory use smaller by eliminating unnecessary columns
    suppressWarnings ({
      cachedBinomialIntervals[, c("method","mean","shape1","shape2","sig", "lower", "upper"):=NULL]
    })
    
    # Calculate the credibility score
    CS = BSnormScoreCache(bsData,
                             cachedBinomialIntervals,
                             confLevel=confLevel, 
                             sdNorm=sdNorm, 
                             samplingRate=samplingRate)
    
  } else {
    
    message("install simplecache to speed up...")
    
    CS = BSnormScore(bsData,
                     confLevel=confLevel, 
                     sdNorm=sdNorm,
                     samplingRate=samplingRate)
  }
  
  data.table::setkey(CS, "chr", "start")
  
  # only keep columns if they exist in input data
  
  keepCols = intersect(colnames(CS), c("chr", "start", "id", "score"))
  
  normIM = CS[, ..keepCols]
  
  normIM
  
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

#' Check to see if bisulfite sequencing data has only one sample
#'
#' @param bsData Bisulfite sequencing data
#'
#'@return A boolean indicating whether or not the data is likely to only
#'represent a single sample
singleSample = function(bsData) {
  
  cond1 = inherits(bsData, "data.table")
  
  cond2 = is.list(bsData) & length(bsData) == 1
  
  any(cond1,cond2)
  
}
