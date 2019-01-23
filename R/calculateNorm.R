#' Given bisulfite sequencing data, calculates the weighted
#' proportion of intermediate methylation sites.
#'
#' @param bsData Bisulfite sequencing data;
#' @param cacheDir If using caching, this argument specifies the directory to
#'     use for storing the cache; defaults to global option for
#'     \code{RESOURCES.RACHE}; if no such option has been specified you must
#'     provide one
#' @param confLevel A decimal indicating the level of confidence to be used
#'     while creating cached the binomial bayes credibility interval; default is
#'     .95 for 95 percent confidence
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
#' @return A single value (numeric vector of length 1) indicating the normalized
#'     proportion of intermediate methylation (normPIM) of an individual sample
#'
#'@examples
#'
#'data("exampleBSDT", package="epihet")
#'
#'if (require(simpleCache)) {
#'normPIM(exampleBSDT, cacheDir="~")
#'
#'simpleCache::setSharedCacheDir("cache")
#'normPIM(exampleBSDT)
#'normPIM(exampleBSDT, cacheDir="~", sdNorm=.2)
#'}
#' @export

normPIM = function(bsData,
               cacheDir=getOption("RESOURCES.RCACHE"),
               confLevel=.95,
               sdNorm=.13,
               samplingRate=100) {

  bsData = bsDataCheck(bsData)

  if(!singleSample(bsData)) {

    stop(strwrap("Your data appears to include more than one sample.
                 Consider reformatting or try using RPIM() to calculate relative
                 proportion of intermediate methylation across all samples.", initial="",
                 prefix=" "))

  }

  imtab = prepNormIM (bsData,
                 cacheDir=cacheDir,
                 confLevel=confLevel,
                 sdNorm=sdNorm,
                 samplingRate=samplingRate)

  sum(imtab$score) / nrow(imtab)

}

#' Helper function to get the relative proportion of flagged sites for a
#' single sample versus all other samples
#'
#' @param sampleBaseline The sample (which should specify a name in the
#'     bisulfite sequencing data) to use as the baseline
#'     the proportion of sites for.
#' @param sampleRelative The sample to compare relatively to the baseline sample for
#'     the proportion of sites.
#' @param bsData Bisulfite sequencing data for multiple samples; a BSDT
#'     (bisulfite data.table) that has been split with splitDataTable
#'     (so, a list of BSDTs); one corresponds to each sample to test.
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
#'@return A vector of the same length as the number of samples being
#' analyzed; each element in the vector represents the the normalized proportion of
#' intermediate methylation relative to the other samples for a single sample
calculateNormRPIM = function(sampleBaseline,
                         sampleRelative,
                         bsData,
                         cacheDir=getOption("RESOURCES.RCACHE"),
                         confLevel=.95,
                         sdNorm=.13,
                         samplingRate=100) {

  message(paste0(sampleBaseline, " to ", sampleRelative))


  sampleBaseline = prepNormIM(bsData[[sampleBaseline]],
                               cacheDir=cacheDir,
                               confLevel=confLevel,
                               sdNorm=sdNorm,
                               samplingRate=samplingRate)


  sampleRelative = prepNormIM(bsData[[sampleRelative]],
                                      cacheDir=cacheDir,
                                      confLevel=confLevel,
                                      sdNorm=sdNorm,
                                      samplingRate=samplingRate)

  result = merge(sampleBaseline, sampleRelative)[,log(sum(score.x/.N)/sum(score.y/.N))]


  return(result)

}

#' Calculate the relative proportion of intermediate methylation (RPIM) score.
#'
#' @param bsData Bisulfite sequencing data for multiple samples. A list of BSDT
#'     (bisulfite data.table), one corresponds to each sample to test. This may
#'     also be a BSseq object
#' @param cacheDir If using caching, this argument specifies the directory to
#'     use for storing the cache; defaults to global option for
#'     \code{RESOURCES.RACHE}. If no such option has been specified you must
#'     provide one
#' @param confLevel A decimal indicating the level of confidence to be used
#'     while creating cached the binomial bayes credibility interval; default is
#'     .95 for 95 percent confidence
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
#' @return A named vector of the same length as the number of samples being
#'     analyzed; each element in the vector represents the the average
#'     proportion of intermediate methylation relative to the other samples; all
#'     samples are represented and elements named accordingly
#'@examples
#'
#'data("BSDTlist", package="epihet")
#'
#'if (require(simpleCache)) {
#'normRPIM(BSDTlist, cacheDir="~")
#'
#'simpleCache::setSharedCacheDir("cache")
#'normRPIM(BSDTlist)
#'normRPIM(BSDTlist, sdNorm=.2)
#'}
#' @importFrom utils combn
#' @export
normRPIM = function(bsData,
                cacheDir=getOption("RESOURCES.RCACHE"),
                confLevel=.95,
                sdNorm=.13,
                samplingRate=100) {

  bsData = bsDataCheck(bsData)

  if(singleSample(bsData)) {

    stop(strwrap("Your data appears to only include one sample. Consider
                 reformatting or try using normPIM() to calculate proportion of intermediate
                 methylation for an individual.", initial="", prefix=" "))

  }

  mysamples = names(bsData)

  allcomb = t(combn(mysamples,2))

  revcomb = t(apply(allcomb,1,rev))

  res = combn(mysamples,
              2,
              simplify=TRUE,
              FUN=function(x) calculateNormRPIM(sampleBaseline=x[1], sampleRelative=x[2], bsData=bsData, cacheDir=cacheDir, confLevel=confLevel,
                                                sdNorm=sdNorm,samplingRate=samplingRate))

  allcomb = cbind(allcomb,res)

  revcomb = cbind(revcomb,-res)

  bothcomb = rbind(allcomb,revcomb)

  ref = expand.grid(V1=mysamples,V2=mysamples, stringsAsFactors=FALSE)

  allres = merge(ref,bothcomb,all=TRUE,stringsAsFactors=FALSE)

  allres = apply(allres, 2, as.character)

  vals = as.numeric(allres[,3])

  dim(vals) = rep(length(mysamples),2)

  colnames(vals) = sort(mysamples)

  rownames(vals) = sort(mysamples)

  colMeans(vals, na.rm=TRUE)[mysamples]

}
