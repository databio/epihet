#' Given bisulfite sequencing data, calculates the proportion
#' of intermediate methylation sites.
#'
#' @param bsData Bisulfite sequencing data;
#' @param cacheDir If using caching, this argument specifies the directory to
#'     use for storing the cache; defaults to global option for
#'     \code{RESOURCES.RACHE}; if no such option has been specified you must
#'     provide one
#' @param imLower The lower boundary for intermediate methylation (IM); if a
#'     site is entirely below this threshold (or if any part of its binomial
#'     credibility interval overlaps this boundary) it is not considered IM;
#'     defaults to .25
#' @param imUpper The upper boundary for intermediate methylation (IM); if a
#'     site is entirely above this threshold (or if any part of its binomial
#'     credibility interval overlaps this boundary) it is not considered IM;
#'     defaults to .75
#' @param confLevel A decimal indicating the level of confidence to be used
#'     while creating cached the binomial bayes credibility interval; default is
#'     .95 for 95 percent confidence
#'
#' @return A single value (numeric vector of length 1) indicating the proportion
#'     of intermediate methylation (PIM) of an individual sample
#'
#'@examples
#'
#'data("exampleBSDT", package="epihet")
#'
#'if (require(simpleCache)) {
#'PIM(exampleBSDT, cacheDir="~")
#'
#'simpleCache::setSharedCacheDir("cache")
#'PIM(exampleBSDT)
#'PIM(exampleBSDT, cacheDir="~", imLower=.2, imUpper=.8)
#'}
#' @export

PIM = function(bsData,
                cacheDir=getOption("RESOURCES.RCACHE"),
                imLower=0.25,
                imUpper=0.75,
                confLevel=.95) {

    bsData = bsDataCheck(bsData)

    if(!singleSample(bsData)) {

        stop(strwrap("Your data appears to include more than one sample.
        Consider reformatting or try using RPIM() to calculate relative
        proportion of intermediate methylation across all samples.", initial="",
        prefix=" ")

    }

    imtab = prepIM(bsData,
                    cacheDir=cacheDir,
                    imLower=imLower,
                    imUpper=imUpper,
                    confLevel=confLevel)

    sum(imtab$IM == TRUE) / nrow(imtab)

}

#' Helper function to get the relative proportion of flagged sites for a
#' single sample versus all other samples
#'
#' @param sampleName The sample (which should specify a name in the
#' bisulfite sequencing data) to use as the baseline
#' the proportion of sites for.
#' @param bsData Bisulfite sequencing data for multiple samples; a BSDT
#' (bisulfite data.table) that has been split with splitDataTable
#' (so, a list of BSDTs); one corresponds to each sample to test.
#' @param cacheDir If using caching, this argument specifies the directory to
#' use for storing the cache;
#' defaults to global option for \code{RESOURCES.RACHE},
#' if no such option has been specified you must provide one
#' @param imLower The lower boundary for intermediate methylation (IM);
#' if a site is entirely below this threshold
#' (or if any part of its binomial credibility interval overlaps this boundary)
#' it is not considered IM;
#' defaults to .25
#' @param imUpper The upper boundary for intermediate methylation (IM);
#' if a site is entirely above this threshold
#' (or if any part of its binomial credibility interval overlaps this boundary)
#' it is not considered IM;
#' defaults to .75
#' @param confLevel A decimal indicating the level of confidence
#' to be used while creating cached the binomial bayes credibility interval;
#' default is .95 for 95 percent confidence
#'
#'@return A vector of the same length as the number of samples being
#' analyzed; each element in the vector represents the the proportion of
#' intermediate methylation relative to the other samples for a single sample
calculateRPIM = function(sampleName,
                            bsData,
                            cacheDir = getOption("RESOURCES.RCACHE"),
                            imLower = .25,
                            imUpper = .75,
                            confLevel = .95) {

    message(sampleName)

    sampleBaseline = prepIM(bsData[[sampleName]],
                            cacheDir = cacheDir,
                            imLower = imLower,
                            imUpper = imUpper,
                            confLevel = confLevel)

    result = vector()

    for (y in names(bsData)) {

        sampleRelative = prepIM(bsData[[y]],
                                cacheDir = cacheDir,
                                imLower = imLower,
                                imUpper = imUpper,
                                confLevel = confLevel)

        result[y] = merge(sampleBaseline, sampleRelative)[,log(sum(IM.x/.N)/sum(IM.y/.N))]

    }

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
#' @param imLower The lower boundary for intermediate methylation (IM); if a
#'     site is entirely below this threshold (or if any part of its binomial
#'     credibility interval overlaps this boundary) it is not considered IM;
#'     defaults to .25
#' @param imUpper The upper boundary for intermediate methylation (IM); if a
#'     site is entirely above this threshold (or if any part of its binomial
#'     credibility interval overlaps this boundary) it is not considered IM;
#'     defaults to .75
#' @param confLevel A decimal indicating the level of confidence to be used
#'     while creating cached the binomial bayes credibility interval; default is
#'     .95 for 95 percent confidence
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
#'RPIM(BSDTlist, cacheDir="~")
#'
#'simpleCache::setSharedCacheDir("cache")
#'RPIM(BSDTlist)
#'RPIM(BSDTlist, imLower=.2, imUpper=.8)
#'}
#' @export
RPIM = function(bsData,
                cacheDir=getOption("RESOURCES.RCACHE"),
                imLower=.25,
                imUpper=.75,
                confLevel=.95) {

    bsData = bsDataCheck(bsData)

    if(singleSample(bsData)) {

        stop(strwrap("Your data appears to only include one sample. Consider
        reformatting or try using PIM() to calculate proportion of intermediate
        methylation for an individual.", initial="", prefix=" ")

    }

    x = sapply(names(bsData),
                calculateRPIM,
                bsData,
                cacheDir,
                imLower,
                imUpper,
                confLevel)

    diag(x) = NA

    colMeans(x, na.rm= RUE)

}
