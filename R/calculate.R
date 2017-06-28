#' Given bisulfite sequencing data, calculates the proportion
#' of intermediate methylation sites.
#'
#' @param bsData Bisulfite sequencing data;
#' @param cacheDir If using caching, this argument specifies the directory
#' to use for storing the cache;
#' defaults to global option for \code{RESOURCES.RACHE};
#' if no such option has been specified you must provide one
#' @param imLower The lower boundary for intermediate methylation (IM);
#' if a site is entirely below this threshold
#' (or if any part of a its binomial credibilty interval overlaps this boundary)
#' it is not considered IM;
#' defaults to .25
#' @param imUpper The upper boundary for intermediate methylation (IM);
#' if a site is entirely above this threshold
#' (or if any part of a its binomial credibilty interval overlaps this boundary) it is not considered IM;
#' defaults to .75
#' @param confLevel A decimal indicating the level of confidence
#' to be used while creating cached the binomial bayes credibility interval;
#' default is .95 for 95 percent confidence
#' @export

PIM = function(bsData,
                cacheDir = getOption("RESOURCES.RCACHE"),
                imLower = 0.25,
                imUpper = 0.75,
                confLevel = .95) {

    bsData = bsDataCheck(bsData)

    imtab = prepIM(bsData,
                    cacheDir = cacheDir,
                    imLower = imLower,
                    imUpper = imUpper,
                    confLevel = confLevel)

    sum(imtab$IM == TRUE) / nrow(imtab)

}

#' Helper function to get the relative proportion of flagged sites for a single sample versus all other samples in a list of bisulfite data tables.
#'
#' @param sampleName The sample (which should specify a name in BSDTlist) to return
#' the proportion of sites for.
#' @param bsData Bisulfite sequencing data for multiple samples; a BSDT
#' (bisulfite data.table) that has been split with splitDataTable
#' (so, a list of BSDTs); one corresponds to each sample to test.
#' @param cacheDir If using caching, this argument specifies the directory to use for storing the cache; defaults to global option for \code{RESOURCES.RACHE}, if no such option has been specified you must provide one
#' @param imLower The lower boundary for intermediate methylation (IM);
#' if a site is entirely below this threshold
#' (or if any part of a its binomial credibilty interval overlaps this boundary)
#' it is not considered IM;
#' defaults to .25
#' @param imUpper The upper boundary for intermediate methylation (IM);
#' if a site is entirely above this threshold
#' (or if any part of a its binomial credibilty interval overlaps this boundary)
#' it is not considered IM;
#' defaults to .75
#' @param confLevel A decimal indicating the level of confidence
#' to be used while creating cached the binomial bayes credibility interval;
#' default is .95 for 95 percent confidence
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

#' Get the relative proportion of flagged sites.
#' @param bsData Bisulfite sequencing data for multiple samples; a BSDT
#' (bisulfite data.table) that has been split with splitDataTable
#' (so, a list of BSDTs); one corresponds to each sample to test.
#' @param cacheDir If using caching, this argument specifies the directory to use
#' for storing the cache; defaults to global option for \code{RESOURCES.RACHE};
#' if no such option has been specified you must provide one
#' @param imLower The lower boundary for intermediate methylation (IM);
#' if a site is entirely below this threshold
#' (or if any part of a its binomial credibilty interval overlaps this boundary)
#' it is not considered IM;
#' defaults to .25
#' @param imUpper The upper boundary for intermediate methylation (IM);
#' if a site is entirely above this threshold (or if any part of a its binomial credibilty interval overlaps this boundary) it is not considered IM;
#' defaults to .75
#' @param confLevel A decimal indicating the level of confidence
#' to be used while creating cached the binomial bayes credibility interval;
#' default is .95 for 95 percent confidence
#' @export
RPIM = function(bsData,
                cacheDir = getOption("RESOURCES.RCACHE"),
                imLower = .25,
                imUpper = .75,
                confLevel = .95) {

    bsData = bsDataCheck(bsData)

    x = sapply(names(bsData),
                calculateRPIM,
                bsData,
                cacheDir,
                imLower,
                imUpper,
                confLevel)

    diag(x) = NA

    colMeans(x, na.rm = TRUE)

}
