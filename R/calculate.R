#' Given bisulfite sequencing data, calculates the proportion
#' of intermediate methylation sites.
#'
#' @param bsData Bisulfite sequencing data;
#' @param cacheDir If using caching, this argument specifies the directory to use for storing the cache; defaults to global option for \code{RESOURCES.RACHE}, if no such option has been specified you must provide one
#' @export

PIM = function(bsData, cacheDir = getOption("RESOURCES.RCACHE")) {

    imtab = prepIM(bsData, cacheDir = cacheDir)

    sum(imtab$IM == TRUE) / nrow(imtab)

}

#' Helper function to get the relative proportion of flagged sites for a single sample versus all other samples in a list of bisulfite data tables.
#'
#' @param sampleName The sample (which should specify a name in BSDTlist) to return
#' the proportion of sites for.
#' @param bsData Bisulfite sequencing data for multiple samples; a BSDT (bisulfite data.table) that has been split with splitDataTable (so, a list of BSDTs); one corresponds to each sample to test.
#' @param cacheDir If using caching, this argument specifies the directory to use for storing the cache; defaults to global option for \code{RESOURCES.RACHE}, if no such option has been specified you must provide one
calculateRPIM = function(sampleName, bsData, cacheDir = getOption("RESOURCES.RCACHE")) {

    message(sampleName)

    result = vector()

    sampleBaseline = prepIM(bsData[[sampleName]], cacheDir = cacheDir)

    for (y in names(bsData)) {

        sampleRelative = prepIM(bsData[[y]], cacheDir = cacheDir)
        result[y] = merge(sampleBaseline, sampleRelative)[,log(sum(IM.x/.N)/sum(IM.y/.N))]

    }
    return(result)
}

#' Get the relative proportion of flagged sites.
#' @export
#' @param bsData Bisulfite sequencing data for multiple samples; a BSDT (bisulfite data.table) that has been split with splitDataTable (so, a list of BSDTs); one corresponds to each sample to test.
#' @param cacheDir If using caching, this argument specifies the directory to use for storing the cache; defaults to global option for \code{RESOURCES.RACHE}, if no such option has been specified you must provide one
RPIM = function(bsData, cacheDir = getOption("RESOURCES.RCACHE")) {

    x = sapply(names(bsData), calculateRPIM, bsData, cacheDir)

    diag(x) = NA

    colMeans(x, na.rm = TRUE)

}
