#' Given bisulfite sequencing data, calculates the proportion
#' of intermediate methylation sites from methylation percenatage
#' histogram and selected weighting function
#'
#' @param bsData Bisulfite sequencing data;
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
#'@examples
#'
#'data("exampleBSDT", package="epihet")
#'
#'
#'histPIM(exampleBSDT)
#'histPIM(exampleBSDT, sdNorm=.2, samplingRate=200)
#'
#' @export

histPIM = function(bsData,
                   sdNorm=.13,
                   samplingRate=100){

    #get histogram of the methylation rates

    histogram = hist(bsData$methylCount / bsData$coverage, breaks = samplingRate, plot = F)

    histDens = histogram$density / 100

    # get normal distribution density function and normalize values to fall into 0-1 interval

    discretize = seq(0, 1, length=samplingRate)

    normDensity = dnorm(discretize, mean=.5, sd=sdNorm)

    normDensity = (normDensity - min(normDensity)) / (max(normDensity) - min(normDensity))

    histPIMscore = sum(histDens * normDensity) / samplingRate

    return(histPIMscore)

}

#' Helper function for the relative PIM (RPIM) score computed from
#' methylation histogram - takes coverage and methylation counts and
#' returns PIM score for these.
#' The function is used in calculateHistRPIM function.
#'
#' @param coverage Vector containing information about CpG coverage
#' @param methylCount Vector containing infomation about methylation counts
#'     at CpG
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
#'@return A PIM score based on methylation histogram and weighting function
#'
prepHistRPIM = function(coverage,
                        methylCount,
                        sdNorm=.13,
                        samplingRate=100){
    #get histogram of the methylation rates

    histogram = hist(methylCount / coverage, breaks = samplingRate, plot = F)

    histDens = histogram$density / 100

    # get normal distribution density function and normalize values to fall into 0-1 interval

    discretize = seq(0, 1, length=samplingRate)

    normDensity = dnorm(discretize, mean=.5, sd=sdNorm)

    normDensity = (normDensity - min(normDensity)) / (max(normDensity) - min(normDensity))

    histPIMscore = sum(histDens * normDensity) / samplingRate
    return(histPIMscore)
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
#' intermediate methylation relative to the other samples for a single sample calculated
#' from methylation histogram and weighting function.

calculateHistRPIM = function(sampleBaseline,
                             sampleRelative,
                             bsData,
                             sdNorm=.13,
                             samplingRate=100) {

    message(paste0(sampleBaseline, " to ", sampleRelative))


    sampleBaseline = bsData[[sampleBaseline]]
    data.table::setkey(sampleBaseline, "chr", "start")

    sampleRelative = bsData[[sampleRelative]]
    data.table::setkey(sampleRelative , "chr", "start")

    result = merge(sampleBaseline, sampleRelative)[,log(prepHistRPIM(coverage.x, methylCount.x, sdNorm = sdNorm, samplingRate = samplingRate)/prepHistRPIM(coverage.y, methylCount.y, sdNorm = sdNorm, samplingRate = samplingRate))]
    return(result)
}

#' Calculate the relative proportion of intermediate methylation (RPIM) score
#' from methylation histogram and weighting function
#'
#' @param bsData Bisulfite sequencing data for multiple samples. A list of BSDT
#'     (bisulfite data.table), one corresponds to each sample to test. This may
#'     also be a BSseq object
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
#'histRPIM(BSDTlist)
#'histRPIM(BSDTlist, sdNorm=.2)
#'
#' @importFrom utils combn
#' @export

histRPIM = function(bsData,
                    sdNorm=.13,
                    samplingRate=100){
    bsData = bsDataCheck(bsData)

    if(singleSample(bsData)) {

        stop(strwrap("Your data appears to only include one sample. Consider
                     reformatting or try using histPIM() to calculate proportion of intermediate
                     methylation for an individual.", initial="", prefix=" "))

    }

    mysamples = names(bsData)

    allcomb = t(combn(mysamples,2))

    revcomb = t(apply(allcomb,1,rev))

    res = combn(mysamples,
                2,
                simplify=TRUE,
                FUN=function(x) calculateHistRPIM(sampleBaseline=x[1], sampleRelative=x[2], bsData=bsData, sdNorm=sdNorm, samplingRate=samplingRate))

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

