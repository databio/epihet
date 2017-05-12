#' Helper function to get the relative proportion of flagged sites for a single sample versus all other samples in a list of bisulfite data tables.
#'
#' @param sampleName The sample (which should specify a name in BSDTlist) to return
#' the proportion of sites for.
#' @param BSDTsplit A BSDT (bisulfite data.table) that has been split with
#' splitDataTable (so, a list of BSDTs); one corresponds to each sample to test.
calculateRPIM = function(sampleName, BSDTsplit) {

  message(sampleName)

  result = vector()

  sampleBaseline = prepIM(BSDTsplit[[sampleName]])

  for (y in names(BSDTsplit)) {

    sampleRelative = prepIM(BSDTsplit[[y]])
    result[y] = merge(sampleBaseline, sampleRelative)[,log(sum(IM.x/.N)/sum(IM.y/.N))]

  }
  return(result)
}

#' Get the relative proportion of flagged sites. This is a general version
#' of a method to get the RPIM (Relative Proportion of Intermediate Methylation).
#' Given a DT with bisulfite reads, and a flag column,
#' and then a huge table with this data for lots of samples, this will calculate
#' the relative proportion of CGs that are flagged, for pairwise comparisons
#' between all samples, subsetted to the CGs present in both samples.
#' @export
#' @param BSDTsplit A BSDT (bisulfite data.table) that has been split with
#' splitDataTable (so, a list of BSDTs); one corresponds to each sample to test.
getRPIM = function(BSDTsplit) {

  x = sapply(names(BSDTsplit, calculateRPIM, BSDTsplit)

  diag(x) = NA

  colMeans(x, na.rm = T)

}
