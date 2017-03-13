
#' Get the relative proportion of flagged sites. This is a general version
#' of a method to get the RPIM (Relative Proportion of Intermediate Methylation)
#' Given a DT with bisulfite reads, and a flag column,
#' and then a huge table with this data for lots of samples, this will calculate
#' the relative proportion of CGs that are flagged, for pairwise comparisons
#' between all samples, subsetted to the CGs present in both samples.
#'
#' @param sampleName The sample (which should specify a name in BSDTlist) to return
#' the proportion of sites for.
#' @param BSDTsplit A BSDT (bisulfite data.table) that has been split with
#' splitDataTable (so, a list of BSDTs); one corresponds to each sample to test.
relativeProportionOfSites = function(sampleName, BSDTsplit) {
  message(sampleName)
  result = vector()
  for (y in names(BSDTsplit)) {
    result[y] = merge(BSDTsplit[[sampleName]], BSDTsplit[[y]])[,log(sum(IM.x/.N)/sum(IM.y/.N))]
  }
  return(result)
}
