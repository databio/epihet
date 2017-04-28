#' Given a Bisulfite data.table (BSDT), prepares the intermediate methylation (IM) table
#' @export
prepIM = function(BSDT) {

  # Grab (or create) the binomial confidence intervals
  simpleCache("cachedBinomialIntervals95", {
    cachedBinomialIntervals95 = cacheBinomConfIntervals(2000, 2000, .95)
  }, cacheDir = getOption("RESOURCES.RCACHE"))
  cachedBinomialIntervals = cachedBinomialIntervals95

  # Make the memory use smaller by eliminating unnecessary columns
  BSDT[,sampleName:=NULL]
  cachedBinomialIntervals[, method:=NULL]
  cachedBinomialIntervals[, mean:=NULL]
  cachedBinomialIntervals[, shape1:=NULL]
  cachedBinomialIntervals[, shape2:=NULL]
  cachedBinomialIntervals[, sig:=NULL]
  cachedBinomialIntervals

  # Calculate the credibility interval
  CI = BScredIntervalCache(BSDT, cachedBinomialIntervals)

  setkey(CI, "chr", "start")

  # Prep matrix for relative PIM calculation
  # We define a site as IM (Intermediate Methylation) if its credibility
  # interval is not completely below .25, or above .75. Other sites are
  # more likely to be either 0 or 1 (or very close).

  # only keep columns if they exist in input data
  IM = CI[, IM := !(upper < .25 | lower > .75) ]
  keepCols = intersect(colnames(CI), c("chr", "start", "id", "IM"))

  IM = CI[, ..keepCols]

  IM
}
