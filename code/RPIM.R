#' Relative Proportion of sites with Intermediate Methylation (RPIM)
#'
#' RPIM is a score that measures the epigenetic heterogeneity in a
#' bisulfite sequencing sample. Under the assumption that a homogeneous
#' sample will have mostly CpGs with either 100% or 0% DNA methylation,
#' it follows that the proportion of sites that differ from these
#' two extremes can be used as a measure of sample heterogeneity.
#'
#' This script (an incipient R package) provides functions
#' for assessing the RPIM given input DNA methylation calls.
#'
#' @docType package
#' @name RPIM
#' @author Nathan Sheffield
#'
#' @references \url{http://github.com/nsheff}
#' @importFrom GenomicRanges GRanges GRangesList elementMetadata strand
#' @import BiocGenerics S4Vectors IRanges
#' @importFrom data.table ":=" setDT data.table setkey fread setnames as.data.table setcolorder melt setkeyv

# These functions will rely on simpleCache (http://github.com/nsheff/simpleCache)

cacheBinomConfIntervals = function(maxHits, maxTotal, confLevel) {
	library(binom);
	allCombinations = cbind(hits=rep(0:maxHits, each=maxTotal), total=rep(1:maxTotal))
	allPossibleCombinations = allCombinations[allCombinations[,"hits"] <= allCombinations[,"total"],]
	#dim(allPossibleCombinations)
	conf = binom.bayes(allPossibleCombinations[,"hits"], allPossibleCombinations[,"total"], conf.level = confLevel, tol=.005, type="central")
	confdt = data.table(conf)
	setkey(confdt, "x", "n")
	return(confdt);
}

BScredIntervalCache = function(BSDT, cachedBinomialIntervals, hitCol="hitCount", readCol="readCount", confLevel=.95) {
	storeKey = key(BSDT);
	if(length(storeKey) != 2) {
		message("Key temporarily set to ", hitCol, " and ", readCol);
		setkeyv(BSDT, c(hitCol, readCol));
	}
	keepCols = colnames(BSDT);
	#Use cache where you can,
	BSDT = cachedBinomialIntervals[BSDT,];
	setnames(BSDT,c("x", "n"), c("hitCount", "readCount"))
	#And otherwise, count directly.
	if(nrow(BSDT[is.na(upper),]) > 0) {
		a= BScredInterval(BSDT[is.na(upper),keepCols, with=FALSE])
		BSDT[is.na(upper),] = a[,colnames(BSDT), with=FALSE]
	}
	BSDT;
}

#' @export
BScredInterval = function(bsdt, hitCol="hitCount", readCol="readCount", confLevel=.95) {
	library(binom);
	conf = binom.bayes(bsdt[,get(hitCol)], bsdt[,get(readCol)], conf.level = confLevel, tol=.005, type="central")
	conf = data.table(conf)
	bsdt = cbind(conf, bsdt)
	bsdt;
}


#' Given a Bisulfite data.table (BSDT), calculates the proportion
#' of intermediate methylation sites.
#'
#' TODO: Make the alpha level a parameter, both for confidence intervals,
#' and also for IM definition
calculatePIM = function(BSDT) {

	# Grab (or create) the binomial confidence intervals
	simpleCache("cachedBinomialIntervals95", {
		cachedBinomialIntervals95 = cacheBinomConfIntervals(2000, 2000, .95)
	}, cacheDir = getOption("RESOURCES.RCACHE"))
	cachedBinomialIntervals = cachedBinomialIntervals95

	# Make the memory use smaller by eliminating unnecessary columns
	allCg[,sampleName:=NULL]
	cachedBinomialIntervals[, method:=NULL]
	cachedBinomialIntervals[, mean:=NULL]
	cachedBinomialIntervals[, shape1:=NULL]
	cachedBinomialIntervals[, shape2:=NULL]
	cachedBinomialIntervals[, sig:=NULL]
	cachedBinomialIntervals

	# Calculate the credibility interval
	CI = BScredIntervalCache(allCg, cachedBinomialIntervals)

	setkey(CI, "chr", "start")

	# Prep matrix for relative PIM calculation
	# We define a site as IM (Intermediate Methylation) if its credibility
	# interval is not completely below .25, or above .75. Other sites are
	# more likely to be either 0 or 1 (or very close).
	IM = CI[, list(chr, start, id, IM = !(upper < .25 | lower > .75)) ]

	# memory hog; clean up!
	rm(CI); rm(BSDT); gc()
	return(IM)
}


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
		result[y] = merge(sdt[[sampleName]], sdt[[y]])[,log(sum(IM.x/.N)/sum(IM.y/.N))]
	}
	return(result)
}
