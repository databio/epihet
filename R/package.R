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
#'
#' @references \url{http://github.com/databio/RPIM}
#' @importFrom data.table ":="

NULL

# Because of some issues with CRAN submission,
# (see here: http://stackoverflow.com/questions/9439256/)
# I have to register stuff used in data.table as non-standard evaluation,
# in order to pass some R CMD check NOTES.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".", "upper", "lower", "..keepCols", ".N", "IM.x", "IM.y"))
}
