#' Assess epigenetic heterogeneity with proportion of intermediate methylation
#'
#' Provides a score called PIM that measures the epigenetic heterogeneity in a
#'     bisulfite sequencing sample. Under the assumption that a homogeneous
#'     sample will have mostly CpGs with either 100% or 0% DNA methylation, it
#'     follows that the proportion of sites that differ from these two extremes
#'     can be used as a measure of sample heterogeneity.
#'
#'
#' @docType package
#' @name RPIM
#'
#' @references \url{http://github.com/databio/epihet}
#' @importFrom data.table ":="
NULL


# register stuff used in data.table as non-standard evaluation
# reference: http://stackoverflow.com/questions/9439256/

if (getRversion() >= "2.15.1") {

    utils::globalVariables(c(".",
                            "upper",
                            "lower",
                            "..keepCols",
                            ".N",
                            "IM.x",
                            "IM.y"))

}
