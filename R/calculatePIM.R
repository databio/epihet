#' Given a Bisulfite data.table (BSDT), calculates the proportion
#' of intermediate methylation sites.
#'
#' TODO: Make the alpha level a parameter, both for confidence intervals,
#' and also for IM definition
#' @export

calculatePIM = function(BSDT) {

  imtab = genIM(BSDT)

  sum(imtab$IM == TRUE) / nrow(imtab)

}
