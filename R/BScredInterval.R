#' @export
BScredInterval = function(BSDT, hitCol="hitCount", readCol="readCount", confLevel=.95) {

  conf = binom::binom.bayes(BSDT[,get(hitCol)], BSDT[,get(readCol)], conf.level = confLevel, tol=.005, type="central")

  conf = data.table(conf)

  BSDT = cbind(conf, BSDT)

  BSDT

}
