#' @export
BScredInterval = function(bsdt, hitCol="hitCount", readCol="readCount", confLevel=.95) {
  library(binom);
  conf = binom.bayes(bsdt[,get(hitCol)], bsdt[,get(readCol)], conf.level = confLevel, tol=.005, type="central")
  conf = data.table(conf)
  bsdt = cbind(conf, bsdt)
  bsdt;
}
