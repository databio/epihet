BScredIntervalCache = function(BSDT, cachedBinomialIntervals, hitCol="hitCount", readCol="readCount", confLevel=.95){
  storeKey = data.table::key(BSDT)

  if(length(storeKey) != 2) {

    message("Key temporarily set to ", hitCol, " and ", readCol)

    data.table::setkeyv(BSDT, c(hitCol, readCol))

  }

  keepCols = colnames(BSDT)

  # Use cache where you can

  BSDT = cachedBinomialIntervals[BSDT,]

  data.table::setnames(BSDT,c("x", "n"), c("hitCount", "readCount"))

  #And otherwise, count directly.

  if(nrow(BSDT[is.na(upper),]) > 0) {

    a= BScredInterval(BSDT[is.na(upper),keepCols, with=FALSE])

    BSDT[is.na(upper),] = a[,colnames(BSDT), with=FALSE]

  }

  BSDT
}

