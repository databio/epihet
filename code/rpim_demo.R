library(data.table)
library(simpleCache)

# from MIRA
# parseBiseq
# BSreadBiSeq
# extractSampleName ... same as tools::file_path_sans_ext()
# lapplyAlias
# setLapplyAlias

# from RPIM
# relativeProportionOfSites
# calculatePIM
# cacheBinomConfInterval
# BScredInterval
# BScredIntervalCache

dat <- BSreadBiSeq("data/RRBS_cpgMethylation_EWS_L10.bed")

# have to setSharedCacheDir()
setSharedCacheDir("cache")

imres <- calculatePIM(dat)

# prop im?
sum(imres$IM == TRUE) / nrow(imres)
