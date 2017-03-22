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




# Rivanna demo:
install.packages("~/code/RPIM", repos=NULL)
library(RPIM)
file = paste0(Sys.getenv("PROCESSED"), "/ews_patients/results_pipeline/EWS_L10/biseq_hg38/RRBS_cpgMethylation_EWS_L10.bed")
#file = "data/RRBS_cpgMethylation_EWS_L10.bed"
dat <- BSreadBiSeq(file)
imres <- calculatePIM(dat)
imres
sum(imres$IM == TRUE) / nrow(imres)