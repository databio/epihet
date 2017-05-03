devtools::install_github("databio/MIRA")
devtools::install_github("databio/simpleCache")

library(RPIM)

dat <- MIRA::BSreadBiSeq("data/RRBS_cpgMethylation_EWS_L10.bed")

# have to setSharedCacheDir()
simpleCache::setSharedCacheDir("cache")

imres <- calculatePIM(dat)

###############################

# relative proportion of sites

dat2 = MIRA::BSreadBiSeq("data/RRBS_cpgMethylation_EWS_T133.bed")
dat3 = MIRA::BSreadBiSeq("data/RRBS_cpgMethylation_EWS_T111.bed")
dat4 = MIRA::BSreadBiSeq("data/RRBS_cpgMethylation_EWS_T120.bed")

alldat = rbind(dat,dat2, dat3, dat4)
allsplitdat = split(alldat, alldat$sampleName)

calculateRPIM("RRBS_cpgMethylation_EWS_L10", allsplitdat)

# Rivanna demo:
# install.packages("~/code/RPIM", repos=NULL)
# library(RPIM)
# file = paste0(Sys.getenv("PROCESSED"), "/ews_patients/results_pipeline/EWS_L10/biseq_hg38/RRBS_cpgMethylation_EWS_L10.bed")
# #file = "data/RRBS_cpgMethylation_EWS_L10.bed"
# dat <- BSreadBiSeq(file)
# imres <- calculatePIM(dat)
# imres
# sum(imres$IM == TRUE) / nrow(imres)
