source("sptranr/R/_Loading.R")
Check_Load_BiocPackages("BayesSpace")
set.seed(1314)

Cluster_BayesSpace <- function(Bayes, platform){
  # Bayes <- tryCatch({scater::logNormCounts(Bayes)}, error=function(e){
  #   clusters <- quickCluster(Bayes, method = "igraph")
  #   Bayes <- computeSumFactors(Bayes)
  #   Bayes <- scran::logNormCounts(Bayes)
  #   return(Bayes)
  # })
  platform <- ifelse(platform == "ST", "ST", "Visium")
  Bayes <- spatialPreprocess(Bayes, platform=platform,
                             n.PCs=7, n.HVGs=2000, log.normalize=T)
  # Bayes <- qTune(Bayes, qs=seq(2, 10), platform=platform, d=10)
  Bayes <- spatialCluster(Bayes, q=10, platform=platform, d=10,
                          init.method="mclust", model="t", gamma=2,
                          nrep=1000, burn.in=100,
                          save.chain=TRUE)
  return(Bayes)
}

Enhance_BayesSpace <- function(Bayes){
  Bayes.enhanced <- spatialEnhance(Bayes, q=10, platform=platform, d=7,
                                   model="t", gamma=2,
                                   jitter_prior=0.3, jitter_scale=3.5,
                                   nrep=1000, burn.in=100,
                                   save.chain=T)
  return(Bayes.enhanced)
}


# clusterPlot(Bayes.enhanced,label = "spatial.cluster" )+ ggtitle("enhanced")
