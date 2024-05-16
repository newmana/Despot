source("sptranr/R/_Loading.R")
source("sptranr/R/transpar.R")
source("sptranr/R/_scRNA-seq.R")
Check_Load_GithubPackages("BASS", "zhengli09/BASS")

Cluster_BASS <- function(BASS){
  listAllHyper(BASS)
  BASS <- BASS.preprocess(BASS)
  BASS <- BASS.run(BASS)
  return(BASS)
}