source("sptranr/R/transpar.R")
source("sptranr/R/_BASS.R")
params <- fromJSON(file = "params.json")
smdFile <- params$smdFile
name <- params$name
imgdir <- paste0(params$dataPath, "/spatial")
platform <- params$platform
if(is.null(platform)){
  platform <- "10X_Visium"  # default using 10X_Visium
}

for(decont in params$Decontamination){
  h5data <- Create_smd_h5data(decont)
  if(platform != "10X_Visium" && h5data == "SpotClean_mat"){
    message("SpotClean only Support 10X Visium data, skip it.")
    next
  }
  sce <- Load_smd_to_Seurat(smdFile, h5data = h5data, platform = platform, imgdir = imgdir)
  cnt<- as.matrix(sce@assays$Spatial@counts)
  X <- list()
  X[[name]] <- cnt
  coord <- as.matrix(sce@images$slice1@coordinates[,c('row', 'col')])
  xy <- list()
  xy[[name]] <- coord

  BASS <- createBASSObject(X, xy, C=15, R=10)
  BASS <- Cluster_BASS(BASS)
  Save_smd_from_BASS(smdFile, BASS, h5data = h5data)
}


