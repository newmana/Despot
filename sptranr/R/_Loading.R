
Check_Load_BiocPackages <- function(pkgName, version=NA){
  if(!require(pkgName, character.only = TRUE)){
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
      library("BiocManager")
    BiocManager::install(pkgName, update=FALSE, version=version)
    library(pkgName, character.only = TRUE)
  }
}

Check_Load_GithubPackages <- function(pkgName, URL){
  if(!require(pkgName, character.only = TRUE)){
    if(!require("devtools", quietly = TRUE))
      install.packages("devtools")
      library("devtools")
    install_github(URL, upgrade = "default", quietly = TRUE, build_vignettes = F)
    library(pkgName, character.only = TRUE)
  }
}

Check_Load_InstallPackages <- function(pkgName){
  if(!require(pkgName, character.only = TRUE)){
    install.packages(pkgName, quiet = TRUE)
    library(pkgName, character.only = TRUE)
  }
}
