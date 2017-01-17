# SOURCE FILE

#=================================================================================*
# ---- SET-UP ----
#=================================================================================*

smartRequire <- function(packages){
  # Get vector of packages listed that aren't on the current computer
  packagesToInstall <- packages[!package %in% rownames(installed.packages())]
  if(length(packagesToInstall > 0)){
    install.packages(packagesToInstall, dependencies = TRUE)
  }
  lapply(packages, require, character.only = TRUE)
}

