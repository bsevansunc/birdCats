# SOURCE FILE

#=================================================================================*
# ---- SET-UP ----
#=================================================================================*

install.packages(tidyverse)

exists(tidyverse)


smartRequire <- function(package){
  if(!package %in% rownames(installed.packages())) {
    install.packages(package)
  }
  require(as.name(package), character.only = TRUE)
}

libraries <- c('tidyverse', 'RCurl')

for(i in 1:length(libraries)) smartRequire(libraries[i])

lapply(list('tidyverse', 'RCurl'), smartRequire)

