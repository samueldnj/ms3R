# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# loadPackages.R
# 
# Checks if required packages are installed. If not, installs them.
# Then loads all required packages.
# 
# 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

cranPackages <- c("coda",
                  "tidyverse",
                  "TMB",
                  "RColorBrewer",
                  "parallel",
                  "stringr",
                  "wesanderson",
                  "scales",
                  "beepr",
                  "tmbstan",
                  "here",
                  "bookdown",
                  "ClimProjDiags",
                  "kableExtra" )

for( pkg in cranPackages )
  while(!require(pkg, character.only = TRUE) )
    install.packages( pkg, repos = "https://cloud.r-project.org/" )


# githubPackages <- c(ggsidekick = "seananderson/ggsidekick",
#                     csasdown = "pbs-assess/csasdown" )

# for( pkgIdx in 1:length(githubPackages) )
#   while(!require(names(githubPackages)[pkgIdx], character.only = TRUE))
#     devtools::install_github(githubPackages[pkgIdx])
