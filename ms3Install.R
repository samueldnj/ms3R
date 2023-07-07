# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# ms3Install.R
# 
# Installs ms3 closed loop simulation package by setting up output
# directories
# 
# 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

if(!dir.exists("Outputs"))
  dir.create("Outputs")

# Load packages
source("loadPackages.R")


# Compile AMs
message("Compiling assessment models\n\n")
compile("hierProd.cpp", flags = "")
message("\n\n Hierarchical Surplus Production model compiled \n\n")


message("\n\n Installation complete \n\n")
beepr::beep()