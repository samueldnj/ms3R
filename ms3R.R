# <><><><><><><><><><><><><><><><><><><><><><><><><><><>
# testSimSCAL.R
#
# Test script for ms3R operating model
#
# Author: SDN Johnson
# Date: July 25, 2019
#
# Last Update: October 8, 2019
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><>

source("loadPackages.R")
source("ms3RrefPts.R")
source("tools.R")
source("simSCAL.R")
source("ms3Rplots.R")
source("applyTMBphase.R")
source("batchTools.R")
source("stats.R")


# TMB::compile("hierProd.cpp")
# dyn.load(dynlib("hierProd"))
