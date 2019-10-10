# <><><><><><><><><><><><><><><><><><><><><><><><><><><>
# simSCAL.R
#
# Simulation model for generating data
# for testing a statistical catch at age 
# and length (SCAL) model, or providing
# an operating model for closed loop simulation.
#
# Author: SDN Johnson
# Date: July 25, 2019
#
# Last Update: October 8, 2019
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><>

source("refPts.R")
source("simSCAL.R")
source("ms3Rplots.R")
load("history/fit_09102019215801/fit_09102019215801.RData")

totRep <- c(reports$repOpt,reports$data)

test <- initMS3pop( totRep, pT = 0 )


diagCondition( repObj = totRep, ms3Obj = test )