#-----------------------------------------------------
# Generate reference point curves for different areas 
# and different MPs for HG Herring fisheries in HG
#
# Created: B Doherty
# Date: Sept 28, 2020
#
#-----------------------------------------------------

# using SISCA refpts for now, since ms3R

source('ms3R.r')

Load ms3R sim with base OM
# load('~/Documents/LANDMARK/projects/2020_Herring_SISCA/ms3R/Outputs/sim_22092020100432/sim_22092020100432.RData')

# function for generating ms3r object with data structure to run calcRefPts
ms3Obj <- function(ctlFile = "refPtCalcs/simCtlFile_RP.txt")
{
  # Load control list
  ctlTable <- .readParFile ( ctlFile )
  # Create control list
  ctlList  <- .createList  ( ctlTable )

  # Load history
  repList               <- .loadFit( ctlList$opMod$histFile )
  ctlList$opMod$histRpt <- c( repList$data, repList$repOpt ) 
  ctlList$opMod$posts   <- repList$posts
  ctlList$opMod$fYear   <- repList$fYear
  ctlList$opMod$species <- repList$species
  ctlList$opMod$stocks  <- repList$stocks
  ctlList$opMod$fleets  <- repList$gearLabs

  # Initialise the data structures to hold simulation states etc.
  simObj <- .initMS3pop(  ctlList )

  # simObj$ctlList$opMod$rep <-1
  simObj <- .condMS3pop_SISCA( simObj )

  return(simObj)

}

# # calculate reference points
simObj <- ms3Obj()

rp <- calcRefPts(simObj)
# Not working, need to look into these calcs....


# Generate equilibrium relationship reference point curves for each each area for:
# a) Yield Per Recruit for increasing fishing mortality (F), 
# b) Spawning biomass (SB) per recruit for increasing fishing mortality (F), 
# c) Yield for increasing fishing mortality (F), 
# d) Spawning biomass for increasing fishing mortality (F), 
# e) Recruitment for different levels of spawning biomass, and 
# f) Yield for different levels of spawning biomass

plotRefPts_p(repList = reports, pIdx=1)
