# <><><><><><><><><><><><><><><><><><><><><><><><><><><>
# stats.R
#
# Functions for calculating statistics from simulations
# produced by the ms3R closed loop simulation platform.
#
# Author: SDN Johnson
# Date: December 6, 2019
#
# Last Update: Dec 6, 2019
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><>


# calcBatchTable()
# Function to calculate a table of stats for 
# a whole batch (usually a given subfolder of
# ./Outputs/)
calcBatchTable <- function( batchFolder = "" )
{
  # set output folder
  simFolder <- here::here("Outputs",batchFolder)
  statsFolder <- file.path(simFolder, "statistics")

  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path=simFolder,full.names = FALSE,
                        recursive=FALSE)
  # Restrict to sim_ folders, count sims
  simList <- dirList[grep(pattern="sim",x=dirList)]
  nSims <- length(simList)

  statsTables <- lapply(  X = 1:nSims, FUN = makeStatTable,
                          folder = batchFolder )

  batchStatTable <- do.call("rbind", statsTables)

  # Now save to batchFolder
  if(!dir.exists(file.path(simFolder,"statistics")))
    dir.create(file.path(simFolder,"statistics"))

  outFile <- file.path(simFolder,"statistics","fullStatTable.csv")

  write.csv( batchStatTable, file = outFile)

  message(" (calcBatchTable) Stats table calculated and saved to ", outFile, "\n")

}


# makeStatTable()
# Wrapper for 
makeStatTable <- function( sims = 1, folder = "" )
{
  # First, load blob
  source("tools.R")

  .loadSim(sims, folder = folder)

  statTable <- .simPerfStats( obj = blob )

  statTable
} # END makeStatTable()

# .simPerfStats()
# Produces a statistics table for a simulation from
# the output produced by a runMS3() call
# inputs:   sim=int indicating which simulation to compute stats for
# outputs:  statTable=data.frame showing conservation and catch performance
# usage:    in lapply to produce stats for a group of simulations
.simPerfStats <- function( obj  )
{
  # Pull sublists
  om        <- obj$om
  opMod     <- obj$ctlList$opMod
  mp        <- obj$mp
  ctlList   <- obj$ctlList
  rp        <- obj$rp[[1]]
  
  # Control info
  nS      <- om$nS
  nP      <- om$nP
  nT      <- om$nT
  tMP     <- om$tMP
  pT      <- opMod$pT
  nReps   <- ctlList$ctl$nReps

  # Species and stock labels
  speciesNames  <- opMod$species
  stockNames    <- opMod$stocks
  fYear         <- opMod$fYear

  # get the replicate numbers for succesful fits (PD Hessians) in
  # MPs
  allConvReps <- obj$goodReps
  pdHess_itsp <- mp$assess$pdHess_itsp[allConvReps,,,]
  

  # Calculat probability of a good replicate (all PD hessians)
  goodReps  <- apply( X = pdHess_itsp, FUN = prod, MARGIN = 1, na.rm = T)
  goodReps  <- which(goodReps == 1)
  nGoodReps <- length(goodReps)
  pGoodReps <- nGoodReps/nReps

  # And the median over replicates of the probability of a 
  # PD Hessian over time
  probPDHess_isp <- apply( X = pdHess_itsp, FUN = mean, MARGIN = c(1,3,4), na.rm = T )
  medProbPDH_sp  <- apply( X = probPDHess_isp, FUN = median, MARGIN = c(2,3), na.rm = T )

  if(is.null(obj$simLabel))
  {
    simLabel <- stringr::str_split(obj$path,"/")[[1]]
    simLabel <- simLabel[length(simLabel)]
  } else simLabel <- obj$simLabel

  # First, create a data.frame of NAs with a row for each of MRE,MARE
  colLabels <- c( "simLabel",
                  "scenario","mp",
                  "species","stock",
                  "projObsErrMult",
                  "pGoodReps", "medProbPDH_sp",
                  "pBtGt.4Bmsy", "PBtGt.8Bmsy",
                  "pCtGtMSY", "pFtGtFmsy",
                  "avgCatch","AAV", "avgTACu" )

  statTable <- matrix( NA,  ncol = length(colLabels),
                            nrow = nS * nP )
  colnames(statTable) <- colLabels

  statTable <- as.data.frame(statTable)

  # Fill in global values
  statTable[,"simLabel"]        <- simLabel
  statTable[,"scenario"]        <- ctlList$ctl$scenarioName
  statTable[,"mp"]              <- ctlList$ctl$mpName
  statTable[,"pGoodReps"]       <- pGoodReps
  statTable[,"projObsErrMult"]  <- opMod$projObsErrMult

  # Need to start layering in performance metrics


  # Pull reference points
  B0_sp     <- opMod$histRpt$B0_sp
  Bmsy_sp   <- rp$FmsyRefPts$BeqFmsy_sp
  Fmsy_sp   <- rp$FmsyRefPts$Fmsy_sp
  Umsy_sp   <- rp$FmsyRefPts$Umsy_sp
  MSY_sp    <- rp$FmsyRefPts$YeqFmsy_sp

  # Calculate depletion wrt Bmsy
  SB_ispt   <- om$SB_ispt[allConvReps,,,,drop = FALSE]
  C_ispt    <- om$C_ispt[allConvReps,,,,drop = FALSE]
  TAC_ispt  <- mp$hcr$TAC_ispt[allConvReps,,,,drop = FALSE]
  TACu_ispt <- C_ispt / TAC_ispt
  F_ispt    <- om$F_ispft[allConvReps,,,2,]
  
  # Calculate probability that Bt above LRP
  pBtGt.4Bmsy_sp <- .calcStatsProportion( TS_ispt = SB_ispt,
                                          ref_sp = Bmsy_sp,
                                          tdx = tMP:nT,
                                          prop = .4,
                                          nS = nS,
                                          nP = nP )
  # In a healthy state (Bt > .8Bmsy)
  pBtGt.8Bmsy_sp <- .calcStatsProportion( TS_ispt = SB_ispt,
                                          ref_sp = Bmsy_sp,
                                          tdx = tMP:nT,
                                          prop = .8,
                                          nS = nS,
                                          nP = nP )

  # Overfishing is occuring (Ft > Fmsy)
  pFtGtFmsy_sp <- .calcStatsProportion( TS_ispt = F_ispt,
                                        ref_sp = Fmsy_sp,
                                        tdx = tMP:nT,
                                        prop = 1,
                                        nS = nS,
                                        nP = nP )

  # Overfishing is occuring (Ct > MSY)
  pCtGtMSY_sp <- .calcStatsProportion(  TS_ispt = C_ispt,
                                        ref_sp = MSY_sp,
                                        tdx = tMP:nT,
                                        prop = 1,
                                        nS = nS,
                                        nP = nP )


  # Average catch and TAC utilisation
  Cbar_sp     <- apply( X = C_ispt[,,,tMP:nT], FUN = mean, MARGIN = c(2,3), na.rm = T)
  TACubar_sp  <- apply( X = TACu_ispt[,,,tMP:nT], FUN = mean, MARGIN = c(2,3), na.rm = T)

  # Catch variability
  AAV_sp      <- .calcStatsAAV( C_ispt = C_ispt,
                                tdx = tMP:nT,
                                qProbs = c(0.5),
                                margin = c(2,3) )

  for( s in 1:nS )
    for(p in 1:nP )
    {
      

      rowIdx <- (s-1) * nP + p

      statTable[rowIdx, "species"]        <- speciesNames[s]
      statTable[rowIdx, "stock"]          <- stockNames[p]
      statTable[rowIdx, "medProbPDH_sp"]  <- medProbPDH_sp[s,p]

      # Conservation performance
      statTable[rowIdx,"pBtGt.4Bmsy"]     <- pBtGt.4Bmsy_sp[s,p]
      statTable[rowIdx,"PBtGt.8Bmsy"]     <- pBtGt.8Bmsy_sp[s,p]
      statTable[rowIdx,"pCtGtMSY"]        <- pCtGtMSY_sp[s,p]
      statTable[rowIdx,"pFtGtFmsy"]       <- pFtGtFmsy_sp[s,p]

      # Catch statistics
      statTable[rowIdx,"avgCatch"]        <- Cbar_sp[s,p]
      statTable[rowIdx,"AAV"]             <- AAV_sp[s,p]
      statTable[rowIdx,"avgTACu"]         <- TACubar_sp[s,p]

    }


  statTable
} # END .simPerfStats()

# .calcStatsProportion()
# Calculates the probability of a time series
# being above a certain proportion of a reference
# level. Used for depletion and overfishing statistics.
.calcStatsProportion <- function( TS_ispt = SB_ispt,
                                  ref_sp = Bmsy_sp,
                                  tdx = tMP:nT,
                                  prop = .4,
                                  nS = 3,
                                  nP = 3 )
{

  # Reduce to time period
  TS_ispt <- TS_ispt[,,,tdx]

  Quotient_ispt <- array(NA, dim = dim(TS_ispt) )

  # First, take quotient by reference
  for(s in 1:nS)
    for(p in 1:nP )
      Quotient_ispt[,s,p,] <- TS_ispt[,s,p,] / ref_sp[s,p] 

  aboveIdx <- which( Quotient_ispt > prop)

  # Set an indicator array
  Ind_ispt <- array(0, dim = dim(Quotient_ispt))
  Ind_ispt[aboveIdx] <- 1

  probGtDep_sp <- apply( X = Ind_ispt, FUN = mean, MARGIN = c(2,3) )
      
  return( probGtDep_sp )
} # END .calcStatsProportion

# .calcStatsAAV()
# Calculate average annual variation
# in catch, measure of the variability
# in landings for each stock/species.
.calcStatsAAV <- function(  C_ispt = C_ispt,
                            tdx = tMP:nT,
                            qProbs = c(0.025,0.5,0.975),
                            margin = c(2,3) )
{
  C_ispt <- C_ispt[,,,tdx]
  # Append margin
  marg <- c(1,margin)
  # Add over margins, in case we are looking
  # complex/species/stock aggregates
  sumC_ispt <- apply( X = C_ispt, FUN = sum, MARGIN = c(marg,4), na.rm = T )


  # Calculate diff
  diffC_ispt        <- aperm(apply( X = sumC_ispt, FUN = diff, MARGIN = marg, na.rm = T ),c(2:4,1))
  absDiffC_ispt     <- abs(diffC_ispt)
  sumAbsDiffC_isp   <- apply( X = absDiffC_ispt, FUN = sum, MARGIN = marg, na.rm = T )
  sumCatch_isp      <- apply( X = C_ispt, FUN = sum, MARGIN = marg, na.rm = T )

  AAV_isp           <- sumAbsDiffC_isp / sumCatch_isp
  AAV_isp[!is.finite(AAV_isp)] <- 0

  AAV_qsp   <- apply(X = AAV_isp, FUN = quantile, probs = qProbs, MARGIN = marg[-1] )

  return(AAV_qsp)
} # END .calcStatsAAV
