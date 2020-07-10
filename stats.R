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


# calcLoss()
# Loads blob and calculates yearly catch and biomass 
# loss values compared to a nominated baseline sim.
# By default, calculates biomass and catch loss
# on relative and absolute scale
calcLoss <- function( sim         = 1,
                      baseline    = "sim_parBatomniRuns_Long3",
                      groupFolder = "diffCV_fixedF_longGrid",
                      lossVars    = c("C_ispt","SB_ispt"),
                      output      = TRUE )
{
  # First, load the sim that we want
  # to calculate loss for
  simFolder <- .loadSim( sim = sim, folder = groupFolder )

  # Save the blob
  lossSim <- blob

  # figure out which reps we want
  goodReps  <- blob$goodReps
  nReps     <- sum(goodReps)
  totReps   <- dim(lossSim$om$SB_ispt)[1]

  # Now load the baseline (omniscient manager)
  .loadSim( sim = baseline, folder = groupFolder )

  baseSim <- blob

  # Model dimensions
  tMP <- lossSim$om$tMP
  nT  <- lossSim$om$nT
  nF  <- lossSim$om$nF
  nS  <- lossSim$om$nS
  nP  <- lossSim$om$nP
  pT  <- lossSim$om$pT

  speciesNames  <- c(lossSim$om$speciesNames,"dataPooled")
  stockNames    <- c(lossSim$om$stockNames,"coastWide")
  fleetNames    <- lossSim$om$fleetNames

  # We want to calculate loss of given OM states
  # Make a list to hold baseline
  baseLineStates        <- vector(  mode = "list", 
                                    length = length(lossVars))
  names(baseLineStates) <- lossVars
  # copy for simulation states and loss
  simStates       <- baseLineStates
  lossRaw         <- baseLineStates
  lossRel         <- baseLineStates


  # Pull control list settings for
  # coastwide and data pooled

  lossCtl         <- lossSim$ctlList


  # Loop over 
  for( varIdx in 1:length(lossVars) )
  {
    # Clean arrays for base and sim states
    baseState_ispt <- array(NA, dim = c(nReps,nS+1,nP+1,nT),
                                dimnames = list(  rep = which(goodReps),
                                                  species = speciesNames,
                                                  stock   = stockNames,
                                                  tdx     = 1:nT ) )
    simState_ispt  <- baseState_ispt

    # Get variable
    var <- lossVars[varIdx]

    # Copy base and sim for individual stocks/species
    simReps <- (1:totReps)[goodReps]
    baseState_ispt[,1:nS,1:nP,] <- baseSim$om[[var]][simReps,,,1:nT]
    simState_ispt[,1:nS,1:nP,]  <- lossSim$om[[var]][simReps,,,1:nT]

    # Now do the aggregates
    # totally aggregated
    baseState_ispt[,nS+1,nP+1,] <- apply(   X = baseSim$om[[var]][simReps,,,1:nT],
                                            FUN = sum,
                                            MARGIN = c(1,4) )

    simState_ispt[,nS+1,nP+1,]  <- apply(   X = lossSim$om[[var]][simReps,,,1:nT],
                                            FUN = sum,
                                            MARGIN = c(1,4) )

    # data-pooled
    baseState_ispt[,nS+1,1:nP,]     <- apply(   X = baseSim$om[[var]][simReps,,,1:nT],
                                            FUN = sum,
                                            MARGIN = c(1,3,4) )

    simState_ispt[,nS+1,1:nP,]      <- apply(   X = lossSim$om[[var]][simReps,,,1:nT],
                                            FUN = sum,
                                            MARGIN = c(1,3,4) )
    
    # coastwide
    baseState_ispt[,1:nS,nP+1,]     <- apply(   X = baseSim$om[[var]][simReps,,,1:nT],
                                            FUN = sum,
                                            MARGIN = c(1,2,4) )

    simState_ispt[,1:nS,nP+1,]      <- apply(   X = lossSim$om[[var]][simReps,,,1:nT],
                                            FUN = sum,
                                            MARGIN = c(1,2,4) )

    baseLineStates[[var]] <- baseState_ispt
    simStates[[var]]      <- simState_ispt

    lossRaw[[var]]        <- baseLineStates[[var]] - simStates[[var]]
    lossRel[[var]]        <- lossRaw[[var]] / baseLineStates[[var]]

   
  }

  outList <- list(  simID         = lossSim$folder,
                    speciesNames  = speciesNames,
                    stockNames    = stockNames,
                    fYear         = lossSim$ctlList$opMod$fYear,
                    simStates     = simStates,
                    baseStates    = baseLineStates,
                    lossRaw       = lossRaw,
                    lossRel       = lossRel,
                    tMP           = tMP,
                    nT            = nT, 
                    nF            = nF, 
                    nS            = nS, 
                    nP            = nP, 
                    pT            = pT,
                    obsCVmult     = lossSim$ctlList$opMod$projObsErrMult,
                    simFolder     = simFolder  )


  # Save loss to sim folder
  save( outList, file = file.path(simFolder,"loss.RData") )

  if(output)
    return(outList)
} # END calcLoss()



makeLossTable <- function(  sim = 1, 
                            groupFolder = "shortGrid",
                            save = TRUE )
{
  # Load loss reports
  objLoss <- .loadLoss(sim = sim, groupFolder)

  # Species/stock names and model dims
  speciesNames    <- objLoss$speciesNames
  stockNames      <- objLoss$stockNames
  tMP             <- objLoss$tMP
  nT              <- objLoss$nT 
  nF              <- objLoss$nF 
  nS              <- objLoss$nS 
  nP              <- objLoss$nP 
  pT              <- objLoss$pT
  simFolder       <- objLoss$simFolder
  obsCVmult       <- objLoss$obsCVmult




  # Here's the loss table
  lossTableColNames <- c( "sim",
                          "scenario",
                          "mp",
                          "species",
                          "stock",
                          "obsCVmult",
                          "relLossCat_short",
                          "relLossCat_med",
                          "relLossCat_long",
                          "absLossCat_short",
                          "absLossCat_med",
                          "absLossCat_long",
                          "relLossBio_short",
                          "relLossBio_med",
                          "relLossBio_long",
                          "absLossBio_short",
                          "absLossBio_med",
                          "absLossBio_long" )

  lossTable <- matrix(NA, nrow = (nS * 1) * (nP + 1), ncol = length(lossTableColNames) )
  colnames(lossTable) <- lossTableColNames
  lossTable <- as.data.frame(lossTable)

  lossTable$sim       <- objLoss$sim
  lossTable$scenario  <- objLoss$scenario
  lossTable$mp        <- objLoss$mp
  lossTable$obsCVmult <- obsCVmult

  
  lossCat_shortTerm  <- calcTotalLossPeriod(objLoss,"C_ispt", period = tMP:(tMP+10))
  lossCat_medTerm    <- calcTotalLossPeriod(objLoss,"C_ispt", period = tMP:(tMP+20))
  lossCat_longTerm   <- calcTotalLossPeriod(objLoss,"C_ispt", period = tMP:nT)
  lossBio_shortTerm  <- calcTotalLossPeriod(objLoss,"SB_ispt", period = tMP:(tMP+10))
  lossBio_medTerm    <- calcTotalLossPeriod(objLoss,"SB_ispt", period = tMP:(tMP+20))
  lossBio_longTerm   <- calcTotalLossPeriod(objLoss,"SB_ispt", period = tMP:nT)

  for( s in 1:(nS+1) )
    for( p in 1:(nP+1) )
    {
      tabRow <- (s - 1)*( nP + 1 ) + p

      lossTable[tabRow,"species"] <- speciesNames[s]
      lossTable[tabRow,"stock"]   <- stockNames[p]
      lossTable[tabRow,"relLossCat_short"]    <- median(lossCat_shortTerm$totRelLoss_isp[,s,p])
      lossTable[tabRow,"relLossCat_med"]      <- median(lossCat_medTerm$totRelLoss_isp[,s,p])
      lossTable[tabRow,"relLossCat_long"]     <- median(lossCat_longTerm$totRelLoss_isp[,s,p])
      lossTable[tabRow,"absLossCat_short"]    <- median(lossCat_shortTerm$totAbsLoss_isp[,s,p])
      lossTable[tabRow,"absLossCat_med"]      <- median(lossCat_medTerm$totAbsLoss_isp[,s,p])
      lossTable[tabRow,"absLossCat_long"]     <- median(lossCat_longTerm$totAbsLoss_isp[,s,p])
      lossTable[tabRow,"relLossBio_short"]    <- median(lossBio_shortTerm$totRelLoss_isp[,s,p])
      lossTable[tabRow,"relLossBio_med"]      <- median(lossBio_medTerm$totRelLoss_isp[,s,p])
      lossTable[tabRow,"relLossBio_long"]     <- median(lossBio_longTerm$totRelLoss_isp[,s,p])
      lossTable[tabRow,"absLossBio_short"]    <- median(lossBio_shortTerm$totAbsLoss_isp[,s,p])
      lossTable[tabRow,"absLossBio_med"]      <- median(lossBio_medTerm$totAbsLoss_isp[,s,p])
      lossTable[tabRow,"absLossBio_long"]     <- median(lossBio_longTerm$totAbsLoss_isp[,s,p])
    
    }

  write.csv( lossTable, file = file.path(simFolder,"lossTable.csv") )

  lossTable
} # END makeLossTable()

# calcTotalLossPeriod()
# Takes a loss object, and variable name,
# and calculates total abs and relative loss
# for a given time period
calcTotalLossPeriod <- function(  obj, 
                                  var = "C_ispt",
                                  period = tMP:nT )
{
  # Pull abs and rel loss
  relLoss_ispt  <- obj$lossRel[[var]][,,,period]
  rawLoss_ispt  <- obj$lossRaw[[var]][,,,period]

  totRelLoss_isp    <- apply( X = abs(relLoss_ispt), FUN = sum, MARGIN = c(1,2,3) )
  totAbsLoss_isp    <- apply( X = abs(rawLoss_ispt), FUN = sum, MARGIN = c(1,2,3) )

  out <- list(  totRelLoss_isp = totRelLoss_isp,
                totAbsLoss_isp = totAbsLoss_isp )

  return(out)
} # END calcTotalLossPeriod



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
.simPerfStats <- function( obj, inclAgg=FALSE, exclNoTAC=FALSE  )
{
  # Pull sublists
  om        <- obj$om
  opMod     <- obj$ctlList$opMod
  mp        <- obj$mp
  ctlList   <- obj$ctlList
  
  # Control info
  nS      <- om$nS
  nP      <- om$nP
  nT      <- om$nT
  nF      <- om$nF
  tMP     <- om$tMP
  pT      <- opMod$pT
  nReps   <- ctlList$ctl$nReps

  # Species and stock labels

  speciesNames  <- opMod$species
  stockNames    <- opMod$stocks
  fYear         <- opMod$fYear

  speciesNames  <- "HG Herring"
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]


  # get the replicate numbers for succesful fits (PD Hessians) in
  # MPs
  allConvReps <- obj$goodReps
  nGoodReps   <- sum(allConvReps)
  pdHess_itsp <- array(NA, dim = c(nGoodReps,pT,nS,nP))
  # Fill
  pdHess_itsp[1:nGoodReps,,,] <- mp$assess$pdHess_itsp[allConvReps,,,]
  

  # Calculat probability of a good replicate (all PD hessians)
  nGoodReps <- sum(allConvReps)
  pGoodReps <- signif(nGoodReps/obj$nSim,2)

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
                  "scenario",
                  "mp",
                  "species",
                  "stock",
                  "projObsErrMult",
                  "nGoodReps", 
                  "medProbPDH_sp",
                  "pGoodReps",
                  "pBtGt.3B0", 
                  # "PBtGt.8Bmsy",
                  # "pBtGtBmsy",
                  # "pCtGtMSY", 
                  # "pFtGtFmsy",
                  "nYrsSOK",
                  "avgLicSOK",
                  "avgPonded_t",
                  "avgSOK_t",
                  "pondAAV",
                  "nYrsComm",
                  "avgCatch_t",
                  "catchAAV"
                  )

  statTable <- matrix( NA,  ncol = length(colLabels),
                            nrow = nS * nP )
  colnames(statTable) <- colLabels

  statTable <- as.data.frame(statTable)

  # Fill in global values
  statTable[,"simLabel"]        <- simLabel
  statTable[,"scenario"]        <- ctlList$ctl$scenarioName
  statTable[,"mp"]              <- ctlList$ctl$mpName
  statTable[,"projObsErrMult"]  <- opMod$projObsErrMult
  statTable[,"nGoodReps"]       <- nGoodReps
  statTable[,"pGoodReps"]       <- pGoodReps
  
  # Need to start layering in performance metrics


  # Pull reference points
  B0_isp <- array(NA, dim=c(nGoodReps,nS,nP))
  for(i in 1:nGoodReps)
  {
    rp          <- obj$rp[[i]]
    B0_isp[i,,] <- rp$B0_sp

  } 
  
  # Bmsy_sp   <- rp$FmsyRefPts$BeqFmsy_sp
  # Fmsy_sp   <- rp$FmsyRefPts$Fmsy_sp
  # Umsy_sp   <- rp$FmsyRefPts$Umsy_sp
  # MSY_sp    <- rp$FmsyRefPts$YeqFmsy_sp

  if( nGoodReps > 0 )
  { 

    # Spawning biomass and fishing mortality
    SB_ispt   <- om$SB_ispt[allConvReps,,,,drop = FALSE]
    F_ispt    <- array(NA, dim = c(nGoodReps,nS,nP,nT))
    F_ispt[,1:nS,,] <- om$F_ispft[allConvReps,,,2,]

    # Catch stats from fisheries
    C_ispft <- om$C_ispft[allConvReps,,,,,drop = FALSE]
    C_ispt  <- apply(C_ispft[,,,1:3,,drop=FALSE], FUN=sum, MARGIN=c(1,2,3,5))

    # Ponded fish
    P_ispft      <- om$P_ispft[allConvReps,,,,,drop = FALSE]
    P_ispt       <- apply(C_ispft[,,,6:7,,drop=FALSE], FUN=sum, MARGIN=c(1,2,3,5))
    
    # SOK licenses
    L_ispft <- mp$hcr$sokEff_ispft [allConvReps,,,,,drop = FALSE]
    L_ispt  <- apply( X = L_ispft[,,,6:7,,drop = FALSE], FUN = sum, MARGIN = c(1,2,3,5))


    # SOK product in kt
    psi_ispft    <- om$psi_ispft[allConvReps,,,,,drop = FALSE]
    SOK_ispft    <- P_ispft*psi_ispft
    SOK_ispt     <- apply(SOK_ispft[,,,6:7,,drop=FALSE], FUN=sum, MARGIN=c(1,2,3,5))
    
    # Calculate probability that Bt above LRP
    pBtGt.3B0_sp <- .calcStatsProportion(   TS_ispt = SB_ispt,
                                            ref_isp = B0_isp,
                                            tdx = tMP:nT,
                                            prop = .3,
                                            nS = nS,
                                            nP = nP )

    # # In a healthy state (Bt > .8Bmsy)
    # pBtGt.8Bmsy_sp <- .calcStatsProportion( TS_ispt = SB_ispt,
    #                                         ref_sp = Bmsy_sp,
    #                                         tdx = tMP:nT,
    #                                         prop = .8,
    #                                         nS = nS,
    #                                         nP = nP )

    # # Above Bmsy
    # pBtGtBmsy_sp <- .calcStatsProportion( TS_ispt = SB_ispt,
    #                                       ref_sp = Bmsy_sp,
    #                                       tdx = tMP:nT,
    #                                       prop = 1.0,
    #                                       nS = nS,
    #                                       nP = nP )  

    # # Overfishing is occuring (Ft > Fmsy)
    # pFtGtFmsy_sp <- .calcStatsProportion( TS_ispt = F_ispt,
    #                                       ref_sp = Fmsy_sp,
    #                                       tdx = tMP:nT,
    #                                       prop = 1,
    #                                       nS = nS,
    #                                       nP = nP )

    # # Overfishing is occuring (Ct > MSY)
    # pCtGtMSY_sp <- .calcStatsProportion(  TS_ispt = C_ispt,
    #                                       ref_sp = MSY_sp,
    #                                       tdx = tMP:nT,
    #                                       prop = 1,
    #                                       nS = nS,
    #                                       nP = nP )

    # Calculate number of years that fisheries are open
    sokYrs_ispt <- P_ispt 
    commYrs_ispt <- C_ispt
    sokYrs_ispt[sokYrs_ispt>0]   <-1
    commYrs_ispt[commYrs_ispt>0] <-1

    sokYrsBar_spt  <- apply( X = sokYrs_ispt[,,,,drop = FALSE], FUN = mean, MARGIN = c(2,3,4))
    sokYrsBar_sp   <- apply( X = sokYrsBar_spt[,,tMP:nT,drop = FALSE], FUN = sum, MARGIN = c(1,2))

    commYrsBar_spt  <- apply( X = commYrs_ispt[,,,,drop = FALSE], FUN = mean, MARGIN = c(2,3,4))
    commYrsBar_sp   <- apply( X = commYrsBar_spt[,,tMP:nT,drop = FALSE], FUN = sum, MARGIN = c(1,2))


    # Average fisheries and SOK stats for years when fisheries are open
    if(exclNoTAC)
    {
      # Assign NAs for years with zero catch to remove from stats
      if( any(C_ispt[,,,tMP:nT] >0) )
        C_ispt[C_ispt==0] <- NA

      if( any(P_ispt[,,,tMP:nT] >0) )
        P_ispt[P_ispt==0] <- NA

      if( any(SOK_ispt[,,,tMP:nT] >0) )
        SOK_ispt[SOK_ispt==0] <- NA
      
      if( any(L_ispt[,,,tMP:nT] >0) )
        L_ispt[L_ispt==0] <- NA

      # calculate average stats across reps
      Cbar_isp      <- apply( X = C_ispt[,,,tMP:nT,drop = FALSE], FUN = mean, MARGIN = c(1,2,3), na.rm = T)
      Pbar_isp      <- apply( X = P_ispt[,,,tMP:nT,drop = FALSE], FUN = mean, MARGIN = c(1,2,3), na.rm = T)
      SOKbar_isp    <- apply( X = SOK_ispt[,,,tMP:nT,drop = FALSE], FUN = mean, MARGIN = c(1,2,3), na.rm = T)
      LBar_isp      <- apply( X = L_ispt[,,,tMP:nT,drop = FALSE], FUN = mean, MARGIN = c(1,2,3), na.rm = T)

      # Now assign zeros for any reps that have avg. catch for full time series
      Cbar_isp[is.na(Cbar_isp)] <- 0
      Pbar_isp[is.na(Pbar_isp)] <- 0
      SOKbar_isp[is.na(SOKbar_isp)] <- 0
      LBar_isp [is.na(LBar_isp )] <- 0

      # catch, live ponded biomass, and SOK product in metric tonnes
      Cbar_sp     <- apply( X = Cbar_isp[,,,drop = FALSE], FUN = mean, MARGIN = c(2,3), na.rm = T)
      Pbar_sp     <- apply( X = Pbar_isp[,,,drop = FALSE], FUN = mean, MARGIN = c(2,3), na.rm = T)
      SOKbar_sp   <- apply( X = SOKbar_isp[,,,drop = FALSE], FUN = mean, MARGIN = c(2,3), na.rm = T)

      # SOK licenses
      Lbar_sp   <- apply( X = LBar_isp[,,,drop = FALSE], FUN = mean, MARGIN = c(2,3), na.rm = T)


      # Catch and ponding variability
      catchAAV_sp      <- .calcStatsAAV_sp( C_ispt = C_ispt,
                                  tdx = tMP:nT,
                                  qProbs = c(0.5),
                                  margin = c(2,3) )

      pondAAV_sp      <- .calcStatsAAV_sp( C_ispt = P_ispt,
                                  tdx = tMP:nT,
                                  qProbs = c(0.5),
                                  margin = c(2,3) )  

    }  



    # Average fisheries and SOK stats for all years
    if(!exclNoTAC)
    {
      # catch, live ponded biomass, and SOK product in metric tonnes
      Cbar_sp     <- apply( X = C_ispt[,,,tMP:nT,drop = FALSE], FUN = mean, MARGIN = c(2,3), na.rm = T)
      Pbar_sp     <- apply( X = P_ispt[,,,tMP:nT,drop = FALSE], FUN = mean, MARGIN = c(2,3), na.rm = T)
      SOKbar_sp   <- apply( X = SOK_ispt[,,,tMP:nT,drop = FALSE], FUN = mean, MARGIN = c(2,3), na.rm = T)

      # SOK licenses for open (o), closed (c) and all ponds (co)
      Lbar_sp <- apply( X = L_ispt[,,,tMP:nT,drop = FALSE], FUN = mean, MARGIN = c(2,3), na.rm = T) # open & closed  
      # cSOKbar_sp  <- apply( X = sokEff_ispft[,,,6,tMP:nT,drop = FALSE], FUN = mean, MARGIN = c(2,3), na.rm = T)
      # oSOKbar_sp  <- apply( X = sokEff_ispft[,,,7,tMP:nT,drop = FALSE], FUN = mean, MARGIN = c(2,3), na.rm = T)
      

      # Catch and ponding variability
      catchAAV_sp      <- .calcStatsAAV_sp( C_ispt = C_ispt,
                                  tdx = tMP:nT,
                                  qProbs = c(0.5),
                                  margin = c(2,3) )

      pondAAV_sp      <- .calcStatsAAV_sp( C_ispt = P_ispt,
                                  tdx = tMP:nT,
                                  qProbs = c(0.5),
                                  margin = c(2,3) )  
    }  





    for( s in 1:nS )
      for(p in 1:nP )
      {
        

        rowIdx <- (s-1) * nP + p

        statTable[rowIdx, "species"]        <- speciesNames[s]
        statTable[rowIdx, "stock"]          <- stockNames[p]
        statTable[rowIdx, "medProbPDH_sp"]  <- medProbPDH_sp[s,p]

        # Conservation performance
        statTable[rowIdx,"pBtGt.3B0"]     <- round(pBtGt.3B0_sp[s,p],2)
        # statTable[rowIdx,"pBtGt.4Bmsy"]     <- round(pBtGt.4Bmsy_sp[s,p],2)
        # statTable[rowIdx,"PBtGt.8Bmsy"]     <- round(pBtGt.8Bmsy_sp[s,p],2)
        # statTable[rowIdx,"pBtGtBmsy"]       <- round(pBtGtBmsy_sp[s,p],2)
        # statTable[rowIdx,"pCtGtMSY"]        <- round(pCtGtMSY_sp[s,p],2)
        # statTable[rowIdx,"pFtGtFmsy"]       <- round(pFtGtFmsy_sp[s,p],2)

        # Catch statistics
        statTable[rowIdx,"nYrsSOK"]         <- round(sokYrsBar_sp[s,p],1)
        statTable[rowIdx,"nYrsComm"]        <- round(commYrsBar_sp[s,p],1)

        statTable[rowIdx,"avgCatch_t"]      <- round(Cbar_sp[s,p]*1e3,1)
        statTable[rowIdx,"catchAAV"]        <- round(catchAAV_sp[s,p]*1e2,1)
        statTable[rowIdx,"avgPonded_t"]     <- round(Pbar_sp[s,p]*1e3,1)
        statTable[rowIdx,"pondAAV"]         <- round(pondAAV_sp[s,p]*1e2,1)
        statTable[rowIdx,"avgSOK_t"]        <- round(SOKbar_sp[s,p]*1e3,1)
        statTable[rowIdx,"avgLicSOK"]       <- round(Lbar_sp[s,p],1)

        # statTable[rowIdx,"avg_cSOKlic"]     <- round(cSOKbar_sp[s,p],1)
        # statTable[rowIdx,"avg_oSOKlic"]     <- round(oSOKbar_sp[s,p],1)

      }


    # Add row for aggregated stats across all species and stocks
    if(inclAgg)
    {
      aggRow <- statTable[1,]

      # Aggregate SOK licenses
      L_ispft <- mp$hcr$sokEff_ispft [allConvReps,,,,,drop = FALSE]

      
      # Calculate probability that SBt above LRP
      SB_it <- apply(SB_ispt, MARGIN=c(1,4), FUN=sum, drop=FALSE)
      B0_i     <- apply(B0_isp, MARGIN=1, FUN=sum)


      pBtGt.3B0_agg <- .calcStatsProportionAgg(   TS_it = SB_it,
                                                  ref_i = B0_i,
                                                  tdx = tMP:nT,
                                                  prop = .3,
                                                  nS = 1,
                                                  nP = 1 )

      # number of years with SOK or commercial fishing at aggregate level
      sokYrs_it                 <- apply(P_ispt, FUN=sum, MARGIN=c(1,4), na.rm=T)
      sokYrs_it[sokYrs_it>0]    <-1
      sokYrs_i                  <- apply( X = sokYrs_it[,tMP:nT], FUN = sum, MARGIN=1)
      sokYrsBar_agg             <- mean(sokYrs_i)

      commYrs_it                <- apply(C_ispt, FUN=sum, MARGIN=c(1,4), na.rm=T)
      commYrs_it[commYrs_it>0]  <-1
      commYrs_i                 <- apply( X = commYrs_it[,tMP:nT], FUN = sum, MARGIN=1)
      commYrsBar_agg            <- mean(commYrs_i)

      # Aggregate fisheries catch, SOK ponded biomass & SOK licenses then calcuate average across iReps and projT
      C_it     <- apply( X = C_ispt[,,,,drop = FALSE], FUN = sum, MARGIN = c(1,4), na.rm = T)
      P_it     <- apply( X = P_ispt[,,,,drop = FALSE], FUN = sum, MARGIN = c(1,4), na.rm = T)
      SOK_it   <- apply( X = SOK_ispt[,,,,drop = FALSE], FUN = sum, MARGIN = c(1,4), na.rm = T)
      L_it     <- apply( X = L_ispt[,,,,drop = FALSE], FUN = sum, MARGIN = c(1,4), na.rm = T)    

      # # Average fisheries and SOK stats for years when fisheries are open
      if(exclNoTAC)
      {
        # Assign NAs for years with zero catch to remove from stats
        if( any(C_it[,tMP:nT] >0) )
          C_it[C_it==0] <- NA

        if( any(P_it[,tMP:nT] >0) )
          P_it[P_it==0] <- NA

        if( any(SOK_it[,tMP:nT] >0) )
          SOK_it[SOK_it==0] <- NA

        if( any(L_it[,tMP:nT] >0) )
          L_it[L_it==0] <- NA

        # calculate average stats across reps
        Cbar_i      <- apply( X = C_it[,tMP:nT,drop = FALSE], FUN = mean, MARGIN = c(1), na.rm = T)
        Pbar_i      <- apply( X = P_it[,tMP:nT,drop = FALSE], FUN = mean, MARGIN = c(1), na.rm = T)
        SOKbar_i    <- apply( X = SOK_it[,tMP:nT,drop = FALSE], FUN = mean, MARGIN = c(1), na.rm = T)
        Lbar_i      <- apply( X = L_it[,tMP:nT,drop = FALSE], FUN = mean, MARGIN = c(1), na.rm = T)

        # Now assign zeros for any reps that have avg. catch of zero for full time series
        Cbar_i[is.na(Cbar_i)] <- 0
        Pbar_i[is.na(Pbar_i)] <- 0
        SOKbar_i[is.na(SOKbar_i)] <- 0
        Lbar_i[is.na(Lbar_i )] <- 0

        # catch, live ponded biomass, and SOK product in metric tonnes
        Cbar_agg    <- mean(Cbar_i)
        Pbar_agg    <- mean(Pbar_i)
        SOKbar_agg  <- mean(SOKbar_i)

        # SOK licenses
        Lbar_agg   <- mean(Lbar_i)

      }  

      if(!exclNoTAC)
      {
        # catch, live ponded biomass, and SOK product in metric tonnes
        Cbar_agg     <- mean(C_it[,tMP:nT], na.rm = T)
        Pbar_agg     <- mean(P_it[,tMP:nT], na.rm = T)
        SOKbar_agg   <- mean(SOK_it[,tMP:nT], na.rm = T)

        # SOK licenses
        Lbar_agg     <- mean(L_it[,tMP:nT], na.rm = T)
      }  
      
      # oSOKbar_agg  <- mean(L_it[,7,tMP:nT], na.rm = T)
      # coSOKbar_agg <- mean(L_it[,7,tMP:nT], na.rm = T)

      # Catch and ponding variability
      catchAAV_agg      <- .calcStatsAAV_agg(C_it = C_it,
                                            tdx = tMP:nT,
                                            qProbs = c(0.5))

      pondAAV_agg      <- .calcStatsAAV_agg(C_it = P_it,
                                            tdx = tMP:nT,
                                            qProbs = c(0.5))
      
      # Fill aggregate row
      aggRow[,"stock"]         <- 'aggregate'
      
      # Conservation
      aggRow[,"pBtGt.3B0"]     <- round(pBtGt.3B0_agg,2)
      
      # commercial catch 
      aggRow[,"nYrsComm"]      <- round(commYrsBar_agg,1)
      aggRow[,"avgCatch_t"]    <- round(Cbar_agg*1e3,1)
      aggRow[,"catchAAV"]      <- round(catchAAV_agg*1e2,2)
      
      # ponded fish
      aggRow[,"nYrsSOK"]       <- round(sokYrsBar_agg,1)
      aggRow[,"avgPonded_t"]   <- round(Pbar_agg*1e3,1)
      aggRow[,"pondAAV"]       <- round(pondAAV_agg*1e2,1)

      # SOK product and licenses
      aggRow[,"avgSOK_t"]      <- round(SOKbar_agg*1e3,1)
      aggRow[,"avgLicSOK"]     <- round(Lbar_agg,1)
      # aggRow[,"avg_cSOKlic"]   <- round(cSOKbar_agg,1)
      # aggRow[,"avg_oSOKlic"]   <- round(oSOKbar_agg,1)


      statTable <- rbind(statTable, aggRow)

    }  
  }

  statTable
} # END .simPerfStats()

# .calcStatsProportionAgg()
# Calculates the probability of a time series
# being above a certain proportion of a reference
# level. Used for depletion and overfishing statistics.
.calcStatsProportionAgg <- function( TS_it = SB_it,
                                  ref_i = B0_i,
                                  tdx = tMP:nT,
                                  prop = .4,
                                  nS = 1,
                                  nP = 1 )
{

  # Reduce to time period
  TS_it <- TS_it[,tdx]

  Quotient_it <- array(NA, dim = dim(TS_it) )

  # First, take quotient by reference
  Quotient_it <- TS_it / ref_i 

  # Set an indicator array
  Ind_it <- Quotient_it > prop

  probGtDep_sp <- mean(Ind_it)
      
  return( probGtDep_sp )
} # END .calcStatsProportionAgg


# .calcStatsProportion()
# Calculates the probability of a time series
# being above a certain proportion of a reference
# level. Used for depletion and overfishing statistics.
.calcStatsProportion <- function( TS_ispt = SB_ispt,
                                  ref_isp = Bmsy_sp,
                                  tdx = tMP:nT,
                                  prop = .4,
                                  nS = 3,
                                  nP = 3 )
{

  # Reduce to time period
  TS_ispt <- TS_ispt[,,,tdx,drop = FALSE]

  Quotient_ispt <- array(NA, dim = dim(TS_ispt) )

  # First, take quotient by reference
  for(s in 1:nS)
    for(p in 1:nP )
      Quotient_ispt[,s,p,] <- TS_ispt[,s,p,] / ref_isp[,s,p] 

  # Set an indicator array
  Ind_ispt <- Quotient_ispt > prop

  probGtDep_sp <- apply( X = Ind_ispt, FUN = mean, MARGIN = c(2,3) )
      
  return( probGtDep_sp )
} # END .calcStatsProportion

# .calcStatsAAV_agg()
# Calculate average annual variation in catch, measure of the variability
# in landings for each aggregate stock/species.
.calcStatsAAV_agg <- function(  C_it = C_it,
                                tdx = tMP:nT,
                                qProbs = c(0.025,0.5,0.975)
                           )
{

  # Calculate diff
  diffC_it        <- aperm(apply( X = C_it[,tdx], FUN = diff, MARGIN = 1, na.rm = T ))
  absDiffC_it     <- abs(diffC_it)
  sumAbsDiffC_i   <- apply( X = absDiffC_it, FUN = sum, MARGIN = 1, na.rm = T )
  sumCatch_i      <- apply( X = C_it[,tdx], FUN = sum, MARGIN = 1, na.rm = T )

  AAV_i           <- sumAbsDiffC_i / sumCatch_i
  AAV_i[!is.finite(AAV_i)] <- 0

  AAV_q   <- quantile(AAV_i,probs = qProbs )

} # END calcStatsAAV_agg


# .calcStatsAAV_sp()
# Calculate average annual variation
# in catch, measure of the variability
# in landings for each stock/species.
.calcStatsAAV_sp <- function(  C_ispt = C_ispt,
                            tdx = tMP:nT,
                            qProbs = c(0.025,0.5,0.975),
                            margin = c(2,3) )
{
  C_ispt <- C_ispt[,,,tdx,drop = FALSE]
  # Append margin
  marg <- c(1,margin)

  # Add over margins, in case we are looking
  # complex/species/stock aggregates
  sumC_ispt <- apply( X = C_ispt, FUN = sum, MARGIN = c(marg,4), na.rm = T )


  # Calculate diff
  diffC_ispt        <- aperm(apply( X = sumC_ispt, FUN = diff, MARGIN = marg, na.rm = T ),c(margin,4,1))
  absDiffC_ispt     <- abs(diffC_ispt)
  sumAbsDiffC_isp   <- apply( X = absDiffC_ispt, FUN = sum, MARGIN = marg, na.rm = T )
  sumCatch_isp      <- apply( X = C_ispt, FUN = sum, MARGIN = marg, na.rm = T )

  AAV_isp           <- sumAbsDiffC_isp / sumCatch_isp
  AAV_isp[!is.finite(AAV_isp)] <- 0

  AAV_qsp   <- apply(X = AAV_isp, FUN = quantile, probs = qProbs, MARGIN = marg[-1] )

  return(AAV_qsp)
} # END .calcStatsAAV_sp


# .getCplxStats()
# A pared down perf stats function for complex level
# quantities (i.e. total catch across all species/stocks etc.)
# Useful for comparing omniscient manager obj fun
# weightings
.getCplxStats <- function( obj )
{
  # Model dims
  tMP <- obj$om$tMP
  nT  <- obj$om$nT
  nS  <- obj$om$nS
  nP  <- obj$om$nP

  # Names of stuff
  scenarioName  <- obj$ctlList$ctl$scenarioName
  mpName        <- obj$ctlList$ctl$mpName

  # Pull simulation label
  if(is.null(obj$simLabel))
  {
    simLabel <- stringr::str_split(obj$path,"/")[[1]]
    simLabel <- simLabel[length(simLabel)]
  } else simLabel <- obj$simLabel

  allConvReps <- obj$goodReps

  projYrs <- tMP:nT

  # E quantities
  C_ispt <- obj$om$C_ispt[allConvReps,,,projYrs, drop=FALSE]
  B_ispt <- obj$om$SB_ispt[allConvReps,,,projYrs, drop=FALSE]

  # Pull Bmsy
  Bmsy_sp <- obj$rp[[1]]$FmsyRefPts$BeqFmsy_sp

  # Turn into complex catch (aggregate across species and stocks)
  C_it <- apply( X = C_ispt, FUN = sum, MARGIN = c(1,4))

  # Calculate AAV at the complex level
  diffC_ti      <- apply( X = C_it, FUN = diff, MARGIN = 1)
  diffC_it      <- abs(t(diffC_ti))
  sumAbsDiff_i  <- apply( X = diffC_it, FUN = sum, MARGIN = 1)
  sumCatch_i    <- apply( X = C_it, FUN = sum, MARGIN = 1)
  AAVC_i        <- sumAbsDiff_i / sumCatch_i
  AAVCplx       <- round(median( AAVC_i),2)
  
  # We want median (over reps) average and total complex catch
  Cbar    <- round(median(apply( X = C_it, FUN = mean, MARGIN = 1 ) ),2)
  totC    <- round(median(apply( X = C_it, FUN = sum,  MARGIN = 1 ) ),2)

  # Now loop over species and stocks, calculate
  # depletion probs (make the following values inputs either
  # as function args or in the ctl file)
  LRP <- 0.4
  USR <- 2

  depBmsy_ispt <- B_ispt
  for( s in 1:nS )
    for( p in 1:nP )
    {
      depBmsy_ispt[,s,p,] <- B_ispt[,s,p,]/Bmsy_sp[s,p]
    }

  depGtLRP_ispt <- depBmsy_ispt > LRP
  depLtUSR_ispt <- depBmsy_ispt < USR

  probDepGtLRP_isp <- apply(X = depGtLRP_ispt, FUN = mean, MARGIN = c(1,2,3))
  probDepLtUSR_isp <- apply(X = depLtUSR_ispt, FUN = mean, MARGIN = c(1,2,3))

  # Get distribution of depletion within desired region 
  # across species, stocks, and replicates
  distProbDepLRP_q <- quantile( probDepGtLRP_isp, probs = c(0.025, 0.5, 0.975) )
  distProbDepUSR_q <- quantile( probDepLtUSR_isp, probs = c(0.025, 0.5, 0.975) )

  # Make above vectors into a distribution as "med (CI limits)"
  distLRP_chr <- paste(distProbDepLRP_q[2], " (", distProbDepLRP_q[1],", ", distProbDepLRP_q[3], ")", sep = "")
  distUSR_chr <- paste(distProbDepUSR_q[2], " (", distProbDepUSR_q[1],", ", distProbDepUSR_q[3], ")", sep = "")

  # Make into a list
  out.df <- data.frame( simLabel    = simLabel,
                        scenario    = scenarioName,
                        mp          = mpName,
                        Cbar        = round(Cbar,2),
                        totC        = round(totC,2),
                        AAV         = round(AAVCplx,2),
                        PGtLRP      = distLRP_chr,
                        PLtUSR      = distUSR_chr )

  out.df
} # END .getCplxStats()


# .getOmniInfo()
# Pulls omniscient manager information from the blob and control
# list.
.getOmniInfo <- function(obj)
{
  # Model dims
  tMP <- obj$om$tMP
  nT  <- obj$om$nT
  nS  <- obj$om$nS
  nP  <- obj$om$nP

  # Names of stuff
  scenarioName  <- obj$ctlList$ctl$scenarioName
  mpName        <- obj$ctlList$ctl$mpName
  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames

  # Pull simulation label
  if(is.null(obj$simLabel))
  {
    simLabel <- stringr::str_split(obj$path,"/")[[1]]
    simLabel <- simLabel[length(simLabel)]
  } else simLabel <- obj$simLabel



  allConvReps <- obj$goodReps

  projYrs <- tMP:nT

  # Get omni sub-ctlList
  omniList    <- obj$ctlList$mp$omni
  omniObjFun  <- obj$omniObjFun

  # Pull obj fun targets
  loDepBmsy     <- omniList$loDepBmsy
  hiDepBmsy     <- omniList$hiDepBmsy
  maxF          <- omniList$maxF
  maxE          <- omniList$maxE
  maxAAV        <- omniList$maxAAV
  maxRelEffDiff <- omniList$maxRelEffDiff
  maxRelCatDiff <- omniList$maxRelCatDiff
  # component weights
  avgCatWt      <- omniList$avgCatWt
  totCatWt      <- omniList$totCatWt
  sumCatWt      <- omniList$sumCatWt
  hiDepBmsyWt   <- omniList$hiDepBmsyWt
  loDepBmsyWt   <- omniList$loDepBmsyWt
  probDepWt     <- omniList$probDepWt
  closedWt      <- omniList$closedWt
  AAVWt         <- omniList$AAVWt
  catDiffWt     <- omniList$catDiffWt
  effDiffWt     <- omniList$effDiffWt
  initCatDiffWt <- omniList$initCatDiffWt
  initEffDiffWt <- omniList$initEffDiffWt

  # Objective function penalty controls
  penType       <- omniList$penType
  linBeta       <- omniList$linBeta

  # Now pull objective function values

  # Complex level
  objFun_i            <- omniObjFun$objFun_i[allConvReps]
  totCbar_i           <- omniObjFun$totCbar_i[allConvReps]
  Csum_i              <- omniObjFun$Csum_i[allConvReps]
  barEffDiff_ip       <- omniObjFun$barEffDiff_ip[allConvReps,]


  # Species/stock
  objFun_isp          <- omniObjFun$objFun_isp[allConvReps,,]
  Cbar_isp            <- omniObjFun$Cbar_isp[allConvReps,,]
  barLoDep_isp        <- omniObjFun$barLoDep_isp[allConvReps,,]
  barHiDep_isp        <- omniObjFun$barHiDep_isp[allConvReps,,]
  barProbLoDep_isp    <- omniObjFun$barProbLoDep_isp[allConvReps,,]
  barProbHiDep_isp    <- omniObjFun$barProbHiDep_isp[allConvReps,,]
  barInitCatDiff_isp  <- omniObjFun$barInitCatDiff_isp[allConvReps,,]
  closedCount_isp     <- omniObjFun$closedCount_isp[allConvReps,,]
  barAAV_isp          <- omniObjFun$barAAV_isp[allConvReps,,]
  barCatDiff_isp      <- omniObjFun$barCatDiff_isp[allConvReps,,]

  # Sum effort diff over stocks
  barEffDiff_i        <- omniObjFun$barEffDiff_i[allConvReps]
  barInitEffDiff_i    <- omniObjFun$barInitEffDiff_i[allConvReps]

  # Now, we want to make a couple of tables
  # make a table of weights as well, this will be good for
  # tuning the weights later...
  wtTableColNames <- c( 'simLabel',
                        'scenario',
                        'mp',
                        'wt_CbarCplx',
                        'wt_CbarStock',
                        'wt_totC',
                        'wt_hiDep',
                        'wt_loDep',
                        'wt_probDep',
                        'wt_closure',
                        'wt_AAV',
                        'wt_catDiff',
                        'wt_effDiff',
                        'wt_initCatDiff',
                        'wt_initEffDiff')

  wtTable <- matrix(NA, nrow = 1, ncol = length(wtTableColNames) )
  colnames(wtTable) <- wtTableColNames

  wtTable[,'simLabel']        <- simLabel
  wtTable[,'scenario']        <- scenarioName
  wtTable[,'mp']              <- mpName
  wtTable[,'wt_CbarCplx']     <- totCatWt
  wtTable[,'wt_CbarStock']    <- avgCatWt
  wtTable[,'wt_totC']         <- sumCatWt
  wtTable[,'wt_hiDep']        <- hiDepBmsyWt
  wtTable[,'wt_loDep']        <- loDepBmsyWt
  wtTable[,'wt_probDep']      <- probDepWt
  wtTable[,'wt_closure']      <- closedWt
  wtTable[,'wt_AAV']          <- AAVWt
  wtTable[,'wt_catDiff']      <- catDiffWt
  wtTable[,'wt_effDiff']      <- effDiffWt
  wtTable[,'wt_initCatDiff']  <- initCatDiffWt
  wtTable[,'wt_initEffDiff']  <- initCatDiffWt


  # First, the complex level quantities
  cplxTableColNames <- c( 'simLabel',
                          'scenario',
                          'mp',
                          'objFun',
                          'wt_CbarCplx',
                          'of_CbarCplx',
                          'wt_totC',
                          'of_totC',
                          'wt_EffDiff',
                          'of_EffDiff',
                          'wt_initEffDiff',
                          'of_initEffDiff' )
  cplxOmniObjFunTable <- matrix(NA, nrow = 1, ncol = length(cplxTableColNames) )
  colnames(cplxOmniObjFunTable) <- cplxTableColNames
  cplxOmniObjFunTable <- as.data.frame(cplxOmniObjFunTable)

  cplxOmniObjFunTable[,'simLabel']    <- simLabel
  cplxOmniObjFunTable[,'scenario']    <- scenarioName
  cplxOmniObjFunTable[,'mp']          <- mpName
  cplxOmniObjFunTable[,'objFun']      <- round(mean(objFun_i),2)
  cplxOmniObjFunTable[,'of_CbarCplx'] <- round(mean(totCbar_i),2)
  cplxOmniObjFunTable[,'wt_CbarCplx'] <- totCatWt
  cplxOmniObjFunTable[,'of_totC']     <- round(mean(Csum_i),2)
  cplxOmniObjFunTable[,'wt_totC']     <- sumCatWt
  cplxOmniObjFunTable[,'of_EffDiff']  <- round(mean(barEffDiff_i),2)
  cplxOmniObjFunTable[,'wt_EffDiff']  <- effDiffWt
  cplxOmniObjFunTable[,'of_EffDiff']  <- round(mean(barInitEffDiff_i),2)
  cplxOmniObjFunTable[,'wt_EffDiff']  <- initEffDiffWt


  # Now a stock/species objective function table
  stockTableColNames <- c(  'simLabel',
                            'scenario',
                            'mp',
                            'species',
                            'stock',
                            'objFun',
                            'wt_CbarStock',
                            'of_CbarStock',
                            'wt_hiDep',
                            'of_hiDep',
                            'wt_loDep',
                            'of_loDep',
                            'wt_probDep',
                            'of_probHiDep',
                            'of_probLoDep',
                            'wt_closure',
                            'of_closure',
                            'wt_AAV',
                            'of_AAV',
                            'wt_catDiff',
                            'of_catDiff',
                            'wt_initCatDiff',
                            'of_initCatDiff' )

  stockObjFunTable <- matrix( NA, ncol = length(stockTableColNames), nrow = nS * nP )
  colnames(stockObjFunTable) <- stockTableColNames
  stockObjFunTable <- as.data.frame(stockObjFunTable)

  stockObjFunTable[,'simLabel']  <- simLabel
  stockObjFunTable[,'scenario']  <- scenarioName
  stockObjFunTable[,'mp']        <- mpName


  for( s in 1:nS )
    for( p in 1:nP )
    {
      rowIdx <- (s - 1) * nP + p
      stockObjFunTable[rowIdx,'species']        <- speciesNames[p]
      stockObjFunTable[rowIdx,'stock']          <- stockNames[p]
      stockObjFunTable[rowIdx,'objFun']         <- round(mean( objFun_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'wt_CbarStock']   <- avgCatWt
      stockObjFunTable[rowIdx,'of_CbarStock']   <- round(mean( Cbar_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'wt_hiDep']       <- hiDepBmsyWt
      stockObjFunTable[rowIdx,'of_hiDep']       <- round(mean( barHiDep_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'wt_loDep']       <- loDepBmsyWt
      stockObjFunTable[rowIdx,'of_loDep']       <- round(mean( barLoDep_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'wt_probDep']     <- probDepWt
      stockObjFunTable[rowIdx,'of_probHiDep']   <- round(mean( barProbHiDep_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'of_probLoDep']   <- round(mean( barProbLoDep_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'wt_closure']     <- closedWt
      stockObjFunTable[rowIdx,'of_closure']     <- round(mean( closedCount_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'wt_AAV']         <- AAVWt
      stockObjFunTable[rowIdx,'of_AAV']         <- round(mean( barAAV_isp[,s,p] ),2)
      stockObjFunTable[rowIdx,'wt_catDiff']     <- catDiffWt
      stockObjFunTable[rowIdx,'of_catDiff']     <- round(mean( barCatDiff_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'wt_initCatDiff'] <- initCatDiffWt
      stockObjFunTable[rowIdx,'of_initCatDiff'] <- round(mean( barInitCatDiff_isp[,s,p]),2)
    }


  outList <- list(  wtTable     = wtTable,
                    cplxTable   = cplxOmniObjFunTable,
                    stockTable  = stockObjFunTable )


  outList
} # END .getOmniInfo()



