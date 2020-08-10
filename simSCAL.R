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
# Last Update: Nov 13, 2019
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><>

runMS3 <- function( ctlFile = "./simCtlFile.txt",
                    outFolder = NULL )
{
  # Load control list
  ctlTable <- .readParFile ( ctlFile )
  # Create control list
  ctlList  <- .createList  ( ctlTable )

  # Load history
  repList               <- .loadFit( ctlList$opMod$histFile )
  ctlList$opMod$histRpt <- c( repList$data, repList$repOpt ) 
  ctlList$opMod$fYear   <- repList$fYear
  ctlList$opMod$species <- repList$species
  ctlList$opMod$stocks  <- repList$stocks
  ctlList$opMod$fleets  <- repList$gearLabs

  # Initialise the data structures
  # to hold simulation states etc.
  simObj <- .initMS3pop(  ctlList )

  # Run mgmtProc
  blob <- .mgmtProc( simObj )

  # Save output
  .saveBlob(blob = blob, ctlTable, outFolder)

  beepr::beep("complete")
}


# calcObsTimes()
# Creates a schedule for survey indices and
# assessments based on control file
# settings
# inputs:   obj: simObj created in mgmtProc
# outputs:  obj: simObj with schedules attached
.calcTimes <- function( obj )
{
  # Pull ctlList
  ctlList <- obj$ctlList

  tMP <- obj$om$tMP
  nT  <- obj$om$nT
  pT  <- ctlList$opMod$pT

  nS <- obj$om$nS
  nP <- obj$om$nP
  nF <- obj$om$nF


  # Get intervals for assessments
  assInt    <- obj$ctlList$mp$assess$assInt
  assTimes  <- seq(from = tMP, by = assInt, length.out = pT)

  assessOn <- c(rep(NA, tMP-1),rep(FALSE,pT))
  assessOn[ assTimes ]<- TRUE

  # Now do the same for indices
  idxOn   <- obj$ctlList$mp$data$idxOn
  obsInt  <- obj$ctlList$mp$data$obsInt_i

  # Start with a FALSE array, add pooled data dimension
  idxOn_spft <- array( FALSE, dim = c(nS+1,nP+1,nF,nT) )
  # Recover historical indices and set 
  # to TRUE
  histdx <- 1:(tMP-1)
  I_spft <- obj$mp$data$I_spft[,,,histdx]

  for( s in 1:nS )
    for( p in 1:nP )
      for( fIdx in 1:nF )
      {
        histOn <- which(!is.na(I_spft[s,p,fIdx,]))
        if( length(histOn) == 0 )
          next

        idxOn_spft[s,p,fIdx,histOn] <- TRUE

        # Now generate new data for the future
        if( fIdx %in% idxOn )
        {
          i <- which( idxOn == fIdx )

          lastObs <- max(histOn)
          if(obsInt[i] == 1)
            lastObs <- tMP

          newObs <- seq( from = lastObs, to = nT, by = obsInt[i] )
          idxOn_spft[s,p,fIdx,newObs] <- TRUE
        }
      }

  # Now fill the pooled data dimensions
  idxOn_spft[s+1,1:nP,,]  <- as.logical( apply( X = idxOn_spft[1:nS,1:nP,,], FUN = sum, MARGIN = c(2,3,4) ) )
  idxOn_spft[1:nS,p+1,,]  <- as.logical( apply( X = idxOn_spft[1:nS,1:nP,,], FUN = sum, MARGIN = c(1,3,4) ) )
  idxOn_spft[s+1,p+1,,]   <- as.logical( apply( X = idxOn_spft[1:nS,1:nP,,], FUN = sum, MARGIN = c(3,4)   ) )

  obj$mp$assess$assessOn_t  <- assessOn
  obj$mp$data$idxOn_spft    <- idxOn_spft

  return(obj)
} # END .calcTimes()


# .updatePop       
# Purpose:        Complete one iteration of feedback system
# Parameters:     obj=complete simulation object up to t-1; t=current time
# Returns:        obj=complete simulation object up to t
# Source:         SDN Johnson, modified from mseR by S.P. Cox
.updatePop <- function( obj, t )
{
  # Flow of updatePop 
  # 1. extract om and mp objects
  # 2. perform stock assessment
  # 3. Determine stock status and set catch limit
  # 4. update operating model state variables and data
  # 5. return updated simulation object

  # 1. Extract operating model and management procedure objects
  ctlList <- obj$ctlList

  rp <- obj$rp      # reference points
  om <- obj$om      # operating model
  mp <- obj$mp      # management procedure

  nT  <- om$nT
  tMP <- om$tMP
  nS  <- om$nS
  nP  <- om$nP

  # 2. Perform stock assessment
  obj  <- .callProcedureAM( obj, t )


  # 3. Determine stock status and set catch limit
  if( ctlList$mp$hcr$type == "ramped")
    obj  <- .applyPrecHCR( obj, t )

  if( ctlList$mp$hcr$type == "conF" )
    obj  <- .applyConstantF( obj, t )

  # Now spread catch among allocation
  for( s in 1:nS )
    for( p in 1:nP )
    {
      obj$mp$hcr$TAC_spft[s,p,,t]    <- obj$mp$hcr$TAC_spt[s,p,t] * obj$om$alloc_spf[s,p,]
    }

  # 4. Update simulated population
  obj <- .ageSexOpMod(obj, t)

  # 5. Return updated sim object
  return(obj)

} # END .updatePop()


# .applyConstantF()
# Applies a harvest control rule to
# the assessment outputs
.applyConstantF <- function( obj, t )
{
  hcr     <- obj$mp$hcr
  ctlList <- obj$ctlList
  assess  <- obj$mp$assess

  # Pull reference points
  refPtList   <- obj$rp
  refPtType   <- ctlList$mp$hcr$refPtType

  # calculate projection time
  tMP <- obj$om$tMP
  pt  <- t - tMP + 1

  # Get HCR quantities
  Fref_sp <- hcr$Fref_spt[,,t]

  # Model dims
  nS  <- obj$om$nS
  nP  <- obj$om$nP

  # Stock status
  projSB_sp   <- assess$retroSB_tspt[pt,,,t] 
  projVB_sp   <- assess$retroVB_tspft[pt,,,2,t] 

  # Pull data
  I_spft <- obj$mp$data$I_spft

  # Use the last two years of data in Synoptic
  # and calculate mean over time to get
  # splitting weights
  rctMeanI_sp <- apply( X = I_spft[1:nS,1:nP,4,(t-2):(t-1)],
                        FUN = mean,
                        MARGIN = c(1,2), na.rm = T )

  # Calculate total recent mean index for a species
  rctMeanI_s <- apply( X = rctMeanI_sp, FUN = sum, MARGIN = 1)

  # if( t >= tMP )
  #   browser()

  
  # Compute target F
  hcr$targetF_spt[,,t] <- Fref_sp
  
  # Now apply "F" as a harvest rate, since
  # we don't have M information
  hcr$TAC_spt[,,t] <- hcr$targetF_spt[,,t] * projVB_sp



  propTAC_sp <- array(1, dim = c(nS,nP) )

  if( ctlList$mp$data$speciesPooling & ctlList$mp$data$spatialPooling )
  {
    # Calc proportion of TAC for each species
    propTAC_s  <- rctMeanI_s / sum(rctMeanI_s)
    
    # Distribute among stocks
    propTAC_sp <- rctMeanI_sp
    for( s in 1:nS )
      propTAC_sp[s,] <- propTAC_sp[s,] / rctMeanI_s[s] * propTAC_s[s]

    hcr$TAC_spt[,,t] <- hcr$TAC_spt[,,t] * propTAC_sp
  }

  if( ctlList$mp$data$spatialPooling & !ctlList$mp$data$speciesPooling )
  {
    # Distribute among stocks
    propTAC_sp <- rctMeanI_sp
    for( s in 1:nS )
      propTAC_sp[s,] <- propTAC_sp[s,] / rctMeanI_s[s]

    hcr$TAC_spt[,,t] <- hcr$TAC_spt[,,t] * propTAC_sp
  }

  if( ctlList$mp$data$speciesPooling & !ctlList$mp$data$spatialPooling )
  {
    # Distribute among species
    propTAC_sp <- rctMeanI_sp
    for( p in 1:nP )
      propTAC_sp[,p] <- propTAC_sp[,p] / sum(propTAC_sp[,p])

    hcr$TAC_spt[,,t] <- hcr$TAC_spt[,,t] * propTAC_sp
  }

  # Apply smoother on TACs
  if( !is.null(ctlList$mp$hcr$maxDeltaTACup) )
  {
    # Determine which TACs are more than maxDeltaTAC
    # above last year
    maxDeltaTACup  <- ctlList$mp$hcr$maxDeltaTACup
    currTAC_sp     <- hcr$TAC_spt[,,t]
    lastTAC_sp     <- hcr$TAC_spt[,,t-1]
    diffTAC_sp     <- currTAC_sp - lastTAC_sp
    DeltaTACup_sp  <- diffTAC_sp / hcr$TAC_spt[,,t-1]
    deltaDiff_sp   <- DeltaTACup_sp - maxDeltaTACup
    smoothIdx <- which( deltaDiff_sp > 0 )

    # Apply smoother
    currTAC_sp[smoothIdx] <- (1 + maxDeltaTACup) * lastTAC_sp[smoothIdx]
    hcr$TAC_spt[,,t] <- currTAC_sp
  }


  # Save proportion of TAC for use in retro SB plots
  hcr$propTAC_spt[,,t] <- propTAC_sp

  # Now put all the lists we modified back into obj
  obj$mp$hcr    <- hcr
  obj$mp$assess <- assess

  return(obj)
} # END .applyConstantF()

# .applyPrecHCR()
# Applies a harvest control rule to
# the assessment outputs
.applyPrecHCR <- function( obj, t )
{
  hcr     <- obj$mp$hcr
  ctlList <- obj$ctlList
  assess  <- obj$mp$assess

  # calculate projection time
  tMP <- obj$om$tMP
  pt  <- t - tMP + 1

  # Get HCR quantities
  Fref_sp <- hcr$Fref_spt[,,t]
  Bref_sp <- hcr$Bref_spt[,,t]

  # Create control points
  hcr$LCP_spt[,,t] <- Bref_sp * ctlList$mp$hcr$LCP
  hcr$UCP_spt[,,t] <- Bref_sp * ctlList$mp$hcr$UCP

  # Calculate stock status
  projSB_sp   <- assess$retroSB_tspt[pt,,,t] 
  projVB_sp   <- assess$retroVB_tspft[pt,,,2,t] 

  nS  <- obj$om$nS
  nP  <- obj$om$nP

  # Pull data
  I_spft <- obj$mp$data$I_spft

  # Use the last two years of data in Synoptic
  # and calculate mean over time to get
  # splitting weights
  rctMeanI_sp <- apply( X = I_spft[1:nS,1:nP,4,(t-2):(t-1)],
                        FUN = mean,
                        MARGIN = c(1,2), na.rm = T )

  # Calculate total recent mean index for a species
  rctMeanI_s <- apply( X = rctMeanI_sp, FUN = sum, MARGIN = 1)

  # if( t >= tMP )
  #   browser()

  # Private function to calculate ramped HCR
  calcRampedHCR <- function( B, LCP, UCP, Fref, lowFmult )
  {
    if( B < LCP )
      F <- lowFmult * Fref
    if( LCP <= B & B < UCP )
      F <- (lowFmult + (1 - lowFmult) * (B - LCP) / (UCP - LCP)) * Fref
    if( UCP <= B )
      F <- Fref

    F
  }

  # Compute target F
  hcr$targetF_spt[,,t] <- mapply( FUN = calcRampedHCR,
                                  B = projSB_sp, 
                                  LCP = hcr$LCP_spt[,,t],
                                  UCP = hcr$UCP_spt[,,t],
                                  Fref = Fref_sp,
                                  MoreArgs = list( lowFmult = .1) )
  
  # Now apply "F" as a harvest rate, since
  # we don't have M information
  hcr$TAC_spt[,,t] <- hcr$targetF_spt[,,t] * projVB_sp

  propTAC_sp <- array(1, dim = c(nS,nP) )

  if( ctlList$mp$data$speciesPooling & ctlList$mp$data$spatialPooling )
  {
    # Calc proportion of TAC for each species
    propTAC_s  <- rctMeanI_s / sum(rctMeanI_s)
    
    # Distribute among stocks
    propTAC_sp <- rctMeanI_sp
    for( s in 1:nS )
      propTAC_sp[s,] <- propTAC_sp[s,] / rctMeanI_s[s] * propTAC_s[s]

    hcr$TAC_spt[,,t] <- hcr$TAC_spt[,,t] * propTAC_sp
  }

  if( ctlList$mp$data$spatialPooling & !ctlList$mp$data$speciesPooling )
  {
    # Distribute among stocks
    propTAC_sp <- rctMeanI_sp
    for( s in 1:nS )
      propTAC_sp[s,] <- propTAC_sp[s,] / rctMeanI_s[s]

    hcr$TAC_spt[,,t] <- hcr$TAC_spt[,,t] * propTAC_sp
  }

  if( ctlList$mp$data$speciesPooling & !ctlList$mp$data$spatialPooling )
  {
    # Distribute among species
    propTAC_sp <- rctMeanI_sp
    for( p in 1:nP )
      propTAC_sp[,p] <- propTAC_sp[,p] / sum(propTAC_sp[,p])

    hcr$TAC_spt[,,t] <- hcr$TAC_spt[,,t] * propTAC_sp
  }

  # Save proportion of TAC for use in retro SB plots
  hcr$propTAC_spt[,,t] <- propTAC_sp

  # Apply smoother on TACs
  if( !is.null(ctlList$mp$hcr$maxDeltaTACup) )
  {
    # Determine which TACs are more than maxDeltaTAC
    # above last year
    maxDeltaTACup  <- ctlList$mp$hcr$maxDeltaTACup
    currTAC_sp     <- hcr$TAC_spt[,,t]
    lastTAC_sp     <- hcr$TAC_spt[,,t-1]
    diffTAC_sp     <- currTAC_sp - lastTAC_sp
    DeltaTACup_sp  <- diffTAC_sp / hcr$TAC_spt[,,t-1]
    deltaDiff_sp   <- DeltaTACup_sp - maxDeltaTACup
    smoothIdx <- which( deltaDiff_sp > 0 )

    # Apply smoother
    currTAC_sp[smoothIdx] <- (1 + maxDeltaTACup) * lastTAC_sp[smoothIdx]
    hcr$TAC_spt[,,t] <- currTAC_sp
  }

  # Now put all the lists we modified back into obj
  obj$mp$hcr    <- hcr
  obj$mp$assess <- assess

  return(obj)
} # END .applyPrecHCR



# .callProcedureAM()
# Generic function to call estimation
# procedure. Uses the do.call function
# to apply the AM using the method name
# on the simObj
.callProcedureAM <- function( obj, t)
{
  mp <- obj$mp
  ctlList <- obj$ctlList

  # Get procedure name
  methodID <- ctlList$mp$assess$method

  # AMfunName <- paste(".AM_",methodID, sep = "")

  if( methodID == "idxBased")
    obj <- .AM_idxBased( obj, t)

  if( methodID == "PerfectInfo" )
    obj <- .AM_perfectInfo(obj,t)

  if( methodID == "hierProd" )
    obj <- .AM_hierProd(obj,t)


  return(obj)
} # END .callProcedureAM()

# solvePTm()
# Function for numerically solving 
# for the m parameter in a Pella Tomlinson
# model given an OM B0 and Bmsy
solvePTm <- function( Bmsy, B0 )
{
  # First, find optFrac
  optFrac <- Bmsy / B0 


  if( optFrac < 0.5 )
    bounds <- c(0,2)

  if( optFrac > .5 )
    bounds <- c(2,1e2)

  m <- seq( from = bounds[1], to = bounds[2], 
            length.out = 1e3)
  m <- m[m != 1]

  mFrac <- 1 / (m^(1/(m-1)))

  mSpline <- splinefun(x = m, y = mFrac - optFrac )

  mSol <- uniroot( mSpline, interval = bounds )$root

  mSol
} # END solvePTm()


# Apply hierarchical production model 
# AM - basically the AM from 
# Johnson and Cox, 2019
.AM_hierProd <- function( obj, t )
{
  # Pull control list
  ctlList     <- obj$ctlList
  refPtList   <- obj$rp

  refPtType   <- ctlList$mp$hcr$refPtType

  # Pull reference points
  refPtsSS <- refPtList$FmsyRefPts
  refPtsMS <- refPtList$EmsyMSRefPts

  # Single species yield and biomass
  YeqSS_sp    <- refPtsSS$YeqFmsy_sp
  BeqSS_sp    <- refPtsSS$BeqFmsy_sp
  expBeqSS_sp <- refPtsSS$expBeqFmsy_sp
  UmsySS_sp   <- YeqSS_sp/expBeqSS_sp

  # Multispecies yield and biomass
  YeqMS_sp    <- refPtsMS$YeqEmsy_sp
  BeqMS_sp    <- refPtsMS$BeqEmsy_sp
  expBeqMS_sp <- refPtsMS$expBeqEmsy_sp

  if( refPtType == "SS" )
  {
    Yeq_sp <- YeqSS_sp
    Beq_sp <- BeqSS_sp
  }

  if( refPtType == "MS" )
  {
    Yeq_sp <- YeqMS_sp
    Beq_sp <- BeqMS_sp
  }

  # Pull B0 for m value calcs
  B0_sp   <- ctlList$opMod$histRpt$B0_sp



  # Get mean lnUmsy from SS ref points, as this
  # reflects the true stock dynamics
  mlnUmsy_SS <- log(sum(YeqSS_sp) / sum(BeqSS_sp))

  # Get complex dims
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nF  <- obj$om$nF
  tMP <- obj$om$tMP

  # PT m parameter for skew in yield curves
  omPTm_sp <- array(2, dim = c(nS,nP))

  if( ctlList$mp$assess$spSkewYieldCurves )
    for( p in 1:nP )
      omPTm_sp[,p]  <- mapply(  FUN = solvePTm, 
                                Bmsy = BeqSS_sp[,p],
                                B0 = B0_sp[,p] )


  # Pull EBSBpars and think about how to use them for
  # data pooling and coastwide - I think we'll just 
  # average them for now... NO - DP and CW JSel pars
  # can be estimated in the ref pts calcs.
  EBSBpars <- refPtList$EBSBpars

  # Calculate projection year
  pt <- t - tMP + 1
  pT <- ctlList$opMod$pT

  # Pull OM devs
  omqDevs_spft <- array(0, dim = c(nS,nP,nF,t-1))
  omqDevs_spft[,,,1:(tMP-1)] <- ctlList$opMod$histRpt$epslnq_spft

  # Pull OM q values
  mq_spf      <- ctlList$opMod$histRpt$q_spf
  mq_sf       <- array(1, dim = c(nS,nF))
  mq_sf[,3:4] <- ctlList$opMod$histRpt$qSurv_sf
  sdlnq_f     <- ctlList$mp$assess$spsdlnq_f

  # Get data for fitting - biomass indices and catch
  spFleetIdx  <- ctlList$mp$assess$spFleets

  # Pull observations and catch
  I_spft      <- obj$mp$data$I_spft[1:(nS+1),1:(nP+1),,1:(t-1)]
  I_spft[,,-spFleetIdx,] <- -1
  C_spft      <- obj$om$C_spft[,,,(1:t-1)]
  C_spt       <- apply( X = C_spft, FUN = sum, 
                        MARGIN = c(1,2,4), na.rm = T )

  nu_spfk <- EBSBpars$stockSpec$nu_spfk
  P1_spf  <- EBSBpars$stockSpec$P1_spf
  P2_spf  <- EBSBpars$stockSpec$P2_spf

  spatialPooling <-  ctlList$mp$data$spatialPooling
  speciesPooling <-  ctlList$mp$data$speciesPooling

  # Spatial Pooling
  if( ctlList$mp$data$spatialPooling & !ctlList$mp$data$speciesPooling )
  {
    I_spft      <- I_spft[1:nS,nP+1,,1:(t-1),drop = FALSE]
    C_spt       <- array(NA, dim = c(nS,1,t-1))
    C_spt[,1,]  <- apply( X = C_spft, FUN = sum, MARGIN = c(1,4) )

    # Need to rearrange the EBSB model
    # fleet indices here.
    # get EBSB pars
    nu_spfk <- EBSBpars$coastWide$nu_spfk
    P1_spf  <- EBSBpars$coastWide$P1_spf
    P2_spf  <- EBSBpars$coastWide$P2_spf

  }
  # Species Pooling
  if( ctlList$mp$data$speciesPooling & !ctlList$mp$data$spatialPooling )
  {
    I_spft      <- I_spft[nS+1,1:nP,,1:(t-1),drop = FALSE]
    C_spt       <- array(NA, dim = c(1,nP,t-1))
    C_spt[1,,]  <- apply( X = C_spft, FUN = sum, MARGIN = c(2,4) )

    nu_spfk <- EBSBpars$dataPooled$nu_spfk
    P1_spf  <- EBSBpars$dataPooled$P1_spf
    P2_spf  <- EBSBpars$dataPooled$P2_spf

  }
  # Total Aggregation
  if( ctlList$mp$data$speciesPooling & ctlList$mp$data$spatialPooling )
  {
    I_spft      <- I_spft[nS+1,nP+1,,1:(t-1),drop = FALSE]
    C_spt       <- array(NA, dim = c(1,1,t-1))
    C_spt[1,1,] <- apply( X = C_spft, FUN = sum, MARGIN = c(4) )

    nu_spfk <- EBSBpars$totalAgg$nu_spfk
    P1_spf  <- EBSBpars$totalAgg$P1_spf
    P2_spf  <- EBSBpars$totalAgg$P2_spf
  }

  nSS <- min(dim(I_spft)[1],nS)
  nPP <- min(dim(I_spft)[2],nP)

  oldStock <- 1:nPP
  oldFleet <- 1:nF


  if( ctlList$mp$assess$spYrFirstIdx > 1 )
    I_spft[,,,1:ctlList$mp$assess$spYrFirstIdx ] <- -1

  # Switches for different model
  # complex structures
  spSingleStock   <- ctlList$mp$assess$spSingleStock
  spFixUmsy       <- ctlList$mp$assess$spFixUmsy
  spOMqDevs       <- ctlList$mp$assess$spOMqDevs


  I_spft[is.na(I_spft)] <- -1


  # Single stock version (no hierarchical method)
  if( spSingleStock | (nSS == 1 & nPP == 1) )
  {

    for( s in 1:nSS)
      for( p in 1:nPP )
      {
        mBmsy_sp  <- refPtList$FmsyRefPts$expBeqFmsy_sp

        omFref_sp <- Yeq_sp / Beq_sp
        omUmsy_sp <- YeqSS_sp / BeqSS_sp

        PTm_sp <- omPTm_sp
        
        if(spSingleStock & !speciesPooling & !spatialPooling)
        {
          mBmsy_sp  <- BeqSS_sp[s,p,drop = FALSE]

          # Calculate true Umsy values for MPs where
          # Umsy is fixed
          omFref_sp <- Yeq_sp[s,p,drop = FALSE] / Beq_sp[s,p,drop = FALSE]
          omUmsy_sp <- UmsySS_sp[s,p,drop = FALSE]

          PTm_sp <- omPTm_sp[s,p,drop = FALSE]
        }

        if( speciesPooling & spatialPooling )
        {
          mBmsy_sp        <- BeqSS_sp[1,1,drop = FALSE]
          mBmsy_sp[1,1]   <- sum(BeqSS_sp)

          omFref_sp       <- array(sum(Yeq_sp) / sum(Beq_sp), dim = c(1,1))
          omUmsy_sp       <- array(sum(YeqSS_sp) / sum(BeqSS_sp), dim = c(1,1))

          PTm_sp          <- omPTm_sp[1,1,drop = FALSE]
          # Solve for total agg yield curve skew
          if( ctlList$mp$assess$spSkewYieldCurves )
            PTm_sp[1,1]     <- solvePTm( Bmsy = sum(BeqSS_sp), B0 = sum(B0_sp) )
        }

        # Add Bmsy across stocks if coastwide
        if( spatialPooling & !speciesPooling & spSingleStock )
        {
          mBmsy_sp          <- BeqSS_sp[1,1,drop = FALSE]
          mBmsy_sp[1,1]     <- sum(BeqSS_sp[s,])

          PTm_sp <- omPTm_sp[,1,drop = FALSE] 

          omFref_sp       <- array(sum(Yeq_sp[s,]) / sum(Beq_sp[s,]), dim = c(1,1))
          omUmsy_sp       <- array(sum(YeqSS_sp[s,]) / sum(BeqSS_sp[s,]), dim = c(1,1))          

          PTm_sp          <- omPTm_sp[1,1,drop = FALSE]
          if( ctlList$mp$assess$spSkewYieldCurves )
            PTm_sp[1,1]    <- solvePTm( Bmsy = sum(BeqSS_sp[s,]), B0 = sum(B0_sp[s,]) )
          
        }

        # Add Bmsy across species if species pooling
        if( speciesPooling & !spatialPooling & spSingleStock )
        {
          mBmsy_sp          <- BeqSS_sp[1,1,drop = FALSE]
          mBmsy_sp[1,1]     <- sum(BeqSS_sp[,p])

          
          PTm_sp <- omPTm_sp[1,1,drop = FALSE] 

          omFref_sp       <- array(sum(Yeq_sp[,p]) / sum(Beq_sp[,p]), dim = c(1,1))
          omUmsy_sp       <- array(sum(YeqSS_sp[,p]) / sum(BeqSS_sp[,p]), dim = c(1,1))          

          if( ctlList$mp$assess$spSkewYieldCurves )
            PTm_sp[1,1]    <- solvePTm( Bmsy = sum(BeqSS_sp[,p]), B0 = sum(B0_sp[,p]) )
        }

        # Make lists for TMB AD function
        tmbLists <- .makeDatParHierProd(  C_spt[s,p,,drop = FALSE], 
                                          I_spft[s,p,,,drop = FALSE],
                                          ctlList,
                                          mlnUmsy = mlnUmsy_SS,
                                          lnUmsy_sp = log(omUmsy_sp),
                                          mBmsy_sp = mBmsy_sp,
                                          PTm_sp = PTm_sp,
                                          nu_spfk = nu_spfk,
                                          P1_spf = P1_spf,
                                          P2_spf = P2_spf,
                                          mq_spf = mq_spf,
                                          mq_sf  = mq_sf,
                                          sdlnq_f = sdlnq_f,
                                          omqDevs_spft = omqDevs_spft,
                                          oldStock = NULL,
                                          oldFleet = NULL,
                                          fixUmsy = spFixUmsy,
                                          importQdevs = spOMqDevs )

        # Set some phases - 
        # if no phase is set then pars not
        # estimated
        phases <- list( lnqShrink_sf    = ctlList$mp$assess$spPhzlnqShrink,
                        lnqFree_spf     = ctlList$mp$assess$spPhzlnqFree,
                        tvlnqDevs_vec   = ctlList$mp$assess$spPhzTVq,
                        lntauspf_vec    = ctlList$mp$assess$spPhztauObs,
                        lnBmsy_sp       = ctlList$mp$assess$spPhzBmsy,
                        lnUmsy          = ctlList$mp$assess$spPhzUmsy,
                        lnsigmaProc     = ctlList$mp$assess$spPhzSigmaProc,
                        zetaspt_vec     = ctlList$mp$assess$spPhzProcErr )

        if( is.null(ctlList$mp$assess$spTVqFleets) | spOMqDevs )
          phases$tvlnqDevs_vec <- -1

        if( ctlList$mp$assess$spSolveInitBio )
          phases$spPhzlnBinit <- -1

        if( spFixUmsy )
          phases$lnUmsy <- -1

        if( !ctlList$mp$assess$spLateStart )
          phases$lnBinit_vec <- -1


        phases$epslnUmsy_s  <- -1
        phases$epslnUmsy_sp <- -1
        phases$deltalnq_spf <- -1

        tmbLists$phases <- phases
        tmbLists$random <- ctlList$mp$assess$spRE


        # Apply AM, save retro biomass
        # and reference points
        amObj <- .applyTMBphase(  tmbLists = tmbLists,
                                  dllName = "hierProd",
                                  optimizer = "nlminb",
                                  silent = !ctlList$ctl$optimOutput,
                                  calcSD = TRUE,
                                  maxPhase = NULL,
                                  base_map = tmbLists$map,
                                  maxEval = ctlList$ctl$maxEval,
                                  maxIter = ctlList$ctl$maxIter,
                                  phaseMsg = ctlList$ctl$phaseMessages )
        # Need to save catchability
        # info for allocating catch
        # under coastwide models

        # Save biomass
        if(spSingleStock & !speciesPooling & !spatialPooling)
        {
          obj$mp$assess$retroSB_tspt[pt,s,p,1:t]    <- amObj$repOpt$B_spt[1,1,]
          for( f in 1:nF )
          {
            obj$mp$assess$retroVB_tspft[pt,s,p,f,1:t]     <- amObj$repOpt$vB_spft[1,1,f,]
            obj$mp$assess$retroq_tspf[pt,s,p,f]           <- amObj$repOpt$qhat_spf[1,1,f]
            obj$mp$assess$retroq_tspft[pt,s,p,f,1:(t-1)]  <- amObj$repOpt$qhat_spft[1,1,f,1:(t-1)]
            obj$mp$assess$retrotauObs_tspf[pt,s,p,f]      <- amObj$repOpt$tau_spf[1,1,f]
          }
          obj$mp$assess$retroUmsy_tsp[pt,s,p] <- amObj$repOpt$Umsy_sp
          obj$mp$assess$retroBmsy_tsp[pt,s,p] <- amObj$repOpt$Bmsy_sp

          if( ctlList$mp$hcr$Fsource == "est")
          {
            if( ctlList$mp$hcr$Fref == "Umsy" )
              obj$mp$hcr$Fref_spt[s,p,t]      <- amObj$repOpt$Umsy_sp

            if( ctlList$mp$hcr$Fref == "Fmsy" )
              obj$mp$hcr$Fref_spt[s,p,t]      <- amObj$repOpt$Umsy_sp
          }

          if( ctlList$mp$hcr$Fsource == "om" )
          {
            obj$mp$hcr$Fref_spt[s,p,t]      <- omFref_sp
          } 

          if( ctlList$mp$hcr$Bref == "Bmsy" )
            obj$mp$hcr$Bref_spt[s,p,t]      <- amObj$repOpt$Bmsy_sp  

          if( ctlList$mp$hcr$Bref == "B0" )
            obj$mp$hcr$Bref_spt[s,p,t]      <- amObj$repOpt$B0_sp


          # Save convergence info
          obj$mp$assess$pdHess_tsp[pt,s,p]  <- amObj$pdHess
          obj$mp$assess$posSDs_tsp[pt,s,p]  <- amObj$posSDs
          obj$mp$assess$maxGrad_tsp[pt,s,p] <- amObj$maxGrad
        }

        if( spSingleStock & speciesPooling & !spatialPooling )
        {
          for( ss in 1:nS)
          {
            obj$mp$assess$retroSB_tspt[pt,ss,p,1:t]    <- amObj$repOpt$B_spt[1,1,]
          
            for( f in 1:nF )
            {
              obj$mp$assess$retroVB_tspft[pt,ss,p,f,1:t]     <- amObj$repOpt$vB_spft[1,1,f,]
              obj$mp$assess$retroq_tspf[pt,ss,p,f]           <- amObj$repOpt$qhat_spf[1,1,f]
              obj$mp$assess$retroq_tspft[pt,ss,p,f,1:(t-1)]  <- amObj$repOpt$qhat_spft[1,1,f,1:(t-1)]
              obj$mp$assess$retrotauObs_tspf[pt,ss,p,f]      <- amObj$repOpt$tau_spf[1,1,f]
            }
          }
          obj$mp$assess$retroUmsy_tsp[pt,,p] <- amObj$repOpt$Umsy_sp
          obj$mp$assess$retroBmsy_tsp[pt,,p] <- amObj$repOpt$Bmsy_sp
        

          if( ctlList$mp$hcr$Fsource == "est")
          {
            if( ctlList$mp$hcr$Fref == "Umsy" )
              obj$mp$hcr$Fref_spt[,p,t]      <- amObj$repOpt$Umsy_sp

            if( ctlList$mp$hcr$Fref == "Fmsy" )
              obj$mp$hcr$Fref_spt[,p,t]      <- amObj$repOpt$Umsy_sp
          }

          if( ctlList$mp$hcr$Fsource == "om" )
          {
            obj$mp$hcr$Fref_spt[,p,t]      <- omFref_sp
          } 

          if( ctlList$mp$hcr$Bref == "Bmsy" )
            obj$mp$hcr$Bref_spt[,p,t]      <- amObj$repOpt$Bmsy_sp  

          if( ctlList$mp$hcr$Bref == "B0" )
            obj$mp$hcr$Bref_spt[,p,t]      <- amObj$repOpt$B0_sp


          # Save convergence info
          obj$mp$assess$pdHess_tsp[pt,,p]  <- amObj$pdHess
          obj$mp$assess$posSDs_tsp[pt,,p]  <- amObj$posSDs
          obj$mp$assess$maxGrad_tsp[pt,,p] <- amObj$maxGrad
          
        }

        if( spatialPooling & !speciesPooling & spSingleStock )
        {
          for( pp in 1:nP)
          {
            obj$mp$assess$retroSB_tspt[pt,s,pp,1:t]    <- amObj$repOpt$B_spt[1,1,]
          
            for( f in 1:nF )
            {
              obj$mp$assess$retroVB_tspft[pt,s,pp,f,1:t]     <- amObj$repOpt$vB_spft[1,1,f,]
              obj$mp$assess$retroq_tspf[pt,s,pp,f]           <- amObj$repOpt$qhat_spf[1,1,f]
              obj$mp$assess$retroq_tspft[pt,s,pp,f,1:(t-1)]  <- amObj$repOpt$qhat_spft[1,1,f,1:(t-1)]
              obj$mp$assess$retrotauObs_tspf[pt,s,pp,f]      <- amObj$repOpt$tau_spf[1,1,f]
            }
            

          }
          obj$mp$assess$retroUmsy_tsp[pt,s,] <- amObj$repOpt$Umsy_sp
          obj$mp$assess$retroBmsy_tsp[pt,s,] <- amObj$repOpt$Bmsy_sp

          if( ctlList$mp$hcr$Fsource == "est")
          {
            if( ctlList$mp$hcr$Fref == "Umsy" )
              obj$mp$hcr$Fref_spt[s,,t]      <- amObj$repOpt$Umsy_sp

            if( ctlList$mp$hcr$Fref == "Fmsy" )
              obj$mp$hcr$Fref_spt[s,,t]      <- amObj$repOpt$Umsy_sp
          }

          if( ctlList$mp$hcr$Fsource == "om" )
          {
            obj$mp$hcr$Fref_spt[s,,t]      <- omFref_sp
          } 

          if( ctlList$mp$hcr$Bref == "Bmsy" )
            obj$mp$hcr$Bref_spt[s,,t]      <- amObj$repOpt$Bmsy_sp  

          if( ctlList$mp$hcr$Bref == "B0" )
            obj$mp$hcr$Bref_spt[s,,t]      <- amObj$repOpt$B0_sp


          # Save convergence info
          obj$mp$assess$pdHess_tsp[pt,s,]  <- amObj$pdHess
          obj$mp$assess$posSDs_tsp[pt,s,]  <- amObj$posSDs
          obj$mp$assess$maxGrad_tsp[pt,s,] <- amObj$maxGrad
          
        }

        if( speciesPooling & spatialPooling )
        {
          for( ss in 1:nS )
          {
            for( pp in 1:nP)
            {
              for( f in 1:nF )
              {
                obj$mp$assess$retroVB_tspft[pt,ss,pp,f,1:t]     <- amObj$repOpt$vB_spft[1,1,f,]
                obj$mp$assess$retroq_tspf[pt,ss,pp,f]           <- amObj$repOpt$qhat_spf[1,1,f] 
                obj$mp$assess$retroq_tspft[pt,ss,pp,f,1:(t-1)]  <- amObj$repOpt$qhat_spft[1,1,f,1:(t-1)]
                obj$mp$assess$retrotauObs_tspf[pt,ss,pp,f]      <- amObj$repOpt$tau_spf[1,1,f]
              }
              # retro SB has to be in a loop because it's a ts, no guarantee
              # that the entries will enter in the right way (OK, there's a guarantee
              # but I never remember the orientation)
            
              obj$mp$assess$retroSB_tspt[pt,ss,pp,1:t]    <- amObj$repOpt$B_spt[1,1,1:t]
            }
              
            
          }
          obj$mp$assess$retroUmsy_tsp[pt,,] <- amObj$repOpt$Umsy_sp
          obj$mp$assess$retroBmsy_tsp[pt,,] <- amObj$repOpt$Bmsy_sp

          # Only one U/B value, so no need for loops here
          if( ctlList$mp$hcr$Fsource == "est")
          {
            if( ctlList$mp$hcr$Fref == "Umsy" )
              obj$mp$hcr$Fref_spt[,,t]      <- amObj$repOpt$Umsy_sp

            if( ctlList$mp$hcr$Fref == "Fmsy" )
              obj$mp$hcr$Fref_spt[,,t]      <- amObj$repOpt$Umsy_sp
          }

          if( ctlList$mp$hcr$Fsource == "om" )
          {
            obj$mp$hcr$Fref_spt[,,t]      <- omFref_sp
          } 


          if( ctlList$mp$hcr$Bref == "Bmsy" )
            obj$mp$hcr$Bref_spt[,,t]      <- amObj$repOpt$Bmsy_sp  

          if( ctlList$mp$hcr$Bref == "B0" )
            obj$mp$hcr$Bref_spt[,,t]      <- amObj$repOpt$B0_sp

          # Save convergence info
          obj$mp$assess$pdHess_tsp[pt,,]  <- amObj$pdHess
          obj$mp$assess$posSDs_tsp[pt,,]  <- amObj$posSDs
          obj$mp$assess$maxGrad_tsp[pt,,] <- amObj$maxGrad 
        }

      }

  }

  if( !spSingleStock & (nSS > 1 | nPP > 1) )
  {
    mBmsy_sp  <- refPtList$FmsyRefPts$expBeqFmsy_sp

    omFref_sp <- Yeq_sp / Beq_sp
    omUmsy_sp <- YeqSS_sp / BeqSS_sp

    PTm_sp <- omPTm_sp

    # Add Bmsy across stocks if coastwide
    if( spatialPooling &  !speciesPooling )
    {
      mBmsy_sp_new      <- mBmsy_sp[,1,drop = FALSE]
      mBmsy_sp_new[,1]  <- apply( X = BeqSS_sp, FUN = sum, MARGIN = 1)

      mBmsy_sp <- mBmsy_sp_new

      PTm_sp <- omPTm_sp[,1,drop = FALSE] 

      for( s in 1:nS )
      {
        omFref_sp[s,] <- sum(Yeq_sp[s,]) / sum(Beq_sp[s,])
        omUmsy_sp[s,] <- sum(YeqSS_sp[s,]) / sum(BeqSS_sp[s,])
        if( ctlList$mp$assess$spSkewYieldCurves )
          PTm_sp[s,]    <- solvePTm( Bmsy = sum(BeqSS_sp[s,]), B0 = sum(B0_sp[s,]) )
      }
    }

    if( speciesPooling & !spatialPooling )
    {
      mBmsy_sp_new <- mBmsy_sp[1,,drop = FALSE]
      mBmsy_sp_new[1,] <- apply( X = BeqSS_sp, FUN = sum, MARGIN = 2 )

      mBmsy_sp <- mBmsy_sp_new
      PTm_sp <- omPTm_sp[1,,drop = FALSE] 

      for( p in 1:nP )
      {
        omFref_sp[,p] <- sum(Yeq_sp[,p]) / sum(Beq_sp[,p])
        omUmsy_sp[,p] <- sum(YeqSS_sp[,p]) / sum(BeqSS_sp[,p])
        if( ctlList$mp$assess$spSkewYieldCurves )
          PTm_sp[,p]    <- solvePTm( Bmsy = sum(BeqSS_sp[,p]), B0 = sum(B0_sp[,p]) )
      }
    }

    tmbLists <- .makeDatParHierProd(  C_spt, 
                                      I_spft[1:nSS,1:nPP,,,drop = FALSE],
                                      ctlList,
                                      mlnUmsy = mlnUmsy_SS,
                                      lnUmsy_sp = log(omUmsy_sp),
                                      mBmsy_sp = mBmsy_sp,
                                      PTm_sp = PTm_sp,
                                      nu_spfk = nu_spfk,
                                      P1_spf = P1_spf,
                                      P2_spf = P2_spf,
                                      mq_spf = mq_spf,
                                      mq_sf  = mq_sf,
                                      sdlnq_f = sdlnq_f,
                                      omqDevs_spft = omqDevs_spft,
                                      oldStock = NULL,
                                      oldFleet = NULL,
                                      fixUmsy = spFixUmsy,
                                      importQdevs = spOMqDevs)

    # Set phases
    # if no phase is set then pars not
    # estimated

    # Need to put in phases that are responsive to changes in
    # the fleet structure...

    phases <- list( lnqShrink_sf    = ctlList$mp$assess$spPhzlnqShrink,
                    lnqFree_spf     = ctlList$mp$assess$spPhzlnqFree,
                    tvlnqDevs_vec   = ctlList$mp$assess$spPhzTVq,
                    lntauspf_vec    = ctlList$mp$assess$spPhztauObs,
                    lnBmsy_sp       = ctlList$mp$assess$spPhzBmsy,
                    lnUmsy          = ctlList$mp$assess$spPhzUmsy,
                    lnsigmaProc     = ctlList$mp$assess$spPhzSigmaProc,
                    zetaspt_vec     = ctlList$mp$assess$spPhzProcErr,
                    lnBinit_vec     = ctlList$mp$assess$spPhzlnBinit )

    if( nSS > 1 )
    {
      phases$epslnUmsy_s      <- ctlList$mp$assess$spPhzepslnUmsy_s
    }
    
    if( nPP > 1 )
    {
      phases$deltalnq_spf  <- ctlList$mp$assess$spPhzdeltalnq_spf
      phases$epslnUmsy_sp  <- ctlList$mp$assess$spPhzepslnUmsy_sp
    }

    # If totally fixing Umsy at OM reference
    if( spFixUmsy )
    {
      phases$lnUmsy       <- -1
      phases$epslnUmsy_s  <- -1
      phases$epslnUmsy_sp <- -1
    }

    if( !ctlList$mp$assess$spLateStart )
      phases$lnBinit_vec <- -1

    if( is.null(ctlList$mp$assess$spTVqFleets) | spOMqDevs )
      phases$tvlnqDevs_vec <- -1

    if( is.null(ctlList$mp$assess$spShrinkqFleets) | nPP == 1 )
      phases$deltalnq_spf <- -1


    if( ctlList$mp$assess$spSolveInitBio )
      phases$spPhzlnBinit <- -1

    # Add checks for 
    #   - condMLEs of comm fleet qs
    #   - fleets for shrinkage of q




    tmbLists$phases <- phases
    tmbLists$random <- ctlList$mp$assess$spRE


    # Apply AM, save retro biomass
    # and reference points
    amObj <- .applyTMBphase(  tmbLists = tmbLists,
                              dllName = "hierProd",
                              optimizer = "nlminb",
                              silent = !ctlList$ctl$optimOutput,
                              calcSD = TRUE,
                              maxPhase = NULL,
                              base_map = tmbLists$map,
                              phaseMsg = ctlList$ctl$phaseMessages,
                              maxEval = ctlList$ctl$maxEval,
                              maxIter = ctlList$ctl$maxIter )

    # browser()

    if( !spatialPooling & !speciesPooling )
    {
      # Save biomass
      obj$mp$assess$retroSB_tspt[pt,1:nSS,1:nPP,1:t]    <- amObj$repOpt$B_spt
      for( f in 1:nF )
      {
        obj$mp$assess$retroVB_tspft[pt,1:nSS,1:nPP,f,1:t]     <- amObj$repOpt$vB_spft[,,f,]
        obj$mp$assess$retroq_tspf[pt,1:nSS,1:nPP,f]           <- amObj$repOpt$qhat_spf[,,f]
        obj$mp$assess$retroq_tspft[pt,1:nSS,1:nPP,f,1:(t-1)]  <- amObj$repOpt$qhat_spft[1:nSS,1:nPP,f,1:(t-1)]
        obj$mp$assess$retrotauObs_tspf[pt,1:nSS,1:nPP,f]      <- amObj$repOpt$tau_spf[,,f]
      }

      obj$mp$assess$retroUmsy_tsp[pt,1:nSS,1:nPP] <- amObj$repOpt$Umsy_sp
      obj$mp$assess$retroBmsy_tsp[pt,1:nSS,1:nPP] <- amObj$repOpt$Bmsy_sp

      # Now reference points
      if( ctlList$mp$hcr$Fsource == "est")
      {
        if( ctlList$mp$hcr$Fref == "Umsy" )
          obj$mp$hcr$Fref_spt[1:nSS,1:nPP,t]      <- amObj$repOpt$Umsy_sp

        if( ctlList$mp$hcr$Fref == "Fmsy" )
          obj$mp$hcr$Fref_spt[1:nSS,1:nPP,t]      <- amObj$repOpt$Umsy_sp
      }

      if( ctlList$mp$hcr$Fsource == "om" )
      {
        obj$mp$hcr$Fref_spt[1:nSS,1:nPP,t]        <- omFref_sp
      }

      if( ctlList$mp$hcr$Bref == "Bmsy" )
        obj$mp$hcr$Bref_spt[1:nSS,1:nPP,t]      <- amObj$repOpt$Bmsy_sp  

      if( ctlList$mp$hcr$Bref == "B0" )
        obj$mp$hcr$Bref_spt[1:nSS,1:nPP,t]      <- 2*amObj$repOpt$Bmsy_sp

      # Save convergence info
      obj$mp$assess$pdHess_tsp[pt,1:nSS,1:nPP]  <- amObj$sdrep$pdHess
      obj$mp$assess$posSDs_tsp[pt,1:nSS,1:nPP]  <- amObj$posSDs
      obj$mp$assess$maxGrad_tsp[pt,1:nSS,1:nPP] <- amObj$maxGrad
    }

    if( spatialPooling )
    {
      # Save retro biomass
      for( p in 1:nP )
      {
        obj$mp$assess$retroSB_tspt[pt,1:nSS,p,1:t]    <- amObj$repOpt$B_spt[1:nSS,1,1:t]
        # Save retro exp. biomass and catchability
        for( f in 1:nF )
        {
          obj$mp$assess$retroVB_tspft[pt,1:nSS,p,f,1:t]     <- amObj$repOpt$vB_spft[1:nSS,1,f,1:t]
          obj$mp$assess$retroq_tspf[pt,1:nSS,p,f]           <- amObj$repOpt$qhat_spf[1:nSS,1,f]
          obj$mp$assess$retroq_tspft[pt,1:nSS,p,f,1:(t-1)]  <- amObj$repOpt$qhat_spft[1:nSS,1,f,1:(t-1)]
          obj$mp$assess$retrotauObs_tspf[pt,1:nSS,p,f]      <- amObj$repOpt$tau_spf[1:nSS,1,f]
        }
      }

      # Now reference points
      for( p in 1:nP )
      {
        obj$mp$assess$retroUmsy_tsp[pt,1:nSS,p] <- amObj$repOpt$Umsy_sp[1:nSS,1]
        obj$mp$assess$retroBmsy_tsp[pt,1:nSS,p] <- amObj$repOpt$Bmsy_sp[1:nSS,1]


        if( ctlList$mp$hcr$Fsource == "est")
        { 
          if( ctlList$mp$hcr$Fref == "Umsy" )
            obj$mp$hcr$Fref_spt[1:nSS,p,t]      <- amObj$repOpt$Umsy_sp[1:nSS,1]

          if( ctlList$mp$hcr$Fref == "Fmsy" )
            obj$mp$hcr$Fref_spt[1:nSS,p,t]      <- amObj$repOpt$Umsy_sp[1:nSS,1]
        }

        if( ctlList$mp$hcr$Fsource == "om")
        {
          obj$mp$hcr$Fref_spt[1:nSS,p,t]      <- omFref_sp[1:nSS,p] 
        }


        if( ctlList$mp$hcr$Bref == "Bmsy" )
          obj$mp$hcr$Bref_spt[1:nSS,p,t]      <- amObj$repOpt$Bmsy_sp[1:nSS,1]

        if( ctlList$mp$hcr$Bref == "B0" )
          obj$mp$hcr$Bref_spt[1:nSS,p,t]      <- 2*amObj$repOpt$Bmsy_sp[1:nSS,1]
      }

      # Save convergence info
      obj$mp$assess$pdHess_tsp[pt,1:nSS,]  <- amObj$sdrep$pdHess
      obj$mp$assess$posSDs_tsp[pt,1:nSS,]  <- amObj$posSDs
      obj$mp$assess$maxGrad_tsp[pt,1:nSS,] <- amObj$maxGrad
    }

    if( speciesPooling )
    {
      # Save retro biomass
      for( s in 1:nS )
        obj$mp$assess$retroSB_tspt[pt,s,1:nP,1:t]    <- amObj$repOpt$B_spt[1,1:nPP,1:t]

      # Save retro exp. biomass and catchability
      for( f in 1:nF )
      {
        for( s in 1:nS )
        {
          obj$mp$assess$retroVB_tspft[pt,s,1:nP,f,1:t]    <- amObj$repOpt$B_spt[1,1:nP,1:t]
          obj$mp$assess$retroq_tspf[pt,s,1:nP,f]          <- amObj$repOpt$qhat_spf[1,1:nP,f]
          obj$mp$assess$retroq_tspft[pt,s,1:nP,f,1:(t-1)] <- amObj$repOpt$qhat_spft[1,1:nP,f,1:(t-1)]
          obj$mp$assess$retrotauObs_tspf[pt,s,1:nP,f]     <- amObj$repOpt$tau_spf[1,1:nP,f]
        }
        
      }

      # Now reference points
      for( s in 1:nS )
      {
        obj$mp$assess$retroUmsy_tsp[pt,s,1:nP] <- amObj$repOpt$Umsy_sp[1,1:nP]
        obj$mp$assess$retroBmsy_tsp[pt,s,1:nP] <- amObj$repOpt$Bmsy_sp[1,1:nP]


        if( ctlList$mp$hcr$Fsource == "est")
        {
          if( ctlList$mp$hcr$Fref == "Umsy" )
            obj$mp$hcr$Fref_spt[s,1:nP,t]      <- amObj$repOpt$Umsy_sp[1,1:nP]

          if( ctlList$mp$hcr$Fref == "Fmsy" )
            obj$mp$hcr$Fref_spt[s,1:nP,t]      <- amObj$repOpt$Umsy_sp[1,1:nP]
        }

        if( ctlList$mp$hcr$Fsource == "om")
          obj$mp$hcr$Fref_spt[s,1:nP,t]      <- omFref_sp[s,1:nP]

        if( ctlList$mp$hcr$Bref == "Bmsy" )
          obj$mp$hcr$Bref_spt[s,1:nP,t]      <- amObj$repOpt$Bmsy_sp[1,1:nP]

        if( ctlList$mp$hcr$Bref == "B0" )
          obj$mp$hcr$Bref_spt[s,1:nP,t]      <- 2*amObj$repOpt$Bmsy_sp[1,1:nP]
      }

      # Save convergence info
      obj$mp$assess$pdHess_tsp[pt,,]  <- amObj$sdrep$pdHess
      obj$mp$assess$posSDs_tsp[pt,,]  <- amObj$posSDs
      obj$mp$assess$maxGrad_tsp[pt,,] <- amObj$maxGrad
    }

  }

  # Return object with assessment complete
  return(obj)

} # END .AM_hierProd()

# makeDatParHierProd() 
# Takes catch and indices set up
# for the hierPrd AM, and makes the
# data and parameter lists for
# running the hierProd AM
.makeDatParHierProd <- function(  C_spt, 
                                  I_spft,
                                  ctlList,
                                  mlnUmsy,
                                  lnUmsy_sp,
                                  mBmsy_sp,
                                  omqDevs_spft,
                                  PTm_sp,
                                  nu_spfk,
                                  P1_spf,
                                  P2_spf,
                                  mq_spf,
                                  mq_sf,
                                  sdlnq_f,
                                  oldFleet    = NULL,
                                  oldStock    = NULL,
                                  spCoastwide = FALSE,
                                  fixUmsy     = FALSE,
                                  importQdevs = FALSE )
{
  # Get model dims
  nSS <- dim(C_spt)[1]
  nPP <- dim(C_spt)[2]
  nT  <- dim(C_spt)[3]
  nF  <- dim(I_spft)[3]

  # Set lnq and lnU prior codes
  if( nSS == 1 & nPP == 1 )
  {
    lnUPriorCode <- 0
    lnqPriorCode <- 0
  } else {
    lnqPriorCode <- 1
    lnUPriorCode <- 1
  }

  # Now make init values for pars
  sumCatch_sp <- apply( X = C_spt, FUN = sum, 
                        MARGIN = c(1,2) )

  # Initial time step of AM and PE devs
  initT_sp  <- array(0, dim = c(nSS,nPP) )
  earliestObs_f <- rep(NA,nF)

  # Switch on solveInitBio if necessary
  solveInitBio_sp <- array(0, dim = c(nSS,nPP))

  # Solve for the first time step
  if( ctlList$mp$assess$spLateStart)
    for( s in 1:nSS)
      for( p in 1:nPP )
      {
        for( f in 1:nF )
          earliestObs_f[f] <- min( which( I_spft[s,p,f,] > 0 ), na.rm = T)
        
        initT_sp[s,p]  <- min(earliestObs_f) - 1

        if( initT_sp[s,p] > 0 & ctlList$mp$assess$spSolveInitBio )
          solveInitBio_sp[s,p] <- which.min(earliestObs_f)
      }



  initPE_sp <- array(ctlList$mp$assess$spInitPE, dim = c(nSS,nPP) )
  tFirstIdx <- ctlList$mp$assess$spYrFirstIdx
  for( s in 1:nSS )
    for( p in 1:nPP)
    {
      if( initPE_sp[s,p] <= initT_sp[s,p] )
        initPE_sp[s,p] <- initT_sp[s,p] + 1

      if( initPE_sp[s,p] < tFirstIdx )
        initPE_sp[s,p] <- tFirstIdx # No need to adjust for zero indexing, since first PE is after first year 

      if( initPE_sp[s,p] <= tFirstIdx & ctlList$mp$assess$spSolveInitBio )
        initPE_sp[s,p] <- tFirstIdx + 1 # Adjust for zero indexing, since the first idx is being used for solution
    }

  initBioCode_sp <-  array(0,dim = c(nSS,nPP))
  initBioCode_sp[initT_sp > 0] <- 1

  lnBinit_vec <- numeric(length = 0)
  if( !ctlList$mp$assess$spSolveInitBio)
    for( s in 1:nSS )
      for( p in 1:nPP )
        if( initBioCode_sp[s,p] == 1)
          lnBinit_vec <- c(lnBinit_vec,log(sum(C_spt[s,p,initT_sp[s,p]:nT])))



  # Now create some vectors to
  # indicate if q is shrunk or not, and if
  # q is conditionally derived or not

  # Add something for when fleets are spread (coastwide)
  condMLEq_f <- rep(0, nF)
  condMLEq_f[ctlList$mp$assess$spCondMLEqFleets] <- 1
  condMLEq_f[-ctlList$mp$assess$spFleets] <- 0

  fleetWeights_f <- ctlList$mp$assess$spFleetWeights

  condMLEobsErr_f <- rep(0,nF)
  condMLEobsErr_f[ctlList$mp$assess$spCondMLEobsErrFleets] <- 1
  condMLEobsErr_f[-ctlList$mp$assess$spFleets] <- 0

  # Shrink q?
  shrinkq_f  <- rep(0, nF) 
  shrinkq_f[ctlList$mp$assess$spShrinkqFleets] <- 1
  shrinkq_f[-ctlList$mp$assess$spFleets] <- 0

  # Time-varying q?
  tvq_f       <- rep(0,nF)
  tvq_f[ctlList$mp$assess$spTVqFleets] <- 1  

  calcIndex_spf <- array(0, dim = c(nSS,nPP,nF) )
  stockq_spf    <- array(0, dim = c(nSS,nPP,nF) )
  speciesq_sf   <- array(0, dim = c(nSS,nF) )
  estObsErr_spf <- array(0, dim = c(nSS,nPP,nF))
  idxLimits_spfj<- array(0, dim = c(nSS,nPP,nF,2))
  fleetq_f      <- rep(0, nF)
  nqDevs        <- 0
  for( s in 1:nSS )
    for( p in 1:nPP )
      for( f in 1:nF )
      {
        if( any( I_spft[s,p,f,] > 0, na.rm = TRUE ) )
          calcIndex_spf[s,p,f] <- 1

        posObs_t <- which(I_spft[s,p,f,] > 0)

        if( length(posObs_t) > 0 )
        {

          if( tvq_f[f] == 1 )
            nqDevs <- nqDevs + (length(posObs_t) - 1)



          idxLimits_spfj[s,p,f,1] <- min(posObs_t - 1)
          idxLimits_spfj[s,p,f,2] <- max(posObs_t - 1)
        }
      }

  nStocks_sf <- apply(X = calcIndex_spf, FUN = sum, MARGIN = c(1,3) )
  for( f in 1:nF )
  {
    for( s in 1:nSS )
    {
      if( nStocks_sf[s,f] > 1 & 
          shrinkq_f[f] == 1 &
          condMLEq_f[f] == 0  )
        stockq_spf[s,,f] <- calcIndex_spf[s,,f]

      if( any(I_spft[s,,f,] > 0, na.rm = T) &
          shrinkq_f[f] == 1 )
        speciesq_sf[s,f] <- 1

      if( condMLEobsErr_f[f] == 0 )
        estObsErr_spf[s,,f] <- calcIndex_spf[s,,f]
    }

    if( any(I_spft[,,f,] > 0, na.rm = T) &
        condMLEq_f[f] == 0 )
      fleetq_f[f] <- 1
  }
     
  # Fix Umsy deviations if asked for
  if( fixUmsy )
  {
    lnUmsyInit_sp <- lnUmsy_sp - mlnUmsy
  } else {
    lnUmsyInit_sp <- array(0, dim = c(nSS,nPP))
  }

  # Time to figure out which fleets have condMLE or estimated
  # catchability, and which have shrinkage priors
  nEstq   <- length(which(condMLEq_f == 0 & fleetq_f == 1))
  nCondq  <- length(which(condMLEq_f == 1))

  nShrinkq    <- sum(speciesq_sf)
  nFreeq      <- length(which(shrinkq_f == 0 & condMLEq_f == 0 & fleetq_f == 1))




  # "Learn to fish" persistent change in q - not currently implemented
  # in this AM
  # logisticq_spf <- array(0, dim = c(nSS,nPP,nF))
  # logisticq_spf[,,(oldFleet %in% ctlList$mp$assess$spLogisticqFleets)] <- 1

  # Replace missing values with -1
  I_spft[is.na(I_spft)] <- -1

  # Prior on observation error variance
  IGtau2alpha   <- ctlList$mp$assess$spIGtau2alpha
  IGtau2beta    <- ctlList$mp$assess$spIGtau2Mode * (ctlList$mp$assess$spIGtau2alpha + 1)


  # Start making data list
  data <- list( I_spft          = I_spft,
                C_spt           = C_spt,
                SigmaPriorCode  = 0,
                condMLEq_f      = condMLEq_f,
                shrinkq_f       = shrinkq_f,
                tvq_f           = tvq_f,
                fleetq_f        = fleetq_f,
                fleetWeights_f  = fleetWeights_f,
                condMLEobsErr_f = condMLEobsErr_f,
                lnqPriorCode    = lnqPriorCode,
                lnUPriorCode    = lnUPriorCode,
                BPriorCode      = 0,
                tauObsPriorCode = ifelse(ctlList$mp$assess$spIGtau2alpha > 0, 1, 0),
                initT_sp        = initT_sp,
                initProcErr_sp  = initPE_sp,
                initIdx_sp      = array(tFirstIdx, dim = c(nSS,nPP)),
                initBioCode_sp  = initBioCode_sp,
                calcIndex_spf   = calcIndex_spf,
                stockq_spf      = stockq_spf,
                speciesq_sf     = speciesq_sf,
                posPenFactor    = ctlList$mp$assess$spPosPenFactor,
                solveInitBio_sp = solveInitBio_sp,
                idxLimits_spfj  = idxLimits_spfj,
                importQdevs     = ifelse(importQdevs,1,0),
                PTm_sp          = PTm_sp,
                nu_spfk         = nu_spfk,
                P1_spf          = P1_spf,
                P2_spf          = P2_spf,
                estExpBio       = ifelse(ctlList$mp$assess$spEstExpBio, 1, 0) )
  

  sumCat_sp <- apply(X = C_spt, FUN = sum, MARGIN = c(1,2) )

  pars <- list( lnBmsy_sp         = log(sumCat_sp),
                lnUmsy            = mlnUmsy,
                lnqFree_spf       = array(0,dim = c(nSS,nPP,nF)),
                lnqShrink_sf      = array(0,dim = c(nSS,nF)),
                lntauspf_vec      = rep(log(ctlList$mp$assess$sptauObs_spf), sum(estObsErr_spf)),
                lnBinit_vec       = lnBinit_vec,
                deltalnq_spf      = array(0, dim = c(nSS,nPP,nF)),
                lntauq_s          = rep(log(ctlList$mp$assess$sptauq),nSS),
                tvlnqDevs_vec     = rep(0,max(1,nqDevs + sum(tvq_f) * nSS * nPP )),
                lntautvqDev       = log(ctlList$mp$assess$sptauqdev),
                # Prior on fleet catchability
                mlnq_sf           = log(mq_sf),
                # Prior on stock specific catchability
                # if freely estimated
                mlnq_spf          = log(mq_spf),
                sdlnq_f           = sdlnq_f,
                epslnUmsy_s       = rep(0,nSS),
                lnsigUmsy         = log(ctlList$mp$assess$spsigU),
                epslnUmsy_sp      = lnUmsyInit_sp,
                lnsigUmsy_s       = rep(log(ctlList$mp$assess$spsigU),nSS),
                mlnUmsy           = mlnUmsy,
                sdlnUmsy          = ctlList$mp$assess$spsdUmsy,
                mBmsy_sp          = mBmsy_sp,
                sdBmsy_sp         = mBmsy_sp*ctlList$mp$assess$spCVBmsy,
                mBinit_sp         = mBmsy_sp/2,
                sdBinit_sp        = (mBmsy_sp/2),
                tau2IGa_f         = rep(IGtau2alpha,nF),
                tau2IGb_f         = rep(IGtau2beta,nF),
                sigma2IG          = c(1,.02),
                deltat            = 1,
                lnsigmaProc       = log(ctlList$mp$assess$spsigmaProc),
                zetaspt_vec       = rep(0,sum(nT - initPE_sp) - nSS * nPP),
                sigmaProcMult_sp  = array(1, dim = c(nSS,nPP)),
                logit_gammaYr     = 0,
                omqDevs_spft      = omqDevs_spft )


  
  qFreeMap <- array( 1:(nSS*nPP*nF), dim = c(nSS,nPP,nF) )
  qFreeMap[,,(shrinkq_f == 1 & condMLEq_f == 0)] <- NA
  # qFreeMap[,,ctlList$mp$assess$spFleetWeights == 0] <- NA


  deltalnqMap <- array( nSS*nPP*nF + 1:(nSS*nPP*nF), dim = c(nSS,nPP,nF) )
  deltalnqMap[stockq_spf == 0] <- NA
  # deltalnqMap[,,ctlList$mp$assess$spFleetWeights == 0] <- NA

  lnqShrinkMap <- array( 2*nSS*nPP*nF + 1:(nSS*nF), dim = c(nSS,nF) )
  lnqShrinkMap[speciesq_sf == 0] <- NA
  # lnqShrinkMap[,ctlList$mp$assess$spFleetWeights == 0] <- NA

  baseMap <- list(  lnqFree_spf   = factor(qFreeMap),
                    lnqShrink_sf  = factor(lnqShrinkMap),
                    deltalnq_spf  = factor(deltalnqMap) )


  # Save
  tmbLists <- list( data = data,
                    pars = pars,
                    map  = baseMap )


  return( tmbLists )

} # END .makeDatParHierProd

# Apply idx based AM - basically
# a smoother of the last n
# survey observations, scaled to
# make a biomass estimate for
# each stock
.AM_idxBased <- function( obj, t )
{
  # Pull control list
  ctlList <- obj$ctlList

  # Get complex dims
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nF  <- obj$om$nF
  tMP <- obj$om$tMP

  # Calculate projection year
  pt <- t - tMP + 1
  pT <- ctlList$opMod$pT

  # Get smoother fleets and number of points
  nPoints   <- ctlList$mp$assess$idxPoints
  idxFleets <- ctlList$mp$assess$idxFleets
  idxScalar <- ctlList$mp$assess$idxScalar

  I_spft    <- obj$mp$data$I_spft
  I_spft[I_spft < 0] <- NA

  sumI_spf  <- apply( X = I_spft, FUN = sum, MARGIN = c(1,2,3), na.rm = T )

  for( s in 1:nS )
    for( p in 1:nP )
    {
      posFleets <- which( sumI_spf[s,p,idxFleets] > 0 )

      # Now get the obs series, and pick out
      # most recent nPoints and calculate the smoother
      fleetObs <- I_spft[s,p,idxFleets[posFleets],1:(t-1)]
      fleetObs <- na.omit(fleetObs)
      nObs <- length(fleetObs)
      fleetObs <- fleetObs[(nObs - nPoints + 1):nObs]

      # Smooth:
      movAvg <- exp(mean(log(fleetObs)))

      # Write the smoothed relative biomass
      # estimate down for all projection years
      # so retro plots look better
      obj$mp$assess$retroSB_tspt[pt:pT,s,p,t] <- movAvg/idxScalar
      obj$mp$assess$retroVB_tspft[pt:pT,s,p,,t] <- movAvg/idxScalar

      
    }

  if( ctlList$mp$hcr$Fref == "Umsy" )
    obj$mp$hcr$Fref_spt[,,t]      <- obj$rp$FmsyRefPts$Umsy_sp

  if( ctlList$mp$hcr$Fref == "Fmsy" )
    obj$mp$hcr$Fref_spt[,,t]      <- obj$rp$FmsyRefPts$Fmsy_sp

  if( ctlList$mp$hcr$Bref == "Bmsy" )
    obj$mp$hcr$Bref_spt[,,t]      <- obj$rp$FmsyRefPts$BeqFmsy_sp  

  if( ctlList$mp$hcr$Bref == "B0" )
    obj$mp$hcr$Bref_spt[,,t]      <- obj$om$B0_sp
  

  return(obj)
} # END .AM_idxBased()

# Apply a perfect info AM - returns
# actual future spawning biomass as the retro
.AM_perfectInfo <- function( obj, t )
{
  # Pull control list
  ctlList <- obj$ctlList

  # Get complex dims
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nF  <- obj$om$nF
  tMP <- obj$om$tMP

  # Calculate projection year
  pt <- t - tMP + 1
  pT <- ctlList$opMod$pT

  obj$mp$hcr$TAC_spt[,,t] <- 0
  tmpObj <- .ageSexOpMod( obj, t )

  for( s in 1:nS )
    for( p in 1:nP )
    {
      obj$mp$assess$retroSB_tspt[pt:pT,s,p,t]     <- tmpObj$om$SB_spt[s,p,t]
      obj$mp$assess$retroVB_tspft[pt:pT,s,p,,t]   <- tmpObj$om$vB_spft[s,p,,t] 
      
    }

  if( ctlList$mp$hcr$Fref == "Umsy" )
    obj$mp$hcr$Fref_spt[,,t]      <- obj$rp$FmsyRefPts$Umsy_sp

  if( ctlList$mp$hcr$Fref == "Fmsy" )
    obj$mp$hcr$Fref_spt[,,t]      <- obj$rp$FmsyRefPts$Fmsy_sp

  if( ctlList$mp$hcr$Bref == "Bmsy" )
    obj$mp$hcr$Bref_spt[,,t]      <- obj$rp$FmsyRefPts$BeqFmsy_sp  

  if( ctlList$mp$hcr$Bref == "B0" )
    obj$mp$hcr$Bref_spt[,,t]      <- obj$om$B0_sp
  

  return(obj)
} # END .AM_perfectInfo()


# .mgmtProc()
# Runs replications of the closed loop
# simulation, applying the simulated mgmt
# procedure to the operating model stocks
.mgmtProc <- function( obj )
{

  ctlList <- obj$ctlList

  # Get model dimensions - refactor
  # so that this is done outside of the
  # initPop function
  nReps <- ctlList$ctl$totReps
  pT    <- ctlList$opMod$pT
  tMP   <- obj$om$tMP
  nT    <- obj$om$nT
  nX    <- obj$om$nX
  nF    <- obj$om$nF
  nA    <- obj$om$nA
  nS    <- obj$om$nS
  nP    <- obj$om$nP
  nL    <- obj$om$nL


  ####---------- THE BLOOOBBBB ----- ####


  # OM lists
  om <- list( iSeed_i       = rep(NA, nReps),
              SB_ispt       = array( NA, dim = c(nReps, nS, nP, nT) ),      # SSB
              B_ispt        = array( NA, dim = c(nReps, nS, nP, nT )),      # Total biomass
              vB_ispft      = array( NA, dim = c(nReps, nS, nP, nF, nT )),  # Vulnerable biomass
              R_ispt        = array( NA, dim = c(nReps, nS, nP, nT) ),      # Rec't
              C_ispft       = array( NA, dim = c(nReps, nS, nP, nF, nT) ),  # Catch by fleet (kt)
              C_ispt        = array( NA, dim = c(nReps, nS, nP, nT) ),      # Total catch (kt)
              F_ispft       = array( NA, dim = c(nReps, nS, nP, nF, nT) ),  # Fishing Mort
              I_ispft       = array( NA, dim = c(nReps, nS+1, nP+1, nF, nT) ),  # Indices (without error)
              E_ipft        = array( NA, dim = c(nReps, nP, nF, nT) ),      # Fishing effort
              q_ispft       = array( NA, dim = c(nReps, nS, nP, nF, nT) ),  # Catchability
              qF_ispft      = array( NA, dim = c(nReps, nS, nP, nF, nT) ),  # F/E
              Rev_ispft     = array( NA, dim = c(nReps, nS, nP, nF, nT ) ), # Revenue 
              landVal_ist   = array( NA, dim = c(nReps, nS, nT) ),          # Landed value
              basePrice_ist = array( NA, dim = c(nReps, nS, nT) ),          # Landed value
              effCost_ipft  = array( NA, dim = c(nReps, nP, nF, nT) )
            )

  om$errors  <- list( omegaR_ispt   = array( NA, dim = c(nReps, nS, nP, nT) ),         # Rec Proc Errors
                      delta_ispft   = array( NA, dim = c(nReps, nS+1, nP+1, nF, nT) ), # Idx obs errors
                      priceDev_ist  = array( NA, dim = c(nReps, nS, nT) )
                    )

  ## ADD AGE/LENGTH SAMPLING ERRORS ## 

  mp <- list( data = NULL, assess = NULL, hcr = NULL )

  mp$data <- list( I_ispft = array( NA, dim = c(nReps, nS+1, nP+1, nF, nT ) ) )    # Biomass indices
  ## ADD AGES AND LENGTH DATA ##

  mp$assess <- list(  retroR_itspt      = array( NA, dim = c(nReps, pT, nS, nP, nT ) ),      # retroR
                      retroSB_itspt     = array( NA, dim = c(nReps, pT, nS, nP, nT ) ),      # retroSB
                      retroVB_itspft    = array( NA, dim = c(nReps, pT, nS, nP, nF, nT ) ),  # retroVB
                      retroq_itspf      = array( NA, dim = c(nReps, pT, nS, nP, nF ) ),      # Retrospective catchability 
                      retroq_itspft     = array( NA, dim = c(nReps, pT, nS, nP, nF, nT ) ),  # Retrospective tv catchability 
                      retrotauObs_itspf = array( NA, dim = c(nReps, pT, nS, nP, nF ) ),      # Retrospective obs error estimates 
                      retroUmsy_itspf   = array( NA, dim = c(nReps, pT, nS, nP ) ),          # Retrospective Umsy 
                      retroBmsy_itspf   = array( NA, dim = c(nReps, pT, nS, nP ) ),          # Retrospective Bmsy 
                      maxGrad_itsp      = array( NA, dim = c(nReps, pT, nS, nP ) ),          # AM max gradient component
                      pdHess_itsp       = array( NA, dim = c(nReps, pT, nS, nP ) ),          # AM pdHessian?
                      posSDs_itsp       = array( NA, dim = c(nReps, pT, nS, nP ) ))          # AM all positive SDs?



  mp$hcr <- list( Fref_ispt     = array( NA, dim = c(nReps, nS, nP, nT ) ),      # Reference F in HCR
                  Bref_ispt     = array( NA, dim = c(nReps, nS, nP, nT ) ),      # Reference biomass (B0 or Bmsy, usually)
                  targetF_ispt  = array( NA, dim = c(nReps, nS, nP, nT) ),       # Target F from HCR
                  LCP_ispt      = array(NA, dim = c(nReps, nS, nP, nT) ),        # Lower control point
                  UCP_ispt      = array(NA, dim = c(nReps, nS, nP, nT) ),        # Upper control point
                  TAC_ispt      = array(NA, dim = c(nReps, nS, nP, nT ) ),       # TAC summed across all fleets
                  propTAC_ispt  = array(NA, dim = c(nReps, nS, nP, nT ) ),       # proportion of TAC for disaggregating pooled TACs
                  TAC_ispft     = array(NA, dim = c(nReps, nS, nP, nF, nT ) ) )  # TAC allocated by fleet


  blob <- list( om = om, mp = mp, ctlList = ctlList,
                rp = vector(mode = "list", length = nReps) )

  blob$goodReps_isp <- array(FALSE, dim = c(nReps,nS,nP) )

  blob$omniObjFun <- list(  objFun_i            = array( NA, dim = nReps),
                            objFun_isp          = array( NA, dim = c(nReps,nS,nP)),
                            Cbar_isp            = array( NA, dim = c(nReps,nS,nP)),
                            totCbar_i           = array( NA, dim = nReps),
                            Csum_i              = array( NA, dim = nReps),
                            barLoDep_isp        = array( NA, dim = c(nReps,nS,nP)),
                            barHiDep_isp        = array( NA, dim = c(nReps,nS,nP)),
                            barProbHiDep_isp    = array( NA, dim = c(nReps,nS,nP)),
                            barProbLoDep_isp    = array( NA, dim = c(nReps,nS,nP)),
                            barInitCatDiff_isp  = array( NA, dim = c(nReps,nS,nP)),
                            closedCount_isp     = array( NA, dim = c(nReps,nS,nP)),
                            barAAV_isp          = array( NA, dim = c(nReps,nS,nP)),
                            barCatDiff_isp      = array( NA, dim = c(nReps,nS,nP)),
                            barEffDiff_i        = array( NA, dim = c(nReps)),
                            barInitEffDiff_i    = array( NA, dim = c(nReps)) )


  ##########################################################
  ######### ------- CLOSED LOOP SIMULATION ------- #########
  ##########################################################

  if(ctlList$opMod$switchCommCatchability)
  {
    set.seed(123)
    newSpeciesIdx <- sample(1:nS,nS,replace = FALSE)

    obj$om$newSpeciesIdx <- newSpeciesIdx
  }

  
  message(" (.mgmtProc) Running feedback loop...\n")

  for( i in 1:nReps )
  {
    simObj <- obj

    # Set rep number and random seed
    simObj$ctlList$opMod$rep <- i
    rSeed <- ctlList$ctl$rSeed + i
    blob$om$iSeed_i[i] <- rSeed
    set.seed( rSeed )

    startTime <- proc.time()[3]    
    t1 <- proc.time()[3]

    # Condition history of simulation model
    # and draw random errors for projection
    simObj <- .condMS3pop( simObj )

    # Calculate effort dynamics parameters
    simObj <- .calcHistEffortDynamics( simObj )
    # Calcualate times for observations
    


    if( ctlList$ctl$omni | ctlList$ctl$perfConF )
    {
      simObj <- .solveProjPop( simObj )
    }
    else if( nT >= tMP )
    {
      for( t in tMP:nT )
      {
        # Try updating the population
        objNext <- try( .updatePop(simObj, t) )

        # Save current rep and time step for
        # later debugging/post-mortem
        blob$maxRep       <- i
        blob$maxTimeStep  <- t

        # Break if objNext fails
        if( class(objNext) == "try-error")
        {
          blob$success <- FALSE
          message( "(.mgmtProc) Error in updatePop for replicate ", i, 
                    " at time step ", t, "\n; saving progress up to here.\n",sep = "")
          break
        } # END objNext conditional

        simObj <- objNext

        message( "(.mgmtProc) Time step t = ", t ," complete\n", sep = "")

      } # END t loop
      gc()
    } # END feedback sim for replicate i

    # Fill blob for the current replicate

    # Operating model
    blob$om$SB_ispt[i,,,]     <- simObj$om$SB_spt
    blob$om$B_ispt[i,,,]      <- simObj$om$B_spt
    blob$om$vB_ispft[i,,,,]   <- simObj$om$vB_spft
    blob$om$R_ispt[i,,,]      <- simObj$om$R_spt
    blob$om$C_ispft[i,,,,]    <- simObj$om$C_spft
    blob$om$C_ispt[i,,,]      <- simObj$om$C_spt
    blob$om$F_ispft[i,,,,]    <- simObj$om$F_spft
    blob$om$I_ispft[i,,,,]    <- simObj$om$I_spft
    blob$om$E_ipft[i,,,]      <- simObj$om$E_pft
    blob$om$q_ispft[i,,,,]    <- simObj$om$q_spft
    blob$om$qF_ispft[i,,,,]   <- simObj$om$qF_spft
    blob$om$Rev_ispft[i,,,,]  <- simObj$om$fleetRev_spft
    blob$om$landVal_ist[i,,]  <- simObj$om$landVal_st
    blob$om$basePrice_ist[i,,]<- simObj$om$basePrice_st
    blob$om$effCost_ipft[i,,,]<- simObj$om$effCost_pft

    # Errors - maybe update simObj structure to match blob here
    blob$om$errors$omegaR_ispt[i,,,]  <- simObj$errors$omegaR_spt
    blob$om$errors$delta_ispft[i,,,,] <- simObj$errors$delta_spft

    # Data
    blob$mp$data$I_ispft[i,,,,]       <- simObj$mp$data$I_spft
    # Add ages and lengths later

    # HCR elements
    blob$mp$hcr$Fref_ispt[i,,,]       <- simObj$mp$hcr$Fref_spt
    blob$mp$hcr$Bref_ispt[i,,,]       <- simObj$mp$hcr$Bref_spt
    blob$mp$hcr$targetF_ispt[i,,,]    <- simObj$mp$hcr$targetF_spt
    blob$mp$hcr$LCP_ispt[i,,,]        <- simObj$mp$hcr$LCP_spt
    blob$mp$hcr$UCP_ispt[i,,,]        <- simObj$mp$hcr$UCP_spt
    blob$mp$hcr$TAC_ispt[i,,,]        <- simObj$mp$hcr$TAC_spt
    blob$mp$hcr$propTAC_ispt[i,,,]    <- simObj$mp$hcr$propTAC_spt
    blob$mp$hcr$TAC_ispft[i,,,,]      <- simObj$mp$hcr$TAC_spft


    # Retrospective AM outputs
    blob$mp$assess$retroR_itspt[i,,,,]        <- simObj$mp$assess$retroR_tspt
    blob$mp$assess$retroSB_itspt[i,,,,]       <- simObj$mp$assess$retroSB_tspt
    blob$mp$assess$retroVB_itspft[i,,,,,]     <- simObj$mp$assess$retroVB_tspft
    blob$mp$assess$retroq_itspf[i,,,,]        <- simObj$mp$assess$retroq_tspf
    blob$mp$assess$retroq_itspft[i,,,,,]      <- simObj$mp$assess$retroq_tspft
    blob$mp$assess$retrotauObs_itspf[i,,,,]   <- simObj$mp$assess$retrotauObs_tspf
    blob$mp$assess$retroUmsy_itsp[i,,,]       <- simObj$mp$assess$retroUmsy_tsp
    blob$mp$assess$retroBmsy_itsp[i,,,]       <- simObj$mp$assess$retroBmsy_tsp
    blob$mp$assess$maxGrad_itsp[i,,,]         <- simObj$mp$assess$maxGrad_tsp
    blob$mp$assess$pdHess_itsp[i,,,]          <- simObj$mp$assess$pdHess_tsp
    blob$mp$assess$posSDs_itsp[i,,,]          <- simObj$mp$assess$posSDs_tsp

    # Determine if a good replicate
    probPosSD_sp  <- apply( X = simObj$mp$assess$posSDs_tsp, FUN = mean, MARGIN = c(2,3))
    probPDHess_sp <- apply( X = simObj$mp$assess$pdHess_tsp, FUN = mean, MARGIN = c(2,3))

    if(ctlList$mp$assess$method %in% c("idxBased","PerfectInfo"))
    {
      blob$goodReps_isp[i,,] <- TRUE
    } else{
      for( s in 1:nS )
        for( p in 1:nP )
          if( ( probPosSD_sp[s,p] >= .95 & probPDHess_sp[s,p] >= .95 ) )
          {
            blob$goodReps_isp[i,s,p] <- TRUE 
          }
    }

    if( ctlList$ctl$omni | ctlList$ctl$perfConF )
    {
      blob$goodReps_isp[i,,] <- TRUE

      blob$omniObjFun$objFun_isp[i,,]         <- simObj$objFun_sp
      blob$omniObjFun$objFun_i[i]             <- simObj$objFun       
      blob$omniObjFun$Cbar_isp[i,,]           <- simObj$Cbar_sp
      blob$omniObjFun$totCbar_i[i]            <- simObj$totCbar
      blob$omniObjFun$Csum_i[i]               <- simObj$Csum
      blob$omniObjFun$barLoDep_isp[i,,]       <- simObj$barLoDep_sp
      blob$omniObjFun$barHiDep_isp[i,,]       <- simObj$barHiDep_sp
      blob$omniObjFun$barProbLoDep_isp[i,,]   <- simObj$barProbLoDep_sp
      blob$omniObjFun$barProbHiDep_isp[i,,]   <- simObj$barProbHiDep_sp
      blob$omniObjFun$barInitCatDiff_isp[i,,] <- simObj$barInitCatDiff_sp
      blob$omniObjFun$closedCount_isp[i,,]    <- simObj$closedCount_sp

      blob$omniObjFun$barAAV_isp[i,,]         <- simObj$barAAV_sp
      blob$omniObjFun$barCatDiff_isp[i,,]     <- simObj$barCatDiff_sp
      blob$omniObjFun$barEffDiff_i[i]         <- simObj$barEffDiff
      blob$omniObjFun$barInitEffDiff_i[i]     <- simObj$barInitEffDiff

    }

    # Save reference points for this replicate - more necessary when
    # we have multiple conditioning assessment posterior dist draws
    blob$rp[[i]] <- simObj$rp

    message( " (.mgmtProc) Completed replicate ", i, " of ", nReps, ".\n", sep = "")

    doneGoodReps_sp <- apply( X = blob$goodReps_isp, FUN = sum, MARGIN = c(2,3))

    if( all( doneGoodReps_sp >= ctlList$ctl$nGoodReps ) )
    {
      message(  " (.mgmtProc) Completed ",  ctlList$ctl$nGoodReps, 
                " replicates with all convergent AMs, ending simulation.\n", sep = "")
      break
    } else {
      finishedGood <- min(doneGoodReps_sp)
      message( " (.mgmtProc) Completed ", finishedGood, " of ", ctlList$ctl$nGoodReps, 
                " replicates with all convergent AMs.")
    }
    # Clear memory of extraneous garbage
    gc()
    message( " (.mgmtProc) Maximum replicate count reached, ending simulations.\n")
  } # END i loop
  gc()

  blob$nSims <- i

  # Name the blob array dimensions...
  tMP   -> blob$om$tMP
  nT    -> blob$om$nT
  nX    -> blob$om$nX
  nF    -> blob$om$nF
  nA    -> blob$om$nA
  nS    -> blob$om$nS
  nP    -> blob$om$nP
  nL    -> blob$om$nL

  blob$om$speciesNames  <- dimnames(ctlList$opMod$histRpt$SB_spt)[[1]]
  blob$om$stockNames    <- dimnames(ctlList$opMod$histRpt$SB_spt)[[2]]
  blob$om$fleetNames    <- dimnames(ctlList$opMod$histRpt$vB_spft)[[3]]




  return(blob)
} # END .mgmtProc()

# .solveProjPop()
# Solves for the fishing mortality or effort rates
# to maximise yield and keep fisheries open. Optimises
# based on a given objective function.
# inputs: obj = simulation object
# outputs: obj = solved simulation object with all time 
#                 steps filled
# usage:  under omniscient manager sims, or for no fishing.
# NOTE: need some objectives input via simCtlFile
.solveProjPop <- function( obj )
{
  # Pull control list items
  om      <- obj$om
  ctlList <- obj$ctlList
  mp      <- ctlList$mp
  # Model dimensions
  nS      <- om$nS
  nP      <- om$nP
  nT      <- om$nT
  tMP     <- om$tMP
  pT      <- ctlList$opMod$pT

  # runModel()
  # Private function to run operating model and return
  # objective function value based on depletion, catch
  # and closures.
  # inputs: pars = vector of F or E values
  #         obj = simulation object
  # outputs: outList = list of  
  #                      - obj = simulation object
  #                      - f = objective function value
  # Modified from SP Cox's function
  runModel <- function( pars, obj )
  {
    # Pull control list items
    om      <- obj$om
    ctlList <- obj$ctlList
    mp      <- ctlList$mp
    opMod   <- ctlList$opMod

    # Some omni control pars
    maxF        <- mp$omni$maxF
    maxE        <- mp$omni$maxE
    minE        <- mp$omni$minE
    parMult     <- mp$omni$expParMult
    Fmsy_sp     <- obj$rp$FmsyRefPts$Fmsy_sp
    loDepBmsy   <- mp$omni$loDepBmsy
    hiDepBmsy   <- mp$omni$hiDepBmsy
    maxAAV      <- mp$omni$maxAAV
    maxEffDiff  <- mp$omni$maxRelEffDiff
    maxCatDiff  <- mp$omni$maxRelCatDiff
    fixEff      <- mp$omni$fixEff
    avgEffYrs   <- mp$omni$avgEffYrs


    # Obj function weights
    avgCatWt      <- mp$omni$avgCatWt
    totCatWt      <- mp$omni$totCatWt
    hiDepBmsyWt   <- mp$omni$hiDepBmsyWt
    loDepBmsyWt   <- mp$omni$loDepBmsyWt
    closedWt      <- mp$omni$closedWt
    AAVWt         <- mp$omni$AAVWt
    catDiffWt     <- mp$omni$catDiffWt
    sumCatWt      <- mp$omni$sumCatWt
    effDiffWt     <- mp$omni$effDiffWt
    initCatDiffWt <- mp$omni$initCatDiffWt
    initEffDiffWt <- mp$omni$initEffDiffWt
    probDepWt     <- mp$omni$probDepWt
    totProfitWt   <- mp$omni$totProfitWt

    # Model dimensions
    nS      <- om$nS
    nP      <- om$nP
    nT      <- om$nT
    tMP     <- om$tMP
    pT      <- ctlList$opMod$pT

    rp      <- obj$rp


    # Some pars that are global to either effort
    # simulation
    projInt <- nT - tMP + 1
    nKnots  <- mp$omni$nKnots
    space   <- projInt / (nKnots - 1)
    knotPts <- round( seq(from = 1, by = space, length.out = nKnots) )
    knotPts[knotPts > projInt] <- projInt

    # Create x points for interpolation spline
    x <- knotPts

    # Set fishing mortality to 0 for tMP:nT
    # so that operating model does not hang
    obj$om$F_spft[,,,tMP:nT] <- 0

    # Choose whether this is effort or fishing mort.
    # If perfect targeting (all TACs fully utilised)
    if( ctlList$opMod$effortMod == "targeting" )
    {
      # This is a projection model that always uses
      # a multiple of Fmsy
      if( ctlList$ctl$perfConF )
      {
        for( s in 1:nS )
          for( p in 1:nP )
            obj$om$F_spft[s,p,2,tMP:nT] <- parMult*rp$FmsyRefPts$Fmsy_sp[s,p]

        if( !is.null(ctlList$omni$inputF) )
          obj$om$F_spft[,,2,tMP:nT] <- ctlList$omni$inputF
      } else {  
        # Then we are solving for stock/species specific
        # Fs (nS * nP at each time step)
        tmpFmult <- exp( pars )
        knotF_spk <- array( tmpFmult, dim = c(nS,nP,nKnots) )

        # Make spline for each species/stock
        for( s in 1:nS )
          for( p in 1:nP )
          {
            splineF <- spline( x = x, y = knotF_spk[s,p,], n = projInt)$y
            splineF <- splineF * rp$FmsyRefPts$Fmsy_sp[s,p]
            splineF[splineF < 0] <- 0
            splineF[splineF > maxF] <- maxF

            obj$om$F_spft[s,p,2,tMP:nT] <- parMult*splineF
          }
      }
    }



    # If bycatch is simulated (single area effort catches all species)
    if( ctlList$opMod$effortMod %in% c("dynModel","Max") )
    {
      # This is a projection model that always uses
      # Fmsy
      if( ctlList$ctl$perfConF )
      {
        for( p in 1:nP )
          obj$om$E_pft[p,2,tMP:nT] <- parMult*rp$EmsyMSRefPts$EmsyMS_p[p]

        if( !is.null(ctlList$omni$inputE) )
          for( p in 1:nP )
            obj$om$E_pft[p,2,tMP:nT] <- ctlList$omni$inputE_p[p]

        for( s in 1:nS )
          for( p in 1:nP )
            obj$om$F_spft[s,p,2,tMP:nT] <- obj$om$E_pft[p,2,tMP:nT] * om$qF_spft[s,p,2,tMP - 1]

      } else {

        if( fixEff )
        {
          # Get fixed total effort
          totalEff_t  <- apply( X = obj$om$E_pft[,2,(tMP-avgEffYrs):(tMP - 1)],
                             MARGIN = 2, FUN = sum )
          totalEff    <- mean(totalEff_t)

          # Make an array for holding pars
          knotEps_pk <- array( pars, dim = c(nP,nKnots))

          splineEps_pt  <- array(0, dim = c(nP,(nT - tMP + 1) ) )
          splineE_pt    <- array(0, dim = c(nP,(nT - tMP + 1) ) )
          for( p in 1:nP )
            splineEps_pt[p,] <- spline( x = x, y = knotEps_pk[p,], n = projInt)$y

          splineEps_pt <- exp(splineEps_pt)
          for( t in 1:projInt )
            splineEps_pt[,t] <- splineEps_pt[,t] / sum(splineEps_pt[,t])

          obj$om$E_pft[,2,tMP:nT] <- splineEps_pt * totalEff

        }

        if( !fixEff )
        {
          # Then we are solving for area specific
          # Es (nP at each time step)
          tmpEmult    <- exp( pars )
          knotE_pk    <- array( tmpEmult, dim = c(nP,nKnots) )


          # Make spline for each species/stock
          for( p in 1:nP )
          {
            # Make spline, multiply by Emsy
            splineE <- spline( x = x, y = knotE_pk[p,], n = projInt)$y

            if( mp$omni$baseEffort == "MSY")
              baseEff <- obj$rp$EmsyMSRefPts$EmsyMS_p[p]

            if( mp$omni$baseEffort == "MEY")
              baseEff <- obj$rp$EmeyRefPts$Emey_p[p]

            splineE <- splineE * baseEff

            maxE <- mp$omni$maxE * obj$rp$EmsyMSRefPts$EmsyMS_p[p]

            # Correct for under and over efforts, limit
            # lowest effort to be the minimum over the
            # historical period
            if( minE == "hist" )
              splineE[splineE < 0] <- min(obj$om$E_pft[p,2,1:(tMP-1)])/rp$EmsyMSRefPts$EmsyMS_p[p]
            else 
              splineE[splineE < 0] <- minE
            
            splineE[splineE > maxE] <- maxE

            if(any(is.na(splineE)))
              browser()

            # Assign to projection
            obj$om$E_pft[p,2,tMP:nT] <- parMult * splineE


          }
        }

        for( s in 1:nS )
          for( p in 1:nP)
            obj$om$F_spft[s,p,2,tMP:nT] <- obj$om$E_pft[p,2,tMP:nT] * om$qF_spft[s,p,2,tMP:nT]
      }

    }

    

    

    # run model for projection period
    for(t in tMP:nT )
    {
      obj <- .ageSexOpMod(obj, t)
    }

    # Pull projection catch and biomass
    Cproj_spt   <- obj$om$C_spft[,,2,tMP:nT]
    TACproj_spt <- obj$om$TAC_spft[,,2,tMP:nT]
    Bproj_spt   <- obj$om$B_spt[,,tMP:nT]
    Bmsy_sp     <- rp$FmsyRefPts$BeqFmsy_sp
    
    # Make arrays to hold stock specific
    # penalties and intermediate values
    # Lots of options here from development of the
    # best objective function
    barLoDep_sp         <- array(0, dim = c(nS,nP))
    barHiDep_sp         <- array(0, dim = c(nS,nP))
    probLoDep_sp        <- array(0, dim = c(nS,nP))
    probHiDep_sp        <- array(0, dim = c(nS,nP))
    barProbLoDep_sp     <- array(0, dim = c(nS,nP))
    barProbHiDep_sp     <- array(0, dim = c(nS,nP))
    barAAV_sp           <- array(0, dim = c(nS,nP))
    closedCount_sp      <- array(0, dim = c(nS,nP))
    barCatDiff_sp       <- array(0, dim = c(nS,nP))
    barInitCatDiff_sp   <- array(0, dim = c(nS,nP))

    # Exponential penalty functions
    expLoDep_sp         <- array(0, dim = c(nS,nP))
    expHiDep_sp         <- array(0, dim = c(nS,nP))
    expProbLoDep_sp     <- array(0, dim = c(nS,nP))
    expProbHiDep_sp     <- array(0, dim = c(nS,nP))
    expAAV_sp           <- array(0, dim = c(nS,nP))
    expCatDiff_sp       <- array(0, dim = c(nS,nP))
    expInitCatDiff_sp   <- array(0, dim = c(nS,nP))

    # Linear penalty functions
    linLoDep_sp         <- array(0, dim = c(nS,nP))
    linHiDep_sp         <- array(0, dim = c(nS,nP))
    linProbLoDep_sp     <- array(0, dim = c(nS,nP))
    linProbHiDep_sp     <- array(0, dim = c(nS,nP))
    linAAV_sp           <- array(0, dim = c(nS,nP))
    linCatDiff_sp       <- array(0, dim = c(nS,nP))
    linInitCatDiff_sp   <- array(0, dim = c(nS,nP))



    # Get minimum catch in historical period
    histCatch_spt <- apply( X = om$C_spft[,,,1:(tMP-1)], FUN = sum, MARGIN = c(1,2,4))
    Cmin_sp <- apply( X = histCatch_spt, FUN = min, MARGIN = c(1,2),
                      na.rm = T )

    # Private function for calculating AAV
    .calcAAV <- function( C )
    {
      diffC <- diff(C)
      absDiffC <- abs(diffC)

      AAV <- sum(absDiffC) / sum(C)
      AAV
    }


    # Catch difference
    catDiff_tsp <- abs(apply( X = obj$om$C_spft[,,2,(tMP-1):nT],
                              FUN = diff, MARGIN = c(1,2) ) )
    catDiff_spt <- aperm( catDiff_tsp, c(2,3,1))

    catDiffRel_spt <- catDiff_spt / (obj$om$C_spft[,,2,(tMP:nT) - 1])
    catDiffRel_spt[!is.finite(catDiffRel_spt)] <- 1

    initCatDiffRel_sp <- catDiffRel_spt[,,1]

    # Effort difference (prefer to keep total effort similar)
    totEff_t   <- apply( X = obj$om$E_pft[,2,(tMP-1):nT],
                          FUN = sum, MARGIN = 2 )
    effDiff_t  <- abs(diff(totEff_t))

    effDiffRel_t <- effDiff_t / totEff_t[-pT]
    effDiffRel_t[!is.finite(effDiffRel_t)] <- 1

    initEffDiffRel <- effDiffRel_t[1]

    # Now loop over species/stocks and
    # calculate penalties
    for( s in 1:nS )
      for( p in 1:nP )
      {
        # biomass above some depletion level
        barLoDep_sp[s,p] <- combBarrierPen(  x = Bproj_spt[s,p,]/Bmsy_sp[s,p],
                                              eps = loDepBmsy,
                                              alpha = loDepBmsy / 2,
                                              above = TRUE )

        expLoDep_sp[s,p]  <- expPenaltyFunction(  x = Bproj_spt[s,p,]/Bmsy_sp[s,p],
                                                  eps = loDepBmsy,
                                                  alpha = 1,
                                                  above = TRUE )

        linLoDep_sp[s,p]  <- linPenaltyFunction(  x = Bproj_spt[s,p,]/Bmsy_sp[s,p],
                                                  eps = loDepBmsy,
                                                  alpha = 4,
                                                  above = TRUE,
                                                  beta = mp$omni$linBeta )

        barHiDep_sp[s,p] <- combBarrierPen( x = Bproj_spt[s,p,]/Bmsy_sp[s,p],
                                            eps = hiDepBmsy,
                                            alpha = hiDepBmsy / 2,
                                            above = FALSE )

        expHiDep_sp[s,p]  <- expPenaltyFunction(  x = Bproj_spt[s,p,]/Bmsy_sp[s,p],
                                                  eps = hiDepBmsy,
                                                  alpha = 1,
                                                  above = FALSE )

        linHiDep_sp[s,p]  <- linPenaltyFunction(  x = Bproj_spt[s,p,]/Bmsy_sp[s,p],
                                                  eps = hiDepBmsy,
                                                  alpha = .1,
                                                  above = FALSE,
                                                  beta = mp$omni$linBeta )

        # What about the probability of being between
        # lowDep and HiDep?
        probHiDep_sp[s,p]   <- sum(Bproj_spt[s,p,]/Bmsy_sp[s,p] > hiDepBmsy) / pT
        probLoDep_sp[s,p]   <- sum( Bproj_spt[s,p,]/Bmsy_sp[s,p] > loDepBmsy) / pT

        barProbHiDep_sp[s,p] <- combBarrierPen( probHiDep_sp[s,p],
                                                eps = mp$omni$hiDepProb,
                                                alpha = mp$omni$hiDepProb,
                                                above = FALSE )

        barProbLoDep_sp[s,p] <- combBarrierPen( probLoDep_sp[s,p],
                                                eps = mp$omni$loDepProb,
                                                alpha = mp$omni$loDepProb,
                                                above = TRUE )


        expProbLoDep_sp[s,p]  <- expPenaltyFunction(  x = probLoDep_sp[s,p],
                                                      eps = loDepBmsy,
                                                      alpha = 10,
                                                      above = TRUE )        

        expProbHiDep_sp[s,p]  <- expPenaltyFunction(  x = probHiDep_sp[s,p],
                                                      eps = hiDepBmsy,
                                                      alpha = 10,
                                                      above = FALSE )     


        linProbLoDep_sp[s,p]  <- linPenaltyFunction(  x = probLoDep_sp[s,p],
                                                      eps = loDepBmsy,
                                                      alpha = 10,
                                                      above = TRUE,
                                                      beta = mp$omni$linBeta )        

        linProbHiDep_sp[s,p]  <- linPenaltyFunction(  x = probHiDep_sp[s,p],
                                                      eps = hiDepBmsy,
                                                      alpha = 10,
                                                      above = FALSE,
                                                      beta = mp$omni$linBeta )     

        # Catch less than historical minimum (closures)
        closedCount_sp[s,p] <- sum( Cproj_spt[s,p,] < Cmin_sp[s,p] )

        # Penalty on AAV
        AAV <- .calcAAV( Cproj_spt[s,p,] )
        barAAV_sp[s,p] <- combBarrierPen( x = AAV, 
                                          eps = maxAAV, 
                                          alpha = maxAAV/2,
                                          above = FALSE )

        expAAV_sp[s,p] <- expPenaltyFunction( x = AAV,
                                              eps = maxAAV,
                                              alpha = 10,
                                              above = FALSE )    

        linAAV_sp[s,p] <- linPenaltyFunction( x = AAV,
                                              eps = maxAAV,
                                              alpha = 2,
                                              above = FALSE,
                                              beta = mp$omni$linBeta )   


        barCatDiff_sp[s,p] <- combBarrierPen( x = catDiffRel_spt[s,p,],
                                              eps = maxCatDiff, 
                                              alpha = maxCatDiff/2,
                                              above = FALSE )

        expCatDiff_sp[s,p] <- expPenaltyFunction( x = catDiffRel_spt[s,p,],
                                                  eps = maxCatDiff,
                                                  alpha = 1,
                                                  above = FALSE )     

        linCatDiff_sp[s,p] <- linPenaltyFunction( x = catDiffRel_spt[s,p,],
                                                  eps = maxCatDiff,
                                                  alpha = 1,
                                                  above = FALSE,
                                                  beta = mp$omni$linBeta ) 

        
        barInitCatDiff_sp[s,p] <- combBarrierPen( x = initCatDiffRel_sp[s,p],
                                                  eps = maxCatDiff, 
                                                  alpha = 1,
                                                  above = FALSE )

        expInitCatDiff_sp[s,p] <- expPenaltyFunction( x = initCatDiffRel_sp[s,p],
                                                      eps = maxCatDiff,
                                                      alpha = 2, 
                                                      above = FALSE )

        linInitCatDiff_sp[s,p] <- linPenaltyFunction( x = initCatDiffRel_sp[s,p],
                                                      eps = maxCatDiff,
                                                      alpha = 2, 
                                                      above = FALSE,
                                                      beta = mp$omni$linBeta )
      }


    barEffDiff <- combBarrierPen( x       = effDiffRel_t,
                                  eps     = maxEffDiff,
                                  alpha   = maxEffDiff/2, 
                                  above   = FALSE )

    expEffDiff <- expPenaltyFunction( x     = effDiffRel_t,
                                      eps   = maxEffDiff,
                                      alpha = 1,
                                      above = FALSE )

    linEffDiff <- linPenaltyFunction( x     = effDiffRel_t,
                                      eps   = maxEffDiff,
                                      alpha = 2,
                                      beta  = mp$omni$linBeta,
                                      above = FALSE )


    barInitEffDiff <- combBarrierPen( x       = initEffDiffRel,
                                      eps     = maxEffDiff,
                                      alpha   = maxEffDiff/2, 
                                      above   = FALSE )

    expInitEffDiff <- expPenaltyFunction( x     = initEffDiffRel,
                                          eps   = maxEffDiff,
                                          alpha = 1,
                                          above = FALSE )

    linInitEffDiff <- linPenaltyFunction( x     = initEffDiffRel,
                                          eps   = maxEffDiff,
                                          alpha = 2,
                                          beta  = mp$omni$linBeta,
                                          above = FALSE )

    # Calculate average catch
    Csum    <- sum(Cproj_spt,na.rm = T)
    Cbar_sp <- apply( X = Cproj_spt, FUN = mean, MARGIN = c(1,2))
    totCbar <- mean( apply(X = Cproj_spt, FUN = sum, MARGIN = 3 ) )

    # Total obj function for each stock/species
    objFun_sp <- -  avgCatWt * log(1e3*Cbar_sp) + 
                    (closedWt * closedCount_sp)^mp$omni$linBeta

    if( mp$omni$penType == "barrier" )
    {
      objFun_sp <- objFun_sp +
                    loDepBmsyWt * barLoDep_sp + 
                    hiDepBmsyWt * barHiDep_sp +
                    AAVWt * barAAV_sp +  
                    catDiffWt * barCatDiff_sp + 
                    initCatDiffWt * barInitCatDiff_sp + 
                    probDepWt * barProbLoDep_sp +
                    probDepWt * barProbHiDep_sp

    }

    if( mp$omni$penType == "exponential" )
    {
      objFun_sp <- objFun_sp +
                    loDepBmsyWt * expLoDep_sp + 
                    hiDepBmsyWt * expHiDep_sp +
                    AAVWt * expAAV_sp +  
                    catDiffWt * expCatDiff_sp + 
                    initCatDiffWt * expInitCatDiff_sp + 
                    probDepWt * expProbLoDep_sp +
                    probDepWt * expProbHiDep_sp
    }

    if( mp$omni$penType == "linear" )
    {
      objFun_sp <- objFun_sp +
                    loDepBmsyWt * linLoDep_sp + 
                    hiDepBmsyWt * linHiDep_sp +
                    AAVWt * linAAV_sp +  
                    catDiffWt * linCatDiff_sp + 
                    initCatDiffWt * linInitCatDiff_sp + 
                    probDepWt * linProbLoDep_sp +
                    probDepWt * linProbHiDep_sp
    }

    econYield_pt  <- array(0, dim = dim(Cproj_spt)[2:3] )
    RevProj_pt    <- array(0, dim = dim(Cproj_spt)[2:3] )
    effCost_pt    <- array(0, dim = dim(Cproj_spt)[2:3] )
                    
    # Now calculate economic yield
    RevProj_spt <- obj$om$fleetRev_spft[,,2,tMP:nT]
    effCost_pt  <- obj$om$effCost_pft[,2,tMP:nT]
    for( p in 1:nP )
    {
      
      RevProj_pt[p,]    <- apply( X = RevProj_spt[,p,], FUN = sum, MARGIN = 2,
                                  na.rm = TRUE )
      econYield_pt[p,]  <- RevProj_pt[p,] * (1 - opMod$crewShare) - effCost_pt[p,]

      econYield_pt[p,]  <- econYield_pt[p,] * (1 + opMod$discountRate)^(-1 * (tMP:nT - tMP))
    }

    discEconYield <- sum(econYield_pt)

   
    objFun    <-  sum(objFun_sp) -
                  totCatWt * log(1e3*totCbar) -
                  sumCatWt * log(1e3*Csum) -
                  totProfitWt * log(1e3*discEconYield)


    if( mp$omni$penType == "barrier" )
      objFun <- objFun + 
                  effDiffWt * barEffDiff +
                  initEffDiffWt * barInitEffDiff
    
    if( mp$omni$penType == "exponential" )
      objFun <- objFun + 
                  effDiffWt * expEffDiff +
                  initEffDiffWt * expInitEffDiff   

    if( mp$omni$penType == "linear" )
      objFun <- objFun + 
                  effDiffWt * linEffDiff +
                  initEffDiffWt * linInitEffDiff

    if(is.na(objFun))
      browser()

    # Return results
    result                    <- obj
    result$objFun_sp          <- objFun_sp
    result$objFun             <- objFun
    result$Cbar_sp            <- Cbar_sp
    result$totCbar            <- totCbar
    result$Csum               <- Csum   
    result$closedCount_sp     <- closedCount_sp

    if( mp$omni$penType == "barrier")
    {
      result$barLoDep_sp        <- barLoDep_sp
      result$barHiDep_sp        <- barHiDep_sp
      result$barProbHiDep_sp    <- barProbHiDep_sp
      result$barProbLoDep_sp    <- barProbLoDep_sp
      result$barAAV_sp          <- barAAV_sp
      result$barCatDiff_sp      <- barCatDiff_sp
      result$barEffDiff         <- barEffDiff
      result$barInitCatDiff_sp  <- barInitCatDiff_sp
      result$barInitEffDiff     <- barInitEffDiff
    }

    if( mp$omni$penType == "exponential" )
    {
      result$barLoDep_sp        <- expLoDep_sp
      result$barHiDep_sp        <- expHiDep_sp
      result$barProbHiDep_sp    <- expProbHiDep_sp
      result$barProbLoDep_sp    <- expProbLoDep_sp
      result$barAAV_sp          <- expAAV_sp
      result$barCatDiff_sp      <- expCatDiff_sp
      result$barEffDiff         <- expEffDiff  
      result$barInitCatDiff_sp  <- expInitCatDiff_sp 
      result$barInitEffDiff     <- expInitEffDiff
    }
    
    if( mp$omni$penType == "linear" )
    {
      result$barLoDep_sp        <- linLoDep_sp
      result$barHiDep_sp        <- linHiDep_sp
      result$barProbHiDep_sp    <- linProbHiDep_sp
      result$barProbLoDep_sp    <- linProbLoDep_sp
      result$barAAV_sp          <- linAAV_sp
      result$barCatDiff_sp      <- linCatDiff_sp
      result$barEffDiff         <- linEffDiff  
      result$barInitCatDiff_sp  <- linInitCatDiff_sp
      result$barInitEffDiff     <- linInitEffDiff
    }

    result
  } # END runModel

  # getObjFunctionVal
  # Purpose:        Private function to run operating model and extract objective function
  # Parameters:     pars=log-Fs for t=tMP,...nT
  # Returns:        Objective function f=barrier penalty for biomass less than Bmsy plus
  #                 log(cumulative catch) and penalty for deviation from Fmsy. 
  #                 Multiplier 1000 is for desired precision
  # Source:         S.P. Cox
  getObjFunctionVal <- function( pars, obj ){
    # Function to run asOM and return objective function
    val <- runModel( pars, obj )$objFun

    
    val
  }     # END function getObjFunctionVal


  # -------- Optimisation ------- #
  # Create vector of initial parameters - need to determine
  # if it's nP * nKnots or nS * nP * nKnots
  nPars <- nP * ctlList$mp$omni$nKnots
  if( ctlList$opMod$effortMod == "targeting" )
    nPars <- nS * nPars

  # Pull average historical effort from fleet 2 and average Fmsy
  histEbar  <- mean(  obj$om$E_pft[,2,1:(tMP - 1)],
                      na.rm = T )

  Fmsybar   <- mean( obj$rp$FmsyRefPts$Fmsy_sp)


  # Might need to refine this later, but for now we
  # use a random vector
  initPars <- rnorm(n = nPars, sd = .3)

  # # Multiply by the appropriate scalar
  # if( ctlList$opMod$effortMod %in% c("Max","dynModel") )
  #   initPars <- log( histEbar * exp(initPars) )

  # if( ctlList$opMod$effortMod %in% c("targeting") )
  #   initPars <- log( Fmsybar * exp(initPars) )


  # OK, run optimisation
  if( ctlList$ctl$omni )
  {


    message( " (.solveProjPop) Running omnsicient manager optimisation, please wait.\n")

    initObjFunVal <- round(getObjFunctionVal( pars = initPars, obj = obj), 3)

    message(" (.solveProjPop) Initial objective function value f =", initObjFunVal, ".\n")


    # Find pars that give the best catch/AAV and minimise
    # closures/low catch
    opt <- optim( par = initPars, fn = getObjFunctionVal,
                  method = "BFGS", obj = obj,
                  control=list(maxit=10000, reltol=0.01 ) )
                  # control=list(maxit=3000, reltol=0.001,ndeps=c(.01,.01) ) )

    # opt <- optim( par = opt$par, fn = getObjFunctionVal,
    #               method = "Nelder-Mead", obj = obj,
    #               control=list(maxit=3000, reltol=0.001,ndeps=c(.01,.01) ) )

    message( " (.solveProjPop) Optimisation for omniscient manager completed with f = ", 
              round(opt$value,2), ".\n")
    optPars <- opt$par

  }

  if( ctlList$ctl$perfConF )
  {
    # Create a dummy optPars object,
    # and replace all future Fs with Fmsy in
    # runModel
    optPars <- initPars 
  }

  # rerun OM with optimised parameters
  message( " (.solveProjPop) Running OM without estimation method.\n")
  obj <- runModel( pars = optPars, obj = obj )

  return(obj)
} # END .solveProjPop()

# barrierPen
# Purpose:        Imposes an increasingly large penalty as a quantity x 
#                 approaches a pre-determined barrier. Approaches from above
#                 and below are treated differently depending on user preference
#                 e.g., for keeping Bt > Bmsy one would use above=TRUE, but for 
#                 AAV < 50%, one would use below=TRUE
# Parameters:     x - vector (or scalar) quantity to test, xBar is the barrier, above
#                 determines from which direction the penalty is applied
# Returns:        scalar penalty value
# Source:         S.P. Cox
barrierPen <- function( x, xBar, above=TRUE)
{
  tmp <- rep(0,length(x))
  if( above )
  {   
    # Penalty grows quadratically below barrier xBar
    tmp[ x <= xBar ]  <- ( x[x<=xBar]-xBar )^2 
  
    # Penalty grows to Inf as Bt approaches Bmsy from above
    tmp[ x > xBar ]   <- (-1.)*log( x[x>xBar] - xBar )

  }
  else
  {
    # Penalty grows quadratically above barrier xBar
    tmp[ x >= xBar ]  <- ( x[x>=xBar]-xBar )^2 

    # Penalty grows to Inf as Bt approaches Bmsy from below
    tmp[ x < xBar ]   <- (-1.)*log( xBar - x[x<xBar] )
  }
  bar <- sum(tmp)
  return( bar )
}


# expPenaltyFunction()
# A simple exponential penalty function
expPenaltyFunction <- function( x, eps, 
                                alpha,
                                above = TRUE )
{
  tmp <- rep(0, length(x))

  if( above )
  {
    tmp <- exp( alpha * ( -x + eps ) )
  }

  if( !above )
    tmp <- exp( alpha * ( x - eps ) )

  pen <- sum(tmp)
  pen
}

# linPenaltyFunction()
# A simple exponential penalty function
linPenaltyFunction <- function( x, eps, 
                                alpha,
                                beta = 1,
                                above = TRUE )
{
  tmp <- rep(0, length(x))

  if( above )
  {
    tmp[x <= eps] <- alpha * (-x[x <= eps] + eps)
  }

  if( !above )
    tmp[x >= eps] <- alpha * (x[x >= eps] - eps)
  pen <- sum(abs(tmp)^beta)
  pen
}

# combBarrierPen()
# Imposes a convex combination of a 
# logarithmic barrier penalty and a classic
# linear penalty as x approaches a pre-determined
# barrier value epsilon, defining the edge of 
# the "feasible" region. 
# Inputs:     x = Quantity to be penalised
#             eps = penalised value of x
#             above = boolean with TRUE => keep above
#             alpha = coefficient for logarithmic barrier
#                       penalty component PsiB in combination.
#                       Also controls coefficient alpha/eps of
#                       linear component PsiP. Guarantees
#                       continuity in first derivative.
# Notes:      Higher values of alpha will increase the weight
#             of the logarithmic component in the feasible region,
#             while at the same time increasing the weight of the 
#             linear penalty in the infeasible region through 
#             the coeff alpha/eps
# Returns:    Scalar penalty summed over components of x
# Source:     SDN Johnson
# Reference:  Srinivasan et al, 2008 
combBarrierPen <- function( x, eps, 
                            above = TRUE,
                            alpha = .05 )
{
  if(eps == 0)
  {
    eps <- eps + 1
    x   <- x + 1
  }

  tmp   <- rep(0,length(x))
  PsiB  <- rep(0,length(x))
  PsiP  <- rep(0,length(x))
  if( !above )
  {
    # If keeping below, translate down
    # by 2eps and apply penalty to translated x...
    xprime <- x - 2 * eps

    PsiB[ xprime <= - eps] <- -log( - xprime [xprime <= -eps] / eps ) 

    PsiP[ xprime > - eps ] <- xprime[ xprime > - eps ] + eps

    tmp <- alpha * PsiB  + alpha/eps * PsiP
  }



  if( above )
  {

    PsiB[ - x <= -eps] <- -log( x[-x <= -eps] / eps ) 

    PsiP[ - x > - eps ] <- -x[ -x > -eps ] + eps

    tmp <- alpha * PsiB  + alpha/eps * PsiP
  }



  # plot( x = range(x), y = range(tmp,PsiB,PsiP) )
  #   lines( x = x, y = tmp, lty = 1, lwd = 2, col = "darkgreen")
  #   lines( x = x, y = PsiB, col = "red", lwd = 2, lty = 2 )
  #   lines( x = x, y = PsiP, col = "darkblue", lwd = 2, lty = 3 )

  pen <- sum(tmp)
  return(pen)
}




# .solveMaxEffort()
# Uses species catchability to solve for
# the maximum effort allowed in a stock area
# before any species TAC is completely caught.
# inputs: simObj = simulation object
#         t = time step
# outputs: simObj = simulation object with
#                   projected effort
# Source: SDN Johnson
.solveMaxEffort <- function(  TAC_spf   = TAC_spft[,,,t], 
                              vB_spf    = vB_spft[,,,t],
                              qF_spf    = qF_spft[,,,t],
                              vB_axspf  = vB_axspft[,,,,,t],
                              sel_axspf = sel_axspft[,,,,,t],
                              N_axsp    = N_axspt[,,,,t],
                              M_xsp     = M_xsp,
                              A_s       = A_s,
                              lastE_pf  = E_pft[,,t-1],
                              wt_axsp   = meanWtAge_axsp,
                              nS = nS, nP = nP, nX = nX, nF = nF )
{

  TAC_f <- apply(X = TAC_spf, FUN = sum, MARGIN = 3)
  whichFleets <- which(TAC_f > 0)
  # First, solve Baranov equation
  
  # baranovSol <- .solveBaranov(  C_spf       = TAC_spf, 
  #                               vB_axspf    = vB_axspf,
  #                               vB_spf      = vB_spf,
  #                               N_axsp      = N_axsp,
  #                               sel_axspf   = sel_axspf,
  #                               M_xsp       = M_xsp,
  #                               A_s         = A_s,
  #                               wt_axsp     = wt_axsp,
  #                               nS = nS, nP = nP, nX = nX, nF = nF,
  #                               nIter       = 10,
  #                               baranovStep = 0.2  )


  # targetedF_spf <- baranovSol$F_spf
  # targetedE_spf <- targetedF_spf/qF_spf
  # targetedE_spf[targetedF_spf == 0] <- 0

  # targE_pf <- apply( X = targetedE_spf, FUN = min, MARGIN = c(2,3))
  # E_pf <- targE_pf

  # Now loop over areas and slowly increase effort
  # in each area until the TAC of one species is met
  C_spf <- array( 0, dim = c(nS,nP,nF) )
  E_pf  <- array( 0, dim = c(nP,nF) )
  solveE_s <- rep(0,nS)
  for( p in 1:nP )
  {
    for( f in whichFleets )
    {
      catDiff_s <- C_spf[,p,f] - TAC_spf[,p,f]
      catRem_s  <- TAC_spf[,p,f] - C_spf[,p,f]
      propE     <- 0
      nSteps    <- 0 

      checkCatDiff <- all(abs(catDiff_s) > 1e-2)

      while( checkCatDiff )
      {
        # Solve for the effort using the TAC and catchability
        solveE_s <- -1/qF_spf[,p,f] * log( 1 - catRem_s / (vB_spf[,p,f] - C_spf[,p,f]))

        propE <- propE + 0.8*min(solveE_s[solveE_s >= 0],na.rm = T)

        propF_s <- propE * qF_spf[,p,f]

        # Now use this to generate a catch
        catArrays <- .applyBaranovCatchEq(  M_xsp = M_xsp[,,p,drop = FALSE],
                                            sel_axspf = sel_axspf[,,,p,f,drop = FALSE],
                                            vB_axspf = vB_axspf[,,,p,f,drop = FALSE],
                                            F_spf = array(propF_s, dim = c(nS,1,1)),
                                            wt_axsp = wt_axsp[,,,p,drop = FALSE],
                                            nS = nS,
                                            nP = 1,
                                            nF = 1,
                                            nX = nX,
                                            A_s = A_s )

        # Save catch
        C_spf[,p,f] <- catArrays$C_spf

        catDiff_s <- C_spf[,p,f] - TAC_spf[,p,f]
        catRem_s  <- TAC_spf[,p,f] - C_spf[,p,f]

        nSteps <- nSteps + 1
        if( nSteps > 5 )
        {
          break

        }

        while( any(catRem_s < 0 ) )
        {
          # Put in routine to wind back E a little
          propE <- 0.9 * propE

          propF_s <- propE * qF_spf[,p,f]

          # Now use this to generate a catch
          catArrays <- .applyBaranovCatchEq(  M_xsp = M_xsp[,,p,drop = FALSE],
                                              sel_axspf = sel_axspf[,,,p,f,drop = FALSE],
                                              vB_axspf = vB_axspf[,,,p,f,drop = FALSE],
                                              F_spf = array(propF_s, dim = c(nS,1,1)),
                                              wt_axsp = wt_axsp[,,,p,drop = FALSE],
                                              nS = nS,
                                              nP = 1,
                                              nF = 1,
                                              nX = nX,
                                              A_s = A_s )

          # Save catch
          C_spf[,p,f] <- catArrays$C_spf

          catDiff_s <- C_spf[,p,f] - TAC_spf[,p,f]
          catRem_s  <- TAC_spf[,p,f] - C_spf[,p,f]

          break

        }

        

      }
      
      E_pf[p,f] <- propE
    }
  }


  out <- list( E_pf = E_pf )

} # END .solveMaxEffort()




# .solveProjEffortDynamics()
# Uses catchability, price, and a utility
# function to allocate effort among areas.
# inputs: simObj = simulation object
#         t = time step
# outputs: simObj = simulation object with
#                   projected effort
# Source: SDN Johnson
.solveProjEffortDynamics <- function( TAC_spf   = TAC_spft[,,,t], 
                                      vB_spf    = vB_spft[,,,t],
                                      qF_spf    = qF_spft[,,,t],
                                      # vB_axspf  = vB_axspft[,,,,,t],
                                      # sel_axspf = sel_axspft[,,,,,t],
                                      # N_axsp    = N_axspt[,,,,t],
                                      # M_xsp     = M_xsp,
                                      # A_s       = A_s,
                                      # wt_axsp   = meanWtAge_axsp,
                                      price_s   = obj$om$price_s,
                                      alphaU    = obj$om$alphaU,
                                      ut_50     = obj$om$ut_50,
                                      ut_95     = obj$om$ut_95,
                                      w_pf      = obj$om$w_pf,
                                      nS = nS, nP = nP, nX = nX, nF = nF )
{
  # First, we need to get whichFleets
  TAC_f <- apply( X = TAC_spf, FUN = sum, MARGIN = 3 )
  whichFleets <- which(TAC_f > 0)

  # For a given vector of effort, we need to
  # calculate the expected revenue - private fn
  # for this purpose
  calcRev <-  function( logE_p,
                        fleetIdx = 2, 
                        vB_spf,
                        qF_spf,
                        TAC_spf,
                        alphaU,
                        price_s,
                        ut_50,
                        ut_95,
                        w_pf,
                        runOptim = TRUE )
  {
    # Exponentiate effort
    E_p <- exp(logE_p)
    # Model dims
    nS  <- dim(vB_spf)[1]
    nP  <- dim(vB_spf)[2]
    # Now compute expected revenue
    vB_sp   <- vB_spf[,,fleetIdx]
    qF_sp    <- qF_spf[,,fleetIdx]
    TAC_sp  <- TAC_spf[,,fleetIdx]

    # estimate catch and calculate utility
    C_sp        <- array(0, dim = c(nS,nP))
    U_sp        <- array(1, dim = c(nS,nP))
    rev_sp      <- array(0, dim = c(nS,nP))
    TACrem_sp   <- array(0, dim = c(nS,nP))
    # utility function pars
    utAlpha     <- -1 / (log(.6) - log(1.3))
    utBeta      <- log(.6) / (log(.6) - log(1.3))

    for( s in 1:nS )
      for(p in 1:nP )
      {
        C_sp[s,p]    <- vB_sp[s,p] * ( 1 - exp(-qF_sp[s,p] * E_p[p]) )

        TACrem_sp[s,p] <- 1 - C_sp[s,p]/TAC_sp[s,p]

        # U_sp[s,p] <- 1 / (1 + exp( - (TACrem_sp - ut_50)(log(19) * (ut95 - ut50))))
        # U_sp[s,p] <- utAlpha * log( TACrem_sp[s,p] + .3 ) + utBeta

        if( TACrem_sp[s,p] < .2 )
          U_sp[s,p] <- - alphaU / (TACrem_sp[s,p])
        if( TACrem_sp[s,p] < 0 )
          U_sp[s,p] <- - alphaU / 1e-9 
        
        rev_sp[s,p]  <- C_sp[s,p] * price_s[s] * U_sp[s,p] * w_pf[p,1]
      }

    if(runOptim)
    {
      totalRev <- sum(rev_sp)
      return(totalRev)  
    }

    return(rev_sp)
    
  }

  # Now optimise
  initE_p <- c(10,4,10)
  E_pf    <- array(0, dim = c(nP, nF) )
  rev_spf <- array(0, dim = c(nS, nP, nF) )
  
  for( fIdx in whichFleets)
  {

    tmp <- optim( par = log(initE_p), fn = calcRev,
                  control = list( fnscale = -1 ), 
                  method = "Nelder-Mead", 
                  fleetIdx = fIdx,
                  vB_spf = vB_spf,
                  qF_spf = qF_spf,
                  TAC_spf = TAC_spf,
                  ut_50 = ut_50,
                  ut_95 = ut_95,
                  alphaU = alphaU,
                  w_pf  = w_pf,
                  price_s = price_s )


    E_pf[,fIdx]     <- exp(tmp$par)

    rev_spf[,,fIdx] <- calcRev( tmp$par,
                                fleetIdx = fIdx,
                                vB_spf = vB_spf,
                                qF_spf = qF_spf,
                                TAC_spf = TAC_spf,
                                ut_50 = ut_50,
                                ut_95 = ut_95,
                                alphaU = alphaU,
                                price_s = price_s,
                                w_pf  = w_pf,
                                runOptim = FALSE )



  }

  outList <- list(  E_pf = E_pf,
                    rev_spf = rev_spf )

  return(outList)

}

# .calcHistEffortDynamics()
# Calculates stock/area specific effort
# dynamics parameters w_p, which are the
# "attractiveness" weightings for each
# area/stock.
# inputs:   simObj
# outputs:  simObj
# Reference: Walters and Bonfil, 1999 (CJFAS)
# Source: SDN Johnson
.calcHistEffortDynamics <- function( obj )
{
  message(" (.calcEffortDynamics) Calculating fleet effort dynamics model pars from history.\n" )
  # OK, get model dims
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nT  <- obj$om$nT
  nF  <- obj$om$nF
  tMP <- obj$om$tMP

  ctlList <- obj$ctlList

  # Count which fleets are fishing
  allocYears <- ctlList$opMod$allocYears
  allocYearCatch <- obj$om$C_spft[,,1:2,(tMP - allocYears):(tMP - 1)]
  totCatch_f <- apply(X = allocYearCatch, FUN = sum, MARGIN = c(3))
  whichFleets <- which( totCatch_f > 0 )

  w_pf <- array(1,dim = c(nP,length(whichFleets)))

  # Save an array to the om object for area weights
  obj$om$w_pf <- array(0, dim = c(nP,nF) )

  if( nP == 1 )
  {
    obj$om$w_pf[, whichFleets] <- w_pf
    
    return(obj)
  }

  # Get historical effort
  histdx <- 1:(tMP - 1)
  histE_pft <- obj$om$E_pft[,,histdx]
  histE_pft[histE_pft < 0] <- NA

  # We only need effort from the fleets that are
  # actually going to be fishing, i.e. trawl.mod for
  # DERPA runs
  histE_pft   <- histE_pft[,whichFleets,,drop = FALSE]

  # Now calculate mean historical effort - hopefully this
  # stands up to dims of size 1
  Ebar_pf <- apply( X = histE_pft, FUN = mean,
                              MARGIN = c(1,2), na.rm = T )
  totEbar_f <- apply( X = Ebar_pf, FUN = sum, 
                              MARGIN = 2 )

  # Now compute mean historical profitability of each area
  # using historical qs and vulnerable biomasses
  price_s <- obj$om$price_s
  vB_spft <- obj$om$vB_spft[,,,histdx,drop = FALSE]
  qF_spft <- obj$om$qF_spft[,,,histdx,drop = FALSE]
  C_spft  <- obj$om$C_spft[,,,histdx,drop = FALSE]


  # Calculate profitability array
  P_spft  <- array( NA, dim = dim(vB_spft) )
  for( s in 1:nS )
    P_spft[s,,,] <- vB_spft[s,,,] * qF_spft[s,,,] * price_s[s]

  P_pft  <- apply( X = P_spft, FUN = sum, MARGIN = c(2,3,4),
                    na.rm = T )

  Pbar_pf <- apply( X = P_pft, FUN = mean, MARGIN = c(1,2),
                    na.rm = T )
  
  Pbar_pf <- Pbar_pf[,whichFleets,drop = FALSE]


  for( i in 1:5 )
  {
    # Multiply weight and profit
    wP_pf     <- w_pf * Pbar_pf
    w_pf_tmp  <- array(0, dim = dim(w_pf)) 
    # Sum for each fleet's total profit
    wPsum_f <- apply( X = wP_pf, FUN = sum, MARGIN = 2)

    for( fIdx in 1:length(whichFleets))
    {
      for( p in 1:nP)
      {
        w_pf_tmp[p, fIdx] <- ( wPsum_f[fIdx] - wP_pf[p,fIdx]) / 
                              (Pbar_pf[p,fIdx] *( totEbar_f[fIdx] / Ebar_pf[p,fIdx] - 1) )
      }
    }
    w_pf <- 0.5 * ( w_pf + w_pf_tmp )
  }

  obj$om$w_pf[, whichFleets] <- w_pf

  return(obj)
}


# .condMS3pop()
# Conditions OM for historical period from a hierSCAL report
# object.
.condMS3pop <- function( obj )
{

  ctlList <- obj$ctlList

  if(!obj$ctlList$ctl$quiet)
    message(" (.condMS3pop) Conditioning ms3R from hierSCAL report\n")

  # Get model historical period report object
  repObj <- ctlList$opMod$histRpt

  # Get model dims
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nT  <- obj$om$nT
  nF  <- obj$om$nF
  tMP <- obj$om$tMP
  pT  <- ctlList$opMod$pT

  # Now, copy model fit
  histdx <- 1:(tMP - 1)
  obj$om$F_spft[,,,histdx]        <- repObj$F_spft
  obj$om$C_spft[,,,histdx]        <- repObj$predCw_spft
  obj$om$C_spt[,,histdx]          <- apply( X = repObj$predCw_spft, FUN = sum, MARGIN = c(1,2,4) )
  obj$om$E_pft[,1:2,histdx]       <- repObj$E_pft[,1:2,]
  obj$om$sel_axspft[,,,,,histdx]  <- repObj$sel_axspft
  obj$om$sel_lspft[,,,,histdx]    <- repObj$sel_lspft

  obj$mp$hcr$TAC_spft[,,,histdx]  <- obj$om$C_spft[,,,histdx]
  obj$mp$hcr$TAC_spt[,,histdx]    <- obj$om$C_spt[,,histdx]

  # If unstandardised commercial CPUE is provided,
  # then read it in and replace report object
  if( !is.null(ctlList$opMod$commCPUEdata) )
  {
    fileList <- list.files( ctlList$opMod$commCPUEdata)
    fileList <- fileList[grepl("DER.csv",fileList)]
    

    cpueTables <- lapply( X = file.path(ctlList$opMod$commCPUEdata,fileList), 
                          FUN = read.csv, header = TRUE,
                          stringsAsFactors = FALSE )

    areas <- stringr::str_split(fileList, "_")

    pullFirst <- function(x) {x[[1]]}

    areas <- unlist(lapply(X = areas, FUN = pullFirst))

    names(cpueTables) <- areas

    # Need to make sure years match
    fYear <- ctlList$opMod$fYear
    histYears <- seq(fYear,by = 1, length.out = tMP-1)
    
    histFleetYrs <- which(histYears <= 1995)
    modFleetYrs  <- which(histYears > 1995 )

    # Now make a dummy effort by species arra
    E_spft <- array(NA, dim = c(nS,nP,nF,tMP-1))

    speciesNames  <- c("Dover","English","Rock")
    stockNames    <- c("HSHG","QCS","WCVI")

    for( s in 1:nS )
      for( p in 1:nP )
      {
        unstdCPUE <- cpueTables[[stockNames[p]]]

        cpueHist <- unstdCPUE[histFleetYrs,speciesNames[s]]
        cpueMod  <- unstdCPUE[modFleetYrs,speciesNames[s]]

        E_spft[s,p,1,histFleetYrs]  <- repObj$predCw_spft[s,p,1,histFleetYrs] / cpueHist
        E_spft[s,p,2,modFleetYrs]  <- repObj$predCw_spft[s,p,2,modFleetYrs] / cpueMod
      }

    histE_pft <- apply(X = 1e3*E_spft, FUN = mean, MARGIN = c(2,3,4), na.rm = T)
    histE_pft[!is.finite(histE_pft)] <- NA
    histE_pft[histE_pft == 0] <- NA
    obj$om$E_pft[,,histdx] <- histE_pft
  }

  # Now fill in future sel
  obj$om$sel_axspft[,,,,,tMP:nT]  <- repObj$sel_axspft[,,,,,tMP-1]
  obj$om$sel_lspft[,,,,tMP:nT]    <- obj$om$sel_lspft[,,,,tMP-1]

  # Copy index catchability
  obj$om$q_spft[,,,histdx]        <- repObj$q_spft
  obj$om$q_spf                    <- repObj$q_spf

  # Now commercial fishing mortality catchability
  for( s in 1:nS )
    obj$om$qF_spft[s,,,histdx]      <- repObj$F_spft[s,,,] / obj$om$E_pft[,,histdx]

  obj$om$qF_spft[!is.finite(obj$om$qF_spft)] <- 0

  # Fill in future q as fixed at conditioning year's value
  # (OR maybe use mean of historical values? look at TS)
  obj$om$q_spft[,,,tMP:nT]        <- repObj$q_spft[,,,tMP-1]
  
  for( s in 1:nS )
    for( p in 1:nP )
      for( f in 1:nF)
      {
        posIdx <- which(obj$om$qF_spft[s,p,f,1:(tMP-1)] != 0 )
        
        if( length(posIdx) > 0)
        {
          obj$om$qF_spft[s,p,f,tMP:nT]  <- exp(mean(log(obj$om$qF_spft[s,p,f,posIdx])))
        }

      }

  # Switch comm trawl catchability randomly
  if( ctlList$opMod$switchCommCatchability)
  {
    newSpeciesIdx <- obj$om$newSpeciesIdx
    
    for( p in 1:nP )
    {
      obj$om$qF_spft[newSpeciesIdx,p,,tMP:nT] <- obj$om$qF_spft[1:nS,p,,tMP:nT]
    }
  }


  # Now we have enough info to calculate reference points
  stime <- Sys.time()
  message(" (.condMS3pop) Calculating Fmsy and Emsy reference points\n")
  
  repObj$om     <- obj$om
  repObj$opMod  <- ctlList$opMod
  refPtList <- calcRefPts( repObj )
  obj$rp    <- refPtList

  etime <- Sys.time()
  rpTime <- round(etime - stime,2)
  message( paste(" (.condMS3pop) Reference point calculations completed in ", 
                  rpTime, " seconds\n", sep = "" ) )


  # Add historical data
  obj$mp$data$I_spft[1:nS,1:nP,,histdx] <- repObj$I_spft
  obj$om$I_spft[1:nS,1:nP,,histdx]      <- repObj$I_spft
  obj$mp$data$A_axspft[,,,,,histdx]  <- repObj$age_axspft
  obj$mp$data$L_lxspft[,,,,,histdx]  <- repObj$len_lxspft

  # replace negatives with NAs for plotting, can change back 
  # later for TMB
  obj$mp$data$I_spft[obj$mp$data$I_spft<0] <- NA
  obj$mp$data$A_axspft[obj$mp$data$A_axspft<0] <- NA
  obj$mp$data$L_lxspft[obj$mp$data$L_lxspft<0] <- NA

  # Now calculate total catch and allocation
  tdxRecent <- (tMP - ctlList$opMod$allocYears):(tMP - 1)
  recentCatch_spf <- apply( X = obj$om$C_spft[,,,tdxRecent], FUN = sum, MARGIN = c(1,2,3) )

  commGears <- ctlList$opMod$commGears


  # Now, check if we're using correlated
  # errors - should refactor the following two routines to
  # avoid so much repetition - uncorrelated devs can be drawn
  # with an ident. mtx as sigma
  if( ctlList$opMod$corrRecDevs )
  {
    # Put historical rec devs into a tMP-1 x (nS x nP) matrix
    histRecDevs_spt <- repObj$omegaR_spt[,,1:(tMP - 1)] 
    histRecDevMat   <- array( 0, dim = c(tMP-1,nS * nP) )
    for( s in 1:nS )
      histRecDevMat[,(s-1)*nP + 1:nP] <- t(histRecDevs_spt[s,,])

    

    histRecDevMat[histRecDevMat==0] <- NA

    histRecDevCorr <- cor(histRecDevMat, use = "pairwise")

    # Now draw correlated deviates
    corrRecDevsMat <- mvtnorm::rmvnorm( n = nT, mean = rep(0,nS*nP),
                                        sigma = histRecDevCorr,
                                        method = "chol")

    for( s in 1:nS )
      for( p in 1:nP)
        obj$errors$omegaR_spt[s,p,] <- corrRecDevsMat[,(s-1)*nP + p ]
  }

  # Loop and fill errors and allocation
  for( s in 1:nS )
    for( p in 1:nP )
    {
      obj$errors$delta_spft[s,p,,] <- rnorm(nF*nT)

      if( !ctlList$ctl$noProcErr & !ctlList$opMod$corrRecDevs )
        obj$errors$omegaR_spt[s,p,] <- rnorm(nT)

      # Save historical proc errors, but use simulated recruitments after
      # last estimated recruitment
      lastDevIdx <- max(which(repObj$omegaR_spt[s,p,] != 0) )

      obj$errors$omegaR_spt[s,p,histdx[1:lastDevIdx]]   <- repObj$omegaR_spt[s,p,1:lastDevIdx] 


      
      if( !ctlList$ctl$noProcErr )
        obj$errors$omegaR_spt[s,p,histdx[1:lastDevIdx]] + 0.5*repObj$sigmaR_sp[s,p]  # rec devs    


      # Might need a different function for this. Also, skip zeroes.
      obj$om$alloc_spf[s,p,commGears] <- recentCatch_spf[s,p,commGears] / sum( recentCatch_spf[s,p,commGears])
    }


  obj$errors$omegaRinit_asp         <- repObj$omegaRinit_asp # Initialisation errors
  
  # Separate obs errors for combined indices - might need a switch...
  obj$errors$delta_spft[nS+1,1:nP,,] <- rnorm(nT * nP * nF)
  obj$errors$delta_spft[1:(nS+1),nP+1,,] <- rnorm(nT * (nS + 1)*nF)

  # Save historical errors
  if( ctlList$mp$data$source == "cond")
    obj$errors$delta_spft[1:nS,1:nP,,histdx]  <- repObj$residCPUE_spft # obs errors
   
  
  # Then add obsErrorMultiplier for different precision scenarios
  obj$errors$obsErrMult_spft        <- array(1, dim = c(nS,nP,nF,nT))

  # Adjust obs error multiplier if 
  # projObsErrMult != 1
  if( ctlList$opMod$projObsErrMult != 1.0 )
  {
    nPhaseIn        <- min(ctlList$opMod$phaseInObsErrMult,pT-1)

    projObsErrMult  <- ctlList$opMod$projObsErrMult
    for( k in 0:nPhaseIn)
    {
      # Calculate this step's adjustment multiplier
      adjMult <- 1.0 + k / nPhaseIn * (projObsErrMult - 1.0)
      obj$errors$obsErrMult_spft[,,,tMP+k] <-  adjMult * obj$errors$obsErrMult_spft[,,,tMP+k]
    }
    obj$errors$obsErrMult_spft[,,,(tMP+nPhaseIn):nT] <- projObsErrMult
  }

  # Get dimension names and economic/effort model pars
  obj$om$speciesNames   <- dimnames(ctlList$opMod$histRpt$SB_spt)[[1]]
  obj$om$stockNames     <- dimnames(ctlList$opMod$histRpt$SB_spt)[[2]]
  obj$om$fleetNames     <- dimnames(ctlList$opMod$histRpt$vB_spft)[[3]]

  obj$om$price_s        <- ctlList$opMod$price_s[obj$om$speciesNames]
  obj$om$alphaU         <- ctlList$opMod$alphaU
  obj$om$ut_50          <- ctlList$opMod$ut_50
  obj$om$ut_95          <- ctlList$opMod$ut_95

  # Fill in economic parameters
  obj$om$basePrice_st[,tMP] <- ctlList$opMod$price_s[obj$om$speciesNames]

  priceDevCorrMat <- diag(1,nS)
  if(ctlList$opMod$corrPriceDevs)
  {
    priceDevCorrMat <- matrix(ctlList$opMod$priceCorr, nrow = nS, ncol = nS)
    diag(priceDevCorrMat) <- 1
  }

  obj$errors$priceDev_st[1:nS,1:nT] <- t( mvtnorm::rmvnorm( n = nT, mean = rep(0,nS),
                                                            sigma = priceDevCorrMat,
                                                            method = "chol") ) * ctlList$opMod$priceSD

  obj <- .calcTimes( obj )

  message(" (.condMS3pop) Running OM for historical period.\n")
  
  # Now, initialise the population
  for( t in 1:(tMP - 1) )
  {
    obj <- .ageSexOpMod( obj, t )
  }

  return(obj)


} # END condMS3pop()


# .initMS3pop()
# Initialises operating model data structures
# for holding model states, observation and
# process errors, and simulated MP data.
.initMS3pop <- function( ctlList )
{
  repObj <- ctlList$opMod$histRpt
  pT     <- ctlList$opMod$pT

  # make OM list
  om <- list()

  message(" (.initMS3pop) Initialising MS3 operating model.\n")

  # Get model dimensions
  tMP   <- repObj$nT + 1   # start of projections
  nT    <- tMP + pT - 1    # Total time period of sims
  nX    <- repObj$nX       # Number of sex classes
  nF    <- repObj$nF       # Number of fleets
  nA    <- repObj$nA       # Number of Age classs
  nS    <- repObj$nS       # Number of species
  nP    <- repObj$nP       # Number of stocks
  nL    <- repObj$nL       # Number of length classes

  # Get species specific largest age/length classes
  om$A_s <- repObj$A_s
  om$L_s <- repObj$L_s

  # Get length bin aggregation info - might refine later.
  om$lenBinWidth   <- repObj$lenBinWidth   # width of length bins
  om$lenBinMids_l  <- repObj$lenBinMids_l  # Mid points
 

  # Set up arrays to hold simulated states
  om$B_axspt      <- array(0,  dim = c(nA,nX,nS,nP,nT) )      # Biomass at age (total)
  om$N_axspt      <- array(0,  dim = c(nA,nX,nS,nP,nT) )      # Numbers at age
  om$SB_spt       <- array(0,  dim = c(nS,nP,nT) )            # Spawning biomass
  om$B_spt        <- array(0,  dim = c(nS,nP,nT) )            # Spawning biomass
  om$R_spt        <- array(0,  dim = c(nS,nP,nT) )            # Recruitment
  om$Z_axspt      <- array(0,  dim = c(nA,nX,nS,nP,nT) )      # Total mortality

  # Catch, fishing mort, vuln bio, selectivity
  om$C_spft       <- array(NA,  dim = c(nS,nP,nF,nT) )        # Total catch by fleet in kt
  om$C_spt        <- array(NA,  dim = c(nS,nP,nT) )           # Total catch in kt
  om$C_axspft     <- array(0,   dim = c(nA,nX,nS,nP,nF,nT) )  # Numbers
  om$Cw_axspft    <- array(0,   dim = c(nA,nX,nS,nP,nF,nT) )  # Weight
  om$vB_axspft    <- array(0,   dim = c(nA,nX,nS,nP,nF,nT) )  # vuln Bio (granular)
  om$vB_spft      <- array(0,   dim = c(nS,nP,nF,nT) )        # vuln Bio (aggregate)
  om$sel_axspft   <- array(0,   dim = c(nA,nX,nS,nP,nF,nT) )  # sel at age (tv)
  om$sel_lspft    <- array(0,   dim = c(nL,nS,nP,nF,nT) )     # sel at len (tv)
  om$sel_axspf    <- array(0,   dim = c(nA,nX,nS,nP,nF) )     # sel at age
  om$sel_lspf     <- array(0,   dim = c(nL,nS,nP,nF) )        # sel at len
  om$F_spft       <- array(NA,  dim = c(nS,nP,nF,nT) )        # Fishing mortality
  om$I_spft       <- array(NA,  dim = c(nS+1,nP+1,nF,nT) )    # Observations without error
  om$E_pft        <- array(NA,  dim = c(nP,nF,nT) )           # Fishing effort
  om$alloc_spf    <- array(0,   dim = c(nS,nP,nF))            # Catch allocation

  # Observation model pars
  om$q_spft       <- array(0,  dim = c(nS,nP,nF,nT) )         # catchability (tv)
  om$q_spf        <- array(0,  dim = c(nS,nP,nF) )            # catchability

  # Commercial catchability scalar between effort and F
  om$qF_spft      <- array(0,  dim = c(nS,nP,nF,nT) )         # catchability (tv)

  # Economic sub-model pars
  om$landVal_st     <- array(0,   dim = c(nS,nT))               # Landed value (after price flexibility) ($/kg)
  om$basePrice_st   <- array(0,   dim = c(nS,nT))               # Base ex-vessel price ($/kg)
  om$fleetRev_spft  <- array(0,   dim = c(nS,nP,nF,nT) )        # Fleet revenue by area/fleet
  om$effCost_pft    <- array(0,   dim = c(nP,nF,nT) )           # Cost of effort by area/fleet
  om$profit_pft     <- array(0,   dim = c(nP,nF,nT) )           # undiscounted profit



  # Variance parameters
  om$tauObs_spf <- repObj$tauObs_spf * ctlList$opMod$obsErrMult
  om$sigmaR_sp  <- repObj$sigmaR_sp

  # Leading bio pars
  om$B0_sp      <- repObj$B0_sp
  om$M_xsp      <- repObj$M_xsp
  om$h_sp       <- repObj$h_sp

  # Growth model pars
  om$L1_s       <- repObj$L1_s
  om$A1_s       <- repObj$A1_s
  om$L2_xsp     <- repObj$L2_xsp
  om$A2_s       <- repObj$A2_s
  om$vonK_xsp   <- repObj$vonK_xsp
  om$LWa_s      <- repObj$LWa_s
  om$LWb_s      <- repObj$LWb_s

  # Later: write calcSchedules function
  om$matAge_asp         <- repObj$matAge_asp
  om$meanWtAge_axsp     <- repObj$meanWtAge_axsp
  om$Wlen_ls            <- repObj$Wlen_ls
  om$probLenAge_laxsp   <- repObj$probLenAge_laxsp
  om$lenAge_axsp        <- repObj$lenAge_axsp

  # Make data list - anything else?
  # mp list - mostly for retrospective
  # analyses and assessment model settings
  mp <- list( assess = NULL, data = NULL, hcr = NULL )
  mp$assess <- list()
  mp$hcr    <- list()
  mp$data <- list()

  # Data list
  mp$data$I_spft    <- array( NA, dim = c(nS+1,nP+1,nF,nT) )      # Add a pooled data dimension
  mp$data$A_axspft  <- array( NA, dim = c(nA,nX,nS,nP,nF,nT) )    
  mp$data$L_lxspft  <- array( NA, dim = c(nL,nX,nS,nP,nF,nT) ) 

  # assessment list
  mp$assess$retroSB_tspt      <- array( NA, dim = c(pT, nS, nP, nT) )
  mp$assess$retroR_tspt       <- array( NA, dim = c(pT, nS, nP, nT) )
  mp$assess$retroVB_tspft     <- array( NA, dim = c(pT, nS, nP, nF, nT) )
  mp$assess$retroq_tspf       <- array( NA, dim = c(pT, nS, nP, nF) )
  mp$assess$retroq_tspft      <- array( NA, dim = c(pT, nS, nP, nF, nT) )
  mp$assess$retrotauObs_tspf  <- array( NA, dim = c(pT, nS, nP, nF) )
  mp$assess$retroUmsy_tsp     <- array( NA, dim = c(pT, nS, nP ) )
  mp$assess$retroBmsy_tsp     <- array( NA, dim = c(pT, nS, nP ) )
  mp$assess$maxGrad_tsp       <- array( NA, dim = c(pT, nS, nP ) )
  mp$assess$pdHess_tsp        <- array( FALSE, dim = c(pT, nS, nP) )
  mp$assess$posSDs_tsp        <- array( FALSE, dim = c(pT, nS, nP) )

  mp$hcr$LCP_spt              <- array( NA, dim = c(nS, nP, nT) )
  mp$hcr$UCP_spt              <- array( NA, dim = c(nS, nP, nT) )
  mp$hcr$Bref_spt             <- array( NA, dim = c(nS, nP, nT) )
  mp$hcr$Fref_spt             <- array( NA, dim = c(nS, nP, nT) )
  mp$hcr$targetF_spt          <- array( NA, dim = c(nS, nP, nT) )
  mp$hcr$TAC_spft             <- array(0,  dim = c(nS,nP,nF,nT) )     # MP TAC by fleet in kt
  mp$hcr$TAC_spt              <- array(0,  dim = c(nS,nP,nT) )        # MP TAC in kt
  mp$hcr$propTAC_spt          <- array(1,  dim = c(nS,nP,nT) )        # MP TAC in kt


  # Process and observation errors
  errors <- list()

  # Arrays for holding errors
  errors$omegaR_spt       <- array(0, dim = c(nS,nP,nT) )
  errors$omegaRinit_asp   <- array(0, dim = c(nA,nS,nP) )
  errors$delta_spft       <- array(0, dim = c(nS+1,nP+1,nF,nT) )
  errors$priceDev_st      <- array(0, dim = c(nS,nT))




  # Save model dims
  om$tMP   <- tMP
  om$nT    <- nT
  om$nX    <- nX
  om$nF    <- nF
  om$nA    <- nA
  om$nS    <- nS
  om$nP    <- nP
  om$nL    <- nL

  # Return initialised object
  obj <- list(  om = om, 
                mp = mp, 
                errors = errors,
                ctlList = ctlList )

  return(obj)
} # END initMS3pop()


# .applyBaranovCatchEq()
# Refactored to apply Baranov eq. in both the
# effort solver and the operating model.
.applyBaranovCatchEq <- function( M_xsp,
                                  sel_axspf,
                                  vB_axspf,
                                  F_spf,
                                  wt_axsp,
                                  nS,
                                  nP,
                                  nF,
                                  nX,
                                  A_s )
{ 
  # Make arrays to hold catch
  nA <- max(A_s)
  Cw_axspf <- array(0, dim = c(nA,nX,nS,nP,nF) )
  C_axspf  <- array(0, dim = c(nA,nX,nS,nP,nF) )
  C_spf    <- array(0, dim = c(nS,nP,nF) )
  Z_axsp   <- array(0, dim = c(nA,nX,nS,nP))

  # Now generate Z
  for( s in 1:nS )
    for( p in 1:nP )
    {
      # Fill Z with Ms
      Z_axsp[1:A_s[s],1,s,p] <- M_xsp[1,s,p]
      Z_axsp[1:A_s[s],2,s,p] <- M_xsp[2,s,p]

      # Now loop over fleets, ages, and sexes, add fishing mortality
      # rates and generate catch
      for( f in 1:nF )
      {
        Z_axsp[,,s,p] <- Z_axsp[,,s,p] + sel_axspf[,,s,p,f] * F_spf[s,p,f]

        # Now generate catch
        for( a in 1:A_s[s] )
        {
          for( x in 1:nX )
          {
            Cw_axspf[a,x,s,p,f]  <- (1 - exp( - Z_axsp[a,x,s,p])) * vB_axspf[a,x,s,p,f] * F_spf[s,p,f] / Z_axsp[a,x,s,p]
            C_axspf[a,x,s,p,f]   <- Cw_axspf[a,x,s,p,f] / wt_axsp[a,x,s,p]
            C_spf[s,p,f]         <- C_spf[s,p,f] + sum(Cw_axspf[a,x,s,p,f],na.rm = T)
          }
        }
        
      } # END f loop        
    } # END p loop
    # END s loop

  out <- list(  Cw_axspf  = Cw_axspf,
                C_axspf   = C_axspf,
                C_spf     = C_spf )

  return(out)
} # END .applyBaranovCatchEq()

# Operating model function (similar to mseR)
.ageSexOpMod <- function( obj, t )
{
  # Pull om
  om      <- obj$om
  rp      <- obj$rp
  data    <- obj$mp$data
  err     <- obj$errors
  hcr     <- obj$mp$hcr
  mp      <- obj$mp
  opMod   <- obj$ctlList$opMod
  ctlList <- obj$ctlList

  # Get model dimensions, and state variables
  tMP               <- om$tMP
  nT                <- om$nT
  nX                <- om$nX
  nF                <- om$nF
  nA                <- om$nA
  nS                <- om$nS
  nP                <- om$nP
  nL                <- om$nL
  A_s               <- om$A_s
  L_s               <- om$L_s

  # Length bin info
  lenBinWidth       <- om$lenBinWidth
  lenBinMids_l      <- om$lenBinMids_l

  # State vars
  B_axspt           <- om$B_axspt
  N_axspt           <- om$N_axspt
  R_spt             <- om$R_spt
  SB_spt            <- om$SB_spt
  B_spt             <- om$B_spt
  Z_axspt           <- om$Z_axspt
  
  # Catch
  C_axspft          <- om$C_axspft
  Cw_axspft         <- om$Cw_axspft
  F_spft            <- om$F_spft
  E_pft             <- om$E_pft
  C_spft            <- om$C_spft
  C_spt             <- om$C_spt
  TAC_spft          <- hcr$TAC_spft
  TAC_spt           <- hcr$TAC_spt 
  Rev_spft          <- obj$om$fleetRev_spft 

  # Vuln bio
  vB_axspft         <- om$vB_axspft
  vB_spft           <- om$vB_spft

  # Selectivity
  sel_axspft        <- om$sel_axspft
  sel_lspft         <- om$sel_lspft
  sel_axspf         <- om$sel_axspf
  sel_lspf          <- om$sel_lspf

  # Catchability 
  q_spft            <- om$q_spft
  qF_spft           <- om$qF_spft
  q_spf             <- om$q_spf

  # Biomass indices
  Ierr_spft         <- data$I_spft
  Iperf_spft        <- om$I_spft

  # Economic model
  v_st              <- om$landVal_st
  basePrice_st      <- om$basePrice_st
  effCost_pft       <- om$effCost_pft

  
  # Biological pars
  B0_sp             <- om$B0_sp
  M_xsp             <- om$M_xsp
  h_sp              <- om$h_sp
  R0_sp             <- rp$R0_sp
  surv_axsp         <- rp$Surv_axsp
  reca_sp           <- rp$rec.a_sp
  recb_sp           <- rp$rec.b_sp  

  # Growth
  L1_s              <- om$L1_s
  A1_s              <- om$A1_s
  L2_spx            <- om$L2_spx
  A2_s              <- om$A2_s
  vonK_spx          <- om$vonK_spx
  LWa_s             <- om$LWa_s
  LWb_s             <- om$LWb_s
  meanWtAge_axsp    <- om$meanWtAge_axsp
  Wlen_ls           <- om$Wlen_ls
  probLenAge_laxsp  <- om$probLenAge_laxsp
  lenAge_axsp       <- om$lenAge_axsp
  matAge_asp        <- om$matAge_asp


  # Initialise population
  if( t == 1 )
  {
    # Pull initial N multiplier in case
    # of non-eq init
    initNmult_asp <- exp(err$omegaRinit_asp)

    for( s in 1:nS )
      for( p in 1:nP)
      {
        N_axspt[,1,s,p,t] <- surv_axsp[,1,s,p] * R0_sp[s,p] * initNmult_asp[,s,p]
        N_axspt[,2,s,p,t] <- surv_axsp[,2,s,p] * R0_sp[s,p] * initNmult_asp[,s,p]
      }

    R_spt[,,t] <- R0_sp


    # Fill TAC history as catch history

  }

  # Otherwise update using total mortality
  if( t > 1 )
  {
    # Loop over species/pops
    for( s in 1:nS )
      for( p in 1:nP )
        for( a in 1:A_s[s] )
        {
          # Recruitment
          if( a == 1 )
          {
            R_spt[s,p,t] <- reca_sp[s,p] * SB_spt[s,p,t-1] / (1 + recb_sp[s,p] * SB_spt[s,p,t-1])
            
            if( !obj$ctlList$ctl$noProcErr )
              R_spt[s,p,t] <- R_spt[s,p,t] * exp( om$sigmaR_sp[s,p] * err$omegaR_spt[s,p,t] - 0.5*om$sigmaR_sp[s,p]^2) 

            N_axspt[a,,s,p,t] <- R_spt[s,p,t]
          }
          # Apply mortality           
          if( a > 1 )
          {
            N_axspt[a,,s,p,t] <- N_axspt[a-1,,s,p,t-1] * exp( - Z_axspt[a-1,,s,p,t-1])
          }
          # Plus group
          if( a == A_s[s] )
          {
            N_axspt[a,,s,p,t] <- N_axspt[a,,s,p,t] + N_axspt[a,,s,p,t-1] * exp( -Z_axspt[a,,s,p,t-1])
          }

        }
  }

  # Now convert into biomass and spawning biomass
  B_axspt[,,,,t]    <- N_axspt[,,,,t] * meanWtAge_axsp
  # Now compute vuln biomass
  for( f in 1:nF )
    vB_axspft[,,,,f,t] <- sel_axspft[,,,,f,t] * B_axspt[,,,,t]

  vB_spft[,,,t] <- apply( X = vB_axspft[,,,,,t], FUN = sum, MARGIN = 3:5, na.rm = T )
  B_spt[,,t]    <- apply( X = B_axspt[,,,,t], FUN = sum, MARGIN = 3:4, na.rm = T )

  # Now calculate spawning biomass
  SB_asp <- matAge_asp * B_axspt[,nX,,,t]
  SB_spt[,,t] <- apply(X = SB_asp, FUN = sum, MARGIN = c(2,3), na.rm = T )

  # update price data
  if( t > tMP )
  {
    i                 <- opMod$interestRate
    basePrice_st[,t]  <- basePrice_st[,t-1] * (1 + i) * exp( err$priceDev_st[,t] )
    
  }

  # Now calculate effort for each area (for closed loop sim)
  if( t >= tMP & all(is.na(F_spft[,,,t])) )
  {
    if( opMod$effortMod == "Max" )
    {
      effortMod <- .solveMaxEffort( TAC_spf   = TAC_spft[,,,t], 
                                    vB_spf    = vB_spft[,,,t],
                                    qF_spf    = qF_spft[,,,t],
                                    vB_axspf  = vB_axspft[,,,,,t],
                                    sel_axspf = sel_axspft[,,,,,t],
                                    N_axsp    = N_axspt[,,,,t],
                                    M_xsp     = M_xsp,
                                    A_s       = A_s,
                                    wt_axsp   = meanWtAge_axsp,
                                    lastE_pf  = E_pft[,,t-1],
                                    nS = nS, nP = nP, nX = nX, nF = nF )

      E_pft[,,t]     <- effortMod$E_pf
      if( any(E_pft[,,t] < 0))
        browser()
      # Rev_spft[,,,t] <- effortMod$rev_spf


      # Convert to F
      for( s in 1:nS )
        F_spft[s,,,t] <- E_pft[,,t] * qF_spft[s,,,t]
    }


    if( opMod$effortMod == "dynModel" )
    {
      effortMod <- .solveProjEffortDynamics(  TAC_spf   = TAC_spft[,,,t], 
                                              vB_spf    = vB_spft[,,,t],
                                              qF_spf    = qF_spft[,,,t],
                                              # vB_axspf  = vB_axspft[,,,,,t],
                                              # sel_axspf = sel_axspft[,,,,,t],
                                              # N_axsp    = N_axspt[,,,,t],
                                              # M_xsp     = M_xsp,
                                              # A_s       = A_s,
                                              # wt_axsp   = meanWtAge_axsp,
                                              price_s   = om$price_s,
                                              alphaU    = om$alphaU,
                                              ut_50     = om$ut_50,
                                              ut_95     = om$ut_95,
                                              w_pf      = om$w_pf,
                                              nS = nS, nP = nP, nX = nX, nF = nF  )
      E_pft[,,t]     <- effortMod$E_pf
      # Rev_spft[,,,t] <- effortMod$rev_spf

      # Convert to F
      for( s in 1:nS )
        F_spft[s,,,t] <- E_pft[,,t] * qF_spft[s,,,t]

    }

    if( opMod$effortMod == "targeting" )
    {
      # F_spft[,,,t] <- .mfPA(  C_spf     = TAC_spft[,,,t], 
      #                         vB_axspf  = vB_axspft[,,,,,t],
      #                         vB_spf    = vB_spft[,,,t],
      #                         N_axsp    = N_axspt[,,,,t],
      #                         sel_axspf = sel_axspft[,,,,,t],
      #                         M_xsp     = M_xsp,
      #                         A_s       = A_s,
      #                         wt_axsp   = meanWtAge_axsp,
      #                         nS = nS, nP = nP, nX = nX, nF = nF  )

      Flist <- .solveBaranov( C_spf       = TAC_spft[,,,t], 
                              vB_axspf    = vB_axspft[,,,,,t],
                              vB_spf      = vB_spft[,,,t],
                              N_axsp      = N_axspt[,,,,t],
                              sel_axspf   = sel_axspft[,,,,,t],
                              M_xsp       = M_xsp,
                              A_s         = A_s,
                              wt_axsp     = meanWtAge_axsp,
                              nS = nS, nP = nP, nX = nX, nF = nF,
                              nIter       = opMod$baranovIter,
                              baranovStep = opMod$baranovStep  )

      F_spft[,,,t] <- Flist$F_spf
    }
   

  } # END F approximation

  # if( t == nT-5 )
  # {

  #   inputF <- obj$ctlList$mp$omni$inputF

  #   whichF <- which.min(abs(rp$refCurves$F - inputF))

  #   rpSurv_axsp <- rp$refCurves$surv_axspf[,,,,whichF]

  #   currentSurv_axsp <- N_axspt[,,,,t]
  #   for( s in 1:nS )
  #     for( p in 1:nP )
  #       currentSurv_axsp[,,s,p] <- currentSurv_axsp[,,s,p] / R_spt[s,p,t]

  # }

  # Now generate Z
  for( s in 1:nS )
    for( p in 1:nP )
    {
      # Fill Z with Ms
      Z_axspt[1:A_s[s],1,s,p,t] <- M_xsp[1,s,p]
      Z_axspt[1:A_s[s],2,s,p,t] <- M_xsp[2,s,p]

      # Now loop over fleets, ages, and sexes, add fishing mortality
      # rates and generate catch
      for( f in 1:nF )
      {
        Z_axspt[,,s,p,t] <- Z_axspt[,,s,p,t] + sel_axspft[,,s,p,f,t] * F_spft[s,p,f,t]

        # Now generate catch
        for( a in 1:A_s[s] )
        {
          for( x in 1:nX )
          {
            Cw_axspft[a,x,s,p,f,t]  <- (1 - exp( - Z_axspt[a,x,s,p,t])) * vB_axspft[a,x,s,p,f,t] * F_spft[s,p,f,t] / Z_axspt[a,x,s,p,t]
            C_axspft[a,x,s,p,f,t]   <- Cw_axspft[a,x,s,p,f,t] / meanWtAge_axsp[a,x,s,p]
            C_spft[s,p,f,t]         <- C_spft[s,p,f,t] + sum(Cw_axspft[a,x,s,p,f,t],na.rm = T)
          }
        }
        
      } # END f loop        
    } # END p loop
    # END s loop

  # Sum realised catch
  C_spft[,,,t] <- apply( X = Cw_axspft[,,,,,t], FUN = sum, MARGIN = c(3,4,5), na.rm = T )
  C_spt[,,t]   <- apply( X = C_spft[,,,t], FUN = sum, MARGIN = c(1,2), na.rm = T )


  # Calculate landed value per kg given price flexibility
  if( t >= tMP )
  {
    pflex   <- opMod$priceFlex
    C_s     <- apply( X = C_spt[,,t], FUN = sum, MARGIN = 1 )
    MSY_sp  <- rp$FmsyRefPts$YeqFmsy_sp
    MSY_s   <- apply(X = MSY_sp, FUN = sum, MARGIN = 1)
    v_st[,t] <- basePrice_st[,t] * ( 1 + pflex * ( C_s/MSY_s - 1 ) )
    
    for( p in 1:nP)
      effCost_pft[p,,t] <- E_pft[p,,t] * rp$EmeyRefPts$effortPrice_p[p]

    if( t > tMP )
      effCost_pft[,,t] <- effCost_pft[,,t] * (1 + i)^(t - tMP)
  }

  # Calculate revenue
  for( p in 1:nP )
    for( f in 1:nF)
      Rev_spft[,p,f,t] <- C_spft[,p,f,t] * v_st[,t]


  # Now generate indices 
  # if( t >= tMP | ctlList$mp$data$source == "OM" )
  # {
  # Get the indices that we're still collecting in the future
  idxOn_spf <- obj$ctlList$mp$data$idxOn_spft[,,,t]

  # Calculate and save
  for( f in 1:nF)
  {
    # First, the species/stock specific observations
    for( s in 1:nS )
      for( p in 1:nP )
      {
        idxOn <- obj$mp$data$idxOn_spft[s,p,f,t]
        if( !idxOn )
          next
        # Compute precision
        tau <- om$tauObs_spf[s,p,f] * err$obsErrMult_spft[s,p,f,t]

        # relative biomass survey
        if( ctlList$mp$data$idxType[f] == 1)
          Iperf_spft[s,p,f,t] <- q_spft[s,p,f,t] * vB_spft[s,p,f,t]

        # CPUE
        if( ctlList$mp$data$idxType[f] == 2 & C_spft[s,p,f,t] > 0 )
          Iperf_spft[s,p,f,t] <- C_spft[s,p,f,t] / E_pft[p,f,t] 

        # Save true observations and those with error
        Ierr_spft[s,p,f,t] <- Iperf_spft[s,p,f,t] * exp(tau * err$delta_spft[s,p,f,t] - 0.5 * tau^2)
      } 
  }
  # }

  # Now make aggregates
  tauObs_spf  <- om$tauObs_spf * err$obsErrMult_spft[s,p,f,t]
  idxOn_spf   <- obj$mp$data$idxOn_spft[,,,t]
  # Species Pooling
  if( ctlList$mp$data$speciesPooling & !ctlList$mp$data$spatialPooling  )
  {
    for( f in 1:nF)
    {
      for(p in 1:nP )
      {
        tauObs_spf[tauObs_spf == 0] <- NA
        tau <- mean(tauObs_spf[,p,f],na.rm = TRUE)
        if( ctlList$mp$data$idxType[f] == 1 )
          Iperf_spft[nS+1,p,f,t] <- sum( q_spft[,p,f,t] * vB_spft[,p,f,t] * idxOn_spf[,p,f], na.rm = T )

        if( ctlList$mp$data$idxType[f] == 2 )
        {
          Iperf_spft[nS+1,p,f,t] <- max(idxOn_spf[,p,f]) * sum( C_spft[,p,f,t] ) 
          if( Iperf_spft[nS+1,p,f,t] > 0 )
            Iperf_spft[nS+1,p,f,t] <- Iperf_spft[nS+1,p,f,t] / max(1e-6,E_pft[p,f,t]) / 1e3
        }

        # Add error
        Ierr_spft[nS+1,p,f,t] <- Iperf_spft[nS+1,p,f,t] * exp( tau * err$delta_spft[nS+1,p,f,t] - 0.5 * tau^2)

      }
    }

    if(any(is.nan(Iperf_spft)))
      browser()
  }

  # Spatial Pooling
  if( ctlList$mp$data$spatialPooling & !ctlList$mp$data$speciesPooling )
  {
    for( f in 1:nF)
    {
      for( s in 1:nS )
      {
        tauObs_spf[tauObs_spf == 0] <- NA
        tau <- mean(tauObs_spf[s,,f],na.rm = TRUE)
        
        if( ctlList$mp$data$idxType[f] == 1 )
          Iperf_spft[s,nP + 1,f,t] <- sum( q_spft[s,,f,t] * vB_spft[s,,f,t] * idxOn_spf[s,,f], na.rm = T )

        if( ctlList$mp$data$idxType[f] == 2 )
        {
          Iperf_spft[s,nP+1,f,t] <- max(idxOn_spf[s,,f]) * sum(  C_spft[s,,f,t] ) 
          if( Iperf_spft[s,nP+1,f,t] > 0 ) 
            Iperf_spft[s,nP+1,f,t] <- Iperf_spft[s,nP+1,f,t] / max(1e-6,sum( E_pft[,f,t], na.rm = TRUE )) / 1e3
        }

        Ierr_spft[s,nP+1,f,t] <- Iperf_spft[s,nP+1,f,t] * exp( tau * err$delta_spft[s,nP+1,f,t] - 0.5 * tau^2)


      }
    }

    if(any(is.nan(Iperf_spft)))
      browser()
  }

  # Total aggregation
  if( ctlList$mp$data$spatialPooling & ctlList$mp$data$speciesPooling )
  {

    for( f in 1:nF)
    {
      tauObs_spf[tauObs_spf == 0] <- NA
      tau <- mean(tauObs_spf[,,f],na.rm = TRUE)

      if( ctlList$mp$data$idxType[f] == 1 )
        Iperf_spft[nS+1,nP + 1,f,t] <- sum( q_spft[,,f,t] * vB_spft[,,f,t] * idxOn_spf[1:nS,1:nP,f], na.rm = T )

      if( ctlList$mp$data$idxType[f] == 2 )
      {
        Iperf_spft[nS+1,nP+1,f,t] <- max(idxOn_spf[,,f]) * sum( C_spft[,,f,t] ) 
        if( Iperf_spft[nS+1,nP+1,f,t] > 0 )
          Iperf_spft[nS+1,nP+1,f,t] <- Iperf_spft[nS+1,nP+1,f,t] / max(1e-6,sum( E_pft[,f,t], na.rm = TRUE )) / 1e3
      }

      Ierr_spft[nS+1,nP+1,f,t] <- Iperf_spft[nS+1,nP+1,f,t] * exp( tau * err$delta_spft[nS+1,nP+1,f,t] - 0.5 * tau^2)
    }

    if(any(is.nan(Iperf_spft)))
      browser()
  }

  Iperf_spft[Iperf_spft == 0] <- NA
  Ierr_spft[Ierr_spft == 0] <- NA


  # put all the state arrays back into OM
  # State vars
  B_axspt           -> obj$om$B_axspt
  N_axspt           -> obj$om$N_axspt
  R_spt             -> obj$om$R_spt
  SB_spt            -> obj$om$SB_spt
  B_spt             -> obj$om$B_spt
  Z_axspt           -> obj$om$Z_axspt
  
  # Catch
  C_axspft          -> obj$om$C_axspft
  Cw_axspft         -> obj$om$Cw_axspft
  F_spft            -> obj$om$F_spft
  C_spft            -> obj$om$C_spft
  C_spt             -> obj$om$C_spt
  E_pft             -> obj$om$E_pft
  Rev_spft          -> obj$om$fleetRev_spft
  
  # Vuln bio
  vB_axspft         -> obj$om$vB_axspft
  vB_spft           -> obj$om$vB_spft

  # Selectivity
  sel_axspft        -> obj$om$sel_axspft
  sel_lspft         -> obj$om$sel_lspft
  sel_axspf         -> obj$om$sel_axspf
  sel_lspf          -> obj$om$sel_lspf

  # Catchability 
  q_spft            -> obj$om$q_spft
  q_spf             -> obj$om$q_spf
  qF_spft           -> obj$om$qF_spft

  # Biomass indices
  Iperf_spft        -> obj$om$I_spft
  Ierr_spft         -> obj$mp$data$I_spft

  # economic model
  v_st              -> obj$om$landVal_st
  basePrice_st      -> obj$om$basePrice_st
  effCost_pft       -> obj$om$effCost_pft


  return( obj )
} # END ageSexOpMod()


# .solveBaranov()
# Generalised Newton-Rhapson solver or the Baranov
# catch equation, given catch, vulnerable biomass
# and other stuff.
.solveBaranov <- function(  C_spf       = newCat_spf, 
                            vB_axspf    = vB_axspft[,,,,,t],
                            vB_spf      = vB_spft[,,,t],
                            N_axsp      = N_axspt[,,,,t],
                            sel_axspf   = sel_axspft[,,,,,t],
                            M_xsp       = M_xsp,
                            A_s         = A_s,
                            wt_axsp     = meanWtAge_axsp,
                            nS = nS, nP = nP, nX = nX, nF = nF,
                            nIter       = 5,
                            baranovStep = 0.5  )
{
  # Make arrays to hold estimated catch with
  # succesive approximations
  appZ_axsp         <- array(0, dim = dim(N_axsp))
  appF_axspf        <- array(0, dim = dim(vB_axspf) )
  appF_spf          <- array(0, dim = dim(C_spf))

  # Approximated catch
  appC_axspf        <- array(0, dim = dim(vB_axspf) )
  appC_spf          <- array(0, dim = dim(C_spf))

  # Jacobian
  J_spf             <- array(0, dim = dim(C_spf))

  # Initialise F at C/B
  appF_spf   <- C_spf / vB_spf

  # Calculate first iteration of Z
  for( s in 1:nS )
    for( p in 1:nP )
    {
      for( x in 1:nX )
      {
        appZ_axsp[1:A_s[s],x,s,p] <- M_xsp[x,s,p]

        for( f in 1:nF )
          appZ_axsp[1:A_s[s],x,s,p] <- appZ_axsp[1:A_s[s],x,s,p] + appF_spf[s,p,f] * sel_axspf[1:A_s[s],x,s,p,f]
      } 
    }

  # if nIter == 0, return F and Z as is.

  if( nIter > 0 & any( C_spf > 0) )
  {
    # Refine F
    for( i in 1:nIter)
    {
      # Reset objective function
      f_spf <- C_spf

      # Calculate implied catch using new Z and F
      for( s in 1:nS )
        for( p in 1:nP )
        {
          whichF <- which( C_spf[s,p,] > 0)
          for( f in whichF )
          {
            
            appC_axspf[1:A_s[s],,s,p,f] <- vB_axspf[1:A_s[s],,s,p,f] * (1 - appZ_axsp[1:A_s[s],,s,p]) * appF_spf[s,p,f]/appZ_axsp[1:A_s[s],,s,p]
          }

          appC_spf[s,p,] <- apply( X = appC_axspf[,,s,p,], FUN = sum, MARGIN = c(3))

          for( x in 1:nX )
            for( a in 1:A_s[s] )
              for( f in 1:whichF )
              {
                # Calculate jacobian
                tmpJ1 <- vB_axspf[a,x,s,p,f] / appZ_axsp[a,x,s,p]
                tmpJ2 <- appF_spf[s,p,f] * sel_axspf[a,x,s,p,f] * exp( - appZ_axsp[a,x,s,p] )
                tmpJ3 <- (1 - exp( -appZ_axsp[a,x,s,p]))
                tmpJ4 <- (appZ_axsp[a,x,s,p] - appF_spf[s,p,f] * sel_axspf[a,x,s,p,f])/appZ_axsp[a,x,s,p]

                J_spf[s,p,f] <- J_spf[s,p,f] -  tmpJ1 * ( tmpJ2 + tmpJ3 * tmpJ4)
              }
          
        }

      for( f in whichF )
      {
        f_spf[,,f] <- f_spf[,,f] - appC_spf[,,f]
        appF_spf[,,f] <- appF_spf[,,f] - baranovStep * f_spf[,,f]/J_spf[,,f]
      }

      # Calculate next iteration of Z
      for( s in 1:nS )
        for( p in 1:nP )
        {
          for( x in 1:nX )
          {
            appZ_axsp[1:A_s[s],x,s,p] <- M_xsp[x,s,p]

            for( f in 1:nF )
              appZ_axsp[1:A_s[s],x,s,p] <- appZ_axsp[1:A_s[s],x,s,p] + appF_spf[s,p,f] * sel_axspf[1:A_s[s],x,s,p,f]
          } 
        }

    } # END i
  } # END nIter > 0

  # outputs for baranov solver
  outList <- list(  F_spf   = appF_spf,
                    Z_axsp  = appZ_axsp )
} # END .solveBaranov

# .mfPA()
# Multi-fleet Pope's approximation - 
# I don't think this is suitable for only one commercial fishery...
.mfPA <- function( C_spf     = newCat_spf, 
                  vB_axspf  = vB_axspft[,,,,,t],
                  vB_spf    = vB_spft[,,,t],
                  N_axsp    = N_axspt[,,,,t],
                  sel_axspf = sel_axspft[,,,,,t],
                  M_xsp     = M_xsp,
                  A_s       = A_s,
                  wt_axsp   = meanWtAge_axsp,
                  nS = nS, nP = nP, nX = nX, nF = nF  )
{
  # Make arrays to hold info
  pvB_axspf         <- array(0, dim = dim(vB_axspf))
  remN_axsp         <- array(0, dim = dim(N_axsp))
  endN_axsp         <- array(0, dim = dim(N_axsp))
  appZ_axsp         <- array(0, dim = dim(N_axsp))
  catNumAge_axspf   <- array(0, dim = dim(vB_axspf))
  appF_axspf        <- array(0, dim = dim(vB_axspf) )
  appF_spf          <- array(0, dim = dim(C_spf))

  for( sIdx in 1:nS )
    for( pIdx in 1:nP )
    {
      for( fIdx in 1:nF )
      { 
        # Skip if catch is zero
        if( C_spf[sIdx,pIdx, fIdx] == 0 )
          next

        # Pull biomass at age, convert to proportions
        pvB_axspf[,,sIdx,pIdx,fIdx] <- vB_axspf[,,sIdx,pIdx,fIdx]/vB_spf[sIdx,pIdx,fIdx]

        # Calculate numbers to remove by converting catch to weight
        catNumAge_axspf[,,sIdx,pIdx,fIdx] <- C_spf[sIdx, pIdx, fIdx] * pvB_axspf[,,sIdx,pIdx,fIdx] / wt_axsp[,,sIdx,pIdx]


        # Now add the catch at age in numbers to the removed fish
        remN_axsp[,,sIdx,pIdx] <- remN_axsp[,,sIdx,pIdx] + catNumAge_axspf[,,sIdx,pIdx,fIdx]

      }
      for( xIdx in 1:nX)
        endN_axsp[,xIdx,sIdx,pIdx] <- N_axsp[,xIdx,sIdx,pIdx] * exp(-M_xsp[xIdx,sIdx,pIdx]) - remN_axsp[,xIdx,sIdx,pIdx]*exp(-M_xsp[xIdx,sIdx,pIdx]/2)

      # Now compute Z approximation
      # if( any(is.nan(log( N_axsp[1:A_s[sIdx],,sIdx,pIdx] / endN_axsp[1:A_s[sIdx],,sIdx,pIdx] ))) ) browser()

      appZ_axsp[1:A_s[sIdx],,sIdx,pIdx] <- log( N_axsp[1:A_s[sIdx],,sIdx,pIdx] / endN_axsp[1:A_s[sIdx],,sIdx,pIdx] )
      # Loop over fIdx again
      for( fIdx in 1:nF )
      {
        if( C_spf[sIdx,pIdx,fIdx] == 0 )
          next
        # This is just checking that the F is the same across age classes
        appF_axspf[1:A_s[sIdx],,sIdx,pIdx,fIdx] <- catNumAge_axspf[1:A_s[sIdx],,sIdx,pIdx,fIdx] * appZ_axsp[1:A_s[sIdx],,sIdx,pIdx]
        appF_axspf[1:A_s[sIdx],,sIdx,pIdx,fIdx] <- appF_axspf[1:A_s[sIdx],,sIdx,pIdx,fIdx]/(N_axsp[1:A_s[sIdx],,sIdx,pIdx]*sel_axspf[1:A_s[sIdx],,sIdx,pIdx,fIdx]*(1 - exp(-appZ_axsp[1:A_s[sIdx],,sIdx,pIdx])))
        appF_axspf[is.nan(appF_axspf)] <- 0

        appF_spf[sIdx,pIdx,fIdx] <- appF_axspf[A_s[sIdx],2,sIdx,pIdx,fIdx]
        # if( fIdx == 2 ) browser()
      }

    }

  appF_spf
} # END .mfPA()


