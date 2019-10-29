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
  .saveBlob(blob, ctlTable, outFolder)

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

  # Start with a FALSE array
  idxOn_spft <- array( FALSE, dim = c(nS,nP,nF,nT) )
  # Recover historical indices and set 
  # to TRUE
  histdx <- 1:(tMP-1)
  I_spft <- obj$mp$data$I_spft[,,,histdx]

  # Now, for the future we want to be a bit more selective
  nIdx <- length(idxOn)
  for( i in 1:nIdx )
  {
    fIdx <- idxOn[i]
    for( s in 1:nS )
      for( p in 1:nP )
      {
        histOn <- which(!is.na(I_spft[s,p,fIdx,]))
        if( length(histOn) == 0 )
          next

        lastObs <- max(histOn)
        newObs <- seq( from = lastObs, to = nT, by = obsInt[i] )
        idxOn_spft[s,p,fIdx,c(histOn,newObs)] <- TRUE
      }
  }

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
  obj  <- .applyPrecHCR( obj, t )

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
  hcr$TAC_spt[,,t] <- hcr$targetF_spt[,,t] * projSB_sp

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


  return(obj)
}

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
      obj$mp$assess$retroVB_tspft[pt:pT,s,p,,t] <- tmpObj$om$vB_spft[s,p,,t] 
      
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
  nReps <- ctlList$ctl$nReps
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
  blob <- list( om = NULL,
                mp = NULL )



  # OM lists
  om <- list( iSeed_i   = rep(NA, nReps),
              SB_ispt   = array( NA, dim = c(nReps, nS, nP, nT) ),      # SSB
              B_ispt    = array( NA, dim = c(nReps, nS, nP, nT )),      # Total biomass
              vB_ispft  = array( NA, dim = c(nReps, nS, nP, nF, nT )),  # Vulnerable biomass
              R_ispt    = array( NA, dim = c(nReps, nS, nP, nT) ),      # Rec't
              C_ispft   = array( NA, dim = c(nReps, nS, nP, nF, nT) ),  # Catch by fleet (kt)
              C_ispt    = array( NA, dim = c(nReps, nS, nP, nT) ),      # Total catch (kt)
              F_ispft   = array( NA, dim = c(nReps, nS, nP, nF, nT) ),  # Fishing Mort
              E_ipft    = array( NA, dim = c(nReps, nP, nF, nT) ),      # Fishing effort
              q_ispft   = array( NA, dim = c(nReps, nS, nP, nF, nT) ),  # Catchability
              qF_ispft  = array( NA, dim = c(nReps, nS, nP, nF, nT) ),  # F/E
              Rev_ispft = array( NA, dim = c(nReps, nS, nP, nF, nT ) )  # Revenue 
            )

  om$errors  <- list( omegaR_ispt = array( NA, dim = c(nReps, nS, nP, nT) ),    # Rec Proc Errors
                      delta_ispft = array( NA, dim = c(nReps, nS, nP, nF, nT) ) # Idx obs errors
                    )

  ## ADD AGE/LENGTH SAMPLING ERRORS ## 

  mp <- list( data = NULL, assess = NULL, hcr = NULL )

  mp$data <- list( I_ispft = array( NA, dim = c(nReps, nS, nP, nF, nT ) ) )    # Biomass indices
  ## ADD AGES AND LENGTH DATA ##

  mp$assess <- list(  retroR_itspt    = array( NA, dim = c(nReps, pT, nS, nP, nT ) ),      # retroR
                      retroSB_itspt   = array( NA, dim = c(nReps, pT, nS, nP, nT ) ),      # retroSB
                      retroVB_itspft  = array( NA, dim = c(nReps, pT, nS, nP, nF, nT ) ),  # retroVB
                      assessMod_i     = vector( mode = "list", length = nReps ) )          # assessment vital stats


  mp$hcr <- list( Fref_ispt = array( NA, dim = c(nReps, nS, nP, nT ) ),       # Reference F in HCR
                  Bref_ispt = array( NA, dim = c(nReps, nS, nP, nT ) ),       # Reference biomass (B0 or Bmsy, usually)
                  targetF_ispt = array( NA, dim = c(nReps, nS, nP, nT) ),     # Target F from HCR
                  LCP_ispt = array(NA, dim = c(nReps, nS, nP, nT) ),          # Lower control point
                  UCP_ispt = array(NA, dim = c(nReps, nS, nP, nT) ),          # Upper control point
                  TAC_ispt = array(NA, dim = c(nReps, nS, nP, nT ) ),         # TAC summed across all fleets
                  TAC_ispft = array(NA, dim = c(nReps, nS, nP, nF, nT ) ) )   # TAC allocated by fleet


  blob <- list( om = om, mp = mp, ctlList = ctlList )


  ##########################################################
  ######### ------- CLOSED LOOP SIMULATION ------- #########
  ##########################################################

  if(!obj$ctlList$ctl$quiet)
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
    simObj <- .calcTimes( simObj )

    if( nT >= tMP )
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
      } # END t loop
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
    blob$om$E_ipft[i,,,]      <- simObj$om$E_pft
    blob$om$q_ispft[i,,,,]    <- simObj$om$q_spft
    blob$om$qF_ispft[i,,,,]   <- simObj$om$qF_spft
    blob$om$Rev_ispft[i,,,,]  <- simObj$om$Rev_spft

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
    blob$mp$hcr$TAC_ispft[i,,,,]      <- simObj$mp$hcr$TAC_spft


    # Retrospective AM outputs
    blob$mp$assess$retroR_itspt[i,,,,]    <- simObj$mp$assess$retroR_tspt
    blob$mp$assess$retroSB_itspt[i,,,,]   <- simObj$mp$assess$retroSB_tspt
    blob$mp$assess$retroVB_itspft[i,,,,,] <- simObj$mp$assess$retroVB_tspft
    blob$mp$assess$assessMod_i[[i]]       <- simObj$mp$assess$assessModOutputs
  


    message( " (.mgmtProc) Completed replicate ", i, " of ", nReps, ".\n", sep = "")
  } # END i loop


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
  whichFleets <- which( ctlList$opMod$alloc > 0 )

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

  # Now, copy model fit
  histdx <- 1:(tMP - 1)
  obj$om$F_spft[,,,histdx]        <- repObj$F_spft
  obj$om$C_spft[,,,histdx]        <- repObj$predCw_spft
  obj$om$C_spt[,,histdx]          <- apply( X = repObj$predCw_spft, FUN = sum, MARGIN = c(1,2,4) )
  obj$om$E_pft[,1:2,histdx]       <- repObj$E_pft[,1:2,]
  obj$om$sel_axspft[,,,,,histdx]  <- repObj$sel_axspft
  obj$om$sel_lspft[,,,,histdx]    <- aperm(repObj$sel_lfspt,c(1,3,4,2,5))

  # Now fill in future sel
  obj$om$sel_axspft[,,,,,tMP:nT]  <- repObj$sel_axspft[,,,,,tMP-1]
  obj$om$sel_lspft[,,,,tMP:nT]    <- obj$om$sel_lspft[,,,,tMP-1]

  # Copy index catchability
  obj$om$q_spft[,,,histdx]        <- repObj$q_spft
  obj$om$q_spf                    <- repObj$q_spf

  # Now commercial fishing mortality catchability
  for( s in 1:nS )
    obj$om$qF_spft[s,,,histdx]      <- repObj$F_spft[s,,,] / repObj$E_pft

  obj$om$qF_spft[!is.finite(obj$om$qF_spft)] <- 0

  # Fill in future q as fixed at conditioning year's value
  obj$om$q_spft[,,,tMP:nT]        <- repObj$q_spft[,,,tMP-1]
  obj$om$qF_spft[,,,tMP:nT]       <- obj$om$qF_spft[,,,tMP-1]


  # Add historical data
  obj$mp$data$I_spft[,,,histdx]      <- repObj$I_spft
  obj$mp$data$A_axspft[,,,,,histdx]  <- aperm( repObj$age_aspftx, c(1,6,2:5) )
  obj$mp$data$L_lxspft[,,,,,histdx]  <- aperm( repObj$len_lspftx, c(1,6,2:5) )

  # replace negatives with NAs for plotting, can change back 
  # later for TMB
  obj$mp$data$I_spft[obj$mp$data$I_spft<0] <- NA
  obj$mp$data$A_axspft[obj$mp$data$A_axspft<0] <- NA
  obj$mp$data$L_lxspft[obj$mp$data$A_axspft<0] <- NA

  # Now calculate total catch and allocation
  tdxRecent <- (tMP - ctlList$opMod$allocYears):(tMP - 1)
  recentCatch_spf <- apply( X = obj$om$C_spft[,,,tdxRecent], FUN = sum, MARGIN = c(1,2,3) )

  commGears <- ctlList$opMod$commGears

  # Loop and fill errors and allocation
  for( s in 1:nS )
    for( p in 1:nP )
    {
      obj$errors$omegaR_spt[s,p,] <- rnorm(nT)
      for( f in 1:nF )
        obj$errors$delta_spft[s,p,f,] <- rnorm(nT)

      obj$om$alloc_spf[s,p,commGears] <- recentCatch_spf[s,p,commGears] / sum( recentCatch_spf[s,p,commGears])

      # Save historical proc errors, but use simulated recruitments after
      # last estimated recruitment
      lastIdx <- max(which(repObj$omegaR_spt[s,p,] != 0) )
      obj$errors$omegaR_spt[s,p,histdx[1:lastIdx]]   <- repObj$omegaR_spt[s,p,1:lastIdx]  # rec devs    
    }


  # Save historical errors
  obj$errors$delta_spft[,,,histdx]  <- repObj$residCPUE_spft # obs errors
  obj$errors$omegaRinit_asp         <- repObj$omegaRinit_asp # Initialisation errors

  message(" (.condMS3pop) Running OM for historical period.\n")

  # Now, initialise the population
  for( t in 1:(tMP - 1) )
  {
    obj <- .ageSexOpMod( obj, t )
  }

  obj$om$speciesNames   <- dimnames(ctlList$opMod$histRpt$SB_spt)[[1]]
  obj$om$stockNames     <- dimnames(ctlList$opMod$histRpt$SB_spt)[[2]]
  obj$om$fleetNames     <- dimnames(ctlList$opMod$histRpt$vB_spft)[[3]]

  obj$om$price_s        <- ctlList$opMod$price_s[obj$om$speciesNames]
  obj$om$alphaU         <- ctlList$opMod$alphaU
  obj$om$ut_50          <- ctlList$opMod$ut_50
  obj$om$ut_95          <- ctlList$opMod$ut_95

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
  om$B_axspt <- array(0,  dim = c(nA,nX,nS,nP,nT) )  # Biomass at age (total)
  om$N_axspt <- array(0,  dim = c(nA,nX,nS,nP,nT) )  # Numbers at age
  om$SB_spt  <- array(0,  dim = c(nS,nP,nT) )        # Spawning biomass
  om$B_spt   <- array(0,  dim = c(nS,nP,nT) )        # Spawning biomass
  om$R_spt   <- array(0,  dim = c(nS,nP,nT) )        # Recruitment
  om$Z_axspt <- array(0,  dim = c(nA,nX,nS,nP,nT) )  # Total mortality

  # Catch, fishing mort, vuln bio, selectivity
  om$C_spft     <- array(NA,  dim = c(nS,nP,nF,nT) )       # Total catch by fleet in kt
  om$C_spt      <- array(NA,  dim = c(nS,nP,nT) )          # Total catch in kt
  om$C_axspft   <- array(0,  dim = c(nA,nX,nS,nP,nF,nT) )  # Numbers
  om$Cw_axspft  <- array(0,  dim = c(nA,nX,nS,nP,nF,nT) )  # Weight
  om$vB_axspft  <- array(0,  dim = c(nA,nX,nS,nP,nF,nT) )  # vuln Bio (granular)
  om$vB_spft    <- array(0,  dim = c(nS,nP,nF,nT) )        # vuln Bio (aggregate)
  om$sel_axspft <- array(0,  dim = c(nA,nX,nS,nP,nF,nT) )  # sel at age (tv)
  om$sel_lspft  <- array(0,  dim = c(nL,nS,nP,nF,nT) )     # sel at len (tv)
  om$sel_axspf  <- array(0,  dim = c(nA,nX,nS,nP,nF) )     # sel at age
  om$sel_lspf   <- array(0,  dim = c(nL,nS,nP,nF) )        # sel at len
  om$F_spft     <- array(NA,  dim = c(nS,nP,nF,nT) )       # Fishing mortality
  om$E_pft      <- array(NA,  dim = c(nP,nF,nT) )          # Fishing effort
  om$alloc_spf  <- array(0,  dim = c(nS,nP,nF))            # Catch allocation
  om$Rev_spft   <- array(NA, dim = c(nS,nP,nF,nT) )        # Fleet revenue

  # Observation model pars
  om$q_spft     <- array(0,  dim = c(nS,nP,nF,nT) )        # catchability (tv)
  om$q_spf      <- array(0,  dim = c(nS,nP,nF) )           # catchability

  # Scalar between effort and F
  om$qF_spft    <- array(0,  dim = c(nS,nP,nF,nT) )        # catchability (tv)

  # Variance parameters
  om$tauObs_spf <- repObj$tauObs_spf
  om$sigmaR_sp  <- repObj$sigmaR_sp

  # Leading bio pars
  om$B0_sp      <- repObj$B0_sp
  om$M_xsp      <- repObj$M_xsp
  om$h_sp       <- repObj$h_sp

  # Growth model pars
  om$L1_s       <- repObj$L1_s
  om$A1_s       <- repObj$A1_s
  om$L2_xsp     <- aperm(repObj$L2_spx,c(3,1,2))
  om$A2_s       <- repObj$A2_s
  om$vonK_xsp   <- aperm(repObj$vonK_spx,c(3,1,2))
  om$LWa_s      <- repObj$LWa_s
  om$LWb_s      <- repObj$LWb_s

  # Later: write calcSchedules function
  om$matAge_asp         <- repObj$matAge_asp
  om$meanWtAge_axsp     <- aperm( repObj$meanWtAge_aspx,c(1,4,2,3) )
  om$Wlen_ls            <- repObj$Wlen_ls
  om$probLenAge_laxsp   <- aperm( repObj$probLenAge_laspx, c(1,2,5,3,4) )
  om$lenAge_axsp        <- aperm( repObj$lenAge_aspx, c(1,4,2,3) )

  # Calculate ref points and get R0 and recruitment pars
  refPtList             <- calcRefPts( repObj )$refPts

  # Make data list - anything else?
  # mp list - mostly for retrospective
  # analyses and assessment model settings
  mp <- list( assess = NULL, data = NULL, hcr = NULL )
  mp$assess <- list()
  mp$hcr    <- list()
  mp$data <- list()

  # Data list
  mp$data$I_spft    <- array( NA, dim = c(nS,nP,nF,nT) )
  mp$data$A_axspft  <- array( NA, dim = c(nA,nX,nS,nP,nF,nT) ) 
  mp$data$L_lxspft  <- array( NA, dim = c(nL,nX+1,nS,nP,nF,nT) ) 

  # assessment list
  mp$assess$retroSB_tspt      <- array( NA, dim = c(pT, nS, nP, nT) )
  mp$assess$retroR_tspt       <- array( NA, dim = c(pT, nS, nP, nT) )
  mp$assess$retroVB_tspft      <- array( NA, dim = c(pT, nS, nP, nF, nT) )
  mp$assess$assessModOutputs  <- data.frame( matrix(NA, nrow = pT, ncol = 25) )

  mp$hcr$LCP_spt              <- array( NA, dim = c(nS, nP, nT) )
  mp$hcr$UCP_spt              <- array( NA, dim = c(nS, nP, nT) )
  mp$hcr$Bref_spt             <- array( NA, dim = c(nS, nP, nT) )
  mp$hcr$Fref_spt             <- array( NA, dim = c(nS, nP, nT) )
  mp$hcr$targetF_spt          <- array( NA, dim = c(nS, nP, nT) )
  mp$hcr$TAC_spft             <- array(0,  dim = c(nS,nP,nF,nT) )     # MP TAC by fleet in kt
  mp$hcr$TAC_spt              <- array(0,  dim = c(nS,nP,nT) )        # MP TAC in kt


  # Process and observation errors
  errors <- list()

  # Arrays for holding errors
  errors$omegaR_spt <- array(NA, dim = c(nS,nP,nT) )
  errors$delta_spft <- array(NA, dim = c(nS,nP,nF,nT) )




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
                rp = refPtList, 
                errors = errors,
                ctlList = ctlList )

  return(obj)
} # END initMS3pop()

# Operating model function (similar to mseR)
.ageSexOpMod <- function( obj, t )
{
  # Pull om
  om    <- obj$om
  rp    <- obj$rp
  data  <- obj$mp$data
  err   <- obj$errors
  hcr   <- obj$mp$hcr

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
  Rev_spft          <- obj$om$Rev_spft 

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
  I_spft            <- data$I_spft
  
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
            R_spt[s,p,t] <- R_spt[s,p,t] * exp( om$sigmaR_sp[s,p] * err$omegaR_spt[s,p,t]) 

            N_axspt[a,,s,p,t] <- R_spt[s,p,t] / nX
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

  # Now calculate effort for each area (for closed loop sim)
  if( t >= tMP )
  {
    effortMod <- .solveProjEffortDynamics( TAC_spf   = TAC_spft[,,,t], 
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
    Rev_spft[,,,t] <- effortMod$rev_spf

    # Convert to F
    for( s in 1:nS )
      F_spft[s,,,t] <- E_pft[,,t] * qF_spft[s,,,t]

  } # END F approximation

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

  # Now generate indices - need to add a schedule later
  if( t >= tMP )
  {
    # Get the indices that we're still collecting
    idxOn_spft <- obj$ctlList$mp$data$idxOn_spft

    # Calculate and save
    for( s in 1:nS )
      for( p in 1:nP )
      {
        idxOn <- obj$mp$data$idxOn_spft[s,p,,t]
        I_spft[s,p,idxOn,t] <- q_spft[s,p,idxOn,t] * vB_spft[s,p,idxOn,t] * exp(om$tauObs_spf[s,p,idxOn] * err$delta_spft[s,p,idxOn,t])
      }
  }


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
  Rev_spft          -> obj$om$Rev_spft
  
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
  I_spft            -> obj$mp$data$I_spft

  return( obj )
} # END ageSexOpMod()

# mfPA()
# Multi-fleet Pope's approximation
mfPA <- function( C_spf     = newCat_spf, 
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
} # END mfPA()


