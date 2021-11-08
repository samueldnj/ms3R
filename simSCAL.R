# <><><><><><><><><><><><><><><><><><><><><><><><><><><>><><><><>><><><><>
# simSCAL.R
#
# Simulation model for generating data for testing a statistical catch 
# at age and length (SCAL) model, or providing an operating model for 
# closed loop simulation.
#
# Author: SDN Johnson
# Date: July 25, 2019
#
# Last Update: SDN Johnson, May 13 2021
#
# TO DO:
# - replace NAs with zeros for histroical Index in Lou
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><>><><><><>><><><><>

runMS3 <- function( ctlFile = "./simCtlFile.txt",
                    outFolder = NULL )
{
  # Load control list
  ctlTable <- .readParFile ( ctlFile )
  # Create control list
  ctlList  <- .createList  ( ctlTable )

  # Load AMs if necessary
  cat("Loading AM for simulation")
  AMid <- ctlList$mp$assess$method
  if(AMid %in% c("hierProd","SISCA"))
    source(paste0("AM_",AMid,".R"))

  # Load history
  repList                   <- .loadFit( ctlList$opMod$histFile )
  repList$repOpt$apF_pgt    <- approxCtsF(repList)$apF_pgt
  repList$repOpt$fYear      <- repList$fYear
  ctlList$opMod$histRpt     <- c( repList$data, repList$repOpt, repList$pars ) 
  ctlList$opMod$posts       <- repList$posts
  ctlList$opMod$fYear       <- repList$fYear
  ctlList$opMod$species     <- repList$species
  ctlList$opMod$stocks      <- repList$stocks
  ctlList$opMod$fleets      <- repList$gearLabs
  ctlList$opMod$histRpt$map <- repList$map
  ctlList$opMod$histRpt$phz <- repList$ctlList$phases

  # Initialise the data structures to hold simulation states etc.
  simObj <- .initMS3pop(  ctlList )

  # Run mgmtProc
  blob <- .mgmtProc( simObj )

  # Save output
  .saveBlob(blob = blob, ctlTable, outFolder)

  # beepr::beep("complete")
}


# calcObsTimes()
# Creates a schedule for survey indices and
# assessments based on control file settings
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
  nX <- obj$om$nX
  nF <- obj$om$nF


  # Get intervals for assessments
  assInt    <- obj$ctlList$mp$assess$assInt
  assTimes  <- seq(from = tMP, by = assInt, length.out = pT)

  assessOn <- c(rep(NA, tMP-1),rep(FALSE,pT))
  assessOn[ assTimes ]<- TRUE

  # Now do the same for indices
  idxOn   <- obj$ctlList$mp$data$idxOn
  obsInt  <- obj$ctlList$mp$data$obsInt_i

  ageOn   <- obj$ctlList$mp$data$ageOn
  ageInt  <- obj$ctlList$mp$data$ageInt

  # Start with a FALSE array, add pooled data dimension
  idxOn_spft <- array( FALSE, dim = c(nS+1,nP+1,nF,nT) )
  ageOn_spft <- array( FALSE, dim = c(nS+1,nP+1,nF,nT) )
  # Recover historical indices and set 
  # to TRUE
  histdx <- 1:(tMP-1)
  I_spft    <- obj$mp$data$I_spft[,,,histdx,drop = FALSE]
  A_axspft  <- obj$mp$data$A_axspft[,,,,,histdx,drop = FALSE]

  # Local helper function
  checkAnyNonNA <- function( data )
  {
    any(!is.na(data))
  }

  for( s in 1:nS )
    for( p in 1:nP )
      for( fIdx in 1:nF )
      {
        # first, do indices
        histOn <- which(!is.na(I_spft[s,p,fIdx,]))
        if( length(histOn) == 0 )
          idxOn_spft[s,p,fIdx,] <- FALSE

        idxOn_spft[s,p,fIdx,histOn] <- TRUE

        # Now generate new data for the future
        if( fIdx %in% idxOn )
        {
          i <- which( idxOn == fIdx )

          if(length(histOn) > 0 )
            lastObs <- max(histOn)
          else
            lastObs <- tMP

          if(obsInt[i] == 1)
            lastObs <- tMP

          newObs <- seq( from = lastObs, to = nT, by = obsInt[i] )
          idxOn_spft[s,p,fIdx,newObs] <- TRUE

        }

        # Then ages

        histOn <- apply(X = A_axspft[2,1:nX,s,p,fIdx,,drop = FALSE], FUN = checkAnyNonNA, MARGIN = 6)
        if(length(histOn) == 0)
          ageOn_spft[s,p,fIdx,] <- FALSE

        if( fIdx %in% ageOn )
        {
          i <- which( ageOn == fIdx)

          if(length(histOn) > 0)
            lastObs <- max(histOn)
          else 
            lastOns <- tMP

          if( ageInt[i] == 1)
            lastObs <- tMP

          newObs <- seq( from  =lastObs, to = nT, by = ageInt[i])
          ageOn_spft[s,p,fIdx,newObs] <- TRUE

        }
      }

  # Now fill the pooled data dimensions
  idxOn_spft[s+1,1:nP,,]  <- as.logical( apply( X = idxOn_spft[1:nS,1:nP,,,drop = FALSE], FUN = sum, MARGIN = c(2,3,4) ) )
  idxOn_spft[1:nS,p+1,,]  <- as.logical( apply( X = idxOn_spft[1:nS,1:nP,,,drop = FALSE], FUN = sum, MARGIN = c(1,3,4) ) )
  idxOn_spft[s+1,p+1,,]   <- as.logical( apply( X = idxOn_spft[1:nS,1:nP,,,drop = FALSE], FUN = sum, MARGIN = c(3,4)   ) )

  ageOn_spft[s+1,1:nP,,]  <- as.logical( apply( X = ageOn_spft[1:nS,1:nP,,,drop = FALSE], FUN = sum, MARGIN = c(2,3,4) ) )
  ageOn_spft[1:nS,p+1,,]  <- as.logical( apply( X = ageOn_spft[1:nS,1:nP,,,drop = FALSE], FUN = sum, MARGIN = c(1,3,4) ) )
  ageOn_spft[s+1,p+1,,]   <- as.logical( apply( X = ageOn_spft[1:nS,1:nP,,,drop = FALSE], FUN = sum, MARGIN = c(3,4)   ) )

  obj$mp$assess$assessOn_t  <- assessOn
  obj$mp$data$idxOn_spft    <- idxOn_spft
  obj$mp$data$ageOn_spft    <- ageOn_spft

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
  # 4. Allocate catch among fleets
  # 5. update operating model state variables and data
  # 6. return updated simulation object

  # 1. Extract operating model and management procedure objects
  ctlList <- obj$ctlList
  
  rp <- obj$rp      # reference points
  om <- obj$om      # operating model
  mp <- obj$mp      # management procedure

  nT      <- om$nT
  tMP     <- om$tMP
  nS      <- om$nS
  nP      <- om$nP
  nF      <- om$nF
  nMICE   <- om$nMICE
  nSpec   <- nS - nMICE

  commGears <- ctlList$opMod$commGears

  # 2. Perform stock assessment
  obj  <- .callProcedureAM( obj, t )

  # 3. Determine stock status and apply HCR to set catch limit
  if(ctlList$ctl$mpName=='NoFish')
    obj$mp$hcr$TAC_spft[1:nSpec,,commGears,t] <- 0
  else 
  {  
    if( ctlList$mp$hcr$type == "ramped")
      obj  <- .applyPrecHCR( obj, t )

    if( ctlList$mp$hcr$type == "conF" )
      obj  <- .applyConstantF( obj, t )

    if( ctlList$mp$hcr$type == "ccRule" )
      obj  <- .calcHCR_ccRule( obj, t )

    if( ctlList$mp$hcr$type == "hgRule" )
      obj  <- .calcHCR_hgRule( obj, t )

    # 4. Allocate catch among fleets
    if( ctlList$mp$hcr$type == "hgRule" )
      obj <- .tacAlloc_hgRule(obj,t)

    # Need to unift the allocation for SOK/HG and SOG. Right now
    if( !ctlList$mp$hcr$type == "hgRule" )
      for( s in 1:nS )
        for( p in 1:nP )
        {
          obj$mp$hcr$TAC_spft[s,p,,t]    <- obj$mp$hcr$TAC_spt[s,p,t] * obj$om$alloc_spf[s,p,]
        }

    # SOK product per license is 16,000 lbs (7.2575 t or 0.0072575 kt) per pond ()
    sokPerLic_kt <- 0.0072575
    # obj <- .tac2SOK(obj, t)
    closedPondAlloc <- 0.0907 # 100 short tons (90.7 t or 0.0907 kt) of TAC allocated to each closed pond,
    openPondAlloc   <- 0.0317 # 35 short tons (31.7t or 0.0317 kt) allocated to each open pond

    histRpt <- ctlList$opMod$histRpt
    sokIdx <- which(histRpt$fleetType_g == 2)

     # TOO PATH SPECIFIC
    if( length(sokIdx) > 0 )
      for( s in 1:nS )
        for( p in 1:nP ) 
        {
          obj$mp$hcr$TAC_spft[s,p,sokIdx,t] <- obj$mp$hcr$TAC_spft[s,p,sokIdx,t]*sokPerLic_kt/closedPondAlloc

          if(nF > histRpt$nG)
            obj$mp$hcr$TAC_spft[s,p,7,t] <- obj$mp$hcr$TAC_spft[s,p,7,t]*sokPerLic_kt/openPondAlloc
        }
  }  

  # 5. Update simulated population
  obj <- .ageSexOpMod(obj, t)

  # 6. Return updated sim object
  return(obj)

} # END .updatePop()


# .mgmtProc()
# Runs replications of the closed loop
# simulation, applying the simulated mgmt
# procedure to the operating model stocks
.mgmtProc <- function( obj )
{

  ctlList <- obj$ctlList
  opMod   <- ctlList$opMod

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


  #### ----- THE BLOOOBBBB ----- ####


  # OM lists
  om <- list( iSeed_i     = rep(NA, nReps),
              SB_ispt     = array( NA, dim = c(nReps, nS, nP, nT) ),          # SSB
              endSB_ispt  = array( NA, dim = c(nReps, nS, nP, nT) ),          # SSB + surviving ponded fish
              B_ispt      = array( NA, dim = c(nReps, nS, nP, nT )),          # Total biomass
              vB_ispft    = array( NA, dim = c(nReps, nS, nP, nF, nT )),      # Vulnerable biomass
              R_ispt      = array( NA, dim = c(nReps, nS, nP, nT) ),          # Rec't
              C_ispft     = array( NA, dim = c(nReps, nS, nP, nF, nT) ),      # Catch by fleet (kt)
              C_ispt      = array( NA, dim = c(nReps, nS, nP, nT) ),          # Total catch (kt)
              F_ispft     = array( NA, dim = c(nReps, nS, nP, nF, nT) ),      # Fishing Mort
              P_ispft     = array( NA, dim = c(nReps, nS, nP, nF, nT ) ),     # Ponded fish in biomass units
              I_ispft     = array( NA, dim = c(nReps, nS+1, nP+1, nF, nT) ),  # Indices (without error)
              rI_ispft    = array( NA, dim = c(nReps, nS+1, nP+1, nF, nT) ),  # proportion of index from each gear
              E_ipft      = array( NA, dim = c(nReps, nP, nF, nT) ),          # Fishing effort
              q_ispft     = array( NA, dim = c(nReps, nS, nP, nF, nT) ),      # Catchability
              qF_ispft    = array( NA, dim = c(nReps, nS, nP, nF, nT) ),      # F/E
              Rev_ispft   = array( NA, dim = c(nReps, nS, nP, nF, nT ) ),     # Revenue 
              psi_ispft   = array( NA, dim = c(nReps, nS, nP, nF, nT) ),      # SOK conversion rate 
              pondM_ift   = array( NA, dim = c(nReps, nF, nT) ),              # post ponding mortality 
              M_iaxspt    = array( NA, dim = c(nReps, nA, nX, nS, nP, nT)),    # Natural mortality             
              N_iaxspt    = array( NA, dim = c(nReps, nA, nX, nS, nP, nT))    # Numbers-at-age

            )

  om$errors  <- list( omegaR_ispt = array( NA, dim = c(nReps, nS, nP, nT) ),    # Rec Proc Errors
                      delta_ispft = array( NA, dim = c(nReps, nS+1, nP+1, nF, nT) ) # Idx obs errors
                    )

  ## ADD AGE/LENGTH SAMPLING ERRORS ## 

  mp <- list( data = NULL, assess = NULL, hcr = NULL )

  mp$data <- list(  I_ispft   = array( NA, dim = c(nReps, nS+1, nP+1, nF, nT ) ),          # Biomass/abundance indices
                    A_iaxspft = array( NA, dim = c(nReps, nA, nX, nS+1, nP+1, nF, nT ) ),  # age comp data
                    L_ilxspft = array( NA, dim = c(nReps, nL, nX, nS+1, nP+1, nF, nT ) )   # len comp data
                  )
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
                  TAC_ispft     = array(NA, dim = c(nReps, nS, nP, nF, nT ) ),   # TAC allocated by fleet
                  sokEff_ispft  = array(NA, dim = c(nReps, nS, nP, nF, nT ) ),   # number of SOK licenses
                  maxTAC_ispft  = array(NA, dim = c(nReps, nS, nP, nF, nT ) )    # maximum TAC before SOK or SR caps are applied
                  )  


  blob <- list( om = om, mp = mp, ctlList = ctlList,
                rp = vector(mode = "list", length = nReps) )

  blob$goodReps       <- rep(FALSE, nReps )

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

  # ---- Posterior Sampling ----

  message(" (.mgmtProc) Posterior sampling for leading parameters...\n")
  if(ctlList$opMod$posteriorSamples)
  {
    
    
    # Get posterior samples
    mcmcPar <- ctlList$opMod$posts
    nMCMC   <- dim(mcmcPar$B0_ip)[1]
    nReps   <- ctlList$ctl$totReps

    if( ctlList$opMod$postSampleType == "quantile" )
    {  
      
      browser(cat('not yet implemented....'))
      samples       <- .quantileStratSample(  seed  = ctlList$opMod$postSampleSeed,
                                              post = mcmcPar,
                                              par1 = "M_iapt", par2 = "B0_ip",
                                              nBreaks = 10 )

    }  

    # if( ctlList$opMod$postSampleType == "full" )
    #   samples <- 1:nrow(mcmcPar)



    if( ctlList$opMod$postSampleType == "random" )
    {
      
      if( nReps > nMCMC )
        setReplace <- TRUE
      else setReplace <- FALSE

      set.seed(ctlList$opMod$postSampleSeed)

      samples <- sample(1:nMCMC, size = nReps, replace = setReplace )
    }

    ctlList$opMod$postDraws_i <- samples

  }



  ##########################################################
  ######### ------- CLOSED LOOP SIMULATION ------- #########
  ##########################################################


  
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

    # replace leading pars
    if(ctlList$opMod$posteriorSamples)
    {
      postDrawIdx <- ctlList$opMod$postDraws_i[i]

      # Replace MLE values for B0 and M with posterior draws
      simObj$om$B0_sp[1:nS,]             <- mcmcPar$B0_ip[postDrawIdx,]
      simObj$om$Rbar_sp[1:nS,]           <- mcmcPar$Rbar_ip[postDrawIdx,]

      # Fix recruitment at values from mcmc draw
      simObj$om$R_spt[1:nS,,1:(tMP-1)]   <- mcmcPar$R_ipt[postDrawIdx,,1:(tMP-1)]

      simObj$om$Rinit_sp[1:nS,1:nP]      <- mcmcPar$Rinit_ip[postDrawIdx,]
      simObj$om$M_axspt[,1,1,,1:(tMP-1)] <- mcmcPar$M_iapt[postDrawIdx,,,1:(tMP-1)]

      # other mcmc pars are set in condMSEpop_SISCA
      # - q, rec devs, selectivity

      simObj$ctlList$opMod$postDrawIdx <- postDrawIdx

    }  

    # Either simulate or condition the model history.

    # If simulated, we cab optimise the removals for a target
    # initial depletion - might need to optimise abs catches for
    # discrete fisheries.
    if(opMod$histType == "sim")
    {
      simObj <- .solveHistPop( simObj )
    }
    # If conditioned, check the 
    if( opMod$histType == "cond")
    {
      # Condition history of simulation model
      # and draw random errors for projection
      if( opMod$condModel == "hierSCAL")
        simObj <- .condMS3pop( simObj )
      if( opMod$condModel == "SISCA" )
        simObj <- .condMS3pop_SISCA( simObj )
    }




    # Calculate effort dynamics parameters
    if( opMod$condModel == "hierSCAL" & opMod$effortMod == "dynModel" )
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
    blob$om$SB_ispt[i,,,]       <- simObj$om$SB_spt
    blob$om$endSB_ispt[i,,,]    <- simObj$om$endSB_spt
    blob$om$B_ispt[i,,,]        <- simObj$om$B_spt
    blob$om$vB_ispft[i,,,,]     <- simObj$om$vB_spft
    blob$om$R_ispt[i,,,]        <- simObj$om$R_spt
    blob$om$C_ispft[i,,,,]      <- simObj$om$C_spft
    blob$om$C_ispt[i,,,]        <- simObj$om$C_spt
    blob$om$F_ispft[i,,,,]      <- simObj$om$F_spft
    blob$om$I_ispft[i,,,,]      <- simObj$om$I_spft
    blob$om$E_ipft[i,,,]        <- simObj$om$E_pft
    blob$om$q_ispft[i,,,,]      <- simObj$om$q_spft
    blob$om$qF_ispft[i,,,,]     <- simObj$om$qF_spft
    blob$om$Rev_ispft[i,,,,]    <- simObj$om$Rev_spft
    blob$om$P_ispft[i,,,,]      <- simObj$om$P_spft
    blob$om$psi_ispft[i,,,,]    <- simObj$om$psi_spft
    blob$om$pondM_ift[i,,]      <- simObj$om$pondM_ft
    blob$om$M_iaxspt[i,,,,,]    <- simObj$om$M_axspt
    blob$om$N_iaxspt[i,,,,,]    <- simObj$om$N_axspt

    # Errors - maybe update simObj structure to match blob here
    blob$om$errors$omegaR_ispt[i,,,]  <- simObj$errors$omegaR_spt
    blob$om$errors$delta_ispft[i,,,,] <- simObj$errors$delta_spft

    # Data
    blob$mp$data$I_ispft[i,,,,]         <- simObj$mp$data$I_spft
    blob$mp$data$A_iaxspft[i,,,,,,]    <- simObj$mp$data$A_axspft
    

    # Add lengths later

    # HCR elements
    blob$mp$hcr$Fref_ispt[i,,,]       <- simObj$mp$hcr$Fref_spt
    blob$mp$hcr$Bref_ispt[i,,,]       <- simObj$mp$hcr$Bref_spt
    blob$mp$hcr$targetF_ispt[i,,,]    <- simObj$mp$hcr$targetF_spt
    blob$mp$hcr$LCP_ispt[i,,,]        <- simObj$mp$hcr$LCP_spt
    blob$mp$hcr$UCP_ispt[i,,,]        <- simObj$mp$hcr$UCP_spt
    blob$mp$hcr$TAC_ispt[i,,,]        <- simObj$mp$hcr$TAC_spt
    blob$mp$hcr$propTAC_ispt[i,,,]    <- simObj$mp$hcr$propTAC_spt
    blob$mp$hcr$TAC_ispft[i,,,,]      <- simObj$mp$hcr$TAC_spft
    blob$mp$hcr$sokEff_ispft[i,,,,]   <- simObj$mp$hcr$sokEff_spft
    # blob$mp$hcr$maxTAC_ispft[i,,,,]   <- simObj$mp$hcr$maxTAC_spft

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


    if( prod(simObj$mp$assess$pdHess_tsp) == 1 | 
        prod(simObj$mp$assess$posSDs_tsp) == 1 | 
        ctlList$mp$assess$method %in% c("ItCt","idxBased","perfectInfo") )
      blob$goodReps[i] <- TRUE

    if( ctlList$ctl$omni | ctlList$ctl$perfConF )
    {
      blob$goodReps[i] <- TRUE

      if(ctlList$ctl$omni)
      {
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

    }

    # Save reference points for this replicate - more necessary when
    # we have multiple conditioning assessment posterior dist draws
    blob$rp[[i]] <- simObj$rp

    message( " (.mgmtProc) Completed replicate ", i, " of ", nReps, ".\n", sep = "")

    if( sum(blob$goodReps) == ctlList$ctl$nGoodReps )
    {
      message(  " (.mgmtProc) Completed ",  ctlList$ctl$nGoodReps, 
                " replicates with all convergent AMs, ending simulation.\n", sep = "")
      break
    } else {
      finishedGood <- sum(blob$goodReps)
      message( " (.mgmtProc) Completed ", finishedGood, " of ", ctlList$ctl$nGoodReps, 
                " replicates with all convergent AMs.")
    }
    # Clear memory of extraneous garbage
    gc()
  } # END i loop

  blob$nSims <- i

  # save initial model times for plotting
  blob$om$tInit_sp <- simObj$om$tInit_sp

  # Name the blob array dimensions...
  tMP   -> blob$om$tMP
  nT    -> blob$om$nT
  nX    -> blob$om$nX
  nF    -> blob$om$nF
  nA    -> blob$om$nA
  nS    -> blob$om$nS
  nP    -> blob$om$nP
  nL    -> blob$om$nL

  blob$om$fleetType_f <- simObj$om$fleetType_f

  blob$om$speciesNames  <- dimnames(ctlList$opMod$histRpt$SB_spt)[[1]]
  blob$om$stockNames    <- dimnames(ctlList$opMod$histRpt$SB_spt)[[2]]
  blob$om$fleetNames    <- dimnames(ctlList$opMod$histRpt$vB_spft)[[3]]




  return(blob)
} # END .mgmtProc()


# .initMS3pop()
# Initialises operating model data structures
# for holding model states, observation and
# process errors, and simulated MP data.
.initMS3pop <- function( ctlList )
{
  repObj <- ctlList$opMod$histRpt
  pT     <- ctlList$opMod$pT
  opMod  <- ctlList$opMod

  # make OM list
  om <- list()

  message(" (.initMS3pop) Initialising MS3 operating model.\n")

  # Get model dimensions
  if(opMod$histType == "sim")
    tMP <- opMod$histT + 1
  if(opMod$histType == "cond")
    tMP   <- repObj$nT + 1   # start of projections
  
  nT    <- tMP + pT - 1    # Total time period of sims

  nA    <- repObj$nA       # Number of Age classs
  nP    <- repObj$nP       # Number of stocks

  if(ctlList$opMod$condModel == "hierSCAL")
  {
    nX    <- repObj$nX       # Number of sex classes
    nF    <- repObj$nF       # Number of fleets
    nS    <- repObj$nS       # Number of species
    nL    <- repObj$nL       # Number of length classes

    # Get species specific largest age/length classes
    om$A_s <- repObj$A_s
    om$L_s <- repObj$L_s

    # Get length bin aggregation info - might refine later.
    om$lenBinWidth   <- repObj$lenBinWidth   # width of length bins
    om$lenBinMids_l  <- repObj$lenBinMids_l  # Mid points
  }

  if( ctlList$opMod$condModel == "SISCA" )
  {
    nP    <- repObj$nP       # Number of stocks
    nX    <- 1               # Number of sex classes
    nF    <- repObj$nG + opMod$addProjFleets    # Number of fleets from SISCA + open pond SOK
    nS    <- 1               # Number of species
    nL    <- 1               # Number of length classes 

    om$A_s <- nA
    om$L_s <- 1

    om$lenBinWidth <- 2
    om$lenBinMids_l <- 1
  }

  # Bad prototype MICE structure
  nMICE <- 0
  isMICE  <- FALSE
  if(ctlList$opMod$MICE > 0)
  {
    isMICE  <- TRUE
    nMICE   <- ctlList$opMod$MICE
    nSpec   <- nS
    nS      <- nS + nMICE
    

    MICEspecies <- names(ctlList$opMod$MICEpars)

    MICE_lhPars <- lapply(X = ctlList$opMod$MICEpars, FUN = readRDS)
    names(MICE_lhPars) <- MICEspecies

    om$MICEpars <- MICE_lhPars

    for( mIdx in 1:nMICE )
    {
      sIdx <- nSpec + mIdx
      MICEpars <- om$MICEpars[[mIdx]]
      
      # Age classes
      om$A_s[nSpec + mIdx]        <- MICEpars$nA

    }
    nA <- max(om$A_s)
  }

  
  # Set up arrays to hold simulated states
  om$B_axspt    <- array(0,  dim = c(nA,nX,nS,nP,nT) )  # Biomass at age (total)
  om$N_axspt    <- array(0,  dim = c(nA,nX,nS,nP,nT) )  # Numbers at age
  om$endN_axspt <- array(0,  dim = c(nA,nX,nS,nP,nT) )  # Numbers at age (at end of time step)
  om$SB_spt     <- array(0,  dim = c(nS,nP,nT) )        # Spawning biomass
  om$endSB_spt  <- array(0,  dim = c(nS,nP,nT) )        # Spawning biomass + surviving ponded fish
  om$B_spt      <- array(0,  dim = c(nS,nP,nT) )        # Spawning biomass
  om$R_spt      <- array(0,  dim = c(nS,nP,nT) )        # Recruitment
  om$Z_axspt    <- array(0,  dim = c(nA,nX,nS,nP,nT) )  # Total mortality

  # Catch, fishing mort, vuln bio, selectivity
  om$C_spft     <- array(NA,  dim = c(nS,nP,nF,nT) )       # Total catch by fleet in kt
  om$P_spft     <- array(0,  dim = c(nS,nP,nF,nT) )        # Total ponded fish by fleet in kt
  om$psi_spft   <- array(0,  dim = c(nS,nP,nF,nT) )        # Conversion from ponded biomass to SOK
  om$P_axspft   <- array(0,  dim = c(nA,nX,nS,nP,nF,nT) )  # Total ponded fish at age (numbers)
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
  om$I_spft     <- array(NA,  dim = c(nS+1,nP+1,nF,nT) )   # Observations without error
  om$rI_spft    <- array(NA,  dim = c(nS+1,nP+1,nF,nT) )   # proportions for blended index
  om$E_pft      <- array(NA,  dim = c(nP,nF,nT) )          # Fishing effort
  om$alloc_spf  <- array(0,  dim = c(nS,nP,nF))            # Catch allocation
  om$Rev_spft   <- array(NA, dim = c(nS,nP,nF,nT) )        # Fleet revenue

  # Observation model pars
  om$q_spft     <- array(0,  dim = c(nS,nP,nF,nT) )        # catchability (tv)
  om$q_spf      <- array(0,  dim = c(nS,nP,nF) )           # catchability

  # Scalar between effort and F
  om$qF_spft    <- array(0,  dim = c(nS,nP,nF,nT) )        # catchability (tv)

  # Variance parameters
  om$tauObs_spf <- array(NA, dim = c(nS,nP,nF))
  om$tauAge_spf <- array(NA, dim = c(nS,nP,nF))
  om$sigmaR_sp  <- array(NA, dim = c(nS,nP))

  # Leading bio pars
  om$B0_sp      <- array(NA, dim = c(nS,nP))
  om$Rinit_sp   <- array(NA, dim = c(nS,nP))
  om$M_xsp      <- array(NA, dim = c(nX,nS,nP))
  om$m1_sp      <- array(NA, dim = c(nS,nP))
  om$h_sp       <- array(NA, dim = c(nS,nP))
  om$Rbar_sp    <- array(NA, dim = c(nS,nP))

  # Time varying natural mortality
  om$M_axspt          <- array(NA, dim = c(nA,nX,nS,nP,nT))
  om$initSurv_axsp    <- array(NA, dim = c(nA,nX,nS,nP) )

  # Spawn on kelp post pars
  om$pondM_ft   <- array(NA, dim = c(nF,nT))      # Post-ponding M
  om$gamma_f    <- array(NA, dim = c(nF) )        # SOK conversion factor
  om$pEff_t     <- array(.35, dim = c(nF,nT))     # Proportion effective spawners  


  # Growth model pars
  om$L1_s       <- array(NA, dim = c(nS))
  om$A1_s       <- array(NA, dim = c(nS))
  om$L2_xsp     <- array(NA, dim = c(nX,nS,nP))
  om$A2_s       <- array(NA, dim = c(nS))
  om$vonK_xsp   <- array(NA, dim = c(nX,nS,nP))
  om$LWa_s      <- array(NA, dim = c(nS))
  om$LWb_s      <- array(NA, dim = c(nS))

  # Weight-at-age
  om$W_axspt     <- array(NA, dim = c(nA,nX,nS,nP,nT))
  om$W_axspft    <- array(NA, dim = c(nA,nX,nS,nP,nF,nT))

  # Historical mean weight-at-age
  om$meanWtAge_axsp <- array(NA, dim = c(nA,nX,nS,nP))

  # maturity
  om$matAge_asp <- array(NA, dim = c(nA,nS,nP) )

   # Make data list - anything else?
  # mp list - mostly for retrospective
  # analyses and assessment model settings
  mp <- list( assess = NULL, data = NULL, hcr = NULL )
  mp$assess <- list()
  mp$hcr    <- list()
  mp$data <- list()

  # Data list
  mp$data$I_spft    <- array( NA, dim = c(nS+1,nP+1,nF,nT) )      # Add a pooled data dimension
  mp$data$A_axspft  <- array( NA, dim = c(nA,nX,nS+1,nP+1,nF,nT) )    
  mp$data$L_lxspft  <- array( NA, dim = c(nL,nX+1,nS,nP,nF,nT) ) 

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
  mp$hcr$sokEff_spft          <- array(0,  dim = c(nS,nP,nF,nT) )    # MP SOK license numbers


  # Process and observation errors
  errors <- list()

  # Arrays for holding errors
  errors$omegaR_spt       <- array(0, dim = c(nS,nP,nT) )         # recruitment PE
  errors$omegaRinit_asp   <- array(0, dim = c(nA,nS,nP) )         # intialisation recruitment PEs
  errors$delta_spft       <- array(0, dim = c(nS+1,nP+1,nF,nT) )  # index observation errors
  errors$ageObsErr_axspft <- array(0, dim = c(nA,nX,nS+1,nP+1,nF,nT) )
  errors$obsErrMult_spft  <- array(1, dim = c(nS,nP,nF,nT))



  if( ctlList$opMod$condModel == "hierSCAL" )
  {
    # Variance parameters
    om$tauObs_spf <- repObj$tauObs_spf
    om$sigmaR_sp  <- repObj$sigmaR_sp

    # Leading bio pars
    om$B0_sp      <- repObj$B0_sp
    om$Rinit_sp   <- repObj$R0_sp
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
  }

  if( ctlList$opMod$condModel == "SISCA" )
  {
    A <- om$A_s[1]
    # Variance parameters
    om$tauObs_spf[1,,1:repObj$nG] <- repObj$tauObs_pg
    om$tauAge_spf[1,,1:repObj$nG] <- sqrt(repObj$tau2Age_pg)
    om$sigmaR_sp[1,]              <- repObj$sigmaR

    # SOK pars
    om$gamma_f[1:repObj$nG]     <- repObj$gamma_g[1:repObj$nG]
    om$pEff_t[1:(tMP-1)]        <- repObj$pEff_t
  

    # Leading bio pars
    om$B0_sp[1,]    <- repObj$B0_p
    om$Rinit_sp[1,] <- repObj$Rinit_p
    om$M_xsp[1,1,]  <- repObj$M0_p
    om$m1_sp[1,]    <- repObj$m1_p
    om$h_sp[1,]     <- repObj$rSteepness_p
    om$Rbar_sp[1,]  <- repObj$Rbar_p

    om$densityDepM  <- repObj$densityDepM
    
    # initial survivorship
    om$initSurv_axsp[,1,1,]        <- repObj$initSurv_ap


    # Time-varying, age-structured natural mortality
    om$M_axspt[1:A,1,1,,1:(tMP-1)] <- repObj$M_apt[1:A,,1:(tMP-1)]

    # Growth model pars
    om$L1_s       <- NA
    om$A1_s       <- NA
    om$L2_xsp     <- NA
    om$A2_s       <- NA
    om$vonK_xsp   <- NA
    om$LWa_s      <- NA
    om$LWb_s      <- NA

    # Later: write calcSchedules function
    for( p in 1:nP)
      om$matAge_asp[1:A,1,p]         <- repObj$mat_a

    # Population weight at age etc
    om$W_axspt[1:A,1,1,,1:(tMP-1)]   <- repObj$W_apt[,,1:(tMP-1)]
    om$W_axspt[1:A,1,1,,tMP:nT]      <- repObj$meanWt_ap

    # Fleet weight-at-age etc.
    om$W_axspft[1:A,1,1,,1:repObj$nG,1:(tMP-1)] <- repObj$W_apgt[,,,1:(tMP-1)]
    om$W_axspft[1:A,1,1,,1:repObj$nG,tMP:nT]    <- aperm(repObj$meanWt_agp,c(1,3,2))

    # Fleet selectivity
    om$sel_axspft[1:A,1,1,,1:repObj$nG,1:(tMP-1)] <- repObj$sel_apgt[,,,1:(tMP-1)]
    for(t in tMP:nT)
    {
      om$sel_axspft[1:A,1,1,,1:repObj$nG,t] <- om$sel_axspft[1:A,1,1,,1:repObj$nG,tMP-1]
    }

    if(opMod$forceMeanWtAge)
    {
      for( t in 1:(tMP-1))
      {
        om$W_axspt[1:A,1,1,,t]               <- repObj$meanWt_ap      
        om$W_axspft[1:A,1,1,,1:repObj$nG,t]  <- aperm(repObj$meanWt_agp,c(1,3,2))
      }
    }

    if( nF > repObj$nG)
    {
      om$W_axspft[1:A,1,1,,repObj$nG+1,1:(tMP-1)]   <- om$W_axspft[,1,1,,repObj$nG,1:(tMP-1)]
      om$W_axspft[1:A,1,1,,repObj$nG+1,tMP:nT]      <- om$W_axspft[,1,1,,repObj$nG,tMP:nT]
    }
    # Get mean weight-at-age for reference point calcs
    om$meanWtAge_axsp[1:A,1,1,]      <- repObj$meanWt_ap

    # Get spawn timing and fleet timing
    om$spawnTiming_s              <- ctlList$opMod$spawnTiming_s
    om$spawnTiming_s[1]           <- repObj$spawnTiming
    om$fleetTiming_f              <- c(repObj$fleetTiming_g)
    om$fleetType_f                <- c(repObj$fleetType_g)

    if(opMod$addProjFleets > 0)
    {
      om$fleetTiming_f              <- c(om$fleetTiming_f,opMod$addFleetTiming)
      om$fleetType_f                <- c(om$fleetType_f,opMod$addFleetType)
    }      

    om$Wlen_ls                    <- NA
    om$probLenAge_laxsp           <- NA
    om$lenAge_axsp                <- NA

    # Overwrite data
    mp$data$I_spft[1,1:nP,,1:(tMP-1)] <- repObj$I_pgt[1:nP,,1:tMP-1]
    om$I_spft[1,1:nP,,1:(tMP-1)]      <- repObj$I_pgt[1:nP,,1:tMP-1]
  }

  # Now fill in MICE model species 
  if( isMICE )
  {

    # Loop over MICE model species
    for( mIdx in 1:nMICE )
    {
      sIdx <- nSpec + mIdx
      A <- om$A_s[sIdx]
      MICEpars <- om$MICEpars[[mIdx]]
  
      # Rec proc error
      om$sigmaR_sp[sIdx,]         <- MICEpars$sigmaR

      # Leading bio pars
      om$B0_sp[sIdx,]             <- MICEpars$B0
      om$Rinit_sp[sIdx,]          <- MICEpars$R0
      om$M_xsp[,sIdx,]            <- MICEpars$M
      om$m1_sp[sIdx,]             <- 0
      om$h_sp[sIdx,]              <- MICEpars$h
      om$Rbar_sp[sIdx,]           <- mean(MICEpars$N_at[1,])

      # Mortality
      om$M_axspt[,1,sIdx,,1:nT]   <- MICEpars$M 

      om$initSurv_axsp[1:A,1,sIdx,]  <- MICEpars$N_at[1:A,1]/MICEpars$R0

      for( p in 1:nP)
        om$matAge_asp[1:A,sIdx,p]         <- MICEpars$mat_a

      # Population weight at age etc
      wt_at     <- MICEpars$B_at/MICEpars$N_at
      wt_at     <- wt_at[,-ncol(wt_at)]
      meanWt_a  <- apply(X = wt_at, FUN = mean, MARGIN = 1)

      fYearPred <- MICEpars$tInitModel
      tInit     <- fYearPred - repObj$fYear + 1
      
      # HACK - fix hard-coded idx (t=7) later
      om$W_axspt[1:A,1,sIdx,1,tInit:(tMP-1)]  <- wt_at
      for(t in 1:(tInit-1))
        om$W_axspt[1:A,1,sIdx,1,t] <- meanWt_a
  
      # project forward with mean wt
      for( t in tMP:nT)
        om$W_axspt[1:A,1,sIdx,1,t]        <- meanWt_a
      # Copy wt-age to fleets
      for(f in 1:nF)
        om$W_axspft[,,sIdx,1,f,] <- om$W_axspt[,,sIdx,1,]

      # ref pt calcs
      om$meanWtAge_axsp[1:A,1,sIdx,1]      <- meanWt_a
    }
  }

  # Save model dims
  om$tMP   <- tMP
  om$nT    <- nT
  om$nX    <- nX
  om$nF    <- nF
  om$nA    <- nA
  om$nS    <- nS
  om$nP    <- nP
  om$nL    <- nL
  om$nMICE <- nMICE

  # Return initialised object
  obj <- list(  om = om, 
                mp = mp, 
                errors = errors,
                ctlList = ctlList )

  return(obj)
} # END initMS3pop()


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

  repObj      <- opMod$histRpt
  mcmcPar     <- ctlList$opMod$posts
  postDrawIdx <- obj$ctlList$opMod$postDrawIdx


  # Get model dimensions, and state variables
  tMP               <- om$tMP
  nT                <- om$nT
  nX                <- om$nX
  nF                <- om$nF
  nA                <- om$nA
  nS                <- om$nS
  nMICE             <- om$nMICE
  nP                <- om$nP
  nL                <- om$nL
  A_s               <- om$A_s
  L_s               <- om$L_s

  nSpec             <- nS - nMICE
  MICEpars          <- om$MICEpars

  # Length bin info
  lenBinWidth       <- om$lenBinWidth
  lenBinMids_l      <- om$lenBinMids_l

  # State vars
  B_axspt           <- om$B_axspt
  N_axspt           <- om$N_axspt
  endN_axspt        <- om$endN_axspt
  R_spt             <- om$R_spt
  SB_spt            <- om$SB_spt
  endSB_spt         <- om$endSB_spt
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
  Ierr_spft         <- data$I_spft
  Iperf_spft        <- om$I_spft

  # Age data
  ageObs_axspft     <- data$A_axspft

  # Do we need one unaffected by error? Isn't C_axspft this?
  # What about for surveys with forced age observations??
  # Get to it later - not an issue for Herring

  # Biological pars
  B0_sp             <- om$B0_sp
  Rbar_sp           <- om$Rbar_sp
  Rinit_sp          <- om$Rinit_sp
  M_xsp             <- om$M_xsp
  M_axspt           <- om$M_axspt
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
  W_axspft          <- om$W_axspft
  W_axspt           <- om$W_axspt
  Wlen_ls           <- om$Wlen_ls
  probLenAge_laxsp  <- om$probLenAge_laxsp
  lenAge_axsp       <- om$lenAge_axsp
  matAge_asp        <- om$matAge_asp

  # Fleet and spawn timing
  fleetTiming_f     <- om$fleetTiming_f
  fleetType_f       <- om$fleetType_f
  spawnTiming_s     <- om$spawnTiming_s

  # SOK stuff
  P_spft            <- om$P_spft
  psi_spft          <- om$psi_spft
  pondM_ft          <- om$pondM_ft

  # Compute pulse limit
  pulseMFrq   <- ctlList$opMod$pulseMFrq
  pulseLim_sp <- ctlList$opMod$pulseMfrac * B0_sp

  # Create a container to hold 
  # spawnN_aspx (for better array dim matching later)
  spawnN_axsp       <- array(0, dim = c(nA,nX,nS,nP) )

  # Add an initial time parameter for each SPECIES, not just 
  # stock
  tInit_sp <- om$tInit_sp



  # loop over stocks to calculate numbers for historical period
  for(p in 1:nP)
  {
    for( s in 1:nS )
    {
    # Initialise population
      if( tInit_sp[s,p] == t )
      {
      # Pull initial N multiplier in case of non-eq init
        if( !opMod$initUnfished )
          for( x in 1:nX )
          {
            N_axspt[1:A_s[s],x,s,p,t] <- obj$om$initSurv_axsp[1:A_s[s],x,s,p] * Rinit_sp[s,p] * exp(err$omegaRinit_asp[1:A_s[s],s,p])
          }
        
        if( opMod$initUnfished )
          for( x in 1:nX )
            N_axspt[1:A_s[s],x,s,p,t] <- obj$om$initSurv_axsp[1:A_s[s],x,s,p] * R0_sp[s,p]

        R_spt[s,p,t] <- N_axspt[1,1,s,p,t]

      }  

      # Otherwise update using total mortality
      if( t > tInit_sp[s,p] )
      {
      # Loop over species/pops
        for( a in 1:A_s[s] )
        {
          # Recruitment
          if( a == 1 )
          {
            # Fix Rt at historical values from rep file or posterior draws
            if(repObj$avgRcode_p[p]==1)
              lastIdx <- max(which(repObj$omegaR_pt[p,] != 0) )
            if(repObj$avgRcode_p[p]==0)
              lastIdx <- max(which(repObj$SRdevs_pt[p,] != 0) )

            if(!is.null(opMod$lastIdxOverwrite))
              lastIdx <- opMod$lastIdxOverwrite

            if(t<=lastIdx & opMod$fixCondRecruits )
            {

              # Fix Rt at historical values from rep file or posterior draws
              if( s <= nSpec )
              {
                if(!ctlList$opMod$posteriorSamples)
                  R_spt[s,p,t] <-  repObj$R_pt[p,t]

                if(ctlList$opMod$posteriorSamples)
                  R_spt[s,p,t] <-  mcmcPar$R_ipt[postDrawIdx,p,t]
              }

              if( s > nSpec )
              {
                mIdx          <- s - nSpec              
                R_spt[s,p,t]  <- MICEpars[[mIdx]]$N_at[1,t - tInit_sp[s,p] + 1]
              }

            }  

            # Simulate Beverton-Holt Recruitments for years where not estimated
            if(t>lastIdx | !opMod$fixCondRecruits)
            {  
              if(obj$ctlList$opMod$recType =='bevHolt')
              {
                R_spt[s,p,t] <- reca_sp[s,p] * SB_spt[s,p,t-1] / (1 + recb_sp[s,p] * SB_spt[s,p,t-1])

                if( !obj$ctlList$ctl$noProcErr)
                {
                  R_spt[s,p,t] <- R_spt[s,p,t] * exp( om$sigmaR_sp[s,p] * err$omegaR_spt[s,p,t] - 0.5*om$sigmaR_sp[s,p]^2)
                }

              }  

              # Average Recruitments
              if( !obj$ctlList$ctl$noProcErr & obj$ctlList$opMod$recType =='avgR')
              {
                  R_spt[s,p,t] <- Rbar_sp[s,p]

                if( !obj$ctlList$ctl$noProcErr)
                  R_spt[s,p,t] <- Rbar_sp[s,p] * exp( om$sigmaR_sp[s,p] * err$omegaR_spt[s,p,t] - 0.5*om$sigmaR_sp[s,p]^2) 


              }  
            }                      

            N_axspt[a,,s,p,t] <- R_spt[s,p,t] * 1/nX
          }
          # Apply mortality           
          if( a > 1 )
          {
            N_axspt[a,,s,p,t] <- endN_axspt[a-1,,s,p,t-1]
          }
          # Plus group
          if( a == A_s[s] )
          {
            N_axspt[a,,s,p,t] <- N_axspt[a,,s,p,t] + endN_axspt[a,,s,p,t-1]
          }

        } # END a loop
      } 
    } # END s loop
  } # END p loop


  # Now convert into biomass at the beginning of the time step
  B_axspt[,,,,t]    <- N_axspt[,,,,t] * W_axspt[,,,,t]

  # Compute total age2+ biomass at beginning of time step - updated to use juveMage, currently path specific
  B_spt[,,t]    <- apply( X = B_axspt[2:nA,,,,t,drop = FALSE], FUN = sum, MARGIN = 3:4, na.rm = T )


  # From here, we have enough info to start generating fishing effort
  # for the MICE model predators (i.e., Hake)
  if(t >= tMP )
  {
    predGears <- ctlList$opMod$predGears
    predType  <- ctlList$opMod$predType

    predGears <- predGears[predType == 2]

    for( mIdx in 1:nMICE )
    {
      sIdx <- nSpec + mIdx
      # Pull which fleets we are using (some redundancy here)
      predFleets <- MICEpars[[mIdx]]$predGearIdx

      # Now get bioenergetic model
      bioEnerModel <- MICEpars[[mIdx]]$bioEnergeticModel

      modelB_at <- array(0,dim = c(A_s[sIdx],1))
      modelB_at[,1] <- B_axspt[1:A_s[sIdx],1,sIdx,1,t]

      predList <- convertHakeBioToEffort( modelB_at = modelB_at,
                                          bioModel = bioEnerModel )

      E_pft[1:nP,predGears,t] <- predList$predEffort_lt[,1]

      # TODO: Here we need to update selectivity model with mean 
      # length
    }
  }

  # Density dependent M
  if( om$densityDepM == 1 )
  {
    totB0_sp <- rp$totB0_sp
    for(s in 1:nS)
      for(p in 1:nP)
      { 
        M_axspt[2:nA,1:nX,s,p,t] <- repObj$M_p[p] + exp( - om$m1_sp[s,p] * B_spt[s,p,t] / totB0_sp[s,p])
        M_axspt[1,1:nX,s,p,t] <- repObj$Mjuve_p[p]
      }
  }
  
  # Pulse Mt
  if( om$densityDepM == 0 & t >= tMP & pulseMFrq > 0 )
  {
    for(s in 1:nS)
      for(p in 1:nP)
      {
        if( SB_spt[s,p,t-1] < pulseLim_sp[s,p] ) 
          M_axspt[2:nA,,s,p,t] <- obj$om$pulseM_axspt[2:nA,,s,p,t]
      } 

  }

  # Now compute vuln biomass
  for( f in 1:nF )
    vB_axspft[,,,,f,t] <- sel_axspft[,,,,f,t] * B_axspt[,,,,t]

  vB_spft <- apply(X = vB_axspft, FUN = sum, MARGIN = c(3:6))

  # Now calculate effort for each area (for closed loop sim)
  if( ctlList$opMod$Ftype == "cts" )
  {
    

    vB_spft[,,,t] <- apply( X = vB_axspft[,,,,,t,drop = FALSE], FUN = sum, MARGIN = 3:5, na.rm = T )
    

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
        Rev_spft[,,,t] <- effortMod$rev_spf

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

        # Make temporary arrays to hold stuff
        tmpN_axsp   <- array(NA, dim = c(nA,nX,nS,nP))
        tmpM_axsp   <- array(NA, dim = c(nA,nX,nS,nP))
        tmpvB_axspf <- array(NA, dim = c(nA,nX,nS,nP,nF))
        tmpC_spf    <- array(NA, dim = c(nS,nP,nF))
        tmpvB_spf   <- array(NA, dim = c(nS,nP,nF))
        tmpsel_axspf<- array(NA, dim = c(nA,nX,nS,nP,nF))

        tmpN_axsp[1:nA,,,]    <- N_axspt[1:nA,,,,t] 
        tmpvB_axspf[1:nA,,,,] <- vB_axspft[1:nA,,,,,t]
        tmpC_spf[1:nS,,]      <- TAC_spft[1:nS,,,t]
        tmpvB_spf[1:nS,,]     <- vB_spft[,,,t]   
        tmpsel_axspf[1:nA,,,,]<- sel_axspft[,,,,,t]
        tmpM_axsp[1:nA,,,]    <- M_axspt[,,,,t]



        Flist <- .solveBaranov( C_spf       = tmpC_spf, 
                                vB_axspf    = tmpvB_axspf,
                                vB_spf      = tmpvB_spf,
                                N_axsp      = tmpN_axsp,
                                sel_axspf   = tmpsel_axspf,
                                M_axsp      = tmpM_axsp,
                                A_s         = A_s,
                                wt_axsp     = meanWtAge_axsp,
                                nS = nS, nP = nP, nX = nX, nF = nF,
                                nIter       = opMod$baranovIter,
                                baranovStep = opMod$baranovStep  )

        F_spft[,,,t] <- Flist$F_spf
      }
     

    } # END F approximation




    # Now generate Z
    for( s in 1:nS )
      for( p in 1:nP )
      {
        # Fill Z with Ms
        for( x in 1:nX )
          Z_axspt[1:A_s[s],x,s,p,t] <- M_axspt[,x,s,p,t]

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
    
    # Calculate numbers at-age for the end of the time step
    endN_axspt[,,,,t] <- N_axspt[,,,,t] * exp( - Z_axspt[,,,,t] )

    # Calculate numbers-at-age for spawning time
    for( s in 1:nS)
      spawnN_axsp[1:nA,,s,]  <- N_axspt[,,s,,t] * exp( - spawnTiming_s[s] * Z_axspt[,,s,,t] )

  }


  # Discrete removals
  if( ctlList$opMod$Ftype == "disc" )
  {
    # Need to apply a function that will step through
    # fleet timings, calculate removals of individuals at age
    # by splitting TAC in proportion to vulnerable biomass
    # at age, then removing them.
    tmpM_axsp     <- array(NA, dim = c(nA,nX,nS,nP))
    tmpN_axsp     <- array(NA, dim = c(nA,nX,nS,nP))
    tmpsel_axspf  <- array(NA, dim = c(nA,nX,nS,nP,nF))
    tmpW_axspf    <- array(NA, dim = c(nA,nX,nS,nP,nF))
    tmpP_spf      <- array(NA, dim = c(nS,nP,nF))

    tmpM_axsp[1:nA,,,]      <- M_axspt[,,,,t]
    tmpW_axspf[1:nA,,,,]    <- W_axspft[,,,,,t]
    tmpN_axsp[1:nA,,,]      <- N_axspt[,,,,t]
    tmpsel_axspf[1:nA,,,,]  <- sel_axspft[,,,,,t]
    tmpP_spf[1:nS,,]        <- P_spft[,,,t]

    sokFleets <- which(fleetType_f >= 2)


    # Now we need to calculate the catch for each predator fleet
    if(length(ctlList$opMod$predGears) & t >= tMP )
    {
      predGears <- ctlList$opMod$predGears
      predType  <- ctlList$opMod$predType

      tmpZ_axsp <- array(0, dim = c(nA,nX,nSpec,nP))

      # effortPred
      effortPredGears <- predGears[predType > 0]
      # Loop over predGears and calculate catch
      for(s in 1:nSpec )
        for( p in 1:nP )
        {
          tmpZ_axsp[1:A_s[s],,s,p] <- M_xsp[,s,p]
          # Calculate F
          for( f in effortPredGears )
          {
            F_spft[s,p,f,t] <- qF_spft[s,p,f,t] * E_pft[p,f,t]

            tmpZ_axsp[,,s,p] <- tmpZ_axsp[,,s,p] + sel_axspft[,,s,p,f,t] * F_spft[s,p,f,t]
          }
          
          # Now calc catch
          for(a in 1:A_s[s] )
            for( x in 1:nX )
              for( f in effortPredGears )
                Cw_axspft[a,x,s,p,f,t] <- (1 - exp( - tmpZ_axsp[a,x,s,p])) * vB_axspft[a,x,s,p,f,t] * F_spft[s,p,f,t] / tmpZ_axsp[a,x,s,p]
          

          # Use Baranov to calculate removals, assuming no comm
          # removals and cts fishing (needs work)
          TAC_spft[s,p,effortPredGears,t] <- apply(X = Cw_axspft[,,s,p,effortPredGears,t,drop = FALSE], FUN = sum, MARGIN = 5,na.rm = T)

          # Now calculat GWs TAC in spawn
          gwPredGears <- predGears[predType==0]
          TAC_spft[1,1,gwPredGears,t] <- E_pft[p,gwPredGears,t] * 0.00043 * exp(-3)
        }


    }


    if( t >= tMP & length(sokFleets) > 0 )
    {
      # Convert SOK TAC to ponded fish (if there is a TAC)
      for( fIdx in sokFleets)
      {
        if( all(TAC_spft[,,fIdx,t] == 0) )
          next

        # Temporary arrays to protect against
        # dropped dims
        tmpN_axsp    <- array(NA, dim = c(nA,nX,nS,nP))
        tmpsel_axsp  <- array(NA, dim = c(nA,nX,nS,nP))
        tmpW_axsp    <- array(NA, dim = c(nA,nX,nS,nP))

        tmpN_axsp[1:nA,,,]    <- N_axspt[,,,,t]
        tmpsel_axsp[1:nA,,,]  <- sel_axspft[,,,,fIdx,t]
        tmpW_axsp[1:nA,,,]   <- W_axspft[,,,,fIdx,t]

        SOKconvList <- calcSOKpsi(  N_axsp    = tmpN_axsp,
                                    sel_axsp  = tmpsel_axsp,
                                    mat_asp   = matAge_asp,
                                    pFem = .5,
                                    fec = 200,
                                    initF = 0.01,
                                    pEff = om$pEff_t[t],
                                    gamma = om$gamma_f[fIdx],
                                    W_axsp = tmpW_axsp,
                                    A_s = A_s,
                                    nA = nA,
                                    nX = nX,
                                    nS = nS,
                                    nP = nP  )
      
        psi_sp     <- SOKconvList$psi_sp
        psi_sp[!is.finite(psi_sp)] <- 1
        propMat_sp <- SOKconvList$propMat_sp
        propMat_sp[!is.finite(propMat_sp)] <- 0


        # Convert SOK product that is set in TAC to total ponded fish (some of which will survive)
        P_spft[,,fIdx,t]    <- TAC_spft[,,fIdx,t] / psi_sp
        psi_spft[,,fIdx,t]  <- psi_sp



        # convert arrays to spf dimensions for applyDiscreteFisheries()
        tmpP_spf[,,fIdx] <- P_spft[,,fIdx,t]

      }
    }
      
    

    if(t >= tMP & nMICE > 0 )
    {
      # Need to harvest some predators
      for(mIdx in 1:nMICE )
      {
        sIdx <- nSpec + mIdx
        for( p in 1:nP )
          TAC_spft[sIdx,p,2,t] <- vB_spft[sIdx,p,2,t] * .22

      }
    }



    # Now apply discrete fisheries 
    # repObj <<- repObj
    
    # if(t == tMP)
    # {
    #   debugFlag <<- TRUE
    #   browser()
    # }
    discRemList <- applyDiscreteFisheries(  N_axsp          = tmpN_axsp,
                                            sel_axspf       = sel_axspft[,,,,,t],
                                            W_axspf         = tmpW_axspf,
                                            fleetTiming_f   = fleetTiming_f,
                                            fleetType_f     = fleetType_f,
                                            spawnTiming_s   = spawnTiming_s,
                                            M_axsp          = tmpM_axsp,
                                            TAC_spf         = TAC_spft[,,,t],
                                            P_spf           = tmpP_spf,
                                            pondM_f         = pondM_ft[,t],
                                            nA = nA, nX = nX, nS = nS, nP = nP, nF = nF )

    # discRemList contains:
    #   - vuln biomass/numbers at the fleet timing
    #   - Catch at age in numbers and biomass
    #   - Spawn timing numbers at age
    #   - Exploitation rates, reported as Fs
    #   - end of time-step number at age (might do this for both cts and disc)


    # Pull vulnerable numbers and biomass at age
    vN_axspf          <- discRemList$vN_axspf
    vB_axspft[,,,,,t] <- discRemList$vB_axspf

    # Compute vulnerable biomass from what's returned and calculate exploitation rate
    vB_spft[,,,t] <- apply( X = vB_axspft[,,,,,t,drop = FALSE], FUN = sum, MARGIN = 3:5, na.rm = T )
    F_spft[,,,t]  <- TAC_spft[,,,t] / (vB_spft[,,,t] + 1e-9)

    # if(debugFlag)
    #   browser()

    for( fIdx in sokFleets)
      F_spft[,,fIdx,t] <- tmpP_spf[,,fIdx] * (1 - exp(-pondM_ft[fIdx,t])) / vB_spft[,,fIdx,t]
    
    # Pull spawn timing numbers at age
    spawnN_axsp[1:nA,,,] <- discRemList$spawnN_axsp

    # Pull end of time-step numbers
    endN_axspt[,,,,t] <- discRemList$endN_axsp

    # Pull catch at age in biomass and numbers
    Cw_axspft[,,,,,t] <- discRemList$Cw_axspf
    C_axspft[,,,,,t]  <- discRemList$C_axspf

    # if( t >= tMP )
    #   browser()

  }

  # Now calculate spawning biomass
  SB_asp <- array(NA, dim = c(nA,nS,nP))
  SB_asp[1:nA,1:nS,1:nP] <- matAge_asp[1:nA,1:nS,1:nP] * spawnN_axsp[1:nA,nX,1:nS,1:nP] * W_axspt[1:nA,nX,1:nS,1:nP,t]
  SB_spt[,,t] <- apply(X = SB_asp, FUN = sum, MARGIN = c(2,3), na.rm = T )

  # Calculate spawning biomass + surviving ponded fish after applying post-ponding M
  # survP_sp <- array(NA, dim = c(nS,nP))
  # for(s in 1:nS)
  #   for(p in 1:nP)
  #     survP_sp[s,p] <- sum(tmpP_spf[s,p,sokFleets]* exp(-pondM_ft[sokFleets,t]))
  
  # endSB_spt[,,t] <- SB_spt[,,t] + survP_sp


  endSB_asp <- SB_asp
  endSB_asp[1:nA,1:nS,1:nP] <- matAge_asp[1:nA,1:nS,1:nP] * endN_axspt[1:nA,nX,1:nS,1:nP,t] * W_axspt[1:nA,nX,1:nS,1:nP,t]
  endSB_spt[,,t] <- apply(X = endSB_asp, FUN = sum, MARGIN = c(2,3), na.rm = T )


  # Sum realised catch and dead ponded fish
  C_spft[,,,t] <- apply( X = Cw_axspft[,,,,,t,drop = FALSE], FUN = sum, MARGIN = c(3,4,5), na.rm = T )
  C_spt[,,t]   <- apply( X = C_spft[,,,t,drop = FALSE], FUN = sum, MARGIN = c(1,2), na.rm = T )

  # browser()


  # Now generate indices - if we are in the projection
  # or if historical data is generated by  OM
  tauObs_spf <- array(NA, dim=c(nS,nP,nF))
  
  for(s in 1:nS)
    for(p in 1:nP)
      tauObs_spf[s,p,]  <- om$tauObs_spf[s,p,1:nF] * err$obsErrMult_spft[s,p,1:nF,t]

  idxOn_spf   <- obj$mp$data$idxOn_spft[,,,t]
  if( t >= tMP | ctlList$mp$data$source == "OM" | ctlList$opMod$histType == "sim" )
  {

    if(!is.null(mp$assess$spTVqFleets))
      browser(cat('check time-varying q implemented correctly and q!=0 after tMP'))
    
    if(is.null(mp$assess$spTVqFleets) & is.null(opMod$blendIdx))
      q_spft[1:nS,1:nP,,t] <- q_spf[1:nS,1:nP,]

    # Calculate and save
    for( f in 1:nF)
    {
      # First, the species/stock specific observations
      for( s in 1:nS )
        for( p in 1:nP )
        {

          idxOn <- obj$mp$data$idxOn_spft[s,p,f,t]
          ageOn <- obj$mp$data$ageOn_spft[s,p,f,t]

          if(idxOn)
          {
            # Compute precision
            tau <- om$tauObs_spf[s,p,f] * err$obsErrMult_spft[s,p,f,t]

            # relative biomass survey
            if( ctlList$mp$data$idxType[f] == 1)
              Iperf_spft[s,p,f,t] <- q_spft[s,p,f,t] * vB_spft[s,p,f,t]

            # CPUE
            if( ctlList$mp$data$idxType[f] == 2 & C_spft[s,p,f,t] > 0 )
              Iperf_spft[s,p,f,t] <- C_spft[s,p,f,t] / E_pft[p,f,t] 

            # Spawn survey
            if( ctlList$mp$data$idxType[f] %in% c(3,4))
              Iperf_spft[s,p,f,t] <- q_spft[s,p,f,t] * SB_spt[s,p,t]

            # Save true observations and those with error
            Ierr_spft[s,p,f,t] <- Iperf_spft[s,p,f,t] * exp(tau * err$delta_spft[s,p,f,t] - 0.5 * tau^2)
          }

          if(ageOn)
          {
            # Need to generate age observation errors in the conditioning step
            if( C_spft[s,p,f,t] > 0)
            {
              for( x in 1:nX)
              {
                catAge    <- C_axspft[,x,s,p,f,t]
                zeroIdx   <- which(catAge == 0)
                logCatAge <- log(C_axspft[-zeroIdx,x,s,p,f,t])
                resids <- err$ageObsErr_axspft[-zeroIdx,x,s,p,f,t]

                logObsCatAge <- logCatAge + om$tauAge_spf[s,p,f] * resids
                obsCatAgeProp <- exp(logObsCatAge)/sum(exp(logObsCatAge))
                ageObs_axspft[-zeroIdx,x,s,p,f,t] <- obsCatAgeProp
                ageObs_axspft[zeroIdx,x,s,p,f,t]  <- 0
              }
            }

            # Now what if age observations are forced (i.e. survey with
            # no catch simulated) - then use vuln biomass or generate
            # a small amt of catch...
            if( C_spft[s,p,f,t] == 0 )
            {
              # Write survey age procedure here.
            }
          }
        } 
    }
  }

  # message("t = ",t,"\n")
  # message("hake SB = ", SB_spt[2,1,t], "\n")
  # message("hake vB = ", vB_spft[2,1,2,t], "\n")
  # message("hake Catch = ", C_spft[2,1,2,t], "\n")

  # Now make pooled data
  # Species Pooling
  if( ctlList$mp$data$speciesPooling & !ctlList$mp$data$spatialPooling  )
  {
    browser(cat('blended index not yet implemented for species pooling... \n'))

    for( f in 1:nF)
    {
      for(p in 1:nP )
      {
        tauObs_spf[tauObs_spf == 0] <- NA
        tau <- mean(tauObs_spf[,p,f],na.rm = TRUE)
        # Vuln bio
        if( ctlList$mp$data$idxType[f] == 1 )
          Iperf_spft[nS+1,p,f,t] <- sum( q_spft[,p,f,t] * vB_spft[,p,f,t] * idxOn_spf[,p,f], na.rm = T )

        # CPUE
        if( ctlList$mp$data$idxType[f] == 2 )
        {
          Iperf_spft[nS+1,p,f,t] <- max(idxOn_spf[,p,f]) * sum( C_spft[,p,f,t] ) 
          if( Iperf_spft[nS+1,p,f,t] > 0 )
            Iperf_spft[nS+1,p,f,t] <- Iperf_spft[nS+1,p,f,t] / max(1e-6,E_pft[p,f,t]) / 1e3
        }

        # Spawn survey
        if( ctlList$mp$data$idxType[f] == 3 )
        {
          Iperf_spft[nS+1,p,f,t] <- sum( q_spft[,p,f,t] * SB_spt[,p,t] * idxOn_spf[,p,f], na.rm = T )
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
        tau <- mean(tauObs_spf[s,1:nP,f],na.rm = TRUE)
        
        # Vuln bio
        if( ctlList$mp$data$idxType[f] == 1 )
          Iperf_spft[s,nP + 1,f,t] <- sum( q_spft[s,,f,t] * vB_spft[s,,f,t] * idxOn_spf[s,1:nP,f], na.rm = T )

        # CPUE
        if( ctlList$mp$data$idxType[f] == 2 )
        {
          browser(cat('check effort scaling is correct.... \n'))

          Iperf_spft[s,nP+1,f,t] <- max(idxOn_spf[s,1:nP,f]) * sum(  C_spft[s,,f,t] ) 
          if( Iperf_spft[s,nP+1,f,t] > 0 ) 
            Iperf_spft[s,nP+1,f,t] <- Iperf_spft[s,nP+1,f,t] / max(1e-6,sum( E_pft[,f,t], na.rm = TRUE )) / 1e3
        }

        # Spawn survey
        if( ctlList$mp$data$idxType[f] %in% c(3,4) )
        {
            Iperf_spft[s,nP + 1,f,t] <- sum( q_spft[s,,f,t] * SB_spt[s,,t] * idxOn_spf[s,1:nP,f], na.rm = T )

        }  
        
        # Spawn survey and vuln bio, aggregate spatial indices with error
        if( ctlList$mp$data$idxType[f] %in% c(1,3,4) )
          Ierr_spft[s,nP+1,f,t] <- sum(Ierr_spft[s,1:nP,f,t], na.rm=T)

        
        if( ctlList$mp$data$idxType[f] ==2  )
        {
          # assuming catch and effort are know without error, otherwise use line below:
          Ierr_spft[s,nP+1,f,t] <- Iperf_spft[s,nP+1,f,t] 

          # Ierr_spft[s,nP+1,f,t] <- Iperf_spft[s,nP+1,f,t] * exp( tau * err$delta_spft[s,nP+1,f,t] - 0.5 * tau^2)
        }  
          


      }
    }

    if(any(is.nan(Iperf_spft)))
      browser()
  }

  # Total aggregation
  if( ctlList$mp$data$spatialPooling & ctlList$mp$data$speciesPooling )
  {
    browser(cat('blended index not yet implemented for species & spatial pooling... \n'))
    
    for( f in 1:nF)
    {
      tauObs_spf[tauObs_spf == 0] <- NA
      tau <- mean(tauObs_spf[,,f],na.rm = TRUE)

      # Vuln bio
      if( ctlList$mp$data$idxType[f] == 1 )
        Iperf_spft[nS+1,nP + 1,f,t] <- sum( q_spft[,,f,t] * vB_spft[,,f,t] * idxOn_spf[1:nS,1:nP,f], na.rm = T )

      # CPUE
      if( ctlList$mp$data$idxType[f] == 2 )
      {
        Iperf_spft[nS+1,nP+1,f,t] <- max(idxOn_spf[,,f]) * sum( C_spft[,,f,t] ) 
        if( Iperf_spft[nS+1,nP+1,f,t] > 0 )
          Iperf_spft[nS+1,nP+1,f,t] <- Iperf_spft[nS+1,nP+1,f,t] / max(1e-6,sum( E_pft[,f,t], na.rm = TRUE )) / 1e3
      }

      # Spawn survey
      if( ctlList$mp$data$idxType[f] == 3 )
        Iperf_spft[nS+1,nP + 1,f,t] <- sum( q_spft[,,f,t] * SB_spt[,,t] * idxOn_spf[1:nS,1:nP,f], na.rm = T )

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
  endN_axspt        -> obj$om$endN_axspt
  R_spt             -> obj$om$R_spt
  SB_spt            -> obj$om$SB_spt
  endSB_spt         -> obj$om$endSB_spt
  B_spt             -> obj$om$B_spt
  Z_axspt           -> obj$om$Z_axspt
  M_axspt           -> obj$om$M_axspt

  # Catch
  C_axspft          -> obj$om$C_axspft
  Cw_axspft         -> obj$om$Cw_axspft
  F_spft            -> obj$om$F_spft
  C_spft            -> obj$om$C_spft
  C_spt             -> obj$om$C_spt
  E_pft             -> obj$om$E_pft
  Rev_spft          -> obj$om$Rev_spft
  TAC_spft          -> obj$mp$hcr$TAC_spft
  TAC_spt           -> obj$mp$hcr$TAC_spt
  
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
  ageObs_axspft     -> obj$mp$data$A_axspft

  # SOK stuff
  P_spft            -> obj$om$P_spft
  psi_spft          -> obj$om$psi_spft

  return( obj )
} # END ageSexOpMod()

#-------------------------------------------------------------------------------#
#--                   Harvest Control Rule Functions                          --#
#-------------------------------------------------------------------------------#

# .tacAlloc_hgRule()      
# Purpose:        allocate TAC to spawn on kelp fisheries first and 
#                 re-distribute remaining catch among fleets.   
# Parameters:     obj=list containing all variables and parameters
#                 necessary to compute the quota
# Returns:        obj with modified alloc_spf and TAC_spft
# Source:         B Doherty
.tacAlloc_hgRule <- function( obj, t)
{

  ctlList <- obj$ctlList # ctl list
  om      <- obj$om      # operating model

  TAC_p <- obj$mp$hcr$TAC_spt[,,t] # Maximum TAC for each area

  nS  <- om$nS
  nP  <- om$nP
  nT  <- om$nT
  nF  <- om$nF

  closedPondAlloc <- 0.0907 # 100 short tons (90.7 t or 0.0907 kt) allocated to each closed pond
  openPondAlloc   <- 0.0317 # 35 short tons (31.7t or 0.0317 kt) allocated to each open pond
  minTAC_f        <- ctlList$mp$hcr$minTAC_f # minimum ammount of TAC for each fleet to open
  fishF           <- ctlList$mp$hcr$fishF # fishing fleets

  # Create array to hold the number of licenses for open and closed ponds
  sokEff_pg     <- array(0, dim=c(nP,2)) # g=1 for closed ponds, g=2 for open ponds
  if(nF - 1 > 5)
  {
    sokEff_pg[,1] <- ctlList$mp$hcr$sokClosedEff_p
    sokEff_pg[,2] <- ctlList$mp$hcr$sokOpenEff_p
  }

  # array for fleet allocation of TAC
  tac_pf <- array(NA, dim=c(nP,nF))

  # Calculate min TAC needed for openings in each area
  minTAC_p     <- rep(closedPondAlloc,3)
  if(any(sokEff_pg[,2]>=1) )
      minTAC_p[sokEff_pg[,2]>=1] <- openPondAlloc

  # Initialise counter at 0 for number of licenses fishing
  nSOK_pg  <- array(NA, dim=c(nP,2))
  maxSOK_g <- ctlList$mp$hcr$maxSOK_g # max number of licenses
  
  if( nP > 1)
  {
    for (p in c(1,3,2))
    {

      remTAC <- TAC_p[p] # object to hold unallocated TAC

      # Adjust number of licenses available for SOK JP/S based on licenses allocated to C/S and Lou
      if(p ==2)
        for(g in 1:2)
          sokEff_pg[2,g] <- maxSOK_g[g] - sum(sokEff_pg[c(1,3),g])

      # Cumshewa/Selwyn: Max 1 SOK license and no commercial fisheries
      # Louscoone: Max 1 SOK license and no commercial fisheries
      # For JP/S, licenses not allocated in Cumshewa/Selwyn or Louscoone can be used
      maxTacSOK_g <- c(sokEff_pg[p,1]*closedPondAlloc, sokEff_pg[p,2]*openPondAlloc)

      # If not enough TAC for open or closed ponds, assign zero for all fleets
      if(TAC_p[p] < minTAC_p[p]  )
        tac_pf[p,] <- 0
      
      # If enough TAC, allocate to open pond, closed pond and seine roe, in that order
      if(TAC_p[p] >= minTAC_p[p]  ) 
      { 
        # allocate open ponds first
        if(fishF[7]==1 & sokEff_pg[p,2] >0)
        { 
            
          # if there is enough TAC, allocate maxTAC to all licenses
          if(maxTacSOK_g[2] <= TAC_p[p] )
            tac_pf[p,7] <- maxTacSOK_g[2]

          # if not enough then allocate all TAC for the area
          if(maxTacSOK_g[2] > TAC_p[p] )
            tac_pf[p,7]  <- TAC_p[p]

          # update available TAC
          remTAC <- remTAC - tac_pf[p,7]

        } # END loop allocating TAC to open ponds
        
        # allocate remaining TAC to closed ponds
        if(fishF[6]==1 & sokEff_pg[p,1] >0 )
        {
          # if there is enough TAC, allocate maxTAC to all licenses
          if(maxTacSOK_g[1] <= remTAC )
            tac_pf[p,6] <- maxTacSOK_g[1]

          # if not enough then allocate all remaining TAC
          if(maxTacSOK_g[1] > remTAC & remTAC > closedPondAlloc)
            tac_pf[p,6]  <- remTAC

          # if remaining TAC is not enough for one license allocate 0
          if(maxTacSOK_g[1] > remTAC & remTAC < closedPondAlloc)
            tac_pf[p,6]  <- 0

          # update available TAC
          remTAC <- remTAC - tac_pf[p,6]

        } # END loop allocating TAC to closed ponds

        # allocate remaining TAC to seine roe fleet in JP/S
        if(fishF[2]==1 & remTAC > minTAC_f[2] & p==2) 
        {

          tac_pf[p,2] <- remTAC

          # update available TAC
          remTAC <- remTAC - tac_pf[p,2]

        } # END loop allocating TAC to seine roe fleet

        # assign zeros for any remaining NAs
        tac_pf[p,][is.na(tac_pf[p,])] <- 0

        # Check if tac_pf exceed allowable TAC
        if(sum(tac_pf[p,]) >  TAC_p[p])
          browser(cat('HCR ERROR: more catch allocated to fleets than available... \n'))
          
      } # END loop for allocating TAC  

      # update number of licenses filled by area
      sokEff_pg[p,1] <- tac_pf[p,6]/closedPondAlloc
      sokEff_pg[p,2] <- tac_pf[p,7]/openPondAlloc


    } # END nP loop    

    # update SOK license, tac and fleet allocation
    obj$mp$hcr$sokEff_spft[1,,6:7,t]  <- sokEff_pg
    obj$mp$hcr$TAC_spft[1,,,t]        <- tac_pf
    obj$mp$hcr$TAC_spt[1,,t]          <- apply(tac_pf, FUN=sum, MARGIN=1)
  }

  if( nP == 1 )
  {
    obj$mp$hcr$TAC_spft[1,1,,t]       <- TAC_p * ctlList$opMod$projAlloc_f
    obj$mp$hcr$TAC_spt[1,1,t]         <- sum(TAC_p * ctlList$opMod$projAlloc_f)
  }

  if(any(is.na(TAC_p)))
    browser()

  return(obj)


} # END  .tacAlloc_hgRule

# .tac2SOK()      
# Purpose:        # Converts total TAC allocated to SOK to 
# Parameters:     obj=list containing all variables and parameters
#                 necessary to compute the quota
# Returns:        obj with modified alloc_spf and TAC_spft
# Source:         modfied from .calcHCR_ccRule, S.P. Cox
.tac2SOK <- function( obj, t)
{

  ctlList <- obj$ctlList
  psi     <- ctlList$mp$hcr$psi
  minSOK  <- ctlList$mp$hcr$minSOK

  om <- obj$om      # operating model
  
  nS  <- om$nS
  nP  <- om$nP

  # Check minimum amount of SOK catch is achieved, and if not assign  minSOK and re-distribute remaining catch among other fleets
  # if( sum(obj$mp$hcr$TAC_spft[,,6,t]) < 0.3 & sum(obj$mp$hcr$TAC_spt[,,t]) >= 0.3 )
  # {  
  #   # recalculate allocation  
  #   alloc_f      <- obj$om$alloc_spf[1,1,]
  #   alloc_f[6]   <- minSOK/sum(obj$mp$hcr$TAC_spt[,,t])
  #   alloc_f[1:5] <- alloc_f[1:5]/sum(alloc_f[1:5])*(1-alloc_f[6])

  #   # Use this allocation for all species and stocks
  #   for( s in 1:nS )
  #     for( p in 1:nP )
  #     {  
  #       obj$om$alloc_spf[s,p,] <- alloc_f
  #       obj$mp$hcr$TAC_spft[s,p,,t]    <- obj$mp$hcr$TAC_spt[s,p,t] * obj$om$alloc_spf[s,p,]
  #     }  

  # }  

  # TAC above represents total mortality for each gear type. The OM reads in SOK landed product from the TAC object so we need to convert TAC for spawn on kelp fishery from assumed dead ponded fish, based on DFO assumption of 65% mortality (Schweigert et al. 2018), to live ponded fish and then SOK product
  for( s in 1:nS )
      for( p in 1:nP )
      {      
        livePondedFish <- obj$mp$hcr$TAC_spft[s,p,6,t]/0.65
        browser(cat('line 1970....'))

        obj$mp$hcr$TAC_spft[s,p,6,t] <- livePondedFish*psi # convert to SOK product
      }

  return(obj)

} # END .tac2SOK

# .calcHCR_hgRule()      
# Purpose:        harvest control rule to generate the quota
# Parameters:     obj=list containing all variables and parameters
#                 necessary to compute the quota
# Requires:       global variable Bmin = initial value for lowest biomass
# Returns:        the quota=catch
# Source:         B Doherty, modified from .applyPrecHCR() 
.calcHCR_hgRule <- function( obj, t )
{

  hcr     <- obj$mp$hcr
  ctlList <- obj$ctlList
  assess  <- obj$mp$assess

  # Model dims
  nS  <- obj$om$nS # number of species
  nP  <- obj$om$nP # of stocks/areas
  nF  <- obj$om$nF # number of fleet
  nT  <- obj$om$nT # number of years

  # Pull reference points
  refPtList   <- obj$rp
  refPtType   <- ctlList$mp$hcr$refPtType

  # calculate projection time
  tMP <- obj$om$tMP
  pt  <- t - tMP + 1

  # Get HCR quantities
  Uref_sp <- obj$rp$FmsyRefPts$Umsy_sp
  Uref_sp[1,]  <- ctlList$mp$hcr$Uref_p[1:nP]

  # Create control points
  # Bref_sp <- obj$om$B0_sp # perfr
  Bref_sp          <- ctlList$mp$hcr$B0Wtd_p[1:nP]
  hcr$LCP_spt[,,t] <- Bref_sp * ctlList$mp$hcr$LCP
  hcr$UCP_spt[,,t] <- Bref_sp * ctlList$mp$hcr$UCP


  # Pull data
  I_spft <- obj$mp$data$I_spft
  C_spt  <- obj$om$C_spt
  idxF   <- obj$mp

  # Stock status
  projSB_p  <- assess$retroSB_tspt[pt,,,t] 

  # Use the last two years of data in Spawn surface Survey
  # and calculate mean over time to get splitting weights
  rctMeanI_sp <- apply( X = I_spft[,,5,(t-2):(t-1),drop = FALSE],
                        FUN = mean,
                        MARGIN = c(1,2), na.rm = T )

  # Calculate total recent mean index for a species
  rctMeanI_s <- apply( X = rctMeanI_sp, FUN = sum, MARGIN = 1, na.rm=T)

  # Create propTAC_sp array
  propTAC_sp <- array(NA, dim=c(nS,nP))

  # Private function to calculate ramped HCR
  calcRampedHCR <- function( B, LCP, UCP, Uref, lowUmult )
  {
    if( B < LCP )
      U <- lowUmult * Uref
    if( LCP <= B & B < UCP )
      U <- (lowUmult + (1 - lowUmult) * (B - LCP) / (UCP - LCP)) * Uref
    if( UCP <= B )
      U <- Uref

    U
  }


  # if spatialPooling=TRUE - MP for aggregated stocks (i.e. larger area)
  # if speciesPooling=TRUE - MP for aggregated species (i.e. mixed fishery)

  # Spatial MP for stocks and separate MP for each species
  if( !ctlList$mp$data$spatialPooling & !ctlList$mp$data$speciesPooling )
  {

    for ( s in 1:nS )
    {  
      for (p in 1:nP)
      {
        
        targU <- calcRampedHCR( B    = projSB_p[p] ,
                                LCP  = hcr$LCP_spt[s,p,t],
                                UCP  = hcr$UCP_spt[s,p,t],
                                Uref = Uref_sp[s,p],
                                lowUmult = 0)

        hcr$TAC_spt[s,p,t] <- targU*projSB_p[p]

      }
      
      # Calculate proportion of TAC for saving in hcr object
      propTAC_sp[s,] <- hcr$TAC_spt[s,,t]/sum(hcr$TAC_spt[s,,t])

    }      
  }

  # Aggregate MP for stocks and separate MP for each species
  if( ctlList$mp$data$spatialPooling & !ctlList$mp$data$speciesPooling )
  {

    for ( s in 1:nS )
    {  

        targU <- calcRampedHCR( B    = sum(projSB_p),
                                LCP  = sum(hcr$LCP_spt[s,,t]),
                                UCP  = sum(hcr$UCP_spt[s,,t]),
                                Uref = Uref_sp[s,min(2,nP)],
                                lowUmult = 0 )

        # Total TAC for all areas
        TAC <- targU*sum(projSB_p)

        # Allocated TAC among areas in proportion to spawning biomass
        propTAC_sp[s,] <- projSB_p/sum(projSB_p)

        hcr$TAC_spt[s,,t] <- TAC*propTAC_sp[s,]

    }
            

  }


  # Aggregate MP for stocks and aggregate MP for each species
  if( ctlList$mp$data$spatialPooling & ctlList$mp$data$speciesPooling )
  {
    browser(cat('Aggregate MP for species pooling not yet implemented'))
  }



  # Spatial MP for stocks and aggregate MP for each species
  if( !ctlList$mp$data$spatialPooling & ctlList$mp$data$speciesPooling )
  {
    browser(cat('Spatial MP for species pooling not yet implemented'))
  }

  # Apply smoother on TACs
  if( !is.null(ctlList$mp$hcr$maxDeltaTACup) )
  {
    # Determine which TACs are more than maxDeltaTAC
    # above last year
    maxDeltaTACup  <- ctlList$mp$hcr$maxDeltaTACup
    currTAC_sp     <- hcr$TAC_spt[,,t] %>% matrix(byrow=TRUE, nrow=nS,ncol=nP)
    lastTAC_sp     <- hcr$TAC_spt[,,t-1] %>% matrix(byrow=TRUE, nrow=nS,ncol=nP)


    # if last years TAC =0, then find last non-zero TAC
    if(any(lastTAC_sp==0) )
    { 
      for(s in 1:nS)
      {  
        pIdx <- which(lastTAC_sp[s,]==0)
        for( p in pIdx)
        {
          tac <- hcr$TAC_spt[s,p,1:t-1]
          tac <- tac[tac>0]
          lastTAC_sp[s,p] <- tac[length(tac)]
        }
      }    
    }    

    diffTAC_sp     <- currTAC_sp - lastTAC_sp
    DeltaTACup_sp  <- diffTAC_sp / lastTAC_sp 
    deltaDiff_sp   <- DeltaTACup_sp - maxDeltaTACup
    smoothIdx <- which( deltaDiff_sp > 0 )

    # Apply smoother
    currTAC_sp[smoothIdx] <- (1 + maxDeltaTACup) * lastTAC_sp[smoothIdx]
    hcr$TAC_spt[,,t] <- currTAC_sp
  }

  # Save proportion of TAC for use in retro SB plots
  hcr$propTAC_spt[,,t] <- propTAC_sp

  if(any(is.na(hcr$TAC_spt)))
    browser()

  # Now put all the lists we modified back into obj
  obj$mp$hcr    <- hcr
  obj$mp$assess <- assess

  return(obj)

}  # END .calcHCR_hgRule()


# .calcHCR_ccRule()      
# Purpose:        harvest control rule to generate the quota
# Parameters:     obj=list containing all variables and parameters
#                 necessary to compute the quota
# Requires:       global variable Bmin = initial value for lowest biomass
# Returns:        the quota=catch
# Source:         modfied from .calcHCR_ccRule, S.P. Cox
.calcHCR_ccRule <- function( obj, t )
{

  # Private function to find the minimum biomass over past 10 years
  # from which the stock has recovered.  
  find_min <- function(x)
  {
    # Goal: find the minumum value of x that has been "recovered" from.
    # In other words, find the smallest datapoint x[t] such that x[t] < x[t+1].
    # If NA is returned, then no such point has been found.

    # remove NAs
    x <- na.omit(x)
  
    if( length(x) == 1 ) return(NA)
  
    # propose first datapoint as current minimum
    min_curr <- x[1]
  
    # find another datapoint that is less than the current minimum
    # and also has a subsequent datapoint that is smaller than itself
    if( length(x)>2 ) {
      for( i in 2:(length(x)-1) ) {
        if( x[i]<min_curr && x[i]<x[i+1] ) min_curr <- x[i]
      }
    }
    # if the current minimum is still the first datapoint,
    # check the subsequent point
    if( min_curr==x[1] && x[1]>x[2] ) min_curr <- NA
    if(!is.na(min_curr)){
      if(min_curr==0) {
        min_curr <- 0.05
      }
    }
    min_curr
  }

  hcr     <- obj$mp$hcr
  ctlList <- obj$ctlList
  assess  <- obj$mp$assess

  # Model dims
  nS  <- obj$om$nS # number of species
  nP  <- obj$om$nP # of stocks/areas
  nF  <- obj$om$nF # number of fleet
  nT  <- obj$om$nT # number of years

  # Pull reference points
  refPtList   <- obj$rp
  refPtType   <- ctlList$mp$hcr$refPtType

  # calculate projection time
  tMP <- obj$om$tMP
  pt  <- t - tMP + 1

  # Get HCR quantities
  # Uref_sp <- hcr$Uref_spt[,,t]
  Uref_sp <- obj$rp$FmsyRefPts$Umsy_sp
  Uref_sp[1,]   <- ctlList$mp$hcr$Uref_p[1:nP]

  # calculate projection time
  tMP <- obj$om$tMP
  pt  <- t - tMP + 1

  # Stock status
  projSB_sp   <- assess$retroSB_tspt[pt,,,t] 
  projVB_sp   <- assess$retroVB_tspft[pt,,,2,t] 

  # Pull data
  I_spft <- obj$mp$data$I_spft
  C_spt  <- obj$om$C_spt
  idxF   <- obj$mp

  # Use the last two years of data in Spawn surface Survey
  # and calculate mean over time to get splitting weights
  rctMeanI_sp <- apply( X = I_spft[,,5,(t-2):(t-1),drop = FALSE],
                        FUN = mean,
                        MARGIN = c(1,2), na.rm = T )

  # Calculate total recent mean index for a species
  rctMeanI_s <- apply( X = rctMeanI_sp, FUN = sum, MARGIN = 1, na.rm=T)

  # Biomass over previous ten years
  I_spt <- I_spft[,,5,] # 5 is fleet ID for dive survey
  if( t > 10 ) 
    B10_spt <- I_spt[1:nS,1:nP,(t-10):(t-1),drop=F] + C_spt[,,(t-10):(t-1), drop=F]
  else         
    B10_spt <- I_spt[1:nS,1:nP,1:(t-1),drop=F] + C_spt[,,1:(t-1),drop=F]

  # Proposed Bmin for each species and stock area
  Bmin_sp <- array(dim=c(nS,nP))

  for (s in 1:nS) # species loop
    for (p in 1:nP) # stock loop  
        Bmin_sp[s,p] <- find_min(B10_spt[s,p,])

  # Create propTAC_sp array
  propTAC_sp <- array(NA, dim=c(nS,nP))


  # spatialPooling=TRUE - MP for aggregated stocks (i.e. larger area)
  # speciesPooling=TRUE - MP for aggregated species (i.e. mixed fishery)



  # Spatial MP for stocks and separate MP for each species
  if( !ctlList$mp$data$spatialPooling & !ctlList$mp$data$speciesPooling )
  {

    for ( s in 1:nS )
    {  
      for (p in 1:nP)
      {
        if( I_spt[s,p,t-1] > Bmin_sp[s,p] ) {
            hcr$TAC_spt[s,p,t] <- Uref_sp[s,p]*(I_spt[s,p,t-1]+C_spt[s,p,t-1])
        } else {
            hcr$TAC_spt[s,p,t] <- 0
          }
      }
      
      # Calculate proportion of TAC for saving in hcr object
      propTAC_sp[s,] <- hcr$TAC_spt[s,,t]/sum(hcr$TAC_spt[s,,t])

    }      
  }

  # Aggregate MP for stocks and separate MP for each species
  if( ctlList$mp$data$spatialPooling & !ctlList$mp$data$speciesPooling )
  {
    # If using aggregate MP, assign constant catch according to total Bmin and split based on weighting from indices
    
    Uref_s     <- vector(length=nS, 'numeric')
    for ( s in 1:nS )
    {  
      # Aggregate spawn index biomass, this will need to be adjusted for other survey types (e.g. CPUE)
      B10_t <- apply(X = B10_spt[s,,,drop=F], FUN = sum, MARGIN = 3, na.rm=T)
      Bmin <- find_min(B10_t)
      
      # Calculate proportion of biomass by area from index
      propTAC_sp[s,1:nP] <- rctMeanI_sp[1:nS,1:nP, drop=FALSE] / sum(rctMeanI_sp[s,1:nP])
      Bmin_p <- Bmin*propTAC_sp[s,1:nP]

      # Use index-weighted average for harvest rate
      Uref_s[s] <- sum(propTAC_sp[s,]*Uref_sp[s,])

      if( I_spt[s,nP+1,t-1] > Bmin) {
          aggregateTAC <- Uref_s[s]*( I_spt[s,nP+1,t-1] + sum(C_spt[s,1:nP,t-1]) )
      } else {
          aggregateTAC <- 0
      }

      
      hcr$TAC_spt[s,1:nP,t] <- aggregateTAC*propTAC_sp[s,]
    }  

  }


  # Aggregate MP for stocks and aggregate MP for each species
  if( ctlList$mp$data$spatialPooling & ctlList$mp$data$speciesPooling )
  {
    browser(cat('Aggregate MP for species pooling not yet implemented'))
  }



  # Spatial MP for stocks and aggregate MP for each species
  if( !ctlList$mp$data$spatialPooling & ctlList$mp$data$speciesPooling )
  {
    browser(cat('Spatial MP for species pooling not yet implemented'))
  }

  # Apply smoother on TACs
  if( !is.null(ctlList$mp$hcr$maxDeltaTACup) )
  {
    # Determine which TACs are more than maxDeltaTAC
    # above last year
    maxDeltaTACup  <- ctlList$mp$hcr$maxDeltaTACup
    currTAC_sp     <- hcr$TAC_spt[,,t] %>% matrix(byrow=TRUE, nrow=nS,ncol=nP)
    lastTAC_sp     <- hcr$TAC_spt[,,t-1] %>% matrix(byrow=TRUE, nrow=nS,ncol=nP)


    # if last years TAC =0, then find last non-zero TAC
    if(any(lastTAC_sp==0) )
    { 
      for(s in 1:nS)
      {  
        pIdx <- which(lastTAC_sp[s,]==0)
        for( p in pIdx)
        {
          tac <- hcr$TAC_spt[s,p,1:t-1]
          tac <- tac[tac>0]
          lastTAC_sp[s,p] <- tac[length(tac)]
        }
      }    
    }    

    diffTAC_sp     <- currTAC_sp - lastTAC_sp
    DeltaTACup_sp  <- diffTAC_sp / lastTAC_sp 
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

}  # END .calcHCR_ccRule()


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

  # Use the last two years of data in Spawn Survey
  # and calculate mean over time to get weights for allocating catch
  rctMeanI_sp <- apply( X = I_spft[1:nS,1:nP,5,(t-2):(t-1),drop = FALSE],
                        FUN = mean,
                        MARGIN = c(1,2), na.rm = T )

  # Calculate total recent mean index for a species
  rctMeanI_s <- apply( X = rctMeanI_sp, FUN = sum, MARGIN = 1)

  # if( t >= tMP )
  #   browser()

  
  # Compute target F
  hcr$targetF_spt[,,t] <- Fref_sp

  # if(t >= tMP)
  #   browser()
  
  # Now apply "F" as a harvest rate, since
  # we don't have M information
  hcr$TAC_spt[,,t] <- hcr$targetF_spt[,,t] * projSB_sp


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

  # Use the last two years of data in Spawn survey
  # and calculate mean over time to get
  # splitting weights
  rctMeanI_sp <- apply( X = I_spft[,,5,(t-2):(t-1),drop = FALSE],
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

  # Pooling TAC across multiple species & areas
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

  # Pooling TAC across multiple areas for separate species
  if( ctlList$mp$data$spatialPooling & !ctlList$mp$data$speciesPooling )
  {
    browser(cat('make sure this is not splitting TAC up among species s and aggregated species'))
    # Distribute among stocks
    propTAC_sp <- rctMeanI_sp
    for( s in 1:nS )
      propTAC_sp[s,] <- propTAC_sp[s,] / rctMeanI_s[s]

    hcr$TAC_spt[,,t] <- hcr$TAC_spt[,,t] * propTAC_sp
  }

  # Pooling TAC across multiple species for separate areas 
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



#-------------------------------------------------------------------------------#
#--                   Asessment functions                                     --#
#-------------------------------------------------------------------------------#

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
  if(ctlList$ctl$mpName=='NoFish')
    methodID <- 'perfectInfo'
  else
    methodID <- ctlList$mp$assess$method


  # AMfunName <- paste(".AM_",methodID, sep = "")

  if( methodID == "ItCt")
    obj <- .AM_ItCt( obj, t)

  if( methodID == "idxBased")
    obj <- .AM_idxBased( obj, t)

  if( methodID == "perfectInfo" )
    obj <- .AM_perfectInfo(obj,t)

  if( methodID == "hierProd" )
    obj <- .AM_hierProd(obj,t)

  if( methodID == "ISCAM" )
    obj <- .AM_ISCAM(obj,t)

  if( methodID == "SISCA" )
    obj <- .AM_SISCA(obj,t)  


  return(obj)
} # END .callProcedureAM()


# Apply the Integrated Statistical Catch-Age Model
# (ISCAM) to aggregate biomass (single-stock) data
# SOURCE: mseR-ISCAM for structure. 
.AM_ISCAM <- function( obj, t )
{

} # END .AM_ISCAM()




# Apply ItCt index AM - this uses spawning biomass index + catch in that year
.AM_ItCt <- function( obj, t )
{
  # Pull control list
  ctlList <- obj$ctlList

  # Get complex dims
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nF  <- obj$om$nF
  nT  <- obj$om$nT
  tMP <- obj$om$tMP

  # Calculate projection year
  pt <- t - tMP + 1
  pT <- ctlList$opMod$pT

  # Pull data
  I_spft <- obj$mp$data$I_spft
  q_spft <- obj$om$q_spft
  C_spt  <- obj$om$C_spt

  # HACK: Assign 0 instead of NAs
  I_spft[is.na(I_spft)] <- 0
  
  if(sum(is.na(C_spt[,,1:(t-1)]))>0)
    browser(cat('NAs in catch... \n'))

  # Get weighted q for surface and dive survey
  qWtd_p <- ctlList$mp$hcr$qWtd_p

  # Stock status, using weighted q from operating model grid
  projSB_spft <- I_spft[1:nS,1:nP,5,(tMP-1):(nT-1),drop = FALSE]/qWtd_p + C_spt[,,(tMP-1):(nT-1)]

  # Stock status, assuming perfect info for q
  # projSB_pt <- I_spft[1:nS,1:nP,5,(tMP-1):(nT-1)]/q_spft[1:nS,1:nP,5,(tMP-1):(nT-1)]  + C_spt[,,(tMP-1):(nT-1)]

  for(s in 1:nS)
  {  
    for (p in 1:nP)
    {
      
      obj$mp$assess$retroSB_tspt[pt:pT,s,p,t]    <- projSB_spft[s,p,1,pt:pT]
      obj$mp$assess$retroVB_tspft[pt:pT,s,p,5,t] <- projSB_spft[s,p,1,pt:pT]  

    }  

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
} # END .AM_ItCt()

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
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nF      <- obj$om$nF
  tMP     <- obj$om$tMP
  nMICE   <- obj$om$nMICE
  nSpec   <- nS - nMICE

  # Calculate projection year
  pt      <- t - tMP + 1
  pT      <- ctlList$opMod$pT

  # Get commGears
  commGears <- ctlList$opMod$commGears


  obj$mp$hcr$TAC_spt[1:nSpec,,t] <- 0
  obj$mp$hcr$TAC_spft[1:nSpec,,,t] <- 0

  tmpObj <- .ageSexOpMod( obj, t )


  for( s in 1:nS )
    for( p in 1:nP )
    {
      obj$mp$assess$retroSB_tspt[pt:pT,s,p,t]     <- tmpObj$om$SB_spt[s,p,t]
      obj$mp$assess$retroVB_tspft[pt:pT,s,p,,t]   <- tmpObj$om$vB_spft[s,p,,t] 
      

    }

  if(nMICE > 0)
  {

    obj$mp$hcr$TAC_spt[nSpec + 1:nMICE,,]   <- tmpObj$mp$hcr$TAC_spt[nSpec + 1:nMICE,,]
    obj$mp$hcr$TAC_spft[nSpec + 1:nMICE,,,] <- tmpObj$mp$hcr$TAC_spft[nSpec + 1:nMICE,,,]
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


#------------------------------------------------------------------------------#
#-- OM conditioning                                                 #
#------------------------------------------------------------------------------#

# .condMS3pop()
# Conditions OM for historical period from a hierSCAL or SISCA report
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
  obj$om$F_spft[,,,histdx]        <- repObj$F_spgt
  obj$om$C_spft[,,,histdx]        <- repObj$predCw_spft
  obj$om$C_spt[,,histdx]          <- apply( X = repObj$predCw_spft, FUN = sum, MARGIN = c(1,2,4) )
  obj$om$E_pft[,1:2,histdx]       <- repObj$E_pft[,1:2,]
  obj$om$sel_axspft[,,,,,histdx]  <- repObj$sel_axspft
  obj$om$sel_lspft[,,,,histdx]    <- aperm(repObj$sel_lfspt,c(1,3,4,2,5))

  obj$mp$hcr$TAC_spft[,,,histdx]  <- obj$om$C_spft[,,,histdx]
  obj$mp$hcr$TAC_spt[,,histdx]    <- obj$om$C_spt[,,histdx]

  obj$om$initSurv_axsp            <- repObj$surv_axsp

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
  obj$om$qF_spft[,,,tMP:nT]       <- obj$om$qF_spft[,,,tMP-1]


  # Now we have enough info to calculate reference points
  stime <- Sys.time()
  message(" (.condMS3pop) Calculating Fmsy and Emsy reference points\n")
  
  repObj$om <- obj$om
  repObj$condModel <- ctlList$opMod$condModel
  refPtList <- calcRefPts( repObj )
  obj$rp    <- refPtList

  etime <- Sys.time()
  rpTime <- round(etime - stime,2)
  message( paste(" (.condMS3pop) Reference point calculations completed in ", 
                  rpTime, " seconds\n", sep = "" ) )

  # Add historical data
  obj$mp$data$I_spft[1:nS,1:nP,,histdx] <- repObj$I_spft
  obj$om$I_spft[1:nS,1:nP,,histdx]      <- repObj$I_spft
  obj$mp$data$A_axspft[,,,,,histdx]  <- aperm( repObj$age_aspftx, c(1,6,2:5) )
  obj$mp$data$L_lxspft[,,,,,histdx]  <- aperm( repObj$len_lspftx, c(1,6,2:5) )

  # replace negatives with NAs for plotting, can change back 
  # later for TMB
  obj$mp$data$I_spft[obj$mp$data$I_spft<0] <- NA
  obj$mp$data$A_axspft[obj$mp$data$A_axspft<0] <- NA
  obj$mp$data$L_lxspft[obj$mp$data$L_axspft<0] <- NA



  # Now calculate total catch and allocation
  tdxRecent <- (tMP - ctlList$opMod$allocYears):(tMP - 1)
  recentCatch_spf <- apply( X = obj$om$C_spft[,,,tdxRecent], FUN = sum, MARGIN = c(1,2,3) )

  commGears <- ctlList$opMod$commGears

  # Loop and fill errors and allocation
  for( s in 1:nS )
    for( p in 1:nP )
    {
      obj$errors$delta_spft[s,p,,] <- rnorm(nF*nT)

      if( !ctlList$ctl$noProcErr )
        obj$errors$omegaR_spt[s,p,] <- rnorm(nT)

      # Save historical proc errors, but use simulated recruitments after
      # last estimated recruitment
      lastDevIdx <- max(which(repObj$omegaR_spt[s,p,] != 0) )

      if(!is.null(opMod$lastIdxOverwrite))
        lastDevIdx <- opMod$lastIdxOverwrite


      if(lastDevIdx > 0)
      {

        obj$errors$omegaR_spt[s,p,histdx[1:lastDevIdx]]   <- repObj$omegaR_spt[s,p,1:lastDevIdx] 
      
        if( !ctlList$ctl$noProcErr )
          obj$errors$omegaR_spt[s,p,histdx[1:lastDevIdx]] <- obj$errors$omegaR_spt[s,p,histdx[1:lastDevIdx]] + 0.5*repObj$sigmaR_sp[s,p]  # rec devs    
      }

      obj$om$alloc_spf[s,p,commGears] <- recentCatch_spf[s,p,commGears] / sum( recentCatch_spf[s,p,commGears])


      # First, fill with uncorrelated errors
      obj$errors$ageObsErr_axspft[,,s,p,,] <- rnorm(nA*nX*nF*nT)

      # # Then check if age observation errors are correlated
      # browser(cat('line4148'))

    }


  # Replace NaNs with 0s when the recent history has no catch
  obj$om$alloc_spf[is.nan(obj$om$alloc_spf)] <- 0

  obj$errors$omegaRinit_asp         <- repObj$omegaRinit_asp # Initialisation errors
  if(!is.null(opMod$initDevsMult))
    obj$errors$omegaRinit_asp <- errors$omegaRinit_asp * opMod$initDevsMult
  
  obj$errors$delta_spft[nS+1,1:nP,,] <- rnorm(nT * nP * nF)
  obj$errors$delta_spft[1:(nS+1),nP+1,,] <- rnorm(nT * (nS + 1)*nF)

  # Save historical errors
  if( ctlList$mp$data$source == "cond")
    obj$errors$delta_spft[1:nS,1:nP,,histdx]  <- repObj$residCPUE_spft # obs errors
   
  
  

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

  obj <- .calcTimes( obj )

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


} # END .condMS3pop()

# .condMS3pop_SISCA()
# Conditions OM for historical period from a hierSCAL report
# object.
.condMS3pop_SISCA <- function( obj )
{

  ctlList <- obj$ctlList

  if(!obj$ctlList$ctl$quiet)
    message(" (.condMS3pop_SISCA) Conditioning ms3R from SISCA report\n")

  # Get model historical period report object
  repObj    <- ctlList$opMod$histRpt
  MICEpars  <- obj$om$MICEpars


  # Get model dims
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nA      <- obj$om$nA

  nMICE   <- obj$om$nMICE
  nSpec   <- nS - nMICE

  A_s <- obj$om$A_s

  nX  <- obj$om$nX
  tMP <- obj$om$tMP
  pT  <- ctlList$opMod$pT

  # Now, copy model fit
  histdx <- 1:(tMP - 1)
  histF  <- 1:repObj$nG # number of fleets in historical time series
  projdx <- tMP:nT

  obj$om$F_spft[1,,histF,histdx]        <- repObj$F_pgt
  obj$om$C_spft[1,,histF,histdx]        <- repObj$totC_pgt

  # Identify time step for initializing each population
  obj$om$tInit_sp <- array(1,dim = c(nS,nP))
  # TMB models have 0 indexing
  obj$om$tInit_sp[1,]  <- repObj$tInitModel_p + 1


  # Load predators  with predType == 1
  if(!is.null(ctlList$opMod$predGears))
  {
    predGears <- ctlList$opMod$predGears
    predType  <- ctlList$opMod$predType

    projPred  <- readRDS("history/MMpredMod_postProj.rds")
    postDraw  <- sample(1:1000, 1)

    # First, do predType == 1
    predGears <- predGears[predType == 1]

    Npred_ift <- aperm(projPred$N_itp[,,c(4,3,1:2)],perm = c(1,3,2))

    meanNpred_ft <- apply(X = Npred_ift, FUN = mean, MARGIN = 2:3)

    # First calculate predator qs
    obj$om$F_spft[1,1:nP,predGears,histdx] <- repObj$apF_pgt[1:nP,predGears,histdx]
    obj$om$E_pft[1:nP,predGears,histdx]    <- meanNpred_ft[,histdx]/1e3

    for(s in 1:nSpec )
      for(p in 1:nP )
      {
        qSeries_ft <- obj$om$F_spft[s,p,predGears,histdx] / (Npred_ift[postDraw,,histdx]/1e3)
        qbar_f     <- exp(apply(X = log(qSeries_ft), FUN = quantile, MARGIN = 1, probs = .5))
        # obj$om$qF_spf[s,p,predGears] <- qbar_f
        for( t in tMP:nT)
          obj$om$qF_spft[s,p,predGears,t] <- qbar_f
      }

    obj$om$E_pft[1:nP,predGears,projdx] <- Npred_ift[postDraw,,projdx]/1e3


    # predType == 2 indicates that the predator is a simulated
    # species in the MICE model, so consumption will be generated 
    # in closed loop sims year by year. However, we can calculate the
    # predator q from historical data

    # These two lines mean nothing, might use them as a safety check later
    # so we can turn predation on and off.
    predGears <- ctlList$opMod$predGears
    predGears <- predGears[predType == 2]

    # Now loop over MICE inputs and generate predator qs
    for( mIdx in 1:nMICE )
    {
      # Get gear indices
      gearIdx <- MICEpars[[mIdx]]$predGearIdx
      # Get historical period biomass for predators
      histEffort_lt <- MICEpars[[mIdx]]$bioEnergeticModel$piscB_lt
      # Pad historical effort by copying first element
      tInitModel <- MICEpars[[mIdx]]$tInitModel - 1951 + 1

      obj$om$E_pft[1:nP,gearIdx,tInitModel:(tMP-1)] <- histEffort_lt
      for(t in 1:(tInitModel-1))
        obj$om$E_pft[1:nP,gearIdx,t] <- histEffort_lt[,1]

      obj$om$F_spft[1,1:nP,gearIdx,histdx] <- repObj$apF_pgt[1:nP,predGears,histdx]

      for(s in 1:nSpec )
        for(p in 1:nP )
        {
          qSeries_ft <- obj$om$F_spft[s,p,gearIdx,histdx] / obj$om$E_pft[p,gearIdx,histdx]
          qbar_f     <- exp(apply(X = log(qSeries_ft), FUN = quantile, MARGIN = 1, probs = 0.5))
          # obj$om$qF_spf[s,p,predGears] <- qbar_f
          for( t in tMP:nT)
            obj$om$qF_spft[s,p,gearIdx,t] <- qbar_f
        }

    }

    # predType == 0 is for Gray Whales right now. There has to 
    # be a better way...
    predGears <- ctlList$opMod$predGears
    predGears <- predGears[predType == 0]

    # Load GW projection model
    gwList    <- readRDS("history/GrayWhales/GW_proj.rds")
    gwN_it    <- gwList$projN_it
    nTrials   <- dim(gwN_it)[1]
    nGWprojT  <- dim(gwN_it)[2]

    if(nGWprojT - 5 < ctlList$opMod$pT )
    {
      cat("Error: projection longer than provided GW abundance series")
    }
    drawIdx  <- sample(1:nTrials, 1)
    # Gray whale "effort" should be included in report
    # object...
    
    # Initial GW projection year is 2015
    projTdx <- seq(from = 6, by = 1, length.out = ctlList$opMod$pT)
    # Abundance is effort
    obj$om$E_pft[1:nP,13,(tMP):nT] <- gwN_it[drawIdx,projTdx]

    obj$om$E_pft[1:nP,13,17:(tMP-6)] <- gwList$histN
    obj$om$E_pft[1:nP,13,1:16] <- mean(obj$om$E_pft[1:nP,13,17:19])
    obj$om$E_pft[1:nP,13,(tMP-5):(tMP-1)] <- apply(X = gwN_it[,1:5], FUN = mean, MARGIN = 2)

    obj$om$qF_spft[1,1:nP,13,] <- repObj$qF_g[13]

    # Scale the effort by egg weight consumed by predator, and the 
    # repObj q value to give the total consumption, skipping a relative 
    # F - can we respec the model to make it a relative F in the history??
    obj$mp$hcr$TAC_spft[1,1:nP,13,1:nT] <- obj$om$E_pft[1:nP,13,1:nT] * obj$om$qF_spft[1,1:nP,13,(1:nT)] * 4.3e-4 
    obj$om$C_spft[1,1:nP,13,1:nT] <- obj$om$E_pft[1:nP,13,1:nT] * obj$om$qF_spft[1,1:nP,13,(1:nT)] * 4.3e-4 
    
  }

  # We should do ALL of the sok/GW history stuff in ONE place.
  sokIdx <- which(repObj$fleetType_g %in% 2:3)

  # For the historical period, we need ponded biomass equivalent for the 
  # applyDiscreteFisheries() function
  obj$om$P_spft[1,1:nP,sokIdx,1:(tMP-1)] <- repObj$pondC_pgt[1:nP,sokIdx,1:(tMP-1)]
  obj$om$C_spft[1,1:nP,sokIdx,1:(tMP-1)] <- repObj$C_pgt[1:nP,sokIdx,1:(tMP-1)]

  # TODO: ADD IN GW CATCH FROM MIXED CATCH
  obj$om$C_spft[1,1,13,1:(tMP-1)] <- obj$om$C_spft[1,1:nP,13,1:(tMP-1)] + repObj$mC_gt[13,1:(tMP-1)]

  # Simulate future post-ponding mortality
  for( f in sokIdx )
  {
    # Post pond M
    obj$om$pondM_ft[f,1:(tMP-1)] <- repObj$postPondM_g[f]
    
    # projections
    if(obj$om$fleetType_f[f] == 3)
      obj$om$pondM_ft[f,tMP:nT] <- 0

    # Draw M in projections from log-normal distribution
    if( obj$om$fleetType_f[f] == 2 )
    {
      pondM <- obj$ctlList$opMod$pondMmu_f[f]
      logSD <- obj$ctlList$opMod$pondMlogSD_f[f]

      obj$om$pondM_ft[f,tMP:nT] <- exp(rnorm(pT, log(pondM), logSD))
    }

    for( p in 1:nP)
    {
      obj$om$P_spft[1,p,f,histdx]   <- repObj$pondC_pgt[p,f,histdx]
      
      if(any(repObj$psi_pgt[p,f,histdx] != 0))
        obj$om$psi_spft[1,p,f,histdx] <- repObj$psi_pgt[p,f,histdx]

      if(all(repObj$psi_pgt[p,f,histdx] == 0))
        obj$om$psi_spft[1,p,f,histdx] <- repObj$psi_gt[f,histdx]
    }
  }  

  
  # Total catch
  obj$om$C_spt[1,,histdx]               <- apply( X = repObj$totC_pgt, FUN = sum, MARGIN = c(1,3) )
  
  # Fill selectivity for first species (target)
  for( t in histdx)
    obj$om$sel_axspft[1:A_s[1],1,1,1:nP,histF,t] <- repObj$sel_apgt[1:A_s[1],1:nP,histF,t]

  # HACK: put a control file par in for this, too path specific
  # if(nF > repObj$nG)
  #   obj$om$sel_axspft[,1,1,,7,histdx]     <- obj$om$sel_axspft[,1,1,,6,histdx]

  # Update selectivity for MICE species
  if(nMICE > 0)
    for( mIdx in 1:nMICE)
    {
      # Read in catch selectivity
      sIdx <- nSpec + mIdx
      initYear <- MICEpars[[mIdx]]$tInitModel
      initT    <- initYear - 1951 + 1
      # Fill selectivity
      for( t in initT:tMP )
        obj$om$sel_axspft[,1,sIdx,1,2,t] <- MICEpars[[mIdx]]$sel_at[,t - initT + 1]
      
      for( t in tMP:nT )
        obj$om$sel_axspft[,1,sIdx,1,2,t] <- obj$om$sel_axspft[,1,sIdx,1,2,tMP]

      # Initial time-step
      obj$om$tInit_sp[nSpec + mIdx,] <- initT

      nYrs <- length(initT:(tMP-1))

      obj$om$C_spft[sIdx,,,]  <- 0
      obj$om$C_spt[sIdx,,]    <- 0
      obj$om$C_spft[sIdx,1,2,initT:(tMP-1)] <- MICEpars[[mIdx]]$C_t[1:nYrs]
      obj$om$C_spt[sIdx,1,initT:(tMP-1)] <- MICEpars[[mIdx]]$C_t[1:nYrs]

    }

  obj$mp$hcr$TAC_spft[,,histF,histdx]   <- obj$om$C_spft[,,histF,histdx]

  obj$mp$hcr$TAC_spt[,,histdx]          <- obj$om$C_spt[,,histdx]
  
  # Put spawn removals into TAC object
  obj$mp$hcr$TAC_spft[1,,sokIdx,histdx] <- repObj$C_pgt[,sokIdx,]


  

  
  # Now fill in future sel
  if(!ctlList$opMod$posteriorSamples)
  {  
    # selectivity
    obj$om$sel_axspft[,1,1,,histF,tMP:nT]  <- obj$om$sel_axspft[,1,1,,histF,tMP-1]
    if(nF > repObj$nG)
      obj$om$sel_axspft[,1,1,,7,tMP:nT]      <- obj$om$sel_axspft[,1,1,,6,tMP:nT]

    # catchability
    obj$om$q_spf[1,,histF]    <- repObj$qhat_pg # same as repObj$qComb_pg
    obj$om$q_spft[1,,histF,]  <- repObj$qhat_pg 
    

  }
  
  if(ctlList$opMod$posteriorSamples)
  {
    postDrawIdx <- obj$ctlList$opMod$postDrawIdx
    mcmcPar     <- ctlList$opMod$posts

    # selectivity
    sel50_pf <- mcmcPar$selAlpha_ipg[postDrawIdx,,]
    sel95_pf <- sel50_pf + mcmcPar$selBeta_ipg[postDrawIdx,,]
    tmp_pf      <- log(19.)/( sel95_pf - sel50_pf )

    sel_apf <- array(0, dim=c(nA,nP,length(histF)))
    # Force selectivity to zero for ages 1:(minAge_g-1)
    for( fIdx in histF)
      for(a in repObj$minAge_g[fIdx]:nA)
        sel_apf[a,,fIdx] = 1/( 1 + exp(-tmp_pf[fIdx]*( a - sel50_pf[,fIdx])))

    obj$om$sel_axspft[,1,1,,histF,] <- sel_apf
    if( nF > repObj$nG)
      obj$om$sel_axspft[,1,1,,7,]     <- obj$om$sel_axspft[,1,1,,6,]

    # catchability
    obj$om$q_spf[1,,histF]    <- mcmcPar$q_ipg[postDrawIdx,,]
    obj$om$q_spft[1,,histF,]  <- mcmcPar$q_ipg[postDrawIdx,,]

  }  


  # Now we have enough info to calculate reference points
  stime <- Sys.time()
  message(" (.condMS3pop_SISCA) Calculating Fmsy and Emsy reference points\n")
  
  repObj$om <- obj$om
  repObj$condModel <- ctlList$opMod$condModel
  refPtList <- calcRefPts( repObj )
  obj$rp    <- refPtList

  etime <- Sys.time()
  rpTime <- round(etime - stime,2)
  message( paste(" (.condMS3pop_SISCA) Reference point calculations completed in ", 
                  rpTime, " seconds\n", sep = "" ) )


  # Add historical data
  for(f in 1:repObj$nG)
  {
    obj$mp$data$I_spft[1,1:nP,f,histdx]        <- repObj$I_pgt[1:nP,f,histdx]
    obj$mp$data$A_axspft[,1,1,1:nP,f,histdx]   <- repObj$A_apgt[,1:nP,f,histdx]

    # Check link between age sample size and other data
    Asample_pt <- apply(X =repObj$A_apgt[,,f,histdx,drop = FALSE], FUN = sum, MARGIN = c(2,4), na.rm = T )
    Asample_pt[Asample_pt < 0] <- 0
    C_pt       <- apply(X = repObj$C_pgt[,,histdx,drop = FALSE], FUN = sum, MARGIN = c(1,3))

    

  }



  # Check if model is using a blended index & overwrite dive survey index & q
  if(!is.null(repObj$combI_pt))
    if( any(repObj$combI_pt > 0))
    {  
      
      # overwrite dive survey index with blended index
      obj$mp$data$I_spft[1,1:nP,5,histdx] <- repObj$combI_pt

      # Flag the om rep file is using blended index
      if(!ctlList$opMod$blendIdx)
        browser(cat('OM rep file uses blended index but ctlList does not... \n'))

    }  

  # Simulate proportions of blended index
  if(ctlList$opMod$blendIdx)
  {
    # add historical proportions of surface and survey index
    obj$om$rI_spft[1,1:nP,histF,histdx]  <- repObj$rI_pgt

    # add simulated proportions for dive survey
    blendG      <- 
    blendType_p <- ctlList$opMod$blendType_p
    binom_p     <- ctlList$opMod$binom_p
    beta1_p     <- ctlList$opMod$betaShape1_p
    beta2_p     <- ctlList$opMod$betaShape2_p

    for(p in 1:nP)  
    {
      if(blendType_p[p]=='binom')
        obj$om$rI_spft[1,p,4,projdx] <- rbinom(pT,1,binom_p[p])

      if(blendType_p[p]=='binomXbeta')
      {
        surfIdxOn  <- rbinom(pT,1,binom_p[p])
        surfYrs    <- which(surfIdxOn==1) + projdx[1] -1

        # random draws from beta distribution
        betaDraws <- rbeta(length(surfYrs), shape1=beta1_p[p], shape2=beta2_p[p])

        obj$om$rI_spft[1,p,4,projdx]  <- surfIdxOn
        obj$om$rI_spft[1,p,4,surfYrs] <- betaDraws

      }

      # Add proportions for dive survey
      obj$om$rI_spft[1,p,5,projdx] <- 1 - obj$om$rI_spft[1,p,4,projdx]

      for( tt in projdx )
        obj$om$q_spft[1,p,5,tt] <- sum(repObj$qComb_pg[p,4:5] * obj$om$rI_spft[1,p,4:5,tt])

    }   
  }  


  # replace negatives with NAs for plotting, can change back 
  # later for TMB
  obj$mp$data$I_spft[obj$mp$data$I_spft<0] <- NA
  obj$mp$data$A_axspft[obj$mp$data$A_axspft<0] <- NA
  # obj$mp$data$L_lxspft[obj$mp$data$L_axspft<0] <- NA

  # Now calculate total catch and allocation
  
  # Allocate according to historical catch by area
  # recentCatch_spf <- apply( X = obj$om$C_spft[,,,tdxRecent,drop = FALSE], FUN = sum, MARGIN = c(1,2,3) )

  # Allocate according to historical catch for all areas
  if(!is.null(obj$ctlList$opMod$projAlloc_f))
    for(s in 1:nS)
      for( p in 1:nP )
        obj$om$alloc_spf[s,p,1:nF] <- obj$ctlList$opMod$projAlloc_f[1:nF]

  if(is.null(obj$ctlList$opMod$projAlloc_f))
  {
    
    # Determin last 20 years when fishery was open
    recentCatch_t <- apply( X = obj$om$C_spft, FUN = sum, MARGIN = c(4) )

    openIdx   <- which(recentCatch_t>0)
    tdxRecent <- openIdx[(length(openIdx) - ctlList$opMod$allocYears+1):length(openIdx)] 

    commGears <- ctlList$opMod$commGears

    # Estimate mortality from closed ponded SOK fish and add to historical catch array
    histC_spft <-  obj$om$C_spft[1,,,histdx, drop=F]
    histC_spft[,,6,] <- repObj$pondC_pgt[,6,]*0.65

    # caculate catch for tdxRecent years for allocating catch in projections
    histCatch_spf <- apply( X = histC_spft[,,,tdxRecent,drop=F], FUN = sum, MARGIN = c(1,2,3))

    if( all(histCatch_spf[,,commGears] == 0) )
      histCatch_spf[,,commGears] <- 1

    for(s in 1:nS)
      for( p in 1:nP )
        obj$om$alloc_spf[s,p,commGears] <- histCatch_spf[s,p,commGears] / sum(histCatch_spf[s,p,commGears])

  }

  # Loop and fill errors 
  for( s in 1:nS )
    for( p in 1:nP )
    {
      if( !ctlList$ctl$noProcErr )
        obj$errors$omegaR_spt[s,p,] <- rnorm(nT)

      for( f in 1:nF )
        obj$errors$delta_spft[s,p,f,] <- rnorm(nT)

      # Calculate last year with recruitment estimates
      if(repObj$avgRcode_p[p]==1)
        lastIdx <- max(which(repObj$omegaR_pt[p,] != 0) )
      if(repObj$avgRcode_p[p]==0)
        lastIdx <- max(which(repObj$SRdevs_pt[p,] != 0) )

      if(!is.null(ctlList$opMod$lastIdxOverwrite))
        lastIdx <- ctlList$opMod$lastIdxOverwrite

      # Save historical proc errors, but use simulated recruitments after
      # last estimated recruitment

      # Beverton Holt Recruitment
      # if(repObj$avgRcode_p[p]==0)
      if(ctlList$opMod$recType == 'bevHolt')
      {

      
        if(!ctlList$opMod$posteriorSamples & lastIdx > 0 )
        {

          if(s <= nSpec)
          {
            histRdevs <- repObj$SRdevs_pt[p,1:lastIdx] * obj$om$sigmaR_sp[s,p]
            obj$om$sigmaR_sp[s,p] <- sd( histRdevs )
            # Standardize
            histRdevs <- histRdevs/sd(histRdevs)  
          }

          if( s > nSpec)
          {
            mIdx <- s - nSpec
            histRdevs <- MICEpars[[mIdx]]$resid_t
            obj$om$sigmaR_sp[s,p] <- sd( histRdevs )
            # pad with 0s and standardize
            histRdevs <- c(rep(0,15),histRdevs)/sd(histRdevs)
          }
          
          obj$errors$omegaR_spt[s,p,histdx[1:lastIdx]] <- histRdevs[1:lastIdx]
        }

        if(ctlList$opMod$posteriorSamples)
        {

          # Calculate expected beverton holt recruitments without error
          SB_t    <- mcmcPar$SB_ipt[postDrawIdx,p,histdx]
          R_t     <- mcmcPar$R_ipt[postDrawIdx,p,histdx]
          reca    <- repObj$reca_p[p]
          recb    <- repObj$recb_p[p]
          sigmaR  <- obj$om$sigmaR_sp[s,p]

          bhR_t <- rep(0, length(R_t))

          # Deterimine first year of recruitment generated from bevholt
          firstdx <- min(which(SB_t>0))+ 1

          for(t in firstdx:(tMP-1))
            bhR_t[t]  <- recb*SB_t[t-1]*1000 / (1 + recb*SB_t[t-1])

          # calculate recruitment deviations
          SRdevs_t <- (log(R_t) - log(bhR_t))/sigmaR
          
          if( lastIdx > 0)
            obj$errors$omegaR_spt[s,p,histdx[firstdx:lastIdx]]   <- SRdevs_t[firstdx:lastIdx]
        }          
      }  

      # Average R recruitment
      # if(repObj$avgRcode_p[p]==1)
      if(ctlList$opMod$recType == 'avgR')
      {  
      
        if(!ctlList$opMod$posteriorSamples & lastIdx > 0)
          obj$errors$omegaR_spt[s,p,histdx[1:lastIdx]]   <- repObj$omegaR_pt[p,1:lastIdx]

        if(ctlList$opMod$posteriorSamples)
        {  
          R_t     <- mcmcPar$R_ipt[postDrawIdx,p,histdx]
          Rbar    <- mcmcPar$Rbar_ip[postDrawIdx,p]

          # Deterimine first year of recruitment generated from bevholt
          firstdx <- min(which(R_t>0))
          
          # calculate recruitment deviations
          recDevs_t <- (log(R_t) - log(Rbar))
          if(lastIdx > 0)
            obj$errors$omegaR_spt[s,p,histdx[firstdx:lastIdx]]   <- recDevs_t[firstdx:lastIdx]
        }  

      }
        

      if( !ctlList$ctl$noProcErr & lastIdx > 0 )
        obj$errors$omegaR_spt[s,p,histdx[1:lastIdx]] <- obj$errors$omegaR_spt[s,p,histdx[1:lastIdx]] + 0.5*obj$om$sigmaR_sp[s,p]  # rec devs    


      # Start by generating uncorrelated errors
      obj$errors$ageObsErr_axspft[1:A_s[s],,s,p,,tMP:nT] <- rnorm(A_s[s] * nX * nF * (nT - tMP + 1))
      # Then loop and apply correlation mtx for resids
      Corr_gaa <- repObj$Corr_gaa

      if(s <= nSpec)
        for( fIdx in 1:repObj$nG )
        {
          
          cholesky_aa <- chol(as.matrix(Corr_gaa[fIdx,1:A_s[s],1:A_s[s]]))

          # Check correlation
          if(all(cholesky_aa %in% c(0,1)))
            next
          
          # Loop and correlate
          for( t in tMP:nT )
            for( x in 1:nX )
            {
              errVec <- array(obj$errors$ageObsErr_axspft[1:A_s[s],x,s,p,f,t],dim = c(A_s[s],1))
              obj$errors$ageObsErr_axspft[1:A_s[s],x,s,p,f,t] <- t(cholesky_aa) %*% errVec
            }
        }

    }




  # Replace NaNs with 0s when the recent history has no catch
  obj$om$alloc_spf[is.nan(obj$om$alloc_spf)] <- 0

  # Save historical errors
  obj$errors$delta_spft[1,1:nP,histF,histdx]   <- repObj$z_pgt[1:nP,,histdx,drop = FALSE] # obs errors
  # Get initialization errors
  for( p in 1:nP)
    obj$errors$omegaRinit_asp[1:A_s[1],1,p]  <- repObj$fDevs_ap[1:A_s[1],p] * repObj$sigmaR # Initialisation errors

  if(!is.null(ctlList$opMod$initDevsMult))
    obj$errors$omegaRinit_asp <- obj$errors$omegaRinit_asp * ctlList$opMod$initDevsMult
  
  obj$errors$obsErrMult_spft         <- array(1, dim = c(nS,nP,nF,nT))


  # Random Walk with a trend for M in projections  
  if(ctlList$opMod$projM == 'ranWalkTrend')
  {
      
      # Create pulse M object
      obj$om$pulseM_axspt <- obj$om$M_axspt

      # Natural mortality: lag-1 autocorrelated, log-normal process errors,
      #                    No linear trend, or 50% higher during pulse yrs
      M_pt        <- array(NA, dim=c(nP,nT))
      pulseM_pt   <- array(NA, dim=c(nP,nT))

      for(p in 1:nP)
        M_pt[p,] <- obj$om$M_axspt[2,1,1,p,]

      endM_p      <- apply(M_pt[,1:tMP-1, drop=F], FUN=mean, MARGIN=1)

      # Random walk M scaled to mean==1
      sigmaM    <- ctlList$opMod$sigmaM           # Natural mortality rate CV.
      gammaM    <- ctlList$opMod$gammaM           # Lag-1 autocorr in M
      # Split Mt series among areas
      
      if( ctlList$opMod$projMtype == "identical")
        deltaM_pt <- matrix( rnorm(pT), nrow = nP, ncol = pT, byrow = T)

      if( ctlList$opMod$projMtype == "independent")
        deltaM_pt <- matrix( rnorm(pT*nP), nrow = nP, ncol = pT, byrow = T)      

      if( ctlList$opMod$projMtype == "correlated")
      {
        # Add in correlated M routine here...
        deltaM_pt <- matrix( rnorm(pT*nP), nrow = nP, ncol = pT, byrow = T)      

        # First, we need to get the deviations from the history
        omegaM_pt <- ctlList$opMod$histRpt$omegaM_pt

        # Then calculate a correlation matrix
        corM_p <- cor(t(omegaM_pt))

        # Then generate correlated deviations by left-multiplying
        # delta - need to do some math to make sure these
        # create the right correlated RW devs after fillRanWalk is applied
        cholM_p <- chol(corM_p)
        deltaM_pt <- t(cholM_p) %*% deltaM_pt
      }

      ranM_pt <- matrix(0, nrow = nP, ncol = pT)
      
      for( p in 1:nP )
      {
        ranM_pt[p,]      <- .fillRanWalk( gamma=gammaM, sigma=sigmaM, deltat=deltaM_pt[p,])
        ranM_pt[p,]      <- ranM_pt[p,]/mean( ranM_pt[p,] ) 
      }
      

      # number of years for trend
      if(!is.null(ctlList$opMod$projMtrendT))
        trendT <- ctlList$opMod$projMtrendT
      else
        trendT <- pT

      for(p in 1:nP)
      {  
        # No Trend in M
        # M_pt[p,tMP:nT]  <- M_pt[p,tMP-1] * ranM_pt

        # Trend M for trendT years
        trendM    <- (log(endM_p[p]) - log(M_pt[p,tMP-1]))/ trendT
        M_pt[p,tMP:(tMP+trendT-1)]  <- M_pt[p,tMP-1]*exp( trendM*c(1:trendT) ) * ranM_pt[p,1:trendT]

        # random walks from remaining projection years (i.e., pT-trendT)
        if(trendT < pT)
          M_pt[p,(tMP+trendT):nT] <- endM_p[p]* ranM_pt[p,(trendT+1):pT]

        # PulseM
        if ( ctlList$opMod$pulseMFrq != 0 ) 
          pPulse   <- 1/ctlList$opMod$pulseMFrq                          
        else pPulse <- 0
        nPulse   <- rbinom( n=1, size=pT, prob=pPulse ) # number of events to sample
        pulseYrs <- sample( x=tMP:nT, size=nPulse, replace=FALSE) # sampled years 
        pulseM_pt[p,] <- M_pt[p,]
        pulseM_pt[p,pulseYrs] <- pulseM_pt[p,pulseYrs] * ctlList$opMod$pulseMSize 

        # Output a base version and a pulse version
        for (a in 2:nA)
        {
      
          obj$om$M_axspt[a,1,1,p,tMP:nT]      <- M_pt[p,tMP:nT]
          obj$om$pulseM_axspt[a,1,1,p,tMP:nT] <- pulseM_pt[p,tMP:nT]

        }   
      
      } # END stock area looop 

      # Keep Juvenile M constant
        for( t in tMP:nT )
          obj$om$M_axspt[1,1,1,,t] <- obj$om$M_axspt[1,1,1,,tMP-1]

  }  


  # Otherwise copy M_axsp forward from last year to get running
  if(is.null(ctlList$opMod$projM) | ctlList$opMod$projM == 'densityDepM')  
    for( t in tMP:nT )
      obj$om$M_axspt[,,,,t] <- obj$om$M_axspt[,,,,tMP-1]

  

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

  message(" (.condMS3pop_SISCA) Running OM for historical period.\n")

  obj <- .calcTimes( obj )

  obj$mp$hcr$TAC_spft <- obj$mp$hcr$TAC_spft
  obj$mp$hcr$TAC_spt <- obj$mp$hcr$TAC_spt
  obj$om$F_spft[,,,1:(tMP-1)] <- obj$om$F_spft[,,,1:(tMP-1)] 
  obj$om$C_spft[,,,1:(tMP-1)] <- obj$om$C_spft[,,,1:(tMP-1)] 

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


} # END condMS3pop_SISCA()

# .quantileStratSample()
# Posterior sampling system for 2 leading parameters
.quantileStratSample <- function( seed = NULL,
                                  post = SOGtvMpost, par1 = "m", par2 = "sbo",
                                  nBreaks = 10 )
{
  if( !is.null(seed) )
    set.seed(seed)

  browser(cat('line 4118....'))

  pctiles <- seq(1/nBreaks,1,length = nBreaks)

  par1 <- post[par1]


  par1Quant <- quantile( x = post[par1], probs = pctiles )

  par1Breaks <- c(0,par1Quant)

  samples <- numeric(length = length(pctiles)^2 )

  for( k in 1:nBreaks )
    {
      LB <- par1Breaks[k]
      UB <- par1Breaks[k+1]

      rowIdx1 <- which(post[,par1] > LB & post[,par1] <= UB  )
      par2CondQuants <- quantile( x = post[rowIdx1,par2], probs = pctiles )
      
      par2CondBreaks <- c( 0, par2CondQuants )

      for( j in 1:nBreaks )
      {
        condLB <- par2CondBreaks[j]
        condUB <- par2CondBreaks[j+1]

        rowIdx2 <- which( post[,par2] > condLB & post[,par2] <= condUB )

        rowIdx <- intersect(rowIdx1, rowIdx2)
        if(length(rowIdx) == 0) browser()
        samples[j + (k-1)*nBreaks ] <- sample( x = rowIdx, size = 1 )
      }

    }  

  samples
}


#------------------------------------------------------------------------------#
#-- Natural Mortality functions.                                                    #
#------------------------------------------------------------------------------#

# .fillRanWalk       
# Purpose:        creates a lag-1 autocorrelated error sequence
# Parameters:     gamma=lag-1 autocorrel.; sigma=std error,
#                 deltat=vector of std normals
# Returns:        nT-vector of log-normal error multipliers
# Usage:          tmp[t] = state[t]*error[t]
# Source:         S.P. Cox
.fillRanWalk <- function( gamma=0.5,sigma=0,deltat,init=NULL )
{
  # Generate autocorrelated deviation vector.
  omegat <- vector( mode="numeric",length=length(deltat) )

  # Uncorrelated error std dev and variance
  sigUnc  <- sigma / sqrt( 1.0 - gamma*gamma )
  sig2Unc <- sigma*sigma / ( 1.0 - gamma*gamma )
  
  # Initialize first value in sequence
  if( is.null(init) ) # if not supplied, scale 1st
    omegat[1] <- sigUnc * deltat[1]
  else            # use init
    omegat[1] <- init

  # The random walk...
  for ( t in 2:length(deltat) )
    omegat[t] <- omegat[t-1] * gamma + sigma*deltat[t]*(1 - gamma)
  
  # note these are now multpliers > 0
  omegat <- exp( omegat - sig2Unc/2.0)
  return(omegat)
}     # END function .fillRanWalk

#------------------------------------------------------------------------------#
#-- SISCA helper functions.                                                    #
#------------------------------------------------------------------------------#

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

# ar1Model()
# First-order autoregression model on log-scale, returns
# SD of residuals and lag-1 autocorrelation
ar1Model <- function(dat)
{
  # create dataset
  n  <- length(dat)
  x <- dat[1:(n-1)]
  y <- dat[2:n]

  arMod <- lm(log(y) ~ log(x))

  residSD <- sd(resid(arMod))
  ar1 <- coefficients(arMod)[2]

  arModList <- list(ar1     = ar1,
                    residSD = residSD)

  return(arModList)

}


#------------------------------------------------------------------------------#
#-- Catch, fishing mortality and effort functions.                             #
#------------------------------------------------------------------------------#

# .solveHistPop()
# Solves for the fishing mortality or effort rates
# to achieve an initial depletion target. Optimises
# based on a given objective function.
# inputs: obj = simulation object
# outputs: obj = solved simulation object with all time 
#                 steps filled
# usage:  under omniscient manager sims, or for no fishing.
# NOTE: need some objectives input via simCtlFile
.solveHistPop <- function( obj )
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



  # Now we have enough info to calculate reference points
  stime <- Sys.time()
  message(" (.condMS3pop) Calculating Fmsy and Emsy reference points\n")
  
  repObj <- ctlList$opMod$histRpt
  repObj$om <- obj$om
  repObj$condModel <- ctlList$opMod$condModel
  refPtList <- calcRefPts( repObj )
  obj$rp    <- refPtList

  # Set tInit_p to 0 (0 index in TMB models)
  obj$om$tInit_p <- rep(0,nP)
  # Calculate times for surveys etc.
  obj <- .calcTimes( obj )


  # Set up process errors
  for( s in 1:nS )
    for( p in 1:nP )
      obj$errors$omegaR_spt[s,p,] <- rnorm( n = nT, 
                                            mean = 0, 
                                            sd = 1 )

  # getObjFunctionVal
  # Purpose:        Private function to run operating model and extract objective function
  # Parameters:     pars=log-Fs for t=tMP,...nT
  # Returns:        Objective function f=barrier penalty for biomass less than Bmsy plus
  #                 log(cumulative catch) and penalty for deviation from Fmsy. 
  #                 Multiplier 1000 is for desired precision
  # Source:         S.P. Cox
  getObjFunctionVal <- function( pars, obj )
  {
    # Function to run projections 
    obj <- .runModel( pars, obj, initT = 1, endT = tMP-1 )
    # Return objective function
    finalDep <- obj$om$SB_spt[,,tMP-1]/obj$om$B0_sp

    # browser()

    targDep     <- obj$ctlList$opMod$initDep
    depScalar   <- obj$ctlList$opMod$initOptScalar

    # Penalise large Fs

    val <- sum(log(finalDep/targDep)^2)*depScalar + 0.5*sum(pars^2)/4


    if(val < 1e-3)
      cat("final Dep = ", finalDep, "\n")

    val
  }     # END function getObjFunctionVal


  # -------- Optimisation ------- #
  # Create vector of initial parameters - need to determine
  # if it's nP * nKnots or nS * nP * nKnots
  nPars <- nP * (tMP-1)
  if( ctlList$opMod$effortMod == "targeting" )
    nPars <- nS * nPars

  # Might need to refine this later, but for now we
  # use a random vector
  initPars <- rnorm(n = nPars, sd = 1)

  # OK, run optimisation
  message( " (.solveHistPop) Optimising historical fishing mortality to meet target depletion, please wait.\n")

  initObjFunVal <- round(getObjFunctionVal( pars = initPars, obj = obj), 3)

  message(" (.solveHistPop) Initial objective function value f =", initObjFunVal, ".\n")


  # Find pars that give the best catch/AAV and minimise
  # closures/low catch
  opt <- optim( par = initPars, fn = getObjFunctionVal,
                method = "BFGS", obj = obj,
                control=list(maxit=10000, reltol=0.01 ) )
                # control=list(maxit=3000, reltol=0.001,ndeps=c(.01,.01) ) )

  opt <- optim( par = opt$par, fn = getObjFunctionVal,
                method = "Nelder-Mead", obj = obj,
                control=list(maxit=3000, reltol=0.001,ndeps=c(.01,.01) ) )

  message( " (.solveProjPop) Optimisation of model history completed with f = ", 
            round(opt$value,2), ".\n")
  optPars <- opt$par


  # if( ctlList$ctl$perfConF )
  # {
  #   # Create a dummy optPars object,
  #   # and replace all future Fs with Fmsy in
  #   # runModel
  #   optPars <- initPars 
  # }



  # rerun OM with optimised parameters
  message( " (.solveProjPop) Running OM for historical period with optimised fihsing effort.\n")
  obj <- .runModel( pars = optPars, obj = obj, 
                    initT = 1, endT = tMP-1 )

  return(obj)
} # END .solveProjPop()




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


  # getObjFunctionVal
  # Purpose:        Private function to run operating model and extract objective function
  # Parameters:     pars=log-Fs for t=tMP,...nT
  # Returns:        Objective function f=barrier penalty for biomass less than Bmsy plus
  #                 log(cumulative catch) and penalty for deviation from Fmsy. 
  #                 Multiplier 1000 is for desired precision
  # Source:         S.P. Cox
  getObjFunctionVal <- function( pars, obj )
  {
    # Function to run projections 
    obj <- .runModel( pars, obj, initT = tMP, endT = nT )
    # Return objective function
    val <- calcOmniObjFun(obj)$objFun
    
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
  obj <- .runModel(  pars = optPars, obj = obj, 
                    initT = tMP, endT = nT )

  return(obj)
} # END .solveProjPop()


# .runModel()
# Runs operating model with an input of fishing 
# mortality/effort rates and returns resulting 
# simulation object
# inputs: pars = vector of F or E values
#         obj = simulation object
# outputs: outList = list of  
#                      - obj = simulation object
#                      - f = objective function value
# Modified from SP Cox's function
.runModel <- function( pars, obj, initT = tMP, endT = nT  )
{
  # Pull control list items
  om      <- obj$om
  ctlList <- obj$ctlList
  mp      <- ctlList$mp

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
  nKnots      <- mp$omni$nKnots

  if(ctlList$opMod$histType=="sim" | endT < om$tMP)
  {
    parMult <- 1
    nKnots  <- om$tMP - 1
  }
  

  # Model dimensions
  nS      <- om$nS
  nP      <- om$nP
  nT      <- om$nT
  tMP     <- om$tMP
  pT      <- ctlList$opMod$pT


  # Some pars that are global to either effort
  # simulation
  projInt <- endT - initT + 1
  
  space   <- projInt / (nKnots - 1)
  knotPts <- round( seq(from = 1, by = space, length.out = nKnots) )
  knotPts[knotPts > projInt] <- projInt

  # Create x points for interpolation spline
  x <- knotPts

  # Set fishing mortality to 0 for tMP:nT
  # so that operating model does not hang
  obj$om$F_spft[,,,initT:endT] <- 0

  # Choose whether this is effort or fishing mort.
  # If perfect targeting (all TACs fully utilised)
  if( ctlList$opMod$effortMod == "targeting" )
  {
    # This is a projection model that always uses
    # a multiple of Fmsy
    if( ctlList$ctl$perfConF & endT >= tMP )
    {
      if( !is.null(ctlList$mp$omni$inputF) )
        obj$om$F_spft[,,2,initT:endT] <- ctlList$mp$omni$inputF

      
      if(is.null(ctlList$mp$omni$inputF))
        for( s in 1:nS )
          for( p in 1:nP )
            obj$om$F_spft[s,p,2,initT:endT] <- parMult*obj$rp$FmsyRefPts$Fmsy_sp[s,p]

      
    } else {  
      # Then we are solving for stock/species specific
      # Fs (nS * nP at each time step)
      tmpFmult <- exp( pars )
      knotF_spk <- array( tmpFmult, dim = c(nS,nP,nKnots) )

      # Make spline for each species/stock
      for( s in 1:nS )
        for( p in 1:nP )
        {
          if( nKnots < (endT - initT + 1) )
          {
            splineF <- spline( x = x, y = knotF_spk[s,p,], n = projInt)$y
            splineF <- splineF * obj$rp$FmsyRefPts$Fmsy_sp[s,p]
            
          } else splineF <- knotF_spk[s,p,] * obj$rp$FmsyRefPts$Fmsy_sp[s,p]

          splineF <- splineF * parMult

          splineF[splineF < 0] <- 0
          splineF[splineF > maxF] <- maxF

          obj$om$F_spft[s,p,2,initT:endT] <- splineF
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
        obj$om$E_pft[p,2,initT:endT] <- parMult*obj$rp$EmsyMSRefPts$EmsyMS_p[p]

      if( !is.null(ctlList$omni$inputE) )
        for( p in 1:nP )
          obj$om$E_pft[p,2,initT:endT] <- ctlList$omni$inputE_p[p]

      for( s in 1:nS )
        for( p in 1:nP )
          obj$om$F_spft[s,p,2,initT:endT] <- obj$om$E_pft[p,2,initT:endT] * om$qF_spft[s,p,2,tMP - 1]

    } else {
      # Then we are solving for area specific
      # Es (nP at each time step)
      tmpEmult    <- exp( pars )
      knotE_pk    <- array( tmpEmult, dim = c(nP,nKnots) )

      # Make spline for each species/stock
      for( p in 1:nP )
      {
        # Make spline, multiply by Emsy
        splineE <- spline( x = x, y = knotE_pk[p,], n = projInt)$y
        splineE <- splineE * obj$rp$EmsyMSRefPts$EmsyMS_p[p]

        # Correct for under and over efforts, limit
        # lowest effort to be the minimum over the
        # historical period
        if( minE == "hist" )
          splineE[splineE < 0] <- min(obj$om$E_pft[p,2,1:(tMP-1)])/obj$rp$EmsyMSRefPts$EmsyMS_p[p]
        else 
          splineE[splineE < 0] <- minE
        
        splineE[splineE > maxE] <- maxE

        # Assign to projection
        obj$om$E_pft[p,2,initT:endT] <- parMult * splineE

        for( s in 1:nS )
          obj$om$F_spft[s,p,2,initT:endT] <- obj$om$E_pft[p,2,initT:endT] * om$qF_spft[s,p,2,tMP - 1]
      }
    }

  }

  

  

  # run model for projection period
  for(t in initT:endT )
  {
    obj <- .ageSexOpMod(obj, t)
  }

  return(obj)
} # END .runModel()

calcOmniObjFun <- function( obj )
{
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


  # Model dimensions
  nS      <- om$nS
  nP      <- om$nP
  nT      <- om$nT
  tMP     <- om$tMP
  pT      <- ctlList$opMod$pT


  # Some pars that are global to either effort
  # simulation
  projInt <- nT - tMP + 1
  nKnots  <- mp$omni$nKnots
  space   <- projInt / (nKnots - 1)
  knotPts <- round( seq(from = 1, by = space, length.out = nKnots) )
  knotPts[knotPts > projInt] <- projInt


  # Pull projection catch and biomass
  Cproj_spft    <- obj$om$C_spft[,,2,tMP:nT,drop = FALSE]
  TACproj_spft  <- obj$om$TAC_spft[,,2,tMP:nT,drop = FALSE]
  Bproj_spt     <- obj$om$B_spt[,,tMP:nT,drop = FALSE]
  Bmsy_sp       <- obj$rp$FmsyRefPts$BeqFmsy_sp
  

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
  histCatch_spt <- apply( X = om$C_spft[,,,1:(tMP-1),drop = FALSE], FUN = sum, MARGIN = c(1,2,4))
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


  catDiff_spt <- array(0, dim = c(nS,nP,pT))
  catDiffRel_spt <- catDiff_spt
  for( s in 1:nS )
    for( p in 1:nP )
    {
      catDiff_spt[s,p,] <- abs(diff(obj$om$C_spft[s,p,2,(tMP-1):nT]))
      catDiffRel_spt[s,p,] <- catDiff_spt[s,p,] / (obj$om$C_spft[s,p,2,(tMP:nT) - 1])
    }
  
  catDiffRel_spt[!is.finite(catDiffRel_spt)] <- 1

  initCatDiffRel_sp <- array(0, dim = c(nS,nP))
  initCatDiffRel_sp[1:nS,1:nP] <- catDiffRel_spt[,,1]

  # Effort difference (prefer to keep total effort similar)
  totEff_t   <- apply( X = obj$om$E_pft[,2,(tMP-1):nT,drop = FALSE],
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
      closedCount_sp[s,p] <- sum( Cproj_spft[s,p,1,] < Cmin_sp[s,p] )

      # Penalty on AAV
      AAV <- .calcAAV( Cproj_spft[s,p,1,] )
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
  Csum    <- sum(Cproj_spft,na.rm = T)
  Cbar_sp <- apply( X = Cproj_spft, FUN = mean, MARGIN = c(1,2))
  totCbar <- mean( apply(X = Cproj_spft, FUN = sum, MARGIN = 4 ) )

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

                  
 
  objFun    <-  sum(objFun_sp) -
                totCatWt * log(1e3*totCbar) -
                sumCatWt * log(1e3*Csum) 


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

  # Return results
  
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
} # END calcOmniObjFun()


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

# applyDiscreteFisheries()
# An alternative removals model to the 
# Baranov catch equation. Steps through fisheries
# in order
applyDiscreteFisheries <- function( N_axsp,
                                    sel_axspf,
                                    W_axspf,
                                    fleetTiming_f,
                                    fleetType_f,
                                    spawnTiming_s,
                                    M_axsp,
                                    TAC_spf,
                                    P_spf,
                                    pondM_f,
                                    nA, nX, nS, nP, nF )
{
  # First, order fleetTiming
  fleetOrder <- order(fleetTiming_f)


  tmpN_axsp <- array(0, dim = c(nA,nX,nS,nP) )
  s_axspf   <- array(0, dim = c(nA,nX,nS,nP,nF) )
  C_spf     <- array(0, dim = c(nS,nP,nF))


  # Arrays to hold states at different times
  endN_axsp     <- array( 0, dim = dim(tmpN_axsp) )
  vN_axspf      <- array( 0, dim = dim(s_axspf) )
  vB_axspf      <- array( 0, dim = dim(s_axspf) )
  pvB_axspf     <- array( 0, dim = dim(s_axspf) )
  Cw_axspf      <- array( 0, dim = dim(s_axspf) )
  C_axspf       <- array( 0, dim = dim(s_axspf) )
  spawnN_axsp   <- array( 0, dim = dim(tmpN_axsp) )
  F_spf         <- array( 0, dim = dim(C_spf) )

  # Create a temporary N_axsp to hold successive
  # numbers at each fractional time step
  tmpN_axsp[1:nA,,,]    <- N_axsp
  s_axspf[1:nA,,,,]     <- sel_axspf
  C_spf[1:nS,,]         <- TAC_spf

  # Overwrite TAC with ponded fish for SOK fleets
  C_spf[1:nS,,fleetType_f %in% 2:3] <- P_spf[,,fleetType_f %in% 2:3]

  # if(debugFlag)
  #   browser()
  # Loop over species to avoid spawn timing issues
  for( s in 1:nS )
  {
    # Initialise lastFleetTime at 0
    lastFleetTime <- 0

    # First, check if spawn timing goes first
    if( spawnTiming_s[s] <= fleetTiming_f[fleetOrder[1]] )
    {
      spawnN_axsp[,,s,] <- calcDiscSpawnN(  lastFleetTime = lastFleetTime,
                                            spawnTiming = spawnTiming_s[s],
                                            currN_axsp = tmpN_axsp[,,s,,drop = FALSE],
                                            M_axsp = M_axsp[,,s,,drop = FALSE] )

      lastFleetTime <- spawnTiming_s[s]
      tmpN_axsp[,,s,]   <- spawnN_axsp[,,s,]
    }

    # if(all(C_spf == 0))
    #     browser()

    # Then loop over fleets, check on spawn timing
    # in each loop
    for( oIdx in 1:length(fleetOrder) )
    {
      fIdx <- fleetOrder[oIdx]
      nextFleetTime <- fleetTiming_f[fIdx]

      # if(debugFlag)
      #   browser()

      # Calculate spawning numbers
      if( spawnTiming_s[s] > lastFleetTime & spawnTiming_s[s] <= nextFleetTime )
      {
        spawnN_axsp[,,s,] <- calcDiscSpawnN(  lastFleetTime = lastFleetTime,
                                              spawnTiming = spawnTiming_s[s],
                                              currN_axsp = tmpN_axsp[,,s,],
                                              M_axsp = M_axsp[,,s,] )

        lastFleetTime <- spawnTiming_s[s]
        tmpN_axsp[,,s,]   <- spawnN_axsp[,,s,]
      }

      # Now calculate fishery numbers
      fracM <- nextFleetTime - lastFleetTime
      # Deplete by natural mortality
      tmpN_axsp[,,s,] <- tmpN_axsp[,,s,] * exp(-fracM * M_axsp[,,s,])

      # Calculate vulnerable numbers at age
      vN_axspf[,,s,,fIdx] <- tmpN_axsp[1:nA,1:nX,s,1:nP] * s_axspf[1:nA,1:nX,s,1:nP,fIdx]

      # Calculate vulnerable biomass at age
      vB_axspf[,,s,,fIdx] <- vN_axspf[,,s,,fIdx] * W_axspf[,,s,,fIdx]

      if(fIdx > nF | any(!is.finite(C_spf[s,,])))
      {
        cat("Bad catch \n")
        browser()
      }

      if(is.null(C_spf[s,,]))
        browser()

      
      if( any(C_spf[s,,fIdx] > 0 ) )
      { 
        # Compute proportion of biomass at age
        # for catch distribution across age classes
        for( p in 1:nP)
        {
          pvB_axspf[,,s,p,fIdx] <- vB_axspf[,,s,p,fIdx] / (sum(vB_axspf[,,s,p,fIdx],na.rm = T) + 1e-9)
          Cw_axspf[,,s,p,fIdx]  <- C_spf[s,p,fIdx] * pvB_axspf[,,s,p,fIdx]
        }

        # Compute catch-at-age in numbers
        C_axspf[,,s,,fIdx] <- Cw_axspf[,,s,,fIdx] / W_axspf[,,s,,fIdx]
        
        # Remove from general population, update lastFleetTime
        tmpN_axsp[1:nA,,s,] <- tmpN_axsp[1:nA,,s,] - C_axspf[,,s,,fIdx]
      }
      
      lastFleetTime <- fleetTiming_f[fIdx]

      # Try to mimic the posfun here
      eps <- 1e-3

      # Overwrite NAs which stem from staggered initTs
      ltepsIdx <- tmpN_axsp[,,s,] < eps
      ltepsIdx[is.na(ltepsIdx)] <- FALSE

      # Apply posfun
      tmpN_axsp[,,s,][ltepsIdx] <- eps/(2-tmpN_axsp[,,s,][ltepsIdx]/eps)
    }

    if( spawnTiming_s[s] > lastFleetTime )
    {
      spawnN_axsp[,,s,] <- calcDiscSpawnN(  lastFleetTime = lastFleetTime,
                                            spawnTiming = spawnTiming_s[s],
                                            currN_axsp = tmpN_axsp[,,s,],
                                            M_axsp = M_axsp[,,s,] )

      lastFleetTime <- spawnTiming_s[s]
      tmpN_axsp[,,s,]   <- spawnN_axsp[,,s,] 
    }

    # Then compute endN by applying last bit of M
    fracM     <- 1 - lastFleetTime
    endN_axsp[,,s,] <- tmpN_axsp[,,s,] * exp(-fracM*M_axsp[,,s,])

    # Add ponded fish at age reduced by post-ponding M
    sokFleets <- which(fleetType_f %in% c(2,3))
    for( fIdx in sokFleets )
    {  
      
      # calculate remaining natural mortality and apply to 
      fracM     <- 1 - fleetTiming_f[fIdx]

      endN_axsp[1:nA,,s,] <- endN_axsp[1:nA,,s,] + C_axspf[1:nA,,s,,fIdx] * exp(-fracM*M_axsp[1:nA,,s,] - pondM_f[fIdx] )

      #Apply post-ponding mortality so catch for SOK fleets represents dead fish
      Cw_axspf[1:nA,,s,,fIdx] <- Cw_axspf[1:nA,,s,,fIdx] * (1-exp(-pondM_f[fIdx]))
      C_axspf[1:nA,,s,,fIdx]  <- C_axspf[1:nA,,s,,fIdx] * (1-exp(-pondM_f[fIdx]))

    }  

  }
  # if(debugFlag)
  #   browser()

  if(sum(endN_axsp,na.rm = T) == 0 )
  {
    cat("Error in discrete fishery calc.\n")
    browser()
  }

  # Return model states
  outList <- list(  endN_axsp   = endN_axsp,
                    vN_axspf    = vN_axspf,
                    vB_axspf    = vB_axspf,
                    Cw_axspf    = Cw_axspf,
                    C_axspf     = C_axspf,
                    spawnN_axsp = spawnN_axsp )

  outList

} # END applyDiscreteFisheries

# calcDiscSpawnN()
# Calculates numbers at age at spawn timing
# under the discrete fisheries model
calcDiscSpawnN <- function(   lastFleetTime,
                              spawnTiming,
                              currN_axsp,
                              M_axsp )
{
  fracM <- spawnTiming - lastFleetTime
  spawnN_axsp <- currN_axsp * exp(-fracM * M_axsp)
  
  return(spawnN_axsp)

} # END calcDiscSpawnN()

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
} # END .calcHistEffortDynamics()

# calcSOKpsi()
# Converts landed SOK product to ponded fish for
# accurate fishery removals. 
calcSOKpsi <- function( N_axsp,
                        sel_axsp,
                        mat_asp,
                        pFem = .5,
                        fec = 200,
                        initF = 0.01,
                        pEff = .5,
                        gamma = 0.0024,
                        W_axsp,
                        A_s,
                        nA,nX,nS,nP  )
{
  # First, we need to estimate proportion
  # mature

  # Total mortality
  Z_axsp        <- array(0, dim = dim(sel_axsp))
  pondC_axsp    <- array(0, dim = dim(sel_axsp))

  tmpBmat_axsp  <- array(0, dim = dim(sel_axsp))
  tmpB_axsp     <- array(0, dim = dim(sel_axsp))

  # tmp ponded fish for pMat calc
  Z_axsp[,,,1:nP] <- sel_axsp[,,,1:nP] * initF
  pondC_axsp[,,,1:nP] <- (1 - exp(-Z_axsp[,,,1:nP])) * N_axsp[,,,1:nP] * sel_axsp[,,,1:nP] * initF/Z_axsp[,,,1:nP]


  for( s in 1:nS)
    for( a in 1:A_s[s] )
      for( x in 1:nX)
      {
        tmpBmat_axsp[a,,s,]  <- pondC_axsp[a,,s,] * mat_asp[a,s,] * W_axsp[a,,s,]
        tmpB_axsp[a,,s,]     <- pondC_axsp[a,,s,] * W_axsp[a,,s,]
      }

  tmpBmat_sp <- apply( X = tmpBmat_axsp[,nX,,,drop = FALSE], FUN = sum, MARGIN = c(3,4), na.rm=T)
  tmpB_sp    <- apply( X = tmpB_axsp[,,,,drop =FALSE], FUN = sum, MARGIN = c(3,4), na.rm=T)

  # Finally, proportion mature
  propMat_sp <- tmpBmat_sp / tmpB_sp

  # Then we calculate psi
  psi_sp <- pEff * pFem * gamma * fec * propMat_sp;

  outList <- list(  propMat_sp = propMat_sp,
                    psi_sp     = psi_sp )

  return(outList)
} # END calcSOKpsi()

# .solveBaranov()
# Generalised Newton-Rhapson solver or the Baranov
# catch equation, given catch, vulnerable biomass
# and other stuff.
.solveBaranov <- function(  C_spf       = newCat_spf, 
                            vB_axspf    = vB_axspft[,,,,,t],
                            vB_spf      = vB_spft[,,,t],
                            N_axsp      = N_axspt[,,,,t],
                            sel_axspf   = sel_axspft[,,,,,t],
                            M_axsp      = M_axspt[,,,,t],
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
        appZ_axsp[1:A_s[s],x,s,p] <- M_axsp[1:A_s[s],x,s,p]

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

          appC_spf[s,p,] <- apply( X = appC_axspf[,,s,p,,drop = FALSE], FUN = sum, MARGIN = c(5))

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
            appZ_axsp[1:A_s[s],x,s,p] <- M_axsp[1:A_s[s],x,s,p]

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
  # First, we need to get whichFleets
  TAC_f <- apply( X = TAC_spf, FUN = sum, MARGIN = 3 )
  whichFleets <- which(TAC_f > 0)



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

        # Solve for the effort using the TAC and 
        solveE_s <- -1/qF_spf[,p,f] * log( 1 - catRem_s / (vB_spf[,p,f] - C_spf[,p,f]))
        propE <- propE + 0.8*min(solveE_s,na.rm = T)

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

        if( any(catRem_s < 0 ) )
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






#------------------------------------------------------------------------------#
#-- Penalty functions                                                        --#
#------------------------------------------------------------------------------#


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

