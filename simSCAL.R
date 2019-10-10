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

  # Initialise the data structures
  # to hold simulation states etc.
  simObj <- .initMS3pop(  ctlList )

  # Run mgmtProc
  blob <- .mgmtProc( simObj )

  # Save output
  .saveBlob(blob)

  beepr::beep("complete")
}


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
  tMP <  om$tMP




  # 2. Perform stock assessment
  stockAssessment <- .callProcedureAM( obj, t )


  # 3. Determine stock status and set catch limit
  catchLim_sp     <- .applyHCR( obj, stockAssessment, t )

  # Now spread catch among allocation
  for( s in 1:nS )
    for( p in 1:nP )
      obj$om$C_spft[s,p,,t]    <- catchLim_sp[s,p] * alloc_spf[s,p,]

  # 4. Update simulated population
  obj$om <- om
  obj$mp <- mp

  obj <- .ageSexOpMod(obj, t)

  # 5. Return updated sim object
  return(obj)

} # END .updatePop()




# .mgmtProc()
# Runs replications of the closed loop
# simulation, applying the simulated mgmt
# procedure to the operating model stocks
mgmtProc <- function( obj )
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
  om <- list( iSeed_i = rep(NA, nReps),
              SB_ispt  = array( NA, dim = c(nReps, nS, nP, nT) ),      # SSB
              B_ispt   = array( NA, dim = c(nReps, nS, nP, nT )),      # Total biomass
              vB_ispft = array( NA, dim = c(nReps, nS, nP, nF, nT )),  # Vulnerable biomass
              R_ispt   = array( NA, dim = c(nReps, nS, nP, nT) ),      # Rec't
              C_ispft  = array( NA, dim = c(nReps, nS, nP, nF, nT) ),  # Catch (kt)
              F_ispft  = array( NA, dim = c(nReps, nS, nP, nF, nT) ),  # Fishing Mort
              E_ipft   = array( NA, dim = c(nReps, nP, nF, nT) )       # Fishing effort
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
                  UCP_ispt = array(NA, dim = c(nReps, nS, nP, nT) )       )   # Upper control point


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

    diagCondition( repObj = ctlList$repObj, ms3Obj = test )

    browser()


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
    blob$om$SBt_ispt[i,,,]    <- simObj$om$SB_spt
    blob$om$B_ispt[i,,,]      <- simObj$om$B_spt
    blob$om$vB_ispft[i,,,,]   <- simObj$om$vB_spft
    blob$om$R_ispt[i,,,]      <- simObj$om$R_spt
    blob$om$C_ispft[i,,,,]    <- simObj$om$C_spft
    blob$om$F_ispft[i,,,,]    <- simObj$om$F_spft
    # blob$om$E_ipft[i,,,]    ## FILL WHEN EFFORT DYNAMICS COMPLETE ##

    # Errors - maybe update simObj structure to match blob here
    blob$om$errors$omegaR_ispt[i,,,]  <- simObj$errors$omegaR_spt
    blob$om$errors$delta_ispft[i,,,]  <- simObj$errors$delta_spft

    # Data
    blob$mp$data$I_ispft[i,,,,]       <- simObj$mp$data$I_spft
    # Add ages and lengths later

    # HCR elements
    blob$mp$hcr$Fref_ispt[i,,,]       <- simObj$mp$hcr$Fref_spt
    blob$mp$hcr$Bref_ispt[i,,,]       <- simObj$mp$hcr$Bref_spt
    blob$mp$hcr$targetF_ispt[i,,,]    <- simObj$mp$hcr$targetF_spt
    blob$mp$hcr$LCP_ispt[i,,,]        <- simObj$mp$hcr$LCP_spt
    blob$mp$hcr$UCP_ispt[i,,,]        <- simObj$mp$hcr$UCP_spt


    # Retrospective AM outputs
    blob$mp$assess$retroR_spt[i,,,,]    <- simObj$mp$assess$retroR_spt
    blob$mp$assess$retroSB_spt[i,,,,]   <- simObj$mp$assess$retroSB_spt
    blob$mp$assess$retroVB_spft[i,,,,,] <- simObj$mp$assess$retroVB_spft
    blob$mp$assess$assessMod_i[[i]]     <- simObj$mp$assess$assessModOutputs
  

  } # END i loop


  # Name the blob array dimensions...

  return(blob)
} # END .mgmtProc()



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

  # Now, copy model fit
  histdx <- 1:(repObj$nT)
  obj$om$F_spft[,,,histdx]        <- repObj$F_spft
  obj$om$C_spft[,,,histdx]        <- repObj$predCw_spft
  obj$om$sel_axspft[,,,,,histdx]  <- repObj$sel_axspft
  obj$om$sel_lspft[,,,,histdx]    <- aperm(repObj$sel_lfspt,c(1,3,4,2,5))


  # Copy catchability
  obj$om$q_spft[,,,histdx]        <- repObj$q_spft
  obj$om$q_spf                    <- repObj$q_spf


  # Add historical data
  obj$data$I_spft[,,,histdx]      <- repObj$I_spft
  obj$data$A_axspft[,,,,,histdx]  <- aperm( repObj$age_aspftx, c(1,6,2:5) )
  obj$data$L_lxspft[,,,,,histdx]  <- aperm( repObj$len_lspftx, c(1,6,2:5) )


  for( s in 1:nS )
    for( p in 1:nP )
    {
      obj$errors$omegaR_spt[s,p,] <- rnorm(nT)
      for( f in 1:nF )
        obj$errors$delta_spft[s,p,f,] <- rnorm(nT)
    }

  # Save historical errors
  obj$errors$omegaR_spt[,,histdx]   <- repObj$omegaR_spt     # rec devs
  obj$errors$delta_spft[,,,histdx]  <- repObj$residCPUE_spft # obs errors
  obj$errors$omegaRinit_asp         <- repObj$omegaRinit_asp # Initialisation errors

  message("(.condMS3pop) Running OM for historical period.\n")

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
initMS3pop <- function( repObj = reports$repOpt,
                        pT = 10 )
{
  # make OM list
  om <- list()

  message("(.initMS3pop) Initialising MS3 operating model.\n")

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
  om$C_spft     <- array(0,  dim = c(nS,nP,nF,nT) )        # Total catch in kt
  om$C_axspft   <- array(0,  dim = c(nA,nX,nS,nP,nF,nT) )  # Numbers
  om$Cw_axspft  <- array(0,  dim = c(nA,nX,nS,nP,nF,nT) )  # Weight
  om$vB_axspft  <- array(0,  dim = c(nA,nX,nS,nP,nF,nT) )  # vuln Bio (granular)
  om$vB_spft    <- array(0,  dim = c(nS,nP,nF,nT) )        # vuln Bio (aggregate)
  om$sel_axspft <- array(0,  dim = c(nA,nX,nS,nP,nF,nT) )  # sel at age (tv)
  om$sel_lspft  <- array(0,  dim = c(nL,nS,nP,nF,nT) )     # sel at len (tv)
  om$sel_axspf  <- array(0,  dim = c(nA,nX,nS,nP,nF) )     # sel at age
  om$sel_lspf   <- array(0,  dim = c(nL,nS,nP,nF) )        # sel at len
  om$F_spft     <- array(0,  dim = c(nS,nP,nF,nT) )        # Fishing mortality

  # Observation model pars
  om$q_spft     <- array(0,  dim = c(nS,nP,nF,nT) )        # catchability (tv)
  om$q_spf      <- array(0,  dim = c(nS,nP,nF) )           # catchability

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

  mp$data$I_spft    <- array( NA, dim = c(nS,nP,nF,nT) )
  mp$data$A_axspft  <- array( NA, dim = c(nA,nX,nS,nP,nF,nT) ) 
  mp$data$L_lxspft  <- array( NA, dim = c(nL,nX+1,nS,nP,nF,nT) ) 


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
                errors = errors )

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
  C_spft            <- om$C_spft
  
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

  vB_spft <- apply( X = vB_axspft, FUN = sum, MARGIN = 3:6, na.rm = T )
  B_spt   <- apply( X = B_axspt, FUN = sum, MARGIN = 3:5, na.rm = T )

  # Now calculate spawning biomass
  SB_asp <- matAge_asp * B_axspt[,nX,,,t]
  SB_spt[,,t] <- apply(X = SB_asp, FUN = sum, MARGIN = c(2,3), na.rm = T )

  # Now calculate Fs from TACs (for closed loop sim)
  if( t >= tMP )
  {
    # Apply pope's approx to get us moving here...
    F_spft[,,,t] <- mfPA( C_spf     = C_spft[,,,t], 
                          vB_axspf  = vB_axspft[,,,,,t],
                          N_axsp    = N_axspt[,,,,t],
                          M_xsp     = M_xsp,
                          A_s       = A_s,
                          wt_axsp   = meanWtAge_axsp,
                          nS = nS, nP = nP, nX = nX, nF = nF )

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

  # Sum catch
  C_spft[,,,t] <- apply( X = Cw_axspft[,,,,,t], FUN = sum, MARGIN = c(3,4,5), na.rm = T )

  # Now generate indices - need to add a schedule later
  if( t >= tMP )
  {
    I_spft[,,,t] <- q_spft[,,,t] * vB_spft[,,,t] * exp(om$tauObs_spf * err$delta_spft[,,,t])
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

  # Biomass indices
  I_spft            -> obj$mp$data$I_spft

  return( obj )
} # END ageSexOpMod()

# mfPA()
# Multi-fleet Pope's approximation
mfPA <- function( C_spf     = newCat_spf, 
                  vB_axspf  = vB_axspft[,,,,,t],
                  N_axsp    = N_axspt[,,,,t],
                  M_xsp     = M_xsp,
                  A_s       = A_s,
                  wt_axsp   = meanWtAge_axsp,
                  nS = nS, nP = nP, nX = nX, nF = nF  )
{
  # Make arrays to hold info
  pvB_axspf         <- array(0, dim = dim(vB_axspf))
  remN_axsp         <- array(0, dim = dim(N_axspt))
  endN_axsp         <- array(0, dim = dim(N_axspt))
  appZ_axsp         <- array(0, dim = dim(N_axspt))
  catNumAge_axspf   <- array(0, dim = dim(vB_axspf))
  appF_axspf        <- array(0, dim = dim(vB_axspf) )
  appF_spf          <- array(0, dim = dim(C_spf))


  for( sIdx in 1:nS )
    for( pIdx in 1:nP )
    {
      for( fIdx in 1:nF )
      {
        # Pull biomass at age, convert to proportions
        pvB_axspf[,,sIdx,pIdx,fIdx] <- pvB_axspf[,,sIdx,pIdx,fIdx]/sum(vB_axspf[,,sIdx,pIdx,fIdx])

        # Calculate numbers to remove by converting catch to weight
        catNumAge_axspf[,,sIdx,pIdx,fIdx] <- C_spft[sIdx, pIdx, fIdx] * pvB_axspf[,,sIdx,pIdx,fIdx] / wt_axsp[,,sIdx,pIdx]


        # Now add the catch at age in numbers to the removed fish
        remN_axsp[,,sIdx,pIdx] <- remN_axsp[,,sIdx,pIdx] + catNumAge_axspf[,,sIdx,pIdx,fIdx]

      }
      for( xIdx in 1:nX)
        endN_axsp[,xIdx,sIdx,pIdx] <- N_axspt[,xIdx,sIdx,pIdx] * exp(-M_xsp[sIdx,pIdx,xIdx]) - remN_axspt[,xIdx,sIdx,pIdx]* exp(-M_xsp[sIdx,pIdx,xIdx]/2)

      # Now compute Z approximation
      appZ_axspt[1:A_s[sIdx],,sIdx,pIdx] <- log( N_axspt[1:A_s[sIdx],,sIdx,pIdx] / endN_axspt[1:A_s[sIdx],,sIdx,pIdx])
      # Loop over fIdx again
      for( fIdx in 1:nF )
      {
        # browser()
        # This is just checking that the F is the same across age classes
        appF_axspf[1:A_s[sIdx],,sIdx,pIdx,fIdx] <- catNumAge_axspf[1:A_s[sIdx],,sIdx,pIdx,fIdx] * appZ_axspt[1:A_s[sIdx],,sIdx,pIdx]
        appF_axspf[1:A_s[sIdx],,sIdx,pIdx,fIdx] <- appF_axspf[1:A_s[sIdx],,sIdx,pIdx,fIdx]/N_axspt[1:A_s[sIdx],,sIdx,pIdx]/sel_aspfxt[1:A_s[sIdx],sIdx,pIdx,fIdx,]/(1 - exp(-appZ_axspt[1:A_s[sIdx],,sIdx,pIdx]))
        appF_spf[sIdx,pIdx,fIdx] <- appF_axspf[A_s[sIdx],2,sIdx,pIdx,fIdx]
      }

    }

  appF_spf
} # END mfPA()


