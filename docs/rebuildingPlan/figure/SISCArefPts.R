# --------------------------------------------------------------------------
# refPts.R
# 
# Takes parameters from a report object and calculates
# reference points and curves under an age structured/YPR formulation. 
# Called from within runHierSCAL() after runnings an AM.
# 
# 
# Author: Samuel Johnson
# Date: March 7, 2019
# 
# --------------------------------------------------------------------------

# calcRefPts()
# Calculates biological reference points based on given biological parameters 
# assuming a delay difference biological model.
# inputs:   obj = list of biological parameters
# ouputs:   refPts = list() of reference points
calcRefPts <- function( obj )
{
  # Calculate selectivity
  # obj <- .calcSel_p(obj, fleetIdx = 1)

  # Calculate reference curves
  refCurves <- .calcRefCurves( obj )

  # First, let's just do Fmsy reference points
  FmsyRefPts <- .getFmsy_p( obj = obj, 
                            refCurves = refCurves )
  
  obj$refPts <- list()
  obj$refPts$refCurves    <- refCurves
  obj$refPts$FmsyRefPts   <- FmsyRefPts
  obj$refPts$Surv_ap      <- refCurves$Surv_apf[,,1]

  obj
}

# .calcRefCurves()
# Calculates equilibrium curves of equilbrium biomass, numbers,
# yield and recruitment as a function of input fishing mortality rates
# inputs:   obj = list of biological parameters
# ouputs:   refCurves = list() of reference curves (vectors)
.calcRefCurves <- function( obj, nFs = 1000 )
{
  # First, compute max F (tolerance of 1e-5)
  # maxF <- max( 10*obj$meanM_p )
  maxF <- 10

  # We're going to need to fill each species' ref curves,
  # so labeling and dimensions are needed
  nP          <- obj$nP
  nA          <- obj$nA
  stockNames  <- dimnames(obj$M_p)[[2]]


  f <- seq( from = 0.0, to = maxF, length = nFs )

  # Create matrices to hold Recruitment reference curve, name rows and copy
  # for each of Beq, Neq and Yeq
  Req_pf      <- array( NA, dim = c(nP, nFs),
                            dimnames = list(  stock = stockNames,
                                              F = f ) )
  Surv_apf    <- array( NA, dim = c(nA, nP, nFs),
                            dimnames = list(  age = 1:nA,
                                              stock = stockNames,
                                              F = f ) )

  Beq_pf       <- Req_pf
  tBeq_pf      <- Req_pf
  expBeq_pf    <- Req_pf
  Yeq_pf       <- Req_pf
  ypr_pf       <- Req_pf
  ssbpr_pf     <- Req_pf
  Ueq_pf       <- Req_pf

  # Loop and fill
  for( i in 1:length(f) )
  {
    tmp           <- .calcEquil( f = f[i], obj = obj )
    Req_pf[,i]    <- tmp$Req_p
    Beq_pf[,i]    <- tmp$Beq_p
    tBeq_pf[,i]   <- tmp$tBeq_p
    expBeq_pf[,i] <- tmp$expBeq_p
    Yeq_pf[,i]    <- tmp$Yeq_p
    ypr_pf[,i]    <- tmp$ypr_p
    ssbpr_pf[,i]   <- tmp$ssbpr_p
    Ueq_pf[,i]    <- tmp$Ueq_p
    Surv_apf[,,i] <- tmp$Surv_ap
  }

  refCurves <- list()
    refCurves$F           <- f
    refCurves$ypr_pf      <- ypr_pf
    refCurves$ssbpr_pf    <- ssbpr_pf
    refCurves$Req_pf      <- Req_pf
    refCurves$Beq_pf      <- Beq_pf
    refCurves$tBeq_pf     <- tBeq_pf
    refCurves$expBeq_pf   <- expBeq_pf
    refCurves$Yeq_pf      <- Yeq_pf
    refCurves$Ueq_pf      <- Ueq_pf
    refCurves$Surv_apf    <- Surv_apf


  return( refCurves )
}

# .calcEquil()
# Calculates equilibrium Biomass, Yield and Recruitment 
# with respect to given fishing mortality rates F
# inputs:   f = input fishing mortality rate
#           obj = list of biological parameters
# ouputs:   equil = list() of equilibrium biomass, yield and recruitment
.calcEquil <- function( f = 0, obj )
{
  nP  <- obj$nP
  A   <- obj$nA

  # Recover recruitment pars
  h_p          <- obj$rSteepness_p
  rec.a_p      <- obj$reca_p
  rec.b_p      <- obj$recb_p

  # Pull B0 and R0
  B0_p     <- obj$B0_p
  R0_p     <- obj$R0_p

  # Beverton-Holt a/b parameters
  rec.a_p <- 4.*h_p*R0_p/(B0_p*(1.-h_p))
  rec.b_p <- (5.*h_p-1.)/(B0_p*(1.-h_p))

  # Now calculate eqbm recruitment at given f value
  tmp <- .calcPerRecruit( f = f, obj = obj )
  yprList <- tmp$yprList

  recruits_p <- ( rec.a_p * yprList$ssbpr_p - 1) / (rec.b_p * yprList$ssbpr_p)
  

  equil <- list()
    equil$Req_p     <- recruits_p
    equil$Beq_p     <- recruits_p * yprList$ssbpr_p
    equil$tBeq_p    <- recruits_p * yprList$totbpr_p
    equil$expBeq_p  <- recruits_p * yprList$expbpr_p
    equil$Yeq_p     <- recruits_p * yprList$ypr_p
    equil$ypr_p     <- yprList$ypr_p
    equil$ssbpr_p   <- yprList$ssbpr_p
    equil$Ueq_p     <- equil$Yeq_p / equil$expBeq_p
    equil$Surv_ap   <- tmp$Surv_ap

  return(equil)
}


.calcSel_p <- function( obj, fleetIdx = 1)
{
  # Model dimensions
  nP    <- obj$nP
  nA    <- obj$nA
  age   <- obj$age
  len   <- obj$len


  # Selectivity - this is mean over each fleet's 
  # time period, not time-varying
  # Might be useful to take the group mean for the comm fleet...
  xSelAlpha_p   <- obj$SelAlpha_gp[ , fleetIdx ]
  xSelBeta_p    <- obj$SelBeta_gp[ , fleetIdx ]
  
  # Calculate selectivity as dome/asymp and age/length based
  selType       <- obj$selType_g[fleetIdx]  # 0 == asymptotic, 1 == normal
  selX          <- obj$selX_g[fleetIdx]     # 0 == age, 1 == length


  selAge_ap     <-  array(NA, dim = c(nA,nP) )

  # Loop over species and stocks, calculate
  # selLen so we can get selAge
  for(p in 1:nP )
  {
    # Asymptotic
    if( selType == 0 )
    {
      tmp <- log(19) / xSelBeta_p[p]

      # Age based
      if( selX == 0)
        selAge_ap[,p] <- 1 / (1 + exp(-tmp * (age - xSelAlpha_p[p])))

      # length based
      if( selX == 1)
        selAge_ap[,p] <- 1 / (1 + exp(-tmp * (len - xSelAlpha_p[p])))
    }

    # Dome shaped (normal)
    if( selType == 1 )
    {
      # Age
      if( selX == 0 )
        selAge_ap[,p] <- exp(-1. * ((age - xSelAlpha_p[p])/xSelBeta_p[p]) )

      # length based
      if( selX == 1)
        selAge_ap[,p] <- exp(-1. * ((len - xSelAlpha_p[p])/xSelBeta_p[p]) )

    }
  }
  
  obj$selAge_ap <- selAge_ap

  obj
}

# .calcPerRecruit
# Purpose:     Calculate all equilibrium per-recruit quantities of interest for an
#              input fishing mortality.
# Parameters:  f=scalar input fishing mortality rate; obj=list of all operating
#              model parameters.
# Returns:     a list with equilibrium quantities - (i) spawning stock biomass-per-recruit
#              and (ii) yield-per-recruit (ypr)
# Source:      S.P. Cox, modified for hierSCAL by SDNJ
.calcPerRecruit <- function( f, obj )
{

  # Compute eqbm spawning biomass per recruit for
  # given f and species/stock pars
  nP    <- obj$nP
  nA    <- obj$nA
  M_p   <- obj$meanM_p
  nT    <- obj$nT

  # Life history schedules
  matAge_a          <- obj$mat
  wtAge_ap          <- obj$meanWt_ap
  selAge_ap         <- array( NA, dim = c(nA,nP))
  selAge_ap[,1:nP]  <- obj$sel_apgt[,,2,nT]

  matAge_ap <- cbind(matAge_a,matAge_a,matAge_a)

  # Compute Z_asp
  Z_ap      <- array( NA, dim = c(nA,nP))
  Surv_ap   <- array( 1, dim = c(nA,nP))
  
  for( p in 1:nP )
  {
    Z_ap[,p] <- M_p[p]
    for( a in 1:nA)
    {
      Z_ap[a,p] <- Z_ap[a,p] + selAge_ap[a,p] * f
      if( a > 1 )
        Surv_ap[a,p] <- Surv_ap[a-1,p] * exp( -Z_ap[a-1,p])
      if( a == nA)
        Surv_ap[a,p] <- Surv_ap[a,p] / (1 - exp(-Z_ap[a,p]))
    }
  }

  # Calculate yield-per-recruit
  C_ap    <- Surv_ap * wtAge_ap * selAge_ap * f * (1 - exp(-Z_ap))/Z_ap
  # Replace NAs with 0 (unmodeled ages)
  Z_ap[is.na(Z_ap)] <- 0
  Surv_ap[is.na(Surv_ap)] <- 0
  C_ap[is.na(C_ap)] <- 0
  
  # Compute YPR
  ypr_p    <- apply( X = C_ap, FUN = sum, MARGIN = c(2),na.rm = T)
  
  totbpr_ap <- Surv_ap * wtAge_ap 
  totbpr_p <- apply( X = totbpr_ap, FUN = sum, MARGIN = 2)

  expbpr_ap <- Surv_ap * selAge_ap * wtAge_ap 
  expbpr_p <- apply( X = expbpr_ap, FUN = sum, MARGIN = 2)

  # spawning biomass per recruit
  ssbpr_ap  <- Surv_ap * wtAge_ap * matAge_a * exp( - Z_ap )
  ssbpr_p   <- apply( X = ssbpr_ap, FUN = sum, MARGIN = c(2), na.rm = T )

  # compile output list
  yprList <- list(  ssbpr_p  = ssbpr_p,
                    ypr_p    = ypr_p,
                    totbpr_p = totbpr_p,
                    expbpr_p = expbpr_p )

  obj$yprList <- yprList
  obj$Surv_ap <- Surv_ap

  return(obj)
}

# # Calculates recruitment parameters, and equilibrium unfished
# # numbers and recruitment.
# .calcRecPars <- function( obj )
# {
#   # Calculate eqbm parameters
#   # Survival

    
#   # Beverton-Holt a parameters
#   rec.a <- 4.*obj$rSteepness*R0/(obj$B0*(1.-obj$rSteepness))
#   rec.b <- (5.*obj$rSteepness-1.)/(obj$B0*(1.-obj$rSteepness))

#   # Now return everything in a list
#   recList <- list(  S0 = S0,
#                     wbar0 = wbar0,
#                     N0 = N0,
#                     R0 = R0,
#                     rec.a = rec.a,
#                     rec.b = rec.b  )

#   return(recList)
# }

# .getFmsy     ()
# Purpose:     fit a spline function to f vs yield, then use a root finder to get Fmsy. 
# Parameters:  obj=list of all operating model parameters, schedules, equilibrium functions
# Returns:     a list with all equilibrium quantities for Fmsy
# Source:      S.P. Cox, modified for hierSCAL by SDNJ
.getFmsy_p <- function( obj, refCurves )
{
  Fseq  <- refCurves$F
  nP    <- obj$nP

  .getFmsy <- function( yieldCurve, F = Fseq )
  {
    minF <- 0
    maxF <- max(F)

    # Now check if yieldCurve goes negative anywhere
    if( any(yieldCurve < 0) )
    {
      minF <- F[min(which(yieldCurve >= 0))]
      maxF <- F[max(which(yieldCurve >= 0))]
    }

    fySplineFun <- splinefun( x=F, y=yieldCurve )  

    # Find stat point for Fmsy
    Fmsy <- try( uniroot( f = fySplineFun, interval = c(minF, maxF),
                          deriv = 1 )$root )
    if(class(Fmsy) == "try-error")
    {
      browser(cat("uniroot failed\n\n"))
      Fmsy <- 0
    }

    Fmsy <- min(Fmsy, maxF)

    return(Fmsy)
  }

  # calculate Fmsy for each species/stock
  Fmsy_p <- apply(  X = refCurves$Yeq_pf, FUN = .getFmsy,
                    MARGIN = c(1) )
  # Now create a vector to hold stock values
  # on each curve
  pVec       <- numeric(length = nP)
  names(pVec) <- names(Fmsy_p)

  FmsyRefPts  <- list(  yprFmsy_p     = pVec,
                        ssbprFmsy_p   = pVec,
                        YeqFmsy_p     = pVec,
                        BeqFmsy_p     = pVec,
                        ReqFmsy_p     = pVec,
                        tBeqFmsy_p    = pVec,
                        expBeqFmsy_p  = pVec,
                        Umsy_p        = pVec )

  
  # Calculate ref points
  FmsyRefPts$Fmsy_p    <- Fmsy_p

  # Loop and get each stock's ref pt
  for( p in 1:nP )
  {
    tmp <- .calcEquil( f = Fmsy_p[p], obj = obj )
    
    FmsyRefPts$yprFmsy_p[p]     <- tmp$ypr_p[p]
    FmsyRefPts$ssbprFmsy_p[p]   <- tmp$ssbpr_p[p]
    FmsyRefPts$YeqFmsy_p[p]     <- tmp$Yeq_p[p]
    FmsyRefPts$BeqFmsy_p[p]     <- tmp$Beq_p[p]
    FmsyRefPts$tBeqFmsy_p[p]    <- tmp$tBeq_p[p]
    FmsyRefPts$expBeqFmsy_p[p]  <- tmp$expBeq_p[p]
    FmsyRefPts$NeqFmsy_p[p]     <- tmp$Neq_p[p]
    FmsyRefPts$ReqFmsy_p[p]     <- tmp$Req_p[p]
    FmsyRefPts$Umsy_p[p]        <- tmp$Ueq_p[p]
  }

  FmsyRefPts
}     # END function .getFmsy


