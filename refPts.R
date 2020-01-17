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
  obj <- .calcSel_xsp(obj, fleetIdx = 2)

  # Calculate R0_sp
  h_sp      <- obj$h_sp
  B0_sp     <- obj$B0_sp

  # temporarily use calcPerRecruit to recalc R0
  tmp           <- .calcPerRecruit( f = 0, obj = obj )
  yprList       <- tmp$yprList
  R0_sp         <- B0_sp / yprList$ssbpr_sp

  # Beverton-Holt a/b parameters
  obj$rec.a_sp  <- 4.*h_sp*R0_sp/(B0_sp*(1.-h_sp))
  obj$rec.b_sp  <- (5.*h_sp-1.)/(B0_sp*(1.-h_sp))

  # Calculate generation time
  obj$genTime_sp <- .calcGenTime(obj)

  # Calculate reference curves
  refCurves <- .calcRefCurves( obj )

  # First, let's just do Fmsy reference points
  FmsyRefPts <- .getFmsy_sp(  obj = obj, 
                                refCurves = refCurves )
  
  obj$refPts <- list()
  obj$refPts$refCurves    <- refCurves
  obj$refPts$FmsyRefPts   <- FmsyRefPts
  # Get survivorship
  obj$refPts$Surv_axsp    <- tmp$Surv_axsp
  obj$refPts$ssbpr_sp     <- yprList$ssbpr_sp
  obj$refPts$R0_sp        <- R0_sp  
  obj$refPts$rec.a_sp     <- obj$rec.a_sp
  obj$refPts$rec.b_sp     <- obj$rec.b_sp
  
  return(obj)

}

# .calcRefCurves()
# Calculates equilibrium curves of equilbrium biomass, numbers,
# yield and recruitment as a function of input fishing mortality rates
# inputs:   obj = list of biological parameters
# ouputs:   refCurves = list() of reference curves (vectors)
.calcRefCurves <- function( obj, nFs = 1000 )
{
  # First, compute max F (tolerance of 1e-5)
  maxF <- max( 10*obj$M_xsp )

  # We're going to need to fill each species' ref curves,
  # so labeling and dimensions are needed
  nS          <- obj$nS
  nP          <- obj$nP
  sexNames    <- dimnames(obj$M_xsp)[[1]]
  specNames   <- dimnames(obj$M_xsp)[[2]]
  stockNames  <- dimnames(obj$M_xsp)[[3]]
  


  f <- seq( from = 0.0, to = maxF, length = nFs )

  # Create matrices to hold Recruitment reference curve, name rows and copy
  # for each of Beq, Neq and Yeq
  Req_spf      <- array( NA,  dim = c(nS, nP, nFs),
                              dimnames = list(  species = specNames,
                                                stock = stockNames,
                                                F = f ) )
  Beq_spf       <- Req_spf
  expBeq_spf    <- Req_spf
  Yeq_spf       <- Req_spf
  ypr_spf       <- Req_spf
  ssbpr_spf     <- Req_spf
  Ueq_spf       <- Req_spf

  # Loop and fill
  for( i in 1:length(f) )
  {
    tmp         <- .calcEquil( f = f[i], obj = obj )
    Req_spf[,,i]    <- tmp$Req_sp
    Beq_spf[,,i]    <- tmp$Beq_sp
    expBeq_spf[,,i] <- tmp$expBeq_sp
    Yeq_spf[,,i]    <- tmp$Yeq_sp
    ypr_spf[,,i]    <- tmp$ypr_sp
    ssbpr_spf[,,i]  <- tmp$ssbpr_sp
    Ueq_spf[,,i]    <- tmp$Ueq_sp
  }

  refCurves <- list()
    refCurves$F           <- f
    refCurves$ypr_spf     <- ypr_spf
    refCurves$Req_spf     <- Req_spf
    refCurves$Beq_spf     <- Beq_spf
    refCurves$expBeq_spf  <- expBeq_spf
    refCurves$Yeq_spf     <- Yeq_spf
    refCurves$Ueq_spf     <- Yeq_spf

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
  nS  <- obj$nS
  nP  <- obj$nP
  A_s <- obj$A_s


  # Now calculate eqbm recruitment at given f value
  tmp <- .calcPerRecruit( f = f, obj = obj )
  yprList <- tmp$yprList

  recruits_sp <- ( obj$rec.a_sp * yprList$ssbpr_sp - 1) / (obj$rec.b_sp * yprList$ssbpr_sp)
  

  equil <- list()
    equil$Req_sp     <- recruits_sp
    equil$Beq_sp     <- recruits_sp * yprList$ssbpr_sp
    equil$expBeq_sp  <- recruits_sp * yprList$expbpr_sp
    equil$Yeq_sp     <- recruits_sp * yprList$ypr_sp
    equil$ypr_sp     <- yprList$ypr_sp
    equil$ssbpr_sp   <- yprList$ssbpr_sp
    equil$Ueq_sp     <- equil$Yeq_sp / equil$expBeq_sp

  return(equil)
}


.calcSel_xsp <- function( obj, fleetIdx = 2)
{
  # Model dimensions
  nS    <- obj$nS
  nP    <- obj$nP
  nX    <- obj$nX
  A_s   <- obj$A_s
  nA    <- obj$nA
  nL    <- obj$nL

  # Probability of being a given length-at-age
  probLenAge_laspx <- obj$probLenAge_laspx

  # Selectivity - this is mean over each fleet's 
  # time period, not time-varying
  # Might be useful to take the group mean for the comm fleet...
  xSel50_sp   <- obj$xSel50_spf[ , , fleetIdx, drop = FALSE ]
  xSel95_sp   <- obj$xSel95_spf[ , , fleetIdx, drop = FALSE ]
  xSelStep_sp <- obj$xSelStep_spf[ , , fleetIdx, drop = FALSE ]

  # Harcode for comm.mod for now,
  # and length based selectivity
  selLen_lsp  <- array(NA, dim = c(nL,nS,nP) )
  selAge_axsp <- array(NA, dim = c(nA,nX,nS,nP) )

  # Loop over species and stocks, calculate
  # selLen so we can get selAge
  for( s in 1:nS )
    for(p in 1:nP )
    {
      selLen_lsp[,s,p] <- (1 + exp(-(1:nL - xSel50_sp[s,p,])/xSelStep_sp[s,p,]/log(19)))^(-1)
      for( x in 1:nX)
      {
        for( a in 1:nA )
        {
          selAge_axsp[a,x,s,p] <- sum(probLenAge_laspx[,a,s,p,x] * selLen_lsp[,s,p])
        }
        selAge_axsp[,x,s,p] <- selAge_axsp[,x,s,p] / max(selAge_axsp[,,s,p],na.rm = T)
      }
    }

  obj$selLen_lsp <- selLen_lsp
  obj$selAge_axsp <- selAge_axsp

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
  nS    <- obj$nS
  nP    <- obj$nP
  A_s   <- obj$A_s
  nA    <- obj$nA
  nL    <- obj$nL
  nX    <- obj$nX
  M_xsp <- obj$M_xsp

  # Life history schedules
  matAge_asp        <- obj$matAge_asp
  lenAge_axsp       <- aperm(obj$lenAge_aspx,c(1,4,2,3))
  wtAge_axsp        <- aperm(obj$meanWtAge_aspx,c(1,4,2,3))
  probLenAge_laspx  <- obj$probLenAge_laspx
  selAge_axsp       <- obj$selAge_axsp

  # Double matAge for easy copying later
  matAge_axsp       <- array(NA, dim = c(nA,nX,nS,nP))
  for( x in 1:nX)
    matAge_axsp[,x,,] <- matAge_asp

  # Compute Z_asp
  Z_axsp    <- array( NA, dim = c(nA,nX,nS,nP))
  Surv_axsp <- array( NA, dim = c(nA,nX,nS,nP))
  Surv_axsp[1,,,] <- 1/nX
  for( x in 1:nX )
    for( s in 1:nS )
      for( p in 1:nP )
      {
        Z_axsp[1:A_s[s],x,s,p] <- M_xsp[x,s,p]
        for( a in 1:A_s[s])
        {
          Z_axsp[a,x,s,p] <- Z_axsp[a,x,s,p] + selAge_axsp[a,x,s,p] * f
          if( a > 1 )
            Surv_axsp[a,x,s,p] <- Surv_axsp[a-1,x,s,p] * exp( -Z_axsp[a-1,x,s,p])
          if( a == A_s[s])
            Surv_axsp[a,x,s,p] <- Surv_axsp[a,x,s,p] / (1 - exp(-Z_axsp[a,x,s,p]))
        }
      }

  # Calculate yield-per-recruit
  C_axsp <- array(0, dim = dim(Surv_axsp))
  for( x in 1:nX)
    C_axsp[,x,,]    <- Surv_axsp[,x,,,drop = FALSE] * wtAge_axsp[,x,,,drop = FALSE] * 
                          selAge_axsp[,x,,,drop = FALSE] * f * 
                          (1 - exp(-Z_axsp[,x,,,drop = FALSE]))/Z_axsp[,x,,,drop = FALSE]
  # Replace NAs with 0 (unmodeled ages)
  Z_axsp[is.na(Z_axsp)] <- 0
  Surv_axsp[is.na(Surv_axsp)] <- 0
  C_axsp[is.na(C_axsp)] <- 0
  
  # Compute YPR
  ypr_sp    <- apply( X = C_axsp, FUN = sum, MARGIN = c(3,4),na.rm = T)

  # spawning biomass per recruit
  ssbpr_asp  <- Surv_axsp[,nX,,] * wtAge_axsp[,nX,,] * matAge_axsp[,nX,,]
  ssbpr_sp   <- apply( X = ssbpr_asp, FUN = sum, MARGIN = c(2,3), na.rm = T )

  expbpr_axsp <- Surv_axsp * selAge_axsp * wtAge_axsp
  expbpr_sp   <- apply( X = expbpr_axsp, FUN = sum, MARGIN = c(3,4), na.rm = T )  

  # compile output list
  yprList <- list(  ssbpr_sp  = ssbpr_sp,
                    ypr_sp    = ypr_sp,
                    expbpr_sp = expbpr_sp )

  obj$yprList <- yprList
  obj$Surv_axsp <- Surv_axsp

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
.getFmsy_sp <- function( obj, refCurves )
{
  Fseq  <- refCurves$F
  nS    <- obj$nS
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
  Fmsy_sp <-  apply(  X = refCurves$Yeq_spf, FUN = .getFmsy,
                      MARGIN = c(1,2) )
  # Now create a matrix to hold the species/stock values
  # on each curve
  spMat       <- matrix(0, nrow = nS, ncol = nP)
  dimnames(spMat) <- dimnames(Fmsy_sp)

  FmsyRefPts  <- list(  yprFmsy_sp    = spMat,
                        YeqFmsy_sp    = spMat,
                        BeqFmsy_sp    = spMat,
                        ReqFmsy_sp    = spMat,
                        expBeqFmsy_sp = spMat,
                        Umsy_sp       = spMat )

  
  # Calculate ref points
  FmsyRefPts$Fmsy_sp    <- Fmsy_sp

  # Loop and get each stock/species ref pt
  for( s in 1:nS )
    for( p in 1:nP )
    {
      tmp <- .calcEquil( f = Fmsy_sp[s,p], obj = obj )
      FmsyRefPts$yprFmsy_sp[s,p]    <- tmp$ypr_sp[s,p]
      FmsyRefPts$YeqFmsy_sp[s,p]    <- tmp$Yeq_sp[s,p]
      FmsyRefPts$BeqFmsy_sp[s,p]    <- tmp$Beq_sp[s,p]
      FmsyRefPts$expBeqFmsy_sp[s,p] <- tmp$expBeq_sp[s,p]
      FmsyRefPts$NeqFmsy_sp[s,p]    <- tmp$Neq_sp[s,p]
      FmsyRefPts$ReqFmsy_sp[s,p]    <- tmp$Req_sp[s,p]
      FmsyRefPts$Umsy_sp[s,p]       <- tmp$Ueq_sp[s,p]
    }

  FmsyRefPts
}     # END function .getFmsy


# .calcGenTime
# Purpose:     Calculate generation time as average age of mature stock
# Parameters:  natural mortality and maturity
# Returns:     generation time
# Source:      S.P. Cox
.calcGenTime <- function( obj )
{
  # Compute eqbm spawning biomass per recruit for
  # given f and species/stock pars
  nS    <- obj$nS
  nP    <- obj$nP
  A_s   <- obj$A_s
  nA    <- obj$nA
  nL    <- obj$nL
  nX    <- obj$nX
  M_xsp <- obj$M_xsp

  # Life history schedules
  matAge_asp        <- obj$matAge_asp

  genTime_sp <- array(0, dim = c(nS,nP))

  for(s in 1:nS )
    for( p in 1:nP )
    {
      surv <- rep(1, length = A_s[s])
      a <- c(1:A_s[s])
      surv[a] <- exp( -M_xsp[nX,s,p] * (a - 1) )
      surv[A_s[s]] <- surv[A_s[s]] / (1 - exp(-M_xsp[nX,s,p]))
      genTime_sp[s,p] <- sum( a * matAge_asp[a,s,p] * surv) / sum( surv * matAge_asp[a,s,p])
    }

  genTime_sp

  
}
