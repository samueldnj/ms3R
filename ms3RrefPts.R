# --------------------------------------------------------------------------
# ms3RrefPts.R
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
  # obj <- .calcSel_xsp(obj, fleetIdx = 2)

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

  # EBSBpars
  EBSBpars <- calcJABBASelPars(obj)

  # First, let's just do Fmsy reference points
  FmsyRefPts <- .getFmsy_sp(  obj = obj, 
                              refCurves = refCurves )

  # Calculate economic yield curves from ref curves
  econYieldCurves <- .calcEconomicYieldCurves(  obj, 
                                                refCurves$EffCurves,
                                                FmsyRefPts )


  tmpEmsyRefPts <- .getEmsy_p(  obj = obj,
                                refCurves = refCurves$EffCurves )


  
  refPts <- list()
  refPts$refCurves    <- refCurves
  refPts$FmsyRefPts   <- FmsyRefPts
  refPts$EmsyRefPts   <- tmpEmsyRefPts$EmsyRefPts
  refPts$EmsyMSRefPts <- tmpEmsyRefPts$EmsyMSRefPts
  refPts$EmeyRefPts   <- econYieldCurves
  # Get survivorship
  refPts$Surv_axsp    <- tmp$Surv_axsp
  refPts$ssbpr_sp     <- yprList$ssbpr_sp
  refPts$R0_sp        <- R0_sp  
  refPts$rec.a_sp     <- obj$rec.a_sp
  refPts$rec.b_sp     <- obj$rec.b_sp
  refPts$B0_sp        <- refCurves$Beq_spf[,,1]

  # Get EBSB pars for JABBA-Select application
  refPts$EBSBpars <- EBSBpars

  return(refPts)

} # END calcRefPts()


# solveSpline()
# Solves a spline (or derivatives) for a given value within certain bounds.
solveSpline <- function(  Yvals, Xvals, value = 0, bounds = c(0,10), 
                          deriv = 0 )
{
  Yvals <- Yvals - value

  xySplineFun <- splinefun( x=Xvals, y=Yvals )  

  if(is.null(bounds))
  {
    UBidx <- min(which(Yvals>0))
    bounds <- c(0,Xvals[UBidx])
  }

  # Find stat point for Fmsy
  soln <- try( uniroot( f = xySplineFun, interval = bounds,
                        deriv = deriv )$root )

  solVal <- xySplineFun(soln) + value

  outList <- list( val = solVal, soln = soln)

  return(outList)
} # END solveEff()



# .calcEconomicYieldCurves()
# Name says it all. Uses MS effort/yield curves
# and calculates economic yield as profit, i.e. profit = revenue - cost.
# Right now, revenue is an average price per kg, and cost is
# an average price per 1000 hours of trawling, back calculated
# from Nelson 2009 (Financial Profile of Pacific Fisheries)
.calcEconomicYieldCurves <- function( obj, refCurves, FmsyRefPts )
{
  # Model dims
  nF <- obj$nF
  nS <- obj$nS
  nP <- obj$nP

  # Pull effort/yield curves, then for each
  # area, solve for the effort related to
  # the given landings (average vessel landings
  # in 2009)
  Yeq_spe     <- refCurves$Yeq_spe
  Beq_spe     <- refCurves$Beq_spe
  expBeq_spe  <- refCurves$expBeq_spe
  Yeq_spe[Yeq_spe < 0] <- 0
  Eff     <- refCurves$E

  stockNames <- dimnames(Yeq_spe)[[2]]

  opMod         <- obj$opMod
  crewShare     <- opMod$crewShare
  invlambda_s   <- opMod$invlambda_s
  alpha_s       <- opMod$alpha_s
  adjC_s        <- opMod$adjC_s
  incomeCoeff_s <- opMod$incomeCoeff_s
  bcIncome      <- opMod$bcIncome

  # Calc total system yield
  totYeq_pe <- apply(X = Yeq_spe, FUN = sum, MARGIN = c(2,3), na.rm = T)

  # solveEff_p <- apply(  X = totYeq_pe, FUN = solveSpline, MARGIN = 1,
  #                       Xvals = Eff, value = opMod$landingsForFuelCalcKt,
  #                       bounds = NULL, deriv = 0 )

  eff_p <- obj$om$E_pft[,2,54]
  C_p   <- apply( X = obj$om$C_spt[,,54], FUN = sum, MARGIN = 2 )

  # Now, we need to convert solveEff_p to the cost of
  # fuel, using the landings and fuel cost per kt
  effortPrice_p <- opMod$varCostPerKt * C_p / eff_p

  # Need to make into 2016 prices - maybe include in priceModel list

  # Array for holding effort cost
  effCost_pe <- array(NA, dim       = c(nP, length(Eff)),
                          dimnames  = list( stock = stockNames,
                                            E = Eff ) )

  # Get Fmsy yield
  YeqFmsy_sp  <- FmsyRefPts$YeqFmsy_sp

  econYeq_spe <- array(0,dim = dim(Yeq_spe), dimnames = dimnames(Yeq_spe))
  econRev_spe <- array(0,dim = dim(Yeq_spe), dimnames = dimnames(Yeq_spe))
  landVal_spe <- array(0,dim = dim(Yeq_spe), dimnames = dimnames(Yeq_spe))

  # Yeq_spe[is.na(Yeq_spe)]  <- 0
  specNames   <- dimnames(Yeq_spe)[[1]]
  price_s     <- opMod$price_s[specNames]
  Yeq_se      <- apply(X = Yeq_spe, FUN = sum, MARGIN = c(1,3))
  for( p in 1:nP )
  {
    effCost_pe[p,] <- effortPrice_p[p] * Eff

    for( s in 1:nS )
    {
      lnP_e     <- alpha_s[s] + invlambda_s[s] * log(Yeq_se[s,] + adjC_s[s]) + incomeCoeff_s[s] * log(bcIncome)
      P_e       <- exp(lnP_e)

      landVal_spe[s,p,] <- P_e


      econRev_spe[s,p,] <- Yeq_spe[s,p,] * landVal_spe[s,p,] * (1 - crewShare)
      econYeq_spe[s,p,] <- Yeq_spe[s,p,] * landVal_spe[s,p,] * (1 - crewShare)  - effCost_pe[p,]
    }
    
  }
  # Now sum across species for area specific yield
  econRev_pe <- apply(X = econRev_spe, FUN = sum, MARGIN = c(2,3), na.rm = T)
  econYeq_pe <- econRev_pe  - effCost_pe
  dimnames(econYeq_pe) <- dimnames(effCost_pe)
  dimnames(econRev_pe) <- dimnames(effCost_pe)

  
  # Emey_sp <- apply( X = econYeq_spe, FUN = solveSpline, Xvals = Eff,
  #                   deriv = 1, bounds = c(0,15), MARGIN = c(1,2) )

  meySolList <- apply(  X = econYeq_pe, FUN = solveSpline, Xvals = Eff,
                        deriv = 1, bounds = c(0,50), MARGIN = c(1),
                        value = 0 )

  Emey_p  <- numeric(length = nP)
  MEY_p   <- numeric(length = nP)
  
  # Get species specific vals

  Bmey_sp   <- array(0, dim = c(nS,nP))
  vBmey_sp  <- array(0, dim = c(nS,nP))
  Ymey_sp   <- array(0, dim = c(nS,nP))

  for( p in 1:nP )
  {
    Emey_p[p] <-  meySolList[[p]]$soln
    MEY_p[p]  <-  meySolList[[p]]$val

    # Now calculate Bmey, vBmey etc
    for( s in 1:nS )
    {
      Bmey_sp[s,p]  <- getSplineVal(x = Eff, y = Beq_spe[s,p,], p = Emey_p[p] )
      Ymey_sp[s,p]  <- getSplineVal(x = Eff, y = Yeq_spe[s,p,], p = Emey_p[p] )
      vBmey_sp[s,p] <- getSplineVal(x = Eff, y = expBeq_spe[s,p,], p = Emey_p[p] )
    }
  }

  if(is.null(opMod$cwMEY))
    opMod$cwMEY <- FALSE

  if( opMod$cwMEY )
  {

    # Need to define an objective function for
    # optimising coastwide effort with a joint
    # downward sloping demand curve - can extend
    # to a proportion estimator later to 
    EmeyCW_optObj <- try(optim( par = log(Emey_p),
                                fn = coastWideEconObjFun,
                                method = "BFGS",
                                prop = FALSE,
                                baseE = 1,
                                # control = list(trace = 6),
                                optimise = TRUE,
                                invlambda_s = invlambda_s,
                                adjC_s = adjC_s,
                                alpha_s = alpha_s,
                                incomeCoeff_s = incomeCoeff_s,
                                income = bcIncome,
                                effortPrice_p = effortPrice_p,
                                MSY_sp = YeqFmsy_sp,
                                refCurves = refCurves,
                                crewShare = crewShare,
                                nF = nF, nS = nS, nP = nP ) )
    if( class(EmeyCW_optObj) == "try-error" )
      browser()

    cwMEYlist <- coastWideEconObjFun( lnE_p = EmeyCW_optObj$par,
                                      prop = FALSE,
                                      baseE = 1,
                                      optimise = FALSE,
                                      invlambda_s = invlambda_s,
                                      adjC_s = adjC_s,
                                      alpha_s = alpha_s,
                                      incomeCoeff_s = incomeCoeff_s,
                                      income = bcIncome,
                                      effortPrice_p = effortPrice_p,
                                      MSY_sp = YeqFmsy_sp,
                                      refCurves = refCurves,
                                      crewShare = crewShare,
                                      nF = nF, nS = nS, nP = nP)

    cwMEY_p <- cwMEYlist$Rent_p
    cwEmey_p <- cwMEYlist$E_p

    
    cwRev_spe     <- array(NA, dim = dim(Yeq_spe))
    cwRent_pe     <- array(NA, dim = dim(econYeq_pe))
    cwEffCost_pe  <- array(NA, dim = dim(econYeq_pe))
    cwEff_pe      <- array(NA, dim = dim(econYeq_pe))
    cwRent_e      <- array(NA, dim = length(Eff) )
    cwYeq_spe     <- array(NA, dim = dim(Yeq_spe))

    cwBmey_sp     <- array(NA, dim = dim(Bmey_sp))
    cwYmey_sp     <- array(NA, dim = dim(Bmey_sp))
    cwvBmey_sp    <- array(NA, dim = dim(Bmey_sp))

    cwRev_spe[,,1]    <- 0
    cwRent_pe[,1]     <- 0
    cwRent_e[1]       <- 0
    cwYeq_spe[,,1]    <- 0
    cwEffCost_pe[,1]  <- 0
    cwEff_pe[,1]      <- 0

    for( eIdx in 2:length(Eff) )
    {
      e <- Eff[eIdx]

      # First, optimise allocation of effort
      optObj <- try(optim(  par = rep(0,3),
                        fn = coastWideEconObjFun,
                        method = "BFGS",
                        prop = TRUE,
                        baseE = e,
                        # control = list(trace = 6),
                        optimise = TRUE,
                        invlambda_s = invlambda_s,
                        adjC_s = adjC_s,
                        alpha_s = alpha_s,
                        incomeCoeff_s = incomeCoeff_s,
                        income = bcIncome,
                        effortPrice_p = effortPrice_p,
                        refCurves = refCurves,
                        crewShare = crewShare,
                        nF = nF, nS = nS, nP = nP ))
      if(class(optObj)== "try-error")
        break

      cwMEYlist <- coastWideEconObjFun( lnE_p = EmeyCW_optObj$par,
                                        prop = TRUE,
                                        baseE = e,
                                        optimise = FALSE,
                                        invlambda_s = invlambda_s,
                                        adjC_s = adjC_s,
                                        alpha_s = alpha_s,
                                        incomeCoeff_s = incomeCoeff_s,
                                        income = bcIncome,
                                        effortPrice_p = effortPrice_p,
                                        MSY_sp = YeqFmsy_sp,
                                        refCurves = refCurves,
                                        crewShare = crewShare,
                                        nF = nF, nS = nS, nP = nP)


      cwRev_spe[,,eIdx]     <- cwMEYlist$Rev_sp
      cwRent_pe[,eIdx]      <- cwMEYlist$Rent_p
      cwRent_e[eIdx]        <- cwMEYlist$totRent
      cwYeq_spe[,,eIdx]     <- cwMEYlist$C_sp
      cwEffCost_pe[,eIdx]   <- cwMEYlist$effCost_p
      cwEff_pe[,eIdx]       <- cwMEYlist$E_p

    }
    # Numerically solve for biomass and catch at MEY
    for( p in 1:nP )
    {
      # Now calculate Bmey, vBmey etc
      for( s in 1:nS )
      {
        cwBmey_sp[s,p]  <- getSplineVal(x = Eff, y = Beq_spe[s,p,], p = cwEmey_p[p] )
        cwYmey_sp[s,p]  <- getSplineVal(x = Eff, y = Yeq_spe[s,p,], p = cwEmey_p[p] )
        cwvBmey_sp[s,p] <- getSplineVal(x = Eff, y = expBeq_spe[s,p,], p = cwEmey_p[p] )
      }
    }

    cwEconYieldCurves <- list()
    cwEconYieldCurves$Eff           <- Eff
    cwEconYieldCurves$cwMEY_p       <- cwMEY_p
    cwEconYieldCurves$cwEmey_p      <- cwEmey_p
    cwEconYieldCurves$cwBmey_sp     <- cwBmey_sp
    cwEconYieldCurves$cwYmey_sp     <- cwYmey_sp
    cwEconYieldCurves$cwvBmey_sp    <- cwvBmey_sp
    cwEconYieldCurves$cwRev_spe     <- cwRev_spe
    cwEconYieldCurves$cwRent_pe     <- cwRent_pe
    cwEconYieldCurves$cwRent_e      <- cwRent_e
    cwEconYieldCurves$cwYeq_spe     <- cwYeq_spe
    cwEconYieldCurves$cwEffCost_pe  <- cwEffCost_pe
    cwEconYieldCurves$cwEff_pe      <- cwEff_pe  
  }



  econYieldCurves <- list()
  econYieldCurves$effCost_pe    <- effCost_pe
  econYieldCurves$crewShare_pe  <- econRev_pe * crewShare
  econYieldCurves$econYeq_spe   <- econYeq_spe
  econYieldCurves$econRev_spe   <- econRev_spe
  econYieldCurves$econYeq_pe    <- econYeq_pe
  econYieldCurves$econRev_pe    <- econRev_pe
  econYieldCurves$landVal_spe   <- landVal_spe
  econYieldCurves$Emey_p        <- Emey_p
  econYieldCurves$MEY_p         <- MEY_p
  econYieldCurves$Bmey_sp       <- Bmey_sp
  econYieldCurves$Ymey_sp       <- Ymey_sp
  econYieldCurves$vBmey_sp      <- vBmey_sp
  econYieldCurves$effortPrice_p <- effortPrice_p
  
  if( opMod$cwMEY )
    econYieldCurves$cwEconYieldCurves <- cwEconYieldCurves

  



  return( econYieldCurves )

} # END .calcEconomicYieldCurves()

# coastwideEconObjFun()
# Objective function/economic model for coastwide
# econ yield, assuming a joint demand curve for
# each species
coastWideEconObjFun <- function(  lnE_p = c(0,0,0),
                                  prop = FALSE,
                                  baseE = 10,
                                  optimise = TRUE,
                                  invlambda_s = c(Inf,Inf,Inf),
                                  alpha_s,
                                  incomeCoeff_s,
                                  effortPrice_p,
                                  adjC_s,
                                  income = bcIncome,
                                  demCurves = NULL,
                                  MSY_sp,
                                  refCurves,
                                  crewShare,
                                  nF, nS, nP )
{
  # First, take exponent of effort
  E_p <- exp(lnE_p)

  # Convert to proportions and multiply by base
  if( prop )
    E_p <- baseE * E_p / (sum(E_p))
  
  # Now compute catch, revenue, and econ yield at this
  # effort
  Eff     <- refCurves$E
  Yeq_spe <- refCurves$Yeq_spe

  C_sp    <- array(NA, dim = c(nS,nP))  
  Rev_sp  <- array(NA, dim = c(nS,nP))  

  for( s in 1:nS )
    for(p in 1:nP )
      C_sp[s,p] <- getSplineVal(x = Eff, y = Yeq_spe[s,p,], p = E_p[p] )

  C_sp[C_sp < 0] <- 0

  # Calculate coastwide catch
  C_s   <- apply( X = C_sp,   FUN = sum, MARGIN = 1 )

  # use coastwide catch to figure out species landed value
  lnP_s     <- alpha_s + invlambda_s * log(C_s + adjC_s+1e-9) + incomeCoeff_s * log(income)
  landVal_s <- exp(lnP_s)

  
  for( s in 1:nS )
    for( p in 1:nP )
    {
      Rev_sp[s,p] <-  C_sp[s,p] * landVal_s[s] * (1 - crewShare)
    }


  # Calculate revenue and rent
  Rev_p   <- apply(X = Rev_sp, FUN = sum, MARGIN = 2 )
  Rent_p  <- Rev_p - E_p * effortPrice_p

  totRent <- sum(Rent_p)

  outList <- list(  landVal_s = landVal_s,
                    E_p = E_p,
                    C_sp = C_sp,
                    Rev_sp = Rev_sp,
                    effCost_p = effortPrice_p * E_p,
                    crewShare = crewShare,
                    Rev_p = Rev_p,
                    Rent_p = Rent_p,
                    totRent = totRent )

  if(optimise)
    return(-totRent )


  return( outList )


} # END coastWideEconObjFun()


calcJABBASelPars <- function( obj )
{
  # Loop over all fleets, and calculate
  # ratios of equilibrium spawning
  # and eqbm exploitable biomass for a vector
  # of fishing mortality rates. Then,
  # fit the JABBA select EB/SB model
  # to the data.
  # Model dims
  nF <- obj$nF
  nS <- obj$nS
  nP <- obj$nP

  # Max F to be calculated for
  maxF <- max( 10*obj$M_xsp )
  Fseq <- seq( from = 0, to = maxF, length.out = 100 )


  Beq_spfj    <- array(NA, dim = c(nS,nP,nF,100) )
  expBeq_spfj <- array(NA, dim = c(nS,nP,nF,100) )
  ratioB_spfj <- array(NA, dim = c(nS,nP,nF,100) )

  nu_spfk <- array(NA, dim = c(nS,nP,nF,3) )
  P1_spf  <- array(NA, dim = c(nS,nP,nF) )
  P2_spf  <- array(NA, dim = c(nS,nP,nF) )

  for( f in 1:nF)
  {
    for( j in 1:100 )
    {
      tmp <- .calcEquil(  f = Fseq[j], obj = obj,
                          fleetIdx = f )
      Beq_spfj[,,f,j]     <- tmp$Beq_sp
      expBeq_spfj[,,f,j]  <- tmp$expBeq_sp

    }
  }

  # Replace crashed stock biomasses with NAs
  Beq_spfj[ Beq_spfj < 0 ]    <- 0
  Beq_spfj[ expBeq_spfj < 0 ] <- 0

  expBeq_spfj[ Beq_spfj < 0 ]     <- 0
  expBeq_spfj[ expBeq_spfj < 0 ]  <- 0

  expBeq_spfj[ expBeq_spfj == 0 ] <- NA
  Beq_spfj[ Beq_spfj == 0 ]       <- NA

  # Compute ratio
  ratioB_spfj <- expBeq_spfj / Beq_spfj
  
  # Now apply a function to solve for nu parameters
  for( s in 1:nS )
    for( p in 1:nP )
      for( f in 1:nF )
      {
        initnu_k <- c(0.25, 1, 2 )

        initlognu_k <- log(initnu_k)

        Peq_j <- Beq_spfj[s,p,f,] / obj$B0_sp[s,p]

        P1_spf[s,p,f] <- min(Peq_j,na.rm = T)
        P2_spf[s,p,f] <- max(Peq_j,na.rm = T)

        optOut <- optim(  par = initlognu_k, fn = .EBSBlikelihood,
                          method = "Nelder-Mead",
                          Beq_j = Beq_spfj[s,p,f,],
                          expBeq_j = expBeq_spfj[s,p,f,],
                          B0 = obj$B0_sp[s,p],
                          P1 = P1_spf[s,p,f],
                          P2 = P2_spf[s,p,f] )



        nu_spfk[s,p,f,] <- exp(optOut$par) 
      }

  # Save model parameters
  EBSB_StockSpec <- list( P1_spf  = P1_spf,
                          P2_spf  = P2_spf,
                          nu_spfk = nu_spfk )

  # Now do coastwide.
  # First, borrow the sub-arrays
  nu_sfk <- nu_spfk[,1,,,drop = FALSE]
  P1_sf  <- P1_spf[,1,,drop = FALSE]
  P2_sf  <- P2_spf[,1,,drop = FALSE]

  Beq_sfj     <- Beq_spfj[,1,,,drop = FALSE]
  Beq_sfj[,1,,] <- apply(X = Beq_spfj, FUN = sum, MARGIN = c(1,3,4))
  expBeq_sfj  <- Beq_spfj[,1,,,drop = FALSE]
  expBeq_sfj[,1,,] <- apply(X = expBeq_spfj, FUN = sum, MARGIN = c(1,3,4))
  # Then loop and compute
  for( s in 1:nS )
  {
    for( f in 1:nF )
    {
      initnu_k <- c(0.25, 1, 2 )

      initlognu_k <- log(initnu_k)

      Peq_j <- Beq_sfj[s,1,f,] / sum(obj$B0_sp[s,])
      P1_sf[s,1,f] <- min(Peq_j,na.rm = T)
      P2_sf[s,1,f] <- max(Peq_j,na.rm = T)

      optOut <- optim(  par = initlognu_k, fn = .EBSBlikelihood,
                        method = "Nelder-Mead",
                        Beq_j = Beq_sfj[s,1,f,],
                        expBeq_j = expBeq_sfj[s,1,f,],
                        B0 = sum(obj$B0_sp[s,]),
                        P1 = P1_sf[s,1,f],
                        P2 = P2_sf[s,1,f] )

      nu_sfk[s,1,f,] <- exp(optOut$par) 
    }
  }

  # Save model parameters
  EBSB_CoastWide <- list( P1_spf  = P1_sf,
                          P2_spf  = P2_sf,
                          nu_spfk = nu_sfk )

  # Now data pooled
  # First, borrow the sub-arrays
  nu_pfk <- nu_spfk[1,,,,drop = FALSE]
  P1_pf  <- P1_spf[1,,,drop = FALSE]
  P2_pf  <- P2_spf[1,,,drop = FALSE]

  Beq_pfj     <- Beq_spfj[1,,,,drop = FALSE]
  Beq_pfj[1,,,] <- apply(X = Beq_spfj, FUN = sum, MARGIN = c(2,3,4))
  expBeq_pfj  <- Beq_spfj[1,,,,drop = FALSE]
  expBeq_pfj[1,,,] <- apply(X = expBeq_spfj, FUN = sum, MARGIN = c(2,3,4))
  # Then loop and compute
  for( p in 1:nP )
  {
    for( f in 1:nF )
    {
      initnu_k <- c(0.25, 1, 2 )

      initlognu_k <- log(initnu_k)

      Peq_j <- Beq_pfj[1,p,f,] / sum(obj$B0_sp[,p])
      P1_pf[1,p,f] <- min(Peq_j,na.rm = T)
      P2_pf[1,p,f] <- max(Peq_j,na.rm = T)

      optOut <- optim(  par = initlognu_k, fn = .EBSBlikelihood,
                        method = "Nelder-Mead",
                        Beq_j = Beq_pfj[1,p,f,],
                        expBeq_j = expBeq_pfj[1,p,f,],
                        B0 = sum(obj$B0_sp[,p]),
                        P1 = P1_pf[1,p,f],
                        P2 = P2_pf[1,p,f] )

      

      nu_pfk[1,p,f,] <- exp(optOut$par) 
    }
  }

  # Save model parameters
  EBSB_DataPooled <- list(  P1_spf  = P1_pf,
                            P2_spf  = P2_pf,
                            nu_spfk = nu_pfk )

  # Now Data Pooled
  # First, borrow the sub-arrays
  nu_fk <- nu_spfk[1,1,,,drop = FALSE]
  P1_f  <- P1_spf[1,1,,drop = FALSE]
  P2_f  <- P2_spf[1,1,,drop = FALSE]

  Beq_fj     <- Beq_spfj[1,1,,,drop = FALSE]
  Beq_fj[1,1,,] <- apply(X = Beq_spfj, FUN = sum, MARGIN = c(3,4))
  expBeq_fj  <- Beq_spfj[1,1,,,drop = FALSE]
  expBeq_fj[1,1,,] <- apply(X = expBeq_spfj, FUN = sum, MARGIN = c(3,4))
  # Then loop and compute
  for( f in 1:nF )
  {
    initnu_k <- c(0.25, 1, 2 )

    initlognu_k <- log(initnu_k)

    Peq_j <- Beq_fj[1,1,f,] / sum(obj$B0_sp)
    P1_f[1,1,f] <- min(Peq_j,na.rm = T)
    P2_f[1,1,f] <- max(Peq_j,na.rm = T)

    optOut <- optim(  par = initlognu_k, fn = .EBSBlikelihood,
                      method = "Nelder-Mead",
                      Beq_j = Beq_fj[1,1,f,],
                      expBeq_j = expBeq_fj[1,1,f,],
                      B0 = sum(obj$B0_sp[,p]),
                      P1 = P1_f[1,1,f],
                      P2 = P2_f[1,1,f] )

    

    nu_fk[1,1,f,] <- exp(optOut$par) 
  }

  # Save model parameters
  EBSB_TotalAgg <- list(  P1_spf  = P1_f,
                          P2_spf  = P2_f,
                          nu_spfk = nu_fk )
  


  # Return Exploitable/spawning biomass
  # ratio model pars
  out <- list(  stockSpec   = EBSB_StockSpec,
                dataPooled  = EBSB_DataPooled,
                coastWide   = EBSB_CoastWide,
                totalAgg    = EBSB_TotalAgg )

  out
} # END .calcJABBASelPars()

# .EBSBratio()
# Parametric function to calculate EBSB
# ratio from a given spawning biomass
# depletion level
.EBSBratio <- function( Pf, P1, P2, nu_k )
{
  # Calculate numerator and denominator
  numerator   <- (1 - exp(-nu_k[3] * (Pf - P1)))
  denominator <- (1 - exp(-nu_k[3] * (P2 - P1)))

  # Calculate expected ratio
  EBSB <- nu_k[1] + (nu_k[2] - nu_k[1]) * numerator / denominator

  EBSB
} # END .EBSBratio

# .EBSBoptim()
# Optimises nu parameters for
# the EBSB model
.EBSBlikelihood <- function( lognu_k, B0, Beq_j, expBeq_j, P1, P2 )
{
  # Exponentiate lognu_k
  nu_k <- exp(lognu_k)

  # First, calculate relative depletion levels
  # for the spawning biomass
  Peq_j <- Beq_j / B0
  Peq_j <- Peq_j[!is.na(Peq_j)]
  
  expRatio <- sapply( X = Peq_j, FUN = .EBSBratio,
                        P1 = P1, P2 = P2, nu_k = nu_k )

  obsRatio <- expBeq_j / Beq_j

  resid <- log(obsRatio[!is.na(Peq_j)]) - log(expRatio[!is.na(Peq_j)])


  SSQ <- sum(resid^2, na.rm = T)

  SSQ
} # END .EBSBlikelihood()



# .calcRefCurves()
# Calculates equilibrium curves of equilbrium biomass, numbers,
# yield and recruitment as a function of input fishing mortality rates
# inputs:   obj = list of biological parameters
# ouputs:   refCurves = list() of reference curves (vectors)
.calcRefCurves <- function( obj, nFs = 500, 
                            maxE = 350, 
                            fleetIdx = 2 )
{
  # First, compute max F (tolerance of 1e-5)
  nT   <- dim(obj$om$qF_spft)[4]
  maxF <- max( 5*obj$M_xsp )

  # We're going to need to fill each species' ref curves,
  # so labeling and dimensions are needed
  nS          <- obj$nS
  nP          <- obj$nP
  nA          <- obj$nA
  nX          <- obj$nX
  sexNames    <- dimnames(obj$M_xsp)[[1]]
  specNames   <- dimnames(obj$M_xsp)[[2]]
  stockNames  <- dimnames(obj$M_xsp)[[3]]
  


  f <- seq( from = 0.0, to = maxF, length = nFs )
  e <- seq( from = 0.0, to = maxE, length = nFs )

  # Create matrices to hold Recruitment reference curve, name rows and copy
  # for each of Beq, Neq and Yeq
  Req_spf      <- array( NA,  dim = c(nS, nP, nFs),
                              dimnames = list(  species = specNames,
                                                stock = stockNames,
                                                F = f ) )

  surv_axspf    <- array( NA, dim = c(nA,nX,nS,nP,nFs) )

  Beq_spf       <- Req_spf
  expBeq_spf    <- Req_spf
  totBeq_spf    <- Req_spf
  Yeq_spf       <- Req_spf
  ypr_spf       <- Req_spf
  ssbpr_spf     <- Req_spf
  Ueq_spf       <- Req_spf


  # Loop and fill
  for( i in 1:length(f) )
  {
    tmp               <- .calcEquil( f = f[i], obj = obj )
    Req_spf[,,i]      <- tmp$Req_sp
    Beq_spf[,,i]      <- tmp$Beq_sp
    expBeq_spf[,,i]   <- tmp$expBeq_sp
    totBeq_spf[,,i]   <- tmp$totBeq_sp
    Yeq_spf[,,i]      <- tmp$Yeq_sp
    ypr_spf[,,i]      <- tmp$ypr_sp
    ssbpr_spf[,,i]    <- tmp$ssbpr_sp
    Ueq_spf[,,i]      <- tmp$Ueq_sp
    surv_axspf[,,,,i] <- tmp$surv_axsp
  }

  Req_spe      <- array( NA,  dim = c(nS, nP, nFs),
                              dimnames = list(  species = specNames,
                                                stock = stockNames,
                                                E = e ) )

  surv_axspe    <- array( NA, dim = c(nA,nX,nS,nP,nFs) )

  Beq_spe       <- Req_spe
  totBeq_spe    <- Req_spe
  expBeq_spe    <- Req_spe
  Yeq_spe       <- Req_spe
  ypr_spe       <- Req_spe
  ssbpr_spe     <- Req_spe
  Ueq_spe       <- Req_spe

  # Loop and fill
  for( i in 1:length(e) )
  {
    tmp               <- .calcEquil( f = e[i], obj = obj, type = "effort" )
    Req_spe[,,i]      <- tmp$Req_sp
    Beq_spe[,,i]      <- tmp$Beq_sp
    expBeq_spe[,,i]   <- tmp$expBeq_sp
    totBeq_spe[,,i]   <- tmp$totBeq_sp
    Yeq_spe[,,i]      <- tmp$Yeq_sp
    ypr_spe[,,i]      <- tmp$ypr_sp
    ssbpr_spe[,,i]    <- tmp$ssbpr_sp
    Ueq_spe[,,i]      <- tmp$Ueq_sp
    surv_axspe[,,,,i] <- tmp$surv_axsp
  }

  # Save F based ref points
  refCurves <- list()
    refCurves$F           <- f
    refCurves$ypr_spf     <- ypr_spf
    refCurves$Req_spf     <- Req_spf
    refCurves$Beq_spf     <- Beq_spf
    refCurves$totBeq_spf  <- totBeq_spf
    refCurves$expBeq_spf  <- expBeq_spf
    refCurves$Yeq_spf     <- Yeq_spf
    refCurves$Ueq_spf     <- Yeq_spf
    refCurves$surv_axspf  <- surv_axspf

  # Save E based ref points
  refCurves$EffCurves <- list()
    refCurves$EffCurves$E           <- e
    refCurves$EffCurves$ypr_spe     <- ypr_spe
    refCurves$EffCurves$Req_spe     <- Req_spe
    refCurves$EffCurves$Beq_spe     <- Beq_spe
    refCurves$EffCurves$totBeq_spe  <- totBeq_spe
    refCurves$EffCurves$expBeq_spe  <- expBeq_spe
    refCurves$EffCurves$Yeq_spe     <- Yeq_spe
    refCurves$EffCurves$Ueq_spe     <- Yeq_spe
    refCurves$EffCurves$surv_axspe  <- surv_axspe    

  return( refCurves )
}

# .calcEquil()
# Calculates equilibrium Biomass, Yield and Recruitment 
# with respect to given fishing mortality rates F
# inputs:   f = input fishing mortality rate
#           obj = list of biological parameters
# ouputs:   equil = list() of equilibrium biomass, yield and recruitment
.calcEquil <- function( f = 0, obj, type = "fmort", 
                        fleetIdx = 2 )
{
  nS  <- obj$nS
  nP  <- obj$nP
  A_s <- obj$A_s


  # Now calculate eqbm recruitment at given f value
  tmp <- .calcPerRecruit( f = f, obj = obj, type = type, fleetIdx = fleetIdx )
  yprList <- tmp$yprList

  recruits_sp <- ( obj$rec.a_sp * yprList$ssbpr_sp - 1) / (obj$rec.b_sp * yprList$ssbpr_sp)
  

  equil <- list()
    equil$Req_sp     <- recruits_sp
    equil$Beq_sp     <- recruits_sp * yprList$ssbpr_sp
    equil$expBeq_sp  <- recruits_sp * yprList$expbpr_sp
    equil$totBeq_sp  <- recruits_sp * yprList$totbpr_sp
    equil$Yeq_sp     <- recruits_sp * yprList$ypr_sp
    equil$ypr_sp     <- yprList$ypr_sp
    equil$ssbpr_sp   <- yprList$ssbpr_sp
    equil$Ueq_sp     <- equil$Yeq_sp / equil$Beq_sp
    equil$surv_axsp  <- tmp$Surv_axsp

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
.calcPerRecruit <- function( f, obj, type = "fmort", fleetIdx = 2 )
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
  nT    <- obj$nT
  qF_sp <- obj$om$qF_spft[,,2,nT+1]

  # Life history schedules
  matAge_asp        <- obj$matAge_asp
  lenAge_axsp       <- obj$lenAge_axsp
  wtAge_axsp        <- obj$meanWtAge_axsp
  probLenAge_laspx  <- obj$probLenAge_laspx
  selAge_axsp       <- obj$sel_axspft[,,,,fleetIdx,nT]

  fmort <- array(f, dim =c(nS,nP) )

  # Compute Z_asp
  Z_axsp    <- array( NA, dim = c(nA,nX,nS,nP))
  Surv_axsp <- array( NA, dim = c(nA,nX,nS,nP))
  Surv_axsp[1,,,] <- 1
  for( x in 1:nX )
    for( s in 1:nS )
      for( p in 1:nP )
      {
        if( type == "effort" )
          fmort[s,p] <- f * qF_sp[s,p]

        Z_axsp[1:A_s[s],x,s,p] <- M_xsp[x,s,p]
        for( a in 1:A_s[s])
        {
          Z_axsp[a,x,s,p] <- Z_axsp[a,x,s,p] + selAge_axsp[a,x,s,p] * fmort[s,p]
          if( a > 1 )
            Surv_axsp[a,x,s,p] <- Surv_axsp[a-1,x,s,p] * exp( -Z_axsp[a-1,x,s,p])
          if( a == A_s[s])
            Surv_axsp[a,x,s,p] <- Surv_axsp[a,x,s,p] / (1 - exp(-Z_axsp[a,x,s,p]))
        }
      }

  # Calculate yield-per-recruit, ssb per recruit, and exp biomass per recruit
  # by using survival array
  ssbpr_asp   <- array( NA, dim = c(nA,nS,nP) )
  totbpr_axsp  <- array( NA, dim = c(nA,nX,nS,nP) )
  expbpr_axsp <- array( NA, dim = c(nA,nX,nS,nP) )
  C_axsp      <- array(0, dim = dim(Surv_axsp))

  for( s in 1:nS )
    for( p in 1:nP)
    {
      ssbpr_asp[,s,p]  <- Surv_axsp[,nX,s,p] * wtAge_axsp[,nX,s,p] * matAge_asp[,s,p]
      for( x in 1:nX )
      {
        C_axsp[,x,s,p]    <- Surv_axsp[,x,s,p] * wtAge_axsp[,x,s,p] * 
                              selAge_axsp[,x,s,p] * fmort[s,p] * 
                              (1 - exp(-Z_axsp[,x,s,p]))/Z_axsp[,x,s,p]

        expbpr_axsp[,x,s,p] <- Surv_axsp[,x,s,p] * selAge_axsp[,x,s,p] * wtAge_axsp[,x,s,p]
        totbpr_axsp[,x,s,p] <- Surv_axsp[,x,s,p] * wtAge_axsp[,x,s,p]
      }

    }
  # Replace NAs with 0 (unmodeled ages)
  # Z_axsp[is.na(Z_axsp)] <- 0
  # Surv_axsp[is.na(Surv_axsp)] <- 0
  # C_axsp[is.na(C_axsp)] <- 0
  
  # Compute ratios
  ypr_sp      <- apply( X = C_axsp, FUN = sum, MARGIN = c(3,4),na.rm = T)
  ssbpr_sp    <- apply( X = ssbpr_asp, FUN = sum, MARGIN = c(2,3), na.rm = T )
  expbpr_sp   <- apply( X = expbpr_axsp, FUN = sum, MARGIN = c(3,4), na.rm = T )  
  totbpr_sp   <- apply( X = totbpr_axsp, FUN = sum, MARGIN = c(3,4), na.rm = T )  

  # compile output list
  yprList <- list(  ssbpr_sp  = ssbpr_sp,
                    ypr_sp    = ypr_sp,
                    expbpr_sp = expbpr_sp,
                    totbpr_sp = totbpr_sp )

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

# .getEmsy_p()
# Calculates effort based reference curves
# for each stock and species, then computes
# multispecies optimum yield and effort
# values for each stock area. Returns
# arrays for plotting of yield curves.
.getEmsy_p <- function( obj, refCurves, FmsyRefPts )
{
  nS  <- obj$nS
  nP  <- obj$nP

  Eseq <- refCurves$E

  .getEmsy <- function( yieldCurve, E = Eseq )
  {
    minE <- 0
    maxE <- max(E)

    # Now check if yieldCurve goes negative anywhere
    if( any(yieldCurve < 0) )
    {
      minE <- E[min(which(yieldCurve >= 0))]
      maxE <- E[max(which(yieldCurve >= 0))]
    }

    eySplineFun <- splinefun( x=E, y=yieldCurve )  

    # Find stat point for Fmsy
    Emsy <- try( uniroot( f = eySplineFun, interval = c(minE, maxE),
                          deriv = 1 )$root )
    if(class(Emsy) == "try-error")
    {
      browser(cat("uniroot failed\n\n"))
      Emsy <- 0
    }

    Emsy <- min(Emsy, maxE)

    return(Emsy)
  }

  # calculate Fmsy for each species/stock
  Emsy_sp <-  apply(  X = refCurves$Yeq_spe, FUN = .getEmsy,
                      MARGIN = c(1,2) )

  Yeq_spe <- refCurves$Yeq_spe
  Yeq_spe[Yeq_spe < 0] <- NA

  # Now we want the total stock area yield curve
  Yeq_pe  <- apply( X  = refCurves$Yeq_spe, FUN = sum, MARGIN = c(2,3), 
                    na.rm = T )
  Yeq_pe[Yeq_pe == 0] <- NA
  Yeq_pe[,1] <- 0
  Emsy_p  <- apply( X = Yeq_pe, FUN = .getEmsy, MARGIN = 1)

  # Now create a matrix to hold the species/stock values
  # on each curve
  spMat       <- matrix(0, nrow = nS, ncol = nP)
  dimnames(spMat) <- dimnames(Emsy_sp)

  EmsyRefPts  <- list(  yprEmsy_sp    = spMat,
                        YeqEmsy_sp    = spMat,
                        BeqEmsy_sp    = spMat,
                        ReqEmsy_sp    = spMat,
                        expBeqEmsy_sp = spMat,
                        totBeqEmsy_sp = spMat,
                        Umsy_sp       = spMat )

  EmsyMSRefPts <- list( yprEmsy_sp    = spMat,
                        YeqEmsy_sp    = spMat,
                        BeqEmsy_sp    = spMat,
                        ReqEmsy_sp    = spMat,
                        expBeqEmsy_sp = spMat,
                        totBeqEmsy_sp = spMat,
                        Umsy_sp       = spMat,
                        YeqEmsy_p     = rep(0,nP) )

  
  # Calculate ref points
  EmsyRefPts$Emsy_sp    <- Emsy_sp
  EmsyMSRefPts$EmsyMS_p   <- Emsy_p

  # Loop and get each stock/species ref pt
  for( p in 1:nP )
  {
    for( s in 1:nS )
    {
      tmp <- .calcEquil( f = Emsy_sp[s,p], obj = obj, type = "effort" )
      EmsyRefPts$yprEmsy_sp[s,p]    <- tmp$ypr_sp[s,p]
      EmsyRefPts$YeqEmsy_sp[s,p]    <- tmp$Yeq_sp[s,p]
      EmsyRefPts$BeqEmsy_sp[s,p]    <- tmp$Beq_sp[s,p]
      EmsyRefPts$expBeqEmsy_sp[s,p] <- tmp$expBeq_sp[s,p]
      EmsyRefPts$totBeqEmsy_sp[s,p] <- tmp$totBeq_sp[s,p]
      EmsyRefPts$NeqEmsy_sp[s,p]    <- tmp$Neq_sp[s,p]
      EmsyRefPts$ReqEmsy_sp[s,p]    <- tmp$Req_sp[s,p]
      EmsyRefPts$Umsy_sp[s,p]       <- tmp$Ueq_sp[s,p]

      tmpMS <- .calcEquil( f = Emsy_p[p], obj = obj, type = "effort" )
      EmsyMSRefPts$yprEmsy_sp[s,p]    <- tmpMS$ypr_sp[s,p]
      EmsyMSRefPts$YeqEmsy_sp[s,p]    <- tmpMS$Yeq_sp[s,p]
      EmsyMSRefPts$BeqEmsy_sp[s,p]    <- tmpMS$Beq_sp[s,p]
      EmsyMSRefPts$expBeqEmsy_sp[s,p] <- tmpMS$expBeq_sp[s,p]
      EmsyMSRefPts$totBeqEmsy_sp[s,p] <- tmpMS$totBeq_sp[s,p]
      EmsyMSRefPts$NeqEmsy_sp[s,p]    <- tmpMS$Neq_sp[s,p]
      EmsyMSRefPts$ReqEmsy_sp[s,p]    <- tmpMS$Req_sp[s,p]
      EmsyMSRefPts$Umsy_sp[s,p]       <- tmpMS$Ueq_sp[s,p]

    }
    # Now sum species eq yield for complex opt yield
    EmsyMSRefPts$YeqEmsy_p[p]         <- sum( EmsyMSRefPts$YeqEmsy_sp[,p] )
  }

  outlist <- list(  EmsyRefPts    = EmsyRefPts,
                    EmsyMSRefPts  = EmsyMSRefPts )
}

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
                        totBeqFmsy_sp = spMat,
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
      FmsyRefPts$totBeqFmsy_sp[s,p] <- tmp$totBeq_sp[s,p]
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
