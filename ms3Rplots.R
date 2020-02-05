# <><><><><><><><><><><><><><><><><><><><><><><><><><><>
# ms3Rplots.R
#
# Plots for ms3R.
#
# Author: SDN Johnson
# Date: October 8, 2019
#
# Last Update: October 10, 2019
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><>


# plotEffYieldCurves()
# Function for plotting effort based
# yield curves for each species and
# the complex in a stock area.
plotEffYieldCurves <- function( obj = blob )
{
  # First, pull reference points and curves
  rp          <- obj$rp[[1]]
  refCurves   <- rp$refCurves
  EmsyRefPts  <- rp$EmsyRefPts

  browser()
  

  # Now compute MSY for SS and MS curves
  # MSY is still the same for SS, just need the effort
  # that corresponds to it, which is Fmsy/qF
  EmsySS_sp <- rp$FmsyRefPts$Fmsy_sp / qF_sp
  MSYSS_sp  <- rp$FmsyRefPts$YeqFmsy_sp

  EmsyMS_p  <- rep(0,nP)
  MSYMS_p   <- rep(0,nP)
  MSYSS_sp  <- array(NA, dim = c(nS,nP))

  # Now we need to solve for the stat point
  # of the complex yield curve
  for( p in 1:nP )
  {
    # First, get complex Emsy, then
    # use splines to get species specific
    # yields
    compYieldSpline <- splinefun( x = Eseq, y = Yeq_pe[p,] )

    # Need to restrict to area where first deriv is non-singular
    # (i.e. where all yield curves are +ve), right now
    # uses Emax/3 as a placeholder

    # Solve for stat point

    EmsyMS_p[p] <- uniroot( f = compYieldSpline, 
                            interval = c(0, maxE / 3),
                            deriv = 1 )$root

    MSYMS_p[p] <- compYieldSpline(EmsyMS_p[p])

    # Now loop over species, fit spline
    # and get MSY for each species
    for( s in 1:nS )
    {
      specYieldSpline <- splinefun( x = Eseq, y = Yeq_spe[s,p,] )
      MSYSS_sp[s,p]   <- specYieldSpline(EmsyMS_p[p])
    }
  }



  specCols <- RColorBrewer::brewer.pal(n = nS, "Dark2")

  par( mfrow = c(3,1), mar = c(1,1,1,1), oma = c(3,3,1,1) )
  for( p in 1:nP )
  {
    plot( x = range(Eseq), y = c(0, max(Yeq_pe[p,],na.rm = T) ),
          type = "n", xlab = "", ylab = "", axes = F )
      axis(side = 2, las = 1)
      box()
      grid()

      for( s in 1:nS )
      {
        lines( x = Eseq, y = Yeq_spe[s,p,],
               col = specCols[s], lty = 1, lwd = 2 )
        # lines( x = Eseq, y = Beq_spe[s,p,],
        #        col = specCols[s], lty = 2, lwd = 2 )

      }
      mfg <- par("mfg")
      if( mfg[1] == mfg[3])
        axis(side = 1)

      lines(  x = Eseq, y = Yeq_pe[p,],
              col = "black", lty = 1, lwd = 3 )

      segments( x0 = EmsyMS_p[p], col = "grey40",
                y0 = 0, y1 = MSYMS_p[p], lty = 2  )

      segments( x0 = 0, x1 = EmsyMS_p[p], col = "grey40",
                y0 = c(MSYMS_p[p],MSYSS_sp[,p]), lty = 2  )

      if(  p == 1 )
        legend( x = "topright", bty = "n",
                col = c(specCols,"black"),
                lwd = c(2,2,2,3),
                legend = c(speciesNames,"Complex") )

  }

  mtext( outer = TRUE, side = 1, text = "Total Effort Index", line = 2 )
  mtext( outer = TRUE, side = 2, text = "Equilibrium Yield (kt)", line = 2 )


}

# plotEmpYieldCurves
# Function to plot median simuated yield 
# resulting from a grid of constant fishing 
# mortality rates - used for checking the
# reference point calculations
plotEmpYieldCurves <- function( sims = 1:11 )
{
  nSims <- length(sims)
  blobList <- vector(mode = "list", length = nSims)

  .loadSim(sims[1])

  nS <- blob$om$nS
  nP <- blob$om$nP
  nT <- blob$om$nT

  goodReps <- blob$goodReps


  # Arrays to hold empirical eqbm yields
  C_spk <- array( 0, dim = c(nS,nP,nSims) )
  B_spk <- array( 0, dim = c(nS,nP,nSims) )
  F_spk <- array( 0, dim = c(nS,nP,nSims) )

  for( x in sims )
  {
    .loadSim(x)
    
    C_spk[,,x] <- apply(X = blob$om$C_ispt[goodReps,,,nT,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
    B_spk[,,x] <- apply(X = blob$om$SB_ispt[goodReps,,,nT,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
    F_spk[,,x] <- apply(X = blob$om$F_ispft[goodReps,,,2,nT,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )

    # Clean up
    gc()
  }

  # Pull ref curves from RP object
  refCurves       <- blob$rp[[1]]$refCurves
  Fvec            <- refCurves$F
  BeqRefCurve_spf <- refCurves$Beq_spf
  YeqRefCurve_spf <- refCurves$Yeq_spf


  Fmsy_sp <- array(0, dim = c(nS,nP))
  Bmsy_sp <- array(0, dim = c(nS,nP))
  MSY_sp <- array(0, dim = c(nS,nP))

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(5,5,3,3) )

  for( s in 1:nS )
    for( p in 1:nP )
    {

      plot( x = c(0,max(F_spk[s,p,])), y = c(0,max(B_spk[s,p,])),
            type = "n", axes = FALSE )
        axis( side = 1 )
        axis( side = 2, las = 1 )
        grid()
        box()

        F <- F_spk[s,p,]
        C <- C_spk[s,p,]
        B <- B_spk[s,p,]

        actualOrder <- order(F)

        C <- C[actualOrder]
        B <- B[actualOrder]
        F <- F[actualOrder]

        CFspline <- splinefun(x = F, y = C)
        BFspline <- splinefun(x = F, y = B)
        empFmsy  <- uniroot(  interval = range(F),
                              f = CFspline,
                              deriv = 1 )$root
        empMSY  <- CFspline(empFmsy)
        empBmsy <- BFspline(empFmsy)

        empFmsy <- round(empFmsy,2)
        empBmsy <- round(empBmsy,2)
        empMSY  <- round(empMSY,2)

        lines( x = F, y = C,
                col = "steelblue", lwd = 2, lty = 1 )
        lines( x = F, y = B,
                col = "black", lwd = 2, lty = 1 )        

        lines( x = Fvec, y = YeqRefCurve_spf[s,p,],
                col = "salmon", lty = 2 )

        lines( x = Fvec, y = BeqRefCurve_spf[s,p,],
                col = "black", lty = 2 )


        legend( "topright",
                legend = c( paste("Fmsy = ", empFmsy, sep = "" ),
                            paste(" MSY = ", empMSY, sep = ""),
                            paste("Bmsy = ", empBmsy, sep = "") ) )


        # Plot some guidelines
        segments(  x0 = empFmsy, x1 = empFmsy,
                    y0 = 0, y1 = empBmsy,
                    lty = 2, col = "red" )
        segments(  x0 = 0, x1 = empFmsy,
                    y0 = empMSY, y1 = empMSY,
                    lty = 2, col = "red" )
        segments(  x0 = 0, x1 = empFmsy,
                    y0 = empBmsy, y1 = empBmsy,
                    lty = 2, col = "red" )
        
        Fmsy_sp[s,p] <- empFmsy
        Bmsy_sp[s,p] <- empBmsy
        MSY_sp[s,p]  <- empMSY
    }

  mtext( side = 1, text = "Fishing mortality (/yr)", outer = T, line = 2)
  mtext( side = 2, text = "Eqbm biomass and catch (kt)", outer = T, line = 2 )



  out <- list(  Fmsy_sp = Fmsy_sp,
                Bmsy_sp = Bmsy_sp,
                MSY_sp  = MSY_sp )

  out
} # END plotEmpYieldCurves()

# compareRefPts()
# Compares OM and conditioning fit reference points
# for a sanity check
compareYieldCurves <- function( obj )
{
  refPtsOM <- obj$rp[[1]]
  refPtsAM <- obj$ctlList$opMod$histRpt$refPts

  # Pull model dimensions
  nS <- obj$om$nS
  nP <- obj$om$nP

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames

  # Now, loop and plot reference points
  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(5,5,3,3) )

  for( s in 1:nS )
    for( p in 1:nP )
    {
      refCurvesOM <- refPtsOM$refCurves
      refCurvesAM <- refPtsAM$refCurves
      plot( x = range( refCurvesOM$F ),
            y = c(0, max(refCurvesOM$Yeq_spf[s,p,], refCurvesOM$Beq_spf[s,p,]  )),
            type = "n", axes = FALSE )
        axis( side = 1 )
        axis( side = 2, las = 1 )

        lines(  x = refCurvesOM$F, y = refCurvesOM$Yeq_spf[s,p,],
                col = "black", lwd = 3 )
        lines(  x = refCurvesOM$F, y = refCurvesOM$Beq_spf[s,p,],
                col = "steelblue", lwd = 3 )

        lines(  x = refCurvesAM$F, y = refCurvesAM$Yeq_spf[s,p,],
                col = "black", lty = 2, lwd = 3 )
        lines(  x = refCurvesAM$F, y = refCurvesAM$Beq_spf[s,p,],
                col = "black", lty = 2, lwd = 3 )


        points( x = refPtsOM$FmsyRefPts$Fmsy_sp[s,p],
                y = refPtsOM$FmsyRefPts$YeqFmsy_sp[s,p], 
                col = "red", pch = 16, cex = 1.5 )
        points( x = refPtsAM$FmsyRefPts$Fmsy_sp[s,p],
                y = refPtsAM$FmsyRefPts$YeqFmsy_sp[s,p], 
                col = "steelblue", pch = 21, cex = 1.5 )

    }
} 


# plotCvsB()
# Fishing mortality as a function of biomass,
# used for showing how well an MP meets the 
# proposed HCR, and for determining the HCR
# implied by the omniscient manager sim
plotCvsB <- function( obj )
{
  goodReps <- obj$goodReps

  SB_ispt   <- obj$om$SB_ispt[goodReps,,,,drop = FALSE]
  C_ispt    <- obj$om$C_ispt[goodReps,,,,drop = FALSE]
  F_ispt    <- obj$om$F_ispft[goodReps,,,2,]
  vB_ispt   <- obj$om$vB_ispft[goodReps,,,2,]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nReps   <- sum(goodReps)

  # Now, let's start plotting. We can fix it later
  speciesNames  <- (obj$om$speciesNames)
  stockNames    <- (obj$om$stockNames)


  projYrs <- tMP:nT

  dep_ispt <- SB_ispt
  Bmsy_sp  <- obj$rp[[1]]$FmsyRefPts$BeqFmsy_sp
  Fmsy_sp  <- obj$rp[[1]]$FmsyRefPts$Fmsy_sp
  MSY_sp   <- obj$rp[[1]]$FmsyRefPts$YeqFmsy_sp

  ctlList <- obj$ctlList

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(5,5,3,3) )

  U_ispt <- C_ispt / vB_ispt

  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      # Convert biomass to Bmsy depletion
      dep_ispt[,s,p,] <- SB_ispt[,s,p,] / Bmsy_sp[s,p]

      # Depletion and F
      maxDep <- max(dep_ispt[,s,p,projYrs], na.rm = T)
      maxC   <- max(C_ispt[,s,p,projYrs]/MSY_sp[s,p], na.rm = T)

      # Plot window
      plot( x = c(0,3), y = c(0,2),
            type = "n", axes = F)
        mfg <- par("mfg")
        # axes
        axis( side = 1 )

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 0.2, cex = 1.5 )

        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text( x = corners[2]+0.5, 
                y = mean(corners[3:4]), 
                labels = stockNames[p], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }
        
        axis( side = 2, las = 1 )
        box()

        ptCol <- scales::alpha("grey70", alpha = .3)

        points( x = dep_ispt[,s,p,projYrs], y = C_ispt[,s,p,projYrs]/MSY_sp[s,p],
                col = ptCol, pch = 1 ) 

        for( i in 1:nReps )
        { 
          # Plot a smoother for each replicate
          lineCol <- scales::alpha("red", alpha = .3)
          smoother <- loess.smooth( x = dep_ispt[i,s,p,projYrs], 
                                    y = C_ispt[i,s,p,projYrs]/MSY_sp[s,p] )
          lines( smoother, col = "grey30", lwd = .8 )
          flB <- dep_ispt[i,s,p,projYrs[c(1,length(projYrs))]]
          flC <- C_ispt[i,s,p,projYrs[c(1,length(projYrs))]]/MSY_sp[s,p]
          points( x = flB,
                  y = flC,
                  col = c("blue","red"), cex = .5,
                  pch = 16 )
        }


        abline( v = 1, lty = 2, lwd = .8)
        abline( h = 1, lty = 2, lwd = .8)

    }
  }

  mtext( side = 1, text = expression(B/B[MSY]), outer = TRUE, line = 2, font = 2)
  mtext( side = 2, text = expression(C/MSY), outer = TRUE, line = 2, font = 2)

}


# plotFvsB()
# Fishing mortality as a function of biomass,
# used for showing how well an MP meets the 
# proposed HCR, and for determining the HCR
# implied by the omniscient manager sim
plotFvsB <- function( obj )
{
  goodReps <- obj$goodReps

  SB_ispt   <- obj$om$SB_ispt[goodReps,,,,drop = FALSE]
  C_ispt    <- obj$om$C_ispt[goodReps,,,,drop = FALSE]
  F_ispt    <- obj$om$F_ispft[goodReps,,,2,]
  vB_ispt   <- obj$om$vB_ispft[goodReps,,,2,]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nReps   <- sum(goodReps)

  # Now, let's start plotting. We can fix it later
  speciesNames  <- (obj$om$speciesNames)
  stockNames    <- (obj$om$stockNames)


  projYrs <- tMP:nT

  dep_ispt <- SB_ispt
  Bmsy_sp  <- obj$rp[[1]]$FmsyRefPts$BeqFmsy_sp
  Fmsy_sp  <- obj$rp[[1]]$FmsyRefPts$Fmsy_sp

  ctlList <- obj$ctlList

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,2.5),
        oma = c(5,5,3,3) )

  U_ispt <- C_ispt / vB_ispt

  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      # Convert biomass to Bmsy depletion
      dep_ispt[,s,p,] <- SB_ispt[,s,p,] / Bmsy_sp[s,p]

      # Depletion and F
      maxDep <- max(dep_ispt[,s,p,projYrs], na.rm = T)
      maxF   <- max(F_ispt[,s,p,projYrs]/Fmsy_sp[s,p], na.rm = T)

      # Plot window
      plot( x = c(0,3), y = c(0,2),
            type = "n", axes = F)
        mfg <- par("mfg")
        # axes
        axis(side =1)

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 0.2, cex = 1.5 )

        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text( x = corners[2]+0.5, 
                y = mean(corners[3:4]), 
                labels = stockNames[p], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }
        
        axis( side = 2, las = 1 )
        box()

        ptCol <- scales::alpha("grey70", alpha = .3)

        points( x = dep_ispt[,s,p,projYrs], y = F_ispt[,s,p,projYrs]/Fmsy_sp[s,p],
                col = ptCol, pch = 1 ) 

        for( i in 1:nReps )
        { 
          # Plot a smoother for each replicate
          lineCol <- scales::alpha("red", alpha = .3)
          smoother <- loess.smooth( x = dep_ispt[i,s,p,projYrs], 
                                    y = F_ispt[i,s,p,projYrs]/Fmsy_sp[s,p] )
          lines( smoother, col = "grey30", lwd = .8 )
          flB <- dep_ispt[i,s,p,projYrs[c(1,length(projYrs))]]
          flF <- F_ispt[i,s,p,projYrs[c(1,length(projYrs))]]/Fmsy_sp[s,p]
          points( x = flB,
                  y = flF,
                  col = c("blue","red"), cex = .5,
                  pch = 16 )
        }


        abline( v = 1, lty = 2, lwd = .8)
        abline( h = 1, lty = 2, lwd = .8)

    }
  }

  mtext( side = 1, text = expression(B/B[MSY]), outer = TRUE, line = 2, font = 2)
  mtext( side = 2, text = expression(F/F[MSY]), outer = TRUE, line = 2, font = 2)

}

# plotBatchPerf_sp()
# Plots multipanels (faceted by species/stock)
# of performance statistics on the y axis, with respect
# to the OM grid on the x axis
# CAUTION: currently only works for numeric xAxis and yAxis
plotBatchPerf_sp <- function( batchFolder = "fourthBatch",
                              xAxis = "projObsErrMult",
                              yAxis = "PBtGt.8Bmsy",
                              yRangeIn = NULL )
{
  # First load full stats table from the batch folder
  statsFile <- here::here(  "Outputs",batchFolder,"statistics",
                            "fullStatTable.csv")
  statTable <- read.csv(statsFile, header = TRUE, stringsAsFactors = FALSE )

  # Now, let's start plotting. We can fix it later
  speciesNames  <- unique(statTable$species)
  stockNames    <- unique(statTable$stock)

  scenarios <- unique(statTable$scenario)
  mps       <- unique(statTable$mp)

  nS <- length(speciesNames)
  nP <- length(stockNames)

  nScen <- length(scenarios)
  nMPs  <- length(mps)



  xLabs <- unique(statTable[,xAxis])
  xLabs <- xLabs[order(xLabs)]
  xMax  <- max(xLabs)
  xMin  <- min(xLabs)
  
  cols <- RColorBrewer::brewer.pal(n = nMPs, "Dark2")


  par(  mfcol = c( nP, nS ), 
        mar = c( 0,2,0,2 ), 
        oma = c(5,7,3,4) )

  yRange <- yRangeIn

  mpJitter <- seq( from = -.3, to = .3, length = nMPs )

  xDiff <- mean(diff(xLabs))
  mpJitter <- mpJitter * xDiff


  for( s in 1:nS )
    for( p in 1:nP )
    {
      subTable <- statTable %>%
                  filter( species == speciesNames[s],
                          stock == stockNames[p] )
      
      if(is.null(yRangeIn))            
        yRange <- range(subTable[,yAxis])      
      
      
      plot( x = c(-.5 + xMin,xMax + .5), y = yRange,
            type = "n", axes = FALSE )
        mfg <- par("mfg")
        # axes
        if( mfg[1] == mfg[3])
          axis( side = 1, at = xLabs )

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
        
        axis( side = 2, las = 1 )
        
        box()
        # Add vertical grid lines to split OMs
        grid()

        # Now, plot the stuff
        for( m in 1:nMPs )
        {
          mpTable <- subTable %>% filter( mp == mps[m] )
          mpTable <- mpTable[order(mpTable[,xAxis]),]

          points( x = mpTable[,xAxis] + mpJitter[m],
                  y = mpTable[,yAxis],
                  col = cols[m], pch = 14 + m )
          lines(  x = mpTable[,xAxis] + mpJitter[m],
                  y = mpTable[,yAxis],
                  col = cols[m], lwd = .8 )
        }

        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text( x = corners[2]+0.2, 
                y = mean(corners[3:4]), 
                labels = stockNames[p], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }



    }
  legend( "topleft",
          legend = mps,
          col = cols,
          pch = 14 + 1:nMPs,
          lwd = 1 )

  

  mtext( side = 1, text = xAxis,
         outer = TRUE, cex = 1.5, 
         font = 2, line = 3 )

  mtext( side = 2, text = yAxis,
         outer = TRUE, cex = 1.5, font = 2, line = 2 )



}



# plotHCR()
# Plots a generic hockey-stick
# harvest control rule, with given
# lower, upper control points, and
# low/high Fs. Stock status is calculated
# as a proportion of Bmsy, and fishing
# mortality rate as a proportion of
# Fmsy.
plotHCR <- function(  LCP = .4,
                      UCP = .6,
                      lowF = .1,
                      highF = 1 )
{
  x <- seq(0,1.3, length.out = 100 )
  y <- rep(lowF, length(x))


  par( mar = c(4,5,1,1), oma = c(2,2,1,1) )

  plot( x = range(x), y = c(0,1.2), type = "n",
        las = 1,
        xlab = expression(B/B[MSY]),
        ylab = expression(F/F[MSY]),
        cex.lab = 1.5,
        cex.axis = 1.5 )
    segments( x0 = 0, x1 = LCP,
              y0 = lowF, y1 = lowF,
              col = "grey50", lwd = 3 )
    segments( x0 = LCP, x1 = UCP,
              y0 = lowF, y1 = highF,
              col = "grey50", lwd = 3 )
    segments( x0 = UCP, x1 = max(x),
              y0 = highF, y1 = highF,
              col = "grey50", lwd = 3 )
    abline( v = c(LCP, UCP),
            col = c("red","darkgreen"),
            lty = 2, lwd = 2 )
    abline( h = highF, lty = 2, col = "grey70" )


} # END plotHCR()

plotUtilityFunction <- function( )
{
  x <- seq(1e-1,1, length = 1000)

  y <- (5*x - 1)/4/x

  plot( x = range(x), y = range(y),
        xlab = expression(C[spt]/TAC[spt]),
        ylab = "Utility", type = "n", las = 1 )
    mtext( side = 3, text = "TAC utilisation utility", font = 2)
    lines( x = x, y = y, lwd = 3, col = "grey40" )
    abline(h = c(0,1), lty = 2, col = "black" )
    abline( v = .2, lty = 2, col = "red", lwd = .8 )

}


# plotConvStats()
plotConvStats <- function( obj = blob )
{
  goodReps      <- obj$goodReps

  # Pull max gradient value and hessian indicator
  maxGrad_itsp  <- obj$mp$assess$maxGrad_itsp[goodReps,,,,drop = FALSE]
  pdHess_itsp   <- obj$mp$assess$pdHess_itsp[goodReps,,,,drop = FALSE]

  nReps <- dim(maxGrad_itsp)[1]

  # Now we want to get the mean and SD
  # of these values over the replicates
  quantsMaxGrad_qtsp <- apply(  X = maxGrad_itsp, FUN = quantile,
                                MARGIN = 2:4, probs = c(0.05, 0.5, 0.95),
                                na.rm = T )

  propPDHess_tsp  <- apply( X = pdHess_itsp, FUN = mean,
                            MARGIN = 2:4, na.rm = T )

  pT <- obj$ctlList$opMod$pT
  nS <- obj$om$nS
  nP <- obj$om$nP

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames


  xJitter     <- seq( from = -.3, to = .3, length.out = nP )
  stockCols   <- RColorBrewer::brewer.pal( nP, "Dark2" )
  stockPts    <- seq( from = 21, by = 1, length.out = nP )
  rectWidth   <- .6 / nP

  par( mfcol = c(nS,2), mar = c(0,2,0,1), oma = c(3,3,3,2) )

  # First, plot the maxGrads, using 
  # different colours for each species, and 
  # pch for each stock... Jitter!
  for( s in 1:nS )
  {
    plot( x = c(1,pT), 
          y = range(quantsMaxGrad_qtsp, na.rm = T),
          type = "n",
          axes = FALSE )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3])
        axis( side = 1 )
      if(mfg[1] == 1 )
        mtext( side = 3, text = "Max Gradient Component")
      
      axis( side = 2, las = 1 )
      box()
      grid()

      for( p in 1:nP )
      {
        points( x = 1:pT + xJitter[p], y = quantsMaxGrad_qtsp[2,,s,p],
                col = stockCols[p], pch = stockPts[p], bg = stockCols[p] )
        segments( x0 = 1:pT + xJitter[p], x1 = 1:pT + xJitter[p],
                  y0 = quantsMaxGrad_qtsp[1,,s,p],
                  y1 = quantsMaxGrad_qtsp[3,,s,p],
                  col = stockCols[p], lty = 1 )
      }
      if( s == 1 )
        legend( "topright", bty = "n",
                legend = stockNames,
                col = stockCols,
                pch = stockPts,
                pt.bg = stockCols )
  }

  # now do proportion of PD hessians
  for( s in 1:nS )
  {
    plot( x = c(0,pT+1), y = c(0,1.3),
          axes = FALSE, type = "n" )
      # Axes
      if( mfg[1] == mfg[3])
        axis( side = 1 )
      axis( side = 2, las = 1 )
      mtext( side = 4, text = speciesNames[s], font = 2, line = 2)
      if(mfg[1] == 1 )
        mtext( side = 3, text = "Proportion of PD Hessians")

      box()
      abline( h = 1.0, lty = 2, lwd = 2, col = "grey40" )

      # Plot rectangles of PD Hessians
      for( p in 1:nP )
      {
        rect( xleft = 1:pT + xJitter[p] - rectWidth/2,
              xright = 1:pT + xJitter[p] + rectWidth/2,
              ybottom = 0,
              ytop = propPDHess_tsp[,s,p],
              col = stockCols[p] )
      }
      if( s == 1 )
        legend( "topright", bty = "n",
                legend = stockNames,
                col = stockCols,
                pch = 22, pt.bg = stockCols )

  }

}

plotTulipF <- function( obj = blob, nTrace = 3 )
{
  # Model dimensions
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF

  # Good replicates
  goodReps <- obj$goodReps

  # Pull reference points
  Fmsy_sp <- obj$rp[[1]]$FmsyRefPts$Fmsy_sp
    
  # Fishing mortality series
  F_ispft <- obj$om$F_ispft[goodReps,,,,]

  # Fishing mortality envelopes
  F_qspft <- apply(  X = F_ispft, FUN = quantile,
                    MARGIN = c(2,3,4,5), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  nReps   <- dim(F_ispft)[1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(5,5,3,3) )


  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(F_ispft[,s,p,,], na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      for( f in 1:nF )
      {
        polygon(  x = c(yrs, rev(yrs)),
                  y = c(F_qspft[1,s,p,f,], rev(F_qspft[3,s,p,f,])),
                  col = "grey65", border = NA )
        lines( x = yrs, y = F_qspft[2,s,p,f,], lwd = 3 )
        for( tIdx in traces )
          lines( x = yrs, y = F_ispft[tIdx,s,p,f,], lwd = .8 )
      }

      abline( v = yrs[tMP], col = "grey30", lty = 3 )
      abline( h = Fmsy_sp[s,p], lty = 2, lwd = 1, col = "red")
    }
  }

  mtext( side = 2, outer = TRUE, text = "Fishing mortality (/yr)",
          line = 2, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)
}

# TAC utilisation envelopes
plotTulipTACu <- function( obj = blob, nTrace = 3 )
{

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT

  goodReps <- obj$goodReps
  

  # Get catch for trawl fleet, in projections only
  C_ispft     <- obj$om$C_ispft[goodReps,,,2,tMP:nT,drop = FALSE]
  TAC_ispft   <- obj$mp$hcr$TAC_ispft[goodReps,,,2,tMP:nT,drop = FALSE]

  nReps   <- dim(TAC_ispft)[1]

  TACu_ispft <- C_ispft / TAC_ispft


  TACu_qspt <- apply( X = TACu_ispft, FUN = quantile,
                      MARGIN = c(2,3,5), probs = c(0.025, 0.5, 0.975),
                      na.rm = T )

  TACu_ispt <- array(0, dim = c(nReps,nS,nP,(nT - tMP + 1)))
  TACu_ispt[1:nReps,,,] <- TACu_ispft[1:nReps,,,1,]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)
  yrs <- yrs[tMP:nT]


  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(5,5,3,3) )


  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(TACu_qspt, na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(TACu_qspt[1,s,p,], rev(TACu_qspt[3,s,p,])),
                col = "grey65", border = NA )
      lines( x = yrs, y = TACu_qspt[2,s,p,], lwd = 3 )
      for( tIdx in traces )
        lines( x = yrs, y = TACu_ispt[tIdx,s,p,], lwd = .8 )

      abline( v = yrs[tMP], col = "grey30", lty = 3 )
    }
  }
  mtext( side = 2, outer = TRUE, text = expression(C[spt]/TAC[spt]),
          line = 2, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

}

# Biomass envelopes
plotTulipBt <- function(  obj = blob, nTrace = 3,
                          dep = FALSE,
                          ref = "B0",
                          Ct  = FALSE )
{
  goodReps <- obj$goodReps

  SB_ispt   <- obj$om$SB_ispt[goodReps,,,,drop = FALSE]
  C_ispt    <- obj$om$C_ispt[goodReps,,,,drop = FALSE]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nReps   <- dim(SB_ispt)[1]

  # Get reference points
  B0_sp     <- obj$ctlList$opMod$histRpt$B0_sp
  BmsySS_sp <- obj$rp[[1]]$FmsyRefPts$BeqFmsy_sp
  BmsyMS_sp <- obj$rp[[1]]$EmsyRefPts$BeqEmsyMS_sp

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if( dep )
  {
    if( ref == "B0" )
    {
      for( s in 1:nS )
        for( p in 1:nP )
        {
          SB_ispt[,s,p,] <- SB_ispt[,s,p,] / B0_sp[s,p]
          C_ispt[,s,p,]  <- C_ispt[,s,p,] / B0_sp[s,p]

        }
      BmsySS_sp <- BmsySS_sp / B0_sp
      BmsyMS_sp <- BmsyMS_sp / B0_sp
      B0_sp     <- B0_sp / B0_sp


    }

    if( ref == "Bmsy" )
    {
      for( s in 1:nS )
        for( p in 1:nP )
        {
          SB_ispt[,s,p,] <- SB_ispt[,s,p,] / Bmsy_sp[s,p]
          C_ispt[,s,p,]  <- C_ispt[,s,p,] / Bmsy_sp[s,p]
        }

      B0_sp     <- B0_sp / Bmsy_sp
      BmsySS_sp <- BmsySS_sp / Bmsy_sp
      BmsyMS_sp <- BmsyMS_sp / Bmsy_sp
      
    }
  }

  if( !dep )
    yAxisLab <- "Biomass (kt)"

  if( dep )
  {
    if( ref == "Bmsy" )
      yAxisLab <- expression(B[t]/B[MSY])

    if( ref == "B0" )
      yAxisLab <- expression(B[t]/B[0])
  }

  # Now take quantiles
  SB_qspt <- apply( X = SB_ispt, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  C_qspt <- apply( X = C_ispt, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )


  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(SB_qspt[,s,p,], na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+3, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(SB_qspt[1,s,p,], rev(SB_qspt[3,s,p,])),
                col = "grey65", border = NA )
      lines( x = yrs, y = SB_qspt[2,s,p,], lwd = 3 )
      for( tIdx in traces )
        lines( x = yrs, y = SB_ispt[tIdx,s,p,], lwd = .8 )

      if( Ct )
      {
        rect( xleft = yrs - .3, xright = yrs + .3,
              ybottom = 0, ytop = C_qspt[2,s,p,], col = "grey65",
              border = NA )
        segments( x0 = yrs[tMP:nT], x1 = yrs[tMP:nT],
                  y0 = C_qspt[1,s,p,tMP:nT], y1 = C_qspt[3,s,p,tMP:nT],
                  col = "black" )
      }

      abline( v = yrs[tMP], col = "grey30", lty = 3 )
      abline( h = B0_sp[s,p], lty = 2, col = "grey50", lwd = 2  )
      abline( h = BmsySS_sp[s,p], lty = 3, col = "darkgreen", lwd = 2)
      abline( h = BmsyMS_sp[s,p], lty = 3, col = "steelblue", lwd = 2)

      if( mfg[1] == 1 & mfg[2] == 1 )
        legend( x = "bottomleft", bty = "n",
                legend = c( "Median Spawning Biomass", 
                            "Central 95%",
                            "Replicate Traces",
                            "Unfished",
                            expression(B[MSY,MS]),
                            expression(B[MSY,SS])),
                col = c(  "black", "grey65", "black",
                          "grey50", "darkgreen","steelblue" ),
                pch = c(NA,15, NA, NA, NA,NA),
                lty = c(1, NA, 1, 2, 3, 3),
                lwd = c(3, NA, .8, 2, 2, 2 ) )
    }
  }
  mtext( side = 2, outer = TRUE, text = yAxisLab,
          line = 1.5, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)
}

# plotTulipEffort_p()
# Effort over time gridded
# by stock area - envelopes
plotTulipEffort_p <- function( obj = blob, nTrace = 3 )
{
  goodReps <- obj$goodReps

  E_ipft <- obj$om$E_ipft[goodReps,,,,drop = FALSE ]

  E_ipft[E_ipft == 0] <- NA

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(E_ipft)[1]

  E_qpft <- apply(  X = E_ipft, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  E_qpft[E_qpft == 0] <- NA

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )


  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  fleetCols <- RColorBrewer::brewer.pal( nF, "Dark2" )

  par(  mfcol = c(nP,1), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,4,3,3) )

  for( p in 1:nP )
  {
    plot( x = range(yrs),
          y = c(0,max(E_qpft[,p,,],na.rm = T) ),
          type = "n", axes = F )
      mfg <- par("mfg")
      # Plot axes and facet labels
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      for( f in 1:nF )
      {
        if(f == 2 )
        {
          polygon( x = c(yrs,rev(yrs)), y = c(E_qpft[1,p,f,],rev(E_qpft[3,p,f,])), 
                  col = scales::alpha(fleetCols[f], alpha = .3), border = NA )
          for( tIdx in traces )
            lines( x = yrs, y = E_ipft[tIdx,p,f,], col = fleetCols[f], lwd = .8 )
        }
        lines( x = yrs, y = E_qpft[2,p,f,], col = fleetCols[f], lwd = 3)
        
      }
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
  }
  mtext( side = 2, outer = TRUE, text = "Trawl Effort (fishing hours?)",
          line = 2, font = 2)
} # END plotTulipEffort_p()

# plotEffort_p()
# Effort over time gridded
# by stock area
plotEffort_p <- function( obj = blob,
                          iRep = 1 )
{
  E_pft <- obj$om$E_ipft[iRep,,,]

  E_pft[E_pft == 0] <- NA

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  fleetCols <- RColorBrewer::brewer.pal( nF, "Dark2" )

  par(  mfcol = c(nP,1), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,4,3,3) )

  for( p in 1:nP )
  {
    plot( x = range(yrs),
          y = c(0,max(E_pft[p,,],na.rm = T) ),
          type = "n", axes = F )
      mfg <- par("mfg")
      # Plot axes and facet labels
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      for( f in 1:nF )
        lines( x = yrs, y = E_pft[p,f,], col = fleetCols[f], lwd = 3 )
        
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
  }
  mtext( side = 2, outer = TRUE, text = "Trawl Effort (fishing hours?)",
          line = 2, font = 2)
} # END plotEffort_p()

# plotCtTACt_sp()
# Comparison plot of TAC and realized catch for each 
# species/stock
plotCtTACt_sp <- function(  obj = blob,
                            iRep = 1,
                            fleets = 1:2 )
{
  C_spft   <- obj$om$C_ispft[iRep,,,,]
  TAC_spft <- obj$mp$hcr$TAC_ispft[iRep,,,,]

  C_spft[C_spft == 0] <- NA

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  fleetCols <- RColorBrewer::brewer.pal( nF, "Dark2" )


  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )

  for(s in 1:nS)
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(C_spft[s,p,,], TAC_spft[s,p,,], na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      for( f in fleets )
      {
        rect( xleft = yrs - .3,
              xright = yrs + .3,
              ybottom = 0,
              ytop = TAC_spft[s,p,f,],
              border = NA, col = fleetCols[f] )
        
        lines(  x = yrs[], y = C_spft[s,p,f,],
                lty = 2, lwd = 2, col = "grey50" )
      }
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }
    mtext( outer = TRUE, side = 2, text = "Catch and TAC (kt)" )
} # END plotCtTACt_sp()


# plotBtCt_sp()
# Biomass and catch plots
# by species/stock for the historical
# and projection period
plotBtCt_sp <- function(  obj = blob,
                          iRep = 1 )
{
  SB_spt  <- obj$om$SB_ispt[iRep,,,]
  B_spt   <- obj$om$B_ispt[iRep,,,]
  C_spt   <- obj$om$C_ispt[iRep,,,]
  TAC_spt <- obj$mp$hcr$TAC_ispt[iRep,,,]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )

  for(s in 1:nS)
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(B_spt[s,p,]) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      rect( xleft = yrs - .3, xright = yrs + .3,
            ybottom = 0, ytop = C_spt[s,p,], col = "grey70", border = NA )
      segments( x0 = yrs, x1 = yrs,
                y0 = 0, y1 = TAC_spt[s,p,], col = "black" )
      lines( x = yrs, y = SB_spt[s,p,], col = "red", lwd = 3 )
      lines( x = yrs, y = B_spt[s,p,], col = "black", lwd = 2, lty = 2 )
      

      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )

    }


  mtext( side =1, outer = TRUE, text = "Year", line = 2)
  mtext( side =2, outer = TRUE, text = "Stock biomass and catch (kt)", line = 1.5)

}


# plotRetroSB()
# Retrospective plots of AM fits for a given replicate
plotRetroSB <- function( obj = blob, iRep = 1 )
{
  # Get biomass arrays
  SB_spt        <- obj$om$SB_ispt[iRep,,,]
  VB_spt        <- obj$om$vB_ispft[iRep,,,2,]
  totB_spt      <- obj$om$B_ispt[iRep,,,]
  retroSB_tspt  <- obj$mp$assess$retroSB_itspt[iRep,,,,]

  retroSB_tspt[retroSB_tspt < 0] <- NA

  # Get proportion of TACs for splitting aggregated biomass
  propTAC_spt   <- obj$mp$hcr$propTAC_ispt[iRep,,,]


  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT

  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fYear         <- obj$ctlList$opMod$fYear
  pT            <- obj$ctlList$opMod$pT

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )
  for(s in 1:nS)
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,1.2*max(totB_spt[s,p,],VB_spt[s,p,],SB_spt[s,p,],na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+1.5, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      lines( x = yrs, y = SB_spt[s,p,], col = "red", lwd = 3 )
      lines( x = yrs, y = VB_spt[s,p,], col = "grey40", lwd = 2, lty = 3 )
      lines( x = yrs, y = totB_spt[s,p,], col = "black", lwd = 2 )
      for( tt in 1:pT )
      {
        propTAC <- propTAC_spt[s,p,tMP + tt - 1]
        lines( x = yrs, y = propTAC * retroSB_tspt[tt,s,p,], col = "grey60", lwd = 1 )
      }
  
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }

}

# plotRetroSBagg()
# Retrospective plots of AM fits for a given replicate
# with biomasses aggregated to match the scale of the AM
plotRetroSBagg <- function( obj = blob, iRep = 1, Ct = TRUE )
{
  # Get biomass arrays
  SB_spt        <- obj$om$SB_ispt[iRep,,,]
  VB_spt        <- obj$om$vB_ispft[iRep,,,2,]
  totB_spt      <- obj$om$B_ispt[iRep,,,]
  retroSB_tspt  <- obj$mp$assess$retroSB_itspt[iRep,,,,]

  ctlList <- obj$ctlList

  retroSB_tspt[retroSB_tspt < 0] <- NA

  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT


  C_spft   <- obj$om$C_ispft[iRep,,,,]
  TAC_spft <- obj$mp$hcr$TAC_ispft[iRep,,,,]

  C_spt    <- apply(X = C_spft, FUN = sum, MARGIN = c(1,2,4))
  TAC_spt  <- apply(X = TAC_spft, FUN = sum, MARGIN = c(1,2,4))

  # Aggregate OM biomasses to match AM biomass
  if( ctlList$mp$assess$spCoastwide )
  {
    newSB_spt        <- apply( X = SB_spt, FUN = sum, MARGIN = c(1,3), na.rm = T )
    newVB_spt        <- apply( X = VB_spt, FUN = sum, MARGIN = c(1,3), na.rm = T )
    newtotB_spt      <- apply( X = totB_spt, FUN = sum, MARGIN = c(1,3), na.rm = T )

    SB_spt    <- SB_spt[,1,,drop = FALSE]
    SB_spt[,1,] <- newSB_spt
    VB_spt    <- VB_spt[,1,,drop = FALSE]
    VB_spt[,1,] <- newVB_spt
    totB_spt  <- totB_spt[,1,,drop = FALSE]
    totB_spt[,1,] <- newtotB_spt

    sumC_spt          <- C_spft[,1,1,,drop = FALSE]
    sumC_spt[,1,1,]   <- apply(X = C_spft, FUN = sum, MARGIN = c(1,4))
    newC_spt          <- array(NA, dim = c(nS,1,nT))
    newC_spt[,1,]     <- sumC_spt[,1,1,]
    C_spt             <- newC_spt

    sumTAC_spt        <- C_spft[,1,,,drop = FALSE]
    sumTAC_spt[,1,1,] <- apply(X = TAC_spft, FUN = sum, MARGIN = c(1,4))
    newTAC_spt        <- array(NA, dim = c(nS,1,nT))
    newTAC_spt[,1,]   <- sumTAC_spt[,1,1,]
    TAC_spt           <- newTAC_spt        

  }

  if( ctlList$mp$assess$spDataPooled )
  {
    newSB_spt        <- apply( X = SB_spt, FUN = sum, MARGIN = c(2,3), na.rm = T )
    newVB_spt        <- apply( X = VB_spt, FUN = sum, MARGIN = c(2,3), na.rm = T )
    newtotB_spt      <- apply( X = totB_spt, FUN = sum, MARGIN = c(2,3), na.rm = T )

    SB_spt    <- SB_spt[1,,,drop = FALSE]
    SB_spt[1,,] <- newSB_spt
    VB_spt    <- VB_spt[1,,,drop = FALSE]
    VB_spt[1,,] <- newVB_spt
    totB_spt  <- totB_spt[1,,,drop = FALSE]
    totB_spt[1,,] <- newtotB_spt

    sumC_spt          <- C_spft[1,,1,,drop = FALSE]
    sumC_spt[1,,1,]   <- apply(X = C_spft, FUN = sum, MARGIN = c(2,4))
    newC_spt          <- array(NA, dim = c(1,nP,nT))
    newC_spt[1,,]     <- sumC_spt[1,,1,]
    C_spt             <- newC_spt
    

    sumTAC_spt        <- TAC_spft[1,,1,,drop = FALSE]
    sumTAC_spt[1,,,]  <- apply(X = TAC_spft, FUN = sum, MARGIN = c(2,4))
    newTAC_spt        <- array(NA, dim = c(1,nP,nT))
    newTAC_spt[1,,]   <- sumTAC_spt[1,,1,]
    TAC_spt           <- newTAC_spt
  }

  SB_spt[SB_spt == 0]     <- NA
  VB_spt[VB_spt == 0]     <- NA
  totB_spt[totB_spt == 0] <- NA

  nSS     <- dim( SB_spt)[1]
  nPP     <- dim( SB_spt)[2]

  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fYear         <- obj$ctlList$opMod$fYear
  pT            <- obj$ctlList$opMod$pT

  if( nPP == 1 )
    stockNames <- "Coastwide"

  if( nSS == 1 )
    speciesNames <- "Data Pooled"

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  par(  mfcol = c(nPP,nSS), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )
  for(s in 1:nSS)
    for( p in 1:nPP )
    {
      plot( x = range(yrs),
            y = c(0,1.2*max(totB_spt[s,p,],VB_spt[s,p,],SB_spt[s,p,],na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+1.5, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      lines( x = yrs, y = SB_spt[s,p,], col = "red", lwd = 3 )
      lines( x = yrs, y = VB_spt[s,p,], col = "grey40", lwd = 2, lty = 3 )
      lines( x = yrs, y = totB_spt[s,p,], col = "black", lwd = 2 )
      for( tt in 1:pT )
      {
        lines( x = yrs, y = retroSB_tspt[tt,s,p,], col = "grey60", lwd = 1 )
      }
      if( Ct )
      {
        # Plot actual catch
        rect( xleft = yrs - .3, xright = yrs + .3, 
              ytop = C_spt[s,p,], ybottom = 0, col = "grey40",
              border = NA )
        # Plot a rectangle of TACs
        rect( xleft = yrs[tMP:nT] - .3, xright = yrs[tMP:nT] + .3, 
              ytop = TAC_spt[s,p,tMP:nT], ybottom = 0, col = NA,
              border = "black" )

      }
  
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }

}


.plotDiagCondition <- function( obj = blob, iRep = 1)
{
  repObj <- obj$ctlList$opMod$histRpt

  diagCondition( repObj, obj, iRep )


}

# diagCondition()
# Plots that help diagnose issues with conditioning
# between the fitted OM report, and the conditioned
# ms3R operating model.
diagCondition <- function(  repObj  = totRep,
                            ms3Obj  = test,
                            iRep    = 1 )
{
  par(mfrow = c(3,2), mar = c(1.5,1,1.5,1), oma = c(3,3,3,3) )

  # Biomass RE
  plotRE_spt( repObj = repObj, omObj = ms3Obj$om, 
              AMseries = "SB_spt", 
              OMseries = "SB_ispt", 
              iRep )
  mtext( side = 3, text = "SB_spt", line = 0, font = 2)

  # Recruitment
  plotRE_spt( repObj = repObj, 
              omObj = ms3Obj$om, 
              AMseries = "R_spt",
              OMseries = "R_ispt", 
              iRep )
  mtext( side = 3, text = "R_spt", line = 0, font = 2)  

  # Recruitment errors
  plotRE_spt( repObj = repObj, 
              omObj = ms3Obj$errors, 
              AMseries = "omegaR_spt",
              OMseries = "omegaR_ispt", 
              iRep )
  mtext( side = 3, text = expression(omega[R]), line = 0, font = 2)  

  # Catch
  plotRE_spft(  repObj = repObj, 
                omObj = ms3Obj$om, 
                AMseries = "C_spft",
                OMseries = "C_ispft", 
                iRep )
  mtext( side = 3, text = expression(C[spft]), line = 0, font = 2)    

  # F
  plotRE_spft(  repObj = repObj, 
                omObj = ms3Obj$om, 
                AMseries = "F_spft",
                OMseries = "F_ispft", 
                iRep )
  mtext( side = 3, text = expression(F[spft]), line = 0, font = 2) 

  # # Numbers at age
  # plotRE_axspt( repObj = repObj, omObj = ms3Obj$om, series = "N_iaxspt" )
  # mtext( side = 3, text = expression(N[axspt]), line = 0, font = 2)    

  # # Biomass at age
  # plotRE_axspt( repObj = repObj, omObj = ms3Obj$om, series = "N_axspt" )
  # mtext( side = 3, text = expression(B[axspt]), line = 0, font = 2)    
}

calcREdist <- function( est, true, marg, 
                        qProbs = c(0.025, 0.5, 0.975) )
{
  # Calculate REs
  re <- (est - true)/true

  # calculate distribution over margin
  reQuants <-  apply( X = re, FUN = quantile, 
                      MARGIN = marg, probs = qProbs,
                      na.rm = TRUE )

  reQuants[!is.finite(reQuants)] <- 0

  return(reQuants)

}

plotRE_spt <- function( repObj, omObj, 
                        AMseries = "SB_spt",
                        OMseries = "SB_ispt", iRep )
{
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT



  true_spt <- repObj[[AMseries]]
  est_spt  <- omObj[[OMseries]][iRep,,,]

  re_qspt  <- calcREdist( true = true_spt,
                          est  = est_spt[,,1:nT],
                          marg = c(1,2,3) )

  specCols <- RColorBrewer::brewer.pal(n = nS, "Dark2")
  stockLty <- 1:nP

  yRange <- range(re_qspt, na.rm = TRUE)
  yRange[2] <- max(.01,yRange[2])
  yRange[1] <- min(-.01,yRange[1])

  plot( x = c(1,nT), y = yRange,
        axes = FALSE, type = "n" )
    mfg <- par("mfg")
    if(mfg[1] == mfg[3])
      axis( side = 1 )
    axis( side = 2, las = 1 )
    box()
    grid()
    for( s in 1:nS )
      for( p in 1:nP )
      {
        polyY <- c(re_qspt[1,s,p,],rev(re_qspt[3,s,p,]))
        polygon( x = c(1:nT,nT:1), y = polyY,
                  col = scales::alpha(specCols[s],alpha = .3),
                  lty = stockLty[p] )
        lines( x = 1:nT, y = re_qspt[2,s,p,], 
                col = specCols[s], lty = stockLty[p] )

      }
    abline( h = 0, lty = 3, lwd = .8 )
}

plotRE_spft <- function( repObj, omObj, series = "C_spft", iRep )
{
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT

  true_spft <- repObj[[AMseries]]
  est_spft  <- omObj[[OMseries]][iRep,,,,]

  re_qspt  <- calcREdist( true = true_spft,
                            est  = est_spft[,,,1:nT],
                            marg = c(1,2,4) )

  specCols <- RColorBrewer::brewer.pal(n = nS, "Dark2")
  stockLty <- 1:nP

  yRange <- range(re_qspt, na.rm = TRUE)
  yRange[2] <- max(.01,yRange[2])
  yRange[1] <- min(-.01,yRange[1])

  plot( x = c(1,nT), y = yRange,
        axes = FALSE, type = "n" )
    mfg <- par("mfg")
    if(mfg[1] == mfg[3])
      axis( side = 1 )
    axis( side = 2, las = 1 )
    box()
    grid()

    for( s in 1:nS )
      for( p in 1:nP )
      {
        polyY <- c(re_qspt[1,s,p,],rev(re_qspt[3,s,p,]))
        polygon( x = c(1:nT,nT:1), y = polyY,
                  col = scales::alpha(specCols[s],alpha = .3),
                  lty = stockLty[p] )
        lines( x = 1:nT, y = re_qspt[2,s,p,], 
                col = specCols[s], lty = stockLty[p] )

      }
    abline( h = 0, lty = 3, lwd = .8 )
}

plotRE_axspt <- function( repObj, omObj, series = "N_axspt", iRep )
{
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT

  true_axspt <- repObj[[AMseries]]
  est_axspt  <- omObj[[OMseries]][iRep,,,,,]

  re_qspt  <- calcREdist( true = true_axspt,
                            est  = est_axspt[,,,,1:nT],
                            marg = c(3:5) )

  specCols <- RColorBrewer::brewer.pal(n = nS, "Dark2")
  stockLty <- 1:nP

  yRange <- range(re_qspt, na.rm = TRUE)
  yRange[2] <- max(.01,yRange[2])
  yRange[1] <- min(-.01,yRange[1])

  plot( x = c(1,nT), y = yRange,
        axes = FALSE, type = "n" )
    mfg <- par("mfg")
    if(mfg[1] == mfg[3])
      axis( side = 1 )
    axis( side = 2, las = 1 )
    box()
    grid()
    for( s in 1:nS )
      for( p in 1:nP )
      {
        polyY <- c(re_qspt[1,s,p,],rev(re_qspt[3,s,p,]))
        polygon( x = c(1:nT,nT:1), y = polyY,
                  col = scales::alpha(specCols[s],alpha = .3),
                  lty = stockLty[p] )
        lines( x = 1:nT, y = re_qspt[2,s,p,], 
                col = specCols[s], lty = stockLty[p] )

      }
    abline( h = 0, lty = 3, lwd = .8 )
}