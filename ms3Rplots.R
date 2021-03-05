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

plotFvsEffHist_p <- function( obj = blob, 
                              fIdx = 2,
                              iRep = 1 )
{
  # Get effort and fishing mortality
  E_pft     <- obj$om$E_ipft[iRep,,,1:61]
  F_spft    <- obj$om$F_ispft[iRep,,,,1:61]
  qF_spf    <- obj$om$qF_ispft[iRep,,,,62]

  maxEff      <- max(E_pft[,fIdx,], na.rm = T)
  rangeEff_p  <- apply( X = E_pft[,fIdx,], FUN = range, na.rm = T, MARGIN = 1)




  specCols    <- wes_palette("Rushmore1", 5, type = "discrete")[3:5]
  stockNames  <- obj$om$stockNames
  speciesNames<- obj$om$speciesNames

  nP <- obj$om$nP
  nS <- obj$om$nS


  par( mfrow = c(nP,1), oma = c(4,4,1,3), mar = c(.5,.5,.5,.5) )

  for( p in 1:nP )
  {
    plot( x = c(0,maxEff), y = c(0,max(max(rangeEff_p[,p]) * qF_spf[,,fIdx],F_spft[,p,fIdx,41:61])),
          type = "n", axes = FALSE, xlab = "", ylab = "")
    mfg <- par("mfg")
    axis(side = 2, las = 1)
    if(mfg[1] == mfg[3])
      axis(side = 1)
    grid()
    box()

    effSeq <- seq(from = rangeEff_p[1,p], 
                  to = rangeEff_p[2,p], length.out = 3)

    for( s in 1:nS )
    {
      points( x = E_pft[p,fIdx,41:61], y = F_spft[s,p,fIdx,41:61],
              pch = 20+s, cex = 2, col = specCols[s], lwd = 2 )
      lines( x = effSeq, y = qF_spf[s,p,fIdx]*effSeq,
              lty = s, col = specCols[s], lwd = 2)

    }
    if( p == 1 )
      legend( x = "topright", bty = "n",
              legend = speciesNames,
              pch = 20 + 1:nS,
              lty = 1:s,
              lwd = 2, pt.lwd = 2, cex = 2,
              col = specCols)
    rmtext( outer = TRUE, line = .1, txt = stockNames[p], font = 2, cex = 1.5)
  }
  mtext(line = 2, side = 1, outer = TRUE, text = "Commercial Trawl Effort (x1000 hours)")
  mtext(line = 2, side = 2, outer = TRUE, text = "Fishing Mortality")




} # END plotFvsEffHist_p()




# plotPriceCatchTS()
# Plots price and catch series for the inverse
# demand curve models.
plotPriceCatchTS <- function( priceModel = priceFlexModel,
                              i = 0.02,
                              baseYr = 2016,
                              firstYr = 2006,
                              specNames = c("Dover","English","Rock") )
{
  # Pull catch
  C_st <- priceModel$C_st

  P_st <- priceModel$P_st
  yrCols <- seq(  from = firstYr, 
                  by = 1, 
                  length.out = ncol(P_st) )
  dimnames(P_st)[[2]] <- yrCols

  for(s in 1:nrow(P_st) )
    P_st[s,] <- P_st[s,] * (1 + i)^(baseYr - as.numeric(dimnames(P_st)[[2]]))

  
  specJitter <- c(-.2,0,.2)
  specDens   <- c(15,15,-1)
  specAngles <- c(45,-45,0)

  # First, let's plot the price data
  par( oma = c(2,2,1,3))
  plot(x = c(firstYr,2016), y = c(0,1.2*max(P_st,C_st)), type = "n",
        xlab = "Year", ylab = "Unit Price ($/kg)", las = 1,
        yaxs = "i" )
    for( s in 1:3 )
    {
      rect( xleft = firstYr:2016 + specJitter[s] - .1,
            xright = firstYr:2016 + specJitter[s] + .1,
            ybottom = 0, ytop = C_st[s,], border = "grey50",
            fill = "grey50",
            col = "grey50", lwd = 2, density = specDens[s], angle = specAngles[s] )
    }
    for( s in 1:3 )
    {
      lines(  x = firstYr:2016, y = P_st[s,],
              col = "grey20", lwd = 3, lty = s )
    }
    axis( side = 4, las = 1)
    rmtext( outer = TRUE, txt = "Catch (kt)", line = .05 )
    legend( x = "topright", bty = "n",
            legend = paste(specNames,"catch"),
            col = "grey20",
            fill = c("grey50","grey50","grey50"),
            border = "grey50",
            density = c(15,15,NA),
            angle = c(45,-45,0) )
    legend( x = "topleft", bty = "n",
            legend = paste(specNames, "price"),
            col = "grey20",
            lty = 1:3, lwd = 3 )

} # END plotPriceCatchTS()

# plotCWeffort()
plotCWeffort <- function( obj,
                          maxE = NULL )
{
  rp <- obj$rp[[1]]

  # Bail out of function if info missing
  if( is.null(rp$EmeyRefPts$cwEconYieldCurves))
  {
    cat("Simulation object did not calculate CW econ yield equilibria.")
    return()
  }

  # Otherwise, grab info, follow other yield curve plots, but
  # this time there's an aggregate curve too
  
  opMod <- obj$ctlList$opMod

  refCurves         <- rp$refCurves
  EmsyRefPts        <- rp$EmsyRefPts
  EmsyMSRefPts      <- rp$EmsyMSRefPts
  FmsyRefPts        <- rp$FmsyRefPts
  EmeyRefPts        <- obj$rp[[1]]$EmeyRefPts
  cwEconYieldCurves <- EmeyRefPts$cwEconYieldCurves
  

  nT  <- obj$om$nT
  nP  <- obj$om$nP
  nS  <- obj$om$nS
  tMP <- obj$om$tMP

  # Pull qF
  qF_sp       <- blob$om$qF_ispft[1,,,2,nT]
  
  speciesNames <- blob$om$speciesNames
  stockNames   <- blob$om$stockNames

  # Now compute MSY for SS and MS curves
  # MSY is still the same for SS, just need the effort
  # that corresponds to it, which is Fmsy/qF
  Eseq      <- as.numeric(dimnames(EmeyRefPts$econRev_pe)[[2]])

  EmsyMS_p  <- EmsyMSRefPts$EmsyMS_p
  MSYMS_sp  <- EmsyMSRefPts$YeqEmsy_sp

  EmsySS_p  <- EmsyRefPts$Emsy_sp
  MSYSS_sp  <- EmsyRefPts$YeqEmsy_sp

  eff_p <- obj$om$E_ipft[1,,2,54]
  C_p   <- apply( X = obj$om$C_ispt[1,,,54], FUN = sum, MARGIN = 2 )

  # Now, we need to convert solveEff_p to the cost of
  # fuel, using the landings and fuel cost per kt
  effortPrice_p <- opMod$varCostPerKt * C_p / eff_p

  
  Rev_spe     <- cwEconYieldCurves$cwRev_spe
  Emey_p      <- cwEconYieldCurves$cwEmey_p
  MEY_p       <- cwEconYieldCurves$cwMEY_p
  Yeq_spe     <- cwEconYieldCurves$Yeq_spe
  Rent_pe     <- cwEconYieldCurves$cwRent_pe
  Rent_e      <- cwEconYieldCurves$cwRent_e
  effCost_pe  <- cwEconYieldCurves$cwEffCost_pe
  E_pe        <- array(NA, dim = dim(effCost_pe))
  Eff         <- cwEconYieldCurves$Eff

  for( p in 1:nP )
    E_pe[p,] <- effCost_pe[p,] / effortPrice_p[p]
  
  Eff         <- cwEconYieldCurves$Eff

  propE_pe    <- E_pe
  for( p in 1:nP )
    propE_pe[p,] <- E_pe[p,]/Eff

  propE_pe[,1] <- NA  # Remove 0/0s

  if( is.null(maxE) )
    maxE <- max(Eff)

  # Now, plot effort on X axis, and the area efforts
  # on the y-axis (or the proportions)
  plot( x = c(0,maxE), y = c(0,1),
        type = "n", xlab = "Total trawl effort (x1000 hrs)",
        ylab = "Proportion of total effort" )
    for(p in 1:nP )
      lines( x = Eff, y = propE_pe[p,],
              col = p )
} # END plotCWeffort

# plotCWeconYield()
# Plot function for economic yield curves
# assuming a coastwide demand curve but spatially
# heterogeneous fishing effort/costs of fishing/catch
plotCWeconYield <- function(  obj,
                              maxE = 120,
                              plotCplx = TRUE )
{
  rp <- obj$rp[[1]]

  # Bail out of function if info missing
  if( is.null(rp$EmeyRefPts$cwEconYieldCurves))
  {
    cat("Simulation object did not calculate CW econ yield equilibria.")
    return()
  }

  # Otherwise, grab info, follow other yield curve plots, but
  # this time there's an aggregate curve too
  
  opMod <- obj$ctlList$opMod

  refCurves         <- rp$refCurves
  EmsyRefPts        <- rp$EmsyRefPts
  EmsyMSRefPts      <- rp$EmsyMSRefPts
  FmsyRefPts        <- rp$FmsyRefPts
  EmeyRefPts        <- obj$rp[[1]]$EmeyRefPts
  cwEconYieldCurves <- EmeyRefPts$cwEconYieldCurves
  

  nT  <- obj$om$nT
  nP  <- obj$om$nP
  nS  <- obj$om$nS
  tMP <- obj$om$tMP

  # Pull qF
  qF_sp       <- blob$om$qF_ispft[1,,,2,nT]
  
  speciesNames <- blob$om$speciesNames
  stockNames   <- blob$om$stockNames

  # Now compute MSY for SS and MS curves
  # MSY is still the same for SS, just need the effort
  # that corresponds to it, which is Fmsy/qF
  Eseq      <- as.numeric(dimnames(EmeyRefPts$econRev_pe)[[2]])

  EmsyMS_p  <- EmsyMSRefPts$EmsyMS_p
  MSYMS_sp  <- EmsyMSRefPts$YeqEmsy_sp

  EmsySS_p  <- EmsyRefPts$Emsy_sp
  MSYSS_sp  <- EmsyRefPts$YeqEmsy_sp

  # Emey_p      <- EmeyRefPts$Emey_p
  # MEY_p       <- EmeyRefPts$MEY_p
  # Rev_pe      <- EmeyRefPts$econRev_pe
  # Rev_spe     <- EmeyRefPts$econRev_spe
  # econYeq_pe  <- EmeyRefPts$econYeq_pe
  # effCost_pe  <- EmeyRefPts$effCost_pe
  # Ymey_sp     <- EmeyRefPts$Ymey_sp

  eff_p <- obj$om$E_ipft[1,,2,54]
  C_p   <- apply( X = obj$om$C_ispt[1,,,54], FUN = sum, MARGIN = 2 )

  # Now, we need to convert solveEff_p to the cost of
  # fuel, using the landings and fuel cost per kt
  effortPrice_p <- opMod$varCostPerKt * C_p / eff_p

  Rev_spe     <- cwEconYieldCurves$cwRev_spe
  Emey_p      <- cwEconYieldCurves$cwEmey_p
  MEY_p       <- cwEconYieldCurves$cwMEY_p
  Yeq_spe     <- cwEconYieldCurves$Yeq_spe
  Rent_pe     <- cwEconYieldCurves$cwRent_pe
  Rent_e      <- cwEconYieldCurves$cwRent_e
  effCost_pe  <- cwEconYieldCurves$cwEffCost_pe
  E_pe        <- array(NA, dim = dim(effCost_pe))
  Eff         <- cwEconYieldCurves$Eff

  for( p in 1:nP )
    E_pe[p,] <- effCost_pe[p,] / effortPrice_p[p]

  Rev_pe    <- apply(X = Rev_spe, FUN = sum, MARGIN = c(2,3), na.rm = TRUE )
  Rev_se    <- apply(X = Rev_spe, FUN = sum, MARGIN = c(1,3), na.rm = TRUE )
  Rev_e     <- apply(X = Rev_spe, FUN = sum, MARGIN = c(3), na.rm = TRUE )
  effCost_e <- apply(X = effCost_pe, FUN = sum, MARGIN = 2 )


  # remove zeroes
  Rev_spe[Rev_spe == 0] <- NA
  Rev_spe[,,1] <- 0

  specCols <- wes_palette("Rushmore1", 5, type = "discrete")[3:5]

  if(is.null(maxE))
    maxE <- 10 * max(Emey_p[pIdx])

  par( mfrow = c(nP+1,1), mar = c(.1,1,.1,1), oma = c(3,4,1,1) )
  for( p in 1:nP )
  {
    plot( x = c(0,maxE), y = c(0, max(Rent_pe[p,],Rev_pe[p,],na.rm = T) ),
          type = "n", xlab = "", ylab = "", axes = F, xaxs="i" )
      axis(side = 2, las = 1)
      box()
      grid()

      for( s in 1:nS )
      {
        lines( x = E_pe[p,], y = Rev_spe[s,p,],
               col = specCols[s], lty = 1, lwd = 2 )
        # lines( x = Eseq, y = Beq_spe[s,p,],
        #        col = specCols[s], lty = 2, lwd = 2 )

      }
      mfg <- par("mfg")
      if( mfg[1] == mfg[3])
        axis(side = 1)

      rmtext( txt = stockNames[p], line = 0.02,
              font = 2, cex = 1.5, outer = TRUE )

      if( plotCplx )
      {
        lines(  x = E_pe[p,], y = Rev_pe[p,],
                col = "black", lty = 1, lwd = 3 )

        lines(  x = E_pe[p,], y = Rent_pe[p,],
                col = "black", lty = 2, lwd = 3 )
      

      lines( x = E_pe[p,], y = effCost_pe[p,], col = "red",
              lwd = 2 )

      abline( v = Emey_p[p], lty = 2, col = "steelblue")
      abline( v = EmsyMS_p[p], lty = 2, col = "grey40")
      }



      # segments( x0 = EmsyMS_p[p], col = "grey40",
      #           y0 = 0, y1 = MSYMS_p[p], lty = 2  )

      # segments( x0 = 0, x1 = EmsyMS_p[p], col = "grey40",
      #           y0 = c(MSYMS_p[p],MSYMS_sp[,p]), lty = 2  )

      if(  p == 1 )
        legend( x = "topright", bty = "n",
                col = c(specCols,"black","red","black"),
                lty = c(1,1,1,1,1,2),
                lwd = c(2,2,2,2,2,2),
                legend = c( paste(speciesNames," Revenue",sep = ""),
                            "Complex Revenue",
                            "Variable Costs",
                            "Total Profits") )

  }

  plot( x = c(0,maxE), y = c(0, max(Rent_e,Rev_e,na.rm = T) ),
          type = "n", xlab = "", ylab = "", axes = F, xaxs="i" )
      axis(side = 2, las = 1)
      box()
      grid()

      for( s in 1:nS )
      {
        lines( x = Eff, y = Rev_se[s,],
               col = specCols[s], lty = 1, lwd = 2 )
      }

      mfg <- par("mfg")
      if( mfg[1] == mfg[3])
        axis(side = 1)

      rmtext( txt = "Coast-wide", line = 0.02,
              font = 2, cex = 1.5, outer = TRUE )

      if( plotCplx )
      {
        lines(  x = Eff, y = Rev_e,
                col = "black", lty = 1, lwd = 3 )

        lines(  x = Eff, y = Rent_e,
                col = "black", lty = 2, lwd = 3 )
      

        lines( x = Eff, y = effCost_e, col = "red",
                lwd = 2 )

        abline( v = sum(Emey_p), lty = 2, col = "steelblue")
        abline( v = sum(EmsyMS_p), lty = 2, col = "grey40")
      }
      



  mtext( outer = TRUE, side = 1, text = "Commercial Trawl Effort (1000 hrs)", line = 2 )
  mtext( outer = TRUE, side = 2, text = "Cost, Revenue, and Profit ($ x 1e6)", line = 2 )



}

# plotChokeHCRs()
# Plot function for making HCRs with an example choke
# effect
plotChokeHCRs <- function(  LCP = .1, UCP = .6,
                            refHR = .1,
                            chokeScale = .2 )
{
  BBmsy <- seq(from = 0, to = 1, length.out = 200 )

  targHR <- rep(0,length.out = 100)

  for( bIdx in 1:length(BBmsy) )
  {
    b <- BBmsy[bIdx]
    
    if( b >= LCP & b <= UCP )
      targHR[bIdx] <- refHR * (b - LCP)/(UCP - LCP)

    if( b >= UCP )
      targHR[bIdx] <- refHR
  }

  chokeHR <- targHR * chokeScale

  plot( x = c(0,1), y = c(0,refHR),
        type = "n", xlab = expression(B/B[MSY]),
        ylab = "Harvest Rate", las = 1 )
    lines( x = BBmsy, y = targHR, lwd = 3 )
    lines( x = BBmsy, y = chokeHR, lwd = 2, col = "red" )


} # END plotChokeHCRs()

# plotBioDist()
# Estimates biomass distributions for given years
# and plots them on a stock/species multipanel plot.
plotBioDist <- function(  groupFolders = c("perfInfo_noCorr","simAssErrs_noCorr"),
                          mpFilter = "",
                          distYrs = 2048,
                          qProbs = c(0.05,0.5,0.95),
                          baseBlob = "sim_baseRun",
                          scenOrder = c("noCorr"),
                          hrOrder = c("MSY","noCorr.MSY","MEY","noCorr.MEY"),
                          pch_m = c(21,21,22,22),
                          pt.lwd_m = c(2,0,2,0),
                          col_m = c("red","black","red","black"),
                          bg_m = c(NA,"black",NA,"black") )
{
  nGroups <- length(groupFolders)
  groupList <- list()
  groupID <- c()

  for( groupIdx in 1:nGroups)
  {
    groupFolder <- groupFolders[groupIdx]
    # First, read info files from the relevant
    # sims
    simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                                recursive = FALSE, full.names = FALSE)

    splitCharColumn <- function( str, pattern = "_", 
                                 entry = 1)
    {
      splitz <- unlist(stringr::str_split(str, pattern = pattern))

      splitz[entry]
    }

    info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
                  filter( grepl( mpFilter, mp ) ) %>%
                  arrange(scenario,mp) %>%
                  mutate( HR = sapply(X = mp, FUN = splitCharColumn, entry = 2) ) %>%
                  arrange(match(HR,hrOrder))

    scenLabels  <- unique(info.df$scenario)
    mpLabels    <- info.df$mp

    # Read in the modelStates we extracted
    nModels <- nrow(info.df)
    modelList <- vector(mode = "list", length = nModels)
    names(modelList) <- info.df$simLabel
    for( k in 1:nModels )
    {
      simLabel <- info.df$simLabel[k]
      modelStatePath <- here::here("Outputs",groupFolder,simLabel,paste(simLabel,"_modelStates.RData",sep=""))
      load(modelStatePath)

      modelList[[k]] <- outList
      modelList[[k]]$mp <- info.df$mp[k]

      modelList[[k]]$scenario <- info.df$scenario[k]
      modelList[[k]]$hr <- info.df$HR[k]

    }
    groupList <- c(groupList, modelList)
    groupID <- c(groupID,rep(groupIdx,nModels))
  }
  outList <- NULL

  .loadSim(baseBlob)
  baseBlob <- blob
  blob <- NULL
  gc()

  nSims <- length(groupList)


  # Get reference points
  rp  <- baseBlob$rp[[1]]
  nS  <- modelList[[1]]$nS
  nP  <- modelList[[1]]$nP
  nT  <- modelList[[1]]$nT
  nF  <- modelList[[1]]$nF
  tMP <- modelList[[1]]$tMP
  pT  <- nT - tMP + 1

  speciesNames <- modelList[[1]]$speciesNames
  stockNames <- modelList[[1]]$stockNames

  nReps <- dim(modelList[[1]]$goodReps_isp)[1]
  fYear <- modelList[[1]]$fYear

  distYrIdx <- distYrs - fYear + 1
  nYrs <- length(distYrIdx)

  # Need to put biomasses for the years of interest into an array
  SB_isptM <- array(NA, dim = c(nReps,nS,nP,length(distYrs),nSims) )

  for( m in 1:nSims )
    SB_isptM[,,,1:nYrs,m] <- groupList[[m]]$modelStates$SB_ispt[,,,distYrIdx]

  # Scale by Bmsy
  for( s in 1:nS )
    for( p in 1:nP )
      SB_isptM[,s,p,,] <- SB_isptM[,s,p,,] / rp$FmsyRefPts$BeqFmsy_sp[s,p]
  
  # Quantiles
  SB_qspM <- apply(X = SB_isptM, FUN = quantile, 
                    MARGIN = c(2,3,5), probs = qProbs )

  # 

  par(mfrow = c(nP,nS), oma = c(3,4,2,2), mar = c(.1,.1,.1,.1) )

  for( p in 1:nP )
    for( s in 1:nS )
    {
      plot( x = c(0.5,nSims/2 + 0.5), y = c(0,3), type = "n",
            axes = FALSE )
        mfg <- par("mfg")
        if( mfg[2] == 1)
          axis( side = 2, las = 1)

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 0.5)

        if( mfg[2] == mfg[4] )
          rmtext( txt = stockNames[p], line = 0.05, font = 2, cex =1.5,
                  outer = TRUE)


        abline(h = 1, lwd = 0.8, lty = 3)
        abline(h = .8, lwd = 0.8, lty = 3, col = "orange")
        for( m in 1:nSims )
        {
          if( m %% 2 == 1 )
            xLoc <- ceiling(m/2) - 0.25
          if( m %% 2 == 0 )
            xLoc <- ceiling(m/2) + 0.25

          segments( x0 = xLoc, 
                    y0 = SB_qspM[1,s,p,m],
                    y1 = SB_qspM[3,s,p,m],
                    lty = groupID[m],
                    lwd = 2,
                    col = "grey60" )
          
          hrLabel <- groupList[[m]]$hr
          if( grepl(pattern = "MSY",x = hrLabel ))
            hrpch <- 21
          if( grepl(pattern = "MEY",x = hrLabel ))
            hrpch <- 22



          if( grepl(pattern = "noCorr",x = hrLabel ))
          {
            hrlwd <- 0
            hrbg <- "black"
            hrcol <- "black"
          }
          if( !grepl(pattern = "noCorr",x = hrLabel ))
          {
            hrlwd <- 2
            hrbg <- NA
            hrcol <- "red"
          }

          points( x = xLoc, y = SB_qspM[2,s,p,m],
                  pch = hrpch, lwd = hrlwd, col = hrcol,
                  bg = hrbg, cex = 2  )

          box()

        }

    }
    mtext( side = 2, outer = TRUE, text = expression(B/B[MSY]),
            font = 2, line = 2 )
}

# rmtext()
# Refactored procedure to plot right hand inner
# margin mtext with the bottom towards the middle
# of the plot
rmtext <- function( line = .05, 
                    txt = "Sample", 
                    font = 1,
                    cex = 1,
                    outer = FALSE,
                    yadj = .5)
{
  corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)

  xRange <- corners[2] - corners[1]


  if( outer )
    par(xpd = NA) #Draw outside the figure region
  if( !outer )
    par(xpd = TRUE)
  text( x = corners[2] + line*xRange, 
        y = yadj * sum(corners[3:4]), 
        labels = txt, srt = 270,
        font = font, cex = cex )
  par(xpd = FALSE)
} # END rmtext()


# Refactored procedure to make rotated
# x axis labels on a plot
xaxisRot <- function( at = 1:10, 
                      labs = 1:10,
                      font = 1,
                      cex = 1,
                      line = .05,
                      rot = 60)
{
  corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)

  yRange <- corners[4] - corners[3]

  par(xpd = NA)
  text( x = at, 
        y = corners[3] - line*yRange, 
        labels = labs, srt = rot,
        adj = 1,
        font = font, cex = cex )
  par(xpd = FALSE)
} # END rmtext()


# plotDynRefPoints_sp()
# Plots distributions of stochastically optimised Bmsy
# for a set of simulations
plotDynBmsy_sp <- function( groupFolder = "detRuns_varEff",
                            mpFilter = "omni",
                            scenOrder = c("infElast_det","constElast_det") )
{
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
                filter( grepl( mpFilter, mp ) ) %>%
                arrange(scenario,mp)

  scenLabels <- unique(info.df$scenario)


  # Read in the modelStates we extracted
  nModels <- nrow(info.df)
  modelList <- vector(mode = "list", length = nModels)
  names(modelList) <- info.df$simLabel
  for( k in 1:nModels )
  {
    simLabel <- info.df$simLabel[k]
    modelStatePath <- here::here("Outputs",groupFolder,simLabel,paste(simLabel,"_modelStates.RData",sep=""))
    load(modelStatePath)

    modelList[[k]] <- outList

  }
  outList <- NULL

  # Get reference points
  rp <- modelList[[1]]$rp
  nS <- modelList[[1]]$nS
  nP <- modelList[[1]]$nP
  nT <- modelList[[1]]$nT
  nF <- modelList[[1]]$nF

  speciesNames <- modelList[[1]]$speciesNames
  stockNames <- modelList[[1]]$stockNames

  Bmsy_sp <- rp$FmsyRefPts$BeqFmsy_sp

  Emey_p  <- rp$EmeyRefPts$Emey_p
  Beq_spe <- rp$refCurves$EffCurves$Beq_spe
  E       <- rp$refCurves$EffCurves$E

  mpJitter <- c(  totCat    = .25,
                  totProfit = -.25 )

  par( mfrow = c( nP, nS ),
        oma = c(4,7,2,2),
        mar = c(.1,.1,.1,.1) )

  for( p in 1:nP )
    for( s in 1:nS )
    {
      BmsyMS <- rp$EmsyMSRefPts$BeqEmsy_sp[s,p]
      relBmsyMS <- BmsyMS/Bmsy_sp[s,p]

      # Calculate BmeyMS
      BmeyMS    <- getSplineVal(x = E, y = Beq_spe[s,p,], p = Emey_p[p] )
      relBmeyMS <- BmeyMS/Bmsy_sp[s,p]

      plot( x = c(0,3.5),
            y = c(0.5,4.5), type = "n", axes = FALSE,
            yaxs = "i")
    

      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis( side = 1 )
      if( mfg[2] == 1 )
        axis( side = 2, at = 1:length(scenOrder), labels = scenOrder, las = 1 )
      box(lwd = 2)
      abline( h = 0.5 + 1:3, lwd = 2 )
      abline( v = c(1), lty = 2, lwd = 1)
      abline( v = c(relBmsyMS), lty = 2, lwd = 2, col = "darkgreen")
      abline( v = c(relBmeyMS), lty = 2, lwd = 2, col = "salmon")

      for( k in 1:nModels )
      {
        Bdist <- modelList[[k]]$stateDists$SB_ispt[,s,p]
        Bmean <- modelList[[k]]$stateMeans$SB_ispt[s,p]
        relBdist <- Bdist / Bmsy_sp[s,p]
        relBmean <- Bmean / Bmsy_sp[s,p]

        scenName  <- info.df$scenario[k]
        mpName    <- info.df$mp[k]

        yIdx <- which(scenOrder == scenName)
        if( grepl(pattern = "totCat", x = mpName ) )
        {
          yJitter <- mpJitter["totCat"]
          MPpch   <- 21
        }
        if( grepl(pattern = "totProfit", x = mpName ) )
        {
          yJitter <- mpJitter["totProfit"]
          MPpch   <- 24
        }

        # Plot central 95%
        segments( x0 = relBdist[1], x1 = relBdist[3],
                  y0 = yIdx + yJitter, lwd = 3, col = "grey50" )
        # Plot median
        points( x = relBdist[2], y = yIdx + yJitter,
                pch = MPpch, cex = 2, bg = "black" )
        # Plot mean
        points( x = relBmean, y = yIdx + yJitter,
                pch = MPpch, cex = 2, bg = NA, col = "red" )



      }

      if( s == 1 & p == 1 )
        legend( x = "topleft", bg = "white",
                legend = c( "maxCatch",
                            "maxProfits",
                            "central 95%",
                            "MS Bmsy"),
                pch = c(21,24,NA,NA),
                pt.bg = c("black","black",NA,NA),
                col = c("black","black","grey50","darkgreen"),
                lwd = c(NA,NA,3,2),
                lty = c(NA,NA,1,2))

      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = .5)

      if( mfg[2] == mfg[4] )
        rmtext( txt = stockNames[p], font = 2, line = .05, outer = TRUE,
                cex = 1.5)
    }

    mtext( side = 1, outer = TRUE, text = expression(B/B[MSY]),
            line = 2.5 )

}

# plotDynUmsy_sp()
# Plots distribution of dynamically optimised harvest rates
# for a set of simulations
plotDynEconYield_p <- function( groupFolder = "omni_econYield_constE_Nov6",
                                mpFilter = "freeEff",
                                scenOrder = c("noCorr","corrRecDevs","corrPriceDevs","corrRecPrice") )
{
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
                filter( grepl( mpFilter, mp ) ) %>%
                arrange(scenario,mp)

  scenLabels <- unique(info.df$scenario)


  # Read in the modelStates we extracted
  nModels <- nrow(info.df)
  modelList <- vector(mode = "list", length = nModels)
  names(modelList) <- info.df$simLabel
  for( k in 1:nModels )
  {
    simLabel <- info.df$simLabel[k]
    modelStatePath <- here::here("Outputs",groupFolder,simLabel,paste(simLabel,"_modelStates.RData",sep=""))
    load(modelStatePath)

    modelList[[k]] <- outList

  }
  outList <- NULL

  # Get reference points
  rp <- modelList[[1]]$rp
  nS <- modelList[[1]]$nS
  nP <- modelList[[1]]$nP
  nT <- modelList[[1]]$nT
  nF <- modelList[[1]]$nF

  speciesNames <- modelList[[1]]$speciesNames
  stockNames <- modelList[[1]]$stockNames

  # Pull static eqbria here
  MEY_p           <- rp$EmeyRefPts$MEY_p
  Emsy_p          <- rp$EmsyMSRefPts$EmsyMS_p
  Emey_p          <- rp$EmeyRefPts$Emey_p
  econYeq_pe      <- rp$EmeyRefPts$econYeq_pe
  E               <- rp$refCurves$EffCurves$E

  # Scale by MEY
  relMEY_p        <- rep(1,nP)

  dynEconYield_kp  <- array( NA, dim = c(nModels,nP) )

  econYeqEmsy_p <- rep(0,nP)
  releconYeqEmsy_p <- rep(0,nP)


  mpJitter <- c(  totCat    = .25,
                  totProfit = -.25 )

  par( mfrow = c( nP, 1 ),
        oma = c(4,7,2,2),
        mar = c(.1,.1,.1,.1) )

  for( p in 1:nP )
  {
    econYeqEmsy_p[p]    <- getSplineVal(x = E, y = econYeq_pe[p,], p = Emsy_p[p])
    releconYeqEmsy_p[p] <- econYeqEmsy_p[p]/MEY_p[p]

    plot( x = c(-1,2),
          y = c(0.5,4.5), type = "n", axes = FALSE,
          yaxs = "i")

    mfg <- par("mfg")
    if(mfg[1] == mfg[3])
      axis( side = 1 )
    if( mfg[2] == 1 )
      axis( side = 2, at = 1:4, labels = scenOrder, las = 1 )
    box(lwd = 2)
    abline( h = 0.5 + 1:3, lwd = 2 )
    abline( v = c(1), lty = 2, lwd = 1)
    abline( v = c(1), lty = 2, lwd = 2, col = "salmon")
    abline( v = c(releconYeqEmsy_p[p]), lty = 2, lwd = 2, col = "darkgreen")

    for( k in 1:nModels )
    {

      econYeqdist <- modelList[[k]]$stateDists$profit_ipft[,p,2]
      econYeqmean <- modelList[[k]]$stateMeans$profit_ipft[p,2]

      relYeqdist <- econYeqdist / MEY_p[p]
      relYeqmean <- econYeqmean / MEY_p[p]

      dynEconYield_kp[k,p] <- relYeqdist[2]

      scenName  <- info.df$scenario[k]
      mpName    <- info.df$mp[k]

      yIdx <- which(scenOrder == scenName)
      if( grepl(pattern = "totCat", x = mpName ) )
      {
        yJitter <- mpJitter["totCat"]
        MPpch   <- 21
      }
      if( grepl(pattern = "totProfit", x = mpName ) )
      {
        yJitter <- mpJitter["totProfit"]
        MPpch   <- 24
      }

      # Plot central 95%
      segments( x0 = relYeqdist[1], x1 = relYeqdist[3],
                y0 = yIdx + yJitter, lwd = 3, col = "grey50" )
      # Plot median
      points( x = relYeqdist[2], y = yIdx + yJitter,
              pch = MPpch, cex = 2, bg = "black" )
      # Plot mean
      points( x = relYeqmean, y = yIdx + yJitter,
              pch = MPpch, cex = 2, bg = NA, col = "red" )



    }

    if( p == 1 )
      legend( x = "topright", bg = "white",
              legend = c( "maxCatch",
                          "maxProfits",
                          "central 95%",
                          "P(Emsy)",
                          "MEY"),
              pch = c(21,24,NA,NA,NA),
              pt.bg = c("black","black",NA,NA,NA),
              col = c("black","black","grey50","darkgreen","salmon"),
              lwd = c(NA,NA,3,2,2),
              lty = c(NA,NA,1,2,2))

    if( mfg[2] == mfg[4] )
      rmtext( txt = stockNames[p], font = 2, line = .02, outer = TRUE,
              cex = 1.5)
  }

  mtext( side = 1, outer = TRUE, text = "Producer surplus relative to MEY",
          line = 2.5 )

  # write.csv(Umey_sp, file = "Umey_sp.csv")
  # write.csv(dynUmey_sp, file = "dynUmey_sp.csv")

} # END plotDynUmsy_sp

# plotDynUmsy_sp()
# Plots distribution of dynamically optimised harvest rates
# for a set of simulations
plotDynOptEffort_p <- function( groupFolder = "omni_econYield_splineE_long",
                                mpFilter = "freeEff",
                                scenOrder = c("noCorr","corrRecDevs","corrPriceDevs","corrRecPrice") )
{
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
                filter( grepl( mpFilter, mp ) ) %>%
                arrange(scenario,mp)

  scenLabels <- unique(info.df$scenario)


  # Read in the modelStates we extracted
  nModels <- nrow(info.df)
  modelList <- vector(mode = "list", length = nModels)
  names(modelList) <- info.df$simLabel
  for( k in 1:nModels )
  {
    simLabel <- info.df$simLabel[k]
    modelStatePath <- here::here("Outputs",groupFolder,simLabel,paste(simLabel,"_modelStates.RData",sep=""))
    load(modelStatePath)

    modelList[[k]] <- outList

  }
  outList <- NULL

  # Get reference points
  rp <- modelList[[1]]$rp
  nS <- modelList[[1]]$nS
  nP <- modelList[[1]]$nP
  nT <- modelList[[1]]$nT
  nF <- modelList[[1]]$nF

  speciesNames <- modelList[[1]]$speciesNames
  stockNames <- modelList[[1]]$stockNames

  Emsy_p  <- rp$EmsyMSRefPts$EmsyMS_p
  Emey_p  <- rp$EmeyRefPts$Emey_p

  relEmey_p   <- Emey_p/Emsy_p

  dynOptE_p <- rep(0,nP)


  mpJitter <- c(  totCat    = .25,
                  totProfit = -.25 )

  par( mfrow = c( nP, 1 ),
        oma = c(4,7,2,2),
        mar = c(.1,.1,.1,.1) )

  for( p in 1:nP )
  {
    plot( x = c(0,2),
          y = c(0.5,4.5), type = "n", axes = FALSE,
          yaxs = "i")

    mfg <- par("mfg")
    if(mfg[1] == mfg[3])
      axis( side = 1 )
    if( mfg[2] == 1 )
      axis( side = 2, at = 1:4, labels = scenOrder, las = 1 )
    box(lwd = 2)
    abline( h = 0.5 + 1:3, lwd = 2 )
    abline( v = c(1), lty = 2, lwd = 1)
    abline( v = c(1), lty = 2, lwd = 2, col = "darkgreen")
    abline( v = c(relEmey_p[p]), lty = 2, lwd = 2, col = "salmon")

    for( k in 1:nModels )
    {
      Edist <- modelList[[k]]$stateDists$E_ipft[,p,2]
      Emean <- modelList[[k]]$stateMeans$E_ipft[p,2]

      relEdist <- Edist / Emsy_p[p]
      relEmean <- Emean / Emsy_p[p]

      dynOptE_p[p] <- Edist[2]

      scenName  <- info.df$scenario[k]
      mpName    <- info.df$mp[k]

      yIdx <- which(scenOrder == scenName)
      if( grepl(pattern = "totCat", x = mpName ) )
      {
        yJitter <- mpJitter["totCat"]
        MPpch   <- 21
      }
      if( grepl(pattern = "totProfit", x = mpName ) )
      {
        yJitter <- mpJitter["totProfit"]
        MPpch   <- 24
      }

      # Plot central 95%
      segments( x0 = relEdist[1], x1 = relEdist[3],
                y0 = yIdx + yJitter, lwd = 3, col = "grey50" )
      # Plot median
      points( x = relEdist[2], y = yIdx + yJitter,
              pch = MPpch, cex = 2, bg = "black" )
      # Plot mean
      points( x = relEmean, y = yIdx + yJitter,
              pch = MPpch, cex = 2, bg = NA, col = "red" )



    }

    if( p == 1 )
      legend( x = "topright", bg = "white",
              legend = c( "maxCatch",
                          "maxProfits",
                          "central 95%",
                          "MS Umsy",
                          "Umey"),
              pch = c(21,24,NA,NA,NA),
              pt.bg = c("black","black",NA,NA,NA),
              col = c("black","black","grey50","darkgreen","salmon"),
              lwd = c(NA,NA,3,2,2),
              lty = c(NA,NA,1,2,2))

    if( mfg[2] == mfg[4] )
      rmtext( txt = stockNames[p], font = 2, line = .02, outer = TRUE,
              cex = 1.5)
  }

  mtext( side = 1, outer = TRUE, text = expression(E/E[MSY]),
          line = 2.5 )

  # write.csv(Umey_sp, file = "Umey_sp.csv")
  # write.csv(dynUmey_sp, file = "dynUmey_sp.csv")

} # END plotDynUmsy_sp



# plotDynUmsy_sp()
# Plots distribution of dynamically optimised harvest rates
# for a set of simulations
plotDynUmsy_sp <- function( groupFolder = "omniRuns_econYield_noCorr",
                            mpFilter = "freeEff",
                            baseRP = "sim_baseRun",
                            scenOrder = c("noCorr") )
{
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
                filter( grepl( mpFilter, mp ) ) %>%
                arrange(scenario,mp)

  scenLabels <- unique(info.df$scenario)


  # Read in the modelStates we extracted
  nModels <- nrow(info.df)
  modelList <- vector(mode = "list", length = nModels)
  names(modelList) <- info.df$simLabel
  for( k in 1:nModels )
  {
    simLabel <- info.df$simLabel[k]
    modelStatePath <- here::here("Outputs",groupFolder,simLabel,paste(simLabel,"_modelStates.RData",sep=""))
    load(modelStatePath)

    modelList[[k]] <- outList

  }
  outList <- NULL

  # Get reference points
  rp <- modelList[[1]]$rp
  nS <- modelList[[1]]$nS
  nP <- modelList[[1]]$nP
  nT <- modelList[[1]]$nT
  nF <- modelList[[1]]$nF

  if(!is.null(baseRP))
  {
    .loadSim(baseRP)
    rp <- blob$rp[[1]]
    blob <- NULL
    gc()
  }

  speciesNames <- modelList[[1]]$speciesNames
  stockNames <- modelList[[1]]$stockNames

  Bmsy_sp <- rp$FmsyRefPts$BeqFmsy_sp
  MSY_sp  <- rp$FmsyRefPts$YeqFmsy_sp

  if( is.null(rp$EmeyRefPts$cwEconYieldCurves) )
  {
    Emey_p  <- rp$EmeyRefPts$Emey_p
  } else {
    Emey_p  <- rp$EmeyRefPts$cwEconYieldCurves$cwEmey_p
  }

  Beq_spe <- rp$refCurves$EffCurves$Beq_spe
  Yeq_spe <- rp$refCurves$EffCurves$Yeq_spe
  E       <- rp$refCurves$EffCurves$E

  Umey_sp <- array(0, dim = c(nS,nP))
  dynUmey_sp <- Umey_sp


  Umsy_sp <- MSY_sp / Bmsy_sp

  mpJitter <- c(  totCat    = .25,
                  totProfit = -.25 )



  par( mfrow = c( nP, 1 ),
        oma = c(4,7,2,3),
        mar = c(.1,.1,.1,.1) )


  xJitter <- c(-.25, .25 )

  for( p in 1:nP )
  {
    plot( x = c(0,2.5),
          y = c(0.5,3.5), type = "n", axes = FALSE,
          yaxs = "i" )
      abline( h = 0.5 + 1:3, lwd = 2 )
      mfg <- par("mfg")
      
      axis( side = 2, at = c(1:nS), labels = speciesNames[nS:1],
            las = 1, tick = FALSE, font = 2 )

      if(mfg[1] == mfg[3])
        axis( side = 1 )
      
      # if(mfg[1] == 1)
      #   mtext(  side = 3, text = speciesNames[s],
      #           line = 0, font = 2)
      #   axis( side = 3, at = 1:3, labels = speciesNames, tick = FALSE,
      #         font = 2, line = -1, lab.cex = 2 )
      
      # abline( v = c(1.5,2.5), lwd = 2 )
      abline( v = seq(from = 0, to = 2.5, by = .5), lwd = .8, col = "grey80", lty = 3 )
      abline( v = 1, lty = 2, lwd = 2, col = "grey60")
      box( lwd = 2 )

    for( s in 1:nS )
    {
      BmsyMS  <- rp$EmsyMSRefPts$BeqEmsy_sp[s,p]
      MSYMS   <- rp$EmsyMSRefPts$YeqEmsy_sp[s,p]

      UmsyMS    <- MSYMS/BmsyMS
      relUmsyMS <- UmsyMS/Umsy_sp[s,p]

      # Calculate BmeyMS
      BmeyMS    <- getSplineVal(x = E, y = Beq_spe[s,p,], p = Emey_p[p] )
      relBmeyMS <- BmeyMS/Bmsy_sp[s,p]

      # Calculate YmeyMS
      YmeyMS    <- getSplineVal(x = E, y = Yeq_spe[s,p,], p = Emey_p[p] )
      relYmeyMS <- YmeyMS/MSY_sp[s,p]

      Umey_sp[s,p]  <- YmeyMS/BmeyMS
      relUmeyMS     <- Umey_sp[s,p]/Umsy_sp[s,p]

      # BmsyMS <- rp$EmeyRefPts$BeqEmsyMS_sp[s,p]
      # relBmsyMS <- BmsyMS/Bmsy_sp[s,p]

      for( k in 1:nModels )
      {
        Udist <- modelList[[k]]$stateDists$U_ispt[,s,p]
        Umean <- modelList[[k]]$stateMeans$U_ispt[s,p]

        relUdist <- Udist / Umsy_sp[s,p]
        relUmean <- Umean / Umsy_sp[s,p]

        dynUmey_sp[s,p] <- Udist[2]

        scenName  <- info.df$scenario[k]
        mpName    <- info.df$mp[k]

        if( grepl(pattern = "totCat", x = mpName ) )
        {
          MPpch   <- 21
          ssU     <- relUmsyMS
        }

        if( grepl(pattern = "totProfit", x = mpName ) )
        {
          MPpch   <- 24
          ssU     <- relUmeyMS
        }

        # Plot central 95%
        segments( x0 = relUdist[1], x1 = relUdist[3],
                  y0 = nS -s + 1 + xJitter[k], lwd = 3, col = "grey50" )
        # Plot median
        points( x = relUdist[2], y = nS -s + 1 + xJitter[k],
                pch = MPpch, cex = 2, bg = "black",
                lwd = 0 )
        points( x = ssU, y = nS -s + 1 + xJitter[k],
                pch = MPpch, cex = 2, bg = NA, col = "red",
                lwd = 2)

        # Plot mean
        # points( x = relUmean, y = yIdx + yJitter,
        #         pch = MPpch, cex = 2, bg = NA, col = "red" )



      }

      # if( s == 1 & p == 1 )
      #   legend( x = "topright", bg = "white",
      #           legend = c( "Catch",
      #                       "Pr. Sp.",
      #                       "Central 95%",
      #                       "UmsyMS",
      #                       "Umey"),
      #           pch = c(21,24,NA,NA,NA),
      #           pt.bg = c("black","black",NA,NA,NA),
      #           col = c("black","black","grey50","darkgreen","salmon"),
      #           lwd = c(NA,NA,3,2,2),
      #           lty = c(NA,NA,1,2,2))

      # if( mfg[1] == 1 )
      #   mtext( side = 3, text = speciesNames[s], font = 2, line = .5)

      if( mfg[2] == mfg[4] )
        rmtext( txt = stockNames[p], font = 2, line = .05, outer = TRUE,
                cex = 1.5)
    }

    if( p == 1 )
      legend( x = "topright", bg = "white", bty = "n",
              legend = c( "Umsy*",
                          "Umey*",
                          "Umsy",
                          "Umey"),
              pch = c(21,24,21,24),
              pt.bg = c("black","black",NA,NA),
              col = c("black","black","red","red"),
              pt.lwd = c(0,0,2,2))   
  }

    mtext( side = 1, outer = TRUE, text = expression(U/U[MSY]),
            line = 2.5 )


    write.csv(Umey_sp, file = "Umey_sp.csv")
    write.csv(dynUmey_sp, file = "dynUmey_sp.csv")

} # END plotDynUmsy_sp

# plotDynUmsy_sp()
# Plots distribution of dynamically optimised harvest rates
# for a set of simulations
plotDynMSY_sp <- function( groupFolder = "omni_econYield_constE_Nov6",
                            mpFilter = "freeEff",
                            scenOrder = c("noCorr","corrRecDevs","corrPriceDevs","corrRecPrice") )
{
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
                filter( grepl( mpFilter, mp ) ) %>%
                arrange(scenario,mp)

  scenLabels <- unique(info.df$scenario)


  # Read in the modelStates we extracted
  nModels <- nrow(info.df)
  modelList <- vector(mode = "list", length = nModels)
  names(modelList) <- info.df$simLabel
  for( k in 1:nModels )
  {
    simLabel <- info.df$simLabel[k]
    modelStatePath <- here::here("Outputs",groupFolder,simLabel,paste(simLabel,"_modelStates.RData",sep=""))
    load(modelStatePath)

    modelList[[k]] <- outList

  }
  outList <- NULL

  # Get reference points
  rp <- modelList[[1]]$rp
  nS <- modelList[[1]]$nS
  nP <- modelList[[1]]$nP
  nT <- modelList[[1]]$nT
  nF <- modelList[[1]]$nF

  speciesNames <- modelList[[1]]$speciesNames
  stockNames <- modelList[[1]]$stockNames

  Bmsy_sp <- rp$FmsyRefPts$BeqFmsy_sp
  MSY_sp  <- rp$FmsyRefPts$YeqFmsy_sp

  Emey_p  <- rp$EmeyRefPts$Emey_p
  Beq_spe <- rp$refCurves$EffCurves$Beq_spe
  Yeq_spe <- rp$refCurves$EffCurves$Yeq_spe
  E       <- rp$refCurves$EffCurves$E


  mpJitter <- c(  totCat    = .25,
                  totProfit = -.25 )

  par( mfrow = c( nP, nS ),
        oma = c(4,7,2,2),
        mar = c(.1,.1,.1,.1) )

  for( p in 1:nP )
    for( s in 1:nS )
    {
      plot( x = c(0,2),
            y = c(0.5,4.5), type = "n", axes = FALSE,
            yaxs = "i")

      MSYMS   <- rp$EmsyMSRefPts$YeqEmsy_sp[s,p]

      relMSYMS <- MSYMS/MSY_sp[s,p]

      # Calculate YmeyMS
      YmeyMS    <- getSplineVal(x = E, y = Yeq_spe[s,p,], p = Emey_p[p] )
      relYmeyMS <- YmeyMS/MSY_sp[s,p]

      # BmsyMS <- rp$EmeyRefPts$BeqEmsyMS_sp[s,p]
      # relBmsyMS <- BmsyMS/Bmsy_sp[s,p]

      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis( side = 1 )
      if( mfg[2] == 1 )
        axis( side = 2, at = 1:4, labels = scenOrder, las = 1 )
      box(lwd = 2)
      abline( h = 0.5 + 1:3, lwd = 2 )
      abline( v = c(1), lty = 2, lwd = 1)
      abline( v = c(relMSYMS), lty = 2, lwd = 2, col = "darkgreen")
      abline( v = c(relYmeyMS), lty = 2, lwd = 2, col = "salmon")

      for( k in 1:nModels )
      {
        Cdist <- modelList[[k]]$stateDists$C_ispt[,s,p]
        Cmean <- modelList[[k]]$stateMeans$C_ispt[s,p]

        relCdist <- Cdist / MSY_sp[s,p]
        relCmean <- Cmean / MSY_sp[s,p]

        scenName  <- info.df$scenario[k]
        mpName    <- info.df$mp[k]

        yIdx <- which(scenOrder == scenName)
        if( grepl(pattern = "totCat", x = mpName ) )
        {
          yJitter <- mpJitter["totCat"]
          MPpch   <- 21
        }
        if( grepl(pattern = "totProfit", x = mpName ) )
        {
          yJitter <- mpJitter["totProfit"]
          MPpch   <- 24
        }

        # Plot central 95%
        segments( x0 = relCdist[1], x1 = relCdist[3],
                  y0 = yIdx + yJitter, lwd = 3, col = "grey50" )
        # Plot median
        points( x = relCdist[2], y = yIdx + yJitter,
                pch = MPpch, cex = 2, bg = "black" )
        # Plot mean
        points( x = relCmean, y = yIdx + yJitter,
                pch = MPpch, cex = 2, bg = NA, col = "red" )



      }

      if( s == 1 & p == 1 )
        legend( x = "topright", bg = "white",
                legend = c( "maxCatch",
                            "maxProfits",
                            "central 95%",
                            "MS MSY"),
                pch = c(21,24,NA,NA),
                pt.bg = c("black","black",NA,NA),
                col = c("black","black","grey50","darkgreen"),
                lwd = c(NA,NA,3,2),
                lty = c(NA,NA,1,2))

      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = .5)

      if( mfg[2] == mfg[4] )
        rmtext( txt = stockNames[p], font = 2, line = .05, outer = TRUE,
                cex = 1.5)
    }

    mtext( side = 1, outer = TRUE, text = expression(C/MSY),
            line = 2.5 )

}

# plotBatchCatchBioTradeoff()
# Plots of catch/biomass tradeoffs for a given
# batch
plotBatchCatchBioTradeoff <- function(  groupFolder = "DERTACS_reruns_Oct10",
                                        prefix = "parBat",
                                        period = 73:82,
                                        lossList = NULL,
                                        qProbs = c(.025,.5,.975),
                                        refPts = "SSrefPts",
                                        dim1 = 1:3, # Species
                                        dim2 = 1:3, # Stock
                                        AMlabs = c( SS = "singleStock",
                                                    HMS = "hierMultiStock",
                                                    SpePool = "speciesPooling",
                                                    SpaPool = "spatialPooling",
                                                    TA = "totalAgg" ),
                                        scenLabs = c( Rich  = "DERfit_HcMcAsSsIdx",
                                                      Poor = "DERfit_AsSsIdx" ),
                                        minSampSize = 100  )
{
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
                filter( grepl( prefix, simLabel ) )

  # Break up MP names into factor levels
  # MP labels are AM_Fsrce_eqbm
  splitMP <- function(  mpLab, 
                        breaks = "_",
                        factorNames = c("AM","Fsrce","eqbm") )
  {
    splitMP <- stringr::str_split(mpLab, breaks)[[1]]

    splitMP
  }

  MPfactors <- lapply( X = info.df$mp, FUN = splitMP )
  MPfactors <- do.call(rbind, MPfactors)
  colnames(MPfactors) <- c("AM","Fsrce","eqbm")
  MPfactors <- as.data.frame(MPfactors)
  

  # Read in the performance tables
  mpTable <- cbind( info.df, MPfactors ) %>%
              mutate( perfPath = here::here("Outputs",groupFolder,simLabel,"simPerfStats.csv"),
                      path = here::here("Outputs",groupFolder,simLabel) ) %>%
              filter(scenario %in% scenLabs, AM %in% AMlabs )

  blobFiles <- file.path(mpTable$path,paste(mpTable$simLabel,".RData",sep = ""))
  lossFiles <- file.path(mpTable$path,"loss.RData")

  lossList <- lapply(X = mpTable$simLabel, FUN = .loadLoss, folder=  groupFolder)

  names(lossList) <- mpTable$simLabel
  nSims <- length(lossList)

  # Now get reference points, any of the blobs will do
  .loadSim( mpTable$simLabel[1], groupFolder = groupFolder )

  rp <- blob$rp[[1]]
  BmsySS_sp <- rp$FmsyRefPts$BeqFmsy_sp
  MSYSS_sp  <- rp$FmsyRefPts$YeqFmsy_sp


  BmsyMS_sp <- rp$EmsyMSRefPts$BeqEmsy_sp
  MSYMS_sp  <- rp$EmsyMSRefPts$YeqEmsy_sp
  

  if(refPts == "MSrefPts" )
  {
    Bmsy_sp <- BmsyMS_sp
    MSY_sp  <- MSYMS_sp
  }

  if( refPts == "SSrefPts" )
  {
    Bmsy_sp <- BmsySS_sp
    MSY_sp  <- MSYSS_sp
  }

  scenarios <- unique(mpTable$scenario)
  AMs       <- unique(mpTable$AM)
  species   <- lossList[[1]]$speciesNames
  stock     <- lossList[[1]]$stockNames
  fYear     <- lossList[[1]]$fYear

  nReps_k <- numeric(length = nSims)
  goodRepsList <- vector(mode = "list", length = nSims)

  for( k in 1:nSims )
  {
    simID <- mpTable$simLabel[k]
    simLoss <- lossList[[1]]$simStates$C_ispt
    nReps_k[k] <- dim(simLoss)[1]
    goodRepsList[[k]]$goodReps_isp <- lossList[[k]]$goodReps_isp
    goodRepsList[[k]]$nReps_sp     <- apply(X = lossList[[k]]$goodReps_isp, FUN = sum, MARGIN = c(2,3) )
    names(goodRepsList)[k] <- simID
  }

  # Pull omni biomass and catch from one of the
  # lossLists
  omniC_ispt <- lossList[[1]]$baseStates$C_ispt[1:100,,,period]
  omniB_ispt <- lossList[[1]]$baseStates$SB_ispt[1:100,,,period]

  

  

  nS <- lossList[[1]]$nS
  nP <- lossList[[1]]$nP
  nT <- lossList[[1]]$nT

  nScen <- length(scenarios)
  nAM   <- length(AMs)
  yrs   <- seq( from = fYear, by = 1, length.out = nT)


  
  # Need to create an array to
  # populate with the simStates
  # scaled by MSY/Bmsy
  # Get nReps
  nReps <- max(nReps_k)
  # Array to hold biomass, with Scenario/AM dimensions
  B_SAispt <- array( NA,  dim = c(nScen,nAM,nReps,nS,nP,nT),
                          dimnames = list(  scenario = scenarios,
                                            AM = AMlabs,
                                            rep = 1:nReps,
                                            species = species[1:nS],
                                            stock = stock[1:nP],
                                            year = yrs ) )
  # Copy for catch
  C_SAispt      <- B_SAispt
  nReps_SAsp    <- B_SAispt[,,1,,,1]

  for( k in 1:nSims )
  {
    scenID    <- mpTable$scenario[k]
    amID      <- mpTable$AM[k]

    nReps_SAsp[scenID,amID,,] <- goodRepsList[[k]]$nReps_sp

    for( s in 1:nS )
      for(p in 1:nP )
      {
        goodRepIdx <- which(goodRepsList[[k]]$goodReps_isp[,s,p])

        B_SAispt[scenID,amID,goodRepIdx,s,p,] <- lossList[[k]]$simStates$SB_ispt[goodRepIdx,s,p,] / Bmsy_sp[s,p]
        C_SAispt[scenID,amID,goodRepIdx,s,p,] <- lossList[[k]]$simStates$C_ispt[goodRepIdx,s,p,] / MSY_sp[s,p]

      }

  }

  for( s in 1:nS )
    for( p in 1:nP )
    {
      omniC_ispt[,s,p,] <- omniC_ispt[,s,p,] / MSY_sp[s,p]
      omniB_ispt[,s,p,] <- omniB_ispt[,s,p,] / Bmsy_sp[s,p]
    }

  omniC_qspt <- apply( X = omniC_ispt, FUN = quantile,
                          MARGIN = c(2,3,4), probs = qProbs)
  omniB_qspt <- apply( X = omniB_ispt, FUN = quantile,
                          MARGIN = c(2,3,4), probs = qProbs)

  omniC_qsp <- apply( X = omniC_qspt, FUN = mean, MARGIN = c(1,2,3) )
  omniB_qsp <- apply( X = omniB_qspt, FUN = mean, MARGIN = c(1,2,3) )



  # Now calculate quantiles
  B_qSAsp <- apply( X = B_SAispt[,,,,,period], FUN = quantile,
                    probs = qProbs,
                    MARGIN =c(1,2,4,5), na.rm = T )
  C_qSAsp <- apply( X = C_SAispt[,,,,,period], FUN = quantile,
                    probs = qProbs,
                    MARGIN =c(1,2,4,5), na.rm = T )

  scenLty <- 1:nScen
  AMpch   <- 21:(21 + nAM - 1)

  AMcols  <- RColorBrewer::brewer.pal(nAM, "Dark2")
  names(AMcols) <- AMlabs
  names(AMpch)  <- AMlabs

  nPP <- length(dim2)
  nSS <- length(dim1)

  # Now plot
  par(  mfcol = c(nPP,nSS),
        mar = c(2,2,1,1),
        oma = c(5,3.5,3,3) )

  for( s in dim1 )
    for( p in dim2 )
    {
      plot( x = range(B_qSAsp[2,,,s,p],omniB_qsp[,s,p],omniC_qsp[,s,p]),
            y = range(C_qSAsp[2,,,s,p],omniB_qsp[,s,p],omniC_qsp[,s,p]),
            type = "n", axes = FALSE )
        mfg <- par( "mfg" )
        
          axis( side = 1 )
          axis( side = 2, las = 1 )

        if( mfg[2] == mfg[4] )
          rmtext( txt = stock[p], line = 0.05, font = 2, cex = 1.5,
                  outer = TRUE )
        if( mfg[1] == 1 )
          mtext( side = 3, text = species[s], font = 2)
        box()
        grid()
        # abline( h = 1, lty = 2, lwd = 2)
        abline( v = .4 * BmsySS_sp[s,p]/Bmsy_sp[s,p], 
                lty = 2, lwd = 2, col = "red")

        # Cat range of omniMgr
        segments( x0 = omniB_qsp[2,s,p],
                  x1 = omniB_qsp[2,s,p],
                  y0 = omniC_qsp[1,s,p],
                  y1 = omniC_qsp[3,s,p], lwd = 2,
                  col = "grey30" )

        segments( x0 = omniB_qsp[1,s,p],
                  x1 = omniB_qsp[3,s,p],
                  y0 = omniC_qsp[2,s,p],
                  y1 = omniC_qsp[2,s,p], lwd = 2,
                  col = "grey30" )
        

        points( x = BmsyMS_sp[s,p]/Bmsy_sp[s,p], y = MSYMS_sp[s,p]/MSY_sp[s,p],
                col = "steelblue", pch = 16, cex = 2,
                lwd = 1.5 )      

        points( x = BmsySS_sp[s,p]/Bmsy_sp[s,p], y = MSYSS_sp[s,p]/MSY_sp[s,p],
                col = "darkgreen", pch = 9, cex = 2,
                lwd = 1.5 )


        for( scenIdx in 1:nScen )
        {
          scenID <- scenarios[scenIdx]
          
          Bmat <- B_qSAsp[,scenID,,s,p]
          Cmat <- C_qSAsp[,scenID,,s,p]

          plotOrder <- order(Bmat[2,])

          nConvReps_A <- nReps_SAsp[scenID,,s,p]

          
          bg <- AMcols
          bg[nConvReps_A < minSampSize] <- NA
          
          scenID <- scenarios[scenIdx]
          
          # browser()
          # segments( x0 = Bmat[2,],
          #           x1 = Bmat[2,],
          #           y0 = Cmat[1,],
          #           y1 = Cmat[3,],
          #           col = AMcols[colnames(Bmat)],
          #           lty = 1, lwd = 2 )
          # segments( x0 = Bmat[1,],
          #           x1 = Bmat[3,],
          #           y0 = Cmat[2,],
          #           y1 = Cmat[2,],
          #           col = AMcols[colnames(Bmat)],
          #           lty = 1, lwd = 2 )

          points( x   = Bmat[2,plotOrder],
                  y   = Cmat[2,plotOrder],
                  pch = AMpch[colnames(Bmat)][plotOrder],
                  bg = bg[colnames(Bmat)][plotOrder],
                  col = AMcols[colnames(Bmat)][plotOrder] )

          lines(  x   = Bmat[2,plotOrder],
                  y   = Cmat[2,plotOrder],
                  lty = scenLty[scenIdx],
                  lwd = 1, col = "grey60" )          
        }
    }

  mtext( side = 1, outer = TRUE,
          text = expression(B/B[MSY]), line = 2 )

  mtext( side = 2, outer = TRUE,
          text = expression(C/MSY), line = 2 )

  par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend( x       = "bottom",
          horiz   = TRUE,
          bty     = "n",
          legend  = c(names(AMlabs),names(scenLabs)),
          pch     = c(AMpch[AMlabs],NA,NA,NA,NA),
          pt.bg   = c(AMcols[AMlabs],NA,NA,NA,NA), 
          lty     = c(rep(NA,5),1:4),
          col     = c(AMcols[AMlabs],rep("grey60",4)),
          lwd     = 1 )


}

# plotRetroBio_Scenario()
# A sample plot of assessment performance
# for each AM type within a scenario.
# Requires 5 completed simulation objects.
plotRetroBio_Scenario <- function(  groupFolder = "DERTACS_reruns_sep24",
                                    iRep = 1,
                                    scenName = "DERfit_HcMcAsSsIdx",
                                    prefix = "parBat",
                                    species = "Dover",
                                    stock = "HSHG" )
{
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
                filter( grepl( prefix, simLabel ),
                        scenario == scenName )

  # Break up MP names into factor levels
  # MP labels are AM_Fsrce_eqbm
  splitMP <- function(  mpLab, 
                        breaks = "_",
                        factorNames = c("AM","Fsrce","eqbm") )
  {
    splitMP <- stringr::str_split(mpLab, breaks)[[1]]

    outList <- vector(  mode = "list", 
                        length = length(splitMP) )
    names(outList) <- factorNames
    for( k in 1:length(splitMP))
      outList[[k]] <- splitMP[k]

    outList
  }

  MPfactors <- lapply( X = info.df$mp, FUN = splitMP )
  MPfactors <- do.call(rbind, MPfactors)

  # Read in the performance tables
  mpTable <- cbind( info.df, MPfactors ) %>%
              mutate( perfPath = here::here("Outputs",groupFolder,simLabel,"simPerfStats.csv"),
                      path = here::here("Outputs",groupFolder,simLabel) )

  blobFiles <- file.path(mpTable$path,paste(mpTable$simLabel,".RData",sep = ""))
  lossFiles <- file.path(mpTable$path,"loss.RData")

  lossList <- lapply(X = mpTable$simLabel, FUN = .loadLoss, folder=  groupFolder)

  names(lossList) <- mpTable$AM

  # Now pull dimensions from the blob
  nT  <- lossList[[1]]$nT
  nS  <- lossList[[1]]$nS
  nP  <- lossList[[1]]$nP
  pT  <- dim(lossList[[1]]$retroSB_itspt)[2]
  tMP <- lossList[[1]]$tMP

  speciesNames <- lossList[[1]]$speciesNames[1:3]
  stockNames   <- lossList[[1]]$stockNames[1:3]

  specIdx   <- which(speciesNames == species)
  stockIdx  <- which(stockNames == stock)

  # Years for plotting
  fYear <- lossList[[1]]$fYear
  years <- seq( from = fYear, by = 1, length.out = nT)

  # Now we want to make arrays to hold
  # Rows/cols are: 
  # 1. SS model fits / HSHG area, 
  # 2. MS model fits / HSHG area,
  # 3. Aggregate model fits / Dover SpatPooled, HSHG SpecPooled, TotalAgg

  # 1. OM Biomass
  omB_rct <- array(NA, dim = c(3,3,nT) )

  # 2. Retrospective Bt fits
  retroB_trct <- array(NA, dim = c(pT,3,3,nT) )

  # Now step through blobList and recover the things we want
  omB_rct[1,,]  <- lossList$singleStock$simStates$SB_ispt[iRep,1:nS,stockIdx,]
  omB_rct[2,,]  <- lossList$hierMultiStock$simStates$SB_ispt[iRep,1:nS,stockIdx,]
  omB_rct[3,1,] <- apply( X = lossList$spatialPooling$simStates$SB_ispt[iRep,specIdx,1:nP,], FUN = sum, MARGIN = 2)
  omB_rct[3,2,] <- apply( X = lossList$speciesPooling$simStates$SB_ispt[iRep,1:nS,stockIdx,], FUN = sum, MARGIN = 2)
  omB_rct[3,3,] <- apply( X = lossList$totalAgg$simStates$SB_ispt[iRep,1:nS,1:nP,], FUN = sum, MARGIN = 3)

  retroB_trct[,1,,]  <- lossList$singleStock$retroSB_itspt[iRep,,,stockIdx,]
  retroB_trct[,2,,]  <- lossList$hierMultiStock$retroSB_itspt[iRep,,,stockIdx,]
  retroB_trct[,3,1,] <- lossList$spatialPooling$retroSB_itspt[iRep,,specIdx,1,]
  retroB_trct[,3,2,] <- lossList$speciesPooling$retroSB_itspt[iRep,,1,stockIdx,]
  retroB_trct[,3,3,] <- lossList$totalAgg$retroSB_itspt[iRep,,1,1,]

  retroB_trct[retroB_trct <= 0] <- NA

  blobList <- NULL
  gc()

  
  rowLabs <- c("Single-stock", "Heriarchical\nmulti-stock","Pooled Data")
  rowLines <- c(3,5,3)

  stockSpecRowTitles <- paste(stock,speciesNames)
  aggRowTitles <- c(paste(species,"Sole Pooled"), paste(stock,"Pooled"), "Total Aggregation")


  plotTitles <- rbind( stockSpecRowTitles,stockSpecRowTitles,aggRowTitles)

  par( mfrow = c(3,3), mar = c(.1,2,2,1), oma = c(3,3,3,4) )
  for( rIdx in 1:3 )
    for( cIdx in 1:3 )
    {
      maxBt <- max( omB_rct[rIdx,cIdx,], retroB_trct[,rIdx,cIdx,], na.rm = T)

      plot( x = range(years), y = c(0,maxBt),
            type = "n", axes = FALSE )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )


      if( mfg[2] == mfg[4] )
        rmtext( txt = rowLabs[rIdx], 
                line = .05, 
                font = 2, 
                cex = 1.5,
                outer = TRUE,
                yadj = .5 )

      mtext( side = 3, text = plotTitles[rIdx,cIdx], font = 2)

      axis( side = 2, las = 1 )
      grid()
      box()
      lines( x = years, y = omB_rct[rIdx,cIdx,],
              col = "red", lty = 1, lwd = 3 )
      for( tt in 1:pT )
      {
        lines( x = years, y = retroB_trct[tt,rIdx,cIdx,],
                col = "grey60", lwd = 1 )
      }
      abline(v = years[tMP], lty = 2)
  
    }

  mtext( side = 1, text = "Year", outer = TRUE, line = 2, font = 2)
  mtext( side = 2, text = "Spawning Biomass", outer = TRUE, line = 1.5, font = 2)
} # END plotRetroBio_Scenario()


# plotTulipEcon_sp
# Simulation envelopes of price and revenue at the species/stock
# level. Doesn't include costs of fishing or profit, as these
# are at the area level aggregated over species.
plotTulipEcon_sp <- function( obj = NULL,
                              simNum = 1,
                              nTrace = 3,
                              groupFolder = "omni_econYield_splineE_long_Jan4",
                              price = TRUE,
                              revenue = TRUE )
{
  # Load sim if not passed in
  if( is.null(obj))
  {
    .loadSim( simNum, groupFolder )
    obj <- blob
  }


  goodReps_isp  <- obj$goodReps_isp
  nReps         <- dim(goodReps_isp)[1]

  # Get catch and econ arrays
  C_ispt            <- obj$om$C_ispt
  Rev_ispft         <- obj$om$Rev_ispft
  landVal_ist       <- obj$om$landVal_ist
  basePrice_ist     <- obj$om$basePrice_ist
  effCost_ipft      <- obj$om$effCost_ipft  


  # Get control list
  ctlList <- obj$ctlList

  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  pT      <- obj$ctlList$opMod$pT

  # Labeling info
  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fYear         <- obj$ctlList$opMod$fYear

  years         <- seq( from = fYear, by = 1, length.out = nT )  

  traceIdx <- sample(x = 1:nReps, size = nTrace)

  
  
  basePrice_qst <- apply(X = basePrice_ist, FUN = quantile,
                          probs = c(0.025, 0.5, 0.975),
                          MARGIN = c(2,3), na.rm = T )
  landVal_qst   <- apply(X = landVal_ist, FUN = quantile,
                          probs = c(0.025, 0.5, 0.975),
                          MARGIN = c(2,3), na.rm = T )
  Rev_qspft     <- apply(X = Rev_ispft, FUN = quantile,
                          probs = c(0.025, 0.5, 0.975),
                          MARGIN = c(2,3,4,5), na.rm = T )



  # Make nP x nS array
  # Now plot
  par(  mfcol = c(nP, nS ),
        mar = c(.1,1.5,.1,1.5),
        oma = c(3,3,2,2) )
  for( s in 1:nS )
    for( p in 1:nP )
    {
      maxY <- 0
      if( price )
        maxY <- max(maxY,basePrice_qst[,s,tMP:nT],landVal_qst[,s,tMP:nT],na.rm =T)
      if( revenue )
        maxY <- max(maxY,Rev_qspft[,s,p,2,tMP:nT],na.rm =T)

      plot( x = range(years[tMP:nT]),
            y = c(0,maxY),
            type = "n", xlab = "", ylab = "", axes = FALSE )

        mfg <- par("mfg")
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        
        axis( side = 2, las = 1 )
        if( mfg[2] == mfg[4] )
          rmtext( txt = stockNames[p], line = .1, font = 2, cex = 1.5, outer = TRUE )
        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s] )
        box()
        grid()

        if( price )
        {
          # Price when catch == MSY
          polygon(  x = c(years[(tMP-11):nT],years[nT:(tMP-11)]),
                    y = c(basePrice_qst[1,s,(tMP-11):nT],rev(basePrice_qst[3,s,(tMP-11):nT])),
                    border = NA, 
                    col = scales::alpha("darkgreen",.3) )
          lines( x = years[(tMP-11):nT], y = basePrice_qst[2,s,(tMP-11):nT],
                  lwd = 3, col = "darkgreen" )

          for( i in traceIdx )
            lines( x = years[(tMP-11):nT], y = basePrice_ist[i,s,(tMP-11):nT],
                  lwd = .8, col = "darkgreen" )            


          # Price from elasticity of demand (downward sloping demand curve)
          polygon(  x = c(years[(tMP-11):nT],years[nT:(tMP-11)]),
                    y = c(landVal_qst[1,s,(tMP-11):nT],rev(landVal_qst[3,s,(tMP-11):nT])),
                    border = NA, 
                    col = scales::alpha("salmon",.3) )
          lines( x = years[(tMP-11):nT], y = landVal_qst[2,s,(tMP-11):nT],
                  lwd = 3, col = "salmon" )

          for( i in traceIdx )
            lines( x = years[(tMP-11):nT], y = landVal_ist[i,s,(tMP-11):nT],
                  lwd = .8, col = "salmon" )    


        }

        if( revenue )
        {
          polygon(  x = c(years[(tMP-11):nT],years[nT:(tMP-11)]),
                    y = c(Rev_qspft[1,s,p,2,(tMP-11):nT],rev(Rev_qspft[3,s,p,2,(tMP-11):nT])),
                    border = NA, 
                    col = scales::alpha("grey50",.3) )
          lines( x = years[(tMP-11):nT], y = Rev_qspft[2,s,p,2,(tMP-11):nT],
                  lwd = 3, col = "grey50" )

          for( i in traceIdx )
            lines( x = years[(tMP-11):nT], y = Rev_ispft[i,s,p,2,(tMP-11):nT],
                  lwd = .8, col = "grey50" )            
        }


    }

    mtext( side = 1, text = "Year", outer = TRUE )
    mtext( side = 2, text = "Unit Price ($/kg), Revenue ($m)", outer = TRUE)

} # END plotTulipEcon_sp

# Envelopes of simulated assessment errors
plotTulipAssError <- function(  simNum = 1,
                                obj = NULL,
                                groupFolder = "DLSurveys_.3tau_Long",
                                save = FALSE,
                                proj = TRUE,
                                clearBadReps = FALSE )
{

  if( is.null(obj))
  {
    .loadSim( simNum, groupFolder )
    obj <- blob
  }

  goodReps_isp  <- obj$goodReps_isp
  nReps         <- dim(goodReps_isp)[1]

  # Get biomass arrays
  SB_ispt        <- obj$om$SB_ispt
  VB_ispt        <- obj$om$vB_ispft[,,,2,]
  totB_ispt      <- obj$om$B_ispt
  retroSB_itspt  <- obj$mp$assess$retroSB_itspt


  ctlList <- obj$ctlList

  retroSB_itspt[retroSB_itspt < 0] <- NA
  # nReps     <- sum(goodReps)

  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT


  # Aggregate OM biomasses to match AM biomass
  if( ctlList$mp$data$spatialPooling )
  {
    newSB_ispt        <- apply( X = SB_ispt, FUN = sum, MARGIN = c(1,2,4), na.rm = T )
    newVB_ispt        <- apply( X = VB_ispt, FUN = sum, MARGIN = c(1,2,4), na.rm = T )
    newtotB_ispt      <- apply( X = totB_ispt, FUN = sum, MARGIN = c(1,2,4), na.rm = T )

    SB_ispt    <- SB_ispt[,,1,,drop = FALSE]
    SB_ispt[,,1,] <- newSB_ispt
    VB_ispt    <- VB_ispt[,,1,,drop = FALSE]
    VB_ispt[,,1,] <- newVB_ispt
    totB_ispt  <- totB_ispt[,,1,,drop = FALSE]
    totB_ispt[,,1,] <- newtotB_ispt     

  }

  if( ctlList$mp$data$speciesPooling )
  {
    newSB_ispt        <- apply( X = SB_ispt, FUN = sum, MARGIN = c(1,3,4), na.rm = T )
    newVB_ispt        <- apply( X = VB_ispt, FUN = sum, MARGIN = c(1,3,4), na.rm = T )
    newtotB_ispt      <- apply( X = totB_ispt, FUN = sum, MARGIN = c(1,3,4), na.rm = T )

    SB_ispt    <- SB_ispt[,1,,,drop = FALSE]
    SB_ispt[,1,,] <- newSB_ispt
    VB_ispt    <- VB_ispt[,1,,,drop = FALSE]
    VB_ispt[,1,,] <- newVB_ispt
    totB_ispt  <- totB_ispt[,1,,,drop = FALSE]
    totB_ispt[,1,,] <- newtotB_ispt

  }

  SB_ispt[SB_ispt == 0]     <- NA
  VB_ispt[VB_ispt == 0]     <- NA
  totB_ispt[totB_ispt == 0] <- NA

  nSS     <- dim( SB_ispt)[2]
  nPP     <- dim( SB_ispt)[3]

  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fYear         <- obj$ctlList$opMod$fYear
  pT            <- obj$ctlList$opMod$pT

  if( nPP == 1 )
    stockNames <- "Spatial Pooled"

  if( nSS == 1 )
    speciesNames <- "Species Pooled"

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")

  if( clearBadReps )
    for( s in 1:nS )
      for( p in 1:nP )
      {
        badIdx <- which(!goodReps_isp[,s,p])
        retroSB_itspt[badIdx,,s,p,] <- NA
        SB_ispt[badIdx,s,p,] <- NA
      }


  assErr_ispt <- array(NA, dim = c(nReps,nSS,nPP,nT) )
  for( s in 1:nSS )
    for( p in 1:nPP )
    {
      # Use first year's assessment for the historical assessment error...
      assErr_ispt[,s,p,1:(tMP-1)] <- (retroSB_itspt[1:nReps,1,s,p,1:(tMP-1)] - SB_ispt[1:nReps,s,p,1:(tMP-1)])/SB_ispt[1:nReps,s,p,1:(tMP-1)]
      for( tt in 1:pT )
      {
        # Now loop over projection years
        tIdx <- tMP + tt - 1
        assErr_ispt[,s,p,tIdx] <- (retroSB_itspt[1:nReps,tt,s,p,tIdx] - SB_ispt[1:nReps,s,p,tIdx])/SB_ispt[1:nReps,s,p,tIdx]
      }
    }

  if( save )
  {
    graphics.off()
    filename <- paste("assessmentError.pdf")
    outFile <- here::here("Outputs",groupFolder,obj$simLabel,filename)

    pdf( outFile, width = 11, height = 8 )
  }

  assErr_qspt <- apply( X = assErr_ispt,
                        FUN = quantile,
                        probs = c(0.025, 0.5, 0.975),
                        MARGIN = c(2,3,4),
                        na.rm = TRUE )

  years <- seq(from = fYear, by = 1, length.out = nT )

  if( proj )
    tIdx <- (tMP - 1):nT
  else
    tIdx <- 1:nT

  if( nReps > 3)
    traces <- sample( 1:nReps, 3)
  else traces <- c()

  # Now plot
  par(  mfcol = c(nPP, nSS ),
        mar = c(.1,1.5,.1,1.5),
        oma = c(3,3,2,2) )
  for( s in 1:nSS )
    for( p in 1:nPP )
    {
      plot( x = range(years[tIdx]),
            y = range(assErr_qspt[,s,p,tIdx],na.rm = T ),
            type = "n", xlab = "", ylab = "", axes = FALSE )
        mfg <- par("mfg")
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        
        axis( side = 2, las = 1 )
        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = stockNames[p], line = 1.5)
        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s])
        box()
        grid()

        # Plot baseline
        polygon(  x = c(years,rev(years)),
                  y = c(assErr_qspt[1,s,p,],rev(assErr_qspt[3,s,p,])),
                  border = NA, col = "grey60")
        lines( x = years, y = assErr_qspt[2,s,p,],
                col = "black", lwd = 2)
        for( traceIdx in traces )
          lines(  x = years, y = assErr_ispt[traceIdx,s,p,],
                  col = "black", lwd = .8 )

        abline( h = 0, lty = 2, lwd = 1 )
    }
  mtext(  side = 1, text = "Year", outer = TRUE, line = 2 )
  mtext(  side = 2, text = "Relative assessment error",
          line = 2, outer = TRUE )

  if(save)
    dev.off()

} # END plotTulipAssError

# Envelopes of simulated assessment errors
plotTulipRE_AM <- function( simNum = 1,
                            obj = NULL,
                            groupFolder = "DLSurveys_.3tau_Long",
                            save = FALSE,
                            proj = TRUE,
                            clearBadReps = FALSE,
                            abs = FALSE )
{

  if( is.null(obj))
  {
    .loadSim( simNum, groupFolder )
    obj <- blob
  }

  goodReps_isp  <- obj$goodReps_isp
  nReps         <- dim(goodReps_isp)[1]

  # Get biomass arrays
  SB_ispt        <- obj$om$SB_ispt
  VB_ispt        <- obj$om$vB_ispft[,,,2,]
  totB_ispt      <- obj$om$B_ispt
  retroSB_itspt  <- obj$mp$assess$retroSB_itspt


  ctlList <- obj$ctlList

  retroSB_itspt[retroSB_itspt < 0] <- NA
  # nReps     <- sum(goodReps)

  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT


  # Aggregate OM biomasses to match AM biomass
  if( ctlList$mp$data$spatialPooling )
  {
    newSB_ispt        <- apply( X = SB_ispt, FUN = sum, MARGIN = c(1,2,4), na.rm = T )
    newVB_ispt        <- apply( X = VB_ispt, FUN = sum, MARGIN = c(1,2,4), na.rm = T )
    newtotB_ispt      <- apply( X = totB_ispt, FUN = sum, MARGIN = c(1,2,4), na.rm = T )

    SB_ispt    <- SB_ispt[,,1,,drop = FALSE]
    SB_ispt[,,1,] <- newSB_ispt
    VB_ispt    <- VB_ispt[,,1,,drop = FALSE]
    VB_ispt[,,1,] <- newVB_ispt
    totB_ispt  <- totB_ispt[,,1,,drop = FALSE]
    totB_ispt[,,1,] <- newtotB_ispt     

  }

  if( ctlList$mp$data$speciesPooling )
  {
    newSB_ispt        <- apply( X = SB_ispt, FUN = sum, MARGIN = c(1,3,4), na.rm = T )
    newVB_ispt        <- apply( X = VB_ispt, FUN = sum, MARGIN = c(1,3,4), na.rm = T )
    newtotB_ispt      <- apply( X = totB_ispt, FUN = sum, MARGIN = c(1,3,4), na.rm = T )

    SB_ispt    <- SB_ispt[,1,,,drop = FALSE]
    SB_ispt[,1,,] <- newSB_ispt
    VB_ispt    <- VB_ispt[,1,,,drop = FALSE]
    VB_ispt[,1,,] <- newVB_ispt
    totB_ispt  <- totB_ispt[,1,,,drop = FALSE]
    totB_ispt[,1,,] <- newtotB_ispt

  }

  SB_ispt[SB_ispt == 0]     <- NA
  VB_ispt[VB_ispt == 0]     <- NA
  totB_ispt[totB_ispt == 0] <- NA

  nSS     <- dim( SB_ispt)[2]
  nPP     <- dim( SB_ispt)[3]

  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fYear         <- obj$ctlList$opMod$fYear
  pT            <- obj$ctlList$opMod$pT

  if( nPP == 1 )
    stockNames <- "Spatial Pooled"

  if( nSS == 1 )
    speciesNames <- "Species Pooled"

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")

  if( clearBadReps )
    for( s in 1:nS )
      for( p in 1:nP )
      {
        badIdx <- which(!goodReps_isp[,s,p])
        retroSB_itspt[badIdx,,s,p,] <- NA
        SB_ispt[badIdx,s,p,] <- NA
      }


  assErr_itspt <- array(NA, dim = c(nReps,pT,nSS,nPP,nT) )
  for( s in 1:nSS )
    for( p in 1:nPP )
    {
      for( tt in 1:pT )
      {
        # Now loop over projection years
        tIdx <- tMP + tt - 1
        assErr_itspt[,tt,s,p,] <- (retroSB_itspt[1:nReps,tt,s,p,] - SB_ispt[1:nReps,s,p,])/SB_ispt[1:nReps,s,p,]
      }
    }

  if( save )
  {
    graphics.off()
    filename <- paste("assessmentError.pdf")
    outFile <- here::here("Outputs",groupFolder,obj$simLabel,filename)

    pdf( outFile, width = 11, height = 8 )
  }

  absAssErr_itspt <- abs(assErr_itspt)
  MAREassErr_sp   <- apply( X = absAssErr_itspt, FUN = median, MARGIN = c(3,4), na.rm = T)
  MREassErr_sp    <- apply( X = assErr_itspt, FUN = median, MARGIN = c(3,4), na.rm = T)

  assErr_qspt <- apply( X = assErr_itspt,
                        FUN = quantile,
                        probs = c(0.025, 0.5, 0.975),
                        MARGIN = c(3,4,5),
                        na.rm = TRUE )


  years <- seq(from = fYear, by = 1, length.out = nT )

  if( proj )
    tIdx <- (tMP - 1):nT
  else
    tIdx <- 1:nT

  if( nReps > 3)
    traces <- sample( 1:nReps, 3)
  else traces <- c()

  # Now plot
  par(  mfcol = c(nPP, nSS ),
        mar = c(.1,1.5,.1,1.5),
        oma = c(3,3,2,2) )
  for( s in 1:nSS )
    for( p in 1:nPP )
    {
      plot( x = range(years[tIdx]),
            y = range(assErr_qspt[,s,p,tIdx], MAREassErr_sp[s,p],MREassErr_sp[s,p],na.rm = T ),
            type = "n", xlab = "", ylab = "", axes = FALSE )
        mfg <- par("mfg")
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        
        axis( side = 2, las = 1 )
        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = stockNames[p], line = 1.5)
        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s])
        box()
        grid()

        # Plot baseline
        polygon(  x = c(years,rev(years)),
                  y = c(assErr_qspt[1,s,p,],rev(assErr_qspt[3,s,p,])),
                  border = NA, col = "grey60")
        lines( x = years, y = assErr_qspt[2,s,p,],
                col = "black", lwd = 2)
        for( traceIdx in traces )
          lines(  x = years, y = assErr_itspt[traceIdx,pT,s,p,],
                  col = "black", lwd = .8 )

        abline( h = 0, lty = 2, lwd = 1 )
        abline( h = MAREassErr_sp[s,p], lty = 3, col = "red" )
        abline( h = MREassErr_sp[s,p], lty = 3, col = "salmon" )
    }
  mtext(  side = 1, text = "Year", outer = TRUE, line = 2 )
  mtext(  side = 2, text = "Relative assessment error",
          line = 2, outer = TRUE )

  if(save)
    dev.off()

} # END plotTulipAssError

# plotTulipBtCtBaseSim()
# Overlays the biomass can catch tulips
# from the baseline omniscient sim (black) and 
# stochastic simulation (red) so we can see
# how the different results are made
plotTulipBtCtBaseSim <- function( sim = 1, 
                                  groupFolder = "diffCV_fixedF_longGrid",
                                  var         = "SB_ispt",
                                  save        = FALSE,
                                  fYear       = 1956,
                                  proj        = TRUE,
                                  clearBadReps = FALSE )
{
  # Load loss reports
  objLoss <- .loadLoss(sim = sim, groupFolder)

  # Load blob, save reference points
  .loadSim(sim = sim, groupFolder)

  goodReps_isp <- blob$goodReps_isp

  rp <- blob$rp[[1]]

  scenName  <- blob$ctlList$ctl$scenarioName
  mpName    <- blob$ctlList$ctl$mpName

  stamp <- paste(scenName,":",mpName, sep = "")

  EmsyMSRefPts <- rp$EmsyMSRefPts
  FmsyRefPts   <- rp$FmsyRefPts


  simFolder <- blob$simLabel

  speciesNames    <- objLoss$speciesNames
  stockNames      <- objLoss$stockNames

  baseState_ispt  <- objLoss$baseStates[[var]]
  simState_ispt   <- objLoss$simStates[[var]]

  B0_isp          <- objLoss$baseStates[["SB_ispt"]][,,,1]


  # Pull model dimensions
  tMP <- objLoss$tMP
  nT  <- objLoss$nT
  nF  <- objLoss$nF
  nS  <- objLoss$nS
  nP  <- objLoss$nP
  pT  <- objLoss$pT

  if( proj )
    tIdx <- (tMP-1):nT
  else tIdx <- 1:nT

  years <- seq( from = fYear, by = 1, length.out = nT )

  nReps <- dim(simState_ispt)[1]
  
  traces <- sample(1:nReps, size = 3)

  # BeqFmsySS_sp <- array(NA, dim =c(nS+1, nP+1))
  # BeqEmsyMS_sp <- array(NA, dim =c(nS+1, nP+1))

  # BeqFmsySS_sp[1:nS,1:nP] <- FmsyRefPts$BeqFmsy_sp
  # BeqFmsySSS_sp[1:nS,p+1] <- apply( X = FmsyRefPts$BeqFmsy_sp, 
  #                                   FUN =  )
  
  if( save )
  {
    graphics.off()
    filename <- paste("baseSimTulipOverlay_",var,".pdf",sep = "")
    outFile <- file.path("./Outputs",groupFolder,simFolder,filename)

    pdf( outFile, width = 11, height = 8 )
  }

  # put on depletion scale
  if( var == "SB_ispt")
  {
    maxY_sp <- matrix(1, nrow = nS + 1, ncol = nP +1 )
    yLab = expression(paste("Spawning biomass relative to unfished (", B[t]/B[0], ")",sep = ""))
    for( p in 1:(nP+1) )  
    {
      
      for( s in 1:(nS+1) )
        for( i in 1:nReps)
        {
          baseState_ispt[i,s,p,]  <- baseState_ispt[i,s,p,] / B0_isp[i,s,p]
          simState_ispt[i,s,p,]   <- simState_ispt[i,s,p,] / B0_isp[i,s,p]
        }

      maxY_sp[,p] <- max( baseState_ispt[,,p,tIdx], simState_ispt[,,p,tIdx])
    }

    
  } else {
    yLab    <- "Catch (kt)"
    maxY_sp <- matrix(NA, nrow = nS+1, ncol = nP+1 )

    for( s in 1:(nS+1))
      for( p in 1:(nP+1))
        maxY_sp[s,p] <- max( baseState_ispt[,s,p,tIdx], simState_ispt[,s,p,tIdx], na.rm = T )
  }

  if( clearBadReps )
    for( s in 1:nS )
      for( p in 1:nP )
      {
        badIdx <- which(!goodReps_isp[,s,p])
        baseState_ispt[badIdx,s,p,] <- NA
        simState_ispt[badIdx,s,p,] <- NA
      }

  baseState_qspt <- apply( X = baseState_ispt, FUN = quantile,
                            probs = c(0.025, 0.5, 0.975 ),
                            MARGIN = c(2,3,4), na.rm = T )

  simState_qspt <- apply( X = simState_ispt, FUN = quantile,
                            probs = c(0.025, 0.5, 0.975 ),
                            MARGIN = c(2,3,4), na.rm = T )

  par( mfcol = c((nP+1), nS+1),
        mar = c(0.1,0.1,0.1,0.1),
        oma = c(4,4.5,3,3) )

  if( var == "C_ispt" )
    par( mar = c(.1, 1.5, .1, 1.5) )


  basePolyCol <- "grey50"
  baseLineCol <- "black"
  simPolyCol <- scales::alpha("red",.3)
  simLineCol <- "red"

  for( s in 1:(nS+1) )
    for( p in 1:(nP+1) )
    {
      plot( x = range(years[tIdx]),
            y = c(0, maxY_sp[s,p] ),
            type = "n", xlab = "", ylab = "", axes = FALSE )
        mfg <- par("mfg")
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        if( mfg[2] == 1 )
          axis( side = 2, las = 1 )
        if( var == "C_ispt" & mfg[2] != 1 )
          axis( side = 2, las = 1)
        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = stockNames[p], line = 1.5)
        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s])
        box()
        grid()

        # Plot baseline
        polygon(  x = c(years,rev(years)),
                  y = c(baseState_qspt[1,s,p,],rev(baseState_qspt[3,s,p,])),
                  border = NA, col = basePolyCol)
        lines( x = years, y = baseState_qspt[2,s,p,],
                col = baseLineCol, lwd = 2)
        for( traceIdx in traces )
          lines(  x = years, y = baseState_ispt[traceIdx,s,p,],
                  col = baseLineCol, lwd = .8 )

        # Plot stochastic sim
        polygon(  x = c(years,rev(years)),
                  y = c(simState_qspt[1,s,p,],rev(simState_qspt[3,s,p,])),
                  border = NA, col = simPolyCol)
        lines( x = years, y = simState_qspt[2,s,p,],
                col = simLineCol, lwd = 2)
        for( traceIdx in traces )
          lines(  x = years, y = simState_ispt[traceIdx,s,p,],
                  col = simLineCol, lwd = .8 )

        # Plot some lines
        abline( v = years[tMP], lty = 2, lwd = 1 )
        if( s <= nS & p <= nP & var == "SB_ispt")
        {
          abline( h = FmsyRefPts$BeqFmsy_sp[s,p]/B0_isp[1,s,p], 
                  lty = 3, lwd = 1 )
          abline( h = EmsyMSRefPts$BeqEmsy_sp[s,p]/B0_isp[1,s,p], 
                  lty = 4, lwd = 1 )
        }


    }
  mtext( side = 1, text = "Year", outer = TRUE, line = 2 )
  mtext( side = 2, text = yLab, outer = TRUE, line = 2 ) 

  mtext( side = 1, text = stamp, outer = TRUE, line = 2.5,
          adj = .85, cex = .6, col = "grey75")

  if(save)
    dev.off()

}



# plotBatchConvergenceRate()
# Abandoned performance metric for simulations,
# basically measuring the usefulness of a model
# under lower data quality conditions - can be
# used to set a threshold for including
# results in the paper...
plotBatchConvergenceRate <- function( groupFolder = "DERTACS_reruns_sep24",
                                      prefix = "parBat",
                                      AMlabs = c( SS = "singleStock",
                                                  HMS = "hierMultiStock",
                                                  SpeP = "speciesPooling",
                                                  SpaP = "spatialPooling",
                                                  TP = "totalAgg" ),
                                      scenLabs = c( HMAS  = "DERfit_HcMcAsSsIdx",
                                                    MAS   = "DERfit_McAsSsIdx",
                                                    AS    = "DERfit_AsSsIdx",
                                                    MS    = "DERfit_McSsIdx" ) )
{
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
                filter( grepl( prefix, simLabel ))

  # Break up MP names into factor levels
  # MP labels are AM_Fsrce_eqbm
  splitMP <- function(  mpLab, 
                        breaks = "_",
                        factorNames = c("AM","Fsrce","eqbm") )
  {
    splitMP <- stringr::str_split(mpLab, breaks)[[1]]

    outList <- vector(  mode = "list", 
                        length = length(splitMP) )
    names(outList) <- factorNames
    for( k in 1:length(splitMP))
      outList[[k]] <- splitMP[k]

    outList
  }

  MPfactors <- lapply( X = info.df$mp, FUN = splitMP )
  MPfactors <- do.call(rbind, MPfactors)

  # Read in the performance tables
  mpTable <- cbind( info.df, MPfactors ) %>%
              mutate( perfPath = here::here("Outputs",groupFolder,simLabel,"simPerfStats.csv") )

  
  perfTableList <- lapply(  X = mpTable$perfPath, 
                            FUN = read.csv,
                            header = TRUE,
                            stringsAsFactors = FALSE )

  names( perfTableList ) <- mpTable$simLabel

  # Summarise the perf tables into convergence rates
  summariseTable <- function( table )
  {
    table <- table %>%
              group_by(simLabel) %>%
              summarise(  obsCVmult = mean(projObsErrMult,na.rm = T ),
                          minConv   = min(pGoodReps,na.rm = T),
                          maxConv   = max(pGoodReps,na.rm = T),
                          meanConv  = mean(pGoodReps,na.rm = T) )

    table
  }

  # Summarise and make table
  perfTableSummList <- lapply(  X = perfTableList, 
                                FUN = summariseTable )

  convRateTable <- do.call(rbind, perfTableSummList )


  mpTable <- mpTable %>% left_join( convRateTable, 
                                    by = "simLabel" )

  CVmults <- unique(mpTable$obsCVmult)
  CVmults <- CVmults[order(CVmults)]

  scenarios   <- unique(mpTable$scenario)

  AMs       <- AMlabs
  AMcols    <- RColorBrewer::brewer.pal(length(AMs), "Set2")
  AMjitter  <- seq( from = -.3, to = .3, length.out = length(AMs) )

  AMwidth <- AMjitter[2] - AMjitter[1]

  par( mar = c(3,3,1,1))

  # Now set up plotting area
  plot( x = c(0,length(scenarios) + 1), y = c(0,1),
        xlab = "", ylab = "", axes=FALSE,
        type = "n" )
    axis( side = 1, at = 1:length(scenLabs),
          labels = names(scenLabs) )
    axis( side = 2, las = 1  )
    box()
    for( aIdx in 1:length(AMs) )
    {
      amLab <- AMs[aIdx]
      SStable <- mpTable %>% filter( AM == amLab,
                                      grepl("SS", eqbm) ) %>%
                              dplyr::select( scenario,
                                              meanConv )


      MStable <- mpTable %>% filter( AM == amLab,
                                      grepl("MS", eqbm) ) %>%
                              dplyr::select( scenario,
                                              meanConv )

      for( sIdx in 1:length(scenLabs) )
      {
        if( nrow(SStable) > 0)
          rect( xleft = sIdx + AMjitter[aIdx] - AMwidth/2,
                xright = sIdx + AMjitter[aIdx],
                ybottom = 0,
                ytop = SStable[SStable$scenario == scenLabs[sIdx], ]$meanConv,
                border = "black",
                col = AMcols[aIdx] )
        if( nrow(MStable) > 0 )
          rect( xleft = sIdx + AMjitter[aIdx],
                xright = sIdx + AMjitter[aIdx] + AMwidth/2,
                ybottom = 0,
                ytop = MStable[MStable$scenario == scenLabs[sIdx], ]$meanConv,
                border = "black",
                col = AMcols[aIdx] )
      }
      

      
    }

    legend( x = "topleft", bty = "n",
            legend = names(AMlabs),
            col = "black",
            pch = 22,
            pt.cex = 2,
            pt.bg = AMcols )

} # END plotBatchConvergenceRate()


# plot
plotSensBatchSummary_SAisp <- function( groupFolder = "sensRuns_UmsyCV_Jan5",
                                        prefix = "parBat",
                                        lossType = "abs",
                                        var = "C_ispt",
                                        period = 73:82,
                                        lossList = NULL,
                                        refPts = NULL,
                                        AMlabs = c( SS = "singleStock",
                                                    HMS = "hierMultiStock",
                                                    SpePool = "speciesPooling",
                                                    SpaPool = "spatialPooling",
                                                    TA = "totalAgg" ),
                                        scenLabs = c( High_hierSD0.1  = "High_shrinkSD0.1",
                                                      High_hierSD0.2  = "High_shrinkSD0.2",
                                                      High_hierSD0.5  = "High_shrinkSD0.5"),
                                        clearBadReps = FALSE,
                                        scenSplitString = "shrinkSD",
                                        xlab = "Bmsy CV",
                                        noPar = FALSE,
                                        plotPts = FALSE,
                                        printLeg = TRUE
                                      )
{
  # browser()
  # First, make cumulative loss array
  cumLossList <- makeCumulativeLossArray_SAisp( groupFolder = groupFolder,
                                                  prefix = prefix,
                                                  lossType = lossType,
                                                  var = var,
                                                  period = period,
                                                  lossList = lossList,
                                                  refPts = refPts,
                                                  AMlabs = AMlabs,
                                                  scenLabs = scenLabs,
                                                  clearBadReps = clearBadReps )


  cumLoss_SAisp <- cumLossList$cumLossArray_SAisp
  if(is.null(lossList))
    lossList      <- cumLossList$lossList
  info.df <- cumLossList$info.df


  # Calc median cumulative loss
  medCumLoss_SAsp <- apply( X = cumLoss_SAisp,
                            FUN = median, MARGIN = c(1,2,4,5),
                            na.rm = TRUE )


  # Can I append stuff to the data.frame? It will
  # make regressions easier later.
  meltedArray.df <- reshape2::melt(medCumLoss_SAsp) %>%
                      filter( species %in% c("Dover","English","Rock"),
                              stock %in% c("HSHG","QCS","WCVI"),
                              scenario %in% scenLabs )

  scenCpts  <- lapply(  X = meltedArray.df$scenario, FUN = splitMP,
                        breaks = "_", factorNames = c("dataScen","CV") )
  scenCpts  <- do.call(rbind, scenCpts) %>% as.data.frame()
  scenCpts$scenario <- meltedArray.df$scenario
  CVs       <- lapply(  X = scenCpts$CV, FUN = splitMP,
                        breaks = scenSplitString, 
                        factorNames = c("blank","CV") )
  CVs       <-  do.call( rbind, CVs) %>% as.data.frame()
  CVs$scenario <- scenCpts$scenario
  scenCpts$CV  <- NULL

  meltedArray.df$dataScen <- as.character(scenCpts$dataScen)
  meltedArray.df$CV       <- as.numeric(CVs$CV)

  AMcols          <- RColorBrewer::brewer.pal(length(AMlabs), "Dark2")
  names(AMcols)   <- AMlabs
  AMpch           <- 20 + 1:length(AMlabs)
  names(AMpch)    <- AMlabs
  AMlty           <- 1:length(AMlabs)
  names(AMlty)    <- AMlabs

  dataScens       <- unique(meltedArray.df$dataScen)
  nScen           <- length(dataScens)
  scenLty         <- 1:nScen

  species         <- unique(meltedArray.df$species)
  stock           <- unique(meltedArray.df$stock)

  batchAMs        <- unique(meltedArray.df$AM)


  nS <- length(species)
  nP <- length(stock)

  # Make lists to hold regressions
  scenRegList       <- vector(mode = "list", length = nScen)

  stockSpecRegList  <- vector(mode = "list", length = nS * nP)
  
  # Now, group and subtract mean of each species
  # response values
  meltedArray.df <- meltedArray.df %>%
                    group_by( species, stock, dataScen ) %>%
                    mutate( stdResp = (value - mean(value))/sd(value) ) %>%
                    ungroup()

  plotAMidx <- which(AMlabs %in% batchAMs)

  xJitter <- 0
  if(length(AMlabs) > 1)
    xJitter <- seq(from = -.05, to = .05, length.out = length(AMlabs) )

  if(!noPar)
    par( mar = c(2,2,1,1), oma = c(3,3,1,1))

  yRange <- range(meltedArray.df$stdResp, na.rm = T)

  plot( x = c(0,1.1*max(meltedArray.df$CV, na.rm = T) ),
        y = c(yRange[1],1.4 * yRange[2]),
        xlab = "",
        ylab = "", axes = FALSE, type = "n" )
    mfg <- par("mfg")
    
    if( mfg[1] == mfg[3])
    {
      axis(side = 1)
      mtext( side = 1, text = xlab, outer = FALSE, line = 3)
    }

    axis( side = 2, las = 1)

    grid()
    abline( h = 0, lty = 2, lwd = 1 )

    box()

    for(amIdx in plotAMidx )
    {
      # Subset dataframe
      subDF <- meltedArray.df %>%
                filter( AM == AMlabs[amIdx] )

      # First, plot all the points
      if(plotPts)
        points( x = subDF$CV + xJitter[amIdx], y = subDF$stdResp,
                col = AMcols[AMlabs[amIdx]],
                bg = "NA", pch = AMpch[AMlabs[amIdx]] )
    }


    amRegList <- vector(mode = "list", length = length(plotAMidx))
    for( i in 1:length(plotAMidx))
    {
      amIdx <- plotAMidx[i]
      # Subset dataframe
      subDF <- meltedArray.df %>%
                filter( AM == AMlabs[amIdx] )
      # Now loop over scenarios
      for(scenIdx in 1:nScen)
      {
        scenID <- dataScens[scenIdx]

        scenSubDF <- subDF %>% filter( dataScen == scenID )

        # Make scenario average regressions
        scenRegList[[scenIdx]]$lm <- lm( stdResp ~ CV, data = scenSubDF )

        scenRegList[[scenIdx]]$stockSpecRegList <- stockSpecRegList

        listIdx <- 0
        for(s in 1:nS )
          for(p in 1:nP )
          {
            listIdx <- listIdx + 1
            specID <- species[s]
            stockID <- stock[p]

            stockSpecDF <- scenSubDF %>%
                            filter( stock == stockID,
                                    species == specID )
            
            names(scenRegList[[scenIdx]]$stockSpecRegList)[listIdx] <- paste(specID,"_",stockID,sep = "")
            scenRegList[[scenIdx]]$stockSpecRegList[[listIdx]] <- lm( stdResp ~ CV, data = stockSpecDF )


            # abline( reg = scenRegList[[scenIdx]]$stockSpecRegList[[listIdx]],
            #         col = AMcols[AMlabs[amIdx]],
            #         lty = scenLty[scenIdx], lwd = .5 )
          }

        summSubDF <- scenSubDF %>%
                      group_by(CV) %>%
                      summarise( meanResp = mean(stdResp) ) 

        summSubDF$predResp <- predict(object = scenRegList[[scenIdx]]$lm, newdata = summSubDF)

        lines( x = summSubDF$CV, y = summSubDF$predResp,
                col = AMcols[AMlabs[amIdx]],
                lty = AMlty[AMlabs[amIdx]], lwd = 2)
        points( x = summSubDF$CV + xJitter[amIdx], y = summSubDF$meanResp,
                col = AMcols[AMlabs[amIdx]],
                bg = "NA", pch = AMpch[AMlabs[amIdx]])


        # abline( reg = scenRegList[[scenIdx]]$lm,
        #         col = AMcols[AMlabs[amIdx]],
        #         lty = AMlty[AMlabs[amIdx]], lwd = 2 )

        amRegList[[i]] <- scenRegList

      }

    }
    amRegPval <- numeric(length = length(plotAMidx))
    # Let's put the significance of each change
    for( i in 1:length(plotAMidx) )
    {
      amRegPval[i] <- coef(summary(amRegList[[i]][[1]]$lm))[2,4]
    }

    legTxt <- paste("p = ", round(amRegPval,2), sep = "")

    if( printLeg )
    {
      legend( x = "topleft",
              bty = "n",
              legend = legTxt,
              lty = AMlty[AMlabs[plotAMidx]],
              col = AMcols[AMlabs[plotAMidx]], 
              pch = AMpch[AMlabs[plotAMidx]],
              lwd = 2 )
    }
}

# Break up MP names into factor levels
# MP labels are AM_Fsrce_eqbm
splitMP <- function(  mpLab, 
                      breaks = "_",
                      factorNames = c("AM","Fsrce","eqbm") )
{
  splitMP <- stringr::str_split(mpLab, breaks)[[1]]

  outList <- vector(  mode = "list", 
                      length = length(splitMP) )
  names(outList) <- factorNames
  for( k in 1:length(splitMP))
    outList[[k]] <- splitMP[k]

  outList
} # END splitMP()

# plotRankDists()
# Multi-panel plot of relative and absolute
# loss under a batch of MPs for all stock/species
# combinations, including data-pooled and 
# coast-wide aggregations.
plotRankDists_sp <- function( groupFolder = "DERTACS_reruns_sep24",
                              prefix = "parBat",
                              lossType = "abs",
                              var = "C_ispt",
                              period = 73:82,
                              lossList = NULL,
                              dim1 = 1:3,   # species (D,E,R, DP)
                              dim2 = 1:3,   # stocks (H,Q,W, CW)
                              qProbs = c(.025,.5,.975),
                              refPts = "MSrefPts",
                              AMlabs = rev(c( SS = "singleStock",
                                              HMS = "hierMultiStock",
                                              SpecPool = "speciesPooling",
                                              SpatPool = "spatialPooling",
                                              TotPool = "totalAgg" ) ),
                              scenLabs = c( Rich  = "DERfit_HcMcAsSsIdx",
                                            Mod  = "DERfit_McAsSsIdx",
                                            Poor = "DERfit_AsSsIdx" ),
                              specLabs = c("Dover","English","Rock"),
                              stockLabs = c("HSHG", "QCS", "WCVI"),
                              clearBadReps = TRUE,
                              minSampSize = 100,
                              rotxlabs = 0,
                              xlab = TRUE,
                              vertLines = TRUE  )
{

  # First, calculate the paired loss ranks

  rankArray_SAisp <- calcPairedLossRank(  groupFolder = groupFolder,
                                          lossVar = var,
                                          lossType = lossType,
                                          prefix = prefix,
                                          period = period,
                                          clearBadReps = clearBadReps,
                                          minSampSize = minSampSize )

  goodReps_SAisp <- attr(rankArray_SAisp,"conv")

  nScen <- length(scenLabs)
  nAM   <- length(AMlabs)
  nS    <- length(dim1)
  nP    <- length(dim2)



  xJitter   <- (1:5) * 0.2 - 0.1
  xRadius   <- 0.1 

  # Now, make the plotting region
  par( mfcol = c(nScen, nAM), mar = c(.2,.2,.2,.2), oma = c(5,5,3,3) )

  for( amIdx in 1:nAM )
  {
    for( scenIdx in 1:nScen )
    {
      AMlab <- AMlabs[amIdx]
      scenLab <- scenLabs[scenIdx]
      # Make a plotting window, we want to break it up into nSxnP
      # cells
      plot( x = c(0,nS), y = c(0,nP),
            xaxs = "i", yaxs = "i", type = "n", axes = FALSE )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis( side = 1, at = 0.5 + 1:nS - 1, labels = specLabs  )

      # reverse y axis
      if(mfg[2] == 1)
        axis( side = 2, at = 0.5 + 1:nP - 1, las = 1, labels = rev(stockLabs) )

      if(mfg[2] == mfg[4])
        rmtext( outer = TRUE, line = 0.1, txt = names(scenLabs)[scenIdx], font = 2, cex = 1.5 )

      if(mfg[1] == 1)
        mtext( line = 1, text = names(AMlabs)[amIdx], font = 2, side = 3 )

      # Close up cells
      grid()
      abline( v = 1:nP )
      abline( h = 1:nS )
      box()

      for( s in 1:nS )
        for( p in 1:nP )
        {
          convReps <- goodReps_SAisp[scenLab,AMlab,,s,p]
          nConv <- sum(convReps)
          rank.df <- reshape2::melt(rankArray_SAisp[scenLab,AMlab, convReps, s, p]) %>%
                      mutate( rank = as.character(value) ) %>%
                      dplyr::group_by( rank ) %>%
                      dplyr::summarise( n = n() )

          totReps   <- sum( rank.df$n )
          rankProps <- rank.df$n/totReps
          names(rankProps) <- rank.df$rank
          whichMode <- which.max(rankProps)
          avgRank   <- sum( rankProps * as.numeric(rank.df$rank) )

          cols <- rep("grey60",nAM)
          cols[whichMode] <- "darkred"

          rect( xleft   = xJitter + s - 1 - xRadius,
                xright  = xJitter + s - 1 + xRadius,
                ybottom = 0 + (nP - p + 1) - 1,
                ytop    = (nP - p + 1) - 1 + rankProps,
                col = cols, border = "black" )

          text( x = 0.3 + s - 1, y = 0.9 + (nP - p + 1) - 1,
                labels = nConv,
                font = 2)


        }

    }
  } 

} # END plotRankDists_sp()

# plotRankDists()
# Multi-panel plot of relative and absolute
# loss under a batch of MPs for all stock/species
# combinations, including data-pooled and 
# coast-wide aggregations.
plotRankDists <- function(  groupFolder = "DERTACS_reruns_sep24",
                            prefix = "parBat",
                            lossType = "abs",
                            var = "C_ispt",
                            period = 73:82,
                            lossList = NULL,
                            dim1 = 1:3,   # species (D,E,R, DP)
                            dim2 = 1:3,   # stocks (H,Q,W, CW)
                            qProbs = c(.025,.5,.975),
                            refPts = "MSrefPts",
                            AMlabs = rev(c( SS = "singleStock",
                                            HMS = "hierMultiStock",
                                            SpecPool = "speciesPooling",
                                            SpatPool = "spatialPooling",
                                            TotPool = "totalAgg" ) ),
                            scenLabs = c( Rich  = "DERfit_HcMcAsSsIdx",
                                          Mod  = "DERfit_McAsSsIdx",
                                          Poor = "DERfit_AsSsIdx" ),
                            clearBadReps = TRUE,
                            minSampSize = 100,
                            rotxlabs = 0,
                            xlab = TRUE,
                            vertLines = TRUE,
                            nGoodReps = 100 )
{

  # First, calculate the paired loss ranks

  rankArray_SAisp <- calcPairedLossRank(  groupFolder = groupFolder,
                                          lossVar = var,
                                          lossType = lossType,
                                          prefix = prefix,
                                          period = period,
                                          clearBadReps = clearBadReps,
                                          minSampSize = minSampSize )
  goodReps_SAisp <- attr(rankArray_SAisp,"conv")[scenLabs,,,,,drop = FALSE]
  
  # Cut down to scenarios we want
  rankArray_SAisp <- rankArray_SAisp[scenLabs,,,,,drop = FALSE]
  goodReps_SAisp <- goodReps_SAisp[scenLabs,,,,,drop = FALSE]



  nScen <- length(scenLabs)
  nAM   <- length(AMlabs)
  nS    <- length(dim1)
  nP    <- length(dim2)


  if( clearBadReps )
  {
    nGood_SAsp <- apply( X = goodReps_SAisp, FUN = sum, MARGIN = c(1,2,4,5),
                          na.rm = T )

    badRuns_arr.ind <-  which(nGood_SAsp < minSampSize, arr.ind = TRUE )

    if( nrow(badRuns_arr.ind) > 0 )
      for( j in 1:nrow(badRuns_arr.ind))
      {
        ind <- badRuns_arr.ind[j,]
        # Remove the species and stock results from all AMs - replace goodReps with TRUE
        # as the NAs in rank will be removed later
        goodReps_SAisp[,,,ind[3],ind[4]] <- NA
        rankArray_SAisp[,,,ind[3],ind[4]] <- NA
      }

    # Ok, now  we want to calculate the largest number of reps we have available
    goodReps_SAi <- apply(X = goodReps_SAisp, FUN = prod, MARGIN = c(1,2,3), na.rm = T)
    goodReps_Si <- apply(X = goodReps_SAi, FUN = prod, MARGIN = c(1,3), na.rm = T)

    goodReps_i <- apply(X = goodReps_SAisp, FUN = prod, MARGIN = 3, na.rm = T)

    goodIdx <- which(goodReps_i == 1)

    if( length(goodIdx) > nGoodReps)
      goodIdx <- goodIdx[1:nGoodReps]

    
    goodReps_SAisp  <- goodReps_SAisp[,,goodIdx,,]
    goodReps_SAi    <- goodReps_SAi[,,goodIdx]
    rankArray_SAisp <- rankArray_SAisp[,,goodIdx,,]

  }



  xJitter   <- (1:5)
  xRadius   <- 0.5

  # Now, make the plotting region
  par( mfcol = c(nScen, nAM), mar = c(.4,.2,.25,.2), oma = c(5,5,3,3) )

  for( amIdx in 1:nAM )
  {
    for( scenIdx in 1:nScen )
    {
      AMlab <- AMlabs[amIdx]
      scenLab <- scenLabs[scenIdx]
      # Make a plotting window, we want to break it up into nSxnP
      # cells
      plot( x = c(0,6), y = c(0,1),
            xaxs = "i", yaxs = "i",
            type = "n", axes = FALSE )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3])
        axis( side = 1, at = 1:5  )

      # # reverse y axis on left margin
      if(mfg[2] == 1)
        axis( side = 2, las = 1 )

      if(mfg[2] == mfg[4])
        rmtext( outer = TRUE, line = 0.1, txt = names(scenLabs)[scenIdx], font = 2, cex = 1.5 )

      if(mfg[1] == 1)
        mtext( line = 1, text = names(AMlabs)[amIdx], font = 2, side = 3 )

      # Close up cells
      grid()
      box()

      rank.df <- reshape2::melt(rankArray_SAisp[scenLab,AMlab,, , ]) %>%
                  filter(!is.na(value)) %>%
                  mutate( rank = as.character(value) ) %>%
                  dplyr::group_by( rank ) %>%
                  dplyr::summarise( n = n() )


      totReps   <- sum( rank.df$n )
      rankProps <- rank.df$n/totReps
      names(rankProps) <- rank.df$rank
      whichMode <- which.max(rankProps)
      avgRank   <- sum( rankProps * as.numeric(rank.df$rank) )

      cols <- rep("grey60",nAM)
      cols[whichMode] <- "darkred"

      rect( xleft   = xJitter - xRadius,
            xright  = xJitter + xRadius,
            ybottom = 0,
            ytop    = rankProps,
            col = cols, border = "black" )
      segments(x0 = avgRank, y0 = 0, y1 = rankProps[whichMode], lty = 2, col = "black", lwd = 3)


    }
  } 
  mtext( side = 1, text = "Rank", outer = TRUE, line = 2 )
  mtext( side = 2, text = "Proportion", outer = TRUE, line = 2.5 )

} # END plotRankDists()\

plotSpecStockBoxes <- function( species, stock )
{
  nS <- length(species)
  nP <- length(stock)

  plot( x = c(0.5,nS+.5), y = c(.5,nP+.5), type = "n", axes = FALSE,
        xaxs = "i", yaxs = "i" )
    axis( side = 1, at = 1:nS, labels = species )
    axis( side = 2, at = 1:nP, labels = rev(stock), las = 1 )
    box(lwd = 3)
    for( s in 1:nS )
      for( p in 1:nP )
      {
        rect( xleft = s - .3,
              xright = s + .3,
              ybottom = p - .3,
              ytop = p + .3, col = "grey60" )
        # text( x = s, y = p, label =  )
      }
}

# plotAMschematic()
# Plots a 3x2 of each AM structure
plotAMschematic <- function(  species = c("Dover","English","Rock"),
                              stock = c("HSHG","QCS","WCVI") )
{
  nS <- length(species)
  nP <- length(stock)

  par(mfrow = c(3,2),
      mar = c(2,3,2,1),
      oma = c(4,4,3,3) )

  # Total Pooling
  plotSpecStockBoxes(species,stock)
    mtext( side = 3, text = "Total Pooling")
    # rect( xleft = 0.6, xright = nS+.4,
    #       ybottom = 0.6, ytop = nP+.4,
    #       border = "black")

  # Spatial Pooling
  plotSpecStockBoxes(species,stock)
    mtext( side = 3, text = "Spatial Pooling")
    abline( h = 1:nS + 0.5, lwd = 3)

  # Species Pooling
  plotSpecStockBoxes(species,stock)
    mtext( side = 3, text = "Species Pooling")
    abline( v = 1:nP + 0.5, lwd = 3)

  # Hierarchical Multi-stock
  plotSpecStockBoxes(species,stock)
    mtext( side = 3, text = "Hierarchical Multi-stock")
    abline( h = 1:nS + 0.5, lty = 2, lwd = 3)
    abline( v = 1:nP + 0.5, lty = 2, lwd = 3)

  # Single Stock
  plotSpecStockBoxes(species,stock)
    mtext( side = 3, text = "Single stock")
    abline( v = 1:nP + 0.5, lwd = 3)
    abline( h = 1:nS + 0.5, lwd = 3)

  mtext( outer = TRUE, side = 1, line = 2, text = "Species")
  mtext( outer = TRUE, side = 2, line = 2, text = "Stock Area")

} # END plotAMschematic

# plotDataSummary()
# Data summary for each species and stock, uses the
# hierSCAL report object
plotDataSummary <- function(  obj,
                              fleetLabs = c(  "Historical",
                                              "Modern",
                                              "HS Ass.",
                                              "Synoptic"),
                              dataCells = c(  "Indices",
                                              "Catch") )
{
  reports <- obj$ctlList$opMod$histRpt

  datObj <- reports
  repObj <- reports
  
  # Pull data
  I_spft      <- datObj$I_spft
  age_axspft  <- datObj$age_axspft
  len_lxspft  <- datObj$len_lxspft
  C_spft      <- datObj$C_spft
  aal_table   <- datObj$aal_table


  # Model dims
  nS      <- repObj$nS
  nP      <- repObj$nP
  nF      <- repObj$nF
  nT      <- repObj$nT

  # Year sequence for x axis
  fYear <- blob$ctlList$opMod$fYear
  years <- seq(from = fYear, by = 1, length.out = nT)


  yLabs <- rev(dataCells)
  specLabs <- blob$ctlList$opMod$species
  stockLabs <- blob$ctlList$opMod$stocks



  fleetCols <- RColorBrewer::brewer.pal(nF,"Dark2")

  checkPosObs <- function( arr, out = 16 )
  {

    if(any(arr > 0) )
      return(out)
    else return(NA)
  }

  par( mfcol = c(nS,nP),
        mar = c(.1,.1,.1,.1),
        oma = c(5,6,2,2) )


  # Prepare catch data
  
  C_spt <- apply( X = C_spft, FUN = sum, MARGIN = c(1,2,4), na.rm = T)
  C_sp.max <- apply(C_spt, FUN = max, MARGIN = c(1,2))
  C_spt <- sweep( x = C_spt, MARGIN = c(1,2), C_sp.max, FUN = "/")

  fleetJitter <- seq( from = -.8, to = .8, 
                      length.out = nF )

  nCells <- length(dataCells)
  cellMins <- seq(from = 0, by = 2, length = nCells)
  names(cellMins) <- rev(dataCells)

  # Prepare pch

  for( s in 1:nS )
    for( p in 1:nP )
    {
      plot( x = range(years), y = c(0,2*length(dataCells)),
            type = "n", axes = FALSE, yaxs = "i" )
        # Axes and labels
        mfg <- par( "mfg")
        if( mfg[1] == mfg[3] )
          axis( side = 1 )

        if( mfg[2] == 1 )
          axis( side = 2, las = 1,
                at = cellMins + 1,
                labels = yLabs )

        if( mfg[1] == 1 )
          mtext( side = 3, text = specLabs[s], font = 2)

        if( mfg[2] == mfg[4] )
          rmtext( txt = stockLabs[p], font = 2, line = .05, outer = TRUE,
                  cex = 1.5 )

        grid()
        box()
        
        # Break up window with hz lines
        abline( h = seq(from = 2, by = 2, length = nCells-1), lwd = .8 )

        # Plot catch
        if("Catch" %in% dataCells)
        {
        rect( xleft = years - .3, xright = years + .3,
              ybottom = cellMins["Catch"], ytop = 1.5*C_spt[s,p,],
              col = "grey50", border = NA )
        }

        if("AAL" %in% dataCells)
        {
          # Plot AAL
          # Prepare table
          aal_table.this <- aal_table %>%
                            as.data.frame() %>%
                            filter( species == s, stock == p )%>%
                            group_by( year, fleet ) %>%
                            summarise( posSamp = 1 )
          points( x =  years[aal_table.this$year],
                  y = cellMins["AAL"] + 1 + fleetJitter[aal_table.this$fleet],
                  pch = 16, col = fleetCols[aal_table.this$fleet] )
        }

        # Plot Lengths
        if("Lengths" %in% dataCells)
        {
          len_spft <- apply(X = len_lxspft, FUN = checkPosObs, MARGIN = c(3:6) )
          for(f in 1:nF )
          {
            plotPts <- which(!is.na(len_spft[s,p,f,]))
            nPts <- length(plotPts)
            points( x = years[plotPts], y = rep(cellMins["Lengths"] + 1 + fleetJitter[f],nPts),
                    col = fleetCols[f], pch = len_spft[s,p,f,plotPts] )
          }
        }

        # Plot Ages
        if("Ages" %in% dataCells )
        {    
          age_spft <- apply(X = age_axspft, FUN = checkPosObs, MARGIN = c(3:6) )
          for(f in 1:nF )
            points( x = years, y = rep(cellMins["Ages"] + 1 + fleetJitter[f],nT),
                    col = fleetCols[f], pch = age_spft[s,p,f,] )
        }

        if( "Indices" %in% dataCells )
        {
          # Plot indices
          I_spft[I_spft > 0] <- 16
          I_spft[I_spft <= 0] <- NA
          for( f in 1:nF )
            points( x = years, y = rep(cellMins["Indices"] + 1 + fleetJitter[f],nT),
                    col = fleetCols[f], pch = I_spft[s,p,f,] )
        }

    }

    par(xpd=TRUE, oma = c(0,0,0,0), mfcol = c(1,1))
    # plot( x = c(0,1), y = c(0,1), add = TRUE)
    legend( x = "bottom", bty = "n",
            # inset = c(0,-1), xpd = TRUE,
            legend = fleetLabs, 
            pch = 16, col = fleetCols,
            horiz = TRUE, cex = 1, pt.cex = 1.5)
} # END plotDataSummary()



# plotBatchLossDists_Scenario()
# Multi-panel plot of relative and absolute
# loss under a batch of MPs for all stock/species
# combinations, including data-pooled and 
# coast-wide aggregations.
plotBatchLossDists_Scenario <- function(  groupFolder = "DERTACS_reruns_sep24",
                                          prefix = "parBat",
                                          lossType = "abs",
                                          var = "C_ispt",
                                          period = 73:82,
                                          lossList = NULL,
                                          dim1 = 1:3,   # species (D,E,R, DP)
                                          dim2 = 1:3,   # stocks (H,Q,W, CW)
                                          qProbs = c(.025,.5,.975),
                                          refPts = "MSrefPts",
                                          AMlabs = rev(c( SS = "singleStock",
                                                          HMS = "hierMultiStock",
                                                          SpecPool = "speciesPooling",
                                                          SpatPool = "spatialPooling",
                                                          TA = "totalAgg" )),
                                          scenLabs = c( Rich  = "DERfit_HcMcAsSsIdx",
                                                        Mod  = "DERfit_McAsSsIdx",
                                                        Poor = "DERfit_AsSsIdx" ),
                                          clearBadReps = TRUE,
                                          minSampSize = 100,
                                          rotxlabs = 0,
                                          xlab = TRUE,
                                          vertLines = TRUE  )
{
  # First, read info files from the relevant
  # sims
  
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  simFolderList <- simFolderList[grepl(prefix, simFolderList)]

  infoFiles     <- lapply(  X = here::here("Outputs",groupFolder,simFolderList,"infoFile.txt"),
                            FUN = lisread )


  if( lossType == "rel" )
  {
    lossArrayName <- "totRelLoss_isp"
    yLab <- "Total Relative Loss"
  }

  if( lossType == "abs" )
  {
    lossArrayName <- "totAbsLoss_isp"
    yLab <- "Cumulative Absolute Catch Loss (kt)"
  }

  # Break up MP names into factor levels
  # MP labels are AM_Fsrce_eqbm
  splitMP <- function(  mpLab, 
                        breaks = "_",
                        factorNames = c("AM","Fsrce","eqbm") )
  {
    splitMP <- stringr::str_split(mpLab, breaks)[[1]]

    outList <- vector(  mode = "list", 
                        length = length(splitMP) )
    names(outList) <- factorNames
    for( k in 1:length(splitMP))
      outList[[k]] <- splitMP[k]

    outList
  }



  for( lIdx in 1:length(infoFiles) )
  {
    infoFiles[[lIdx]] <- c(infoFiles[[lIdx]], splitMP(infoFiles[[lIdx]]$mp))
    infoFiles[[lIdx]] <- as.data.frame(infoFiles[[lIdx]])
  }



  info.df <- do.call(rbind, infoFiles) %>% 
              mutate_if(is.factor, as.character) %>%
              filter( scenario %in% scenLabs )


  if( !is.null(refPts) )
    info.df <- info.df %>% filter( eqbm == refPts )


  # Load loss files, place in a list
  simFolderNames  <- info.df$simLabel
  if( is.null(lossList) )
  {
    lossList        <- lapply(  X = simFolderNames,
                                FUN = .loadLoss,
                                folder = groupFolder )
    names(lossList) <- simFolderNames
  }


  # Calculate total loss for variable/period
  totLossList  <- lapply( X = lossList,
                          FUN = calcTotalLossPeriod,
                          var = var, period = period )
  names(totLossList) <- names(lossList)


  nS  <- lossList[[1]]$nS
  nP  <- lossList[[1]]$nP
  nT  <- lossList[[1]]$nT
  tMP <- lossList[[1]]$tMP
  pT  <- lossList[[1]]$pT

  # Get number of MPs
  nSims <- length(lossList)


  AMs       <- unique(info.df$AM)
  Fsources  <- unique(info.df$Fsrce)
  eqbm      <- unique(info.df$eqbm)
  MPs       <- unique(info.df$mp)
  scenarios <- unique(info.df$scenario)

  nAM       <- length(AMs)
  nSrce     <- length(Fsources)
  nEqbm     <- length(eqbm)
  nMP       <- length(MPs)
  nScen     <- length(scenarios)



  # Some axis labels and translations from
  # simple axis label to info file AM name
  AMcodes <- c( SS = 1,
                HMS = 2,
                SpecPool = 3,
                SpatPool = 4,
                TotPool = 5 )

  AMdictionary <- AMlabs

  AMcols          <- RColorBrewer::brewer.pal(length(AMcodes), "Dark2")
  names(AMcols)   <- names(AMcodes)
  AMpch           <- 20 + 1:length(AMcodes)
  names(AMpch)    <- names(AMcodes)



  MPgrid <- list( AMs = names(AMdictionary),
                  Fsrce = Fsources,
                  eqbm = eqbm )

  MPgrid <- expand.grid(MPgrid)


  xJitter <- seq(from = -.4, to = .4, length.out = nAM )
  if(nAM == 1 )
    xJitter <- 0
  nReps_k <- numeric(length = nSims)
  goodRepsList <- vector(mode = "list", length = nSims)

  for( k in 1:nSims )
  {
    simID <- info.df$simLabel[k]
    simLoss <- totLossList[[simID]][[lossArrayName]] 
    nReps_k[k] <- dim(simLoss)[1]
    goodRepsList[[k]]$goodReps_isp <- lossList[[k]]$goodReps_isp
    goodRepsList[[k]]$nReps_sp     <- apply(X = lossList[[k]]$goodReps_isp, FUN = sum, MARGIN = c(2,3) )
    names(goodRepsList)[k] <- simID
  }

  nReps         <- max(nReps_k)
  speciesNames  <- lossList[[1]]$speciesNames
  stockNames    <- lossList[[1]]$stockNames

  speciesNames[nS+1]  <- "SpeciesPooled"
  stockNames[nP+1]    <- "SpatialPooled"

  # Make an array to hold loss function values
  totLossArray_misp <- array(NA,  dim = c(nSims, nReps, nS+1, nP+1 ),
                                  dimnames = list(  simID = info.df$simLabel,
                                                    rep = 1:nReps,
                                                    species = speciesNames,
                                                    stock = stockNames ) )

  for( k in 1:nSims )
  {
    simID <- info.df$simLabel[k]
    simLoss <- totLossList[[simID]][[lossArrayName]] 

    totLossArray_misp[simID,1:dim(simLoss)[1],,] <- simLoss[1:dim(simLoss)[1],,]

    if(clearBadReps)
    {
      goodReps_isp <- goodRepsList[[k]]$goodReps_isp[1:nReps_k[k],,]
      totLossArray_misp[simID,,,][!goodReps_isp] <- NA

      # # remove loss when nReps is below min sample size
      # for(s in 1:nS)
      #   for( p in 1:nP )
      #     if( goodRepsList[[k]]$nReps_sp[s,p] < minSampSize )
      #       totLossArray_misp[simID,,s,p] <- NA            
    }
  }

  totLossQuantiles_qmsp <- apply( X = totLossArray_misp,
                                  FUN = quantile,
                                  probs = qProbs,
                                  MARGIN = c(1,3,4),
                                  na.rm = TRUE )

  # Start plotting windows
  nSpec   <- length(dim1)
  nStock  <- length(dim2) 

  par(  mfcol = c(nStock,nSpec),
        mar = c(0.1,1.5,.1,1.5),
        oma = c(6,4,3,3) )

  for( sIdx in 1:nSpec )
    for( pIdx in 1:nStock )
    {
      s <- dim1[sIdx]
      p <- dim2[pIdx]

      lossRange <- max(abs(range(totLossQuantiles_qmsp[,,s,p],na.rm = TRUE)))

      if( lossRange == Inf)
        lossRange <- 1

      plot( y = c(0,lossRange), x = c(0.5,nScen + .5),
            type = "n", axes = FALSE, xaxs = "i"  )
        
        mfg <- par("mfg")
        if( mfg[1] == mfg[3])
        {
          if( rotxlabs == 0 )
            axis( side = 1, at = 1:nScen, labels = names(scenLabs) )
          if( rotxlabs > 0)
          {
            axis( side = 1, at = 1:nScen, labels = NA )
            xaxisRot( at = 1:nScen, 
                      labs = names(scenLabs),
                      font = 1,
                      cex = 1,
                      line = .05,
                      rot = rotxlabs )
            

          }

        }


        axis( side = 2, las = 1 )

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 1 )

        if( mfg[2] == mfg[4] )
        {
          rmtext( txt = stockNames[p], font = 2, line = .05,
                  outer = TRUE, cex = 1.5 )
        }

        grid()
        box()
        
        jitIdx <- 0
        for( amIdx in 1:nAM )
          for( eqIdx in 1:nEqbm )
          {
            jitIdx  <- jitIdx + 1
            AMname  <- AMlabs[amIdx]
            AMid    <- AMdictionary[amIdx]
            eqbmID  <- eqbm[eqIdx]
            
            subInfo <- info.df %>%
                        filter( AM    == AMid,
                                eqbm  == eqbmID ) 
            
            simLabels <- subInfo$simLabel
            
            FsrceID <- subInfo$Fsrce[1]
            AMcode  <- AMdictionary[amIdx]

            MPlab <- paste( AMid, FsrceID, eqbmID, sep = "_")

            for( scenIdx in 1:nScen )
            {

              thisSimLabel <- subInfo[subInfo$scenario == scenLabs[scenIdx],]$simLabel


              points( x = scenIdx + xJitter[jitIdx],
                      y = totLossQuantiles_qmsp[2,thisSimLabel,s,p],
                      pch = AMpch[names(AMid)],
                      bg = AMcols[names(AMid)],
                      col = AMcols[names(AMid)], cex = 1.5 )
              # lines(  x = 1:3 + xJitter[jitIdx],
              #         y = totLossQuantiles_qmsp[2,thisSimLabel,s,p],
              #         lty = 1, lwd = .8, col = AMcols[AMid] )

              if(goodRepsList[[thisSimLabel]]$nReps_sp[s,p] < minSampSize)
                lineType <- 2
              else lineType <- 1

              segments( y0 = totLossQuantiles_qmsp[1,thisSimLabel,s,p],
                        y1 = totLossQuantiles_qmsp[3,thisSimLabel,s,p],
                        x0 = scenIdx + xJitter[jitIdx],
                        lty = lineType,
                        col = AMcols[names(AMid)], lwd = 2  )
            }

          }
        if( vertLines)
          abline( v = 1:(nScen-1) + .5,
                  col = "black", lty = 1, lwd = 1.2 )


      }

  mtext( side = 2, outer = TRUE, text = yLab, line = 2 )
  
  if(xlab)
    mtext( side = 1, outer = TRUE, text = "Data Scenario", line = 2 )

  par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend( x       = "bottom",
          horiz   = TRUE,
          bty     = "n",
          legend  = c(names(AMdictionary)),
          pch     = AMpch[names(AMdictionary)],
          pt.bg   = AMcols[names(AMdictionary)], 
          lty     = c(rep(1,5),NA),
          col     = AMcols[names(AMdictionary)], 
          lwd     = 2 )
  
}

plotSensRuns_pub <- function()
{
  scenLabs <- list( MSYCVhigh  = c( High_MSYCV0.1  = "High_MSYCV0.1",
                                    High_MSYCV0.5  = "High_MSYCV0.5",
                                    High_MSYCV1.0  = "High_MSYCV1.0"),
                    MSYCVlow  = c( Low_MSYCV0.1  = "Low_MSYCV0.1",
                                    Low_MSYCV0.5  = "Low_MSYCV0.5",
                                    Low_MSYCV1.0  = "Low_MSYCV1.0" ),
                    UmsyCVhigh = c( High_UmsyCV0.1  = "Rich_UmsyCV0.1",
                                    High_UmsyCV0.5  = "Rich_UmsyCV0.5",
                                    High_UmsyCV1.0  = "Rich_UmsyCV1.0"),
                    UmsyCVlow = c( Low_UmsyCV0.1  = "Poor_UmsyCV0.1",
                                    Low_UmsyCV0.5  = "Poor_UmsyCV0.5",
                                    Low_UmsyCV1.0  = "Poor_UmsyCV1.0" ),
                    hierSDhigh = c( High_hierSD0.1  = "High_shrinkSD0.1",
                                    High_hierSD0.2  = "High_shrinkSD0.2",
                                    High_hierSD0.5  = "High_shrinkSD0.5"),
                    hierSDlow = c( Low_hierSD0.1  = "Low_shrinkSD0.1",
                                    Low_hierSD0.2  = "Low_shrinkSD0.2",
                                    Low_hierSD0.5  = "Low_shrinkSD0.5" ),
                    obsErrhigh = c( High_obsErr0.1  = "High_obsErr0.1",
                                    High_obsErr0.5  = "High_obsErr0.5",
                                    High_obsErr1.0  = "High_obsErr1.0"),
                    obsErrlow = c( Low_obsErr0.1  = "Low_obsErr0.1",
                                    Low_obsErr0.5  = "Low_obsErr0.5",
                                    Low_obsErr1.0  = "Low_obsErr1.0" )
                     )

  scenSplitString <- c( MSY     = "MSYCV",
                        Umsy    = "UmsyCV",
                        hierSD  = "shrinkSD",
                        obsErr  = "obsErr" )

  xlab <- c(  Bmsy = expression(CV(MSY)), 
              Umsy = expression( SD(U[MSY])), 
              hierSD = expression(paste( tau[q], " and ", sigma[U[MSY]] ) ),
              obsErr = expression(tau) )

  groupFolders <- c("sensRuns_MSYCV_Jan5", "sensRuns_UmsyCV_Jan5","sensRuns_hierSD_Jan5", "sensRuns_obsErr_Jan5")


  par( mfcol = c(2,4), oma = c(6,3,2,1), mar = c(.1,2,.1,1) )
  for( sensIdx in 1:4 )
  {

    plotSensBatchSummary_SAisp( groupFolder = groupFolders[sensIdx],
                                prefix = "parBat",
                                lossType = "abs",
                                var = "C_ispt",
                                period = 73:82,
                                lossList = NULL,
                                refPts = NULL,
                                AMlabs = rev(c( SS = "singleStock",
                                                HMS = "hierMultiStock",
                                                SpePool = "speciesPooling",
                                                SpaPool = "spatialPooling",
                                                TotPool = "totalAgg" )),
                                scenLabs = scenLabs[[2 * sensIdx - 1]],
                                clearBadReps = FALSE,
                                scenSplitString = scenSplitString[sensIdx],
                                xlab = xlab[sensIdx],
                                noPar = TRUE,
                                printLeg = TRUE
                             )

    if( sensIdx == 4)
      rmtext( txt = "High", font = 2, line = .05, outer = TRUE, cex = 1.5 )

    plotSensBatchSummary_SAisp( groupFolder = groupFolders[sensIdx],
                                prefix = "parBat",
                                lossType = "abs",
                                var = "SB_ispt",
                                period = 73:82,
                                lossList = NULL,
                                refPts = NULL,
                                AMlabs = rev(c( SS = "singleStock",
                                                HMS = "hierMultiStock",
                                                SpePool = "speciesPooling",
                                                SpaPool = "spatialPooling",
                                                TotPool = "totalAgg" ) ),
                                scenLabs = scenLabs[[2 * sensIdx]],
                                clearBadReps = FALSE,
                                scenSplitString = scenSplitString[sensIdx],
                                xlab = xlab[sensIdx],
                                noPar = TRUE,
                                printLeg = TRUE
                             )

    if( sensIdx == 4)
      rmtext( txt = "Low", font = 2, line = .05, outer = TRUE, cex = 1.5 )
    

    # mtext( side = 1, text = xlab[sensIdx], font = 2, line = 3)
  }

  mtext( side = 2, text = "Standardised difference from mean absolute cumulative loss", 
          outer = TRUE, line = 1)
  # Plot legend
  par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend( x       = "bottom",
            horiz   = TRUE,
            bty     = "n",
            legend  = rev(c("Single-species","Hierarchical","Species Pooling","Spatial Pooling","Total Pooling")),
            lty     = 1:5,
            pch     = 21:25,
            bg      = NA,
            col     = RColorBrewer::brewer.pal(5,"Dark2"), 
            lwd     = 2 )
}


# plotBatchLossDists()
# Multi-panel plot of relative and absolute
# loss under a batch of MPs for all stock/species
# combinations, including data-pooled and 
# coast-wide aggregations.
plotBatchLossDists <- function( groupFolder = "diffCV_fixedF_longGrid",
                                prefix = "sim_parBatfixedF",
                                lossType = "rel",
                                var = "C_ispt",
                                period = 62:72,
                                lossList = NULL )
{
  # First, read info files from the relevant
  # sims
  
  simFolderList <- list.dirs(here::here("Outputs",groupFolder))
  simFolderList <- simFolderList[grepl(prefix, simFolderList)]

  infoFiles     <- lapply(  X = file.path(simFolderList,"infoFile.txt"),
                            FUN = lisread )

  if( lossType == "rel" )
  {
    lossArrayName <- "totRelLoss_isp"
    yLab <- "Total Relative Loss"
  }

  if( lossType == "abs" )
  {
    lossArrayName <- "totAbsLoss_isp"
    yLab <- "Total Absolute Loss (kt)"
  }

  # Break up MP names into factor levels
  # MP labels are AM_Fsrce_eqbm
  splitMP <- function(  mpLab, 
                        breaks = "_",
                        factorNames = c("AM","Fsrce","eqbm") )
  {
    splitMP <- stringr::str_split(mpLab, breaks)[[1]]

    outList <- vector(  mode = "list", 
                        length = length(splitMP) )
    names(outList) <- factorNames
    for( k in 1:length(splitMP))
      outList[[k]] <- splitMP[k]

    outList
  }



  for( lIdx in 1:length(infoFiles) )
  {
    infoFiles[[lIdx]] <- c(infoFiles[[lIdx]], splitMP(infoFiles[[lIdx]]$mp))
    infoFiles[[lIdx]] <- as.data.frame(infoFiles[[lIdx]])
  }

  info.df <- do.call(rbind, infoFiles) %>% 
              mutate_if(is.factor, as.character)


  # Load loss files, place in a list
  simFolderNames  <- info.df$simLabel
  if( is.null(lossList) )
  {
    lossList        <- lapply(  X = simFolderNames,
                                FUN = .loadLoss,
                                folder = groupFolder )
    names(lossList) <- simFolderNames
  }
  

  # Calculate total loss for variable/period
  totLossList  <- lapply( X = lossList,
                          FUN = calcTotalLossPeriod,
                          var = var, period = period )
  names(totLossList) <- names(lossList)


  nS  <- lossList[[1]]$nS
  nP  <- lossList[[1]]$nP
  nT  <- lossList[[1]]$nT
  tMP <- lossList[[1]]$tMP
  pT  <- lossList[[1]]$pT

  # Get number of MPs
  nSims <- length(lossList)


  AMs       <- unique(info.df$AM)
  Fsources  <- unique(info.df$Fsrce)
  eqbm      <- unique(info.df$eqbm)

  nAM       <- length(AMs)
  nSrce     <- length(Fsources)
  nEqbm     <- length(eqbm)

  # Some axis labels and translations from
  # simple axis label to info file AM name
  AMcodes <- c( SS = 1,
                MS = 2,
                DP = 3,
                CW = 4,
                TA = 5 )

  AMdictionary <- c(  singleStock     = "SS",
                      hierMultiStock  = "MS",
                      dataPooled      = "DP",
                      coastwide       = "CW",
                      totalAgg        = "TA" )

  AMcols          <- RColorBrewer::brewer.pal(nAM, "Dark2")
  names(AMcols)   <- names(AMdictionary)
  FsrcePCH  <- 15 + 1:nSrce
  names(FsrcePCH) <- Fsources
  eqbmLTY   <- 1:nEqbm
  names(eqbmLTY)  <- eqbm

  nReps         <- dim(lossList[[1]]$lossRel[[var]])[1]
  speciesNames  <- lossList[[1]]$speciesNames
  stockNames    <- lossList[[1]]$stockNames

  # Make an array to hold loss function values
  totLossArray_misp <- array(NA,  dim = c(nSims, nReps, nS+1, nP+1 ),
                                  dimnames = list(  simID = info.df$simLabel,
                                                    rep = 1:nReps,
                                                    species = speciesNames,
                                                    stock = stockNames ) )

  for( k in 1:nSims )
  {
    simID <- info.df$simLabel[k]
    totLossArray_misp[k,,,] <- totLossList[[simID]][[lossArrayName]] 
  }

  totLossQuantiles_qmsp <- apply( X = totLossArray_misp,
                                  FUN = quantile,
                                  probs = c(0.25, 0.5, 0.75),
                                  MARGIN = c(1,3,4) )

  

  # Loop and plot - start by plotting ALL
  # combos
  yJitter <- seq(from = -.3, to = .3, length = nSrce * nEqbm )
  names(yJitter) <- apply(expand.grid(Fsources, eqbm), 1, paste, collapse="_")

  # Plot total loss for given period
  # on grid of Species/stocks
  par(  mfcol = c(nP+1,nS+1), 
        mar = c(.2,1.5,.2,1),
        oma = c(3,3.5,3,3) )

  for( s in 1:(nS+1) )
    for( p in 1:(nP+1) )
    {
      lossRange <- max(abs(range(totLossQuantiles_qmsp[,,s,p])))

      plot( y = c(0,lossRange), x = c(0,6),
            type = "n", axes = FALSE  )
        
        mfg <- par("mfg")
        if( mfg[1] == mfg[3])
          axis( side = 1, at = AMcodes, labels = names(AMcodes) )

        axis( side = 2, las = 1 )

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 1 )

        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = stockNames[p], font = 2, line = 1 )

        grid()
        box()

        for( k in 1:nSims )
        {
          simID   <- info.df$simLabel[k]
          AMid    <- info.df$AM[k]
          FsrceID <- info.df$Fsrce[k]
          eqbmID  <- info.df$eqbm[k]
          AMcode  <- AMdictionary[AMid]

          jitY    <- yJitter[paste(FsrceID,eqbmID,sep = "_")]

          points( y   = totLossQuantiles_qmsp[2,simID,s,p],
                  x   = AMcodes[AMcode] + jitY,
                  pch = FsrcePCH[FsrceID],
                  col = AMcols[AMid], cex = 1.5 )
          segments( y0 = totLossQuantiles_qmsp[1,simID,s,p],
                    y1 = totLossQuantiles_qmsp[3,simID,s,p],
                    x0 = AMcodes[AMcode] + jitY,
                    lty = eqbmLTY[eqbmID],
                    col = AMcols[AMid], lwd = 2  )
        }
    } # END p loop
    # END s loop

  legend( x       = "bottomright",
          bty     = "n",
          legend  = c(Fsources,eqbm),
          pch     = c(FsrcePCH,NA,NA),
          lty     = c(NA,NA,eqbmLTY),
          lwd     = 2 )


  mtext( side = 2, outer = TRUE, text = yLab, line = 2 )
  mtext( side = 1, outer = TRUE, text = "AM Configuration", line = 2 )


} # END plotBatchLossDists()

# plotTotLossDists()
# Plots relative and absolute loss
# for a given simulation. Requires
# loss for a given baseline to have
# been calculated first, and saved
# into the sim folder
plotTotLossDists <- function( sim = 1, 
                              groupFolder = "shortGrid",
                              lossType    = "rel",
                              var         = "SB_ispt",
                              save        = FALSE )
{
  # Load loss reports
  objLoss <- .loadLoss(sim = sim, groupFolder)
  simFolder <- objLoss$simFolder

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
  
  loss_shortTerm  <- calcTotalLossPeriod(objLoss,var, period = tMP:(tMP+10))
  loss_medTerm    <- calcTotalLossPeriod(objLoss,var, period = tMP:(tMP+20))
  loss_longTerm   <- calcTotalLossPeriod(objLoss,var, period = tMP:nT)

  if( lossType == "rel" )
  {
    lossListName  <- "totRelLoss_isp"  
    yLab          <- "Total relative loss (unitless)"
  }

  if( lossType == "abs" )
  {
    lossListName  <- "totAbsLoss_isp"    
    yLab          <- "Total loss (kt)"
  }

  if( save )
    pdf( file = file.path(simFolder,paste(lossType,var,"LossCleveland.pdf",sep = "_")),
          width = 11, height = 8.5 )

  # Plot window - just for loss right now
  par(  mfcol = c(nP + 1, nS + 1),
        mar = c(.5,.5,.5,.5),
        oma = c(3,3,3,3) )

  for( s in 1:(nS+1) )
    for( p in 1:(nP + 1) )
    {
      lossShort_i <- loss_shortTerm[[lossListName]][,s,p]
      lossMed_i   <- loss_medTerm[[lossListName]][,s,p]
      lossLong_i  <- loss_longTerm[[lossListName]][,s,p]

      lossShort_q <- quantile(lossShort_i, probs = c(0.025, 0.5, 0.975) )
      lossMed_q   <- quantile(lossMed_i, probs = c(0.025, 0.5, 0.975) )
      lossLong_q  <- quantile(lossLong_i, probs = c(0.025, 0.5, 0.975) )

      maxLoss     <- max( abs(lossShort_q[is.finite(lossShort_q)]),
                          abs(lossMed_q[is.finite(lossMed_q)]),
                          abs(lossLong_q[is.finite(lossLong_q)]), na.rm = T )
      
      plot( x = c(0,4), y = c(0,maxLoss), type = "n",
            axes = FALSE )
        mfg <- par("mfg")
        if( mfg[1] == mfg[3] )
          axis( side = 1, at = 1:3, labels = c("S", "M", "L") )
        if( mfg[2] == 1 )
          axis( side = 2, las = 1)
        if( mfg[1] == 1 )
        {
          mtext( side = 3, text = speciesNames[s] )
        }
        if( mfg[2] == mfg[4] )
        {
          mtext( side = 4, text = stockNames[p] )
        }
        grid()
        box()

        segments( x0 = 1, y0 = lossShort_q[1], y1 = lossShort_q[3]  )
        points( x = 1, y = lossShort_q[2]  )
        segments( x0 = 2, y0 = lossMed_q[1], y1 = lossMed_q[3]  )
        points( x = 2, y = lossMed_q[2]  )
        segments( x0 = 3, y0 = lossLong_q[1], y1 = lossLong_q[3]  )
        points( x = 3, y = lossLong_q[2]  )

    }

  mtext( side = 1, text = "Time Period", outer = T, line = 2)
  mtext( side = 2, text = yLab, outer = T, line = 2)

  if( save )
    dev.off()
} # END plotLossDists()

# plotLossTulips()
# Plots relative and absolute loss
# for a given simulation. Requires
# loss for a given baseline to have
# been calculated first, and saved
# into the sim folder
plotLossTulip <- function(  sim = 2, 
                            groupFolder   = "DERTACS_Reruns_Jan5",
                            lossType      = "abs",
                            var           = "C_ispt",
                            dim1          = 1:2,
                            dim2          = 1:3,
                            save          = FALSE,
                            clearBadReps  = FALSE )
{
  # Load loss reports
  objLoss <- .loadLoss(sim = sim, groupFolder)
  simFolderName <- basename(objLoss$simFolder)
  simFolder <- here::here("Outputs",groupFolder,simFolderName)

  goodReps_isp <- objLoss$goodReps_isp

  if(lossType == "rel" )
  {
    loss_ispt <- objLoss$lossRel[[var]]
    yLab      <- "Relative Loss"
  }

  if(lossType == "abs" )
  {
    loss_ispt <- objLoss$lossRaw[[var]]
    yLab      <- "Absolute loss (kt)"
  }


  # Pull model dimensions
  tMP <- objLoss$tMP
  nT  <- objLoss$nT
  nF  <- objLoss$nF
  nS  <- objLoss$nS
  nP  <- objLoss$nP
  pT  <- objLoss$pT

  if( clearBadReps )
    for( s in 1:nS )
      for( p in 1:nP )
      {
        badIdx <- which(!goodReps_isp[,s,p])
        loss_ispt[badIdx,s,p,] <- NA
      }

  nReps <- dim(loss_ispt)[1]

  # Compute max number of reps
  traces <- c()
  if( nReps >= 3)
    traces <- sample(1:nReps, size = 3)


  # Time steps
  fYear   <- objLoss$fYear
  tdxPlot <- tMP:nT
  years   <- seq(from = fYear, by = 1, length.out = nT)

  if( save )
  { 
    graphics.off()
    pdf( file = file.path(simFolder,paste(lossType,var,"LossEnvelopes.pdf",sep = "_")),
          width = 11, height = 8.5 )
  }

  # Plot window - just for loss right now
  par(  mfcol = c(length(dim2), length(dim1)),
        mar = c(.5,2,.5,1),
        oma = c(3,3,3,3) )

  speciesNames  <- c(objLoss$speciesNames,"data-pooled")
  stockNames    <- c(objLoss$stockNames,"coast-wide")

  loss_qspt <- apply( X = loss_ispt, FUN = quantile,
                      probs = c(0.025, 0.5, 0.975),
                      MARGIN = c(2,3,4), na.rm = TRUE )

  plotYrs <- years[tdxPlot]


  # Loop and plot
  for( s in dim1 )
    for( p in dim2 )
    {
      maxLoss <- max(abs(loss_qspt[,s,p,]),na.rm = T)
      

      plot( x = range(years[tdxPlot]), y = c(-maxLoss,maxLoss),
            type = "n", axes = FALSE  )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      axis( side = 2, las = 1)
      if( mfg[1] == 1 )
      {
        mtext( side = 3, text = speciesNames[s], font = 2 )
      }
      if( mfg[2] == mfg[4] )
      {
        rmtext( outer =TRUE, line = 0.05, txt = stockNames[p], font = 2, cex = 1.5 )
      }
      grid()
      box()

      abline( v = years[tMP], lty = 2, lwd = 3 )

      polygon(  x = c(years,rev(years)),
                y = c(loss_qspt[1,s,p,],rev(loss_qspt[3,s,p,])),
                border = NA, col = "grey65" )
      lines(  x = years, y = loss_qspt[2,s,p,],
              col = "black", lwd = 3 )
      for( traceIdx in traces )
        lines( x = years, y = loss_ispt[traceIdx,s,p,],
                col = "black", lwd = 1 )
      grid()

      abline( h = 0, lty = 2, lwd = 1 )


    }
    mtext( side = 1, text = "Year", line = 2, outer = TRUE )
    mtext( side = 2, text = yLab, line = 1.5, outer = TRUE )

  if( save )
    dev.off()

} # END plotLoss()

plotSSvsMSrefPts <- function( obj = blob )
{
  # Pull reference points and curves
  rp            <- obj$rp[[1]]
  refCurves     <- rp$refCurves
  EmsyRefPts    <- rp$EmsyRefPts
  EmsyMSRefPts  <- rp$EmsyMSRefPts
  FmsyRefPts    <- rp$FmsyRefPts

  nP  <- obj$om$nP
  nS  <- obj$om$nS
  tMP <- obj$om$tMP
  nT  <- obj$om$nT

  # Pull qF
  qF_sp       <- blob$om$qF_ispft[1,,,2,nT]
  
  speciesNames <- blob$om$speciesNames
  stockNames   <- blob$om$stockNames

  # Want to compare Umsy and Bmsy between
  # multi-stock and single-stock models
  BmsySS_sp <- FmsyRefPts$BeqFmsy_sp
  BmsyMS_sp <- EmsyMSRefPts$BeqEmsy_sp

  UmsySS_sp <- FmsyRefPts$Umsy_sp
  UmsyMS_sp <- EmsyMSRefPts$Umsy_sp  

  B0_sp     <- obj$ctlList$opMod$histRpt$B0_sp

  speciesCols <- RColorBrewer::brewer.pal(n = nS, "Set2")

  cols_sp   <- matrix(  speciesCols, nrow = nS, ncol = nP,
                        byrow = FALSE )
  pch_sp    <- matrix(  21:23, nrow = nS, ncol = nP, 
                        byrow = TRUE )


  par(  mfrow = c(1,2),
        mar = c(1.5,4,0,0),
        oma = c(2,1,1,1) )

  # First plot Umsy
  plot( x = c(0,max(UmsySS_sp, UmsyMS_sp) ),
        y = c(0,max(UmsySS_sp, UmsyMS_sp) ),
        type = "n", 
        xlab = "", 
        ylab = "", las = 1 )
    mtext(side =1, text = "Single-species Umsy", line = 2 )
    mtext(side =2, text = "Multi-species Umsy", line = 3 )
    points( x = UmsySS_sp, y = UmsyMS_sp,
            cex = 3 * B0_sp / max(B0_sp),
            bg = cols_sp, pch = pch_sp )
    abline( a = 0, b = 1, lty = 2 )

  # First plot Umsy
  plot( x = c(0,max(BmsySS_sp, BmsyMS_sp) ),
        y = c(0,max(BmsySS_sp, BmsyMS_sp) ),
        type = "n", las = 1,
        xlab = "", 
        ylab = "" )
    mtext(side =1, text = "Single-species Bmsy", line = 2 )
    mtext(side =2, text = "Multi-species Bmsy", line = 2 )
    points( x = BmsySS_sp, y = BmsyMS_sp,
            cex = 3 * B0_sp / max(B0_sp),
            bg = cols_sp, pch = pch_sp )
    abline( a = 0, b = 1, lty = 2 )


} # END plotSSvsMSrefPts

plotInvDemandCurves <- function(  obj = blob,
                                  maxC = 3.5,
                                  maxP = 4,
                                  minQ = 0.05,
                                  sIdx = 1, 
                                  lambda_s = NULL )
{
  if(is.null(lambda_s))
    lambda_s  <- blob$ctlList$opMod$lambda_s

  price_s   <- blob$ctlList$opMod$price_s
  Yeq_sp    <- blob$rp[[1]]$FmsyRefPts$YeqFmsy_sp
  Yeq_s     <- apply( X = Yeq_sp, FUN = sum, MARGIN = 1)


  nT  <- obj$om$nT
  nP  <- obj$om$nP
  nS  <- obj$om$nS
  tMP <- obj$om$tMP
  specCols <- wes_palette("Rushmore1", 5, type = "discrete")[3:5]
  speciesNames <- blob$om$speciesNames

  if(is.null(maxC))
    maxC <- 3*max(Yeq_sp)

  C <- seq(from = minQ, to = maxC, length.out = 100)
  P_sc  <- array(0, dim = c(length(sIdx),100))

  for(s in sIdx)
    P_sc[s,] <- price_s[s] * ( C / Yeq_s[s] )^(-1/lambda_s[s])

  plot( x = c(0,maxC), y = c(0,maxP), 
        xlab = "Catch (kt)", ylab = "Unit price ($/kg)",
        las = 1, type = "n", yaxs = "i", xaxs = "i")
    
    for( s in sIdx )
      lines( x = C, y = P_sc[s,], lty = 1, lwd = 3, col = specCols[s] )

    segments( x0 = Yeq_s[sIdx], y0 = -1, y1 = price_s[sIdx], lty = 2, col = specCols[sIdx])
    segments( x0 = -1, x1 = Yeq_s[sIdx], y0 = price_s[sIdx], lty = 2, col = specCols[sIdx])


    legend( x = "topleft", bty = "n",
            col = c(specCols[sIdx]),
            lty = rep(1,length(sIdx)),
            lwd = rep(3,length(sIdx)),
            legend = paste0(speciesNames[sIdx],", lambda = ",round(lambda_s[sIdx],2) ))


} # END plotInvDemandCurves()

# plotUnitPriceTS()
# 
plotUnitPriceTS <- function( priceModel = priceFlexModel )
{
  refPrice_st <- priceModel$refPrice_st
  P_st        <- priceModel$P_st

  yrs <- 2006:2016

  specJitter <- c(-.2,0,.2)
  nS <- 3

  specCols <- wes_palette("Rushmore1", 5, type = "discrete")[3:5]
  speciesNames <- c("Dover","English","Rock")

  plot( x = range(yrs), y = c(0,max(P_st,refPrice_st)),
        type = "n", las = 1, xlab = "", ylab = "" )
    grid()
    mtext( side = 1, text = "Year", line = 2.5 )
    mtext( side = 2, text = "Unit price ($/kg)", line = 2.5 )
    for( s in 1:3 )
    {
      lines(  x = yrs + specJitter[s], y = refPrice_st[s,],
              col = "grey70", lwd = 1, lty = s  )
      points( x = yrs + specJitter[s], y = refPrice_st[s,],
              pch = 21, cex = .8, col = specCols[s],
              bg = "white" )
      segments( x0 = yrs + specJitter[s],
                y0 = refPrice_st[s,],
                y1 = P_st[s,],
                lty = 1, lwd = .8, col = "grey70" )
      points( x = yrs + specJitter[s], y = P_st[s,],
              pch = 16, cex = .8, col = specCols[s] )
    }
} # END plotUnitPriceTS()

# plotPriceDevTS()
# plots the time series of reference price RW jumps
plotPriceDevTS <- function( priceModel = priceFlexModel )
{
  eps_st <- priceModel$eps_st

  meanDev_s <- round(apply(X = eps_st, FUN = mean, MARGIN = 1),2)

  yrs <- 2006:2016

  specJitter <- c(-.2,0,.2)
  nS <- 3

  specCols <- wes_palette("Rushmore1", 5, type = "discrete")[3:5]
  speciesNames <- c("Dover","English","Rock")

  plot( x = range(yrs), y = 2.5*c(-1,1),
        type = "n", las = 1, xlab = "", ylab = "" )
    grid()
    abline(h = 0, lty = 2, lwd = 1.5 )
    for( s in 1:nS )
    {
      points( x = yrs + specJitter[s],
              y = eps_st[s,],
              col = specCols[s], pch = s, cex = 1.5 )
      # Add regression lines?
    }

    legend( x = "topleft", bty = "n",
            legend = paste0(speciesNames,", m = ", meanDev_s),
            pch = 1:3,
            col = specCols)

} # END plotPriceDevTS

# plotEconModelDash
# A dashboard of the economic model fit
plotEconModelDash <- function( obj = blob,
                                priceModel = priceFlexModel )
{
  par(mfrow = c(2,2), mar = c(3,3,2,1), oma = c(2,2,2,2) )

  # Inverse demand curves
  plotInvDemandCurves(obj = obj,sIdx = 1:3)
  mtext( side = 3, "2016 demand curves", font = 2)
  mtext( side = 1, text = "Catch (kt)", line = 2.5 )
  mtext( side = 2, text = "Unit price ($/kg)", line = 2.5)

  # Time series of reference prices
  plotUnitPriceTS(priceModel)
  mtext(side = 3, text = "Price time series", font = 2)

  # Deviations
  plotPriceDevTS(priceModel)
    mtext( side = 1, text= "Year", line = 2.5 )
    mtext( side = 2, text= expression(delta[s][t]), line = 2.5 )
    mtext( side = 3, text = "RW deviations", font = 2)

}

plotFYieldCurves <- function( obj = blob,
                              maxF = NULL,
                              pIdx = 1, sIdx = 1 )
{
  # First, pull reference points and curves
  rp            <- obj$rp[[1]]
  refCurves     <- rp$refCurves
  FmsyRefPts    <- rp$FmsyRefPts

  nT  <- obj$om$nT
  nP  <- obj$om$nP
  nS  <- obj$om$nS
  tMP <- obj$om$tMP

  speciesNames <- blob$om$speciesNames
  stockNames   <- blob$om$stockNames

  specCols <- wes_palette("Rushmore1", 5, type = "discrete")[3:5]

  # Now compute MSY for SS and MS curves
  # MSY is still the same for SS, just need the effort
  # that corresponds to it, which is Fmsy/qF

  # get SS ref points
  Fmsy_sp <- rp$FmsyRefPts$Fmsy_sp
  MSY_sp  <- rp$FmsyRefPts$YeqFmsy_sp
  Yeq_spf <- refCurves$Yeq_spf
  F       <- refCurves$F

  if(is.null(maxF))
    maxF <- F[max(which(Yeq_spf[sIdx,pIdx,] > 0) )]

  maxY <- max(Yeq_spf[sIdx,pIdx,])

  Yeq_spf[Yeq_spf < 0] <- NA

  plot( x = c(0,maxF), y = c(0,maxY), type = "n",
        xlab = "Fishing Mortality (/yr)",
        ylab = "Equilibrium Yield (kt)", las = 1)
    lines( x = F, y = Yeq_spf[sIdx,pIdx,], lwd = 3,
            col = specCols[sIdx] )
    mtext( side = 3, text = speciesNames[sIdx], font = 2 )


}

# plotEffYieldCurves()
# Function for plotting effort based
# yield curves for each species and
# the complex in a stock area.
plotEffYieldCurves <- function( obj = blob, 
                                maxE = NULL,
                                pIdx = 1:3,
                                sIdx = 1:3,
                                plotCplx = TRUE )
{
  # First, pull reference points and curves
  rp            <- obj$rp[[1]]
  refCurves     <- rp$refCurves
  EmsyRefPts    <- rp$EmsyRefPts
  EmsyMSRefPts  <- rp$EmsyMSRefPts
  FmsyRefPts    <- rp$FmsyRefPts

  nT  <- obj$om$nT
  nP  <- obj$om$nP
  nS  <- obj$om$nS
  tMP <- obj$om$tMP

  # Pull qF
  qF_sp       <- blob$om$qF_ispft[1,,,2,nT]
  
  speciesNames <- blob$om$speciesNames
  stockNames   <- blob$om$stockNames

  # Now compute MSY for SS and MS curves
  # MSY is still the same for SS, just need the effort
  # that corresponds to it, which is Fmsy/qF

  # get SS ref points
  EmsySS_sp <- rp$FmsyRefPts$Fmsy_sp / qF_sp
  MSYSS_sp  <- rp$FmsyRefPts$YeqFmsy_sp

  # Now get MS ref points
  Yeq_spe   <- rp$refCurves$EffCurves$Yeq_spe
  Yeq_spe[Yeq_spe < 0] <- NA
  Yeq_pe    <- apply(X = Yeq_spe, FUN = sum, MARGIN = c(2,3), na.rm = T )
  Yeq_pe[Yeq_pe == 0] <- NA
  Yeq_pe[,1] <- 0
  Eseq      <- rp$refCurves$EffCurves$E

  EmsyMS_p  <- EmsyMSRefPts$EmsyMS_p
  MSYMS_sp  <- EmsyMSRefPts$YeqEmsy_sp
  MSYMS_p   <- EmsyMSRefPts$YeqEmsy_p
  BmsyMS_sp <- EmsyMSRefPts$BeqEmsy_sp

  Yeq_e <- apply( X = Yeq_pe, FUN = sum, na.rm = T,
                  MARGIN = 2 )
  maxEval <- max(which(Yeq_e > 0))

  specCols <- wes_palette("Rushmore1", 5, type = "discrete")[3:5]

  if(is.null(maxE))
    maxE <- Eseq[maxEval]

  par( mfrow = c(length(pIdx),1), mar = c(.1,1,.1,1), oma = c(3,4,1,1) )
  for( p in pIdx )
  {
    if( plotCplx | length(sIdx) == nS )
      maxY <- 1.05 * max(Yeq_pe[p,],na.rm = TRUE )
    else
      maxY <- max(Yeq_spe[,p,], na.rm = TRUE)

    plot( x = c(0,maxE), y = c(0, maxY ),
          type = "n", xlab = "", ylab = "", axes = F, xaxs="i",
          yaxs = "i" )
      axis(side = 2, las = 1)
      box()
      grid()

      for( s in sIdx )
      {
        lines( x = Eseq, y = Yeq_spe[s,p,],
               col = specCols[s], lty = 1, lwd = 2 )
        # lines( x = Eseq, y = Beq_spe[s,p,],
        #        col = specCols[s], lty = 2, lwd = 2 )

      }
      mfg <- par("mfg")
      if( mfg[1] == mfg[3])
        axis(side = 1)

      rmtext( txt = stockNames[p], line = 0.02,
              font = 2, cex = 1.5, outer = TRUE )

      if(plotCplx)
      {
        lines(  x = Eseq, y = Yeq_pe[p,],
              col = "black", lty = 1, lwd = 3 )

        segments( x0 = EmsyMS_p[p], col = "grey40",
                  y0 = 0, y1 = MSYMS_p[p], lty = 2  )
 
        segments( x0 = 0, x1 = EmsyMS_p[p], col = "grey40",
                  y0 = c(MSYMS_p[p],MSYMS_sp[,p]), lty = 2  )
      }

      if(  p == 1)
        legend( x = "topright", bty = "n",
                col = c(specCols,"black"),
                lwd = c(2,2,2,3),
                legend = c(speciesNames,"Complex") )

  }

  mtext( outer = TRUE, side = 1, text = "Commercial Trawl Effort (1000 hrs)", line = 2 )
  mtext( outer = TRUE, side = 2, text = "Equilibrium Yield (kt)", line = 2 )


  outList <- list(  EmsyMS_p = EmsyMS_p,
                    EmsySS_sp = EmsySS_sp,
                    MSYMS_p = MSYMS_p,
                    MSYSS_sp = MSYSS_sp,
                    MSYMS_sp = MSYMS_sp )
} # END plotEffYieldCurves()

# plotEffYieldCurves()
# Function for plotting effort based
# yield curves for each species and
# the complex in a stock area.
plotEconYieldCurves <- function(  obj = blob, 
                                  maxE = 60,
                                  pIdx = 1:3,
                                  sIdx = 1:3,
                                  plotCplx = TRUE )
{
  # First, pull reference points and curves
  rp            <- obj$rp[[1]]
  
  refCurves     <- rp$refCurves
  EmsyRefPts    <- rp$EmsyRefPts
  EmsyMSRefPts  <- rp$EmsyMSRefPts
  FmsyRefPts    <- rp$FmsyRefPts
  EmeyRefPts    <- rp$EmeyRefPts

  nT  <- obj$om$nT
  nP  <- obj$om$nP
  nS  <- obj$om$nS
  tMP <- obj$om$tMP

  crewShare <- obj$ctlList$opMod$crewShare

  # Pull qF
  qF_sp       <- blob$om$qF_ispft[1,,,2,nT]
  
  speciesNames <- blob$om$speciesNames
  stockNames   <- blob$om$stockNames

  # Now compute MSY for SS and MS curves
  # MSY is still the same for SS, just need the effort
  # that corresponds to it, which is Fmsy/qF
  Eseq      <- as.numeric(dimnames(EmeyRefPts$econRev_pe)[[2]])

  EmsyMS_p  <- EmsyMSRefPts$EmsyMS_p
  MSYMS_sp  <- EmsyMSRefPts$YeqEmsy_sp

  EmsySS_p  <- EmsyRefPts$Emsy_sp
  MSYSS_sp  <- EmsyRefPts$YeqEmsy_sp

  Emey_p      <- EmeyRefPts$Emey_p
  MEY_p       <- EmeyRefPts$MEY_p
  Rev_pe      <- EmeyRefPts$econRev_pe
  Rev_spe     <- EmeyRefPts$econRev_spe
  econYeq_pe  <- EmeyRefPts$econYeq_pe
  effCost_pe  <- EmeyRefPts$effCost_pe
  Ymey_sp     <- EmeyRefPts$Ymey_sp


  # remove zeroes
  Rev_spe[Rev_spe == 0] <- NA
  Rev_spe[,,1] <- 0

  specCols <- wes_palette("Rushmore1", 5, type = "discrete")[3:5]

  if(is.null(maxE))
    maxE <- 10 * max(Emey_p[pIdx])

  par( mfrow = c(length(pIdx),1), mar = c(.1,1,.1,1), oma = c(3,4,1,1) )
  for( p in pIdx )
  {
    plot( x = c(0,maxE), y = c(0, max(econYeq_pe[p,],Rev_pe[p,],na.rm = T) ),
          type = "n", xlab = "", ylab = "", axes = F, xaxs="i" )
      axis(side = 2, las = 1)
      box()
      grid()

      for( s in sIdx )
      {
        lines( x = Eseq, y = Rev_spe[s,p,],
               col = specCols[s], lty = 1, lwd = 2 )
        # lines( x = Eseq, y = Beq_spe[s,p,],
        #        col = specCols[s], lty = 2, lwd = 2 )

      }
      mfg <- par("mfg")
      if( mfg[1] == mfg[3])
        axis(side = 1)

      rmtext( txt = stockNames[p], line = 0.02,
              font = 2, cex = 1.5, outer = TRUE )

      if(plotCplx)
      {
        lines(  x = Eseq, y = Rev_pe[p,],
                col = "black", lty = 1, lwd = 3 )

        lines(  x = Eseq, y = econYeq_pe[p,],
                col = "black", lty = 2, lwd = 3 )
      

      lines( x = Eseq, y = effCost_pe[p,], col = "red",
              lwd = 2 )

      abline( v = Emey_p[p], lty = 2, col = "steelblue")
      abline( v = EmsyMS_p[p], lty = 2, col = "grey40")
      }



      # segments( x0 = EmsyMS_p[p], col = "grey40",
      #           y0 = 0, y1 = MSYMS_p[p], lty = 2  )

      # segments( x0 = 0, x1 = EmsyMS_p[p], col = "grey40",
      #           y0 = c(MSYMS_p[p],MSYMS_sp[,p]), lty = 2  )

      if(  p == 1 & length(sIdx) == nS )
        legend( x = "topright", bty = "n",
                col = c(specCols,"black","red","black"),
                lty = c(1,1,1,1,1,2),
                lwd = c(2,2,2,2,2,2),
                legend = c( paste(speciesNames," Revenue",sep = ""),
                            "Complex Revenue",
                            "Variable Costs",
                            "Total Profits") )

  }

  mtext( outer = TRUE, side = 1, text = "Commercial Trawl Effort (1000 hrs)", line = 2 )
  mtext( outer = TRUE, side = 2, text = "Cost, Revenue, and Profit ($ x 1e6)", line = 2 )


  # outList <- list(  EmsyMS_p = EmsyMS_p,
  #                   EmsySS_sp = EmsySS_sp,
  #                   MSYMS_p = MSYMS_p,
  #                   MSYSS_sp = MSYSS_sp,
  #                   MSYMS_sp = MSYMS_sp )
} # END plotEffYieldCurves()

# plotEmpYieldCurves
# Function to plot median simuated yield 
# resulting from a grid of constant fishing 
# mortality rates - used for checking the
# reference point calculations
plotEmpYieldCurves <- function( sims = 1:101, 
                                folder = "EmsyTesy",
                                indepVar = "E",
                                redoEmpRefCurves = FALSE )
{
  nSims <- length(sims)
  blobList <- vector(mode = "list", length = nSims)

  .loadSim(sims[1], folder = folder)

  nS <- blob$om$nS
  nP <- blob$om$nP
  nT <- blob$om$nT

  goodReps <- 1:dim(blob$goodReps_isp)[1]


  # Arrays to hold empirical eqbm yields
  C_spk <- array( 0, dim = c(nS,nP,nSims) )
  B_spk <- array( 0, dim = c(nS,nP,nSims) )
  F_spk <- array( 0, dim = c(nS,nP,nSims) )
  E_spk <- array( 0, dim = c(nS,nP,nSims) )
  E_pk  <- array( 0, dim = c(nP,nSims) )

  qF_sp <- blob$om$qF_ispft[1,,,2,nT]

  if(!file.exists(here::here("Outputs",folder,"empRefCurves.RData")) |
      redoEmpRefCurves )
  {
    for( x in sims )
    {
      .loadSim(x, folder = folder)
      C_spk[,,x]  <- apply(X = blob$om$C_ispt[goodReps,,,nT,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
      B_spk[,,x]  <- apply(X = blob$om$SB_ispt[goodReps,,,nT,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
      F_spk[,,x]  <- apply(X = blob$om$F_ispft[goodReps,,,2,nT,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
      E_pk[,x]    <- apply(X = blob$om$E_ipft[goodReps,,2,nT,drop = FALSE], FUN = median, MARGIN = c(2), na.rm = T )
      E_spk[,,x]  <- F_spk[,,x] / qF_sp

      # Clean up
      gc()
    }

    saveEmpRefCurves <- list( C_spk = C_spk,
                              B_spk = B_spk,
                              F_spk = F_spk,
                              E_spk = E_spk,
                              E_pk  = E_pk )

    save( saveEmpRefCurves, file = here::here("Outputs",folder,"empRefCurves.RData"))
  } else {
    # Load ref curve list
    load(here::here("Outputs",folder,"empRefCurves.RData"))

    C_spk <- saveEmpRefCurves$C_spk
    B_spk <- saveEmpRefCurves$B_spk
    F_spk <- saveEmpRefCurves$F_spk
    E_spk <- saveEmpRefCurves$E_spk
    E_pk  <- saveEmpRefCurves$E_pk 
  }



  

  # Pull F based ref curves from RP object
  refCurves       <- blob$rp[[1]]$refCurves
  Fvec            <- refCurves$F
  BeqRefCurve_spf <- refCurves$Beq_spf
  YeqRefCurve_spf <- refCurves$Yeq_spf

  # Pull effort based ref points
  EmsyRefPts      <- blob$rp[[1]]$EmsyRefPts
  Evec            <- refCurves$EffCurves$E

  YeqRefCurve_spe <- refCurves$EffCurves$Yeq_spe
  BeqRefCurve_spe <- refCurves$EffCurves$Beq_spe

  Xmsy_sp <- array(0, dim = c(nS,nP))
  Bmsy_sp <- array(0, dim = c(nS,nP))
  MSY_sp <- array(0, dim = c(nS,nP))

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(5,5,3,3) )

  for( s in 1:nS )
    for( p in 1:nP )
    {
      if( indepVar == "F" )
        maxX <- max(F_spk[s,p,])
      if( indepVar == "E" )
        maxX <- max(E_pk[p,],E_spk[s,p,])

      plot( x = c(0,maxX), y = c(0,max(B_spk[s,p,])),
            type = "n", axes = FALSE )
        axis( side = 1 )
        axis( side = 2, las = 1 )
        grid()
        box()

        F <- F_spk[s,p,]
        C <- C_spk[s,p,]
        B <- B_spk[s,p,]
        E <- E_spk[s,p,]

        actualOrder <- order(F)

        C <- C[actualOrder]
        B <- B[actualOrder]
        F <- F[actualOrder]
        E <- E[actualOrder]

        if( indepVar == "F" )
        {
          X     <- F
          Xvec  <- Fvec
          Yeq   <- YeqRefCurve_spf[s,p,]
          Beq   <- BeqRefCurve_spf[s,p,]

          xLab  <- "Fishing Mortality"
        }

        if( indepVar == "E" )
        {
          
          X     <- E
          Xvec  <- Evec
          Yeq   <- YeqRefCurve_spe[s,p,]
          Beq   <- BeqRefCurve_spe[s,p,]

          xLab  <- "Fishing Effort"
        }


        ubX <- max(which(Yeq >= 0) )

        CXspline <- splinefun(x = X, y = C)
        BXspline <- splinefun(x = X, y = B)

        empXmsy  <- try(uniroot(  interval = range(X),
                              f = CXspline,
                              deriv = 1 )$root)
        if( class(empXmsy) == "try-error")
          empXmsy <- 0
        
        empMSY  <- CXspline(empXmsy)
        empBmsy <- BXspline(empXmsy)

        empXmsy <- round(empXmsy,2)
        empBmsy <- round(empBmsy,2)
        empMSY  <- round(empMSY,2)


        

        lines( x = X, y = C,
                col = "steelblue", lwd = 2, lty = 1 )
        lines( x = X, y = B,
                col = "black", lwd = 2, lty = 1 )        

        lines( x = Xvec, y = Yeq,
                col = "salmon", lty = 2 )

        lines( x = Xvec, y = Beq,
                col = "black", lty = 2 )


        legend( "topright", bty = "n",
                legend = c( paste(indepVar,"msy = ", empXmsy, sep = "" ),
                            paste(" MSY = ", empMSY, sep = ""),
                            paste("Bmsy = ", empBmsy, sep = "") ) )


        # Plot some guidelines
        segments(  x0 = empXmsy, x1 = empXmsy,
                    y0 = 0, y1 = empBmsy,
                    lty = 2, col = "red" )
        segments(  x0 = 0, x1 = empXmsy,
                    y0 = empMSY, y1 = empMSY,
                    lty = 2, col = "red" )
        segments(  x0 = 0, x1 = empXmsy,
                    y0 = empBmsy, y1 = empBmsy,
                    lty = 2, col = "red" )
        
        Xmsy_sp[s,p] <- empXmsy
        Bmsy_sp[s,p] <- empBmsy
        MSY_sp[s,p]  <- empMSY
    }

  # if( indepVar == "F" )

  mtext( side = 1, text = xLab, outer = T, line = 2)
  mtext( side = 2, text = "Eqbm biomass and catch (kt)", outer = T, line = 2 )



  out <- list(  Xmsy_sp = Xmsy_sp,
                Bmsy_sp = Bmsy_sp,
                MSY_sp  = MSY_sp )

  out
} # END plotEmpYieldCurves()

# plotEmpYieldCurves
# Function to plot median simulated economic yield 
# resulting from a grid of constant fishing 
# effort rates - used for checking the
# reference point calculations
plotEmpEconYieldCurves <- function( sims = 1:101, 
                                    folder = "EmeyTesy",
                                    indepVar = "E",
                                    redoEmpRefCurves = FALSE )
{
  nSims <- length(sims)
  blobList <- vector(mode = "list", length = nSims)

  .loadSim(sims[1], folder = folder)

  nS <- blob$om$nS
  nP <- blob$om$nP
  nT <- blob$om$nT

  goodReps <- 1:dim(blob$goodReps_isp)[1]


  # Arrays to hold empirical eqbm yields
  C_spk <- array( 0, dim = c(nS,nP,nSims) )
  B_spk <- array( 0, dim = c(nS,nP,nSims) )
  F_spk <- array( 0, dim = c(nS,nP,nSims) )
  E_spk <- array( 0, dim = c(nS,nP,nSims) )
  E_pk  <- array( 0, dim = c(nP,nSims) )

  qF_sp <- blob$om$qF_ispft[1,,,2,nT]

  if(!file.exists(here::here("Outputs",folder,"empRefCurves.RData")) |
      redoEmpRefCurves )
  {
    for( x in sims )
    {
      .loadSim(x, folder = folder)
      C_spk[,,x]  <- apply(X = blob$om$C_ispt[goodReps,,,nT,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
      B_spk[,,x]  <- apply(X = blob$om$SB_ispt[goodReps,,,nT,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
      F_spk[,,x]  <- apply(X = blob$om$F_ispft[goodReps,,,2,nT,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
      E_pk[,x]    <- apply(X = blob$om$E_ipft[goodReps,,2,nT,drop = FALSE], FUN = median, MARGIN = c(2), na.rm = T )
      E_spk[,,x]  <- F_spk[,,x] / qF_sp

      # Clean up
      gc()
    }

    saveEmpRefCurves <- list( C_spk = C_spk,
                              B_spk = B_spk,
                              F_spk = F_spk,
                              E_spk = E_spk,
                              E_pk  = E_pk )

    save( saveEmpRefCurves, file = here::here("Outputs",folder,"empRefCurves.RData"))
  } else {
    # Load ref curve list
    load(here::here("Outputs",folder,"empRefCurves.RData"))

    C_spk <- saveEmpRefCurves$C_spk
    B_spk <- saveEmpRefCurves$B_spk
    F_spk <- saveEmpRefCurves$F_spk
    E_spk <- saveEmpRefCurves$E_spk
    E_pk  <- saveEmpRefCurves$E_pk 
  }



  

  # Pull F based ref curves from RP object
  refCurves       <- blob$rp[[1]]$refCurves
  Fvec            <- refCurves$F
  BeqRefCurve_spf <- refCurves$Beq_spf
  YeqRefCurve_spf <- refCurves$Yeq_spf

  # Pull effort based ref points
  EmsyRefPts      <- blob$rp[[1]]$EmsyRefPts
  Evec            <- refCurves$EffCurves$E

  YeqRefCurve_spe <- refCurves$EffCurves$Yeq_spe
  BeqRefCurve_spe <- refCurves$EffCurves$Beq_spe

  Xmsy_sp <- array(0, dim = c(nS,nP))
  Bmsy_sp <- array(0, dim = c(nS,nP))
  MSY_sp <- array(0, dim = c(nS,nP))

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(5,5,3,3) )

  for( s in 1:nS )
    for( p in 1:nP )
    {
      if( indepVar == "F" )
        maxX <- max(F_spk[s,p,])
      if( indepVar == "E" )
        maxX <- max(E_pk[p,],E_spk[s,p,])

      plot( x = c(0,maxX), y = c(0,max(B_spk[s,p,])),
            type = "n", axes = FALSE )
        axis( side = 1 )
        axis( side = 2, las = 1 )
        grid()
        box()

        F <- F_spk[s,p,]
        C <- C_spk[s,p,]
        B <- B_spk[s,p,]
        E <- E_spk[s,p,]

        actualOrder <- order(F)

        C <- C[actualOrder]
        B <- B[actualOrder]
        F <- F[actualOrder]
        E <- E[actualOrder]

        if( indepVar == "F" )
        {
          X     <- F
          Xvec  <- Fvec
          Yeq   <- YeqRefCurve_spf[s,p,]
          Beq   <- BeqRefCurve_spf[s,p,]

          xLab  <- "Fishing Mortality"
        }

        if( indepVar == "E" )
        {
          
          X     <- E
          Xvec  <- Evec
          Yeq   <- YeqRefCurve_spe[s,p,]
          Beq   <- BeqRefCurve_spe[s,p,]

          xLab  <- "Fishing Effort"
        }


        ubX <- max(which(Yeq >= 0) )

        CXspline <- splinefun(x = X, y = C)
        BXspline <- splinefun(x = X, y = B)

        empXmsy  <- try(uniroot(  interval = range(X),
                              f = CXspline,
                              deriv = 1 )$root)
        if( class(empXmsy) == "try-error")
          empXmsy <- 0
        
        empMSY  <- CXspline(empXmsy)
        empBmsy <- BXspline(empXmsy)

        empXmsy <- round(empXmsy,2)
        empBmsy <- round(empBmsy,2)
        empMSY  <- round(empMSY,2)


        

        lines( x = X, y = C,
                col = "steelblue", lwd = 2, lty = 1 )
        lines( x = X, y = B,
                col = "black", lwd = 2, lty = 1 )        

        lines( x = Xvec, y = Yeq,
                col = "salmon", lty = 2 )

        lines( x = Xvec, y = Beq,
                col = "black", lty = 2 )


        legend( "topright", bty = "n",
                legend = c( paste(indepVar,"msy = ", empXmsy, sep = "" ),
                            paste(" MSY = ", empMSY, sep = ""),
                            paste("Bmsy = ", empBmsy, sep = "") ) )


        # Plot some guidelines
        segments(  x0 = empXmsy, x1 = empXmsy,
                    y0 = 0, y1 = empBmsy,
                    lty = 2, col = "red" )
        segments(  x0 = 0, x1 = empXmsy,
                    y0 = empMSY, y1 = empMSY,
                    lty = 2, col = "red" )
        segments(  x0 = 0, x1 = empXmsy,
                    y0 = empBmsy, y1 = empBmsy,
                    lty = 2, col = "red" )
        
        Xmsy_sp[s,p] <- empXmsy
        Bmsy_sp[s,p] <- empBmsy
        MSY_sp[s,p]  <- empMSY
    }

  # if( indepVar == "F" )

  mtext( side = 1, text = xLab, outer = T, line = 2)
  mtext( side = 2, text = "Eqbm biomass and catch (kt)", outer = T, line = 2 )



  out <- list(  Xmsy_sp = Xmsy_sp,
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
  goodReps_isp <- obj$goodReps

  goodReps <- apply( X = goodReps_isp, FUN = prod, MARGIN = 1)

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
  goodReps_isp <- obj$goodReps

  goodReps <- as.logical(apply(X = goodReps_isp, FUN = prod, MARGIN = 1))

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
plotConvStats <- function( obj = blob, clearBadReps = FALSE )
{
  nSims         <- obj$nSims
  goodReps_isp  <- obj$goodReps_isp[1:nSims,,]

  # Pull max gradient value and hessian indicator
  maxGrad_itsp  <- obj$mp$assess$maxGrad_itsp[1:nSims,,,,drop = FALSE]
  pdHess_itsp   <- obj$mp$assess$pdHess_itsp[1:nSims,,,,drop = FALSE]

  nReps <- dim(maxGrad_itsp)[1]

  if( clearBadReps )
  {
    for( s in 1:nS )
      for( p in 1:nP )
      {
        badIdx <- which(!goodReps_isp[,s,p])
        maxGrad_itsp[badIdx,,s,p] <- NA
        pdHess_itsp[badIdx,,s,p] <- NA
      }
  }


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

plotTulipF <- function( obj = blob, nTrace = 3, clearBadReps = FALSE )
{
  # Model dimensions
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF

  # Good replicates
  nSims         <- obj$nSims
  goodReps_isp  <- obj$goodReps_isp[1:nSims,,]

  # Pull reference points
  Fmsy_sp <- obj$rp[[1]]$FmsyRefPts$Fmsy_sp
    
  # Fishing mortality series
  F_ispft <- obj$om$F_ispft[1:nSims,,,,]

  if(clearBadReps)
  {
    for( s in 1:nS )
      for(p in 1:nP )
      {
        badIdx <- which(!goodReps_isp[,s,p])
        F_ispft[badIdx,s,p,,] <- NA
      }
  }

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
plotTulipTACu <- function(  obj = blob, 
                            nTrace = 3,
                            clearBadReps = FALSE )
{

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT

  nSims         <- obj$nSims
  goodReps_isp  <- obj$goodReps_isp[1:nSims,,]
  

  # Get catch for trawl fleet, in projections only
  C_ispft     <- obj$om$C_ispft[1:nSims,,,2,tMP:nT,drop = FALSE]
  TAC_ispft   <- obj$mp$hcr$TAC_ispft[1:nSims,,,2,tMP:nT,drop = FALSE]

  nReps   <- dim(TAC_ispft)[1]

  TACu_ispft <- C_ispft / TAC_ispft

  if( clearBadReps )
  {
    for( s in 1:nS )
      for(p in 1:nP )
      {
        badIdx <- which(!goodReps_isp[,s,p])
        TACu_ispft[badIdx,s,p,,] <- NA
        C_ispft[badIdx,s,p,,] <- NA
        TAC_ispft[badIdx,s,p,,] <- NA
      }
  }



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
                          Ct  = FALSE,
                          leg = FALSE,
                          printStamp = TRUE,
                          clearBadReps = FALSE,
                          objPtOffset = 3.6 )
{
  nSims        <- obj$nSims
  goodReps_isp <- obj$goodReps_isp[1:nSims,,]
  
  SB_ispt   <- obj$om$SB_ispt[1:nSims,,,,drop = FALSE]
  C_ispt    <- obj$om$C_ispt[1:nSims,,,,drop = FALSE]


  # Pull ref pts
  rp <- obj$rp[[1]]

  # Pull Bmsy
  BmsySS_sp <- rp$FmsyRefPts$BeqFmsy_sp
  BmsyMS_sp <- rp$EmsyMSRefPts$BeqEmsy_sp


  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nReps   <- dim(SB_ispt)[1]

  # Get reference points
  B0_sp     <- obj$ctlList$opMod$histRpt$B0_sp
  BmsySS_sp <- obj$rp[[1]]$FmsyRefPts$BeqFmsy_sp
  BmsyMS_sp <- obj$rp[[1]]$EmsyMSRefPts$BeqEmsy_sp

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

  if( clearBadReps )
  {
    for( s in 1:nS )
      for(p in 1:nP )
      {
        badIdx <- which(!goodReps_isp[,s,p])
        SB_ispt[badIdx,s,p,] <- NA
        C_ispt[badIdx,s,p,] <- NA
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

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")

  par(  mfcol = c(nP,nS), 
        mar = c(.1,1.5,.1,1.5),
        oma = c(4,3,3,3) )


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
        rmtext( txt = stockNames[p], outer = TRUE,
                font = 2, line = 0.05, cex = 1.5 )
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
      par(xpd = TRUE)
        points( x = yrs[length(yrs)]+objPtOffset,
                y = BmsyMS_sp[s,p], col = "steelblue",
                pch = 16, cex = 1.5, lwd = 2 )
        points( x = yrs[length(yrs)]+objPtOffset,
                y = BmsySS_sp[s,p], col = "darkgreen",
                pch = 1, cex = 1.5, lwd = 2)
      par(xpd = FALSE)

      abline( v = yrs[tMP], col = "grey30", lty = 3 )
      abline( h = B0_sp[s,p], lty = 2, col = "grey50", lwd = 2  )

      if( mfg[1] == 1 & mfg[2] == 1 & leg )
        legend( x = "topright", bty = "n",
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
  if( printStamp )
    mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
              text = stamp, col = "grey60" )
}

plotTulipRt <- function(  obj = blob,
                          nTrace = 3, 
                          clearBadReps = FALSE )
{

  nSims         <- obj$nSims
  goodReps_isp  <- obj$goodReps_isp[1:nSims,,]
  
  R_ispt        <- 2*obj$om$R_ispt[1:nSims,,,,drop = FALSE]

  R_qspt        <- apply( X = R_ispt, FUN = quantile, MARGIN = c(2,3,4),
                          probs = c(0.025, 0.5, 0.975))

  # Pull ref pts
  rp <- obj$rp[[1]]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nReps   <- dim(R_ispt)[1]

  # Get reference points
  R0_sp     <- obj$ctlList$opMod$histRpt$R0_sp
  
  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

   par( mfcol = c(nP,nS), 
        mar = c(.1,1.5,.1,1.5),
        oma = c(4,4,3,3) )

   for( s in 1:nS )
    for( p in 1:nP )
    {
      maxR <- max(R_qspt[,s,p,],na.rm = T) 
      plot( x = range(yrs), y = c(0,maxR), type = "n",
            axes = FALSE )
        mfg <- par("mfg")
        if(mfg[1] == mfg[3])
          axis( side = 1)
        axis( side = 2, las = 1)
        grid()
        box()
        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2)

        if( mfg[2] == mfg[4] )
          rmtext( txt = stockNames[p], outer = TRUE, line = 0.05,
                  font = 2, cex = 1.5)
        

        polygon(  x = c(yrs,rev(yrs)), 
                  y = c(R_qspt[1,s,p,],rev(R_qspt[3,s,p,]) ),
                  border = NA, col = "grey60" )
        lines( x = yrs, y = R_qspt[2,s,p,], lwd = 2)
        points( x = yrs, y = R_qspt[2,s,p,], pch = 16, col = "grey15" )

        abline( h = R0_sp[s,p], lty = 2, lwd = 2)
        abline( v = yrs[tMP]-.5, lty = 3, lwd = 1)
        

    }
} # END plotTulipRt

plotTulipEconYield <- function( obj = blob,
                                nTrace = 3, 
                                fIdx = 2,
                                clearBadReps = FALSE,
                                colAlpha = 0.5 )
{

  nSims         <- obj$nSims
  goodReps_isp  <- obj$goodReps_isp[1:nSims,,]

  nReps         <- sum(apply(X = goodReps_isp, FUN=prod, MARGIN = c(1)))
  
  Rev_ispt      <- obj$om$Rev_ispft[1:nSims,,,2,]
  effCost_ipt   <- obj$om$effCost_ipft[1:nSims,,2,]
  crewShare     <- obj$ctlList$opMod$crewShare

  Rev_ipt       <- apply( X = Rev_ispt, FUN = sum, MARGIN = c(1,3,4), 
                          na.rm = T )
  Prof_ipt      <- Rev_ipt - effCost_ipt

  # Make Polygons
  Rev_qpt       <- apply( X = Rev_ipt, FUN = quantile, MARGIN = c(2,3),
                          probs = c(0.025, 0.5, 0.975), na.rm = T )
  effCost_qpt   <- apply( X = effCost_ipt, FUN = quantile, MARGIN = c(2,3),
                          probs = c(0.025, 0.5, 0.975), na.rm = T )
  Prof_qpt      <- apply( X = Prof_ipt, FUN = quantile, MARGIN = c(2,3),
                          probs = c(0.025, 0.5, 0.975), na.rm = T )




  # Pull steady state ref pts
  rp      <- obj$rp[[1]]
  Emey_p  <- rp$EmeyRefPts$Emey_p
  MEY_p   <- rp$EmeyRefPts$MEY_p

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nReps   <- dim(Rev_ipt)[1]
  
  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  par(  mfcol = c(nP,1), 
        mar = c(.1,1,.1,1),
        oma = c(4,4,3,3) )

  
  revCol  <- scales::alpha("grey30",alpha = colAlpha)
  costCol <- scales::alpha("red",alpha = colAlpha)
  profCol <- scales::alpha("steelblue",alpha = colAlpha)

  for( p in 1:nP )
  {
    maxRev <- max(Rev_qpt[,p,],na.rm = T) 
    plot( x = range(yrs[tMP:nT]), y = c(0,maxRev), type = "n",
          axes = FALSE )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis( side = 1)
      axis( side = 2, las = 1)
      grid()
      box()

      if( mfg[2] == mfg[4] )
        rmtext( txt = stockNames[p], outer = TRUE, line = 0.05,
                font = 2, cex = 1.5)
      

      polygon(  x = c(yrs,rev(yrs)), 
                y = c(Rev_qpt[1,p,],rev(Rev_qpt[3,p,]) ),
                border = NA, col = revCol )
      lines( x = yrs, y = Rev_qpt[2,p,], col = "black", lwd = 2)
      
      polygon(  x = c(yrs,rev(yrs)), 
                y = c(effCost_qpt[1,p,],rev(effCost_qpt[3,p,]) ),
                border = NA, col = costCol )
      lines( x = yrs, y = effCost_qpt[2,p,], col = "red", lwd = 2)
      
      polygon(  x = c(yrs,rev(yrs)), 
                y = c(Prof_qpt[1,p,],rev(Prof_qpt[3,p,]) ),
                border = NA, col = profCol )
      lines( x = yrs, y = Prof_qpt[2,p,], col = "steelblue", lwd = 2)


      abline( h = MEY_p[p], lty = 2, lwd = 2)
      abline( v = yrs[tMP]-.5, lty = 3, lwd = 1)
      

  }
} # END plotTulipRt

# plotTulipEffort_p()
# Effort over time gridded
# by stock area - envelopes
plotTulipEffort_p <- function(  obj = blob, 
                                nTrace = 3, 
                                fIdx = 1:2,
                                combineFisheries = TRUE,
                                envelopeFleet = 2,
                                clearBadReps = FALSE,
                                plotEmsy = TRUE )
{
  nSims         <- obj$nSims
  goodReps_isp  <- obj$goodReps_isp[1:nSims,,]

  E_ipft <- obj$om$E_ipft[1:nSims,,fIdx,,drop = FALSE ]

  

  if( clearBadReps )
  {
    for(p in 1:nP )
    {
      badIdx <- which(!goodReps_isp[,s,p])
      E_ipft[badIdx,p,,] <- NA
    }
  }

  # Pull ref pts
  rp <- obj$rp[[1]]

  # Pull Bmsy
  EmsyMS_p <- rp$EmsyMSRefPts$EmsyMS_p


  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(E_ipft)[1]

  if( combineFisheries )
  {
    newE_ipft       <- array( NA, dim = c(nSims,nP,1,nT))
    sumEff          <- apply(X = E_ipft, FUN = sum, MARGIN = c(1,2,4), na.rm = T )
    newE_ipft[,,1,] <- sumEff
    E_ipft <- newE_ipft
  }
  E_ipft[E_ipft == 0] <- NA

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

  if( envelopeFleet > 0 & combineFisheries )
    envelopeFleet <- 1

  if( combineFisheries )
    fleetCols <- "gray65"


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
        rmtext( txt = stockNames[p], line = 0.02,
                outer = TRUE, font = 2, cex = 1.5 )
      }
      box()
      grid()
      for( f in 1:dim(E_qpft)[3] )
      {
        if(f == envelopeFleet )
        {
          polygon( x = c(yrs,rev(yrs)), y = c(E_qpft[1,p,f,],rev(E_qpft[3,p,f,])), 
                  col = fleetCols[f], border = NA )
          for( tIdx in traces )
            lines( x = yrs, y = E_ipft[tIdx,p,f,], col = "black", lwd = .8 )
        }
        lines( x = yrs, y = E_qpft[2,p,f,], col = "black", lwd = 3)
        
      }
      par(xpd = TRUE)
        points( x = yrs[length(yrs)]+5.5,
                y = EmsyMS_p[p], col = "steelblue",
                pch = 16, cex = 1.5, lwd = 2 )
      par(xpd = FALSE)
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
  }
  mtext( side = 2, outer = TRUE, text = "Commercial Trawl Effort (1000 hrs)",
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
plotRetroSB <- function(  obj = blob, iRep = 1, 
                          totB = FALSE,
                          TAC = FALSE )
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
      if( totB )
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
plotRetroSBagg <- function( obj = blob, iRep = 1, 
                            Ct = TRUE,
                            totB = FALSE,
                            TAC = TRUE,
                            plotStamp = TRUE )
{
  # Get biomass arrays
  SB_spt        <- obj$om$SB_ispt[iRep,,,]
  VB_spt        <- obj$om$vB_ispft[iRep,,,2,]
  totB_spt      <- obj$om$B_ispt[iRep,,,]
  retroSB_tspt  <- obj$mp$assess$retroSB_itspt[iRep,,,,]

  pdHess_tsp    <- obj$mp$assess$pdHess_itsp[iRep,,,]

  ctlList <- obj$ctlList

  # Pull ref pts
  rp <- obj$rp[[1]]

  # Pull Bmsy
  BmsySS_sp <- rp$FmsyRefPts$BeqFmsy_sp
  BmsyMS_sp <- rp$EmsyMSRefPts$BeqEmsy_sp

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
  if( ctlList$mp$data$spatialPooling )
  {
    newSB_spt        <- apply( X = SB_spt, FUN = sum, MARGIN = c(1,3), na.rm = T )
    newVB_spt        <- apply( X = VB_spt, FUN = sum, MARGIN = c(1,3), na.rm = T )
    newtotB_spt      <- apply( X = totB_spt, FUN = sum, MARGIN = c(1,3), na.rm = T )

    newBmsySS_sp     <- apply( X = BmsySS_sp, FUN = sum, MARGIN = c(1),na.rm = T )
    newBmsyMS_sp     <- apply( X = BmsyMS_sp, FUN = sum, MARGIN = c(1),na.rm = T )

    BmsySS_sp  <- BmsySS_sp[,1,drop = FALSE]
    BmsySS_sp[,1]  <- newBmsySS_sp
    BmsyMS_sp  <- BmsyMS_sp[,1,drop = FALSE]
    BmsyMS_sp[,1]  <- newBmsyMS_sp


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

  if( ctlList$mp$data$speciesPooling )
  {
    newSB_spt        <- apply( X = SB_spt, FUN = sum, MARGIN = c(2,3), na.rm = T )
    newVB_spt        <- apply( X = VB_spt, FUN = sum, MARGIN = c(2,3), na.rm = T )
    newtotB_spt      <- apply( X = totB_spt, FUN = sum, MARGIN = c(2,3), na.rm = T )

    newBmsySS_sp     <- apply( X = BmsySS_sp, FUN = sum, MARGIN = c(2),na.rm = T )
    newBmsyMS_sp     <- apply( X = BmsyMS_sp, FUN = sum, MARGIN = c(2),na.rm = T )

    BmsySS_sp  <- BmsySS_sp[1,,drop = FALSE]
    BmsySS_sp[1,]  <- newBmsySS_sp
    BmsyMS_sp  <- BmsyMS_sp[1,,drop = FALSE]
    BmsyMS_sp[1,]  <- newBmsyMS_sp

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
    stockNames <- "Spatial Pooled"

  if( nSS == 1 )
    speciesNames <- "Species Pooled"

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,":",iRep,sep = "")

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  par(  mfcol = c(nPP,nSS), 
        mar = c(1,1.5,1,1.5),
        oma = c(4,4,3,3) )
  for(s in 1:nSS)
    for( p in 1:nPP )
    {
      plot( x = range(yrs),
            y = c(0,max(VB_spt[s,p,],SB_spt[s,p,],retroSB_tspt[,s,p,],na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
        rmtext( line = .05, txt = stockNames[p], cex = 1.5, font = 2)
      box()
      grid()
      if( Ct )
      {
        # Plot actual catch
        rect( xleft = yrs - .3, xright = yrs + .3, 
              ytop = C_spt[s,p,], ybottom = 0, col = "grey40",
              border = NA )
        # Plot a rectangle of TACs
        if( TAC )
          rect( xleft = yrs[tMP:nT] - .3, xright = yrs[tMP:nT] + .3, 
                ytop = TAC_spt[s,p,tMP:nT], ybottom = 0, col = NA,
                lwd = .5,
                border = "black" )

      }
      lines( x = yrs, y = SB_spt[s,p,], col = "red", lwd = 3 )
      lines( x = yrs, y = VB_spt[s,p,], col = "grey40", lwd = 2, lty = 3 )
      if( totB )
        lines( x = yrs, y = totB_spt[s,p,], col = "black", lwd = 2 )
      for( tt in 1:pT )
      {
        if( pdHess_tsp[tt,s,p] )
          lineCol <- "grey60"
        if( !pdHess_tsp[tt,s,p] )
          lineCol <- "purple"
        lines( x = yrs, y = retroSB_tspt[tt,s,p,], col = lineCol, lwd = 1 )
      }
      par(xpd = TRUE)
        points( x = yrs[length(yrs)]+3.6,
                y = BmsyMS_sp[s,p], col = "steelblue",
                pch = 16, cex = 1.5, lwd = 2 )
        points( x = yrs[length(yrs)]+3.6,
                y = BmsySS_sp[s,p], col = "darkgreen",
                pch = 1, cex = 1.5, lwd = 2)
      par(xpd = FALSE)
  
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }
  mtext( side = 1, text = "Year", outer = TRUE, line = 2, font = 2)
  mtext( side = 2, text = "Spawning Biomass, TAC, and Catch (kt)", 
          outer = TRUE, line = 1.5, font = 2)
  if( plotStamp )
    mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
              text = stamp, col = "grey60" )

}

# plotEBSBratio()
# Plots of ratio between vulnerable
# biomass and exploitable biomass
# for given depletion levels.
plotEBSBratio <- function( obj = blob )
{
  # Pull model dims
  nS <- obj$om$nS
  nP <- obj$om$nP
  nF <- obj$om$nF

  fleetNames <- obj$om$fleetNames

  # Pull reference points
  rp <- obj$rp[[1]]

  nu_spfk <- array(NA, dim = c(nS+1,nP+1,nF,3))
  P1_spf  <- array(NA, dim = c(nS+1,nP+1,nF))
  P2_spf  <- array(NA, dim = c(nS+1,nP+1,nF))

  # Fill!
  # nu pars
  nu_spfk[1:nS,1:nP,1:nF,] <- rp$EBSBpars$stockSpec$nu_spfk
  nu_spfk[nS+1,1:nP,1:nF,] <- rp$EBSBpars$coastWide$nu_spfk
  nu_spfk[1:nS,1+nP,1:nF,] <- rp$EBSBpars$dataPooled$nu_spfk
  nu_spfk[1+nS,1+nP,1:nF,] <- rp$EBSBpars$totalAgg$nu_spfk

  # P1
  P1_spf[1:nS,1:nP,1:nF] <- rp$EBSBpars$stockSpec$P1_spf
  P1_spf[nS+1,1:nP,1:nF] <- rp$EBSBpars$coastWide$P1_spf
  P1_spf[1:nS,1+nP,1:nF] <- rp$EBSBpars$dataPooled$P1_spf
  P1_spf[1+nS,1+nP,1:nF] <- rp$EBSBpars$totalAgg$P1_spf

  # P2
  P2_spf[1:nS,1:nP,1:nF] <- rp$EBSBpars$stockSpec$P2_spf
  P2_spf[nS+1,1:nP,1:nF] <- rp$EBSBpars$coastWide$P2_spf
  P2_spf[1:nS,1+nP,1:nF] <- rp$EBSBpars$dataPooled$P2_spf
  P2_spf[1+nS,1+nP,1:nF] <- rp$EBSBpars$totalAgg$P2_spf

  depSeq <- seq(0,1,length.out = 100 )

  par( mfcol = c(nP+1,nS+1),
        mar = c(.5,.5,.5,.5),
        oma = c(3,3,2,2) )


  fleetCols <- RColorBrewer::brewer.pal(nF,"Dark2")

  for( s in 1:(nS+1) )
    for( p in 1:(nP+1) )
    {
      plot( x = c(0,1), y = c(0,2), axes = FALSE,
            type = "n" )
        mfg <- par( "mfg" )

        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        if( mfg[2] == 1 )
          axis( side = 2, las = 1 )

        box()

        abline( h = 1, lwd = .8, lty = 2 )

        for( f in 1:nF )
        {
          numerator   <- (1 - exp( -nu_spfk[s,p,f,3]*( depSeq - P1_spf[s,p,f]))) 
          denominator <- (1 - exp( -nu_spfk[s,p,f,3]*( P2_spf[s,p,f] - P1_spf[s,p,f]))) 

          ratioSeq    <- nu_spfk[s,p,f,1] + (nu_spfk[s,p,f,2] - nu_spfk[s,p,f,1]) * numerator/denominator


          lines( x = depSeq, y = ratioSeq, lwd = 2, col = fleetCols[f] )
        }
        
    }

  legend( x = "topleft",
          legend = fleetNames,
          col = fleetCols,
          lwd = 2, bty = "n" )

} # END plotEBSBratio()

# plotRetroCatchability()
# Plot of retrospective catchability estimates for each 
# fleet compared to the mean taken from the OM
plotRetroCatchability <- function(  obj = blob,
                                    iRep = 1 )
{
  # Get model dims
  tMP <- obj$om$tMP
  nT  <- obj$om$nT
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nF  <- obj$om$nF

  # Get model state arrays and AM fits
  SB_spt        <- obj$om$SB_ispt[iRep,,,1:t]
  VB_spt        <- obj$om$vB_ispft[iRep,,,2,1:t]
  totB_spt      <- obj$om$B_ispt[iRep,,,1:t]
  fitSB_spt     <- obj$mp$assess$retroSB_itspt[iRep,projt,,,1:t]
  fitVB_spft    <- obj$mp$assess$retroVB_itspft[iRep,projt,,,,1:t]
  fitq_spf      <- obj$mp$assess$retroq_itspf[iRep,,,,]
  fitq_spft     <- obj$mp$assess$retroq_itspft[iRep,,,,,]

  # Get ctlList
  ctlList <- obj$ctlList
  pT      <- ctlList$opMod$pT

  # Get mean catchability (OM fits)
  mq_spf      <- ctlList$opMod$histRpt$q_spf
  mq_sf       <- array(1, dim = c(nS,nF))
  mq_sf[,3:4] <- ctlList$opMod$histRpt$qSurv_sf
  sdlnq_f     <- ctlList$mp$assess$spsdlnq_f


}



# plotScaledIndices()
# Scaled indices with model fits for a given
# replicate and year. 
plotScaledIndices <- function(  obj = blob, 
                                iRep = 1, 
                                Ct = TRUE,
                                t = blob$om$tMP,
                                totB = FALSE )
{
  # Get model dims
  tMP <- obj$om$tMP
  nT  <- obj$om$nT
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nF  <- obj$om$nF

  # Calculate the projection time step
  projt <- t - tMP + 1

  # Get biomass arrays
  SB_spt        <- obj$om$SB_ispt[iRep,,,1:t]
  VB_spt        <- obj$om$vB_ispft[iRep,,,2,1:t]
  totB_spt      <- obj$om$B_ispt[iRep,,,1:t]
  fitSB_spt     <- obj$mp$assess$retroSB_itspt[iRep,projt,,,1:t]
  fitVB_spft    <- obj$mp$assess$retroVB_itspft[iRep,projt,,,,1:t]
  fitq_spf      <- obj$mp$assess$retroq_itspf[iRep,projt,,,]
  fitq_spft     <- obj$mp$assess$retroq_itspft[iRep,projt,,,,1:t]

  ctlList <- obj$ctlList
  spFleets <- ctlList$mp$assess$spFleets


  fitSB_spt[fitSB_spt < 0] <- NA

  C_spft   <- obj$om$C_ispft[iRep,,,,1:t]

  C_spt    <- apply(X = C_spft, FUN = sum, MARGIN = c(1,2,4))

  # Now pull indices
  I_spft <- obj$mp$data$I_ispft[iRep,,,,1:t]
  I_spft[,,-spFleets,] <- -1
  I_spft[I_spft < 0] <- NA



  # Aggregate OM biomasses to match AM biomass
  if( ctlList$mp$data$spatialPooling & !ctlList$mp$data$speciesPooling )
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

    I_spft <- I_spft[1:nS,1+nP,,,drop = FALSE]

    sumC_spt          <- apply(X = C_spft, FUN = sum, MARGIN = c(1,4))
    newC_spt          <- array(NA, dim = c(nS,1,t))
    
    newC_spt[,1,]     <- sumC_spt
    C_spt             <- newC_spt

  }

  if( ctlList$mp$data$speciesPooling &  !ctlList$mp$data$spatialPooling )
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

    I_spft <- I_spft[1+nS,1:nP,,, drop= FALSE]

    sumC_spt          <- apply(X = C_spft, FUN = sum, MARGIN = c(2,4))
    newC_spt          <- array(NA, dim = c(1,nP,t))
    newC_spt[1,,]     <- sumC_spt
    C_spt             <- newC_spt

  }

  if( ctlList$mp$data$speciesPooling &  ctlList$mp$data$spatialPooling )
  {
    newSB_spt        <- apply( X = SB_spt, FUN = sum, MARGIN = c(3), na.rm = T )
    newVB_spt        <- apply( X = VB_spt, FUN = sum, MARGIN = c(3), na.rm = T )
    newtotB_spt      <- apply( X = totB_spt, FUN = sum, MARGIN = c(3), na.rm = T )

    SB_spt    <- SB_spt[1,1,,drop = FALSE]
    SB_spt[1,1,] <- newSB_spt
    VB_spt    <- VB_spt[1,1,,drop = FALSE]
    VB_spt[1,1,] <- newVB_spt
    totB_spt  <- totB_spt[1,1,,drop = FALSE]
    totB_spt[1,1,] <- newtotB_spt

    I_spft <- I_spft[1+nS,1+nP,,, drop= FALSE]

    sumC_spt          <- apply(X = C_spft, FUN = sum, MARGIN = c(4))
    newC_spt          <- array(NA, dim = c(1,1,t))
    newC_spt[1,1,]    <- sumC_spt
    C_spt             <- newC_spt

  }

  tFirstIdx <- ctlList$mp$assess$spYrFirstIdx

  if( tFirstIdx > 1 )
    I_spft[,,,1:tFirstIdx] <- NA


  
  fleetPCH  <- 20 + 1:nF
  fleetBG   <- RColorBrewer::brewer.pal(nF, "Set1")
  stockCol  <- RColorBrewer::brewer.pal(nP, "Spectral")


  spTVqFleets <- ctlList$mp$assess$spTVqFleets

  SB_spt[SB_spt == 0]     <- NA
  VB_spt[VB_spt == 0]     <- NA
  totB_spt[totB_spt == 0] <- NA

  nSS     <- dim( SB_spt)[1]
  nPP     <- dim( SB_spt)[2]

  scaledIdx_spft <- array( NA, dim = c(nSS, nP, nF, t ) )
  for( s in 1:nSS )
    for( p in 1:nPP )
      for( f in 1:nF )
      {
        if( f %in% spTVqFleets)
          scaledIdx_spft[s,p,f,] <- I_spft[s,p,f,] / fitq_spft[s,p,f,] 
        else
          scaledIdx_spft[s,p,f,] <- I_spft[s,p,f,] / fitq_spf[s,p,f]

        # Scale by ratio of SB and VB
        scaledIdx_spft[s,p,f,] <- scaledIdx_spft[s,p,f,] * fitSB_spt[s,p,] / fitVB_spft[s,p,f,]
      }

  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear
  pT            <- obj$ctlList$opMod$pT

  if( nPP == 1 )
    stockNames <- "Spatially Pooled"

  if( nSS == 1 )
    speciesNames <- "DER Complex"

  yrs <- seq( from = fYear, by = 1, length.out = t)
  ppJitter <- seq(from = -.3, to = .3, length.out = nP )

  par(  mfcol = c(nPP,nSS), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )
  for(s in 1:nSS)
    for( p in 1:nPP )
    {
      plot( x = range(yrs),
            y = c(0,max(totB_spt[s,p,],VB_spt[s,p,],SB_spt[s,p,],na.rm = T) ),
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
      if( totB )
        lines( x = yrs, y = totB_spt[s,p,], col = "black", lwd = 2 )
      
      lines( x = yrs, y = fitSB_spt[s,p,], col = "grey60", lwd = 1 )

      for( f in 1:nF )
        points( x = yrs, y = scaledIdx_spft[s,p,f,],
                pch = fleetPCH[f], bg = fleetBG[f],
                cex = 1.3, col = stockCol[p] )
      
      if( Ct )
      {
        # Plot actual catch
        rect( xleft = yrs - .3, xright = yrs + .3, 
              ytop = C_spt[s,p,], ybottom = 0, col = "grey40",
              border = NA )

      }
  
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }
    legend( x= "topright",
            legend = fleetNames,
            pch = fleetPCH,
            pt.bg = fleetBG )

} # END plotScaledIndices()

# plotAMIdxResids()
# Plots of standardised AM residuals
plotAMIdxResids <- function(  obj = blob, 
                              iRep = 1, 
                              Ct = TRUE,
                              t = blob$om$tMP )
{
  # Get model dims
  tMP <- obj$om$tMP
  nT  <- obj$om$nT
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nF  <- obj$om$nF

  # Calculate the projection time step
  projt <- t - tMP + 1

  # Get biomass arrays
  SB_spt        <- obj$om$SB_ispt[iRep,,,1:t]
  VB_spt        <- obj$om$vB_ispft[iRep,,,2,1:t]
  totB_spt      <- obj$om$B_ispt[iRep,,,1:t]
  fitSB_spt     <- obj$mp$assess$retroSB_itspt[iRep,projt,,,1:t]
  fitVB_spft    <- obj$mp$assess$retroVB_itspft[iRep,projt,,,,1:t]
  fitq_spf      <- obj$mp$assess$retroq_itspf[iRep,projt,,,]
  fitq_spft     <- obj$mp$assess$retroq_itspft[iRep,projt,,,,1:t]
  tauObs_spf    <- obj$mp$assess$retrotauObs_itspf[iRep,projt,,,]

  tauObs_spf[is.na(tauObs_spf)] <- 0

  ctlList <- obj$ctlList

  fitSB_spt[fitSB_spt < 0] <- NA

  firstIdx <- ctlList$mp$assess$spYrFirstIdx

  spTVqFleets <- ctlList$mp$assess$spTVqFleets
  spFleets    <- ctlList$mp$assess$spFleets

  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT

  # Now pull indices
  I_spft <- obj$mp$data$I_ispft[iRep,,,,1:t]
  
  I_spft[,,-spFleets,] <- -1
  
  if(firstIdx > 1)
    I_spft[,,,1:firstIdx] <- -1
  
  I_spft[I_spft <= 0] <- NA

  nSS <- nS
  nPP <- nP

  # Aggregate OM biomasses to match AM biomass
  if( ctlList$mp$data$spatialPooling & !ctlList$mp$data$speciesPooling )
  {
    nPP    <- 1
    I_spft <- I_spft[1:nS,nP+1,,,drop = FALSE]
  }

  if( ctlList$mp$data$speciesPooling & !ctlList$mp$data$spatialPooling  )
  {
    nSS <- 1
    I_spft <- I_spft[1+nS,1:nP,,,drop = FALSE] 
  }

  if( ctlList$mp$data$speciesPooling & ctlList$mp$data$spatialPooling )
  {
    nSS <- 1
    nPP <- 1
    I_spft <- I_spft[1+nS,1+nP,,,drop = FALSE] 
  }


  stdResids_spft <- array( NA, dim = c(nSS, nPP, nF, t) )

  for( s in 1:nSS )
    for( p in 1:nPP )
      for( f in spFleets )
      {
        if( tauObs_spf[s,p,f] > 0 )
        {

          if( f %in% spTVqFleets )
            stdResids_spft[s,p,f,] <- (-log( I_spft[s,p,f,]/fitq_spft[s,p,f,] ) + log(fitVB_spft[s,p,f,]))/tauObs_spf[s,p,f]
          if( !f %in% spTVqFleets ) 
            stdResids_spft[s,p,f,] <- (-log( I_spft[s,p,f,]/fitq_spf[s,p,f] ) + log(fitVB_spft[s,p,f,]))/tauObs_spf[s,p,f]


        }
      }

  
  fleetPCH <- 20 + 1:nF
  fleetBG <- RColorBrewer::brewer.pal(nF, "Set1")



  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear
  pT            <- obj$ctlList$opMod$pT


  if( nSS == 1 )
    speciesNames <- "Data Pooled"

  yrs <- seq( from = fYear, by = 1, length.out = t)

  ppJitter <- seq( from = -.3, to = .3, length.out = nP )

  par(  mfcol = c(nPP,nSS), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )
  for(s in 1:nSS)
    for( p in 1:nPP )
    {
      maxResid <- max(abs(stdResids_spft[s,p,,]),na.rm = T)
      if(!is.finite(maxResid))
        maxResid <- 1
      plot( x = range(yrs),
            y = range(-maxResid,maxResid),
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

      for( f in 1:nF )
      {
        points( x = yrs, y = stdResids_spft[s,p,f,],
                pch = fleetPCH[f], bg = fleetBG[f],
                cex = 1.3 )
        nonNA <- which(!is.na(stdResids_spft[s,p,f,]))
        if( length(nonNA) > 0 )
        {
          yVal <- stdResids_spft[s,p,f,nonNA]
          xVal <- yrs[nonNA]



          dat <- data.frame(x = xVal, y = yVal )
          regLine <- lm( y~x, data = dat )

          pVal <- round(summary(regLine)$coefficients[2,4],3)

          pLabel <- paste("p = ", pVal, sep = "")

          dat <- dat %>%
                  mutate( regLine = predict.lm(regLine, newdata = dat) )


          lines( x = dat$x, y = dat$regLine, col = fleetBG[f], lwd = 2 )
          text( x = dat$x[1], y = 1.2*dat$y[1], label = pLabel, col = fleetBG[f], font = 2 )

        }
      }
      abline( h = 0, lty = 2)
      
      
    }
    
    legend( x= "bottomleft", bty = "n",
            legend = c(fleetNames),
            pch = fleetPCH,
            pt.bg = fleetBG,
            cex = 1.3 )



} # END plotAMIdxResids()


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


  library(wesanderson)
  specCols <- wes_palette("Rushmore1", 5, type = "discrete")[3:5]

  true_spt <- repObj[[AMseries]]
  est_spt  <- omObj[[OMseries]][iRep,,,]

  re_qspt  <- calcREdist( true = true_spt,
                          est  = est_spt[,,1:nT],
                          marg = c(1,2,3) )

  # specCols <- wes_palette("Rushmore1", 5, type = "discrete")[3:5]
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

  specCols <- wes_palette("Rushmore1", 5, type = "discrete")[3:5]
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

  specCols <- wes_palette("Rushmore1", 5, type = "discrete")[3:5]
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


