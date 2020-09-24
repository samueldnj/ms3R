# Plotting functions for assessCA TMB model
# assessCAplots.R

library(reshape2)

# rmtext()
# Refactored procedure to plot right hand inner
# margin mtext with the bottom towards the middle
# of the plot
rmtext <- function( line = 1, 
                    txt = "Sample", 
                    font = 1,
                    cex = 1,
                    outer = FALSE,
                    yadj = .5)
{
  corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
  if( outer )
    par(xpd = NA) #Draw outside the figure region
  if( !outer )
    par(xpd = TRUE)

  xRange <- corners[2] - corners[1]


  text( x = corners[2]+line*xRange, 
        y = yadj * sum(corners[3:4]), 
        labels = txt, srt = 270,
        font = font, cex = cex )
  par(xpd = FALSE)
} # END rmtext()

# Plot multi-panel time series plots
plotMulti <- function(  ts = c("Nt","Bt","Ft"),
                        repList = reports,
                        labcex = 0.8, heading = "dimLabs" )
{

  nP <- repList$repOpt$nP
  par(mfcol = c(length(ts),nP), oma = c(3,4,3,4), mar = c(0.2,3,0.2,2),
        cex.lab = labcex )

  argList <- list(  repList = repList,  
                    noPar = TRUE, 
                    labcex = labcex )

  for( pIdx in 1:nP)  
  {
    argList$pIdx = pIdx
    for( tsName in ts )
    {
      if( tsName %in% c("SBt","SBtIdx") )
        plotArgs <- c(argList, list(heading = heading) )
      else plotArgs <- argList
      plotCall <- paste("plot",tsName,sep = "")
      do.call(  what = eval(plotCall), 
                args = plotArgs, 
                envir = .GlobalEnv )
    }
  }
}



# Plot biomass
plotSBt <- function(  repList = reports,  
                      noPar = FALSE, 
                      pIdx = 1,
                      labcex = .8,
                      heading = NULL )
{
  report    <- repList$repOpt
  initYear  <- repList$fYear


  # Pull stuff from report
  SBt       <- report$SB_pt[pIdx,]
  Cgt       <- report$totC_pgt[pIdx,,]

  B0        <- signif(report$B0_p[pIdx],3)
  M         <- signif(report$M_p[pIdx],3)

  plotCI <- FALSE

  if(!is.null(repList$posts))
  {
    plotCI <- TRUE
    tInitModel_p <- repList$repOpt$tInitModel_p
    SB_it <- repList$posts$SB_ipt[,pIdx,]

    SB_qt <- apply(X = SB_it, FUN = quantile,
                    MARGIN = c(2), probs = c(0.025, 0.5, 0.975),
                    na.rm = TRUE )

    SBt <- apply(X = SB_it, FUN = mean, MARGIN = c(2))

    B0 <- signif(mean(repList$posts$B0_ip[,pIdx]),3)
    M  <- signif(mean(repList$posts$M_iapt[,1,pIdx,tInitModel_p[pIdx]+ 1]),3)

  } else {
    SB_qt <- 0
  }

  SBt[SBt == 0] <- NA

  # Sum catch
  Ct        <- apply(X = Cgt, FUN = sum, MARGIN = 2)

  stockLabs <- dimnames(report$SB_pt)[[1]][pIdx]

  if( !is.null(heading))
  {
    if( heading == "dimLabs")
      stockLabel <- stockLabs
    else
      stockLabel <- try(heading[stockLabs])

    if( class(stockLabel) == "try-error" ) 
      stockLabel <- stockLabs
  }

  
  # Get number of time steps
  nT    <- report$nT
  years <- seq(from = initYear, length = nT+1, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  depRange  <- range(SBt/B0, na.rm = T)
  maxDep    <- ceiling(depRange[2])

  depAxis   <- round(seq(0, maxDep, length.out = 5),3)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot recruitments
  plot( x = range(years), y = c(0,max(SBt,SB_qt, na.rm =T) ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    mfg <- par("mfg")
    # Add SSB axis
    axis(side = 2, las =1 )
    # Add depletion axis
    axis( side = 4, las = 1, labels = depAxis, at = B0 * depAxis)

    if( mfg[1] == mfg[3] )
      axis(side = 1 )

    box()
    grid()
    rect( xleft = years[1:nT] - .3, xright = years[1:nT] + .3,
          ybottom = 0, ytop = Ct, border = NA, col = "grey60" )
    # Plot Biomass
    if( plotCI )
    {
      polyCol <- scales::alpha("red", .5)
      polygon(  x = c(years[1:nT],rev(years[1:nT])),
                y = c(SB_qt[1,1:nT], rev(SB_qt[3,1:nT])),
                col = polyCol, border = NA )
    }
    lines( x = years[1:nT], y = SBt[1:nT], lwd = 2, col = "red" )
    # abline( v = vertLines, lwd = .8, lty = 2, col = "grey80")
    abline( h = B0, lty = 3, lwd = .8, col = "red" )
    points(x = years[(nT+1)], y = SBt[(nT+1)], col = "black", bg = "red", pch = 21)
    
    if( mfg[2] == 1 )
      mtext(side = 2, text = "Spawning \n Biomass (kt)", line = 3.5, cex = labcex)
    if( mfg[2] == mfg[4] )
      mtext(side = 4, text = "Depletion", line = 3, cex = labcex)
    if( mfg[1] == 1 & !is.null(heading) )
      mtext(side = 3, text = stockLabel, line = 1, font = 2, cex = labcex )

    panLab( x = 0.7, y = 0.8, txt = paste( "B0 = ",  B0, sep = "") )
    panLab( x = 0.7, y = 0.7, txt = paste( "M = ",  M, sep = "") )
} # END plotSBt()

# Plot biomass
plotSBtIdx <- function( repList = reports,  
                        noPar = FALSE, 
                        pIdx = 1,
                        labcex = .8,
                        heading = NULL,
                        plotCt = FALSE )
{
  initYear <- repList$fYear
  report   <- repList$repOpt

  # Pull stuff from report
  SBt       <- report$SB_pt[pIdx,]
  Cgt       <- report$totC_pgt[pIdx,,]
  Igt       <- report$I_pgt[pIdx,,]
  qg        <- report$qhat_pg[pIdx,]
  rI_gt     <- repList$data$rI_pgt[pIdx,,]
  vulnBgt   <- report$vulnB_pgt[pIdx,,]
  vulnNgt   <- report$vulnN_pgt[pIdx,,]

  combIdx <- FALSE

  if( any(repList$data$whichCombIdx_g > 0) )
  {
    combIdx <- TRUE
    # Replace q values with qComb_g
    qg <- report$qComb_pg[pIdx,]
  }

  
  vulnBgt[ vulnBgt == 0 ] <- NA
  vulnNgt[ vulnNgt == 0 ] <- NA



  # Pull prob positive idx.
  probPosIdx_gt <- report$probPosIdx_pgt[pIdx,,]

  # Get survey type and whether
  # idx is relative or absolute
  surveyType  <- report$survType_g
  indexType   <- report$indexType_g
  calcIndex   <- report$calcIndex_g
  surveyGears <- which(calcIndex == 1)

  # Get number of time steps/gear types
  nT    <- report$nT
  nG    <- report$nG

  B0        <- signif(report$B0_p[pIdx],3)
  M         <- signif(report$M_p[pIdx],3)

  plotCI <- FALSE

  if(!is.null(repList$posts))
  {
    plotCI <- TRUE
    tInitModel_p <- repList$repOpt$tInitModel_p
    SB_it <- repList$posts$SB_ipt[,pIdx,]

    SB_qt <- apply(X = SB_it, FUN = quantile,
                    MARGIN = c(2), probs = c(0.025, 0.5, 0.975),
                    na.rm = TRUE )

    SBt <- apply(X = SB_it, FUN = mean, MARGIN = c(2))

    B0 <- signif(mean(repList$posts$B0_ip[,pIdx]),3)
    M  <- signif(mean(repList$posts$M_iapt[,1,pIdx,tInitModel_p[pIdx]+ 1]),3)

    q_ig  <- repList$posts$q_ipg[,pIdx,]
    qg   <- apply( X = q_ig, FUN = mean, MARGIN = 2)

  } else {
    SB_qt <- array(0, dim = c(3,nT+1))
  }



  SBt[ SBt == 0 ] <- NA
  SB_qt[SB_qt == 0] <- NA

  # Sum catch
  Ct        <- apply(X = Cgt, FUN = sum, MARGIN = 2)

  stockLabs <- dimnames(report$SB_pt)[[1]][pIdx]

  if( !is.null(heading))
  {
    if( heading == "dimLabs")
      stockLabel <- stockLabs
    else
      stockLabel <- try(heading[stockLabs])

    if( class(stockLabel) == "try-error" ) 
      stockLabel <- stockLabs
  }


  # if( stockLabs == "Aggregate" )
  # {
  #   # We need to pull the aggregated stock
  #   # indices in at this point
  #   mItg    <- report$mI_tg
  #   mqg     <- report$qhat_g

  #   mItg[mItg < 0] <- 0

  #   for( g in surveyGears )
  #   {
  #     if( all(Itg[,g] < 0) & any(mItg[,g] > 0))
  #     {
  #       Itg[,g] <- mItg[,g]
  #       qg[g]   <- mqg[g]
  #     }
  #   }
  # }

  # browser()

  Igt[Igt<0] <- NA

  scaledIndices <- Igt

  for( g in surveyGears )
  {
    posIdx <- which(Igt[g,1:nT] > 0)

    if( surveyType[g] == 1 )
      vulnTS <- rI_gt[g,1:nT] * vulnNgt[g,1:nT]
    
    if( surveyType[g] == 0 )
      vulnTS <- rI_gt[g,1:nT] * vulnBgt[g,1:nT]

    if( surveyType[g] == 2 )
      vulnTS <- rI_gt[g,1:nT] * SBt[1:nT]
    
    # Compute scaled indices
    scaledIndices[g,posIdx] <- Igt[g,posIdx] / qg[g]  * SBt[posIdx] / vulnTS[posIdx] / probPosIdx_gt[g,posIdx] 

  }

  tInitModel <- report$tInitModel_p[pIdx] + 1

  scaledCombIdx <- repList$data$combI_pt[pIdx,]
  scaledCombIdx[scaledCombIdx < 0 ] <- NA

  if( combIdx )
  {
    # # Need to use posteriors here
    # if(!is.null(reports$posts))
    # {
    #   probPosIdx_ipgt <- reports$posts$probPosCombIdx_ipgt
    #   nPosts <- dim(probPosIdx_ipgt)[1]
    #   probPosIdx_pgt  <- apply( X = probPosIdx_ipgt, FUN = mean, MARGIN = c(2,3,4))

    #   meanProbPosIdx_ipg  <- reports$posts$meanProbPosIdx_ipg
    #   SDProbPosIdx_ipg    <- reports$posts$SDProbPosIdx_ipg

      

    #   meanProbPosIdx_pg <- apply( X = reports$posts$meanProbPosIdx_ipg, FUN = mean, MARGIN = c(2,3))
    #   SDProbPosIdx_pg   <- apply( X = reports$posts$SDProbPosIdx_ipg  , FUN = mean, MARGIN = c(2,3))

    #   modelProbPos_igb <- array( 0, dim = c(nPosts,nG,100) )
    #   for(i in 1:nPosts )
    #   {
    #     for( g in 1:nG )
    #       modelProbPos_igb[i,g,] <- 1 / (1 + exp( - (meanProbPosIdx_ipg[i,1,g] + SDProbPosIdx_ipg[i,1,g] * BtSeq)))
    #   }

    #   modelProbPos_qgb <- apply( X = modelProbPos_igb, FUN = quantile, na.rm = T,
    #                               MARGIN = c(2,3), probs = c(0.025, 0.5, 0.975) )

    #   meanModelProb_gb <- apply( X = modelProbPos_igb, FUN = mean,
    #                               MARGIN = c(2,3), na.rm = T )
    # } 

    probPosCombIdx <- report$probPosCombIdx_pt[pIdx,]
    qComb_t <- report$qComb_pt[pIdx,]
    scaledCombIdx <- scaledCombIdx / qComb_t / probPosCombIdx
  }


  # Create x axis label vector and vertical lines
  # for easy multipanel plotting
  years <- seq(from = initYear, length = nT+1, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  cols <- brewer.pal(nG,"Dark2")

  depRange  <- range(SBt/B0, na.rm = T)
  maxDep    <- ceiling(depRange[2])

  depAxis   <- round(seq(0,maxDep, length.out = 5),2)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  yMax <- max(SBt,scaledIndices, na.rm =T)

  # Now plot recruitments
  plot( x = range(years), y = c(0,yMax ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    mfg <- par("mfg")
    # Add SSB axis
    axis(side = 2, las =1 )
    # Add depletion axis
    axis( side = 4, las = 1, labels = depAxis, at = B0 * depAxis)

    if( mfg[1] == mfg[3] )
      axis(side = 1 )

    box()
    grid()
    # Plot data
    if( combIdx )
    {
      posIdx  <- which(scaledCombIdx > 0)
      zeroIdx <- which(scaledCombIdx == 0)
      zeroIdx <- zeroIdx[zeroIdx >= tInitModel ]
      points( x = years[posIdx], y = scaledCombIdx[posIdx],
              col = "grey50", pch = 16)

      # points( x = years[zeroIdx], y = scaledCombIdx[zeroIdx],
      #         col = "grey40", pch = 1)
      axis( side = 3, at = years[zeroIdx], labels = FALSE,
            col.ticks = "grey40", tck = .02, lwd.ticks = 3 )
    }

    # Plot Biomass
    if( plotCI )
    {
      polyCol <- scales::alpha("red", .5)
      polygon(  x = c(years[1:nT],rev(years[1:nT])),
                y = c(SB_qt[1,1:nT], rev(SB_qt[3,1:nT])),
                col = polyCol, border = NA )
    }
    abline( h = B0, lty = 3, lwd = .8, col = "red" )
    lines( x = years[1:nT], y = SBt[1:nT], lwd = 2, col = "red" )
    if(plotCt)
      rect( xleft = years[1:nT] - .3, xright = years[1:nT] + .3,
            ybottom = 0, ytop = Ct, border = NA, col = "grey60" )
    # points(x = years[(nT+1)], y = SBt[(nT+1)], col = "black", bg = "red", pch = 21)
    
    if( mfg[2] == 1 )
      mtext(side = 2, text = "Spawning\nBiomass (kt)", line = 3.5, cex = labcex)
    if( mfg[2] == mfg[4] )
      rmtext(txt = "Depletion", line = .2, cex = 1.5, outer = TRUE)
    if( mfg[1] == 1 & !is.null(heading) )
      mtext(side = 3, text = stockLabel, line = 1, font = 2, cex = labcex )

    for( g in surveyGears )
    {
      posIdx  <- which(scaledIndices[g,] > 0)
      zeroIdx <- which(scaledIndices[g,] == 0)
      points( x = years[posIdx], y = scaledIndices[g,posIdx],
            col = alpha(cols[g],.5), pch = 16 )
      points( x = years[zeroIdx], y = scaledIndices[g,zeroIdx],
            col = alpha(cols[g],.5), pch = 1 )
    }

    
    text( x = years[nT - 10], y = c(0.9,.8,.7,.6,.5) * yMax,
          labels = c( paste( "B0 = ",  B0, sep = ""),
                      paste( "Mbar = ",  M, sep = ""),
                      paste( "qs = ",  round(qg[4],2), sep = ""),
                      paste( "qd = ",  round(qg[5],2), sep = ""),
                      "| Miss Idx"),

          cex = .9 )
} # END plotSBtIdx



# Plot fishing mortality
plotNt <- function( report = repFE, 
                    initYear = fYear, 
                    noPar = FALSE, pIdx = 1 )
{
  # Pull stuff from report
  Nat       <- report$N_atp[,,pIdx]
  Nt        <- colSums(Nat, na.rm = T)
  
  # Get number of time steps
  nT    <- report$nT
  years <- seq(from = initYear, length = nT+1, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot recruitments
  plot( x = range(years), y = c(0,max(Nt, na.rm =T) ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    mfg <- par("mfg")
    axis(side = 2, las =1 )
    if( mfg[1] == mfg[3] )
      axis(side = 1 )
    box()
    # Plot recruitment
    abline( v = vertLines, lwd = .8, lty = 2, col = "grey80")
    lines( x = years, y = Nt, lwd = 2, col = "grey40" )
    if(mfg[2] == 1)
      mtext(side = 2, text = "Numbers (1e6)", line = 3)
}


# Plot fishing mortality
plotFtg <- function(  repList = reports,  
                      noPar = FALSE, 
                      pIdx = 1,
                      labcex = 1 )
{
  report    <- repList$repOpt
  initYear  <- repList$fYear

  # Pull stuff from report
  Fgt     <- report$F_pgt[pIdx,,,drop = FALSE]
  gearIDs <- dimnames(report$F_pgt)[[2]]
  nG      <- report$nG
  nT      <- report$nT

  # Pull catch
  C_gt   <- report$totC_pgt[pIdx,,]
  C_t    <- apply( X = C_gt, FUN = sum, MARGIN = 2)

  if(!is.null(repList$posts))
  {
    tInitModel_p <- repList$repOpt$tInitModel_p
    F_igt <- repList$posts$F_ipgt[,pIdx,,,drop = FALSE]
    Fgt   <- apply( X = F_igt, FUN = mean,
                    MARGIN = c(2,3,4), na.rm = T)
  } 

  commGears <- c(1:3,6)

  minTime_g <- rep(0,nG)
  maxTime_g <- rep(nT,nG)
  for( g in commGears )
  {
    minTime_g[g] <- min(which(Fgt[1,g,] > 0),na.rm = T)
    maxTime_g[g] <- max(which(Fgt[1,g,] > 0),na.rm = T)
  }

  cols    <- brewer.pal( nG, "Dark2" )

  years <- seq(from = initYear, length.out = nT+1, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot F series
  plot( x = range(years), y = c(0,1 ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE, yaxs = "i" )
    mfg <- par("mfg")
    axis(side = 2, las =1 )
    if( mfg[1] == mfg[3] )
      axis(side = 1 )

    # Add catch axis
    maxC <- max(C_t,na.rm = T)
    maxCaxis <- 5 * maxC / 4
    CaxisLabs <- round(seq(from = 0, to = maxCaxis, length.out = 6),1)
    CaxisTicks <- seq(from = 0, to = 1, length.out = 6)
    axis( side = 4, at = CaxisTicks, labels = rev(CaxisLabs), las = 1 )
    box()
    # Plot recruitment
    grid()
    C_t <- C_t / maxC * 0.8
    rect( xleft = years - .3,
          xright = years + .3,
          ytop = 1, ybottom = 1 - C_t,
          col = "grey50", border = NA )
    for(g in commGears)
    {
      gYrs <-  max(minTime_g[g]-1,1):min(maxTime_g[g]+1,nT)
      lines( x = years[gYrs], y = Fgt[,g,gYrs], lwd = 2, col = cols[g] )
    }
    # Rescale catch so that it goes from 0 to .8?
    
    if(mfg[2] == 1)
    {
      legend( x = "topright", col = cols[c(1:3,6)],
              legend = gearIDs[c(1:3,6)], lwd = 2, bty = "n")
      mtext(side = 2, text = "Harvest\nRate", line = 3.5, cex = labcex )
    }
    if( mfg[2] == mfg[4] )
      rmtext( txt = "Catch (kt)", outer = TRUE, line = .2, cex = 1.5)
}


# Plot natural mortality
plotMt <- function( repList = reports, 
                    noPar = FALSE, 
                    pIdx = 1,
                    labcex = .8 )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  # Pull stuff from report
  Mat       <- report$M_apt[,pIdx,,drop = FALSE]
  Mbar      <- report$M
  M         <- report$M_p[pIdx]
  nT        <- report$nT
  Mjuve     <- report$Mjuve_p[pIdx]
  juveMage  <- report$juveMage


  plotCI <- FALSE

  if(!is.null(repList$posts))
  {
    plotCI <- TRUE
    tInitModel_p  <- repList$repOpt$tInitModel_p
    M_iapt        <- repList$posts$M_iapt[,,pIdx,,drop = FALSE]
    Mjuve_i       <- repList$posts$M_iapt[,juveMage,pIdx,tInitModel_p[pIdx]+1]

    M_qapt <- apply(  X = M_iapt, FUN = quantile,
                      MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                      na.rm = TRUE )

    Mat <- apply(X = M_iapt, FUN = mean, MARGIN = c(2,3,4))

    Mjuve <- mean(Mjuve_i)


  } else {
    M_qapt <- 0
  }




  years     <- seq(from = initYear, length.out = nT+1, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot F series
  plot( x = range(years), y = c(0,max(Mat,M_qapt, na.rm =T) ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    mfg <- par("mfg")
    axis(side = 2, las =1 )
    if( mfg[1] == mfg[3] )
      axis(side = 1 )
    box()
    grid()
    # Plot recruitment
    polyCol <- scales::alpha("salmon", .5)
    if( plotCI )
    {
      polygon(  x = c(years[1:nT],rev(years[1:nT])),
                y = c(M_qapt[1,juveMage+1,1,1:nT], rev(M_qapt[3,juveMage+1,1,1:nT])),
                col = polyCol, border = NA )
    }
    
    # abline( h = M, lty = 2, lwd = 2, col = "salmon")
    abline( h = Mbar, lty = 2, lwd = 2, col = "grey50")
    abline( h = Mjuve, lty = 4, lwd = 2, col = "salmon" )
    lines( x = years[1:nT], y = Mat[juveMage+1,1,1:nT], col = "salmon",
            lwd = 3  )
    # points( x = years[nT+1], y = Mat[juveMage+1,1,nT+1], col = "black",
    #         pch = 21, bg = "salmon" )
    if(mfg[2] == 1)
    {
      mtext(side = 2, text = "Natural\nMortality (/yr)", line = 3.5, cex = labcex)
      legend( x = "topleft", 
              lwd = c(2,2,2,NA),
              lty = c(1,4,2,NA),
              pch = c(22,NA,NA),
              pt.bg = c(polyCol,NA,NA), 
              pt.lwd = c(0,NA,NA),
              pt.cex = c(1.5,NA,NA),
              cex = .8,
              col = c("salmon","salmon","grey50"),
              legend = c("Mat M","Juve M","Cplx M0"), bty = "n")
    }
}


# Selectivity at age for each gear
plotSagComp <- function(  repList = list(repFE),
                          gIdx  = 1:3,
                          noPar = TRUE,
                          legend = FALSE,
                          gearLabels = NULL )
{
  SagList <- vector(mode = "list", length = length(repList) )

  nReps <- length(repList)

  # model dimensions
  nA      <- repList[[1]]$repOpt$nA
  nG      <- min(repList$repOpt$nG,length(gIdx))

  legNames <- character(length = 0)

  for( rIdx in 1:nReps )
  {
    SagList[[rIdx]] <- repList[[rIdx]]$repOpt$sel_ag[,gIdx, drop = FALSE]
    legNames <- c(legNames,repList[[rIdx]]$ctlList$ctrl$dataScenarioName)
  }

  gearLabs <- dimnames(SagList[[1]])[[2]]

  if(is.null(gearLabels))
    names(gearLabs) <- gearLabs
  else names(gearLabs) <- gearLabels

  cols <- alpha(brewer.pal( n = nReps, "Set1" ), .8)


  # Make plotting window
  if(!noPar)
    par(  mfrow = c(nG,1), 
          oma = c(3,4,.1,.1), 
          mar = c(0,1,0,1) )

  # loop over gears
  for( g in 1:nG )
  {
    plot( x = c(1,nA), y = c(0,1), type = "n",
          axes = FALSE, xlab = "", ylab = "" )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis(side = 1)
      if( mfg[2] == 1)
        axis( side = 2, las =1 )
      box()
      for( rIdx in 1:nReps)
      {
        lines( x = 1:nA, y = SagList[[rIdx]][,g], col = cols[rIdx],
                lwd = 3 )
        points( x = 1:nA, y = SagList[[rIdx]][,g], col = cols[rIdx],
                pch = 16 )
      }
      if( mfg[2] == mfg[4])
        mtext( side = 4, text = names(gearLabs)[g], line = 2 )
  }
  if( legend )
    legend( x = "bottomright", bty = "n",
            legend = legNames, col = cols, lwd = 2 )

  if( !noPar )
  {
    mtext( side = 1, text = "Age", outer = T, line = 1 )
    mtext( side = 2, text = "Selectivity", outer = T, line = 2 )
  }

}


# plotC_pgt()
# Breaks catch out into stacked bars
# by gear type. Shows SOK catch as total ponded
# fish, with a red border for dead ponded
# fish.
plotC_pgt <- function(  repList = reports,
                        yrRange = 1951:2019 )
{
  repObj    <- repList$repOpt 
  fYear     <- repList$fYear
  # Pull stuff from report
  C_pgt     <- repObj$totC_pgt
  gearLabs  <- dimnames(repObj$F_pgt)[[2]]
  nG        <- repObj$nG
  nT        <- repObj$nT
  nP        <- repObj$nP

  stockLabs <- dimnames(repObj$F_pgt)[[1]]

  gearCols  <- brewer.pal( nG, "Dark2" )

  fleetType <- repObj$fleetType_g
  postPondM <- repObj$postPondM

  gears <- which(fleetType > 0)
  nComm <- length(gears)

  years <- seq(from = fYear, length = nT, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  plotdx <- which(years %in% yrRange)

  Csum_pt <- apply( X = C_pgt, FUN = sum, MARGIN = c(1,3), na.rm = T )

  # Now loop over p
  par( mfrow = c(nP,1), mar = c(0.5,.5,0,0), oma = c(3,3,3,3))

  for( p in 1:nP )
  {
    plot( x = range(years[plotdx]), y = c(0,max(Csum_pt[p,plotdx])),
          type = "n", axes = FALSE )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      axis( side = 2, las = 1)
      box()
      grid()
      mtext( side = 4, text = stockLabs[p], font = 2 )

      for( g in gears )
      {
        if( g == 1 )
          baseC <- rep(0,nT)
        if( g > 1 )
          baseC <- apply( X = C_pgt[p,1:(g-1),,drop = FALSE], FUN = sum,
                          MARGIN = 3 )

        # if SOK, add shading
        if( fleetType[g] == 2 )
        {
          catVec <- (1 - exp(-postPondM))*C_pgt[p,g,]
        } else catVec <- C_pgt[p,g,]

        rect( xleft = years - .3, xright = years + .3,
              ybottom = baseC, ytop = baseC + catVec,
              col = gearCols[g], border = NA )

        
        

        if( p == 1 )
          legend( x = "topright", bty = "n",
                  pt.bg = c(gearCols[gears]),
                  pch = 22,
                  col = c(rep(NA,nComm)),
                  pt.lwd = c(rep(0,nComm)),
                  legend = c(gearLabs[gears]))

      }
  }

  mtext( side = 1, outer = TRUE, text = "Year", line = 2)
  mtext( side = 2, outer = TRUE, text = "Catch (kt)", line = 2 )

} # END plotC_pgt()

# Selectivity at age for each gear
plotSag <- function(  repList = reports, fleets = NULL )
{
  report <- repList$repOpt
  # pull selectivity from report
  Sag       <- report$sel_ag
  Sapgt     <- report$sel_apgt

  gearLabs <- dimnames(Sag)[[2]]

  stockLabels <- dimnames(Sapgt)[[2]]

  # model dimensions
  nA      <- report$nA
  nG      <- report$nG
  nT      <- report$nT
  nP      <- report$nP

  if(is.null(fleets))
    fleets <- 1:nG

  gearCols <- RColorBrewer::brewer.pal(nG, "Dark2")

  # Make plotting window
  par(  mfrow = c(length(fleets),1), 
        oma = c(3,4,.1,.1), 
        mar = c(0.1,1,0.1,1) )
  # loop over gears
  for( g in fleets )
  {
    plot( x = c(1,nA), y = c(0,1), type = "n",
          axes = FALSE, xlab = "", ylab = "" )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis(side = 1)
      axis( side = 2, las =1 )
      box()
      lines( x = 1:nA, y = Sag[,g], col = "grey50",
              lwd = 3 )
      panLab( x = 0.1, y = 0.9, txt = gearLabs[g], font = 2 )
      for( p in 1:nP )
        lines(  x = 1:nA, y = Sapgt[,p,g,1], 
                col = gearCols[g], lwd = 2,
                lty = p + 1 )

      if( g == 1 )
      {
        legend( x = "bottomright", bty = "n",
                legend = stockLabels, col = "black", lwd = 2,
                lty = 1:nP + 1 )
      }

  }
  
  mtext( side = 1, text = "Age", outer = T, line = 1 )
  mtext( side = 2, text = "Selectivity", outer = T, line = 2 )

}

# Vulnerable biomass
plotIdxFits <- function(  repList = reports,
                          labcex = .8 )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  # Pull vulnerable biomass and indices
  vulnBpgt      <- report$vulnB_pgt
  vulnNpgt      <- report$vulnN_pgt
  SBpt          <- report$SB_pt
  Ipgt          <- report$I_pgt
  qpg           <- report$qhat_pg
  taupg         <- report$tauObs_pg

  # Mixed indices
  mIgt          <- report$mI_gt
  qg            <- report$qhat_g
  taug          <- report$tauObs_g


  # Get survey type and whether
  # idx is relative or absolute
  surveyType  <- report$survType_g
  indexType   <- report$indexType_g
  calcIndex   <- report$calcIndex_g

  # Model dimensions
  nP          <- report$nP
  nG          <- report$nG
  nT          <- report$nT
  surveyGears <- which(calcIndex == 1)

  gearLabels  <- dimnames(vulnBpgt)[[2]]
  stockLabels <- dimnames(vulnBpgt)[[1]]

  # replace missing indices with NAs
  Ipgt[ Ipgt < 0 ] <- NA
  mIgt[ mIgt < 0 ] <- NA

  # Can we combine Ipgt and mItg into one, with 
  # an extra slice for the mixed stock if nP > 1,
  # or just overwrite the single stock gear
  # entries if nP == 1

  if( nP > 1 )
  {
    # Combine the indices
    combinedI_pgt <- array( NA, dim = c(nP+1,nG,nT))
    combinedI_pgt[1:nP,,] <- Ipgt
    combinedI_pgt[nP+1,,] <- mIgt

    Ipgt <- combinedI_pgt

    # now combine the qs
    comb_qpg <- array(NA, dim = c(nP+1,nG))
    comb_qpg[1:nP,] <- qpg
    comb_qpg[nP+1,] <- qg

    qpg <- comb_qpg

    # And the taus
    comb_taupg <- array(NA, dim = c(nP+1,nG))
    comb_taupg[1:nP,] <- taupg
    comb_taupg[nP+1,] <- taug

    taupg <- comb_taupg

    stockLabels <- c( stockLabels, "Mixed stock" )
  }

  if(nP == 1)
  {
    # Loop over gear types, merge the two together
    for( g in surveyGears )
    {
      if( all( is.na(Ipgt[,g,]) ) & any(!is.na(mIgt[g,])) )
      {
        Ipgt[,g,] <- mIgt[g,]
        qpg[,g]   <- qg[g]
        taupg[,g] <- taug[g]
      }
    }    

  }

  cols <- brewer.pal(nG,"Dark2")

  if( nP > 1 )
    nSeries <- nP + length(surveyGears) - 1
  else nSeries <- 1

  if(length(surveyGears) > 1)
    par(  mfcol = c(length(surveyGears),nSeries), 
          oma = c(3,4,1,1),
          mar = c(0,1.5,0,1.5) )
  if( length(surveyGears) == 1 )
    par(  mfrow = c(nSeries,1), 
          oma = c(3,4,1,1),
          mar = c(0,1.5,0,1.5) )

  years <- seq(initYear, by = 1, length = nT + 1 )
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  for( sIdx in 1:nSeries)
  {
    if( sIdx == nP + 1)
      pIdx <- 1:nP
    else pIdx <- sIdx
    for( g in surveyGears )
    {
      if( surveyType[g] == 1 )
      {
        vulnTS <- vulnNpgt[pIdx,g,,drop = FALSE]
        vulnTS <- apply( X = vulnTS, FUN = sum, MARGIN = 1)
        yLabel  <- "Vulnerable Numbers (1e6)"
      }
      if( surveyType[g] == 0 )
      {
        vulnTS <- vulnBpgt[pIdx,g,,drop = FALSE]
        vulnTS <- apply( X = vulnTS, FUN = sum, MARGIN = 1)
        yLabel  <- "Vulnerable Biomass (kt)"
      }
      if( surveyType[g] == 2 )
      {
        vulnTS <- SBpt[pIdx,1:nT,drop = FALSE]
        vulnTS <- apply( X = vulnTS, FUN = sum, MARGIN = 2)
        yLabel  <- "Spawning Biomass (kt)"
      }

      yRange <- range(vulnTS,Ipgt[sIdx,g,]/qpg[sIdx,g], na.rm = T )

      plot( x = range(years), y = yRange,
            xlab = "", ylab = "", type = "n",
            axes = F )
        mfg <- par("mfg")
        if(mfg[1] == mfg[3])
          axis(side = 1)
        if(mfg[1] == 1)
          mtext( side = 3, text = stockLabels[sIdx], font = 2,
                  line = 1, cex = labcex )
        axis(side = 2, las = 1)
        box()
        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = gearLabels[g], font = 2, cex = 1,
                  line = 1 )
        # Only plot lines/points if there is data
        if( all(is.na(Ipgt[sIdx,g,]) ) )
          next
        points( x = years[1:nT], y = Ipgt[sIdx,g,]/qpg[sIdx,g], 
                col = cols[g], pch = 16 )
        lines( x = years[1:nT], y = vulnTS, lwd = 2, col = "grey40" )
        panLab( x = 0.8, y = 0.8, 
                txt = paste( "q = ", signif(qpg[sIdx,g],2), sep = "") )
        panLab( x = 0.8, y = 0.7, 
                txt = paste( "tau = ", signif(taupg[sIdx,g],2), sep = "") )
        panLab( x = 0.1, y = 0.9, txt = stockLabels[sIdx] )
        
    }
  }
  mtext(  side = 1, text = "Year", outer = T, line = 2 )
  mtext(  side = 2, text = yLabel, outer = T, 
          line = 2, cex = labcex)
}

# Vulnerable biomass
plotMixedIdxFits <- function( repList = reports )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  # Pull vulnerable biomass and indices
  vulnBtgp      <- report$vulnB_tgp
  vulnNtgp      <- report$vulnN_tgp
  SBpt          <- report$SB_pt
  Itgp          <- report$I_tgp
  qpg           <- report$qhat_gp
  taugp         <- report$tauObs_gp

  # Mixed indices
  mItg          <- report$mI_tg
  qg            <- report$qhat_g
  taug          <- report$tauObs_g

  # Get survey type and whether
  # idx is relative or absolute
  surveyType  <- report$survType_g
  indexType   <- report$indexType_g
  calcIndex   <- report$calcIndex_g

  # Model dimensions
  nP          <- report$nP
  nG          <- report$nG
  nT          <- report$nT
  surveyGears <- which(calcIndex == 1)

  gearLabels  <- dimnames(vulnBtgp)[[2]]
  stockLabels <- dimnames(vulnBtgp)[[3]]

  # replace missing indices with NAs
  mItg[ mItg < 0 ] <- NA
  cols <- brewer.pal(nG,"Dark2")

  if(length(surveyGears) > 1)
    par(  mfcol = c(length(surveyGears),1), 
          oma = c(3,4,1,1),
          mar = c(0,1.5,0,1.5) )
  if( length(surveyGears) == 1 )
    par(  mfrow = c(1,1), 
          oma = c(3,4,1,1),
          mar = c(0,1.5,0,1.5) )

  years <- seq(initYear, by = 1, length = nT + 1 )
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  for( g in surveyGears )
  {
    if( surveyType[g] == 1 )
    {
      vulnTS  <- vulnNtgp[,g,1:nP,drop = FALSE]
      vulnTS  <- apply( X = vulnTS, FUN = sum, MARGIN = 1 )
      yLabel  <- "Vulnerable Numbers (1e6)"
    }
    if( surveyType[g] == 0 )
    {
      vulnTS  <- vulnBtgp[,g,1:nP,drop = FALSE]
      vulnTS  <- apply( X = vulnTS, FUN = sum, MARGIN = 1 )
      yLabel  <- "Vulnerable Biomass (kt)"
    }
    if( surveyType[g] == 2 )
    {
      vulnTS  <- SBpt[1:nP,1:nT,drop = FALSE]
      vulnTS  <- apply( X = vulnTS, FUN = sum, MARGIN = 2 )
      yLabel  <- "Spawning Biomass (kt)"
    }

    yRange <- range(vulnTS,mItg[,g]/qg[g], na.rm = T )

    plot( x = range(years), y = yRange,
          xlab = "", ylab = "", type = "n",
          axes = F )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis(side = 1)
      if(mfg[1] == 1)
        mtext( side = 3, text = "Mixed", font = 2, cex = 1,
                line = 1 )
      axis(side = 2, las = 1)
      box()
      if( mfg[2] == mfg[4] )
        mtext( side = 4, text = gearLabels[g], font = 2, cex = 1,
                line = 1 )
      # Only plot lines/points if there is data
      if( all(is.na(mItg[,g]) ) )
        next
      points( x = years[1:nT], y = mItg[,g]/qg[g], 
              col = cols[g], pch = 16 )
      lines( x = years[1:nT], y = vulnTS, lwd = 2, col = "grey40" )
      panLab( x = 0.8, y = 0.8, 
              txt = paste( "q = ", signif(qg[g],2), sep = "") )
      panLab( x = 0.8, y = 0.7, 
              txt = paste( "tau = ", signif(taug[g],2), sep = "") )
      panLab( x = 0.1, y = 0.9, txt = "Mixed stock" )
      
  }
  mtext(  side = 1, text = "Year", outer = T, line = 2 )
  mtext(  side = 2, text = yLabel, outer = T, 
          line = 2)
}

# Plot stock-recruitment curve
plotSR <- function( repList = reports )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  # Pull stuff from report
  Rpt       <- report$R_pt
  Bpt       <- report$SB_pt
  omegaRtp  <- report$omegaR_tp
  sigmaR    <- report$sigmaR
  R0_p      <- report$R0_p
  B0_p      <- report$B0_p
  reca_p    <- report$reca_p
  recb_p    <- report$recb_p
  h_p       <- round(report$rSteepness_p,2)
  h         <- round(report$rSteepness,2)
  phi_p     <- report$phi_p
  
  # Get model dimensions
  nT      <- report$nT
  nP      <- report$nP

  # Get max depletion level
  D_pt    <- Bpt
  for( p in 1:nP )
    D_pt[p,] <- D_pt[p,] / B0_p[p]

  # stock names
  stockNames <- dimnames(Rpt)[[1]]

  cols    <- brewer.pal(n = max(nP,3), "Set1") 

  SB <- seq(0,1.2*max(Bpt,B0_p),length = 1000 )
  # Recruitment
  R_p <-  matrix(0, nrow = nP, ncol = length(SB) )
  for( p in 1:nP) 
    R_p[p,] <- reca_p[p] * SB / (1 + recb_p[p]*SB)
  # 20% level
  B20_p <- 0.2*B0_p
  R20_p <- reca_p * B20_p / (1 + recb_p*B20_p)

  par( mfrow = c(nP,1), mar = c(1,1,1.5,1), oma = c(3,3,1,1))

  for( p in 1:nP )
  {
    plot( x = c(0,max(Bpt[p,])), y = c(0,1.2 * max(R_p[p,],Rpt[p,],na.rm = T ) ),
          type = "n", las = 1, xlab = "",
          ylab = "")
      mfg <- par("mfg")
      lines( x = SB, y = R_p[p,], lwd = 3, col = cols[p] )
      points( x = Bpt[p,1:nT], y = Rpt[p,2:(nT+1)], pch = 16, col = cols[p] )
      mtext( side = 3, text = stockNames[p], font = 2)
      # Plot B0,R0
      segments( x0 = B0_p[p], y0 = 0, y1 = R0_p[p],
                lty = 2, lwd = 2, col = cols[p] )
      segments( x0 = 0, x1 = B0_p[p], y0 = R0_p[p], y1 = R0_p[p],
                lty = 2, lwd = 2, col = cols[p] )
      # Plot B20,R20
      segments( x0 = B20_p[p], x1 = B20_p[p], y0 = 0, y1 = R20_p[p],
                lty = 2, lwd = 2, col = cols[p] )
      segments( x0 = 0, x1 = B20_p[p], y0 = R20_p[p], y1 = R20_p[p],
                lty = 2, lwd = 2, col = cols[p] )
      # Label with steepness
      text( x = 1.1, y = R0_p[p]*1.3, labels = paste("h = ", h_p[p], sep = "") )
  }
  


  mtext( side = 1, text = "Spawning Stock Biomass", line = 1.5, outer = TRUE )
  mtext( side = 2, text = "Age-1 Recruits (1e6)", line = 2, outer =  TRUE )

}

# recruitments
plotRt <- function( repList = reports,
                    noPar = FALSE, 
                    pIdx = 1, 
                    labcex = .8,
                    plotLog = TRUE,
                    heading = "dimLabs" )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  # Pull stuff from report
  Rt      <- report$R_pt[pIdx,]
  omegaRt <- report$omegaR_pt[pIdx,]
  sigmaR  <- report$sigmaR
  R0      <- report$R0_p[pIdx]

  ageData_agt <- repList$data$A_apgt[,pIdx,,]

  colRec_t <- rep("grey40", length(Rt))

  checkPos <- function( x )
  {
    if(any(x > 0))
      return(1)
    else return(0)
  }

  posAgeIdx_t <- apply( X = ageData_agt, FUN = checkPos, MARGIN = 3 )

  colRec_t[posAgeIdx_t == 0] <- "white"
  


  plotCI <- FALSE

  if(!is.null(repList$posts))
  {
    plotCI <- TRUE
    tInitModel_p <- repList$repOpt$tInitModel_p
    R_it <- repList$posts$R_ipt[,pIdx,]

    R_qt <- apply(X = R_it, FUN = quantile,
                    MARGIN = c(2), probs = c(0.025, 0.5, 0.975),
                    na.rm = TRUE )

    Rt <- apply(X = R_it, FUN = mean, MARGIN = c(2))

    R0 <- round(mean(repList$posts$R0_ip[,pIdx]),2)

  } else {
    R_qt <- 0
  }

  R0lab <- round(R0,2)

  Rt[Rt == 0] <- NA
  R_qt[R_qt == 0] <- NA

  if(plotLog)
  {
    R0 <- log(R0, base = 10)
    R_qt <- log(R_qt, base = 10)
    Rt <- log(Rt, base = 10)

  }

  # stock labels
  stockLabs <- dimnames(report$R_pt)[[1]]


  if( !is.null(heading))
  {
    if( heading == "dimLabs")
      stockLabel <- stockLabs
    else
      stockLabel <- try(heading[stockLabs])

    if( class(stockLabel) == "try-error" ) 
      stockLabel <- stockLabs
  }

  # Get number of time steps
  nT    <- report$nT
  years <- seq(from = initYear, length = nT+1, by = 1)

  vertLines <- seq(from = 1950, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot recruitments
  plot( x = range(years), y = c(0,max(Rt,R_qt, na.rm =T) ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    mfg <- par("mfg")
    if(plotLog)
    {
      maxY <- ceiling(max(Rt,R_qt, na.rm =T))
      yTicks <- 0:maxY
      yLabs <- 10^yTicks
    }
    if( mfg[1] == mfg[3])
      axis(side = 1)
    if(!plotLog)
      axis(side = 2, las = 1  )
    if(plotLog)
      axis(side = 2, las = 1, at = yTicks, labels = yLabs )
    box()
    grid()
    # Plot recruitment
    # lines( x = years[1:nT], y = Rt[1:nT], lwd = 2, col = "grey40" )
    abline( h = R0, lty = 2, lwd = 2)
    if( plotCI )
      segments( x0 = years[1:nT],
                y0 = R_qt[1,1:nT],
                y1 = R_qt[3,1:nT],
                lwd = 1, col = "grey40" )
    points( x = years[1:nT], y = Rt[1:nT], lwd = 2, bg = "grey40",
            pch = 21, col = "grey40" )

    axis( side = 3, at = years[posAgeIdx_t == 0], labels = FALSE,
          col.ticks = "grey40", tck = .02, lwd.ticks = 3 )
    
    panLab( x = 0.8, y = 0.1, txt = paste("R0 = ", R0lab, sep = "") )
    panLab( x = 0.8, y = 0.2, txt = "|  Miss Ages")
    if( mfg[1] == mfg[3])
      mtext( side = 1, text = "Year", line = 2)
    if( mfg[2] == 1 )
      mtext(side = 2, text = "Recruits\n(millions)", line = 3.5, cex = labcex)
    if( mfg[1] == 1 & !is.null(heading) )
      mtext(side = 3, text = stockLabel[pIdx], line = 1, font = 2, cex = labcex )

}

# recruitments, adjusted for brood year
plotRtResids <- function( repList = reports,
                          noPar = FALSE, 
                          pIdx = 1, labcex = .8 )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  # Pull stuff from report
  Rt      <- report$R_pt[pIdx,]
  omegaRt <- report$omegaR_pt[pIdx,]
  sigmaR  <- report$sigmaR
  R0      <- report$R0_p[pIdx]

  omegaRt <- omegaRt / sigmaR

  # Get number of time steps
  nT    <- report$nT
  years <- seq(from = initYear, length = nT+1, by = 1)

  # Adjust by brood year
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot recruitment resids
  plot( x = range(years), y = range(omegaRt,omegaRt-sigmaR^2/2, na.rm =T),
        type = "n", xlab = "", ylab = "", axes = FALSE,
        las = 1 )
    # Add axis label if on the left margin
    mfg <- par("mfg")
    if( mfg[2] == 1 )
      mtext(side = 2, text = "Recruitment std.\nlog-residuals", line = 3, cex = labcex)  
    if( mfg[1] == mfg[3] )
    {
      axis( side = 1 )
      mtext( side = 1, text = "Year", line = 2 )
    }
    # Add y axis
    axis( side = 2, las = 1 )
    box()
    # Plot recruitment - update to include SEs later
    abline( h = 0, lty = 2, lwd = 1)
    grid()
    abline( h = mean(omegaRt,na.rm =T), lty = 3, lwd = 2, col = "red")
    points( x = years[1:(nT)], y = omegaRt[1:nT], col = "grey40", pch = 16 )
    panLab( x = 0.8, y = 0.8, txt = paste("sigmaR = ", round(sigmaR,2), sep = "" ) )

    
    panLegend(  x = 0.1, y = 0.95, 
                bty = "n",
                legTxt = c("Mean Resid"),
                lty = c(3),
                lwd = c(2),
                col = c("red") )
}

# # Plot age fits
# plotAgeFitYrs <- function(  report = repFE,
#                             initYear = fYear,
#                             gearLabels = gearLabs,
#                             pIdx = 1 )
# {
#   # Pull predicted and observed ages
#   predAge <- report$predPA_apgt[,,,pIdx]
#   obsAge  <- report$A_apgt[,,,pIdx]

#   # Pull model dims
#   nG      <- report$nG
#   nT      <- report$nT
#   nA      <- report$nA

  

#   # Make years vector
#   years   <- seq(initYear, length = nT+1, by = 1)

#   # Now, we want to loop over gear types now, 
#   # and record the gIdxes for which there are
#   # observations
#   ageGears <- c()
#   gearTimes <- vector(mode = "list", length = nG)
#   for( gIdx in 1:nG )
#   {
#     if( any(obsAge[1,,gIdx] >= 0) )
#     {
#       ageGears <- c(ageGears,gIdx)
#       gearTimes[[gIdx]] <- which(obsAge[1,,gIdx] >= 0)
#     }
#   }

#   # Make colours vector
#   cols    <- brewer.pal( n = length(ageGears), "Dark2" )

#   # ok, ageGears are the ones we want to plot,
#   # and gearTimes is the list of time indices
#   for( aIdx in 1:length(ageGears) )
#   {
#     gIdx <- ageGears[aIdx]
#     dev.new()
#     times <- gearTimes[[gIdx]]
#     # Count the number of age observations
#     # there are, and make the plotting window
#     nObs <- length(times)
#     nCols <- round(sqrt(nObs))
#     nRows <- ceiling(nObs/nCols)

#     par(  mfcol = c(nRows,nCols), 
#           mar = c(1,1,1,1),
#           oma = c(3,3,3,3) )

#     ageObs <- obsAge[,,gIdx]
#     # ageObs <- sweep( x = ageObs, FUN = "/", MARGIN = c(2), STATS = "sum")

#     for( tIdx in times )
#     { 
#       # get age obs and preds
#       ageObs.t  <- ageObs[,tIdx]/sum(ageObs[,tIdx])
#       agePred <- predAge[,tIdx,gIdx]

#       plot( x = c(1,nA), y = c(0,max(ageObs.t,agePred,na.rm = T) ),
#             xlab = "", ylab = "", type = "n", las = 1 )
#         rect( xleft = 1:nA - .3, xright = 1:nA + .3,
#               ybottom = 0, ytop = ageObs.t,
#               col = "grey40", border = NA )
#         lines(  x = 1:nA, y = agePred, lwd = 1,
#                 col = cols[aIdx] )
#         points(  x = 1:nA, y = agePred,
#                 col = cols[aIdx], pch = 21 )
#         panLab( x=.5, y = .95, txt = years[tIdx] )

#     }
#     mtext( side = 1, outer = T, text = "Age", line = 2 )
#     mtext( side = 2, outer = T, text = "Proportion", line = 2 )
#   }

# }

# Plot age fits averaged over time and stock
plotAgeFitAggAvg <- function( repList = reports )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  # Pull predicted and observed ages
  predAge <- report$predPA_apgt
  obsAge  <- report$A_apgt

  minAge <- repList$data$minAge

  gearLabs  <- dimnames(obsAge)[[3]]

  # replace missing entries with NAs
  predAge[predAge < 0] <- NA
  obsAge[obsAge < 0] <- NA

  # Get time steps
  nT <- dim(predAge)[4]
  nG <- dim(predAge)[3]
  nP <- dim(predAge)[2]

  # Multiply the yearly predicted proportions by the
  # yearly sample sizes
  for( g in 1:nG )
    for( p in 1:nP )
      for( t in 1:nT )
      {
        thisSamp <- sum( obsAge[,p,g,t], na.rm = T )
        if( thisSamp < 0 )
          thisSamp <- 0

        predAge[,p,g,t] <- thisSamp * predAge[,p,g,t]
      }

  
  
  # Average over time
  predAge <- apply( X = predAge, FUN = sum, MARGIN = c(1,3), na.rm = T )
  obsAge <- apply( X = obsAge, FUN = sum, MARGIN = c(1,3), na.rm = T )    

  
  # Pull model dims
  nG      <- report$nG
  nA      <- report$nA
  nT      <- report$nT
  nP      <- report$nP
  minPA   <- report$minPropAge


  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  ageGears <- c()
  gearTimes <- vector(mode = "list", length = nG)
  for( gIdx in 1:nG )
    if( sum(obsAge[,gIdx]) > 0  )
      ageGears <- c(ageGears,gIdx)

  # Make colours vector
  cols    <- brewer.pal( n = nG, "Dark2" )


  par(  mfcol = c(length(ageGears),1), 
        mar = c(1,2,1,2),
        oma = c(3,3,3,3) )

  # ok, ageGears are the ones we want to plot,
  for( aIdx in 1:length(ageGears) )
  { 
    gIdx <- ageGears[aIdx]
    # get age obs and preds
    ageObs  <- obsAge[,gIdx]
    nSamp   <- sum(ageObs)
    ageObs  <- ageObs / sum(ageObs)
    agePred <- predAge[,gIdx]
    agePred <- agePred / sum( agePred )

    plot( x = c(1,nA), y = c(0,1.1*max(ageObs,agePred,na.rm = T) ),
          xlab = "", ylab = "", type = "n", las = 1 )
      legend( x = "topright",
              bty = "n",
              legend = paste("N = ", nSamp, sep = "" ) )
      rect( xleft = minAge:nA - .3, xright = minAge:nA + .3,
            ybottom = 0, ytop = ageObs[minAge:nA],
            col = "grey40", border = NA )
      lines(  x = minAge:nA, y = agePred[minAge:nA], lwd = 3,
              col = cols[gIdx] )
      points(  x = minAge:nA, y = agePred[minAge:nA],
              col = cols[gIdx], pch = 16, cex = 1.5 )
      abline( h = minPA, lty = 2, lwd = .8 )
      
      
      # panLab( x=.5, y = .95, txt = gearLabels[gIdx] )
      mfg <- par("mfg")
      if(mfg[1] == 1 )
        mtext( side = 3, text = "Aggregate SAR", line = 1, font = 2)

      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text( x = corners[2]+0.5, 
              y = mean(corners[3:4]), 
              labels = gearLabs[gIdx], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
  }
  
  mtext( side = 1, outer = T, text = "Age", line = 2 )
  mtext( side = 2, outer = T, text = "Proportion", line = 2 )

}

# Plot age fits averaged over time
plotAgeFitAvg <- function(  repList = reports )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs

  # Pull predicted and observed ages
  predAge <- report$predPA_apgt
  obsAge  <- report$A_apgt

  minAge <- repList$data$minAge

  gearLabs  <- dimnames(obsAge)[[3]]
  stockLabs <- dimnames(obsAge)[[2]]


  # replace missing entries with NAs
  predAge[predAge < 0] <- NA
  obsAge[obsAge < 0] <- NA


  # Get time steps
  nT <- dim(predAge)[4]
  nG <- dim(predAge)[3]
  nP <- dim(predAge)[2]

  # Multiply the yearly predicted proportions by the
  # yearly sample sizes
  for( g in 1:nG )
    for( p in 1:nP )
      for( t in 1:nT )
      {
        thisSamp <- sum( obsAge[,p,g,t], na.rm = T )
        if( thisSamp < 0 )
          thisSamp <- 0

        predAge[,p,g,t] <- thisSamp * predAge[,p,g,t]
      }


  
  # Average over time
  predAge <- apply( X = predAge, FUN = sum, MARGIN = c(1,2,3), na.rm = T )
  obsAge <- apply( X = obsAge, FUN = sum, MARGIN = c(1,2,3), na.rm = T )    

  # Pull model dims
  nG      <- report$nG
  nA      <- report$nA
  nT      <- report$nT
  nP      <- report$nP
  minPA   <- report$minPropAge

  # Make years vector
  years   <- seq(initYear, length = nT+1, by = 1)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  ageGears <- c()
  gearTimes <- vector(mode = "list", length = nG)
  for( gIdx in 1:nG )
    if( sum(obsAge[,,gIdx]) > 0  )
      ageGears <- c(ageGears,gIdx)

  # Make colours vector
  cols    <- brewer.pal( n = nG, "Dark2" )


  par(  mfcol = c(length(ageGears),nP), 
        mar = c(1,2,1,2),
        oma = c(3,3,3,3) )

  # ok, ageGears are the ones we want to plot,
  for( pIdx in 1:nP)
    for( aIdx in 1:length(ageGears) )
    { 
      gIdx <- ageGears[aIdx]
      # get age obs and preds
      ageObs  <- obsAge[,pIdx,gIdx]
      nSamp   <- sum(ageObs)
      ageObs  <- ageObs / sum(ageObs)
      agePred <- predAge[,pIdx,gIdx]
      agePred <- agePred / sum( agePred )

      plot( x = c(1,nA), y = c(0,1.1*max(ageObs,agePred,na.rm = T) ),
            xlab = "", ylab = "", type = "n", las = 1 )
        legend( x = "topright",
                bty = "n",
                legend = paste("N = ", nSamp, sep = "" ) )
        rect( xleft = minAge:nA - .3, xright = minAge:nA + .3,
              ybottom = 0, ytop = ageObs[minAge:nA],
              col = "grey40", border = NA )
        lines(  x = minAge:nA, y = agePred[minAge:nA], lwd = 3,
                col = cols[gIdx] )
        points(  x = minAge:nA, y = agePred[minAge:nA],
                col = cols[gIdx], pch = 16, cex = 1.5 )
        abline( h = minPA, lty = 2, lwd = .8 )
        
        
        # panLab( x=.5, y = .95, txt = gearLabels[gIdx] )
        mfg <- par("mfg")
        if(mfg[1] == 1 )
          mtext( side = 3, text = stockLabs[pIdx], line = 1, font = 2)

        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text( x = corners[2]+0.5, 
                y = mean(corners[3:4]), 
                labels = gearLabs[gIdx], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }
    }
  
  mtext( side = 1, outer = T, text = "Age", line = 2 )
  mtext( side = 2, outer = T, text = "Proportion", line = 2 )

}

# Plot age fits averaged over time
plotMixedAgeFitAvg <- function( repList = reports )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs

  # Pull predicted and observed ages
  predAge <- report$predMixedPA_atg
  obsAge  <- report$mA_atg

  gearLabs  <- dimnames(obsAge)[[3]]

  # replace missing entries with NAs
  predAge[predAge < 0] <- NA
  obsAge[obsAge < 0] <- NA
  
  # Average over time
  predAge <- apply( X = predAge, FUN = mean, MARGIN = c(1,3), na.rm = T )
  obsAge <- apply( X = obsAge, FUN = mean, MARGIN = c(1,3), na.rm = T )    

  # Pull model dims
  nG      <- report$nG
  nA      <- report$nA
  nT      <- report$nT
  minPA   <- report$minPropAge

  # Make years vector
  years   <- seq(initYear, length = nT+1, by = 1)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  ageGears <- c()
  gearTimes <- vector(mode = "list", length = nG)
  for( gIdx in 1:nG )
    if( any(!is.na(obsAge[1,gIdx]) ) )
      ageGears <- c(ageGears,gIdx)

  # Make colours vector
  cols    <- brewer.pal( n = nG, "Dark2" )


  par(  mfcol = c(length(ageGears),1), 
        mar = c(1,1,1,1),
        oma = c(3,3,3,3) )

  # ok, ageGears are the ones we want to plot,
  for( aIdx in 1:length(ageGears) )
  { 
    gIdx <- ageGears[aIdx]
    # get age obs and preds
    ageObs  <- obsAge[,gIdx]
    ageObs  <- ageObs / sum(ageObs)
    agePred <- predAge[,gIdx]

    plot( x = c(1,nA), y = c(0,max(ageObs,agePred,na.rm = T) ),
          xlab = "", ylab = "", type = "n", las = 1 )
      rect( xleft = 1:nA - .3, xright = 1:nA + .3,
            ybottom = 0, ytop = ageObs,
            col = "grey40", border = NA )
      lines(  x = 1:nA, y = agePred, lwd = 3,
              col = cols[aIdx] )
      points(  x = 1:nA, y = agePred,
              col = cols[aIdx], pch = 16, cex = 1.5 )
      abline( h = minPA, lty = 2, lwd = .8 )
      panLab( x=.5, y = .95, txt = gearLabels[gIdx] )
      mfg <- par("mfg")
      if(mfg[1] == 1)
        mtext( side = 3, text = "Mixed stock", line = 1, font = 2)
  }
  
  mtext( side = 1, outer = T, text = "Age", line = 2 )
  mtext( side = 2, outer = T, text = "Proportion", line = 2 )

}


# Plot Biomass, catch and indices on the same
# axes
plotItResids <- function( repList =reports,
                          noPar = FALSE,
                          pIdx = 1, labcex = .8 )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs
  # Pull stuff from report
  Bt     <- report$SB_pt[pIdx,]
  Igt    <- report$I_pgt[pIdx,,]
  qg     <- report$qhat_pg[pIdx,]
  rIgt   <- repList$data$rI_pgt[pIdx,,]

  # Need to pull vuln biomass/numbers
  vulnBgt <- report$vulnB_pgt[pIdx,,]
  vulnNgt <- report$vulnN_pgt[pIdx,,]

  # Pull probability of positive indices
  probPosIdx_gt <- report$probPosIdx_pgt[pIdx,,]

  # Get survey type and whether
  # idx is relative or absolute
  surveyType  <- report$survType_g      # abs or relative?
  indexType   <- report$indexType_g     # biomass or numbers
  calcIndex   <- report$calcIndex_g     # calc index at all?

  # Model dimensions
  nG          <- report$nG
  nT          <- report$nT
  surveyGears <- which(calcIndex == 1)

  # Get survey names
  surveyNames <- dimnames(vulnBgt)[[1]][surveyGears]
  stockLabel  <- dimnames(report$I_pgt)[[1]]

  if( stockLabel == "Aggregate" )
  {
    # We need to pull the aggregated stock
    # indices in at this point
    mIgt    <- report$mI_gt
    mqg     <- report$qhat_g

    mIgt[mIgt < 0] <- 0

    for( g in surveyGears )
    {
      if( all(Igt[g,] < 0) & any(mIgt[g,] > 0))
      {
        Igt[g,] <- mIgt[g,]
        qg[g]   <- mqg[g]
      }
    }
  }

  Igt[Igt<=0] <- NA

  if(!is.null(repList$posts))
  {
    tInitModel_p <- repList$repOpt$tInitModel_p
    SB_it <- repList$posts$SB_ipt[,pIdx,]

    Bt <- apply(X = SB_it, FUN = mean, MARGIN = c(2))

    q_ig  <- repList$posts$q_ipg[,pIdx,]
    qg   <- apply( X = q_ig, FUN = mean, MARGIN = 2)

    probPosIdx_igt <- repList$posts$probPosIdx_ipgt[,pIdx,,]
    probPosIdx_gt  <- apply( X = probPosIdx_igt, FUN = mean, MARGIN = c(2,3))

  } 


  # Scale indices by q
  IgtScaled <- Igt
  for( g in 1:nG)
  {
    IgtScaled[g,] <- Igt[g,] / qg[g] / probPosIdx_gt[g,]
  }

  IgtScaled[IgtScaled < 0] <- NA

  cols <- brewer.pal(nG, name = "Dark2")

  # Get number of time steps
  years <- seq(from = initYear, length = nT, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  tauObs_g <- report$tauObs_pg[pIdx,]

  resids <- IgtScaled

  # Calc residuals
  for( g in 1:nG )
  {
    posIdx <- which(IgtScaled[g,] > 0)

    # Check index type, calcs resids
    if( surveyType[g] == 0 )
      resids[g,posIdx] <- log(IgtScaled[g,posIdx] / vulnBgt[g,posIdx] ) / tauObs_g[g]
    if( surveyType[g] == 1 )
      resids[g,posIdx] <- log(IgtScaled[g,posIdx] / vulnNgt[g,posIdx] ) / tauObs_g[g]
    if( surveyType[g] == 2 )
      resids[g,posIdx] <- log(IgtScaled[g,posIdx] / Bt[posIdx] ) / tauObs_g[g]

  }   

  combResids <- rep(NA,nT)
  if( any(repList$data$whichCombIdx_g == 1) )
    combResids <- -1 * report$zComb_pt[pIdx,]


  # Now set up plotting window
  if( !noPar )
    par( mfrow = c(1, 1), mar = c(2,3,1,3), oma = c(3,3,1,1) )

  plot( x = range(years), y = range(resids,combResids,na.rm=T), xlab = "",
          ylab = "", type = "n", las = 1 )
    mfg <- par("mfg")
    # show 0 line
    abline( h = 0, lty = 3, lwd = .8 )
    abline( v = vertLines, lwd = .8, lty = 2, col = "grey80")
    # Now show the resids
    for(g in surveyGears)
    {
      residTimes <- which(!is.na(resids[g,]))
      if(length(residTimes) > 0)
      {
        points( x = years[residTimes], y = resids[g,residTimes], pch = 16,
                col = cols[g] )
        # Fit a regression
        dat <- data.frame(x = years[residTimes], y = resids[g,residTimes])
        residModel <- lm( y~x, data = dat )
        dat$pred <- predict.lm(object = residModel, newdata = dat)

        pVal <- round(summary(residModel)$coefficients[2,4],2)

        lines( x = dat$x, y = dat$pred,
                col = cols[g], lwd = 2 )
        text( x = dat$x[4], y = min(dat$y) + .3, label = paste("p = ", pVal, sep = ""),
              col = cols[g], font = 2 )

        meanResid <- mean(resids[g,residTimes])
        segments(  x0 = years[residTimes[1]],
                  x1 = years[residTimes[length(residTimes)]], 
                  y0 = meanResid, col = cols[g], lty = 2, lwd = 2 ) 
      }
    }
    combResids[combResids == 0] <- NA
    combResidTimes <- which(!is.na(combResids))

    if(length(combResidTimes) > 0)
    {
      points( x = years[combResidTimes], y = combResids[combResidTimes], pch = 16,
                col = "grey40" )
        # Fit a regression
        dat <- data.frame(x = years[combResidTimes], y = combResids[combResidTimes])
        residModel <- lm( y~x, data = dat )
        dat$pred <- predict.lm(object = residModel, newdata = dat)

        pVal <- round(summary(residModel)$coefficients[2,4],2)

        lines( x = dat$x, y = dat$pred,
                col = "grey40", lwd = 2 )
        text( x = dat$x[7], y = min(dat$y + .3), label = paste("p = ", pVal, sep = ""),
              col = "grey40", font = 2 )

        meanResid <- mean(combResids[combResidTimes])
        segments(  x0 = years[combResidTimes[1]],
                  x1 = years[combResidTimes[length(combResidTimes)]], 
                  y0 = meanResid, col = cols[g], lty = 2, lwd = 2 )

        surveyNames <- c(surveyNames,"Spawn Index")
    }
    if(mfg[2] == 1 )
      legend( x = "topleft",
              legend = surveyNames,
              pch = 16,
              col = c(cols[which(calcIndex == 1)],"grey40") )
    if( mfg[2] == 1)
      mtext(  side = 2, outer = F,
              text = paste("Std. log-residuals", sep = ""),
              line = 3, cex = labcex )
}

# Plot standardised residuals for each gear
plotCompResids <- function( repList = reports,
                            comps = "age",
                            pIdx = 1 )
{
  repObj    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs

  # Pull age data and predictions
  nA            <- repObj$nA
  resids_agt    <- repObj$ageResids_apgt[,pIdx,,]
  yLab          <- "Age"
  tau2Age_g     <- repObj$tau2Age_pg[pIdx,]

  # Take sqrt for stdizing reisds
  tauAge_g      <- sqrt(tau2Age_g)

  
  # Pull gear names
  gearLabs  <- dimnames( repObj$predPA_apgt )[[3]]
  stockName <- dimnames( repObj$predPA_apgt )[[2]][pIdx]
  nG        <- dim(resids_agt)[2]
  nT        <- dim(resids_agt)[3]

  yrs <- seq( initYear, by = 1, length.out = nT )

  ageGears  <- c()

  for( g in 1:nG )
    if( any( resids_agt[,g,] != 0 ) )
      ageGears <- c(ageGears,g)

  # Count age gears
  nAgeGears <- length(ageGears)

  par(  mfrow = c(nAgeGears,1), 
        mar = c(1,.5,1,.5),
        oma = c(3,3,3,2) )

  resids_agt <- resids_agt / max(abs(range(resids_agt)))

  for( g in ageGears )
  {
    # Plot window
    plot( x = range(yrs), y = c(0,nA+1),
          axes = FALSE, type = "n", 
          xlab = "", ylab = "" )
      axis( side = 2, las = 1, at = 1:nA )
      grid()
      box()
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      
      # Gear label
      mtext( side = 3, text = paste(stockName," ", gearLabs[g], sep = ""), 
              font = 2, line = 0.5 )

      # Loop over years, plot circles
      for( t in 1:nT )
      {
        if( all( resids_agt[,g,t] == 0 ) )
          next

        res_a <- resids_agt[,g,t]/tauAge_g[g]

        cols  <- rep( "black", nA )
        cols[res_a < 0] <- "red"
        cols[res_a == 0 ] <- NA

        radians <- seq(0,2*pi, length.out = 100)
        for( a in 1:nA)
          if( res_a[a] != 0 )
            lines(  x = yrs[t] + abs(res_a[a]) * cos(radians),
                    y = a + abs(res_a[a]) * sin(radians),
                    col = cols[a], lwd = 1 )


      }
      legend( x = "topright", bty = "n",
              legend = paste("tau = ", round(tauAge_g[g],2), sep = "") )


    
  }
  mtext( side = 2, text = "Age", outer = T, line = 2)
  mtext( side = 1, text = "Year", outer = T, line = 2)

} #END plotCompResids()

# Plot comp fits
plotCompFitYrs <- function( repList = reports,
                            comps = "age",
                            save = FALSE,
                            savePath = "plotFitYrs",
                            pIdx = 1,
                            tc = FALSE )
{
  repObj    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs



  # Pull predicted and observed ages
  if( comps == "age" )
  {
    max       <- repObj$nA
    pred_xgt  <- repObj$predPA_apgt[1:max,pIdx,,]
    obs_xgt   <- repObj$A_apgt[1:max,pIdx,,]  
    xLab      <- "Age"

    minX      <- repList$data$minAge

    # browser()

    if( tc )
    {
      pred_xgt  <- repObj$tcPred_apgt[1:max,pIdx,,]
      obs_xgt   <- repObj$tcComps_apgt[1:max,pIdx,,]
    }
  }
  
  dimNames    <- dimnames(pred_xgt)
  compNames   <- dimNames[[1]]
  yearNames   <- dimNames[[3]]
  gearNames   <- dimNames[[2]]
  stockName   <- dimnames(repObj$A_apgt)[[2]][pIdx]

  # Pull model dims
  nG      <- repObj$nG
  nT      <- repObj$nT  
  minPA   <- repObj$minPropAge

  # Make colours vector
  cols    <- brewer.pal( n = nG, "Dark2" )

  # Make years vector
  years   <- seq(initYear, length = nT+1, by = 1)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  obsGears <- c()
  gearTimes <- vector(mode = "list", length = nG)
  for( gIdx in 1:nG )
  {
    if( any(obs_xgt[1,gIdx,] >= 0) )
    {
      obsGears <- c(obsGears,gIdx)
      gearTimes[[gIdx]] <- which(obs_xgt[1,gIdx,] >= 0)
    }
  }

  # ok, obsGears are the ones we want to plot,
  # and gearTimes is the list of time indices
  for( gIdx in obsGears )
  {
    if(save)
    {
      gearPath <- paste(savePath,gearNames[gIdx],".png",sep = "")
      png(  gearPath, 
            width = 11, height = 8.5,
            units = "in", res = 300 )
    }

    times <- gearTimes[[gIdx]]
    # Count the number of age observations
    # there are, and make the plotting window
    nObs <- length(times)
    nCols <- round(sqrt(nObs))
    nRows <- ceiling(nObs/nCols)

    par(  mfcol = c(nRows,nCols), 
          mar = c(0,0.5,0,0.5),
          oma = c(3,3,3,3) )

    yMax <- max(pred_xgt, na.rm = T ) 

    for( tIdx in times )
    { 
      Nobs <- sum(obs_xgt[,gIdx,tIdx])
      # get age obs and preds
      compObsProp_x  <- obs_xgt[,gIdx,tIdx]/sum(obs_xgt[,gIdx,tIdx])
      compPred_x     <- pred_xgt[,gIdx,tIdx]



      plot( x = c(1,max), y = c(0,yMax ),
            xlab = "", ylab = "", type = "n", las = 1,
            axes = FALSE )
        box()
        grid()
        mfg <- par("mfg")
        if( mfg[1] == mfg[3] )
          axis( side = 1 )

        if( mfg[2] == 1 )
          axis( side = 2, las = 1 )

        rect( xleft = minX:max - .3, xright = minX:max + .3,
              ybottom = 0, ytop = compObsProp_x[minX:max],
              col = "grey40", border = NA )
        lines(  x = minX:max, y = compPred_x[minX:max], lwd = 1,
                col = cols[gIdx] )
        points(  x = minX:max, y = compPred_x[minX:max],
                col = cols[gIdx], pch = 16 )
        abline( h = minPA, lty = 2, lwd = .8 )
        panLab( x=.5, y = .95, txt = years[tIdx], font = 2 )
        legend( x = "topright",
                legend = paste("N = ", Nobs, sep = "" ),
                bty = "n" )

    }
    mtext( side = 3, outer = T, text = paste(stockName, " ", gearNames[gIdx]), line = 2, font = 2)
    mtext(  side = 1, outer = T, text = xLab, line = 2 )
    mtext(  side = 2, outer = T, text = paste("Proportion-at-",xLab,sep=""), 
            line = 2 )
    if(save)
      dev.off()
  }

} # END plotCompFitYrs()

plotCompLikelihoodError <- function(  repList = reports)
{
  repObj    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs
  # Pull tail compressed compositions and
  # predictions
  tcComps_apgt  <- repObj$tcComps_apgt
  tcPred_apgt   <- repObj$tcPred_apgt 

  nP <- repObj$nP
  nG <- repObj$nG
  nT <- repObj$nT

  stockLabs <- dimnames(tcComps_apgt)[[2]]
  gearLabs <- dimnames(tcComps_apgt)[[3]]

  # Pull annual weights and geometric means
  gmObs_pgt   <- repObj$gmObs_pgt
  gmPred_pgt  <- repObj$gmPred_pgt
  ageWt_pgt   <- repObj$ageWt_pgt

  # Pull correlation matrices
  Corr_gaa    <- repObj$Corr_gaa

  # Pull internally calculated component
  intrnlAgeLikCpt_pg <- repObj$intrnlAgeLikCpt_pg

  # Pull etaSumSq
  modelEtaSumSq_pg  <- repObj$etaSumSq_pg

  # Pull model tau2Age
  tau2Age_pg        <- repObj$tau2Age_pg

  # And full likelihood function value
  ageObsNLL_pg <- repObj$ageObsNLL_pg

  # Calculate mean sample size
  meanA_apgt <- repObj$tcComps_apgt
  meanA_apgt[meanA_apgt < 0] <- NA
  meanA_pgt <- apply( X = meanA_apgt, FUN = sum, MARGIN = c(2,3,4), na.rm = T)
  meanA_pgt[meanA_pgt == 0] <- NA
  meanA_pg <- apply( X = meanA_pgt, FUN = mean, MARGIN = c(1,2), na.rm = T )
  meanA_pg[is.nan(meanA_pg)] <- 0

  # Now, we need to loop over stock, gear and time,
  # plot the errors in the annual weights, 
  # and annotate with the difference between
  # the etaSumSq values and the internally
  # calculated likelihood component.
  wtErr_pgt     <- array( NA, dim = c(nP, nG, nT) )  # On the log scale
  nBins_pgt     <- array( 0,  dim = c(nP, nG, nT) )  # number of bins
  gmObsErr_pgt  <- array( NA, dim = c(nP, nG, nT) )  # error in geometric mean calc for obs
  gmPredErr_pgt <- array( NA, dim = c(nP, nG, nT) )  # error in geometric mean calc for preds
  etaSumSq_pgt  <- array( 0,  dim = c(nP, nG, nT) )  # sum of squared resids
  intNLLCpt_pgt <- array( 0,  dim = c(nP, nG, nT) )  # Internal likelihood component

  # Differences in likelihood values
  diffIntCpt_pg   <- array( 0, dim = c(nP, nG) )
  diffEtaSumSq_pg <- array( 0, dim = c(nP, nG) )
  diffTotLike_pg  <- array( 0, dim = c(nP, nG) )
  tau_pg          <- array( 0, dim = c(nP, nG) )
  difftau_pg      <- array( 0, dim = c(nP, nG) )


  for( g in 1:nG )
    for( p in 1:nP )
    {
      if( any( tcComps_apgt[,p,g,] > 0 ) )
      {
        for( t in 1:nT )
        {
          if( all(tcComps_apgt[,p,g,t] < 0 ) )
            next

          thisSampSize <- sum( tcComps_apgt[,p,g,t] )
          thisWt <- sqrt(meanA_pg[p,g] / thisSampSize)

          wtErr_pgt[p,g,t] <- log(thisWt) - log(ageWt_pgt[p,g,t])
          # Count the number of bin
          posBins <- which(tcComps_apgt[,p,g,t] > 0)
          nBins   <- length(posBins)
          nBins_pgt[p,g,t] <- nBins - 1

          # Calculate the covariance matrix
          Kmat <- matrix( 0, nrow = nBins - 1, ncol = nBins )
          Fmat <- matrix( 0, nrow = nBins - 1, ncol = nBins )
          Hmat <- matrix( 1, nrow = nBins - 1, ncol = nBins - 1)
          Kmat[1:(nBins - 1),1:(nBins-1)] <- diag(1, nBins - 1 )
          Fmat[1:(nBins - 1),1:(nBins-1)] <- diag(1, nBins - 1 )
          Kmat[,nBins] <- -1
          Fmat[,nBins] <- -1
          diag(Hmat) <- 2

          Corr    <- Corr_gaa[g,posBins,posBins]
          Vcheck  <- Kmat %*% Corr %*% t(Kmat)
          Vinv    <- ginv(Vcheck)
          logdetV <- log(det(Vcheck))
          Hinv    <- ginv(Hmat)

          Gamma   <- t(Fmat) %*% Hinv %*% Vcheck %*% Hmat %*% Fmat


          # Get obs and pred vecs without zeroes
          obsVec  <- tcComps_apgt[posBins,p,g,t]
          obsVec  <- obsVec / sum(obsVec)
          predVec <- tcPred_apgt[posBins,p,g,t]

          gmObsErr_pgt[p,g,t]   <-  (prod(obsVec)^(1/nBins))  - gmObs_pgt[p,g,t]
          gmPredErr_pgt[p,g,t]  <-  (prod(predVec)^(1/nBins)) - gmPred_pgt[p,g,t]

          obsY  <- obsVec[-nBins] / obsVec[nBins]
          predY <- predVec[-nBins] / predVec[nBins]

          # Calculate SSR
          etaSumSq_pgt[p,g,t] <- (t( log(obsY) - log(predY) ) %*% Vinv %*% ( log(obsY) - log(predY) )) / thisWt^2

          # Calculate internal component 
          intNLLCpt_pgt[p,g,t] <- 0.5 * logdetV + (nBins - 1) * log(thisWt);
        } 

        totBins <- sum(nBins_pgt[p,g,])

        # Compute internal likelihood component
        internalLike        <- sum( intNLLCpt_pgt[p,g,] )
        diffIntCpt_pg[p,g]  <- intrnlAgeLikCpt_pg[p,g] - internalLike
        
        # Compute sigma
        tau_pg[p,g]         <- sqrt(sum( etaSumSq_pgt[p,g,] ) / totBins )
        difftau_pg[p,g]     <- tau_pg[p,g] - sqrt(tau2Age_pg[p,g])

        # total Likelihood
        totLike             <- log(tau_pg[p,g]) * totBins + 
                                internalLike + 
                                0.5 * sum(etaSumSq_pgt[p,g,]) / (tau_pg[p,g]^2)

        diffEtaSumSq_pg[p,g] <- sum( etaSumSq_pgt[p,g,] ) - modelEtaSumSq_pg[p,g]
        diffTotLike_pg[p,g]  <- totLike - ageObsNLL_pg[p,g]

      }

    }
  
  years <- fYear:(fYear + nT - 1)

  errorCols <- RColorBrewer::brewer.pal(3, "Set2" )

  par(  mfrow = c(3,nP),
        mar = c(1,1,1,1),
        oma = c(3,4,2,3) )

  for( g in 1:3 )
    for( p in 1:nP )
    {
      plot( x = range(years), 
            y = range(wtErr_pgt,gmObsErr_pgt,gmPredErr_pgt,na.rm = T),
            type = "n", xlab = "", ylab = "",
            axes = FALSE )
        mfg <- par("mfg")
        if( mfg[1] == mfg[3] )
          axis(side = 1 )
        if( mfg[2] == 1 )
          axis( side = 2, las = 1 )
        if( mfg[1] == 1 )
          mtext( side = 3, text = stockLabs[p] )

        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = gearLabs[g] )
        box()
        # Plot points for each error seriess
        points( x = years, y = wtErr_pgt[p,g,],
                pch = 1, col = errorCols[1] )
        points( x = years, y = gmObsErr_pgt[p,g,],
                pch = 1, col = errorCols[2] )
        points( x = years, y = gmPredErr_pgt[p,g,],
                pch = 1, col = errorCols[3] )

        text( x = years[nT - 10], y = c(.5,.6,.7),
              labels = c( paste("diffNLL = ", round(diffTotLike_pg[p,g],2)),
                          paste("diffSSR = ", round(diffEtaSumSq_pg[p,g],2)),
                          paste("diffInt = ", round(diffIntCpt_pg[p,g],2)) ) )
    }

  mtext( side = 1, outer = TRUE, text = "Year" )
  mtext( side = 2, outer = TRUE, text = "Error" )

}

# Plot comp fits
plotMixedCompFitYrs <- function(  repList = reports,
                                  comps = "age",
                                  save = FALSE,
                                  savePath = "plotFitYrs" )
{
  repObj    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs

  # Pull predicted and observed ages
  if( comps == "age" )
  {
    max       <- repObj$nA
    pred_xtg  <- repObj$predMixedPA_atg[1:max,,, drop = FALSE]
    obs_xtg   <- repObj$mA_atg[1:max,,, drop = FALSE]  
    xLab      <- "Age"
  }
  

  dimNames  <- dimnames(obs_xtg)
  compNames <- dimNames[[1]]
  yearNames <- dimNames[[2]]
  gearNames <- dimNames[[3]]

  # Pull model dims
  nG      <- repObj$nG
  nT      <- repObj$nT  
  minPA   <- repObj$minPropAge

  # Make colours vector
  cols    <- brewer.pal( n = nG, "Dark2" )

  # Make years vector
  years   <- seq(initYear, length = nT+1, by = 1)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  obsGears <- c()
  gearTimes <- vector(mode = "list", length = nG)
  for( gIdx in 1:nG )
  {
    if( any(obs_xtg[1,,gIdx] >= 0) )
    {
      obsGears <- c(obsGears,gIdx)
      gearTimes[[gIdx]] <- which(obs_xtg[1,,gIdx] >= 0)
    }
  }

  # ok, obsGears are the ones we want to plot,
  # and gearTimes is the list of time indices
  for( gIdx in obsGears )
  {
    if(!save)
      dev.new()

    if(save)
    {
      gearPath <- paste(savePath,gearNames[gIdx],".png",sep = "")
      png(  gearPath, 
            width = 11, height = 8.5,
            units = "in", res = 300 )
    }

    times <- gearTimes[[gIdx]]
    # Count the number of age observations
    # there are, and make the plotting window
    nObs <- length(times)
    nCols <- round(sqrt(nObs))
    nRows <- ceiling(nObs/nCols)

    par(  mfcol = c(nRows,nCols), 
          mar = c(1,2,1,2),
          oma = c(3,3,3,3) )

    for( tIdx in times )
    { 
      Nobs <- sum(obs_xtg[,tIdx,gIdx])
      # get age obs and preds
      compObsProp_x  <- obs_xtg[,tIdx,gIdx]/sum(obs_xtg[,tIdx,gIdx])
      compPred_x     <- pred_xtg[,tIdx,gIdx]



      plot( x = c(1,max), y = c(0,max(compObsProp_x,compPred_x,na.rm = T) ),
            xlab = "", ylab = "", type = "n", las = 1 )
        rect( xleft = 1:max - .3, xright = 1:max + .3,
              ybottom = 0, ytop = compObsProp_x,
              col = "grey40", border = NA )
        lines(  x = 1:max, y = compPred_x, lwd = 1,
                col = cols[gIdx] )
        points(  x = 1:max, y = compPred_x,
                col = cols[gIdx], pch = 16 )
        abline( h = minPA, lty = 2, lwd = .8 )
        panLab( x=.5, y = .95, txt = years[tIdx] )
        panLab( x=.1, y = .95, txt = paste( "N = ", Nobs, sep = "" ) )

    }
    mtext( side = 3, outer = T, text = gearNames[gIdx], line = 2, font = 2)
    mtext(  side = 1, outer = T, text = xLab, line = 2 )
    mtext(  side = 2, outer = T, text = paste("Proportion-at-",xLab,sep=""), 
            line = 2 )
    if(save)
      dev.off()
  }

} # END plotCompFitYrs()


# plotRefpts_p
# 6 panel reference points plot from mseR
plotRefPts_p <- function(repList = reports, pIdx=1)
{
  repObj <- repList$repOpt
  stockNames <- dimnames(repObj$SB_pt)[[1]][pIdx]

  par(mfrow=c(3,2), mgp=c(2.1,0.4,0), mar=c(3,3.2,1,1), tck=-0.015,
      oma=c(0,0,1.5,0), cex.axis=1.3, cex.lab=1.3)

  .plotYprF_p(repList=repList, pIdx=pIdx )
  mtext('a)', adj=0)

  .plotSbPerRecF_p(repList=repList, pIdx=pIdx )
  mtext('b)', adj=0)

  .plotYieldF_p(repList=repList, pIdx=pIdx )
  mtext('c)', adj=0)

  .plotSbF_p(repList=repList, pIdx=pIdx ) 
  mtext('d)', adj=0)

  .plotRecSb_p(repList=repList, pIdx=pIdx )
  mtext('e)', adj=0)

  .plotYieldSb_p(repList=repList, pIdx=pIdx )
  mtext('f)', adj=0)

  mtext( side=3, line=0, cex=1, outer=TRUE, stockNames )


} # END plotRefpts_p()

# .plotYprF_p
# Plot equilibrium yield per recruit as a function of F
.plotYprF_p <- function( repList = reports, pIdx = 1 )
{
  repObj <- repList$repOpt
  repObj <- calcRefPts(repObj) 
  
  # Pull reference points objects
  refPoints <- repObj$refPts
  refCurves <- refPoints$refCurves
  FmsyRefPts<- refPoints$FmsyRefPts


  # Pull eqbm yields
  F         <- refCurves$F
  ypr_f     <- refCurves$ypr_pf[pIdx,]*1000
  Fmsy      <- FmsyRefPts$Fmsy_p[pIdx]
  yprFmsy   <- FmsyRefPts$yprFmsy_p[pIdx]*1000

  # extract values for SB>0
  ssb_f   <- refCurves$Beq_pf[pIdx,]
  fIdx    <- which(ssb_f>0)

  plot(x=F[fIdx], y=ypr_f[fIdx], type='l', ylab='Yield per recruit (t/million)', 
       xlab='F', las=1, lwd=2)
  points(x=Fmsy, y=yprFmsy, pch=21, bg='green', cex=2)


} # END .plotYprF_p()

# .plotSbPerRec_p
# Plot equilibrium spawning biomass per recruit as a function of F
.plotSbPerRecF_p <- function( repList = reports, pIdx = 1 )
{
  repObj <- repList$repOpt
  
  # recalculate to include ssbpr, since it isn't in older refpts calcs used
  repObj <- calcRefPts(repObj) 
  
  # Pull reference points objects
  refPoints <- repObj$refPts 
  refCurves <- refPoints$refCurves
  FmsyRefPts<- refPoints$FmsyRefPts


  # Pull eqbm yields
  F         <- refCurves$F
  ssbpr_f   <- refCurves$ssbpr_pf[pIdx,]*1000
  Fmsy      <- FmsyRefPts$Fmsy_p[pIdx]
  ssbprFmsy <- FmsyRefPts$ssbprFmsy_p[pIdx]*1000

  # extract values for SB>0
  ssb_f   <- refCurves$Beq_pf[pIdx,]
  fIdx    <- which(ssb_f>0)

  plot(x=F[fIdx], y=ssbpr_f[fIdx], type='l', ylab='SB per recruit (t/million)', 
       xlab='F', las=1, lwd=2)
  points(x=Fmsy, y=ssbprFmsy, pch=21, bg='green', cex=2)


} # END .plotSbPerRec_p()

# .plotYieldF_p
# Plot equilibrium yield reference curve as a function of F
.plotYieldF_p <- function( repList = reports, pIdx = 1 )
{
  repObj <- repList$repOpt
  repObj <- calcRefPts(repObj) 
  
  # Pull reference points objects
  refPoints <- repObj$refPts
  refCurves <- refPoints$refCurves
  FmsyRefPts<- refPoints$FmsyRefPts

  # Pull eqbm yields
  F         <- refCurves$F
  YeqF_f    <- refCurves$Yeq_pf[pIdx,]
  Fmsy      <- FmsyRefPts$Fmsy_p[pIdx]
  YeqFmsy   <- FmsyRefPts$YeqFmsy_p[pIdx]

  # extract values for SB>0
  ssb_f   <- refCurves$Beq_pf[pIdx,]
  fIdx    <- which(ssb_f>0)

  plot(x=F[fIdx], y=YeqF_f[fIdx], type='l', ylab='Yield (kt)', 
       xlab='F', las=1, lwd=2)
  points(x=Fmsy, y=YeqFmsy, pch=21, bg='green', cex=2)


} # END .plotYieldF_p()

# .plotSbF
# Plot equilibrium spawning biomass as a function of F
.plotSbF_p <- function( repList = reports, pIdx = 1 )
{
  repObj <- repList$repOpt
  repObj <- calcRefPts(repObj) 
  
  # Pull reference points objects
  refPoints <- repObj$refPts 
  refCurves <- refPoints$refCurves
  FmsyRefPts<- refPoints$FmsyRefPts


  # Pull eqbm yields
  F         <- refCurves$F
  ssb_f     <- refCurves$Beq_pf[pIdx,]
  Fmsy      <- FmsyRefPts$Fmsy_p[pIdx]
  ssbFmsy   <- FmsyRefPts$BeqFmsy_p[pIdx]

  # extract values for SB>0
  fIdx <- which(ssb_f>0)

  plot(x=F[fIdx], y=ssb_f[fIdx], type='l', ylab='Spawning biomass (kt)', 
       xlab='F', las=1, lwd=2)
  points(x=Fmsy, y=ssbFmsy, pch=21, bg='green', cex=2)


} # END .plotSbF()

# .plotRecSb
# Plot equilibrium recruits vs SB
.plotRecSb_p <- function( repList = reports, pIdx = 1 )
{
  repObj <- repList$repOpt
  repObj <- calcRefPts(repObj) 
  
  # Pull reference points objects
  refPoints <- repObj$refPts 
  refCurves <- refPoints$refCurves
  FmsyRefPts<- refPoints$FmsyRefPts

  # Pull eqbm yields
  R_f       <- refCurves$Req_pf[pIdx,]
  ssb_f     <- refCurves$Beq_pf[pIdx,]
  ReqFmsy   <- FmsyRefPts$ReqFmsy_p[pIdx]
  ssbFmsy   <- FmsyRefPts$BeqFmsy_p[pIdx]

  # extract values for SB>0
  fIdx <- which(ssb_f>0)

  plot(x=ssb_f[fIdx], y=R_f[fIdx], type='l', ylab='Recruits (millions)', 
       xlab='Spawning biomass (kt)', las=1, lwd=2)
  points(x=ssbFmsy, y=ReqFmsy,  pch=21, bg='green', cex=2)


} # END .plotRecSb_p()

# .plotYieldSb_p
# Plot equilibrium Yield vs SB
.plotYieldSb_p <- function( repList = reports, pIdx = 1 )
{
  repObj <- repList$repOpt
  repObj <- calcRefPts(repObj) 
  
  # Pull reference points objects
  refPoints <- repObj$refPts 
  refCurves <- refPoints$refCurves
  FmsyRefPts<- refPoints$FmsyRefPts

  # Pull eqbm yields
  yield_f   <- refCurves$Yeq_pf[pIdx,]
  ssb_f     <- refCurves$Beq_pf[pIdx,]
  YeqFmsy   <- FmsyRefPts$YeqFmsy_p[pIdx]
  ssbFmsy   <- FmsyRefPts$BeqFmsy_p[pIdx]

  # extract values for SB>0
  fIdx <- which(ssb_f>0)

  plot(x=ssb_f[fIdx], y=yield_f[fIdx], type='l', ylab='Yield (kt)', 
       xlab='Spawning biomass (kt)', las=1, lwd=2)
  points(x=ssbFmsy, y=YeqFmsy,  pch=21, bg='green', cex=2)


} # END .plotYieldSb_p()




# plotYeq
# Plot equilibrium yield reference curve as a function
# of F
plotYeqF <- function( repList = reports,
                      pIdx = 1 )
{
  repObj <- repList$repOpt
  # Pull reference points object
  refPoints <- repObj$refPts
  refCurves <- refPoints$refCurves

  stockNames <- dimnames(repObj$SB_pt)[[1]]

  # Pull eqbm yields
  YeqF_pf     <- refCurves$Yeq_p[pIdx,,drop = FALSE]
  Fmsy_p      <- refPoints$FmsyRefPts$Fmsy_p[pIdx, drop = FALSE]
  YeqFmsy_p   <- refPoints$FmsyRefPts$YeqFmsy_p[pIdx, drop = FALSE]

  dimnames(YeqF_pf)[[1]] <- stockNames[pIdx]
  names(Fmsy_p) <- stockNames[pIdx]
  names(YeqFmsy_p) <- stockNames[pIdx]

  # Now melt
  YeqF.df     <- melt(YeqF_pf) %>%
                  rename(yield = value) %>%
                  filter( yield >= 0 )
  Fmsy.df     <- melt(Fmsy_p) %>%
                  rename( Fmsy = value) %>%
                  mutate( stock = stockNames[pIdx])
  YeqFmsy.df  <- melt(YeqFmsy_p) %>%
                  rename( yield = value) %>%
                  mutate( stock = stockNames[pIdx]) %>%
                  left_join(Fmsy.df, by = "stock" )

  # Figure out how to stack multiple YPR/SPR ref point
  # tables together, with a column for the ref point
  # label, for easier plotting of multiple ref points

  tmp <-  ggplot(data = YeqF.df, aes(x=F, y=yield)) + 
          geom_line() +
          facet_grid( stock ~ ., scales = "free_y") +
          theme_sleek() +
          geom_point( data = YeqFmsy.df, inherit.aes = FALSE,
                      mapping = aes( x = Fmsy, y = yield ),
                      col = "red", size = 1.5 )

  print(tmp)

} # END plotYeqF()

# plotAvgAgeFitComparison()
# Comparison of age fits, averaged over time, 
# for ISCAM and SISCA.
plotAvgAgeFitComparison <- function(  repObj      = reports$repOpt,
                                      iscamRep    = ISCAMrep,
                                      fYearSISCA  = reports$fYear,
                                      lYearSISCA  = reports$lYear )
{
  # Pull predicted and observed ages
  predAge_apgt <- repObj$predPA_apgt
  obsAge_apgt  <- repObj$A_apgt

  gearLabs  <- dimnames(obsAge_apgt)[[3]]

  # replace missing entries with NAs
  predAge_apgt[predAge_apgt < 0] <- NA
  obsAge_apgt[obsAge_apgt < 0] <- NA

  fYearISCAM <- iscamRep$syr
  lYearISCAM <- iscamRep$nyr


  # Get time steps
  nT <- dim(predAge_apgt)[4]
  nG <- dim(predAge_apgt)[3]
  nP <- dim(predAge_apgt)[2]


  # Multiply the yearly predicted proportions by the
  # yearly sample sizes
  for( g in 1:nG )
    for( p in 1:nP )
      for( t in 1:nT )
      {
        thisSamp <- sum( obsAge_apgt[,p,g,t], na.rm = T )
        if( thisSamp < 0 )
          thisSamp <- 0

        predAge_apgt[,p,g,t] <- thisSamp * predAge_apgt[,p,g,t]
      }

  
  
  # Average over time
  predAge_ag <- apply( X = predAge_apgt, FUN = sum, MARGIN = c(1,3), na.rm = T )
  obsAge_ag <- apply( X = obsAge_apgt, FUN = sum, MARGIN = c(1,3), na.rm = T )    

  
  # Pull model dims
  nG      <- repObj$nG
  nA      <- repObj$nA
  nT      <- repObj$nT
  nP      <- repObj$nP
  minPA   <- repObj$minPropAge

  # Now pull ISCAM ages
  iscamObsAge   <- iscamRep$d3_A
  iscamPredAge  <- iscamRep$A_hat

  # Gotta rearrange
  iscamObsAge_agt   <- array( NA, dim = c(nA,nG,nT) )
  iscamPredAge_agt  <- array( NA, dim = c(nA,nG,nT) )

  for( k in 1:nrow(iscamObsAge) )
  {
    iscamYr <- iscamObsAge[k,1]

    if( iscamYr < fYearSISCA | iscamYr > lYearSISCA )
      next

    gearIdx <- iscamObsAge[k,2]
    yrIdx   <- iscamYr - fYearSISCA + 1

    ageObs <- iscamObsAge[k,6:14]
    sampSize <- sum(ageObs)

    # browser()

    iscamObsAge_agt[2:10,gearIdx,yrIdx] <- ageObs
    iscamPredAge_agt[2:10,gearIdx,yrIdx] <- sampSize * iscamPredAge[k,]

  }

  iscamObsAge_ag  <- apply( X = iscamObsAge_agt, FUN = sum, MARGIN = c(1,2), na.rm = T)
  iscamPredAge_ag <- apply( X = iscamPredAge_agt, FUN = sum, MARGIN = c(1,2), na.rm = T)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  ageGears <- c()
  gearTimes <- vector(mode = "list", length = nG)
  for( gIdx in 1:nG )
    if( sum(obsAge_ag[,gIdx]) > 0  )
      ageGears <- c(ageGears,gIdx)

  # Make colours vector
  cols    <- RColorBrewer::brewer.pal( n = nG, "Dark2" )


  par(  mfrow = c(length(ageGears),2), 
        mar = c(1,2,1,2),
        oma = c(3,3,3,3) )

  # ok, ageGears are the ones we want to plot,
  for( aIdx in 1:length(ageGears) )
  { 
    gIdx <- ageGears[aIdx]
    # get age obs and preds for SISCA
    ageObs  <- obsAge_ag[,gIdx]
    nSamp   <- sum(ageObs)
    ageObs  <- ageObs / sum(ageObs)
    agePred <- predAge_ag[,gIdx]
    agePred <- agePred / sum( agePred )

    plot( x = c(1,nA), y = c(0,1.1*max(ageObs,agePred,na.rm = T) ),
          xlab = "", ylab = "", type = "n", las = 1 )
      legend( x = "topright",
              bty = "n",
              legend = paste("N = ", nSamp, sep = "" ) )
      rect( xleft = 3:nA - .3, xright = 3:nA + .3,
            ybottom = 0, ytop = ageObs[3:nA],
            col = "grey40", border = NA )
      lines(  x = 3:nA, y = agePred[3:nA], lwd = 3,
              col = cols[gIdx] )
      points(  x = 3:nA, y = agePred[3:nA],
              col = cols[gIdx], pch = 16, cex = 1.5 )
      abline( h = minPA, lty = 2, lwd = .8 )
      
      
      # panLab( x=.5, y = .95, txt = gearLabels[gIdx] )
      mfg <- par("mfg")
      if(mfg[1] == 1 )
        mtext( side = 3, text = "SISCAH", line = 1, font = 2)

    # get age obs and preds for ISCAM
    ageObs  <- iscamObsAge_ag[,gIdx]
    nSamp   <- sum(ageObs)
    ageObs  <- ageObs / sum(ageObs)
    agePred <- iscamPredAge_ag[,gIdx]
    agePred <- agePred / sum( agePred )

    plot( x = c(1,nA), y = c(0,1.1*max(ageObs,agePred,na.rm = T) ),
          xlab = "", ylab = "", type = "n", las = 1 )
      legend( x = "topright",
              bty = "n",
              legend = paste("N = ", nSamp, sep = "" ) )
      rect( xleft = 2:nA - .3, xright = 2:nA + .3,
            ybottom = 0, ytop = ageObs[2:nA],
            col = "grey40", border = NA )
      lines(  x = 2:nA, y = agePred[2:nA], lwd = 3,
              col = cols[gIdx] )
      points(  x = 2:nA, y = agePred[2:nA],
              col = cols[gIdx], pch = 16, cex = 1.5 )
      abline( h = 0.02, lty = 2, lwd = .8 )
      
      
      # panLab( x=.5, y = .95, txt = gearLabels[gIdx] )
      mfg <- par("mfg")
      if(mfg[1] == 1 )
        mtext( side = 3, text = "ISCAM", line = 1, font = 2)


      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text( x = corners[2]+0.5, 
              y = mean(corners[3:4]), 
              labels = gearLabs[gIdx], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
  }
  
  mtext( side = 1, outer = T, text = "Age", line = 2 )
  mtext( side = 2, outer = T, text = "Proportion", line = 2 )

} # END plotAvgAgeFitComparison()


# plotRetroSBt()
# Takes in a report list, and 
# plots the spawning biomass retrospective
# plots
plotRetroSBt <- function( reportList = repList,
                          dep = FALSE,
                          relative = FALSE,
                          refIdx = 1,
                          legendText = NULL  )
{

  nFits     <- length(reportList)
  nTs       <- numeric( length = nFits )
  fYears    <- numeric( length = nFits )
  lYears    <- numeric( length = nFits )


  # Get number of stocks
  nP <- reportList[[1]]$repOpt$nP

  # Calculate Mohn's rho
  # Number of peels is nFits - 1
  Bpeel_pf  <- array( NA, dim = c(nP,nFits) )
  Bref_pf   <- array( NA, dim = c(nP,nFits) )

  stockNames <- reportList[[1]]$stock

  # Get nTs and model start/end years
  for( k in 1:nFits )
  {
    fYears[k] <-  reportList[[k]]$fYear
    lYears[k] <-  reportList[[k]]$lYear
  }
  nTs    <- lYears - fYears + 1
  maxT <- max(nTs)

  minfYear <- min(fYears)
  maxlYear <- max(lYears)

  if( all(fYears == minfYear ))
    retroType <- "end"
  if( all(lYears == maxlYear ))
    retroType <- "start"

  # Make an array for holding biomass
  bio_kpt <- array(NA, dim = c(nFits,nP,maxT) )

  for( k in 1:nFits )
  {
    thisfYear <- fYears[k]
    thislYear <- lYears[k]

    tdx <- thisfYear:thislYear - minfYear + 1

    bio_kpt[k,,tdx] <- reportList[[k]]$repOpt$SB_pt[,1:nTs[k]]

    # if(retroType == "end")
    # {
    #   Bpeel_pf[,k] <- reportList[[k]]$repOpt$SB_pt[,nTs[k]]
    #   Bref_pf[,k]  <- reportList[[refIdx]]$repOpt$SB_pt[,nTs[k]]
    # }

    # if(retroType == "start")
    # {
    #   Bpeel_pf[,k] <- reportList[[k]]$repOpt$SB_pt[,1]
    #   Bref_pf[,k]  <- reportList[[refIdx]]$repOpt$SB_pt[,tdx[1]]
    # }

    # Make depletion if asked for
    if( dep )
    {
      for( p in 1:nP)
        bio_kpt[k,p,tdx] <- bio_kpt[k,p,tdx] / reportList[[k]]$repOpt$B0_p[p]
    }
  }

  # mohnBerr_pf <- (Bpeel_pf - Bref_pf) / Bref_pf
  # mohnBerr_pf <- mohnBerr_pf[,-refIdx]

  # rho_p       <- round(apply( X = mohnBerr_pf, FUN = mean, MARGIN = 1, na.rm = T),3)

  # Replace zeroes with NAs
  bio_kpt[bio_kpt == 0 ] <- NA

  # Pull reference out
  refB_pt <- array(NA, dim = c(nP,maxT))
  refB_pt[1:nP,] <- bio_kpt[refIdx,1:nP,]

  if( relative )
  {
    for( k in 1:nFits )
      bio_kpt[k,,] <- bio_kpt[k,,] / refB_pt
  
    refB_pt[1:nP,] <- refB_pt[1:nP,] / refB_pt[1:nP,]
  }


  # And remove from biomasss
  bio_kpt <- bio_kpt[-refIdx,,,drop = FALSE]

  # Pull out colours

  caseCols <- scales::viridis_pal(alpha = 1, 
                                  begin = 0,
                                  end = 1,
                                  direction = -1)(nFits-1)

  par( mfrow = c( nP,1 ), mar = c(.1,1,.1,2),
        oma = c(3,4,2,2) )

  for( p in 1:nP )
  {
    plot( x = c(minfYear,maxlYear),
          y = c(0,max(1,bio_kpt[,p,],refB_pt[p,],na.rm = T)),
          type = "n", axes = FALSE )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )

      if( mfg[2] == mfg[4] )
      {
        rmtext( txt = stockNames[p], outer = TRUE,
                font = 2, line = .02, cex = 1.5)
      }

      axis( side = 2, las = 1 )
      box()
      grid()


      # rhoLab <- bquote(rho == .(rho_p[p]) )

      # text( x = 1960, y = 1 * max(bio_kpt[,p,],na.rm = T),
      #       labels = rhoLab )

      # Now do the reference model
      lines(  x = minfYear:maxlYear, y = refB_pt[p,],
              col = "black", lwd = 3 )

      # And retro lines
      for( k in 1:(nFits-1) )
      {
        lines( x = minfYear:maxlYear, y = bio_kpt[k,p,],
                col = caseCols[k], lwd = 2 )
      }

      
  }

  if(!is.null(legendText))
  {
    legOrder <- order(legendText)
    legCols <- character(length = nFits)
    legCols[refIdx] <- "black"
    legCols[-refIdx] <- caseCols

    legend( x = "topright", bty = "n",
            legend = legendText[legOrder],
            lty = 1, lwd = 2,
            col = legCols )

  }

  yLab <- "Spawning Biomass (kt)"
  if( dep )
    yLab <- expression( paste("Spawning Biomass Depletion (", SB[t]/SB[0], ")", sep = "") )

  mtext( outer = TRUE, line = 2, "Year", side = 1 )
  mtext( outer = TRUE, line = 2, yLab,   side = 2 )


} # END plotRetroSBt()


# plotRetroRt()
# Takes in a report list, and 
# plots the spawning biomass retrospective
# plots
plotRetroRt <- function(  reportList = repList,
                          dep = FALSE,
                          relative = FALSE,
                          refIdx = 1,
                          legendText = NULL )
{

  nFits     <- length(reportList)
  nTs       <- numeric( length = nFits )
  fYears    <- numeric( length = nFits )
  lYears    <- numeric( length = nFits )

  # Get number of stocks
  nP <- reportList[[1]]$repOpt$nP

  stockNames <- reportList[[1]]$stock


  # Get nTs and model start/end years
  for( k in 1:nFits )
  {
    fYears[k] <-  reportList[[k]]$fYear
    lYears[k] <-  reportList[[k]]$lYear
  }
  nTs    <- lYears - fYears + 1
  maxT <- max(nTs)

  minfYear <- min(fYears)
  maxlYear <- max(lYears)

  # Make an array for holding biomass
  R_kpt <- array(NA, dim = c(nFits,nP,maxT) )

  for( k in 1:nFits )
  {
    thisfYear <- fYears[k]
    thislYear <- lYears[k]

    tdx <- thisfYear:thislYear - minfYear + 1


    R_kpt[k,,tdx] <- reportList[[k]]$repOpt$R_pt[,1:nTs[k]]

    # Make depletion if asked for
    if( dep )
    {
      for( p in 1:nP)
        R_kpt[k,p,tdx] <- R_kpt[k,p,tdx] / reportList[[k]]$repOpt$R0_p[p]
    }
  }

  R_kpt[R_kpt == 0 ] <- NA

  refR_pt <- array(NA, dim = c(nP,maxT))
  refR_pt[1:nP,] <- R_kpt[refIdx,1:nP,]

  if( relative )
  {
    

    for( k in 1:nFits )
      R_kpt[k,,] <- R_kpt[k,,] / refR_pt

  }
  R_kpt <- R_kpt[-refIdx,,,drop = FALSE]

  # Pull out colours
  caseCols <- RColorBrewer::brewer.pal(nFits - 1,"Set3")

  par( mfrow = c( nP,1 ), mar = c(.1,1,.1,2),
        oma = c(3,4,2,2) )

  for( p in 1:nP )
  {
    plot( x = c(minfYear,maxlYear),
          y = c(0,max(1,R_kpt[,p,],refR_pt[p,],na.rm = TRUE)),
          type = "n", axes = FALSE )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text( x = corners[2]+1, 
              y = mean(corners[3:4]), 
              labels = stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }

      axis( side = 2, las = 1 )
      box()
      grid()

      # And retro lines
      for( k in 1:(nFits-1) ) 
        lines( x = minfYear:maxlYear, y = R_kpt[k,p,],
                col = caseCols[k], lwd = 2 )

      # Now do the reference case
      lines(  x = minfYear:maxlYear, y = R_kpt[refIdx,p,],
              col = "black", lwd = 3 )
      

  }

  if(!is.null(legendText))
  {
    legOrder <- order(legendText)
    legend( x = "topright", bty = "n",
            legend = legendText[legOrder],
            lty = 1, lwd = 2,
            col = c("black",caseCols)[legOrder] )
  }

  yLab <- "Recruitment (1e6)"
  if( dep )
    yLab <- expression( paste("Relative Recruitment  (", R[t]/R[0], ")", sep = "") )
  if( relative )
    yLab <- expression(paste( "Recruitment relative to reference model, (", R[t]/Rref[t] ,")" ))

  mtext( outer = TRUE, line = 2, "Year", side = 1 )
  mtext( outer = TRUE, line = 2, yLab,   side = 2 )


} # END plotRetroRt

# plotACF()
# Plots autocorrelation estimates for
# an array of time series.
plotACF <- function(arr_pt, max.lag = 10, maxy = .5)
{
  nP <- dim(arr_pt)[1]

  # Now plot indicators of +ve or 0 idx
  stockCols <- scales::viridis_pal(option = "C",alpha = 1,
                                    begin = 0, end = .5)(3)

  calcACF <- function(idx, arr_pt)
  {
    xVec <- arr_pt[idx,]

    ACF <- acf(xVec, na.action = na.pass, plot = FALSE)
    ACF
  }

  acfList <- lapply( X = 1:nP, FUN = calcACF, arr_pt = arr_pt)
  xJitter <- seq( from = -.15, to = .15, length.out = nP)

  plot( x = c(0, max.lag), y = c(-maxy, maxy),
        xlab = "", ylab = "", axes = FALSE, type = "n" )
  mfg <- par("mfg")
  if( mfg[1] == mfg[3])
    axis( side = 1, at = 0:max.lag)

  if( mfg[2] == 1)
    axis( side = 2, las = 1 )

  if( mfg[1] == 1 )
    legend( x = "topright",
            legend = c("C/S","JP/S","Lou"),
            lty = 1, lwd = 2,
            col = stockCols, bty = "n")

  box()
  grid()
  
  abline(h = 0, lty = 2, col = "grey30")

  for( p in 1:nP )
  {
    segments( x0 = acfList[[p]]$lag + xJitter[p],
              y0 = 0, y1 = acfList[[p]]$acf,
              col = stockCols[p], lwd = 2 )
  }
  # abline(v = .5 + c(0:10), lty = 1, lwd = 1)
} # END plotACF()

plotDeltaLNFits <- function( reports = reports )
{
  # Get reports and fYear
  repObj <- reports$repOpt
  fYear  <- reports$fYear
  datObj <- reports$data

  # Get dimension
  nP <- repObj$nP
  nG <- repObj$nG
  nT <- repObj$nT

  stock <- reports$stock



  # Get indicator vector for delta model
  deltaIdx_pg <- datObj$deltaIdx_pg
  deltaIdx_g  <- apply( X = deltaIdx_pg, FUN = sum, MARGIN = 2)
  deltaIndVar <- datObj$deltaIndVar

  # Pull biomass, probability of positive indices,
  B0_p              <- repObj$B0_p
  SB_pt             <- repObj$SB_pt
  probPosIdx_pgt    <- repObj$probPosIdx_pgt
  I_pgt             <- repObj$I_pgt

  probPosCombIdx_pt <- repObj$probPosCombIdx_pt


  maxBt <- max(SB_pt[,38:69], na.rm = T)
  BtSeq <- seq(from = 0, to = maxBt, length.out = 100)



  # Compute all probability lines for the posterior
  indI_pgt <- I_pgt
  indI_pgt[I_pgt > 0] <- 1
  indI_pgt[I_pgt < 0] <- NA


  if( any(datObj$combI_pt >= 0) )
    indI_pgt[indI_pgt == 0] <- NA

  indCombI_pt <- datObj$combI_pt
  indCombI_pt[indCombI_pt > 0] <- 1
  indCombI_pt[indCombI_pt < 0] <- NA

  # get model parameters too
  meanProbPosIdx_pg <- repObj$meanProbPosIdx_pg
  SDProbPosIdx_pg   <- repObj$SDProbPosIdx_pg

  if(!is.null(reports$posts))
  {
    probPosIdx_ipgt <- reports$posts$probPosIdx_ipgt
    nPosts <- dim(probPosIdx_ipgt)[1]
    probPosIdx_pgt  <- apply( X = probPosIdx_ipgt, FUN = mean, MARGIN = c(2,3,4))

    meanProbPosIdx_ipg  <- reports$posts$meanProbPosIdx_ipg
    SDProbPosIdx_ipg    <- reports$posts$SDProbPosIdx_ipg

    

    meanProbPosIdx_pg <- apply( X = reports$posts$meanProbPosIdx_ipg, FUN = mean, MARGIN = c(2,3))
    SDProbPosIdx_pg   <- apply( X = reports$posts$SDProbPosIdx_ipg  , FUN = mean, MARGIN = c(2,3))

    modelProbPos_igb <- array( 0, dim = c(nPosts,nG,100) )
    for(i in 1:nPosts )
    {
      for( g in 1:nG )
        modelProbPos_igb[i,g,] <- 1 / (1 + exp( - (meanProbPosIdx_ipg[i,1,g] + SDProbPosIdx_ipg[i,1,g] * BtSeq)))
    }

    modelProbPos_qgb <- apply( X = modelProbPos_igb, FUN = quantile, na.rm = T,
                                MARGIN = c(2,3), probs = c(0.025, 0.5, 0.975) )

    meanModelProb_gb <- apply( X = modelProbPos_igb, FUN = mean,
                                MARGIN = c(2,3), na.rm = T )
  } 

  probPosIdx_pgt[I_pgt < 0] <- NA




  par(  mar = c(.1,.1,.1,.1),
        oma = c(3,5.5,1,3) )


  fleetCols <- RColorBrewer::brewer.pal( nG, "Dark2" )
  fleetPch  <- 1:nG

  # Now plot indicators of +ve or 0 idx
  stockCols <- scales::viridis_pal(option = "C",alpha = 1,
                                    begin = 0, end = .5)(nP)

  plot( x = c(0,30),y = c(0,1),
        type = "n", axes = FALSE, xaxs = "i" ) 
    mfg <- par("mfg")
    if(mfg[1] == mfg[3])
    {
      axis( side = 1 )
      mtext(  side = 1, text = "Biomass",
              line = 2)
    }

    axis( side = 2, las = 1)

    box()
    grid()
    deltaFleets <- which(deltaIdx_g > 0)
    nDeltaFleets <- length(deltaFleets)

    # Loop over gears, 
    for( g in deltaFleets)
    {
      # plot model of detection prob
      if(is.null(reports$posts))
      {
        indVar <- BtSeq
        if(deltaIndVar == 2)
          indVar <- BtSeq / B0_p[p]
        modelProb <- 1 / (1 + exp( - (meanProbPosIdx_pg[1,g] + SDProbPosIdx_pg[1,g] * indVar)))
      } else {
        modelProb <- meanModelProb_gb[g,]
        polygon(  x = c(BtSeq,rev(BtSeq)),
                  y = c( modelProbPos_qgb[1,g,], rev(modelProbPos_qgb[3,g,]) ),
                  col = "grey70",
                  border = NA )
      }

      lines(  x = BtSeq, y = modelProb,
              col = fleetCols[g], lwd = 2 )

      for( p in 1:nP )
      {
        # Get ticks if using split indices
        if( any(!is.na(indI_pgt)))
        {
          zeroTicks <- SB_pt[p,which(indI_pgt[p,g,] == 0)] 
          oneTicks <- SB_pt[p,which(indI_pgt[p,g,] == 1)] 
        }

        if( any(!is.na(indCombI_pt)))
        {
          zeroTicks <- SB_pt[p,which(indCombI_pt[p,] == 0)] 
          oneTicks <- SB_pt[p,which(indCombI_pt[p,] == 1)]  
        }

        axis( side = 1, at = zeroTicks, labels = FALSE,
              tck = .02, col.ticks = stockCols[p])
        axis( side = 3, at = oneTicks, labels = FALSE,
              tck = .02, col.ticks = stockCols[p])
       
      }

      
    }

      
    legend( x = "bottomright",bty = "n",
            legend = c(stock,"Dive"),
            col = c(stockCols,fleetCols[deltaFleets]),
            pch = c(NA,NA,NA,22),
            pt.bg = "grey70",
            pt.lwd = 0,
            lty = 1,
            lwd = 1,
            pt.cex = 1.5 )

  mtext( side = 2, text = "Probability of\ndetecting spawn",
         outer = TRUE, line = 3, font = 2 )

}

plotIratio <- function( repList = reports,
                        gIdx = c(4,5) )
{
  datObj <- repList$data

  I_pgt     <- datObj$I_pgt
  rI_pgt    <- datObj$rI_pgt
  combI_pt  <- datObj$combI_pt

  stockLabs <- dimnames(I_pgt)[[1]]
  gearLabs  <- dimnames(I_pgt)[[2]]
  nT        <- dim(I_pgt)[3]
  nP        <-  dim(I_pgt)[1]

  yrs <- 1951:2019

  combI_pt[combI_pt <= 0] <- NA

  fleetCols <- brewer.pal(n = 6, "Dark2")

  par( mfrow = c(nP,1), oma = c(4,4,2,3), mar = c(.1,1,.1,1) )

  for( p in 1:nP)
  {
    maxIdx <- max(combI_pt[p,],na.rm = T)

    ratTicks <- seq(from = 0, to = maxIdx, length.out = 5)
    ratLabs  <- (0:4)/4
    plotRat <- rI_pgt[p,4,] * maxIdx
    plot( x = range(yrs), y = c(0,maxIdx),
          axes = FALSE, type = "n" )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      axis( side = 2, las = 1 )
      box()
      grid()
      # Surface survey
      rect( xleft = yrs - .3,
            xright = yrs + .3,
            ybottom = 0,
            ytop = combI_pt[p,] * rI_pgt[p,4,],
            col = fleetCols[4], border = NA )
      # Dive survey
      rect( xleft = yrs - .3,
            xright = yrs + .3,
            ybottom = combI_pt[p,] * rI_pgt[p,4,],
            ytop = combI_pt[p,],
            col = fleetCols[5], border = NA )
      if( p == 1 )
        legend( x = "topright", bty = "n",
                legend = c("Surface","Dive"),
                col = fleetCols[4:5],
                pch = 22,
                pt.bg = fleetCols[4:5],
                cex = 2)

    rmtext( txt = stockLabs[p], font = 2, line = 0.02, cex = 1.5, outer = TRUE )
  }
  mtext( side = 1, outer = TRUE, text =  "Year", line = 2)
  mtext( side = 2, outer = TRUE, text =  "Spawn Index (kt)", line = 2)
}

# plotJuveMLikeProfile()
plotJuveMLikeProfile <- function( groupFolder = "juveM_likeProf",
                                  prefix = "juveM" )
{ 
  # Load batch reports
  batchReports <- .loadBatch( groupFolder = groupFolder,
                            baseDir = "Outputs/fits",
                            prefix = prefix )
  # Load the base OM
  .loadFit("fit_MSbase")

  baseJuveM <- reports$repOpt$M_apt[1,2,1]


  baseModelLike <-  sum(reports$repOpt$obsIdxNLL_pg) +
                    sum(reports$repOpt$obsIdxDeltaNLL_pg) +
                    sum(reports$repOpt$obsCombIdxNLL_p) +
                    sum(reports$repOpt$obsCombIdxDeltaNLL_p) +
                    sum(reports$repOpt$ageObsNLL_pg) +
                    sum(reports$repOpt$obsMixedIdxNLL_g) +
                    sum(reports$repOpt$ageObsNLL_g)



  reports$repOpt$objFun
  baseB0    <- sum(reports$repOpt$B0_p)

  nModels <- length(batchReports)
  # Now get x value and likelihood function values
  juveM_m     <- numeric(length = nModels)
  modelLike_m <- numeric(length = nModels)
  B0_m        <- numeric(length = nModels)

  for( k in 1:nModels)
  {
    rep_k <- batchReports[[k]]$repOpt
    juveM_m[k]      <-  rep_k$M_apt[1,2,1]
    modelLike_m[k]  <-  sum(rep_k$obsIdxNLL_pg) +
                        sum(rep_k$obsIdxDeltaNLL_pg) +
                        sum(rep_k$obsCombIdxNLL_p) +
                        sum(rep_k$obsCombIdxDeltaNLL_p) +
                        sum(rep_k$ageObsNLL_pg) +
                        sum(rep_k$obsMixedIdxNLL_g) +
                        sum(rep_k$ageObsNLL_g)
    B0_m[k]         <-  sum(rep_k$B0_p)
  }

  modelLike_m <- modelLike_m[order(juveM_m)]
  B0_m <- B0_m[order(juveM_m)]
  juveM_m <- juveM_m[order(juveM_m)]


  minLikeIdx <- which.min(modelLike_m)

  likeDiff <- modelLike_m - baseModelLike

  # Fit a spline to likeDiff
  likeDiffSpline <- splinefun(  x = juveM_m, y = modelLike_m - baseModelLike )

  bounds <- c(juveM_m[minLikeIdx],2)
  zeroLikeDiff  <- uniroot(f = likeDiffSpline, interval = bounds)$root
  minLikeM      <- uniroot( f = likeDiffSpline, interval = c(0.6,2), deriv = 1)$root

  browser()
  
  # Now make a spline for plotting
  likeSpline <- splinefun(  x = juveM_m, y = modelLike_m )
  plotM <- seq(from = 0.6, to = 2, by = 0.01)
  plotLike <- likeSpline(plotM)

  par(mfrow = c(2,1), mar = c(.1,.1,.1,.1), oma = c(3,5,2,2))
  plot( x = c(0.6,2),
        y = range(modelLike_m),
        las = 1, type = "n", xlab = "Age-1 M (/yr)",
        ylab = "Model minimum negative log-posterior",
        axes = FALSE)
    lines( x = plotM, y = plotLike, lty = 1, lwd = 3)
    points( x = baseJuveM, y  = baseModelLike, pch = 24, bg = "green", col = "black", cex = 2)
    points( x = minLikeM, y = likeSpline(minLikeM), 
            pch = 21, bg = "red", cex = 2 )
    abline( h = baseModelLike, lty = 2, lwd = 2, col = "steelblue")
    axis( side = 2, las = 1)
    box()
    points( x = zeroLikeDiff, y = likeSpline(zeroLikeDiff), 
            pch = 21, bg = "steelblue", lwd = 1, cex = 2 )
    mtext( side = 2, text = "Model minimum negative log-likelihood", line = 3.5)

    legend( x = "topleft", bty = "n",
            legend = c( "Base OM, Age-1 M = 0.72",
                        paste("Min NLL, Age-1M =",round(minLikeM,2)),
                        paste("Age-1 M =", round(zeroLikeDiff,2)),
                        "Base OM NLL"),
            pch = c(24,21,21, NA),
            lwd = c(NA,NA,NA,2),
            lty = c(NA,NA,NA,2),
            col = c("black","black","black","steelblue"),
            pt.bg = c("green","red","steelblue",NA),
            pt.lwd = 1,
            cex = 1,
            pt.cex = 2 )

  B0spline <- splinefun(  x = juveM_m, y = B0_m )

  plot( x = c(0.6,2),
        y = range(B0_m),
        las = 1, type = "n", xlab = "Age-1 M (/yr)",
        ylab = "Aggregate unfished biomass (kt)", axes = FALSE)
    axis( side = 1 )
    axis( side = 2, las = 1)
    box()
    lines( x = juveM_m, y = B0_m, lty = 1, lwd = 3)
    points( x = baseJuveM, y  = baseB0, pch = 24, bg = "green", cex = 2)
    points( x = minLikeM, y = B0spline(minLikeM), pch = 21, bg = "red", cex = 2 )
    points( x = zeroLikeDiff, y = B0spline(zeroLikeDiff), pch = 21, 
            bg = "steelblue", lwd = 1, cex = 2 )
    mtext( side = 2, text = "Aggregate unfished biomass (kt)", line = 3.5)
    mtext( side = 1, text = "Age-1 Natural Mortality (/yr)", line = 2)



} # END plotJuveMLikeProfile

# PlotAggBtRtComparison()
# Comparison of SISCA and ISCAM Bt and age-2 Rt at the major stock
# level.
plotAggBtRtComparison <- function(  repObj      = reports$repOpt,
                                    datObj      = reports$data,
                                    iscamRep    = ISCAMrep,
                                    fYearSISCA  = reports$fYear,
                                    lYearSISCA  = reports$lYear,
                                    simpleBt    = FALSE )
{
  # Get ISCAM years
  fYearISCAM <- iscamRep$syr
  lYearISCAM <- iscamRep$nyr

  nP <- repObj$nP
  nT <- repObj$nT
  nG <- repObj$nG

  if( any( datObj$combI_pt > 0 ) )
  {
    combIdx <- TRUE
    qComb_pt <- repObj$qComb_pt
  } else combIdx <- FALSE


  # Get SBt, sum
  SB_pt <- repObj$SB_pt[1:nP,,drop = FALSE]
  
  SB_t  <- apply( X = SB_pt, FUN = sum, MARGIN = 2, na.rm = T)



  # Get mortality-at-age
  M_apt <- repObj$M_apt[,,1:nT]
  Mbar_ap <- apply(X = M_apt, FUN = mean, MARGIN = c(1,2), na.rm = T)

  iscamB0   <- iscamRep$sbo
  siscaB0   <- sum(repObj$B0_p)

  iscamR0   <- iscamRep$ro
  siscaR0   <- sum(repObj$R0_p * exp(-Mbar_ap[2,]) )

  # Get age-2 numbers at age
  R_pt        <- array(NA, dim = c(nP,nT))
  R_pt[1:nP,] <- repObj$N_apt[2,1:nP,1:nT]

  R_pt[R_pt == 0] <- NA
  R_t   <- apply(X = R_pt, FUN = sum, MARGIN = 2, na.rm = T)


  # Get iscam SBt and age-2 recruitment
  iscamSB_t <- iscamRep$sbt
  iscamR_t  <- iscamRep$N[,1]

  # Scaled indices
  qSISCA_pg <- repObj$qhat_pg
  qISCAM_g  <- c(NA,NA,NA,iscamRep$q,NA)


  # reduce ISCAM to the subset that matches SISCA years
  siscaYrs <- fYearSISCA:lYearSISCA - fYearISCAM  + 1
  iscamSB_t <- iscamSB_t[siscaYrs]

  yrs <- fYearSISCA:lYearSISCA
  tdx <- yrs - fYearSISCA + 1

  # Pull indices
  siscaI_pgt  <- repObj$I_pgt
  siscaI_pgt[siscaI_pgt < 0] <- NA

  iscamIdata  <- iscamRep$d3_survey_data

  # Remove years of data that aren't relevant to this fit
  keepIdx <- which(iscamIdata[,1] %in% fYearSISCA:lYearSISCA)
  iscamIdata <- iscamIdata[keepIdx,]

  # Sum and remove zeroes
  scaledI_pgt <- siscaI_pgt
  if(!combIdx)
  {
    for( p in 1:nP )
    {
      scaledI_pgt[p,4,] <- siscaI_pgt[p,4,] / qSISCA_pg[p,4]
      scaledI_pgt[p,5,] <- siscaI_pgt[p,5,] / qSISCA_pg[p,5]
    }

  }

  scaledI_gt   <- apply( X = scaledI_pgt, FUN = sum, MARGIN = c(2,3), na.rm = T )
  scaledI_gt[ scaledI_gt == 0 ] <- NA

  if( combIdx )
  { 
    scaledI_pt <- datObj$combI_pt/repObj$qComb_pt
    scaledI_pt[scaledI_pt < 0] <- NA

    scaledI_t  <- apply( X = scaledI_pt, FUN = sum, MARGIN = 2, na.rm = T)
    scaledI_t[ scaledI_t == 0 ] <- NA

    scaledI_gt[scaledI_gt > 0] <- NA

  }

  iscamI_gt <- array(NA, dim = c(nG,nT))
  surveyFleets <- unique(iscamIdata[,3])

  iscamIdata[,1] <- iscamIdata[,1] - fYearISCAM + 1
  for( k in 1:nrow(iscamIdata))
    iscamI_gt[iscamIdata[k,3],iscamIdata[k,1]] <- iscamIdata[k,2]/qISCAM_g[iscamIdata[k,3]]

  iscamI_gt <- iscamI_gt[,siscaYrs]

  iscamC_gt <- array(NA, dim = c(nG,nT))
  iscamCtable <- iscamRep$dCatchData
  iscamCtable[,1] <- iscamCtable[,1] - fYearISCAM + 1
  for( k in 1:nrow(iscamCtable))
    iscamC_gt[iscamCtable[k,2], iscamCtable[k,1]] <- iscamCtable[k,7]

  # Pull catch and calculate HR
  iscamC_t <- apply( X = iscamC_gt, FUN = sum, MARGIN = 2 )

  stockCols <- brewer.pal(n = 3, "Set2")

  iscamResids_gt <- array(NA, dim = c(2,nT) )
  SISCAHresids_t <- rep(NA,nT)

  iscamResids_gt[1,] <- log(iscamI_gt[4,] / iscamSB_t[1:nT]) / iscamRep$sig
  iscamResids_gt[2,] <- log(iscamI_gt[5,] / iscamSB_t[1:nT]) / (iscamRep$sig / 1.166)

  SISCAHresids_t[1:nT] <- log(scaledI_t / SB_t[1:nT] )
  SISCAHresids_t[1:nT] <- SISCAHresids_t[1:nT] / sd(SISCAHresids_t[1:nT])

  # Now plot
  par(  mfrow = c(3,1), mar = c(0.1,1,.1,1),
        oma = c(3,4,1,1) )
  # First, SB
  plot( x = c(fYearSISCA,lYearSISCA),
        y = c(0, max(SB_t, iscamSB_t)), type = "n",
        axes = FALSE, xlab = "", ylab = "" )
    axis( side = 2, las = 1)
    box()
    grid()
    mtext( side = 2, text=  "Spawning Biomass (kt)", line = 3)

    # Plot SISCAH biomass
    if(simpleBt)
    {
      lines( x = yrs, y = SB_t[1:nT], col = "steelblue", lwd = 3)
    }
    if(!simpleBt)
    {
      bottom_t <- rep(0,nT)
      for( p in 1:nP )
      {
        rect( xleft = yrs - .3,
              xright = yrs + .3,
              ybottom = bottom_t,
              ytop = bottom_t + SB_pt[p,1:nT],
              col = stockCols[p], border = NA )

        bottom_t <- bottom_t + SB_pt[p,1:nT]
      }
    }

    lines( x = yrs, y = iscamSB_t,
            lwd = 3, col = "red" )

    if( simpleBt )
    {
      legText <- c( "ISCAM Bt",
                    "ISCAM B0",
                    "SISCAH Bt",
                    "SISCAH B0",
                    "Surface Survey Index",
                    "Dive Survey Index",
                    "Combined Spawn Index")

      legend( x = "topright", bty = "n",
              legend = legText,
              col = c("red","red","steelblue","steelblue","black","black","black"),
              pt.bg = c(NA,NA,NA,NA,NA,NA,"grey40"),
              lwd = c(3,3,3,3,NA,NA,NA),
              pch = c(NA,NA,NA,NA,4,5,21),
              lty = c(1,2,1,3,NA,NA,NA),
              pt.cex = 1.5)

      abline( h = c(iscamB0,siscaB0),
            lty = c(2,3), lwd = 2, col = c("red","steelblue") )
    } else {

      legText <- c( "ISCAM Bt",
                    "ISCAM B0",
                    "C/S Bt",
                    "JP/S Bt",
                    "Lou Bt",
                    "SISCAH B0",
                    "Surface Survey Index",
                    "Dive Survey Index",
                    "Combined Spawn Index")

      legend( x = "topright", bty = "n",
              legend = legText,
              col = c("red","red",stockCols,"grey40","black","black","black"),
              pt.bg = c(NA,NA,stockCols,NA,NA,NA,"grey40"),
              lwd = c(3,3,NA,NA,NA,3,NA,NA,NA),
              pch = c(NA,NA,22,22,22,NA,4,5,21),
              lty = c(1,2,NA,NA,NA,3,NA,NA,NA),
              pt.cex = 1.5)

      abline( h = c(iscamB0,siscaB0),
            lty = c(2,3), lwd = 2, col = c("red","grey40") )
    }
    
    for( g in surveyFleets)
    {
      points( x = yrs, y = iscamI_gt[g,],
              col = "black", pch = g, cex = 1.5 )
    }
    if( combIdx )
      points( x = yrs, y = scaledI_t,
              bg = "grey60", pch = 21, col = "black", cex = 1.5)

  plot( x = c(fYearSISCA,lYearSISCA),
        y = c(-3,4), type = "n",
        axes = FALSE, xlab = "", ylab = "" )
    box()
    grid()
    axis( side = 2, las = 1 )
    mtext( side = 2, text = "Std. log-residuals", line = 3)
    abline( h = 0, lty = 2, lwd = 1)
    points( x = yrs, y = iscamResids_gt[1,], pch = 4, cex = 1.5 )
    points( x = yrs, y = iscamResids_gt[2,], pch = 5, cex = 1.5 )
    points( x = yrs, y = SISCAHresids_t, pch = 21, col = "black", bg = "grey40", cex = 1.5 )

    meanISCAMresids <- apply( X = iscamResids_gt, FUN = mean, MARGIN = 1, na.rm = T)

    legText <- c( paste("Surface Survey, mean ", round(meanISCAMresids[1],2), sep = "" ),
                  paste("Dive Survey, mean ", round(meanISCAMresids[2],2), sep = "" ),
                  paste("Combined Spawn Index, mean ", round(mean(SISCAHresids_t),2), sep = "" ) )
    legend( x = "topleft", bty = "n",
            col = c("black"),
            pch = c(4,5,21),
            pt.bg = c(NA,NA,"grey40"),
            cex = 1.5, legend = legText )


  plot( x = c(fYearSISCA,lYearSISCA),
        y = c(0, max(R_t, iscamR_t)), type = "n",
        axes = FALSE, xlab = "", ylab = "" )
    axis( side = 2, las = 1)
    axis( side = 1 )
    mtext( side = 2, text=  "Age-2 Recruitment (1e6)", line = 3)
    mtext( side = 1, text=  "Year", line = 2)
    box()
    grid()
    legend( x = "topright", bty = "n",
            legend = c( "SISCAH Rt", "SISCAH R0", 
                        "ISCAM Rt", "ISCAM R0"),
            col = c("steelblue","steelblue","red","red"),
            lwd = c(3,3,3,3),
            lty = c(1,3,1,2))
    lines( x = yrs, y = R_t[tdx],
            lwd = 3, col = "steelblue" )
    lines( x = yrs, y = iscamR_t[tdx],
            lwd = 3, col = "red" )
    abline( h = c(iscamR0,siscaR0),
            lty = c(2,3), lwd = 2, col = c("red","steelblue") )

} # END plotAggBtRtComparison()

# plotPerCapSP
# Per Capita surplus production plotted against
# biomass depletion
plotPerCapSP <- function( repList = reports,
                          minAge = NULL,
                          maxY = 2,
                          minYear = 1951  )
{
  repObj <- repList$repOpt

  # Pull biomass and catch
  B_apt   <- repObj$B_apt
  C_pgt   <- repObj$C_pgt
  nA      <- repObj$nA
  nP      <- repObj$nP
  nT      <- repObj$nT
  initT_p <- repObj$tInitModel_p
  B0_p    <- repObj$B0_p
  R0_p    <- repObj$R0_p

  fYear <- repList$fYear
  tInit <- minYear - fYear + 1


  if( !is.null(minAge))
  {
    B0minAge_p <- B0_p

    for( p in 1:nP )
    {
      eqbmN_a <- R0_p[p] * repObj$initSurv_ap[,p]
      B0minAge_p[p] <- sum(eqbmN_a[minAge:nA] * repObj$meanWt_ap[minAge:nA,p])
    }

    B0_p <- B0minAge_p
  }

  # Get biomass and catch
  B_pt    <- repObj$SB_pt[,1:nT]
  C_pt    <- apply( X = C_pgt, FUN = sum, MARGIN = c(1,3) )

  if( !is.null(minAge))
    B_pt    <- apply( X = B_apt[minAge:nA,,1:nT], FUN = sum, MARGIN = c(2,3) )
  
  # Calculate surplus production
  P_pt     <- array(NA, dim = c(nP,nT))

  for( p in 1:nP)
    for( t in (initT_p[p]+1):(nT-1) )
      P_pt[p,t] <- B_pt[p,t+1] - B_pt[p,t] + C_pt[p,t]



  # Now do the aggregate
  # B_t <- apply( X = B_apt[minAge:nA,,1:nT], FUN = sum, MARGIN = c(3))
  B_t <- apply( X = B_pt[,1:nT], FUN = sum, MARGIN = c(2))
  C_t <- apply( X = C_pgt, FUN = sum, MARGIN = c(3))

  P_t   <- array(NA, dim = c(nT))
  P_t[1:(nT-1)] <- B_t[2:nT] - B_t[1:(nT - 1)] + C_t[1:(nT-1)]



  ptCols <- scales::viridis_pal(option = "C",alpha = 1)(nT - tInit+1)

  names(ptCols) <- minYear:2019

  stockNames <- repList$stock

  pcP_pt <- P_pt/B_pt
  pcP_t  <- P_t/B_t

  depB_pt <- B_pt
  for( p in 1:nP )
    depB_pt[p,] <- B_pt[p,]/B0_p[p]
  depB_t  <- B_t/sum(B0_p)

  if( is.null(maxY) )
    maxY <- max(pcP_pt,pcP_t,na.rm =T)

  minY <- min(pcP_pt,pcP_t,na.rm =T)

  par(  mfrow = c(2,2),
        mar = c(1,.1,1,.1),
        oma = c(4,4,1,1) )
  for( p in 1:nP )
  {
    plot( x = c(0,max(B_pt[p,tInit:nT],na.rm = T)),
          y = c(minY,maxY),
          axes = FALSE, type = "n", las = 1 )
      mfg <- par("mfg")
      # if( mfg[1] == mfg[3])
        axis(side = 1)
      if( mfg[2] == 1)
        axis(side = 2, las = 1)
      box()
      abline(h = 0, lty = 2 )
      
      arrows( x0 = B_pt[p,tInit:(nT-1)],
              x1 = B_pt[p,(tInit+1):(nT)],
              y0 = P_pt[p,tInit:(nT-1)]/B_pt[p,tInit:(nT-1)],
              y1 = P_pt[p,(tInit+1):(nT)]/B_pt[p,(tInit+1):(nT)],
              col = "grey75", lwd = .8, length = .1 )
      points( x = B_pt[p,tInit:nT],
              y = P_pt[p,tInit:nT]/B_pt[p,tInit:nT],
              bg = ptCols, pch = 21 )
      abline( v = .3*B0_p[p], lty = 3, lwd = .8)
      abline( v = B0_p[p], lty = 3, lwd = .8, col = "red")

      legend( x = "topleft",
              legend = stockNames[p], bty = "n")

  }

  plot( x = c(0,max(B_t[tInit:nT],na.rm = T)),
          y = c(minY,maxY),
          axes = FALSE, type = "n", las = 1 )
    mfg <- par("mfg")
    # if( mfg[1] == mfg[3])
      axis(side = 1)
    if( mfg[2] == 1)
      axis(side = 2, las = 1)
    box()
    abline(h = 0, lty = 2 )
    arrows( x0 = B_t[tInit:(nT-1)],
            x1 = B_t[(tInit+1):(nT)],
            y0 = P_t[tInit:(nT-1)]/B_t[tInit:(nT-1)],
            y1 = P_t[(tInit+1):(nT)]/B_t[(tInit+1):(nT)],
            col = "grey75", lwd = .8, length = .1 )
    points( x = B_t[tInit:nT],
            y = P_t[tInit:nT]/B_t[tInit:nT],
            bg = ptCols, pch = 21 )

    abline( v = .3 * sum(B0_p), lty = 3, lwd = .8)
    abline( v = sum(B0_p), lty = 3, lwd = .8, col = "red")

    legend( x = "topright",
            bty = "n",
            legend = c( minYear,
                        "2019"),
            pch = 21, pt.bg = ptCols[c(1,nT-tInit+1)] )
    legend( x = "topleft",
            legend = "Aggregate",
            bty = "n")

    xLab <- "Spawning Biomass (kt)"

    if( !is.null(minAge))
      xLab <- paste( minAge,"+ Biomass (kt)",sep = "")

    mtext(side = 1, outer = TRUE, text = xLab,
          line = 2.5)
    mtext(side = 2, outer = TRUE, text = "Per Capita Surplus Production",
          line = 2.5)
} # END plotPerCapSP()

# plotRecPerSpawner()
# Per capita recruits per spawner, as a function
# of biomass
plotRecPerSpawner <- function(  repList = reports,
                                maxY = NULL  )
{
  repObj <- repList$repOpt

  # Pull biomass and catch
  SB_pt   <- repObj$SB_pt
  R_pt    <- repObj$R_pt
  C_pgt   <- repObj$C_pgt
  nA      <- repObj$nA
  nP      <- repObj$nP
  nT      <- repObj$nT
  initT_p <- repObj$tInitModel_p
  B0_p    <- repObj$B0_p

  B_pt    <- repObj$SB_pt[,1:nT]

  # Calculate surplus production
  RPS_pt     <- array(NA, dim = c(nP,nT))

  for( p in 1:nP)
    for( t in (initT_p[p]+1):(nT-1) )
      RPS_pt[p,t] <- R_pt[p,t+1]/B_pt[p,t]



  # Now do the aggregate
  # B_t <- apply( X = B_apt[minAge:nA,,1:nT], FUN = sum, MARGIN = c(3))
  B_t <- apply( X = B_pt[,1:nT], FUN = sum, MARGIN = c(2))
  R_t <- apply( X = R_pt, FUN = sum, MARGIN = c(2))

  RPS_t   <- array(NA, dim = c(nT))
  RPS_t[1:(nT-1)] <- R_t[2:nT]/ B_t[1:(nT - 1)]

  DecadeCols <- scales::viridis_pal(option = "C",alpha = 1)(7)

  ptCols <- c(  rep(DecadeCols[1],10),
                rep(DecadeCols[2],10),
                rep(DecadeCols[3],10),
                rep(DecadeCols[4],10),
                rep(DecadeCols[5],10),
                rep(DecadeCols[6],10),
                rep(DecadeCols[7],9) )

  names(ptCols) <- 1951:2019

  stockNames <- repList$stock

  if( is.null(maxY) )
    maxY <- max(RPS_pt,RPS_t,na.rm =T)

  minY <- 0

  par(  mfrow = c(2,2),
        mar = c(1,1.5,1,1),
        oma = c(4,3,1,1) )
  for( p in 1:nP )
  {
    plot( x = c(0,max(B_pt[p,],na.rm = T)),
          y = c(0,max(RPS_pt[p,],na.rm = T)),
          axes = FALSE, type = "n", las = 1 )
      mfg <- par("mfg")
      # if( mfg[1] == mfg[3])
        axis(side = 1)
      # if( mfg[2] == 1)
        axis(side = 2, las = 1)
      box()
      points( x = B_pt[p,],
              y = RPS_pt[p,],
              bg = ptCols, pch = 21 )
      abline( v = .3*B0_p[p], lty = 3, lwd = .8)
      abline( v = B0_p[p], lty = 3, lwd = .8, col = "red")

      legend( x = "topleft",
              legend = stockNames[p], bty = "n")

  }

  plot( x = c(0,max(B_t,na.rm = T)),
          y = c(0,max(RPS_t,na.rm = T)),
          axes = FALSE, type = "n", las = 1 )
    mfg <- par("mfg")
    # if( mfg[1] == mfg[3])
      axis(side = 1)
    # if( mfg[2] == 1)
      axis(side = 2, las = 1)
    box()
    points( x = B_t,
            y = RPS_t,
            bg = ptCols, pch = 21 )

    abline( v = .3 * sum(B0_p), lty = 3, lwd = .8)
    abline( v = sum(B0_p), lty = 3, lwd = .8, col = "red")

    legend( x = "topright",
            bty = "n",
            legend = c( "1951 - 1960",
                        "1961 - 1970",
                        "1971 - 1980",
                        "1981 - 1990",
                        "1991 - 2000",
                        "2001 - 2010",
                        "2011 - 2019"),
            pch = 21, pt.bg = DecadeCols )
    legend( x = "topleft",
            legend = "Aggregate",
            bty = "n")

    mtext(side = 1, outer = TRUE, text = "Spawning Biomass (kt)",
          line = 1.5)
    mtext(side = 2, outer = TRUE, text = "Recruits per Spawning Stock Biomass (1e6/kt)",
          line = 1.5)
} # END plotRecPerSpawner()

# plotMvsB()
# Plots M as a function of beginning year biomass in each 
# year
plotMvsB <- function( repList=reports,
                      minAge = 3 )
{
  repObj <- repList$repOpt

  # Pull biomass and catch
  B_apt   <- repObj$B_apt
  M_pt    <- repObj$M_apt[minAge,,]
  nA      <- repObj$nA
  nP      <- repObj$nP
  nT      <- repObj$nT
  initT_p <- repObj$tInitModel_p
  B0_p    <- repObj$B0_p
  R0_p    <- repObj$R0_p

  stockNames <- repList$stock

  B_pt    <- apply( X = B_apt[minAge:nA,,], FUN = sum, MARGIN = c(2,3) )

  B_t     <- apply( X = B_pt, FUN = sum, MARGIN = 2 )

  B0minAge_p <- B0_p

  for( p in 1:nP )
  {
    eqbmN_a <- R0_p[p] * repObj$initSurv_ap[,p]
    B0minAge_p[p] <- sum(eqbmN_a[minAge:nA] * repObj$meanWt_ap[minAge:nA,p])
  }

  B0_p <- B0minAge_p

  M_t     <- numeric( length = nT )
  for( t in 1:nT )
    M_t[t] <- sum(M_pt[,t] * B0_p)/sum(B0_p)

  depB_pt <- B_pt
  for( p in 1:nP )
    depB_pt[p,] <- depB_pt[p,] / B0_p[p]
  
  depB_t  <- B_t / sum(B0_p)

  depB_pt[depB_pt == 0] <- NA
  M_pt[M_pt == 0] <- NA

  stockCols <- brewer.pal(n = 4, "Set2")

  plot( x = c(0,max(depB_pt,depB_t,na.rm = T)),
        y = c(0,max(M_pt,M_t,na.rm = T)),
        axes = FALSE, type = "n", las = 1,
        xlab = "", ylab = "" )
    mfg <- par("mfg")
    axis(side = 1)
    axis(side = 2, las = 1)
    box()
    for(p in 1:nP )
    {
      points( x = depB_pt[p,],
              y = M_pt[p,],
              bg = alpha(stockCols[p],alpha = .4), pch = 21, col = NA )
      lines(loess.smooth( x = depB_pt[p,],
                          y = M_pt[p,], span = 1),
            col = stockCols[p], lwd = 2)
    }
    points( x = depB_t[1:nT],
              y = M_t[1:nT],
              bg = alpha(stockCols[4],alpha = .4), pch = 21, col = NA )
    lines(  loess.smooth( x = depB_t[1:nT],
                          y = M_t[1:nT], span = 1),
            col = stockCols[4], lwd = 2)


    abline( v = .3, lty = 3, lwd = .8)
    abline( v = 1, lty = 3, lwd = .8, col = "red")


    legend( x = "topright",
            bty = "n",
            legend = c(stockNames,"Aggregate"),
            pch = 21, pt.bg = stockCols, pt.lwd = 0,
            lty = 1, col = stockCols, lwd = 2 )


    mtext(side = 1, outer = FALSE, text = "3+ Biomass Depletion",
          line = 2)
    mtext(side = 2, outer = FALSE, text = "Natural Mortality (yr)",
          line = 2.5)

}# END plotMvsB()



# plotStockDataSummary()
# A visual summary of the data
# available for a given stock
plotStockDataSummary <- function( dataList,
                                  fYear = 1951 )
{
  # Pull data
  I_pgt     <- dataList$I_pgt
  A_apgt    <- dataList$A_apgt
  C_pgt     <- dataList$C_pgt
  combI_pt  <- dataList$combI_pt
  # Mixed catch
  mC_gt  <- dataList$mC_gt

  # First, get dimension
  nT <- dim(I_pgt)[3]
  nP <- dim(I_pgt)[1]
  yrs <- seq( from = fYear, by = 1, length.out = nT+4)

  # Replace -1s with NA
  A_apgt[A_apgt < 0] <- NA
  I_pgt[I_pgt < 0] <- 0
  combI_pt[combI_pt < 0] <- NA

  gearCols <- RColorBrewer::brewer.pal(6, "Dark2")

  A_pgt <- apply( X = A_apgt, FUN = sum, MARGIN = c(2,3,4),
                  na.rm = TRUE )
  A_gt <- apply( X = A_pgt, FUN = sum, MARGIN = c(2,3),
                  na.rm = TRUE )
  sampSize_pg <- apply( X = A_pgt, FUN = sum, MARGIN = c(1,2) )
  sampSize_g  <- apply( X = A_gt, FUN = sum, MARGIN = c(1) )
  
  pchA_pgt <- A_pgt
  pchA_pgt[A_pgt < 200 ] <- 1
  pchA_pgt[A_pgt >= 200 ] <- 16
  pchA_pgt[A_pgt == 0 ] <- NA

  pchA_gt <- A_gt
  pchA_gt[A_gt < 200 ] <- 1
  pchA_gt[A_gt >= 200 ] <- 16
  pchA_gt[A_gt == 0 ] <- NA

  pchI_pgt <- I_pgt
  pchI_pgt[ I_pgt > 0 ] <- 16
  pchI_pgt[ I_pgt == 0 ] <- NA

  pchC_pgt <- C_pgt
  pchC_pgt[ C_pgt > 0 ] <- 16
  pchC_pgt[ C_pgt == 0 ] <- NA

  pchCombI_pt <- combI_pt
  pchCombI_pt[ combI_pt > 0] <- 16
  pchCombI_pt[ combI_pt == 0] <- 1


  I_gt <- apply( X = I_pgt, FUN = sum, MARGIN = c(2,3))
  pchI_gt <- I_gt
  pchI_gt[ I_gt > 0 ] <- 16
  pchI_gt[ I_gt == 0 ] <- NA  

  C_gt <- apply( X = C_pgt, FUN = sum, MARGIN = c(2,3))
  C_gt <- C_gt + mC_gt
  pchC_gt <- C_gt
  pchC_gt[ C_gt > 0 ] <- 16
  pchC_gt[ C_gt == 0 ] <- NA 

  combI_t <- apply( X = combI_pt, FUN = sum, MARGIN = c(2), na.rm = T)
  pchCombI_t <- combI_t
  pchCombI_t[combI_t > 0] <- 16


  ageGears <- 1:3
  catGears <- c(1:3,6)
  idxGears <- 4:5

  # Make a plotting area

  stockLab <- c("C/S","JP/S","Lou")

  yLabs <- c("Ages", "Catch", "Spawn Idx" )

  gearLabs <- c("Reduction","SeineRoe","Gillnet","Surface","Dive","SOK")

  if( any(!is.na(combI_pt)))
    gearLabs <- c("Reduction","SeineRoe","Gillnet","Spawn Idx","SOK")

  par(  mfrow = c(nP+1,1), mar = c(0.5,1,.5,1),
        oma = c(7,6,3,4), xpd = FALSE )
  # Now plot aggregate/mixed data
  plot( x = range(yrs), y = c(0.5,3.5),
          type = "n", axes = FALSE )
    axis( side = 2, las = 1, at = 1:3,
          labels = yLabs)
    mtext( side = 3, text = "SISCA Data Summary", font = 2)

    box()
    grid()  
    abline(h = c(1:2) + .5, lwd = 2)
    rmtext( txt = "Aggregate", cex = 1.5,
            line = .05, outer = TRUE, font = 2 )
    # Plot age comps
    points( x = yrs[1:nT], y = rep(1,nT)-.2,
            pch = pchA_gt[1,],
            col = gearCols[1] )
    points( x = yrs[1:nT], y = rep(1,nT),
            pch = pchA_gt[2,],
            col = gearCols[2] )
    points( x = yrs[1:nT], y = rep(1,nT)+.2,
            pch = pchA_gt[3,],
            col = gearCols[3] )
    # Plot Catch
    points( x = yrs[1:nT], y = rep(2,nT)-.3,
            pch = pchC_gt[1,],
            col = gearCols[1] )
    points( x = yrs[1:nT], y = rep(2,nT)-.1,
            pch = pchC_gt[2,],
            col = gearCols[2] )
    points( x = yrs[1:nT], y = rep(2,nT)+.1,
            pch = pchC_gt[3,],
            col = gearCols[3] )
    points( x = yrs[1:nT], y = rep(2,nT)+.3,
            pch = pchC_gt[6,],
            col = gearCols[6] )


    # Plot Indices
    points( x = yrs[1:nT], y = rep(3,nT)-.2,
            pch = pchI_gt[4,],
            col = gearCols[4] )
    points( x = yrs[1:nT], y = rep(3,nT)+.2,
            pch = pchI_gt[5,],
            col = gearCols[5] )
    points( x = yrs[1:nT], y = rep(3,nT),
            pch = pchCombI_t,
            col = "grey40" )



    text( x = yrs[nT]+4, y = 1 + c(-.2,0,.2), 
          labels = sampSize_g[1:3], cex = .8 )
  for( p in 1:nP )
  {
    plot( x = range(yrs), y = c(0.5,3.5),
          type = "n", axes = FALSE )

      axis( side = 2, las = 1, at = 1:3,
            labels = yLabs)
      box()
      grid()
      rmtext( txt = stockLab[p], cex = 1.5,
              font = 2, line = .05, outer = TRUE )

      # Plot age comps
      points( x = yrs[1:nT], y = rep(1,nT)-.2,
              pch = pchA_pgt[p,1,],
              col = gearCols[1] )
      points( x = yrs[1:nT], y = rep(1,nT),
              pch = pchA_pgt[p,2,],
              col = gearCols[2] )
      points( x = yrs[1:nT], y = rep(1,nT) + .2,
              pch = pchA_pgt[p,3,],
              col = gearCols[3] )
      # Plot catch
      points( x = yrs[1:nT], y = rep(2,nT)-.3,
              pch = pchC_pgt[p,1,],
              col = gearCols[1] )
      points( x = yrs[1:nT], y = rep(2,nT)-.1,
              pch = pchC_pgt[p,2,],
              col = gearCols[2] )
      points( x = yrs[1:nT], y = rep(2,nT)+.1,
              pch = pchC_pgt[p,3,],
              col = gearCols[3] )
      points( x = yrs[1:nT], y = rep(2,nT)+.3,
              pch = pchC_pgt[p,6,],
              col = gearCols[6] )

      # Plot indices
      points( x = yrs[1:nT], y = rep(3,nT)-.2,
              pch = pchI_pgt[p,4,],
              col = gearCols[4] )
      points( x = yrs[1:nT], y = rep(3,nT)+.2,
              pch = pchI_pgt[p,5,],
              col = gearCols[5] )
      points( x = yrs[1:nT], y = rep(3,nT),
              pch = pchCombI_pt[p,],
              col = "grey40" )
      abline(h = c(1:2)+.5, lwd = 2)

  
      text( x = yrs[nT]+4, y = 1 + c(-.2,0,.2), 
            labels = sampSize_pg[p,1:3], cex = .8 )

  }
    axis( side = 1)
    mtext( side = 1, text = "Year", font = 2, outer = TRUE,
            line = 2)

    if( any(!is.na(combI_pt)))
      gearCols <- c(gearCols[1:3],"grey40",gearCols[6])

    par(xpd=TRUE, oma = c(0,0,0,0))
    # plot( x = c(0,1), y = c(0,1), add = TRUE)
    legend( x = "bottom", bty = "n",
            # inset = c(0,-1), xpd = TRUE,
            legend = gearLabs, 
            pch = 16, col = gearCols,
            horiz = TRUE, cex = 1.5, pt.cex = 2)
}


plotSimilarityMetrics <- function(  sensTable, 
                                    winX = c(.3,.6),
                                    winY = c(0,.32) )
{
  AREB0range      <-  c(0,max(sensTable$AREB0))
  AREstatusrange  <-  c(0,max(sensTable$AREstatus))
  MAREBtrange     <-  c(0,max(sensTable$MAREBt))

  nSens <- length(unique(sensTable$Sensitivity))

  cols <- scales::viridis_pal(option = "C",alpha = .8)(nSens)


  sensNames <-  sensTable %>%
                group_by(Sensitivity) %>%
                summarise( nModels = n() ) %>%
                ungroup() %>%
                mutate( sensCol = cols )

  sensTable <-  sensTable %>%
                left_join( sensNames, by = "Sensitivity" ) %>%
                mutate( cex = .8 + 4*( AREB0 - min(AREB0) )/( max(AREB0) - min(AREB0)) ,
                        pch = ifelse(pdHess, 21, 22) )

  simMetrics <- as.matrix(sensTable[,c("AREstatus","AREB0","MAREBt")])
  simMetricCorrelation <- cor(simMetrics)


  par( mfrow = c(2,1), mar = c(1,2,1,2), oma =c(3,3,1,1) )
  
  plot( x = MAREBtrange, y = AREstatusrange,
        type = "n", las = 1 )
    
    points( x = sensTable$MAREBt,
            y = sensTable$AREstatus,
            cex = sensTable$cex,
            bg = sensTable$sensCol,
            pch = sensTable$pch,
            col = "black" )
    points( x = baseCase$MAREBt, y = baseCase$AREstatus,
            col = "black", pch = 2, cex = 1.5)
    # Make a window for the zoomed in plot
    segments( x0 = winX,
              y0 = rep(winY[1],2), y1 = rep(winY[2],2),
              col = "black", lwd = 1)

    segments( x0 = rep(winX[1],2), x1 = rep(winX[2],2),
              y0 = winY,
              col = "black", lwd = 1)
                
    plot( x = winX, y = winY,
          type = "n", las = 1 )
      # mtext( side = 2, text = "ARE(B2019/B0)")
      points( x = sensTable$MAREBt,
              y = sensTable$AREstatus,
              cex = sensTable$cex,
              bg = sensTable$sensCol,
              pch = sensTable$pch,
              col = "black" )
      points( x = baseCase$MAREBt, y = baseCase$AREstatus,
              col = "black", pch = 2, cex = 1.5)
      legend( x = "topright", bty = "n",
            legend = c(sensNames$Sensitivity[-nSens],"Base Model"),
            pt.bg = c(sensNames$sensCol[-nSens],NA),
            pch = c(rep(21,nSens-1),2),
            col = c(rep("black",nSens-1),"black"),
            pt.cex = 1.5 )
  
  mtext( side = 2, text = expression(rho(B[2019]/B[0])), outer = TRUE,
          line = 1.5) 
  mtext( side = 1, text = expression(rho(B[t])),
          outer = TRUE, line = 2 )
}

plotSimMetricCorrelation <- function(  sensTable )
{
  require(corrplot)
  simMetrics <- as.matrix(sensTable[,c("AREstatus","AREB0","MAREBt")])
  colnames(simMetrics) <- c("rho(Status)","rho(B0)","rho(Bt)")
  # colnames(simMetrics) <- c(  expression( rho(B[2019]/B[0])) ,
  #                             expression( rho(B[0]) ),
  #                             expression( rho(B[t]) ) )
  simMetricCorrelation <- cor(simMetrics)

  Dmetrics <- sensTable %>%
              dplyr::select(  Equal,
                              Status = StatusDom,
                              B0 = B0dom,
                              Bt = BtDom,
                              ageSE = meanAgeSE,
                              idxSE = meanIdxSE ) %>%
              as.matrix()


  DmetricCorrelation <- cor(Dmetrics)

  par( mfrow = c(2,1), mar = c(0,0,0,0), oma = c(0.,0.,0,0) )

  corrplot.mixed( simMetricCorrelation, lower.col = "black", number.cex = 1 )
  corrplot.mixed( DmetricCorrelation, lower.col = "black", number.cex = 1 )
}

# plotSpawnIdxSP()
# Plot of spawn indices and surplus production,
# reproducing Kronlund et al Figure 22 - 29
plotSpawnIdxSP <- function( repList = reports,
                            stock = "Agg",
                            lowSSBquant = .3 )
{
  # First, make SP tables
  SPtables <- makeSPtables(repList, save = FALSE)

  # Pull table for stock
  spTable <- SPtables[[stock]]

  repObj <- repList$repOpt

  # Pull biomass and catch
  nA      <- repObj$nA
  nP      <- repObj$nP
  nT      <- repObj$nT

  fYear   <- repList$fYear
  yrs     <- spTable$Year

  SSBquant  <- quantile(spTable$SSB, probs = lowSSBquant, na.rm = T)
  spTable <- spTable %>% mutate( ptCol = ifelse( SSB <= SSBquant, "grey40","white") )


  par(mfrow = c(3,1), oma = c(3,4,2,2), mar = c(.1,.1,.1,.1))

  # Plot spawn indices
  plot( x = range(yrs), y = range(spTable$spawnIdx, na.rm = T),
        type = "n", axes = FALSE )
    axis( side = 2, las = 1 )
    grid()
    box()
    mtext( side = 2, line = 2.5, text = "Spawn Index (kt)" )
    lines( x = yrs, y = spTable$spawnIdx, lwd = 2 )
    points( x = yrs, y = spTable$spawnIdx, pch = 21,
            bg = spTable$ptCol, cex = 1.5 )
    lines( x = yrs[-(1:2)], y = pracma::movavg(x = spTable$spawnIdx, n = 3)[-(1:2)],
            col = "steelblue", lwd = 3 )
    abline(v = 1987.5, lty = 2, lwd = 1)
    legend( x = "topright",
            legend = c( paste("SSB < ", 100*lowSSBquant, "th percentile", sep = ""),
                        "3 Year Moving Avg"),
            pch = c(21,NA),
            col = c("black","steelblue"),
            pt.bg = c("grey40",NA),
            lty = c(NA,1),
            lwd = c(NA,3) )

  plot( x = range(yrs), y = range(spTable$SurplusProd, na.rm = T),
        type = "n", axes = FALSE )
    axis( side = 2, las = 1 )
    grid()
    box()
    mtext( side = 2, line = 2.5, text = "Surplus Production (kt)" )
    abline( h = 0, lty = 2, lwd = .8)
    lines( x = yrs, y = spTable$SurplusProd, lwd = 2 )
    points( x = yrs, y = spTable$SurplusProd, pch = 21,
            bg = spTable$ptCol, cex = 1.5 )
    lines( x = yrs[-(1:2)], y = pracma::movavg(x = spTable$SurplusProd, n = 3)[-(1:2)],
            col = "steelblue", lwd = 3 )
    abline(v = 1987.5, lty = 2, lwd = 1)

  plot( x = range(yrs), y = range(spTable$perCapSurplusProd, na.rm = T),
        type = "n", axes = FALSE )
    axis( side = 2, las = 1 )
    axis( side = 1)
    grid()
    box()
    abline( h = 0, lty = 2, lwd = .8)
    mtext( side = 2, line = 2.5, text = "Surplus Production Rate" )
    lines( x = yrs, y = spTable$perCapSurplusProd, lwd = 2 )
    points( x = yrs, y = spTable$perCapSurplusProd, pch = 21,
            bg = spTable$ptCol, cex = 1.5 )
    lines(  x = yrs[-(1:2)], 
            y = pracma::movavg(x = spTable$perCapSurplusProd, n = 3)[-(1:2)],
            col = "steelblue", lwd = 3 )
    abline(v = 1987.5, lty = 2, lwd = 1)


} # END plotSpawnIdxSP

# plotSPphasePlots()
# Plots SP phase plots similar to Kronlund et al
plotSPphasePlots <- function( repList = reports,
                              stock = "Agg",
                              depLines = c(.1,.25,.3,.5,.6,1.0),
                              depCols  = c("red","grey30") )
{
  # First, make SP tables
  SPtables <- makeSPtables(repList, save = FALSE)

  # Pull table for stock
  spTable <- SPtables[[stock]]

  repObj <- repList$repOpt

  # Pull biomass and catch
  nA      <- repObj$nA
  nP      <- repObj$nP
  nT      <- repObj$nT

  if( stock == "Agg" )
    B0 <- sum(repObj$B0_p)
  else
    B0 <- repObj$B0_p[which(names(SPtables == stock))]

  fYear   <- repList$fYear
  yrs     <- spTable$Year

  surfYrs <- 1951:1987
  diveYrs <- 1988:2019

  nSurf <- length(surfYrs)
  nDive <- length(diveYrs)

  surfTable <- spTable %>% filter(Year %in% surfYrs)
  diveTable <- spTable %>% filter(Year %in% diveYrs)

  depLinePalette  <- colorRampPalette(depCols)
  depLineCols     <- depLinePalette(length(depLines)) 

  ptCols <- scales::viridis_pal(option = "C",alpha = 1)(nT)
  names(ptCols) <- fYear:2019

  # 2x2 plot, 
  # left column is 1951-1987, right is 1988+
  # top row is abs SP, bottom is per cap
  par(  mfrow = c(2,2),
        mar = c(.1,1.5,.1,1.5),
        oma = c(4,2,1,1) )
  
  plot( x = c(0,max(surfTable$SSB)),
        y = range(surfTable$SurplusProd),
        axes = FALSE, type = "n", las = 1 )
    axis(side = 2, las = 1)
    box()
    abline(h = 0, lty = 2 )
    abline(v = depLines * B0, col = depLineCols, lty = 2 )
    mtext( side = 2, text = "Surplus Production (kt)", line = 2)
    # Put in arrows first

    arrows( x0 = surfTable$SSB[1:(nSurf-1)],
            x1 = surfTable$SSB[2:(nSurf)],
            y0 = surfTable$SurplusProd[1:(nSurf-1)],
            y1 = surfTable$SurplusProd[2:(nSurf)],
            col = "grey75", lwd = .8, length = .1 )
    points( x = surfTable$SSB,
            y = surfTable$SurplusProd,
            bg = ptCols[1:nSurf], pch = 21 )
    text( x = surfTable$SSB[c(1,nSurf)] ,
          y = surfTable$SurplusProd[c(1,nSurf)]+ c(-5,5),
          labels = surfTable$Year[c(1,nSurf)] )

    legend( x = "topright",
            bty = "n",
            legend = c( "1951",
                        "1987"),
            pch = 21, pt.bg = ptCols[c(1,nSurf)] )
    mtext( side = 3, text="1951-1987", font = 2)

  plot( x = c(0,max(diveTable$SSB)),
        y = range(diveTable$SurplusProd, na.rm = T),
        axes = FALSE, type = "n", las = 1 )
    axis(side = 2, las = 1)
    box()
    abline(h = 0, lty = 2 )
    abline(v = depLines * B0, col = depLineCols, lty = 2 )
    # Put in arrows first

    arrows( x0 = diveTable$SSB[1:(nDive-1)],
            x1 = diveTable$SSB[2:(nDive)],
            y0 = diveTable$SurplusProd[1:(nDive-1)],
            y1 = diveTable$SurplusProd[2:(nDive)],
            col = "grey75", lwd = .8, length = .1 )
    points( x = diveTable$SSB,
            y = diveTable$SurplusProd,
            bg = ptCols[(nSurf+1):nT], pch = 21 )
    text( x = diveTable$SSB[c(1,nDive-1)] + c(-1,0),
          y = diveTable$SurplusProd[c(1,nDive-1)] + c(0,0),
          labels = diveTable$Year[c(1,nDive-1)] )

    legend( x = "topright",
            bty = "n",
            legend = c( "1988",
                        "2018"),
            pch = 21, pt.bg = ptCols[c(nSurf+1,nT-1)] )
    mtext( side = 3, text="1988-2019", font = 2)


  plot( x = c(0,max(surfTable$SSB)),
        y = range(surfTable$perCapSurplusProd, na.rm = T),
        axes = FALSE, type = "n", las = 1 )
    axis(side = 1)
    axis(side = 2, las = 1)
    box()
    abline(h = 0, lty = 2 )
    abline(v = depLines * B0, col = depLineCols, lty = 2 )
    mtext( side = 2, text = "Surplus Production Rate", line = 2)
    # Put in arrows first

    arrows( x0 = surfTable$SSB[1:(nSurf-1)],
            x1 = surfTable$SSB[2:(nSurf)],
            y0 = surfTable$perCapSurplusProd[1:(nSurf-1)],
            y1 = surfTable$perCapSurplusProd[2:(nSurf)],
            col = "grey75", lwd = .8, length = .1 )
    points( x = surfTable$SSB,
            y = surfTable$perCapSurplusProd,
            bg = ptCols[1:nSurf], pch = 21 )
    text( x = surfTable$SSB[c(1,nSurf)],
          y = surfTable$perCapSurplusProd[c(1,nSurf)] + c(.1,.1),
          labels = surfTable$Year[c(1,nSurf)] )

    legend( x = "topright",
            bty = "n",
            legend = paste(depLines,"B0",sep = ""),
            lty = 2, col = depLineCols )

  plot( x = c(0,max(diveTable$SSB)),
        y = range(diveTable$perCapSurplusProd, na.rm = T),
        axes = FALSE, type = "n", las = 1 )
    axis(side = 1)
    axis(side = 2, las = 1)
    box()
    abline(h = 0, lty = 2 )
    abline(v = depLines * B0, col = depLineCols, lty = 2 )
    # Put in arrows first

    arrows( x0 = diveTable$SSB[1:(nDive-1)],
            x1 = diveTable$SSB[2:(nDive)],
            y0 = diveTable$perCapSurplusProd[1:(nDive-1)],
            y1 = diveTable$perCapSurplusProd[2:(nDive)],
            col = "grey75", lwd = .8, length = .1 )
    points( x = diveTable$SSB,
            y = diveTable$perCapSurplusProd,
            bg = ptCols[(nSurf+1):nT], pch = 21 )
    text( x = diveTable$SSB[c(1,nDive-1)],
          y = diveTable$perCapSurplusProd[c(1,nDive-1)] + .1,
          labels = diveTable$Year[c(1,nDive-1)] )



    mtext(side = 1, outer = TRUE, text = "Spawning Stock Biomass (kt)",
          line = 2.5)

}
