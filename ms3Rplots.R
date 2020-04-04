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

plotBatchConvergenceRate <- function( groupFolder = "diffCV_newObsCV_short",
                                      prefix = "diffCV",
                                      AMlabs = c( singleStock = "singleStock",
                                                  hierMultiStock = "hierMultiStock",
                                                  dataPooled = "dataPooled",
                                                  coastwide = "coastwide",
                                                  totalAgg = "totalAgg" ) )
{
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs(here::here("Outputs",groupFolder))
  simFolderList <- simFolderList[grepl(prefix, simFolderList)]


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
              summarise(  obsCVmult = mean(projObsErrMult),
                          minConv   = min(pGoodReps),
                          maxConv   = max(pGoodReps),
                          meanConv  = mean(pGoodReps) )

    table
  }

  # Summarise and make table
  perfTableSummList <- lapply(  X = perfTableList, 
                                FUN = summariseTable )

  convRateTable <- do.call(rbind, perfTableSummList )

  browser()

  mpTable <- mpTable %>% left_join( convRateTable, 
                                    by = "simLabel" )

  CVmults <- unique(mpTable$obsCVmult)
  CVmults <- CVmults[order(CVmults)]

  AMs       <- names(AMlabs)
  AMcols    <- RColorBrewer::brewer.pal(length(AMs), "Set2")
  AMjitter  <- seq( from = -.3, to = .3, length.out = length(AMs) )

  AMwidth <- AMjitter[2] - AMjitter[1]

  # Now set up plotting area
  plot( x = c(0,length(CVmults) + 1), y = c(0,1),
        xlab = "", ylab = "", axes=FALSE,
        type = "n" )
    axis( side = 1, at = 1:length(CVmults),
          labels = CVmults )
    axis( side = 2, las = 1  )
    box()
    for( aIdx in 1:length(AMs) )
    {
      amLab <- AMs[aIdx]
      SStable <- mpTable %>% filter( AM == amLab,
                                      grepl("SS", eqbm) ) %>%
                              dplyr::select( obsCVmult,
                                              meanConv )

      SStable <- SStable[order(SStable$obsCVmult),]

      MStable <- mpTable %>% filter( AM == amLab,
                                      grepl("MS", eqbm) ) %>%
                              dplyr::select( obsCVmult,
                                              meanConv )
      MStable <- MStable[order(MStable$obsCVmult),]

      rect( xleft = 1:3 + AMjitter[aIdx] - AMwidth/2,
            xright = 1:3 + AMjitter[aIdx],
            ybottom = 0,
            ytop = SStable$meanConv,
            border = "black",
            col = AMcols[aIdx] )

      rect( xleft = 1:3 + AMjitter[aIdx],
            xright = 1:3 + AMjitter[aIdx] + AMwidth/2,
            ybottom = 0,
            ytop = MStable$meanConv,
            border = "black",
            col = AMcols[aIdx] )
      

      
    }



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
  par(  mfcol = c(nS+1,nP+1), 
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
  par(  mfcol = c(nS + 1, nP + 1),
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
plotLossTulip <- function(  sim = 1, 
                            groupFolder = "shortGrid",
                            lossType    = "rel",
                            var         = "SB_ispt",
                            save        = FALSE )
{
  # Load loss reports
  objLoss <- .loadLoss(sim = sim, groupFolder)
  simFolder <- objLoss$simFolder

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


  nReps <- dim(loss_ispt)[1]
  
  traces <- sample(1:nReps, size = 3)

  # Time steps
  fYear   <- objLoss$fYear
  tdxPlot <- tMP:nT
  years   <- seq(from = fYear, by = 1, length.out = nT)

  if( save )
    pdf( file = file.path(simFolder,paste(lossType,var,"LossEnvelopes.pdf",sep = "_")),
          width = 11, height = 8.5 )

  # Plot window - just for loss right now
  par(  mfcol = c(nS + 1, nP + 1),
        mar = c(.5,.5,.5,.5),
        oma = c(3,3,3,3) )

  speciesNames  <- c(blob$om$speciesNames,"data-pooled")
  stockNames    <- c(blob$om$stockNames,"coast-wide")

  loss_qspt <- apply( X = loss_ispt, FUN = quantile,
                      probs = c(0.025, 0.5, 0.975),
                      MARGIN = c(2,3,4) )

  plotYrs <- years[tdxPlot]


  # Loop and plot
  for( s in 1:(nS+1) )
    for( p in 1:(nP+1) )
    {
      maxLoss <- 1
      

      plot( x = range(years[tdxPlot]), y = c(-maxLoss,maxLoss),
            type = "n", axes = FALSE  )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
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

      abline( v = years[tMP], lty = 2, lwd = 3 )

      polygon(  x = c(years,rev(years)),
                y = c(loss_qspt[1,s,p,],rev(loss_qspt[3,s,p,])),
                border = NA, col = scales::alpha("grey10",.4) )
      lines(  x = years, y = loss_qspt[2,s,p,],
              col = "black", lwd = 3 )
      for( traceIdx in traces )
        lines( x = years, y = loss_ispt[traceIdx,s,p,],
                col = "black", lwd = 1 )
      grid()

      abline( h = 0, lty = 2, lwd = 1 )


    }
    mtext( side = 1, text = "Year", line = 2, outer = TRUE )
    mtext( side = 2, text = yLab, line = 2, outer = TRUE )

  if( save )
    dev.off()

} # END plotLoss()



# plotEffYieldCurves()
# Function for plotting effort based
# yield curves for each species and
# the complex in a stock area.
plotEffYieldCurves <- function( obj = blob )
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

  specCols <- RColorBrewer::brewer.pal(n = nS, "Dark2")

  par( mfrow = c(3,1), mar = c(1,1,1,1), oma = c(3,3,1,1) )
  for( p in 1:nP )
  {

    plot( x = c(0,Eseq[maxEval]), y = c(0, max(Yeq_pe[p,],na.rm = T) ),
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
                y0 = c(MSYMS_p[p],MSYMS_sp[,p]), lty = 2  )

      if(  p == 1 )
        legend( x = "topright", bty = "n",
                col = c(specCols,"black"),
                lwd = c(2,2,2,3),
                legend = c(speciesNames,"Complex") )

  }

  mtext( outer = TRUE, side = 1, text = "Total Effort Index", line = 2 )
  mtext( outer = TRUE, side = 2, text = "Equilibrium Yield (kt)", line = 2 )


  outList <- list(  EmsyMS_p = EmsyMS_p,
                    EmsySS_sp = EmsySS_sp,
                    MSYMS_p = MSYMS_p,
                    MSYSS_sp = MSYSS_sp,
                    MSYMS_sp = MSYMS_sp )
}

# plotEmpYieldCurves
# Function to plot median simuated yield 
# resulting from a grid of constant fishing 
# mortality rates - used for checking the
# reference point calculations
plotEmpYieldCurves <- function( sims = 1:11, 
                                folder = "",
                                indepVar = "F",
                                redoEmpRefCurves = FALSE )
{
  nSims <- length(sims)
  blobList <- vector(mode = "list", length = nSims)

  .loadSim(sims[1], folder = folder)

  nS <- blob$om$nS
  nP <- blob$om$nP
  nT <- blob$om$nT

  goodReps <- blob$goodReps


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

  if( indepVar == "F" )

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

  par( mfcol = c(nS+1,nP+1),
        mar = c(.5,.5,.5,.5),
        oma = c(3,3,2,2) )


  fleetCols <- RColorBrewer::brewer.pal(nF,"Set1")

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

  ctlList <- obj$ctlList
  


  fitSB_spt[fitSB_spt < 0] <- NA

  C_spft   <- obj$om$C_ispft[iRep,,,,1:t]

  C_spt    <- apply(X = C_spft, FUN = sum, MARGIN = c(1,2,4))

  # Now pull indices
  I_spft <- obj$mp$data$I_ispft[iRep,,,,1:t]
  I_spft[I_spft < 0] <- NA

  tFirstIdx <- ctlList$mp$assess$spYrFirstIdx

  if( tFirstIdx > 1 )
    I_spft[,,,tFirstIdx] <- NA

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
    newC_spt          <- array(NA, dim = c(nS,1,t))
    
    newC_spt[,1,]     <- sumC_spt[,1,1,]
    C_spt             <- newC_spt

  }

  if( ctlList$mp$assess$spDataPooled )
  {
    newSB_spt        <- apply( X = SB_spt, FUN = sum, MARGIN = c(2,3), na.rm = T )
    newVB_spt        <- apply( X = VB_spt, FUN = sum, MARGIN = c(2,3), na.rm = T )
    newtotB_spt      <- apply( X = totB_spt, FUN = sum, MARGIN = c(2,3), na.rm = T )

    newI_spft        <- apply( X = I_spft, FUN = sum, MARGIN = c(2,3,4), na.rm = T )

    SB_spt    <- SB_spt[1,,,drop = FALSE]
    SB_spt[1,,] <- newSB_spt
    VB_spt    <- VB_spt[1,,,drop = FALSE]
    VB_spt[1,,] <- newVB_spt
    totB_spt  <- totB_spt[1,,,drop = FALSE]
    totB_spt[1,,] <- newtotB_spt

    I_spft <- I_spft[1, , , , drop= FALSE]
    I_spft[1,,,] <- newI_spft
    I_spft[I_spft <= 0] <- NA

    sumC_spt          <- C_spft[1,,1,,drop = FALSE]
    sumC_spt[1,,1,]   <- apply(X = C_spft, FUN = sum, MARGIN = c(2,4))
    newC_spt          <- array(NA, dim = c(1,nP,t))
    newC_spt[1,,]     <- sumC_spt[1,,1,]
    C_spt             <- newC_spt

  }

  
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
    for( p in 1:nP )
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
    stockNames <- "Coastwide"

  if( nSS == 1 )
    speciesNames <- "Data Pooled"

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
      lines( x = yrs, y = totB_spt[s,p,], col = "black", lwd = 2 )
      lines( x = yrs, y = fitSB_spt[s,p,], col = "grey60", lwd = 1 )

      for( f in 1:nF )
      {
        # Plot scaled indices
        if( nP > nPP )
          for( pp in 1:nP )
          {
            if( f > 1 )
            points( x = yrs+ppJitter[pp], y = scaledIdx_spft[s,pp,f,],
                    pch = fleetPCH[f], bg = fleetBG[f],
                    cex = .8, col = stockCol[pp] )
          }
        if( nP == nPP) 
        {
          points( x = yrs, y = scaledIdx_spft[s,p,f,],
                  pch = fleetPCH[f], bg = fleetBG[f],
                  cex = 1.3, col = stockCol[p] )
        }
        
      }
      
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

  ctlList <- obj$ctlList

  fitSB_spt[fitSB_spt < 0] <- NA

  spTVqFleets <- ctlList$mp$assess$spTVqFleets

  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT

  # Now pull indices
  I_spft <- obj$mp$data$I_ispft[iRep,,,,1:t]
  I_spft[I_spft < 0] <- NA

  nSS <- nS
  nPP <- nP

  # Aggregate OM biomasses to match AM biomass
  if( ctlList$mp$assess$spCoastwide )
  {
    nPP <- 1
  }

  if( ctlList$mp$assess$spDataPooled )
  {
    nSS <- 1
    newI_spft        <- apply( X = I_spft, FUN = sum, MARGIN = c(2,3,4), na.rm = T )

    I_spft <- I_spft[1, , , , drop= FALSE]
    I_spft[1,,,] <- newI_spft
    I_spft[I_spft <= 0] <- NA

  }

  stdResids_spft <- array( NA, dim = c(nSS, nP, nF, t) )

  for( s in 1:nSS )
    for( p in 1:nP )
      for( f in 1:nF )
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

  par(  mfcol = c(nP,nSS), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )
  for(s in 1:nSS)
    for( p in 1:nP )
    {

      maxResid <- max(abs(stdResids_spft[s,p,,]),na.rm = T)
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
          lines(  loess.smooth(  x = yrs[nonNA], y = stdResids_spft[s,p,f,nonNA]), 
                  lwd = 2, col = fleetBG[f] )
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