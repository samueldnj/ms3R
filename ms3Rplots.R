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
        lines( x = yrs, y = C_spft[s,p,f,], col = fleetCols[f], lwd = 2, lty = 1)
        lines( x = yrs, y = TAC_spft[s,p,f,], col = fleetCols[f], lwd = 2, lty = 2)
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

plotRetroSB <- function( obj = blob, iRep = 1 )
{
  # Get biomass arrays
  SB_spt        <- obj$om$SB_ispt[iRep,,,]
  retroSB_tspt  <- obj$mp$assess$retroSB_itspt[iRep,,,,]


  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT

  speciesNames  <- blob$ctlList$opMod$species
  stockNames    <- blob$ctlList$opMod$stock
  fYear         <- blob$ctlList$opMod$fYear
  pT            <- blob$ctlList$opMod$pT

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  par(  mfcol = c(nS,nP), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )
  for(s in 1:nS)
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(SB_spt[s,p,],retroSB_tspt[,s,p,],na.rm = T) ),
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
      for( tt in 1:pT )
        lines( x = yrs, y = retroSB_tspt[tt,s,p,], col = "grey60", lwd = 1 )
  
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }

}


# diagCondition()
# Plots that help diagnose issues with conditioning
# between the fitted OM report, and the conditioned
# ms3R operating model.
diagCondition <- function(  repObj  = totRep,
                            ms3Obj   = test  )
{
  par(mfrow = c(3,3), mar = c(1.5,1,1.5,1), oma = c(3,3,3,3) )

  # Biomass RE
  plotRE_spt( repObj = repObj, omObj = ms3Obj$om, series = "SB_spt" )
  mtext( side = 3, text = "SB_spt", line = 0, font = 2)

  # Recruitment
  plotRE_spt( repObj = repObj, omObj = ms3Obj$om, series = "R_spt" )
  mtext( side = 3, text = "R_spt", line = 0, font = 2)  

  # Recruitment errors
  plotRE_spt( repObj = repObj, omObj = ms3Obj$errors, series = "omegaR_spt" )
  mtext( side = 3, text = expression(omega[R]), line = 0, font = 2)  

  # Catch
  plotRE_spft( repObj = repObj, omObj = ms3Obj$om, series = "C_spft" )
  mtext( side = 3, text = expression(C[spft]), line = 0, font = 2)    

  # F
  plotRE_spft( repObj = repObj, omObj = ms3Obj$om, series = "F_spft" )
  mtext( side = 3, text = expression(F[spft]), line = 0, font = 2) 

  # Numbers at age
  plotRE_axspt( repObj = repObj, omObj = ms3Obj$om, series = "N_axspt" )
  mtext( side = 3, text = expression(N[axspt]), line = 0, font = 2)    

  # Biomass at age
  plotRE_axspt( repObj = repObj, omObj = ms3Obj$om, series = "N_axspt" )
  mtext( side = 3, text = expression(B[axspt]), line = 0, font = 2)    
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

plotRE_spt <- function( repObj, omObj, series = "SB_spt" )
{
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT

  true_spt <- repObj[[series]]
  est_spt  <- omObj[[series]]

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

plotRE_spft <- function( repObj, omObj, series = "C_spft" )
{
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT

  true_spft <- repObj[[series]]
  est_spft  <- omObj[[series]]

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

plotRE_axspt <- function( repObj, omObj, series = "N_axspt" )
{
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT

  true_axspt <- repObj[[series]]
  est_axspt  <- omObj[[series]]

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