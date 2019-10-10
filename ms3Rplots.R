# <><><><><><><><><><><><><><><><><><><><><><><><><><><>
# ms3Rplots.R
#
# Plots for ms3R.
# 
# 
#
# Author: SDN Johnson
# Date: July 25, 2019
#
# Last Update: October 8, 2019
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><>


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
                            est  = est_spt,
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
                            est  = est_spft,
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
                            est  = est_axspt,
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