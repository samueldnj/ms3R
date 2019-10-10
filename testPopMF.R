# source("control.R")
.loadFit(16)

dat   <- reports$data
fit  <- reports$repOpt

C_spft      <- dat$C_spft
N_axspt     <- fit$N_axspt
vB_axspft   <- fit$vB_axspft
sel_aspfxt  <- fit$sel_aspfxt
A_s         <- fit$A_s
M_spx       <- fit$M_spx

vB_spft     <- fit$vB_spft

wt_axsp     <- fit$meanWtAge_axsp



remN_axspt  <- array(0, dim = dim(N_axspt))
endN_axspt  <- array(0, dim = dim(N_axspt))
appZ_axspt  <- array(0, dim = dim(N_axspt))

appF_spft   <- array(0, dim = dim(C_spft) )
appF_axspft <- array(0, dim = dim(vB_axspft) )
pvB_axspft  <- array(0, dim = dim(vB_axspft) )



nS <- fit$nS
nP <- fit$nP
nX <- fit$nX
nA <- fit$nA
nF <- fit$nF
nT <- fit$nT

vBtmp_ax          <- matrix(0, nrow = nA, ncol = nX)
catNumAge_axspft  <- array(0, dim = dim(vB_axspft))

for( t in 1:(nT-1) )
  for( sIdx in 1:nS )
    for( pIdx in 1:nP )
    {
      for( fIdx in 1:nF )
      {
        # Pull biomass at age, convert to proportions
        pvB_axspft[,,sIdx,pIdx,fIdx,t] <- vB_axspft[,,sIdx,pIdx,fIdx,t];
        pvB_axspft[,,sIdx,pIdx,fIdx,t] <- pvB_axspft[,,sIdx,pIdx,fIdx,t]/sum(vB_axspft[,,sIdx,pIdx,fIdx,t])

        # Calculate numbers to remove by converting catch to weight
        catNumAge_axspft[,,sIdx,pIdx,fIdx,t] <- C_spft[sIdx, pIdx, fIdx,t] * pvB_axspft[,,sIdx,pIdx,fIdx,t] / wt_axsp[,,sIdx,pIdx]


        # Now add the catch at age in numbers to the removed fish
        remN_axspt[,,sIdx,pIdx,t] <- remN_axspt[,,sIdx,pIdx,t] + catNumAge_axspft[,,sIdx,pIdx,fIdx,t]

      }
      for( xIdx in 1:nX)
        endN_axspt[,xIdx,sIdx,pIdx,t] <- N_axspt[,xIdx,sIdx,pIdx,t] * exp(-M_spx[sIdx,pIdx,xIdx]) - remN_axspt[,xIdx,sIdx,pIdx,t]* exp(-M_spx[sIdx,pIdx,xIdx]/2)

      # Now compute Z approximation
      appZ_axspt[1:A_s[sIdx],,sIdx,pIdx,t] <- log( N_axspt[1:A_s[sIdx],,sIdx,pIdx,t] / endN_axspt[1:A_s[sIdx],,sIdx,pIdx,t])
      # Loop over fIdx again
      for( fIdx in 1:nF )
      {
        # browser()
        # This is just checking that the F is the same across age classes
        appF_axspft[1:A_s[sIdx],,sIdx,pIdx,fIdx,t] <- catNumAge_axspft[1:A_s[sIdx],,sIdx,pIdx,fIdx,t] * appZ_axspt[1:A_s[sIdx],,sIdx,pIdx,t]
        appF_axspft[1:A_s[sIdx],,sIdx,pIdx,fIdx,t] <- appF_axspft[1:A_s[sIdx],,sIdx,pIdx,fIdx,t]/N_axspt[1:A_s[sIdx],,sIdx,pIdx,t]/sel_aspfxt[1:A_s[sIdx],sIdx,pIdx,fIdx,,t]/(1 - exp(-appZ_axspt[1:A_s[sIdx],,sIdx,pIdx,t]))
        appF_spft[sIdx,pIdx,fIdx,t] <- appF_axspft[A_s[sIdx],2,sIdx,pIdx,fIdx,t]
      }

    }

Fratio <- appF_spft / fit$F_spft
Fratio[is.nan(Fratio)] <- NA

hist(Fratio)

plot( fit$F_spft, appF_spft, xlim = range(fit$F_spft), ylim = range(appF_spft,na.rm=T) )
abline( a = 0, b = 1)