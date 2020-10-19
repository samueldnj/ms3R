#-----------------------------------------------------
# Generate reference point curves for different areas 
# and different MPs for HG Herring fisheries in HG
#
# Created: B Doherty
# Date: Sept 28, 2020
#
#-----------------------------------------------------

library(here)


# load reference pt calcs and plotting functions from SISCA
source('~/Documents/LANDMARK/projects/2020_Herring_SISCA/siscaom/SISCAplots.R')
source('~/Documents/LANDMARK/projects/2020_Herring_SISCA/siscaom/refPts.R')

# load base OM from SISCA
load(here::here('history/fit_MSbase_MCMC/fit_MSbase_MCMC.RData'))


# Generate equilibrium relationship reference point curves for each each area for:
# a) Yield Per Recruit for increasing fishing mortality (F), 
# b) Spawning biomass (SB) per recruit for increasing fishing mortality (F), 
# c) Yield for increasing fishing mortality (F), 
# d) Spawning biomass for increasing fishing mortality (F), 
# e) Recruitment for different levels of spawning biomass, and 
# f) Yield for different levels of spawning biomass

# ref curves plot directory
plotDir <- here('refPtCalcs')

# Calculate reference points:
# Seine Roe fishery
repObjSR <- calcRefPts(reports$repOpt, fleetIdx=2) 

# Spawn on kelp fishery
repObjSOK <- calcRefPts(reports$repOpt, fleetIdx=6)

# Spawn on kelp fishery with Seine Roe fishery in JP/S for quota > 907 t
repObjSOKSR <- calcRefPts(reports$repOpt, fleetIdx=6, 
						  sokCatchLim_p=c(1e6,0.907,1e6))

# get estimate psi = SOK product/ponded fish
psi_pt <- reports$repOpt$psi_pt
psi_pt[psi_pt==0] <- NA #assign NAs for non-estimated years
psi_p <- apply(psi_pt, FUN=mean, MAR=1, na.rm=T) 


# arguments for reference points plots
fxListEggs = list(	pIdx=NULL, 
            		yieldType='Eggs',
            		gamma=0.0024)

fxListCatch = list(	pIdx=NULL, 
            		yieldType='Catch',
            		gamma=0.0024)


fxListSOK = list(	pIdx=NULL, 
            		yieldType='sokProduct',
            		gamma=0.0024)

# Seine Roe reference curves
pdf(file=file.path(plotDir,'refPtsF_SR_CatchYield.pdf'))
plotRefPtsF_p(repObj=repObjSR, fxList=fxListCatch)
dev.off()

pdf(file=file.path(plotDir,'refPtsF_SR_EggYield.pdf'))
plotRefPtsF_p(repObj=repObjSR, fxList=fxListEggs)
dev.off()

pdf(file=file.path(plotDir,'refPtsU_SR_CatchYield.pdf'))
plotRefPtsU_p(repObj=repObjSR, fxList=fxListCatch)
dev.off()

pdf(file=file.path(plotDir,'refPtsU_SR_EggYield.pdf'))
plotRefPtsU_p(repObj=repObjSR, fxList=fxListEggs)
dev.off()


pdf(file=file.path(plotDir,'refPtsU_SR_sokProduct.pdf'))
plotRefPtsU_p(repObj=repObjSR, fxList=fxListSOK)
dev.off()


# Spawn on kelp reference curves
pdf(file=file.path(plotDir,'refPtsF_SOK_CatchYield.pdf'))
plotRefPtsF_p(repObj=repObjSOK, fxList=fxListCatch)
dev.off()

pdf(file=file.path(plotDir,'refPtsF_SOK_EggYield.pdf'))
plotRefPtsF_p(repObj=repObjSOK, fxList=fxListEggs)
dev.off()

pdf(file=file.path(plotDir,'refPtsU_SOK_CatchYield.pdf'))
plotRefPtsU_p(repObj=repObjSOK, fxList=fxListCatch)
dev.off()

pdf(file=file.path(plotDir,'refPtsU_SOK_EggYield.pdf'))
plotRefPtsU_p(repObj=repObjSOK, fxList=fxListEggs)
dev.off()

pdf(file=file.path(plotDir,'refPtsU_SOK_sokProduct.pdf'))
plotRefPtsU_p(repObj=repObjSOK, fxList=fxListSOK)
dev.off()

pdf(file=file.path(plotDir,'sokRefPtsU_SOK.pdf'))
plotSokRefPtsU_p (repObj=repObjSOK)
dev.off()


# Harvest rate reference curves for SOK and Seine Roe fleet
pdf(file=file.path(plotDir,'refPtsF_SOKSR_CatchYield.pdf'))
plotRefPtsF_p(repObj=repObjSOKSR, fxList=fxListCatch)
dev.off()

pdf(file=file.path(plotDir,'refPtsF_SOKSR_EggYield.pdf'))
plotRefPtsF_p(repObj=repObjSOKSR, fxList=fxListEggs)
dev.off()

pdf(file=file.path(plotDir,'refPtsU_SOKSR_CatchYield.pdf'))
plotRefPtsU_p(repObj=repObjSOKSR, fxList=fxListCatch)
dev.off()

pdf(file=file.path(plotDir,'refPtsU_SOKSR_EggYield.pdf'))
plotRefPtsU_p(repObj=repObjSOKSR, fxList=fxListEggs)
dev.off()

# plot U vs F for SOK and Seine Roe fleet in JP/S
refCurves <- repObjSOKSR$refPts$refCurves

# add lines for Fmsy and Umsy under SR and SOK only
FmsySOK	<- repObjSOK$refPts$FmsyRefPts$Fmsy_p[2]
UmsySOK <- repObjSOK$refPts$FmsyRefPts$Umsy_p[2]
FmsySR	<- repObjSR$refPts$FmsyRefPts$Fmsy_p[2]
UmsySR  <- repObjSR$refPts$FmsyRefPts$Umsy_p[2]

# Identify point where yield is greater than 907t
F          	<- refCurves$F
U 	   		<- refCurves$Ueq_pf[2,]
Udead  		<- refCurves$Udead_pf[2,]

# SOK F
Fsok       <- repObjSOK$refPts$refCurves$F
UdeadSOK   <- repObjSOK$refPts$refCurves$Udead_pf[2,]

Yeq_f      <- repObjSOKSR$refPts$refCurves$Yeq_pf[2,]

pdf(file=file.path(plotDir,'UvsF_SOKSR_JPS.pdf'))
plot(x=refCurves$F,y=refCurves$Ueq_pf[2,], xlab='Fishing mortality',
	ylab='Harvest Rate', type='l', lwd=2)
lines(x=refCurves$F,y=UdeadSOKSR, col='blue', lty=1)
lines(x=Fsok,y=UdeadSOK, col='green', lty=3)

abline(h=c(UmsySOK, UmsySR), lty=c(2,3))
abline(v=FmsySOK, lty=1)

legend('bottomright', bty='n',
		legend=c('SOKSR HR','SOKSR deadHR','SOK deadHR','Umsy','Fmsy'),
		lwd=c(2,1,1,1,NA),
		lty=c(1,1,3,3,NA),
		pch=c(NA,NA,NA,NA,'|'),
		col=c('black','blue','green','black','black') )
dev.off()





