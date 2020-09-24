#################################################

# Plots for blended surface and dive survery index
# for historical period and simulation

#################################################

setwd('~/Documents/LANDMARK/projects/2020_Herring_SISCA/ms3R/blendIdxCalcs')

# Load in historical data from SISCA OM fits
load('~/Documents/LANDMARK/projects/2020_Herring_SISCA/ms3R/history/fit_parBatblendIdx2/fit_parBatblendIdx2.RData')

# extract proportions of survey index from 1988:2019
rI_pgt   <- reports$data$rI_pgt[,,38:69]
combI_pt <- reports$data$combI_pt[,38:69]

areaNames<- c('Cumshewa/Selwyn','Juan Perez/Skincuttle','Louscoone')

# Calculate probability that there is any surface survey
binom_p <- c(NA,NA,NA)

for(p in 1:3)
  binom_p[p] <- length(which(rI_pgt[p,4,]>0))/length(rI_pgt[p,4,])



# Plot proportions
pdf(file.path('propSurfIndex.pdf'))
par(mfrow=c(3,1), mar=c(3,4,0.2,0.2))
for(p in 1:3)
{ 
  plot(y=rI_pgt[p,4,],x=1988:2019, ylab='surface proportion of SI (%) ',
     xlab='',las=1, col='blue',pch=16)
  legend('left',bty='n',
        legend=paste('Proportion of years with surface index = ',
                              round(binom_p[p]*100,1),'%',sep=''))

}   
dev.off()

# Plot years with non-zero proportions for surface survey
pdf(file.path('nonZeroSurfaceProps.pdf'))
par(mfrow=c(3,1), mar=c(3,4,0.2,0.2))
for(p in 1:3)
{
    yrs      <- which(rI_pgt[p,4,] > 0)
    propSurf_t <- rI_pgt[p,4,yrs]

    dat <- data.frame(propSurf_t    = propSurf_t[2:length(yrs)],
                      propSurf_tPrev  = propSurf_t[1:(length(yrs)-1)]
                      )

  plot(y=dat[,1],x=dat[,2], 
    ylab='PropSurf_t ', xlab='PropSurf_t-1',las=1, col='blue',pch=16)


} 
  
dev.off()


# Plot years with non-zero proportions for surface survey and estimate parameters for beta distribution
pdf(file.path('betaFits.pdf'))
par(mfrow=c(2,2), mar=c(3,4,2,0.2))
for(p in 1:2)
{
  yrs      <- which(rI_pgt[p,4,] > 0)
  propSurf_t <- rI_pgt[p,4,yrs]

  betaFits <- ebeta(propSurf_t)
  alpha <- betaFits$parameters[1]
  beta  <- betaFits$parameters[2]

  pdfX <-seq(0.01,0.99,0.01)
  pdf <- dbeta(x = pdfX, shape1 = alpha, shape2 = beta)
  meanBeta <- alpha/(alpha+beta) 
  
  hist(propSurf_t, freq=FALSE, breaks=seq(0,1,0.1),
       xlab='% surface survey', main=areaNames[p])
  abline(v=mean(propSurf_t), lty=3)

  plot(x=pdfX, y=pdf, type='l',ylab='density',
       main=paste(areaNames[p], 'Beta Pars'))
  abline(v=meanBeta, lty=3)
  

  legend('topright',bty='n', cex=0.8,
            legend=c(paste('alpha =', round(alpha,5)),
                     paste('beta =', round(beta,5))))


} 
  
dev.off()


# Plot logit fits to proportions for C/S and JP/S, including zero years
pdf(file.path('logitFits.pdf'))
par(mfrow=c(2,1), mar=c(3,4,0.2,0.2))
for (p in 1:2)
  fitlogitRW(p=p, rI_pgt=rI_pgt)
dev.off()

# Plot logit fits to proportions for C/S and JP/S, excluding zero and 1 years
pdf(file.path('logitFits_exclude01.pdf'))
par(mfrow=c(2,1), mar=c(3,4,0.2,0.2))
for (p in 1:2)
  fitlogitRW(p=p, rI_pgt=rI_pgt, exclude01=TRUE)
dev.off()



