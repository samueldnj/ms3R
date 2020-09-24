#################################################

# Simulate blended index using different approaches
# and compare results to historical period

#################################################

setwd('~/Documents/LANDMARK/projects/2020_Herring_SISCA/ms3R/blendIdxCalcs')
source('blendIdxFuns.R')

# Run Sims for logit random walk
rw <- simIdx('simCtlFile_logitRW.txt')
save(rw, file='sims/rw.Rdata')
pdf(file='sims/simDix_rw_rep1.pdf')
plotSimIdx_i(rw, rep=1)
dev.off()

# Run Sims for logit random walk with higher SD, sd=1
rwSD1 <- simIdx('simCtlFile_logitRWsd1.txt')
save(rwSD1, file='sims/rwSD1.Rdata')
pdf(file='sims/simDix_rwSD1_rep1.pdf')
plotSimIdx_i(rwSD1, rep=1)
dev.off()

# Run Sims for binomial draw with logit random walk
binRw <- simIdx('simCtlFile_binomXlogitRW.txt')
save(binRw, file='sims/binRw.Rdata')
pdf(file='sims/simDix_binRw_rep1.pdf')
plotSimIdx_i(binRw, rep=1)
dev.off()

# Run Sims for binomial draw with logit random walk
binRwSD1 <- simIdx('simCtlFile_binomXlogitRWsd1.txt')
save(binRwSD1, file='sims/binRwSD1.Rdata')
pdf(file='sims/simDix_binRwSD1_rep1.pdf')
plotSimIdx_i(binRwSD1, rep=1)
dev.off()

# Run Sims for binomial draw with random beta distribution
beta <- simIdx('simCtlFile_binomXbeta.txt')
save(beta, file='sims/beta.Rdata')
pdf(file='sims/simDix_beta_rep1.pdf')
plotSimIdx_i(beta, rep=1)
dev.off()

projList <- list( rw    = rw,
                  rwSD1 = rwSD1,
                  binRw = binRw,
                  binRwSD1 = binRwSD1,
                  beta     = beta)


stats <- data.frame()

for(s in 1:length(projList))
  stats <- rbind(stats, statTable(projList[[s]]))

write.csv(stats,'sims/statTable.cvs',row.names=F)


# dot plots
pdf(file='sims/dotPlots.pdf')
par(mfrow=c(3,2), mar=c(4,2,1,1))
for (p in 1:3)
{
    pStats <- subset(stats, area==p)

    if(p==3)
    {
      pStats <- pStats[1,]
      pStats$simName <- 'binom'

    }  
      

    # Get historical values
    histP_p <- rw$hist$histP_pt[p,]
    
    meanHistProp  <- mean(histP_p)
    histPropZero  <- length(histP_p[histP_p==0])/length(histP_p)

    # plot mean proportions
    dotchart(x=pStats$meanPropQ50, labels=pStats$simName,
             xlim=c(0,1), xlab='surface proportion of SI')
    segments(x0=pStats$meanPropQ10, x1=pStats$meanPropQ90,
            y1=1:5, y0=1:5)
    abline(v=meanHistProp,lty=3)

    # plot mean proportion zero
    dotchart(x=pStats$meanPropZeroQ50, labels=pStats$simName,
             xlim=c(0,1), xlab='proportion of years with only dive index')
    segments(x0=pStats$meanPropZeroQ10, x1=pStats$meanPropZeroQ90,
            y1=1:5, y0=1:5)
    abline(v=histPropZero,lty=3)

}
dev.off()





 


