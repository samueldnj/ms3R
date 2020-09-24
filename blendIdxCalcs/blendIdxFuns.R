#################################################

# Plot blended index data and figure out approach 
# to simulate for projections

#################################################


# -- sim functions --
simIdx <- function(ctlFileName)
{
  #--------INITIALIZATION-----------------#

  ctrl <- initCtlFile(ctlFileName)
  hist <- ctrl$hist # info from historical blended index
  proj <- ctrl$proj # info for projections of blended index

  # Get Historical index data
  omFitFolder <- '~/Documents/LANDMARK/projects/2020_Herring_SISCA/ms3R/history'
  load(file.path(omFitFolder,hist$omFitName,paste(hist$omFitName,'.Rdata',sep='')))

  # extract proportion of index composed of surface survey for blended index years
  rI_pgt     <- reports$data$rI_pgt
  fYearBlend <- hist$fYearBlend
  lYearBlend <- hist$lYearBlend

  yrChar   <- dimnames(rI_pgt)[[3]]
  years    <- as.integer(yrChar)
  blendYrsIdx <- which(years>=fYearBlend)

  histP_pt <- rI_pgt[,4,blendYrsIdx]

  # ---- SIMULATED INDEX ----#

  nSims  <- proj$nSims
  nP     <- dim(rI_pgt)[1]
  nProjT <- proj$nProjT

  # array for holding projected proprotions
  projP_ipt      <- array(NA, dim=c(nSims, nP, nProjT))
  
  # pars for simulation
  simType_p       <- proj$type
  binomP_p        <- proj$binomP # probability for binomial draws
  logitSlope_p    <- proj$logitSlope
  logitSlopeSD_p  <- proj$logitSlopeSD
  beta1_p         <- proj$betaShape1
  beta2_p         <- proj$betaShape2


  # run sims
  for(p in 1:nP)
  {
    set.seed(hist$seed)

    if(simType_p[p]=='logitRW')
    {
      for(i in 1:nSims)
      {
        for(t in 1:nProjT)
        {
          # draw std normally distributed errors on logit-scale
          rwDevs <- rnorm(nProjT, 0, logitSlopeSD_p[p])

          if(t==1)
            projP_ipt[i,p,t] <- logit(histP_pt[p,as.character(lYearBlend)])*logitSlope_p[p] +rwDevs[t]

          if(t>1)
            projP_ipt[i,p,t] <- projP_ipt[i,p,t-1]*logitSlope_p[p] + rwDevs[t]

        }  
      }

      # transfrom to natural scale
      projP_ipt[,p,] <- expit(projP_ipt[,p,])

    }

    if(simType_p[p]=='binom')
      for(i in 1:nSims)
          projP_ipt[i,p,] <- rbinom(nProjT,1,binomP_p[p])

    if(simType_p[p]=='binomXlogitRW')
      for(i in 1:nSims)
      {
          projP_ipt[i,p,] <- rbinom(nProjT,1,binomP_p[p])
          surfaceYrs      <- which(projP_ipt[i,p,]==1)
          projSurfYrs_tt  <- vector('numeric',sum(projP_ipt[i,p,]))
          nTT <- length(projSurfYrs_tt)

          # draw std normally distributed errors on logit-scale
          rwDevs <- rnorm(nTT, 0, logitSlopeSD_p[p])

          for(tt in 1:nTT)
          {
            if(tt==1)
              projSurfYrs_tt[tt] <- logit(histP_pt[p,as.character(lYearBlend)])*logitSlope_p[p] +rwDevs[tt]

            if(tt>1)
               projSurfYrs_tt[tt] <- projSurfYrs_tt[tt-1]*logitSlope_p[p] + rwDevs[tt]
          }

          
          projP_ipt[i,p,surfaceYrs] <- expit(projSurfYrs_tt)

      } 

    if(simType_p[p]=='binomXbeta')
      for(i in 1:nSims)
      {
          projP_ipt[i,p,] <- rbinom(nProjT,1,binomP_p[p])
          surfaceYrs      <- which(projP_ipt[i,p,]==1)
          projSurfYrs_tt  <- vector('numeric',sum(projP_ipt[i,p,]))
          nTT <- length(projSurfYrs_tt)

          # random draws from beta distribution
          betaDraws <- rbeta(nTT, shape1=beta1_p[p], shape2=beta2_p[p])
          
          projP_ipt[i,p,surfaceYrs] <- betaDraws 

      }   

  }  

  hist$histP_pt   <- histP_pt
  hist$years      <- fYearBlend:lYearBlend

  proj$projP_ipt  <- projP_ipt
  proj$years      <- (lYearBlend+1):(lYearBlend+nProjT)

  report <- list(hist = hist,
                 proj = proj,
                 ctrlFileName = ctlFileName)

  return(report)
} 

initCtlFile <- function(ctlFileName='simCtlFile.txt')
{
  
  # Read in control file
  controlTable  <- .readParFile(ctlFileName )
  
  # Create control list
  ctrl <- .createList( controlTable )

  return(ctrl)

}  

# -- fit functions --

# Fit linear regression on logit-scale to time-varying proportions
fitlogitRW <- function(p, rI_pgt, fYearBlend=1988, exclude01=FALSE,
                       plot=TRUE, plotLogitScale=TRUE)
{

  # extract proportion of index composed of surface survey
  yrChar   <- dimnames(rI_pgt)[[3]]
  years    <- as.integer(yrChar)
  blendYrsIdx <- which(years>=fYearBlend)
  P_t <- rI_pgt[p,4,blendYrsIdx]

  # add 0.0001 and substract 0.0001 for 1s and zeros
  if(!exclude01)
  {
    P_t[P_t==1] <- 1-0.0001
    P_t[P_t==0] <- 0+0.0001
  }  


  if(exclude01)
  {
    blendYrsIdx <- which(P_t>0 & P_t<1)
    P_t <- P_t[blendYrsIdx]
  }  
    

  # Fit regression to estimate random walk coefficients
  logitP_t <- logit(P_t)
  dat <- data.frame(yrs = years[blendYrsIdx],
                    logitP_t = logitP_t)
  dat$logitP_tPrev <- NA
  dat$logitP_tPrev[2:nrow(dat)] <- dat$logitP_t[1:(nrow(dat)-1)]
  dat$P_t     <- P_t
  dat$P_tPrev <- expit(dat$logitP_t)

  # This should probablly have zero interecept
  m1 <- lm(logitP_t ~ 0 + logitP_tPrev, data=dat)
  m2 <- lm(logitP_t ~ logitP_tPrev, data=dat)

  
  # calculate auto-correlation coefficient
  dat <- na.omit(dat)
  autoCorP <- cor(dat$P_t, dat$P_tPrev)
  autoCorLogitP <- cor(dat$logitP_t, dat$logitP_tPrev)

  slope  <- coefficients(m1)['logitP_tPrev']
  slopeSD <- sqrt(diag(vcov(m1)))['logitP_tPrev']

  # plot
  if(plot)
  {  
    if(plotLogitScale)
    {  
      plot(x=dat$logitP_tPrev, y=dat$logitP_t, las=1,
           ylab = expression('logit propSurf'[t]),
           xlab = expression('logit propSurf'[t-1]))
      text(x=dat$logitP_tPrev+0.2, y=dat$logitP_t+0.2, labels=dat$yrs, cex=0.6)
      abline(m1)

      legend('topleft',bty='n', cex=0.8,
            legend=c(paste('slope =', round(slope,2)),
                     paste('slopeSD =', round(slopeSD,2)),
                     paste('logit propSurf acc =', round(autoCorLogitP,2))
                      )
            )
    }

    if(!plotLogitScale)
    {  
      
      plot(x=dat$p_tPrev, y=dat$p_t, las=1, 
           ylim=c(0,1), xlim=c(0,1),
           ylab = expression('lambda'[t]),
           xlab = expression('lambda'[t-1]),
           main=paste(omName,om$stock[pIdx], sep=' - ' ))
      text(x=dat$p_tPrev+0.02, y=dat$p_t+0.02, labels=dat$yrs, cex=0.6)

      # Generate predictions
      xPred <- seq(0.001,0.999,0.01) # vector of pt-1
      
      newDat <- data.frame(logitP_Prev=logit(xPred))
      u2Pred <- predict.lm(m1, newData=newDat)

      u2Pred <- coefficients(m1)['logitP_Prev']*logit(xPred)
      p2Pred <- expit(u2Pred)
      lines(x=xPred, y=p2Pred, col='red')
      
      # u2Pred <- coefficients(m2)['logitP_Prev']*logit(xPred)
      # p2Pred <- expit(u2Pred)
      # lines(x=xPred, y=p2Pred, col='red', lty=3)

      legend('topleft',bty='n', cex=0.8,
            legend=c(paste('logit slope =', round(slope,2)),
                     paste('logit slopeSD =', round(slopeSD,2)),
                     paste('propSurf acc =', round(autoCorLogitP,2))
                      )
            )

    }
  }  
}


# -- Plot functions --

plotSimIdx_i <- function(projObj, rep=1,   
                        areaNames= c('Cumshewa/Selwyn','Juan Perez/Skincuttle','Louscoone'))
{


  histP_pt <- projObj$hist$histP_pt
  histYrs  <- projObj$hist$years

  projP_pt <- projObj$proj$projP_ipt[rep,,]
  projYrs  <- projObj$proj$years

  nP <- dim(projP_pt)[1]

  par(mfrow=c(nP,1), mar=c(3,4,1,0.2))
  for(p in 1:3)
  { 
    plot(x=c(histYrs,projYrs), y=c(histP_pt[p,],projP_pt[p,]) ,
         ylab='surface proportion of SI (%) ', type='n',
         xlab='',las=1, col='blue',pch=16,
         main=areaNames[p])
    points(x=histYrs, y=histP_pt[p,], col='black', pch=16)
    points(x=projYrs, y=projP_pt[p,], col='blue', pch=16)

    # add lines for mean proportions
    abline(v=max(histYrs),lty=3)
    abline(h=c(mean(histP_pt[p,]),mean(projP_pt[p,])), col=c('black','blue'))

    # add lines for proportions of zeros
    histPzero <- length(histP_pt[p,][histP_pt[p,]==0])/length(histP_pt[p,])
    projPzero <- length(projP_pt[p,][projP_pt[p,]==0])/length(projP_pt[p,])
    abline(h=c(mean(histPzero),mean(projPzero)), col=c('black','blue'), lty=3)
  }  

}

# --- Stat functions ----


statTable <- function(simObj)
{
  projP_ipt <- simObj$proj$projP_ipt
  nP <- dim(projP_ipt)[2]

   # First, create a data.frame of NAs with a row for each nP
  colLabels <- c( "simName",
                  "type",
                  "area",
                  "meanPropQ50",
                  "meanPropQ90",
                  "meanPropQ10",
                  "medianPropQ50",
                  "medianPropQ90",
                  "medianPropQ10",
                  "meanPropZeroQ50",
                  "meanPropZeroQ90",
                  "meanPropZeroQ10"
                  )

  stats <- matrix( NA,  ncol = length(colLabels),
                            nrow = nP )
  
  colnames(stats) <- colLabels
  stats <- as.data.frame(stats)

  # sim name and stock area
  stats$simName <- simObj$proj$name
  stats$type    <- simObj$proj$type

  for(p in 1:nP)
  {
    stats$area[p] <- p

    projP_it <- projP_ipt[,p,]
    
    # calculate mean and median proportion
    mean_i   <- apply(projP_it, FUN=mean, MAR=c(1))
    median_i <- apply(projP_it, FUN=median, MAR=c(1))

    # calculate proportion of 0s and 1s
    nReps <- dim(projP_it)[1]
    propZero_i <- vector('numeric', length=nReps)
    for(i in 1:nReps)
    {
      projP_t    <- projP_it[i,]

      nZero         <- length(projP_t[projP_t==0])
      propZero_i[i] <- nZero/length(projP_t)
    }  

    stats[p,"meanPropQ50"] <- median(mean_i)
    stats[p,"meanPropQ90"] <- quantile(mean_i, 0.90)
    stats[p,"meanPropQ10"] <- quantile(mean_i, 0.10)

    stats[p,"medianPropQ50"] <- median(median_i)
    stats[p,"medianPropQ90"] <- quantile(median_i, 0.90)
    stats[p,"medianPropQ10"] <- quantile(median_i, 0.10)

    stats[p,"meanPropZeroQ50"] <- median(propZero_i)
    stats[p,"meanPropZeroQ90"] <- quantile(propZero_i, 0.90)
    stats[p,"meanPropZeroQ10"] <- quantile(propZero_i, 0.10)
  }  

  return(stats)
}


# -- tools --

expit <- function(x)
{
exp(x)/(1+exp(x))
}


logit <- function(x)
{
log(x/(1-x))
}

# .readParFile   (reads an ASCII file with 1 comment line, header, data frame)
# Purpose:      Reads an ASCII file: 1 comment, 1 header, space-delimited
#               data frame usually containing columns "parameter" and "value".
# Parameters:   parFile is a character string indicating the input file.
# Returns:      result, a data frame.
# Source:       A.R. Kronlund
.readParFile <- function( parFile="inputFile.par" )
{
  # Read the file and store as a dataframe.
  result <- read.table( file=parFile, as.is=TRUE, header=TRUE, skip=1,
                        quote="",sep=" " )
  result
}

.createList <- function( obj )
{
  # Input  a data frame with columns "parameter" and "value".
  # Output a list with elements named as parameter and values in "value".

  result <- list()

  # Shut off whining, then coerce to numeric to let NA indicate non-numerics.
  options( warn=-1 )
  numericVal <- as.numeric( obj[,"value"] )
  options( warn=0 )

  for ( i in 1:nrow(obj) )
  {
    # Value is numeric, build the parse string.
    if ( !is.na(numericVal[i]) )
      listText <- paste( "result$",obj[i,"parameter"],"=",
                    obj[i,"value"],sep="" )
    # Value is character, build the parse string.
    else
      listText <- paste( "result$",obj[i,"parameter"],"=",
                  obj[i,"value"], sep="" )

    # ARK: At one point I had this code, presumably to handle a different
    #      input format convention, perhaps assuming "value" was all character.
    #                   sQuote(obj[i,"value"]),sep="" )
    
    # Evaluate the parsed string.
    eval( parse( text=listText ) )
  }
  result
}








