####################################################################

# Calculate values for Harvest Control Rules for HG Herring
# - weighted q across OM grid for empirical assessment rule
# - weighted B0 across OM grid for estimating stock status for HCR

###################################################################


omNames <- paste("fit_parBatHGherringMCMC",1:18,sep = "")


calcWtdPars <- function( omFitName )
{

  # Get Historical index data
  omFitFolder <- '../history'
  load(file.path(omFitFolder,omFitName,paste(omFitName,'.Rdata',sep='')))

  # --- weighted q ----

  # extract proportion of index composed of surface survey for blended index years
  rI_pgt     <- reports$data$rI_pgt
  fYearBlend <- 1988
  lYearBlend <- 2019

  yrChar   <- dimnames(rI_pgt)[[3]]
  years    <- as.integer(yrChar)
  blendYrsIdx <- which(years>=fYearBlend)

  # proportion of index composed by surface survey
  histPs_pt <- rI_pgt[,4,blendYrsIdx]
  meanPs_p <- apply(histPs_pt, FUN=mean, MAR=c(1))

  # calculate weighted q using mean posterior q values
  q_ipg    <- reports$posts$q_ipg
  meanq_pg <- apply(q_ipg, FUN=mean, MAR=c(2,3))

  qWtd_p <- rep(NA,3)
  for(p in 1:3)
    qWtd_p[p] <- meanq_pg[p,4]*meanPs_p[p] + meanq_pg[p,5]*(1-meanPs_p[p])

  # --- weighted B0 ----
  B0_ip <- reports$posts$B0_ip
  wtdB0_p <- apply(B0_ip, FUN=mean, MAR=c(2))

  outMat <- matrix(NA, nrow = 3, ncol=2)
  colnames(outMat) <- c("qWtd","B0")
  outMat[,1] <- qWtd_p
  outMat[,2] <- wtdB0_p

  outMat
}

outMats_pPm <- vapply( X = omNames, FUN = calcWtdPars, FUN.VALUE = array(0, dim = c(3,2)) )

meanPars <- apply(X = outMats_pPm, FUN = mean, MARGIN = c(1,2))

