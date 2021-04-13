# <><><><><><><><><><><><><><><><><><><><><><><><><><><>
# stats.R
#
# Functions for calculating statistics from simulations
# produced by the ms3R closed loop simulation platform.
#
# Author: SDN Johnson
# Date: December 6, 2019
#
# Last Update: Dec 6, 2019
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><>

# makeDynEqbriaTab( )
# Uses reference points object from a simulation
# to create a table of multispecies catch and economic 
# reference points.
makeDynEqbriaTab <- function( folder = "./Outputs/omniRuns_priceDevPED",
                              scenario = "noCorr" )
{
  csvFiles <- list.files(folder)
  csvFiles <- csvFiles[grepl(x = csvFiles, pattern = ".csv")]

  nFiles <- length(csvFiles)
  tabLabs <- unlist(stringr::str_split(csvFiles,".csv"))[2*(1:nFiles)-1]

  csvPath <- file.path(folder,csvFiles)
  tabs <- lapply(X = csvPath, FUN = read.csv, header = TRUE)
  names(tabs) <- tabLabs

  stockNames    <- unique(tabs$dynEqHR$stock)[1:3]
  speciesNames  <- unique(tabs$dynEqHR$species)[1:3]


  nS <- length(speciesNames)
  nP <- length(stockNames)
  
  colNames <- c("Stock",
                "Species",
                "Emsy",
                "RentMSY",
                "MSY",
                "Bmsy",
                "Umsy",
                "Emey",
                "MEY",
                "Cmey",
                "Bmey",
                "Umey")

  tableFrame <- matrix( NA, nrow = 12, ncol = length(colNames))
  colnames(tableFrame) <- colNames

  tableFrame <- as.data.frame(tableFrame)
  msyCol <- paste0(scenario,".totCat")
  meyCol <- paste0(scenario,".totProfit")

  for( p in 1:nP )
  {
    stockRow  <- (p-1) * (nS+1) + 1  

    tableFrame$Stock[stockRow]    <- stockNames[p]

    tableFrame$Emsy[stockRow]     <- round(tabs$dynEqEff[p,msyCol],2)
    tableFrame$RentMSY[stockRow]  <- round(tabs$dynEqPrSp[p,msyCol],2)

    tableFrame$Emey[stockRow]     <- round(tabs$dynEqEff[p,meyCol],2)
    tableFrame$MEY[stockRow]      <- round(tabs$dynEqPrSp[p,meyCol],2)

    for( s in 1:nS )
    {
      specRow   <- (p-1) * (nS+1) + s + 1

      tableFrame$Species[specRow]   <- speciesNames[s]

      tabRow <- (p - 1) * nS + s

      tableFrame$MSY[specRow]       <- round(tabs$dynEqCat[tabRow,msyCol],2)
      tableFrame$Bmsy[specRow]      <- round(tabs$dynEqBio[tabRow,msyCol],2)
      tableFrame$Umsy[specRow]      <- round(tabs$dynEqHR[tabRow,msyCol],3)

      tableFrame$Cmey[specRow]      <- round(tabs$dynEqCat[tabRow,meyCol],2)
      tableFrame$Bmey[specRow]      <- round(tabs$dynEqBio[tabRow,meyCol],2)
      tableFrame$Umey[specRow]      <- round(tabs$dynEqHR[tabRow,meyCol],3)
    }


  }

  tableFrame

} # END makeDynEqbriaTab()


makeSpeciesDynEqTable <- function( file = "./Outputs/omniRuns_priceDevPED/dynEqHR.csv" )
{
  specTable <- read.csv(file, stringsAsFactors = FALSE, header = TRUE, row.names = 1)

  specTable <- t(specTable)
  colnames(specTable) <- paste(specTable[1,],specTable[2,],sep = ".")
  specTable <- specTable[-c(1,2,nrow(specTable)),]
  specTable <- specTable[,-ncol(specTable)]


  scen <- c("Steady State","Steady State",
            "No Correlation","No Correlation",
            "Recruitment Correlation","Recruitment Correlation",
            "Price Correlation","Price Correlation",
            "Recruitment and Price Correlation","Recruitment and Price Correlation")

  objFun <- rep(c("Catch","Rent"),5)
              
  scenOrder <- scen[c(1,3,5,7,9)]

  specTable <- as.data.frame(specTable)
  specTable$scen <- scen
  specTable$objFun <- objFun

  specTable <- specTable %>%
                arrange( match(objFun,c("Catch","Rent")),
                         match(scen,scenOrder) ) %>%
                relocate(scen) %>%
                relocate(objFun)


  specTable
}

# makeStatMSEqbriaTab()
# Uses reference points object from a simulation
# to create a table of multispecies catch and economic 
# reference points.
makeStatMSEqbriaTab <- function( obj = blob )
{
  rp <- obj$rp[[1]]

  nS <- obj$om$nS
  nP <- obj$om$nP
  nT <- obj$om$nT
  nF <- obj$om$nF
  tMP<- obj$om$tMP

  stockNames    <- obj$om$stockNames
  speciesNames  <- obj$om$speciesNames

  EmeyRefPts    <- rp$EmeyRefPts
  EmsyMSRefPts  <- rp$EmsyMSRefPts

  # Get CW ref pts
  if(!is.null(EmeyRefPts$cwEconYieldCurves))
  {
    Emey_p      <- EmeyRefPts$cwEconYieldCurves$cwEmey_p
    MEY_p       <- EmeyRefPts$cwEconYieldCurves$cwMEY_p
    Bmey_sp     <- EmeyRefPts$cwEconYieldCurves$cwBmey_sp
    Ymey_sp     <- EmeyRefPts$cwEconYieldCurves$cwYmey_sp    
  } else {
    Emey_p  <- EmeyRefPts$Emey_p
    MEY_p   <- EmeyRefPts$MEY_p
    Bmey_sp <- EmeyRefPts$Bmey_sp
    Ymey_sp <- EmeyRefPts$Ymey_sp  
  }

  Emsy_p  <- EmsyMSRefPts$EmsyMS_p
  MSY_sp  <- EmsyMSRefPts$YeqEmsy_sp
  Bmsy_sp <- EmsyMSRefPts$BeqEmsy_sp

  qF_spf     <- obj$om$qF_ispft[1,,,,tMP]

  # gonna have to solve for rent at Emsy_p
  E           <- rp$refCurves$EffCurves$E
  econYmsy_p  <- array(NA, dim = c(nP))
  econYeq_pe  <- rp$EmeyRefPts$econYeq_pe
  for( p in 1:nP )
    econYmsy_p[p] <- getSplineVal( x = E, y = econYeq_pe[p,], p = Emsy_p[p])


  colNames <- c("Stock",
                "Species",
                "qComm",
                "Emsy",
                "RentMSY",
                "MSY",
                "Bmsy",
                "Umsy",
                "Emey",
                "MEY",
                "Cmey",
                "Bmey",
                "Umey")

  tableFrame <- matrix( NA, nrow = 12, ncol = length(colNames))
  colnames(tableFrame) <- colNames

  tableFrame <- as.data.frame(tableFrame)

  for( p in 1:nP )
  {
    stockRow  <- (p-1) * (nS+1) + 1

    tableFrame$Stock[stockRow]    <- stockNames[p]

    tableFrame$Emsy[stockRow]     <- round(Emsy_p[p],2)
    tableFrame$RentMSY[stockRow]  <- round(econYmsy_p[p],2)

    tableFrame$Emey[stockRow]     <- round(Emey_p[p],2)
    tableFrame$MEY[stockRow]      <- round(MEY_p[p],2)

    for( s in 1:nS )
    {
      specRow   <- (p-1) * (nS+1) + s + 1

      tableFrame$Species[specRow]   <- speciesNames[s]
      tableFrame$qComm[specRow]     <- signif(qF_spf[s,p,2],2)

      tableFrame$MSY[specRow]        <- round(MSY_sp[s,p],2)
      tableFrame$Bmsy[specRow]      <- round(Bmsy_sp[s,p],2)
      tableFrame$Umsy[specRow]      <- round(MSY_sp[s,p]/Bmsy_sp[s,p],3)

      tableFrame$Cmey[specRow]      <- round(Ymey_sp[s,p],2)
      tableFrame$Bmey[specRow]      <- round(Bmey_sp[s,p],2)
      tableFrame$Umey[specRow]      <- round(Ymey_sp[s,p]/Bmey_sp[s,p],3)
    }


  }

  return(tableFrame)
} # END makeStatMSEqbriaTab()

# makeParEstTable
makeParEstTable <- function( obj )
{
  repObj  <- obj$ctlList$opMod$histRpt
  rp      <- obj$rp[[1]]
  # Pull number of species/stocks from 
  # the blob
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT

  species <- obj$om$speciesNames
  stock   <- obj$om$stockNames

  FmsyRefPts <- rp$FmsyRefPts

  # need to combine these all into a table - also want B0/R0 etc
  B0_sp       <- round(repObj$B0_sp,2)
  R0_sp       <- round(repObj$R0_sp,2)
  h_sp        <- round(repObj$h_sp,2)
  M_xsp       <- round(repObj$M_xsp,2)
  Bmsy_sp     <- round(FmsyRefPts$BeqFmsy_sp,2)
  Fmsy_sp     <- round(FmsyRefPts$Fmsy_sp,2)
  MSY_sp      <- round(FmsyRefPts$YeqFmsy_sp,2)
  Umsy_sp     <- round(MSY_sp/Bmsy_sp,2)
  SSBT_sp     <- round(repObj$SB_spt[,,nT,drop = FALSE],2)
  DT_sp       <- round(SSBT_sp[,,1]/B0_sp,2)
  DmsyT_sp    <- round(SSBT_sp[,,1]/Bmsy_sp,2)
  
  nTabRow <- nS * nP
  
  
  tabColNames <- c( "Species",
                    "Stock",
                    "B0",
                    "R0",
                    "M_m",
                    "M_f",
                    "h",
                    "Bmsy",
                    "Umsy",
                    "MSY",
                    "SSB_T",
                    "D_T",
                    "Dmsy_T" )
  
  refPtsTab <- matrix(NA, nrow = nTabRow, ncol = length(tabColNames) )
  colnames(refPtsTab) <- tabColNames

  refPtsTab <- as.data.frame(refPtsTab)
  for( s in 1:nS )
    for(p in 1:nP )
    {
      rowIdx <- (s - 1) * nP + p
      refPtsTab$Species[rowIdx]   <- species[s]
      refPtsTab$Stock[rowIdx]     <- stock[p]
      refPtsTab$B0[rowIdx]        <- B0_sp[s,p]
      refPtsTab$R0[rowIdx]        <- R0_sp[s,p]
      refPtsTab$M_m[rowIdx]       <- M_xsp[1,s,p]
      refPtsTab$M_f[rowIdx]       <- M_xsp[2,s,p]
      refPtsTab$h[rowIdx]         <- h_sp[s,p]
      refPtsTab$Bmsy[rowIdx]      <- Bmsy_sp[s,p]
      refPtsTab$Umsy[rowIdx]      <- Umsy_sp[s,p]
      refPtsTab$MSY[rowIdx]       <- MSY_sp[s,p]
      refPtsTab$SSB_T[rowIdx]     <- SSBT_sp[s,p,1]
      refPtsTab$D_T[rowIdx]       <- DT_sp[s,p]
      refPtsTab$Dmsy_T[rowIdx]    <- DmsyT_sp[s,p]
    }


  refPtsTab
} # END makeParEstTab

# dynEqbriaTab()
# Collects output from omniscient manager simulations
# and turns them into tables of dynamic equilibria
dynEqbriaTab <- function( groupFolder = "omniRuns_priceDevPED",
                          mpFilter = "",
                          distYrs = 2040:2060,
                          qProbs = c(0.05,0.5,0.95),
                          econYieldFile = "./Outputs/sim_baseRun_invDem/sim_baseRun_invDem.Rdata",
                          scenOrder = c("noCorr_infPED","noCorr_infPED_fixGDP"))
{
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
                filter( grepl( mpFilter, mp ) ) %>%
                arrange(scenario,mp)

  scenLabels <- unique(info.df$scenario)

  # Read in the modelStates we extracted
  nModels <- nrow(info.df)
  modelList <- vector(mode = "list", length = nModels)
  names(modelList) <- info.df$simLabel
  for( k in 1:nModels )
  {
    simLabel <- info.df$simLabel[k]
    modelStatePath <- here::here("Outputs",groupFolder,simLabel,paste(simLabel,"_modelStates.RData",sep=""))
    load(modelStatePath)

    modelList[[k]] <- outList

  }
  outList <- NULL

  load(econYieldFile)
  baseBlob <- blob

  discRate <- blob$ctlList$opMod$discountRate

  # Get reference points
  rp  <- baseBlob$rp[[1]]
  nS  <- modelList[[1]]$nS
  nP  <- modelList[[1]]$nP
  nT  <- modelList[[1]]$nT
  nF  <- modelList[[1]]$nF
  tMP <- modelList[[1]]$tMP
  pT  <- nT - tMP + 1

  speciesNames <- modelList[[1]]$speciesNames
  stockNames <- modelList[[1]]$stockNames

  fYear <- modelList[[1]]$fYear
  distYearIdx <- distYrs - fYear + 1

  objFuns <- c("totCat","totProfit")

  egHeaders <- expand.grid(scen = scenOrder,objFun = objFuns) %>%
                mutate_all(as.character) %>%
                arrange(match(scen, scenOrder))

  colHeaders <- apply(X = egHeaders, FUN = paste, collapse = ".", MARGIN = 1)


  colHeaders <- c(  "species",
                    "stock",
                    "Xmsy",
                    "Xmey",
                    colHeaders,
                    "var" )

  specEqTab <- matrix(NA, nrow = nS * nP + 1, ncol = length(colHeaders) )
  areaEqTab <- matrix(NA, nrow = nP + 1, ncol = length(colHeaders) )

  colnames(specEqTab) <- colHeaders
  colnames(areaEqTab) <- colHeaders

  specEqTab <- as.data.frame(specEqTab)
  areaEqTab <- as.data.frame(areaEqTab)

  areaEqTab$species <- "DER"
  areaEqTab$stock   <- c(stockNames,"BC")

  specEqTab[nS*nP+1,"species"] <- "DER"
  specEqTab[nS*nP+1,"stock"]   <- "BC"



  
  # equlibria:
  # - Species specific
  #     Biomass 
  dynEqBio  <- specEqTab
  dynEqBio$var <- "SSB"
  #     Catch
  dynEqCat  <- specEqTab
  dynEqCat$var <- "Ct"
  #     HR
  dynEqHR   <- specEqTab
  dynEqHR$var <- "Ut"
  # - Area sum
  #     Effort
  dynEqEff  <- areaEqTab
  dynEqEff$var <- "Et"
  #     Producer Surplus
  dynEqPrSp <- areaEqTab
  dynEqPrSp$var <- "ProdSurplus"

  # Discounted total rent
  dynEqDiscRent     <- areaEqTab
  dynEqDiscRent$var <- "discRent"

  # Percentile tabs
  # Producer surplus
  dynEqPrSpDist     <- areaEqTab
  dynEqPrSpDist$var <- "ProdSurplus"

  #     Catch
  dynEqCatDist      <- specEqTab
  dynEqCatDist$var  <- "Ct"

  # Harvest rates
  dynEqHRDist       <- specEqTab
  dynEqHRDist$var   <- "Ut"

  # Biomass
  dynEqBioDist      <- specEqTab
  dynEqBioDist$var  <- "SSB"

  dynEqDiscRentDist <- areaEqTab
  dynEqDiscRentDist$var <- "discRent"


  # Pull static RPs
  # Species
  Bmsy_sp     <- rp$EmsyMSRefPts$BeqEmsy_sp
  # Bmey_sp     <- rp$EmeyRefPts$Bmey_sp
  # Ymey_sp     <- rp$EmeyRefPts$Ymey_sp
  MSY_sp      <- rp$EmsyMSRefPts$YeqEmsy_sp

  # Area
  Emey_p      <- rp$EmeyRefPts$cwEconYieldCurves$cwEmey_p
  Emsy_p      <- rp$EmsyMSRefPts$EmsyMS_p
  MEY_p       <- rp$EmeyRefPts$cwEconYieldCurves$cwMEY_p


  # Solve for Bmey with coastwide demand curve
  Beq_spe     <- rp$refCurves$EffCurves$Beq_spe
  Yeq_spe     <- rp$refCurves$EffCurves$Yeq_spe
  # Want to solve for econ yield at MSY
  econYeq_pe  <- rp$EmeyRefPts$econYeq_pe

  E           <- rp$refCurves$EffCurves$E

  econYmsy_p  <- array(NA, dim = c(nP))
  Bmey_sp     <- array(NA, dim = c(nS,nP))
  Ymey_sp     <- array(NA, dim = c(nS,nP))
  Umey_sp     <- array(NA, dim = c(nS,nP))

  Umsy_sp     <- MSY_sp / Bmsy_sp
  

  
  # Loop over areas, then species
  for( k in 1:nModels)
  {
    # Get scenario name and objective function
    scenName  <- info.df$scenario[k]
    mpName    <- info.df$mp[k]
    objFun <- str_split(mpName,"_")[[1]][2]
    objFun <- str_split(objFun,"1")[[1]][1]

    colName <- paste(scenName,objFun,sep = ".")

    # make discounted rent - wasn't done with pullModelStates
    dynDiscRent_ipt <- modelList[[k]]$modelStates$profit_ipft[,,2,tMP:nT]
    for(p in 1:nP)
      for( t in 1:pT)
        dynDiscRent_ipt[,p,t] <- dynDiscRent_ipt[,p,t] * (( 1 + discRate)^(-t))

    # By area
    dynDiscRent_ip <- apply( X = dynDiscRent_ipt, FUN = sum, MARGIN = c(1,2) )
    dynDiscRent_qp <-  round(apply(X = dynDiscRent_ip, FUN = quantile, probs = qProbs, MARGIN = 2),2)

    for( p in 1:nP )
    {
      # Do species tables first
      for( s in 1:nS )
      {
        specRowIdx <- s + (p-1) * nS

        # Solve for biomass/catch at Emey
        Bmey_sp[s,p] <- getSplineVal( x = E, y = Beq_spe[s,p,], p = Emey_p[p] )
        Ymey_sp[s,p] <- getSplineVal( x = E, y = Yeq_spe[s,p,], p = Emey_p[p] )

        Umey_sp[s,p] <- Ymey_sp[s,p] / Bmey_sp[s,p]

        dynEqBio[specRowIdx,"species"] <- speciesNames[s]
        dynEqCat[specRowIdx,"species"] <- speciesNames[s]
        dynEqHR[specRowIdx,"species"] <- speciesNames[s]
        dynEqHRDist[specRowIdx,"species"] <- speciesNames[s]
        dynEqCatDist[specRowIdx,"species"] <- speciesNames[s]
        dynEqBioDist[specRowIdx,"species"] <- speciesNames[s]

        dynEqBio[specRowIdx,"stock"] <- stockNames[p]
        dynEqCat[specRowIdx,"stock"] <- stockNames[p]
        dynEqHR[specRowIdx,"stock"] <- stockNames[p]
        dynEqHRDist[specRowIdx,"stock"] <- stockNames[p]
        dynEqBioDist[specRowIdx,"stock"] <- stockNames[p]
        dynEqCatDist[specRowIdx,"stock"] <- stockNames[p]

        # Fill static equilibria
        dynEqBio[specRowIdx,"Xmsy"] <- round(Bmsy_sp[s,p],2)
        dynEqBio[specRowIdx,"Xmey"] <- round(Bmey_sp[s,p],2)

        dynEqCat[specRowIdx,"Xmsy"] <- round(MSY_sp[s,p],2)
        dynEqCat[specRowIdx,"Xmey"] <- round(Ymey_sp[s,p],2)

        dynEqHR[specRowIdx,"Xmsy"]  <- round(Umsy_sp[s,p],3)
        dynEqHR[specRowIdx,"Xmey"]  <- round(Umey_sp[s,p],3)


        # Now pull dynamic eqbria
        dynB <- modelList[[k]]$stateDists$SB_ispt[2,s,p]
        dynU <- modelList[[k]]$stateDists$U_ispt[2,s,p]
        dynC <- modelList[[k]]$stateDists$C_ispt[2,s,p]

        dynBpctiles <- modelList[[k]]$stateDists$SB_ispt[c(1,3),s,p]
        dynUpctiles <- modelList[[k]]$stateDists$U_ispt[c(1,3),s,p]
        dynCpctiles <- modelList[[k]]$stateDists$C_ispt[c(1,3),s,p]

        dynBcntrl95 <- paste(round(dynBpctiles,2),collapse = ", ")
        dynCcntrl95 <- paste(round(dynCpctiles,2),collapse = ", ")
        dynUcntrl95 <- paste(round(dynUpctiles,3),collapse = ", ")

        dynEqBio[specRowIdx,colName]  <- round(dynB,2)
        dynEqCat[specRowIdx,colName]  <- round(dynC,2)
        dynEqHR[specRowIdx,colName]   <- round(dynU,3)

        dynEqBioDist[specRowIdx,colName] <- paste0(round(dynB,2), " (", dynBcntrl95, ")")
        dynEqCatDist[specRowIdx,colName] <- paste0(round(dynC,2), " (", dynCcntrl95, ")")
        dynEqHRDist[specRowIdx,colName] <- paste0(round(dynU,3), " (", dynUcntrl95, ")")


      }

      # Then do area tables 
      econYmsy_p[p] <- getSplineVal( x = E, y = econYeq_pe[p,], p = Emsy_p[p])

      # Effort table
      dynEqEff[p,"Xmsy"] <- Emsy_p[p]
      dynEqEff[p,"Xmey"] <- Emey_p[p]

      # Producer surplus
      dynEqPrSp[p,"Xmsy"] <- econYmsy_p[p]
      dynEqPrSp[p,"Xmey"] <- MEY_p[p]

      # Now get dynamic values
      dynPrSp       <- modelList[[k]]$stateDists$profit_ipft[2,p,2]
      dynPrSpctiles <- modelList[[k]]$stateDists$profit_ipft[c(1,3),p,2]
      dynEff        <- modelList[[k]]$stateDists$E_ipft[2,p,2]
      
      dynEqPrSp[p,colName]  <- round(dynPrSp,2)
      dynEqEff[p,colName]   <- round(dynEff,2)

      dynPrSpCntrl95 <- paste(round(dynPrSpctiles,2),collapse = ", ")
      dynEqPrSpDist[p,colName] <- paste0(round(dynPrSp,2), " (", dynPrSpCntrl95, ")")

      # Median and distribution of discounted rent
      dynDiscRent <- dynDiscRent_qp[2,p]
      dynDiscRentPctiles <- dynDiscRent_qp[c(1,3),p]
      dynDiscRentCntrl95 <- paste(dynDiscRentPctiles,collapse = ", ")
      dynEqDiscRent[p,colName] <- dynDiscRent
      dynEqDiscRentDist[p,colName] <- paste0(dynDiscRent, " (", dynDiscRentCntrl95, ")")


    }

    # Now calculate coastwide values for VAR5 stuff
    dynB_ispt   <- modelList[[k]]$modelStates$SB_ispt[,,,distYearIdx]
    dynC_ispt   <- modelList[[k]]$modelStates$C_ispt[,,,distYearIdx]
    dynE_ipt    <- modelList[[k]]$modelStates$E_ipft[,,2,distYearIdx]
    dynRent_ipt <- modelList[[k]]$modelStates$profit_ipft[,,2,distYearIdx]

    dynBcw_it     <- apply(X = dynB_ispt, FUN = sum, MARGIN =c(1,4) )
    dynCcw_it     <- apply(X = dynC_ispt, FUN = sum, MARGIN =c(1,4) )
    dynEcw_it     <- apply(X = dynE_ipt, FUN = sum, MARGIN =c(1,3) )
    dynRentcw_it  <- apply(X = dynRent_ipt, FUN = sum, MARGIN =c(1,3) )

    dynDiscRentcw_i <- apply( X = dynDiscRent_ip, FUN = sum, MARGIN = 1)




    # Now take quantiles
    dynBcw_q        <- round(quantile( dynBcw_it, probs = qProbs ),2)
    dynCcw_q        <- round(quantile( dynCcw_it, probs = qProbs ),2)
    dynEcw_q        <- round(quantile( dynEcw_it, probs = qProbs ),2)
    dynRentcw_q     <- round(quantile( dynRentcw_it, probs = qProbs ),2)
    dynDiscRentcw_q <- round(quantile( dynDiscRentcw_i, probs = qProbs ),2)


    

    # Add medians to median table
    dynEqBio[nS*nP+1,colName] <- dynBcw_q[2]
    dynEqCat[nS*nP+1,colName] <- dynCcw_q[2]

    dynEqEff[nP+1,colName]      <- dynEcw_q[2]
    dynEqPrSp[nP+1,colName]     <- dynRentcw_q[2]
    dynEqDiscRent[nP+1,colName] <- dynDiscRentcw_q[2]

    # Now do the dists again
    dynEqBioDist[nS*nP+1,colName]   <- paste0(dynBcw_q[2]," (", dynBcw_q[1],", ",dynBcw_q[3],")")
    dynEqCatDist[nS*nP+1,colName]   <- paste0(dynCcw_q[2]," (", dynCcw_q[1],", ",dynCcw_q[3],")")
    dynEqPrSpDist[nP+1,colName]     <- paste0(dynRentcw_q[2]," (", dynRentcw_q[1],", ",dynRentcw_q[3],")")
    dynEqDiscRentDist[nP+1,colName] <- paste0(dynDiscRentcw_q[2]," (", dynDiscRentcw_q[1],", ",dynDiscRentcw_q[3],")")

  }

  # Save tables
  outFolder <- here::here("Outputs",groupFolder)
  # Bio
  write.csv(dynEqBio, file = file.path(outFolder,"dynEqBio.csv"))
  write.csv(dynEqCat, file = file.path(outFolder,"dynEqCat.csv"))
  write.csv(dynEqHR, file = file.path(outFolder,"dynEqHR.csv"))
  write.csv(dynEqEff, file = file.path(outFolder,"dynEqEff.csv"))
  write.csv(dynEqPrSp, file = file.path(outFolder,"dynEqPrSp.csv"))
  write.csv(dynEqDiscRent, file = file.path(outFolder,"dynEqDiscRent.csv"))

  write.csv(dynEqBioDist, file = file.path(outFolder,"dynEqBioDist.csv"))
  write.csv(dynEqCatDist, file = file.path(outFolder,"dynEqCatDist.csv"))
  write.csv(dynEqHRDist, file = file.path(outFolder,"dynEqHRDist.csv"))
  write.csv(dynEqPrSpDist, file = file.path(outFolder,"dynEqPrSpDist.csv"))
  write.csv(dynEqDiscRentDist, file = file.path(outFolder,"dynEqDiscRentDist.csv"))

  outList <- list(  dynEqBio  = dynEqBio,
                    dynEqCat  = dynEqCat,
                    dynEqHR   = dynEqHR,
                    dynEqEff  = dynEqEff,
                    dynEqPrSp = dynEqPrSp,
                    dynEqDiscRent = dynEqDiscRent,
                    dynEqBioDist = dynEqBioDist,
                    dynEqCatDist = dynEqCatDist,
                    dynEqHRDist = dynEqHRDist,
                    dynEqPrSpDist = dynEqPrSpDist, 
                    dynEqDiscRentDist = dynEqDiscRentDist)
  
  
  outList
}


# pullModelStates()
# Loads blob pulls catch, effort and biomass 
# time series, and calculates distributions
# of each inside a nominated time period
pullModelStates <- function(  sim         = 1,
                              groupFolder = "omniRuns_priceDevPED",
                              stateVars   = c("C_ispt","SB_ispt","E_ipft"),
                              output      = FALSE,
                              distPeriod  = 2040:2060 )
{
  # First, load the sim that we want
  # to calculate loss for
  simFolder <- .loadSim( sim = sim, folder = groupFolder )

  # Save the blob
  simObj <- blob

  # figure out which reps we want
  goodReps_isp    <- simObj$goodReps_isp
  maxRep_sp       <- apply(X = goodReps_isp, FUN = maxWhich, MARGIN = c(2,3))
  totReps         <- max(maxRep_sp)


  blob <- NULL
  gc()

  # Calculate profit
  rev_ispft      <- simObj$om$Rev_ispft[1:totReps,,,,]
  rev_ipft       <- apply( X = rev_ispft, FUN = sum, MARGIN = c(1,3,4,5))
  effCost_ipft   <- simObj$om$effCost_ipft[1:totReps,,,]
  crewShare      <- simObj$ctlList$opMod$crewShare
  discRate       <- simObj$ctlList$opMod$discountRate

  C_ispt <- simObj$om$C_ispt
  v_ist  <- simObj$om$landVal_ist
  nS <- simObj$om$nS
  nP <- simObj$om$nP

  profit_ipft    <- array(0, dim = dim(rev_ipft))
  profit_ipft    <- rev_ipft - effCost_ipft
  

  fYear <- simObj$ctlList$opMod$fYear

  distPerIdx <- distPeriod - fYear + 1

  # Model dimensions
  tMP <- simObj$om$tMP
  nT  <- simObj$om$nT
  nF  <- simObj$om$nF
  nS  <- simObj$om$nS
  nP  <- simObj$om$nP
  pT  <- simObj$om$pT

  speciesNames  <- c(simObj$om$speciesNames)
  stockNames    <- c(simObj$om$stockNames)
  fleetNames    <- simObj$om$fleetNames

  # We want to calculate loss of given OM states
  # Make a list to hold baseline
  modelStates        <- vector( mode = "list", 
                                length = length(stateVars))
  names(modelStates) <- stateVars

  # Make a list of distributions too
  stateDists         <- modelStates
  stateMeans         <- modelStates
  

  # Pull control list settings 
  ctlList         <- simObj$ctlList

  for( var in stateVars )
  {
    modelStates[[var]] <- simObj$om[[var]]

    # Now calculate distribution
    stateDists[[var]] <- apply( X = simObj$om[[var]][,,,distPerIdx],
                                FUN = quantile,
                                probs = c(0.025, 0.5, 0.975),
                                MARGIN = c(2,3), na.rm = TRUE )

    stateMeans[[var]] <- apply( simObj$om[[var]][,,,distPerIdx],
                                FUN = mean,
                                MARGIN = c(2,3), na.rm = TRUE )
  }

  # Calculate harvest rates
  modelStates$U_ispt <- modelStates$C_ispt / modelStates$SB_ispt
  stateDists$U_ispt <- apply(   X = modelStates$U_ispt[,,,distPerIdx],
                                FUN = quantile,
                                probs = c(0.025, 0.5, 0.975),
                                MARGIN = c(2,3), na.rm = TRUE )

  stateMeans$U_ispt <- apply(   X = modelStates$U_ispt[,,,distPerIdx],
                                FUN = mean,
                                MARGIN = c(2,3), na.rm = TRUE )

  modelStates$profit_ipft <- profit_ipft
  stateDists$profit_ipft  <- apply( X = modelStates$profit_ipft[,,,distPerIdx],
                                    FUN = quantile,
                                    probs = c(0.025, 0.5, 0.975),
                                    MARGIN = c(2,3), na.rm = TRUE )
  stateMeans$profit_ipft <- apply(  X = modelStates$profit_ipft[,,,distPerIdx],
                                    FUN = mean,
                                    MARGIN = c(2,3), na.rm = TRUE )


  outList <- list(  simID         = simObj$folder,
                    goodReps_isp  = goodReps_isp,
                    speciesNames  = speciesNames,
                    stockNames    = stockNames,
                    fYear         = fYear,
                    modelStates   = modelStates,
                    stateDists    = stateDists,
                    stateMeans    = stateMeans,
                    B0_sp         = simObj$om$B0_sp,
                    rp            = simObj$rp[[1]],
                    tMP           = tMP,
                    nT            = nT, 
                    nF            = nF, 
                    nS            = nS, 
                    nP            = nP, 
                    pT            = pT,
                    simFolder     = simFolder  )


  # Save loss to sim folder
  save( outList, file = file.path(simFolder,paste(simObj$simLabel,"_modelStates.RData",sep = "") ) )

  if(output)
    return(outList)
} # END calcLoss()


# calcProbOverfished()
# Calculates mean and range of probability of being
# overfished wrt single species reference points
calcProbOverfished <- function( groupFolder="DERTACS_reruns_Oct10",
                                prefix = "parBat" )
{
  # First get info files so we can load the right 
  # loss objects
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  simFolderList <- simFolderList[grepl(prefix, simFolderList)]

  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
              filter(grepl(prefix,simLabel)) %>%
              mutate( perfPath = here::here("Outputs",groupFolder,simLabel,"simPerfStats.csv"),
                      path = here::here("Outputs",groupFolder,simLabel),
                      pBtLt.4Bmsy = NA,
                      minProbBtLt.4Bmsy = NA,
                      maxProbBtLt.4Bmsy = NA,
                      pBtLt.8Bmsy = NA,
                      minProbBtLt.8Bmsy = NA,
                      maxProbBtLt.8Bmsy = NA,
                      pFtGtFmsy = NA,
                      minProbFtGtFmsy = NA,
                      maxProbFtGtFmsy = NA )

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

    unlist(outList)
  }

  MPfactors <- lapply( X = info.df$mp, FUN = splitMP )
  MPfactors <- do.call(rbind, MPfactors) %>% as.data.frame()

  # Append MP factors
  info.df <- cbind( info.df, MPfactors )

  # Read in the performance tables
  
  perfTables <- lapply(X = info.df$perfPath, FUN = read.csv)
  
  for( k in 1:length(perfTables) )
  {
    perfTables[[k]] <- perfTables[[k]] %>% filter( ! species %in% "Rock" )
    info.df$pBtLt.4Bmsy[k]        <- round(mean(1 - perfTables[[k]]$pBtGt.4Bmsy ),2)
    info.df$minProbBtLt.4Bmsy[k]  <- round(min(1 - perfTables[[k]]$pBtGt.4Bmsy ),2)
    info.df$maxProbBtLt.4Bmsy[k]  <- round(max(1 - perfTables[[k]]$pBtGt.4Bmsy ),2)
    info.df$pBtLt.8Bmsy[k]        <- round(mean(1 - perfTables[[k]]$PBtGt.8Bmsy ),2)
    info.df$minProbBtLt.8Bmsy[k]  <- round(min(1 - perfTables[[k]]$PBtGt.8Bmsy ),2)
    info.df$maxProbBtLt.8Bmsy[k]  <- round(max(1 - perfTables[[k]]$PBtGt.8Bmsy ),2)
    info.df$pFtGtFmsy[k]          <- round(mean(perfTables[[k]]$pFtGtFmsy ),2)
    info.df$minProbFtGtFmsy[k]    <- round(min(perfTables[[k]]$pFtGtFmsy ),2)
    info.df$maxProbFtGtFmsy[k]    <- round(max(perfTables[[k]]$pFtGtFmsy ),2)
  }

  probOverfishedTable <- info.df %>%
                          dplyr::select(  scenario, 
                                          AM, 
                                          pBtLt.4Bmsy,
                                          minProbBtLt.4Bmsy,
                                          maxProbBtLt.4Bmsy, 
                                          pBtLt.8Bmsy,
                                          minProbBtLt.8Bmsy,
                                          maxProbBtLt.8Bmsy,
                                          pFtGtFmsy,
                                          minProbFtGtFmsy,
                                          maxProbFtGtFmsy )


  probOverfishedTable
} # END calcProbOverfished


calc_AssErrDist <- function(  groupFolder="DERTACs_reruns_Jan5",
                              prefix = "parBat",
                              scenName = "DERfit_AsSsIdx",
                              mpName = "hierMultiStock_omF_MSrefPts" )
{
  # First get info files so we can load the right 
  # loss objects
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  simFolderList <- simFolderList[grepl(prefix, simFolderList)]

  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
              filter( grepl(prefix,simLabel),
                      scenario %in% scenName,
                      mp %in% mpName )

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

    unlist(outList)
  }

  MPfactors <- lapply( X = info.df$mp, FUN = splitMP )
  MPfactors <- do.call(rbind, MPfactors) %>% as.data.frame()

  # Read in the performance tables
  mpTable <- cbind( info.df, MPfactors ) %>%
              mutate( perfPath = here::here("Outputs",groupFolder,simLabel,"simPerfStats.csv"),
                      path = here::here("Outputs",groupFolder,simLabel) )

  blobFiles <- file.path(mpTable$path,paste(mpTable$simLabel,".RData",sep = ""))
  names(blobFiles) <- mpTable$simLabel
  
  nSims <- length(blobFiles)
  retroEstList <- vector(mode ="list", length = nSims)
  names(retroEstList) <- mpTable$simLabel

  errList <- vector(mode = "list", length = length(retroEstList))
  names(errList) <- names(retroEstList)



  for( i in 1:nSims )
  {
    # load blob
    simID <- mpTable$simLabel[i]
    load(blobFiles[[simID]])

    # Save info
    retroEstList[[simID]]$SB_ispt        <- blob$om$SB_ispt
    retroEstList[[simID]]$retroSB_itspt  <- blob$mp$assess$retroSB_itspt


    if( i == 1 )
    {
      pT      <- blob$ctlList$opMod$pT
      nS      <- blob$om$nS
      nP      <- blob$om$nP
      nT      <- blob$om$nT
      nF      <- blob$om$nF
      tMP     <- blob$om$tMP
      nReps   <- dim(retroEstList[[i]]$SB_ispt)[1]
      species <- blob$om$speciesNames
      stock   <- blob$om$stockNames

      assErr_ispt <- array(NA, dim = c(nReps,nS,nP,nT) )
    }

    blob <- NULL
    gc()

    errList[[simID]]$assErr_ispt <- assErr_ispt

    
    for( t in tMP:nT )
    {
      pTdx <- t - tMP + 1
      errList[[simID]]$assErr_ispt[,,,t] <-  log(retroEstList[[simID]]$retroSB_itspt[,pTdx,,,t]/retroEstList[[simID]]$SB_ispt[,,,t])
    }

    # Now, we want to remove the nSxnP structure so we can calculate correlation
    # and covariance
    errList[[simID]]$assErr_itSP <- array(NA, dim = c(nReps,pT, nS*nP) )
    errList[[simID]]$cov_i       <- array(NA, dim = c(nReps,nS*nP,nS*nP))
    errList[[simID]]$cor_i       <- array(NA, dim = c(nReps,nS*nP,nS*nP))
    errList[[simID]]$mu_i        <- array(NA, dim = c(nReps,nS*nP))
    errList[[simID]]$sd_i        <- array(NA, dim = c(nReps,nS*nP))
    for( i in 1:nReps )
    {
      for( s in 1:nS )
        for( p in 1:nP )
        {
          spIdx <- (s-1) * nP + p
          errList[[simID]]$assErr_itSP[i,1:pT,spIdx] <- errList[[simID]]$assErr_ispt[i,s,p,tMP:nT]
        }

      errList[[simID]]$cov_i[i,,] <- cov(errList[[simID]]$assErr_itSP[i,,])
      errList[[simID]]$cor_i[i,,] <- cor(errList[[simID]]$assErr_itSP[i,,])
      errList[[simID]]$mu_i[i,]   <- apply(X = errList[[simID]]$assErr_itSP[i,,], FUN = mean, MARGIN = 2)
      errList[[simID]]$sd_i[i,]   <- apply(X = errList[[simID]]$assErr_itSP[i,,], FUN = sd, MARGIN = 2)

    }

    # Now look at distribution of correlations, covariance and mean errors
    errList[[simID]]$cov_q <- apply(  X = errList[[simID]]$cov_i, FUN = quantile, 
                                      probs = c(0.025, 0.5, 0.975),
                                      MARGIN = c(2,3) )
    errList[[simID]]$cor_q <- apply(  X = errList[[simID]]$cor_i, FUN = quantile, 
                                      probs = c(0.025, 0.5, 0.975),
                                      MARGIN = c(2,3) )
    errList[[simID]]$mu_q <- apply(  X = errList[[simID]]$mu_i, FUN = quantile, 
                                      probs = c(0.025, 0.5, 0.975),
                                      MARGIN = c(2) )
    errList[[simID]]$sd_q <- apply(  X = errList[[simID]]$sd_i, FUN = quantile, 
                                      probs = c(0.025, 0.5, 0.975),
                                      MARGIN = c(2) )
  }

  # Save assErr list
  save(errList, file = "assErrDists.RData")

  errList

}

calcMSE_AMestimates <- function(  groupFolder="DERTACS_Reruns_Jan5",
                                  prefix = "parBat" )
{
  # First get info files so we can load the right 
  # loss objects
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  simFolderList <- simFolderList[grepl(prefix, simFolderList)]

  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
              filter(grepl(prefix,simLabel))

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

    unlist(outList)
  }

  MPfactors <- lapply( X = info.df$mp, FUN = splitMP )
  MPfactors <- do.call(rbind, MPfactors) %>% as.data.frame()

  # Read in the performance tables
  mpTable <- cbind( info.df, MPfactors ) %>%
              mutate( perfPath = here::here("Outputs",groupFolder,simLabel,"simPerfStats.csv"),
                      path = here::here("Outputs",groupFolder,simLabel) )

  blobFiles <- file.path(mpTable$path,paste(mpTable$simLabel,".RData",sep = ""))
  names(blobFiles) <- mpTable$simLabel
  
  nSims <- length(blobFiles)
  retroEstList <- vector(mode ="list", length = nSims)
  names(retroEstList) <- mpTable$simLabel

  for( i in 1:nSims )
  {
    # load blob
    simID <- mpTable$simLabel[i]
    load(blobFiles[[simID]])

    # Save info
    retroEstList[[simID]]$SB_ispt        <- blob$om$SB_ispt
    retroEstList[[simID]]$retroSB_itspt  <- blob$mp$assess$retroSB_itspt
    retroEstList[[simID]]$retroUmsy_itsp <- blob$mp$assess$retroUmsy_itsp
    retroEstList[[simID]]$retroBmsy_itsp <- blob$mp$assess$retroBmsy_itsp
    retroEstList[[simID]]$retroMSY_itsp  <- blob$mp$assess$retroBmsy_itsp * blob$mp$assess$retroUmsy_itsp
    retroEstList[[simID]]$Bmsy_sp        <- blob$rp[[1]]$FmsyRefPts$BeqFmsy_sp
    retroEstList[[simID]]$YeqFmsy_sp     <- blob$rp[[1]]$FmsyRefPts$YeqFmsy_sp
    retroEstList[[simID]]$Umsy_sp        <- retroEstList[[simID]]$YeqFmsy_sp / retroEstList[[simID]]$Bmsy_sp


    if( i == 1 )
    {
      pT      <- blob$ctlList$opMod$pT
      nS      <- blob$om$nS
      nP      <- blob$om$nP
      nT      <- blob$om$nT
      nF      <- blob$om$nF
      tMP     <- blob$om$tMP
      nReps   <- dim(retroEstList[[i]]$SB_ispt)[1]
      species <- blob$om$speciesNames
      stock   <- blob$om$stockNames
    }

    blob <- NULL
    gc()
  }


  
  # No go through the list and calculate errors
  errList <- vector(mode = "list", length = nSims )
  names(errList) <- names(retroEstList)
  tabColNames <- c( "simID",
                    "Scenario","AM",
                    "Species","Stock",
                    "MARE_Bt","MARE_Bmsy","MARE_Umsy","MARE_MSY")
  
  mareTable <- matrix(NA, nrow = nSims * nS * nP, ncol = length(tabColNames))
  
  colnames(mareTable) <- tabColNames
  
  mareTable <- as.data.frame(mareTable)

  # browser()

  rowIdx <- 1
  for( i in 1:nSims )
  {
    simID <- mpTable$simLabel[i]
    errList[[simID]]$mseSB_itsp   <- array(NA, dim = c(nReps,pT,nS,nP))
    errList[[simID]]$mseBmsy_itsp <- array(NA, dim = c(nReps,pT,nS,nP))
    errList[[simID]]$mseUmsy_itsp <- array(NA, dim = c(nReps,pT,nS,nP))
    errList[[simID]]$mseMSY_itsp  <- array(NA, dim = c(nReps,pT,nS,nP))

    errList[[simID]]$mseSB_sp   <- array(NA, dim = c(nS,nP))
    errList[[simID]]$mseBmsy_sp <- array(NA, dim = c(nS,nP))
    errList[[simID]]$mseUmsy_sp <- array(NA, dim = c(nS,nP))
    errList[[simID]]$mseMSY_sp  <- array(NA, dim = c(nS,nP))



    retroSB_itspt   <- retroEstList[[simID]]$retroSB_itspt
    omSB_ispt       <- retroEstList[[simID]]$SB_ispt
    retroUmsy_itsp  <- retroEstList[[simID]]$retroUmsy_itsp
    retroBmsy_itsp  <- retroEstList[[simID]]$retroBmsy_itsp
    retroMSY_itsp   <- retroEstList[[simID]]$retroMSY_itsp
    omUmsy_sp       <- retroEstList[[simID]]$Umsy_sp
    omBmsy_sp       <- retroEstList[[simID]]$Bmsy_sp
    omYeqFmsy_sp    <- retroEstList[[simID]]$YeqFmsy_sp

    nSS <-  nS
    nPP <-  nP
    for( tt in 1:pT )
    {
      errBmsy_isp    <- array(NA, dim = c(nReps,nS,nP))
      errUmsy_isp    <- array(NA, dim = c(nReps,nS,nP))
      errMSY_isp    <- array(NA, dim = c(nReps,nS,nP))
      errSB_ispt     <- array(NA, dim = c(nReps,nS,nP,nT))

      if( mpTable[i,"AM"] %in% c("speciesPooling","totalAgg") )
      {
        omSB_ispt_sum     <- apply( X = omSB_ispt, FUN = sum, MARGIN = c(1,3,4) )
        omBmsy_sp_sum     <- apply( X = omBmsy_sp, FUN = sum, MARGIN = c(2) )
        omYeqFmsy_sp_sum  <- apply( X = omYeqFmsy_sp, FUN = sum, MARGIN = c(2) )
        omUmsy_sp_sum     <- omYeqFmsy_sp_sum / omBmsy_sp_sum

        omSB_ispt <- omSB_ispt[,1,,,drop = FALSE]
        omBmsy_sp <- omBmsy_sp[1,,drop = FALSE]
        omUmsy_sp <- omUmsy_sp[1,,drop = FALSE]
        omYeqFmsy_sp <- omYeqFmsy_sp[1,,drop = FALSE]

        omSB_ispt[,1,,] <- omSB_ispt_sum
        omBmsy_sp[1,]   <- omBmsy_sp_sum
        omUmsy_sp[1,]   <- omUmsy_sp_sum
        omYeqFmsy_sp[1,]<- omYeqFmsy_sp_sum

        
        nSS <- 1
      }

      if( mpTable[i,"AM"] %in% c("spatialPooling","totalAgg") )
      {

        omSB_ispt_sum     <- apply( X = omSB_ispt, FUN = sum, MARGIN = c(1,2,4) )
        omBmsy_sp_sum     <- apply( X = omBmsy_sp, FUN = sum, MARGIN = c(1) )
        omYeqFmsy_sp_sum  <- apply( X = omYeqFmsy_sp, FUN = sum, MARGIN = c(1) )
        omUmsy_sp_sum     <- omYeqFmsy_sp_sum / omBmsy_sp_sum

        omSB_ispt <- omSB_ispt[,,1,,drop = FALSE]
        omBmsy_sp <- omBmsy_sp[,1,drop = FALSE]
        omUmsy_sp <- omUmsy_sp[,1,drop = FALSE]
        omYeqFmsy_sp <- omYeqFmsy_sp[,1,drop = FALSE]

        omSB_ispt[,,1,] <- omSB_ispt_sum
        omBmsy_sp[,1]   <- omBmsy_sp_sum
        omUmsy_sp[,1]   <- omUmsy_sp_sum
        omYeqFmsy_sp[,1]<- omYeqFmsy_sp_sum

        nPP <- 1
      }

      errSB_ispt[,1:nSS,1:nPP,]    <- abs(retroSB_itspt[,tt,1:nSS,1:nPP,] - omSB_ispt[,1:nSS,1:nPP,])/omSB_ispt[,1:nSS,1:nPP,]

      # browser()

      # Loop and calculate error
      for( j in 1:nReps )
      {

        errBmsy_isp[j,1:nSS,1:nPP]    <- abs(retroBmsy_itsp[j,tt,1:nSS,1:nPP] - omBmsy_sp[1:nSS,1:nPP])/omBmsy_sp[1:nSS,1:nPP]
        errUmsy_isp[j,1:nSS,1:nPP]    <- abs(retroUmsy_itsp[j,tt,1:nSS,1:nPP] - omUmsy_sp[1:nSS,1:nPP])/omUmsy_sp[1:nSS,1:nPP]
        errMSY_isp[j,1:nSS,1:nPP]     <- abs(retroMSY_itsp[j,tt,1:nSS,1:nPP] - omYeqFmsy_sp[1:nSS,1:nPP])/omYeqFmsy_sp[1:nSS,1:nPP]
      }

      errList[[simID]]$mseSB_itsp[,tt,1:nSS,1:nPP] <- apply(X = errSB_ispt[,1:nSS,1:nPP,,drop = FALSE], FUN = median, MARGIN = c(1,2,3),na.rm = T )
      errList[[simID]]$mseBmsy_itsp[,tt,1:nSS,1:nPP] <- errBmsy_isp[,1:nSS,1:nPP]
      errList[[simID]]$mseUmsy_itsp[,tt,1:nSS,1:nPP] <- errUmsy_isp[,1:nSS,1:nPP]
      errList[[simID]]$mseMSY_itsp[,tt,1:nSS,1:nPP] <- errMSY_isp[,1:nSS,1:nPP]

    }

    errList[[simID]]$mseSB_sp[1:nSS,1:nPP]   <- apply( X = errList[[simID]]$mseSB_itsp[,,1:nSS,1:nPP,drop = FALSE], FUN = median, MARGIN = c(3,4),na.rm = T)
    errList[[simID]]$mseBmsy_sp[1:nSS,1:nPP] <- apply( X = errList[[simID]]$mseBmsy_itsp[,,1:nSS,1:nPP,drop = FALSE], FUN = median, MARGIN = c(3,4),na.rm = T)
    errList[[simID]]$mseUmsy_sp[1:nSS,1:nPP] <- apply( X = errList[[simID]]$mseUmsy_itsp[,,1:nSS,1:nPP,drop = FALSE], FUN = median, MARGIN = c(3,4),na.rm = T)
    errList[[simID]]$mseMSY_sp[1:nSS,1:nPP]  <- apply( X = errList[[simID]]$mseMSY_itsp[,,1:nSS,1:nPP,drop = FALSE], FUN = median, MARGIN = c(3,4),na.rm = T)

    

    for( s in 1:nSS )
    {
      specLab   <- species[s]

      if( mpTable[i,"AM"] %in% c("speciesPooling","totalAgg") )
        specLab <- "SpeciesPooled"
      
      for( p in 1:nPP )
      {
        mareTable[rowIdx,"simID"]     <- simID
        mareTable[rowIdx,"Scenario"]  <- mpTable[i,"scenario"]
        mareTable[rowIdx,"AM"]        <- mpTable[i,"AM"]

        stockLab <- stock[p]

        if( mpTable[i,"AM"] %in% c("spatialPooling","totalAgg") )
          stockLab <- "spatialPooled"

        mareTable[rowIdx,"Species"]         <- specLab
        mareTable[rowIdx,"Stock"]           <- stockLab
        mareTable[rowIdx,"MARE_Bt"]         <- round(errList[[simID]]$mseSB_sp[s,p],3)
        mareTable[rowIdx,"MARE_Bmsy"]       <- round(errList[[simID]]$mseBmsy_sp[s,p],3)
        mareTable[rowIdx,"MARE_Umsy"]       <- round(errList[[simID]]$mseUmsy_sp[s,p],3)
        mareTable[rowIdx,"MARE_MSY"]        <- round(errList[[simID]]$mseMSY_sp[s,p],3)

        rowIdx <- rowIdx + 1
      }
    }
    
  }
  # reduce to only the rows we used
  mareTable <- mareTable[1:(rowIdx-1),] %>%
                dplyr::arrange( Scenario, AM, Species, Stock )
  # Write to groupFolder
  write.csv(mareTable, file = here::here("Outputs",groupFolder,"mareTable.csv"))

  mareTable

}

maxWhich <- function( vector )
{
  idx <- max(which(vector))

  idx
}

# makePairedLossRankTable()
# Runs calcLossRank for abs and rel bio and catch
# then summarises the mean (min, max) rank for all
# four into a single table
makePairedLossRankTable <- function(  groupFolder = "DERTACS_reruns_sep24",
                                      prefix = "parBat",
                                      period = 73:82,
                                      clearBadReps = TRUE,
                                      minSampSize = 100,
                                      nGoodReps = 100,
                                      dim1 = 1:3,
                                      dim2 = 1:3,
                                      var = "C_ispt",
                                      lossType = "abs",
                                      AMlabs = rev(c( SS = "singleStock",
                                                      HMS = "hierMultiStock",
                                                      SpecPool = "speciesPooling",
                                                      SpatPool = "spatialPooling",
                                                      TotPool = "totalAgg" ) ),
                                      scenLabs = c( High  = "DERfit_HcMcAsSsIdx",
                                                    Mod  = "DERfit_McAsSsIdx",
                                                    Low = "DERfit_AsSsIdx" ) )
{
  # First, calculate the paired loss ranks
  rankArray_SAisp <- calcPairedLossRank(  groupFolder = groupFolder,
                                          lossVar = var,
                                          lossType = lossType,
                                          prefix = prefix,
                                          period = period,
                                          clearBadReps = clearBadReps,
                                          minSampSize = minSampSize )
  goodReps_SAisp <- attr(rankArray_SAisp,"conv")[scenLabs,,,,,drop = FALSE]
  
  # Cut down to scenarios we want
  rankArray_SAisp <- rankArray_SAisp[scenLabs,,,,,drop = FALSE]
  goodReps_SAisp <- goodReps_SAisp[scenLabs,,,,,drop = FALSE]

  nScen <- length(scenLabs)
  nAM   <- length(AMlabs)
  nS    <- length(dim1)
  nP    <- length(dim2)


  if( clearBadReps )
  {
    nGood_SAsp <- apply( X = goodReps_SAisp, FUN = sum, MARGIN = c(1,2,4,5),
                          na.rm = T )

    badRuns_arr.ind <-  which(nGood_SAsp < minSampSize, arr.ind = TRUE )

    if( nrow(badRuns_arr.ind) > 0 )
      for( j in 1:nrow(badRuns_arr.ind))
      {
        ind <- badRuns_arr.ind[j,]
        # Remove the species and stock results from all AMs - replace goodReps with TRUE
        # as the NAs in rank will be removed later
        goodReps_SAisp[,,,ind[3],ind[4]] <- TRUE
        rankArray_SAisp[,,,ind[3],ind[4]] <- NA
      }

    # Ok, now  we want to calculate the largest number of reps we have available
    goodReps_SAi <- apply(X = goodReps_SAisp, FUN = prod, MARGIN = c(1,2,3), na.rm = T)
    goodReps_Si <- apply(X = goodReps_SAi, FUN = prod, MARGIN = c(1,3))

    goodReps_i <- apply(X = goodReps_SAisp, FUN = prod, MARGIN = 3)

    goodIdx <- which(goodReps_i == 1)

    if( length(goodIdx) > nGoodReps)
      goodIdx <- goodIdx[1:nGoodReps]

    
    goodReps_SAisp  <- goodReps_SAisp[,,goodIdx,,]
    goodReps_SAi    <- goodReps_SAi[,,goodIdx]
    rankArray_SAisp <- rankArray_SAisp[,,goodIdx,,]

  }

  dimnames(rankArray_SAisp)[[1]] <- names(scenLabs)


  rank.df <- reshape2::melt(rankArray_SAisp) %>%
              filter(!is.na(value)) %>%
              rename(rank = value )
              

  rankSummary.df <- rank.df %>%
                    group_by(scenario, AM, rank) %>%
                    summarise( nRuns = n() ) %>%
                    mutate( rankProp = nRuns / sum(nRuns) ) %>%
                    summarise( modalRank = which.max(rankProp),
                                avgRank = sum(rankProp * rank) ) %>%
                    ungroup()

  tabList <- list(  rank.df = rank.df,
                    rankSummary.df = rankSummary.df )

  fullRankFilename <- paste("fullRankTable_",lossType,var,".csv", sep = "")
  summRankFilename <- paste("summRankTable_",lossType,var,".csv", sep = "")
  write.csv( rank.df, file = file.path(here::here("Outputs",groupFolder,fullRankFilename)))
  write.csv( rankSummary.df, file = file.path(here::here("Outputs",groupFolder,summRankFilename)))

  return(tabList)
} # END makePairedLossRankTable()

# makeLossRankTable()
# Runs calcLossRank for abs and rel bio and catch
# then summarises the mean (min, max) rank for all
# four into a single table
makeLossRankTable <- function(  groupFolder = "DERTACS_reruns_sep24",
                                prefix = "parBat",
                                period = 73:82,
                                clearBadReps = TRUE,
                                minSampSize = 10 )
{

  # Make absC rank tables
  absCrank <- calcLossRank( groupFolder = groupFolder,
                            lossVar = "C_ispt",
                            lossType = "abs",
                            prefix = prefix,
                            period = period,
                            clearBadReps = clearBadReps,
                            minSampSize = minSampSize )$summRank %>%
                mutate( absCrank = paste( meanRank, " (",minRank, ", ", maxRank, ")", sep = "")) 

  # # Rel C
  relCrank <- calcLossRank( groupFolder = groupFolder,
                            lossVar = "C_ispt",
                            lossType = "rel",
                            prefix = prefix,
                            period = period,
                            clearBadReps = clearBadReps,
                            minSampSize = minSampSize )$summRank %>%
                mutate( relCrank = paste( meanRank, " (",minRank, ", ", maxRank, ")", sep = "")) 

  # Make absC rank tables
  absBrank <- calcLossRank( groupFolder = groupFolder,
                            lossVar = "SB_ispt",
                            lossType = "abs",
                            prefix = prefix,
                            period = period,
                            clearBadReps = clearBadReps,
                            minSampSize = minSampSize )$summRank %>%
                mutate( absBrank = paste( meanRank, " (",minRank, ", ", maxRank, ")", sep = ""))

  # Rel C
  relBrank <- calcLossRank( groupFolder = groupFolder,
                            lossVar = "SB_ispt",
                            lossType = "rel",
                            prefix = prefix,
                            period = period,
                            clearBadReps = clearBadReps,
                            minSampSize = minSampSize )$summRank %>%
                mutate( relBrank = paste( meanRank, " (",minRank, ", ", maxRank, ")", sep = "")) 


  # Join into one table
  rankTable <- absCrank

  write.csv(absCrank, file = file.path(here::here("Outputs",groupFolder,"absCrankTable.csv")))
  write.csv(relCrank, file = file.path(here::here("Outputs",groupFolder,"relCrankTable.csv")))
  write.csv(absBrank, file = file.path(here::here("Outputs",groupFolder,"absBrankTable.csv")))
  write.csv(relBrank, file = file.path(here::here("Outputs",groupFolder,"relBrankTable.csv")))

  return(rankTable)
} # END makeLossRankTable()

# makeCumulativeLossArray_SAisp()
# Refactored from other functions where I kept copying
# it. Loads loss objects from single sims
# and arranges it into an array for better plotting.
makeCumulativeLossArray_SAisp <- function(  groupFolder = "DERTACs_reruns_.75tau_.05PE_180reps",
                                            prefix = "parBat",
                                            lossType = "abs",
                                            var = "C_ispt",
                                            period = 73:82,
                                            lossList = NULL,
                                            refPts = NULL,
                                            AMlabs = c( SS = "singleStock",
                                                        HMS = "hierMultiStock",
                                                        SpePool = "speciesPooling",
                                                        SpaPool = "spatialPooling",
                                                        TA = "totalAgg" ),
                                            scenLabs = c( Comm1  = "DERfit_HcMcAsSsIdx",
                                                          Comm2  = "DERfit_McAsSsIdx",
                                                          Comm3  = "DERfit_McSsIdx",
                                                          Surv = "DERfit_AsSsIdx" ),
                                            clearBadReps = FALSE
                                        )
{
  # First get info files so we can load the right 
  # loss objects
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  simFolderList <- simFolderList[grepl(prefix, simFolderList)]

  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
              filter(grepl(prefix,simLabel))

  if( lossType == "rel" )
  {
    lossArrayName <- "totRelLoss_isp"
  }

  if( lossType == "abs" )
  {
    lossArrayName <- "totAbsLoss_isp"
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

    unlist(outList)
  }

  MPfactors <- lapply( X = info.df$mp, FUN = splitMP )
  MPfactors <- do.call(rbind, MPfactors) %>% as.data.frame()

  # Read in the performance tables
  mpTable <- cbind( info.df, MPfactors ) %>%
              mutate( perfPath = here::here("Outputs",groupFolder,simLabel,"simPerfStats.csv"),
                      path = here::here("Outputs",groupFolder,simLabel) )

  blobFiles <- file.path(mpTable$path,paste(mpTable$simLabel,".RData",sep = ""))
  lossFiles <- file.path(mpTable$path,"loss.RData")

  lossList <- lapply(X = mpTable$simLabel, FUN = .loadLoss, folder=  groupFolder)

  names(lossList) <- mpTable$simLabel

  # Calculate total loss for variable/period
  totLossList  <- lapply( X = lossList,
                          FUN = calcTotalLossPeriod,
                          var = var, period = period )
  names(totLossList) <- names(lossList)

  # Now pull dimensions from the blob
  nT  <- lossList[[1]]$nT
  nS  <- lossList[[1]]$nS
  nP  <- lossList[[1]]$nP
  pT  <- dim(lossList[[1]]$retroSB_itspt)[2]
  tMP <- lossList[[1]]$tMP

  speciesNames <- lossList[[1]]$speciesNames[1:3]
  stockNames   <- lossList[[1]]$stockNames[1:3]


  # Get number of MPs
  nSims <- length(lossList)

  AMs       <- unique(mpTable$AM)
  Fsources  <- unique(mpTable$Fsrce)
  eqbm      <- unique(mpTable$eqbm)
  MPs       <- unique(info.df$mp)
  scenarios <- unique(info.df$scenario)

  nAM       <- length(AMs)
  nSrce     <- length(Fsources)
  nEqbm     <- length(eqbm)
  nMP       <- length(MPs)
  nScen     <- length(scenarios)

  nReps_k <- numeric(length = nSims)
  goodRepsList <- vector(mode = "list", length = nSims)

  for( k in 1:nSims )
  {
    simID <- info.df$simLabel[k]
    simLoss <- totLossList[[simID]][[lossArrayName]] 
    nReps_k[k] <- dim(simLoss)[1]
    goodRepsList[[k]]$goodReps_isp <- lossList[[k]]$goodReps_isp
    goodRepsList[[k]]$nReps_sp     <- apply(X = lossList[[k]]$goodReps_isp, FUN = sum, MARGIN = c(2,3) )
    names(goodRepsList)[k] <- simID
  }

  nReps         <- max(nReps_k)
  speciesNames  <- lossList[[1]]$speciesNames
  stockNames    <- lossList[[1]]$stockNames

  # Make an array to hold loss function values
  totLossArray_SAisp <- array(NA,  dim = c(nScen, nAM, nReps, nS+1, nP+1 ),
                                  dimnames = list(  scenario = scenarios,
                                                    AM = AMs,
                                                    rep = 1:nReps,
                                                    species = speciesNames,
                                                    stock = stockNames ) )

  for( k in 1:nSims )
  {
    simID <- mpTable$simLabel[k]
    simLoss <- totLossList[[simID]][[lossArrayName]] 

    scenID  <- mpTable$scenario[k]
    amID    <- mpTable$AM[k]


    totLossArray_SAisp[scenID,amID,1:dim(simLoss)[1],,] <- simLoss[1:dim(simLoss)[1],,]

    if(clearBadReps)
    {
      goodReps_isp <- goodRepsList[[k]]$goodReps_isp[1:nReps_k[k],,]
      totLossArray_SAisp[scenID,amID,,,][!goodReps_isp] <- NA

      # remove loss when nReps is below min sample size
      for(s in 1:nS)
        for( p in 1:nP )
          if( goodRepsList[[k]]$nReps_sp[s,p] < minSampSize )
            totLossArray_SAisp[scenID,amID,,s,p] <- NA            
    }
  }

  # Return stuff
  outList <- list(  cumLossArray_SAisp = totLossArray_SAisp,
                    info.df = mpTable,
                    lossList = lossList )
  return( outList )
} # END makeCumulativeLoss_SAisp()


# calcLossRank()
# Calculates each AMs rank with respect
# to abs/rel biomass/catch loss. Outputs
# a table with scenarios as rows and
# AMs as columns
calcLossRank <- function( groupFolder = "DERTACS_reruns_sep24",
                          lossVar = "SB_ispt",
                          lossType = "rel",
                          prefix = "parBat",
                          period = c(73:82),
                          clearBadReps = TRUE,
                          minSampSize = 10 )
{
  # First get info files so we can load the right 
  # loss objects
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  simFolderList <- simFolderList[grepl(prefix, simFolderList)]

  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
              filter(grepl(prefix,simLabel))

  if( lossType == "rel" )
  {
    lossArrayName <- "totRelLoss_isp"
  }

  if( lossType == "abs" )
  {
    lossArrayName <- "totAbsLoss_isp"
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

    unlist(outList)
  }

  MPfactors <- lapply( X = info.df$mp, FUN = splitMP )
  MPfactors <- do.call(rbind, MPfactors) %>% as.data.frame()

  # Read in the performance tables
  mpTable <- cbind( info.df, MPfactors ) %>%
              mutate( perfPath = here::here("Outputs",groupFolder,simLabel,"simPerfStats.csv"),
                      path = here::here("Outputs",groupFolder,simLabel) )

  blobFiles <- file.path(mpTable$path,paste(mpTable$simLabel,".RData",sep = ""))
  lossFiles <- file.path(mpTable$path,"loss.RData")

  lossList <- lapply(X = mpTable$simLabel, FUN = .loadLoss, folder=  groupFolder)

  names(lossList) <- mpTable$simLabel

  # Calculate total loss for variable/period
  totLossList  <- lapply( X = lossList,
                          FUN = calcTotalLossPeriod,
                          var = lossVar, period = period )
  names(totLossList) <- names(lossList)

  # Now pull dimensions from the blob
  nT  <- lossList[[1]]$nT
  nS  <- lossList[[1]]$nS
  nP  <- lossList[[1]]$nP
  pT  <- dim(lossList[[1]]$retroSB_itspt)[2]
  tMP <- lossList[[1]]$tMP

  speciesNames <- lossList[[1]]$speciesNames[1:3]
  stockNames   <- lossList[[1]]$stockNames[1:3]


  # Get number of MPs
  nSims <- length(lossList)

  AMs       <- unique(mpTable$AM)
  Fsources  <- unique(mpTable$Fsrce)
  eqbm      <- unique(mpTable$eqbm)
  MPs       <- unique(info.df$mp)
  scenarios <- unique(info.df$scenario)

  nAM       <- length(AMs)
  nSrce     <- length(Fsources)
  nEqbm     <- length(eqbm)
  nMP       <- length(MPs)
  nScen     <- length(scenarios)

  nReps_k <- numeric(length = nSims)
  goodRepsList <- vector(mode = "list", length = nSims)

  for( k in 1:nSims )
  {
    simID <- info.df$simLabel[k]
    simLoss <- totLossList[[simID]][[lossArrayName]] 
    nReps_k[k] <- dim(simLoss)[1]
    goodRepsList[[k]]$goodReps_isp <- lossList[[k]]$goodReps_isp
    goodRepsList[[k]]$nReps_sp     <- apply(X = lossList[[k]]$goodReps_isp, FUN = sum, MARGIN = c(2,3) )
    names(goodRepsList)[k] <- simID
  }

  nReps         <- max(nReps_k)
  speciesNames  <- lossList[[1]]$speciesNames
  stockNames    <- lossList[[1]]$stockNames

  # Make an array to hold loss function values
  totLossArray_SAisp <- array(NA,  dim = c(nScen, nAM, nReps, nS+1, nP+1 ),
                                  dimnames = list(  scenario = scenarios,
                                                    AM = AMs,
                                                    rep = 1:nReps,
                                                    species = speciesNames,
                                                    stock = stockNames ) )

  for( k in 1:nSims )
  {
    simID <- mpTable$simLabel[k]
    simLoss <- totLossList[[simID]][[lossArrayName]] 

    scenID  <- mpTable$scenario[k]
    amID    <- mpTable$AM[k]


    totLossArray_SAisp[scenID,amID,1:dim(simLoss)[1],,] <- simLoss[1:dim(simLoss)[1],,]

    if(clearBadReps)
    {
      goodReps_isp <- goodRepsList[[k]]$goodReps_isp[1:nReps_k[k],,]
      totLossArray_SAisp[scenID,amID,,,][!goodReps_isp] <- NA

      # remove loss when nReps is below min sample size
      for(s in 1:nS)
        for( p in 1:nP )
          if( goodRepsList[[k]]$nReps_sp[s,p] < minSampSize )
            totLossArray_SAisp[scenID,amID,,s,p] <- NA            
    }
  }

  # Find median total loss
  medTotLoss_SAsp <- apply( X = totLossArray_SAisp,
                            FUN = median,
                            MARGIN = c(1,2,4,5),
                            na.rm = TRUE )

  # If there were NAs from non-convergent
  # runs then they will make a 
  medTotLoss_SAsp[medTotLoss_SAsp == 0] <- Inf
  medTotLoss_SAsp[is.na(medTotLoss_SAsp)] <- Inf

  # Now construct the table. I think we need 9 rows (spec/stock),
  # then 20 columns (scenario * AM) plus a couple for
  # labels
  tabColNames <- c( "Species",
                    "Stock",
                    "Scenario",
                    "AM",
                    "Rank" )
  rankTable <- matrix(NA, nrow = nS*nP*nScen*nAM, ncol = length(tabColNames) )
  colnames(rankTable) <- tabColNames

  rankTable <- as.data.frame(rankTable)

  rankStatTable <- matrix( NA, nrow = nScen, ncol = nAM )
  colnames(rankStatTable) <- AMs
  rownames(rankStatTable) <- scenarios

  rowIdx <- 1
  for( scenIdx in 1:length(scenarios))
  {
    
    scenID  <- scenarios[scenIdx]
    
    for( s in 1:nS )
      for( p in 1:nP )
      {
        medLossThisScen <- medTotLoss_SAsp[scenID,AMs,s,p]

        rankOrder <- order(medLossThisScen)
        rankVec <- numeric(length = nAM)
        rankVec[rankOrder] <- 1:5
        names(rankVec) <- AMs

        for( amIdx in 1:nAM )
        {
          amID    <- AMs[amIdx]
          # Now loop through scenarios, rank AMs
          # and put into table
          rankTable[rowIdx,"Species"]   <- speciesNames[s]
          rankTable[rowIdx,"Stock"]     <- stockNames[p]
          rankTable[rowIdx,"Scenario"]  <- scenID
          rankTable[rowIdx,"AM"]        <- amID
          rankTable[rowIdx,"Rank"]      <- rankVec[amID]

          rowIdx <- rowIdx + 1
        }

      }
  }

  summaryRankTable <- rankTable %>%
                      group_by(Scenario,AM) %>%
                      summarise(  minRank = min(Rank),
                                  maxRank = max(Rank),
                                  meanRank = sprintf("%.2f",round(mean(Rank),2)) )

  tabFileName <- paste(lossVar,lossType, sep = "_")
  rankTableFile <- here::here("Outputs",groupFolder,paste(tabFileName,"_rankTableFull.csv",sep = ""))
  summaryRankTableFile <- here::here("Outputs",groupFolder,paste(tabFileName,"_rankTableSumm.csv",sep = ""))

  write.csv(summaryRankTable, file = summaryRankTableFile)
  write.csv(rankTable, file = rankTableFile)

  outList <- list(  summRank = summaryRankTable,
                    fullRank = rankTable )

  outList
} # END calcLossRank()

# calcPairedLossRank()
# Calculates each AMs rank with respect
# to abs/rel biomass/catch loss WITHIN A REPLICATE. Outputs
# a table with scenarios as rows and
# AMs as columns
calcPairedLossRank <- function( groupFolder = "DERTACS_reruns_Jan5",
                                lossVar = "C_ispt",
                                lossType = "abs",
                                prefix = "parBat",
                                period = c(73:82),
                                clearBadReps = TRUE,
                                minSampSize = 100 )
{
  # First get info files so we can load the right 
  # loss objects
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs( here::here("Outputs",groupFolder),
                              recursive = FALSE, full.names = FALSE)

  simFolderList <- simFolderList[grepl(prefix, simFolderList)]

  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
              filter(grepl(prefix,simLabel))

  if( lossType == "rel" )
  {
    lossArrayName <- "totRelLoss_isp"
  }

  if( lossType == "abs" )
  {
    lossArrayName <- "totAbsLoss_isp"
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

    unlist(outList)
  }

  MPfactors <- lapply( X = info.df$mp, FUN = splitMP )
  MPfactors <- do.call(rbind, MPfactors) %>% as.data.frame()

  # Read in the performance tables
  mpTable <- cbind( info.df, MPfactors ) %>%
              mutate( perfPath = here::here("Outputs",groupFolder,simLabel,"simPerfStats.csv"),
                      path = here::here("Outputs",groupFolder,simLabel) )

  blobFiles <- file.path(mpTable$path,paste(mpTable$simLabel,".RData",sep = ""))
  lossFiles <- file.path(mpTable$path,"loss.RData")

  lossList <- lapply(X = mpTable$simLabel, FUN = .loadLoss, folder=  groupFolder)

  names(lossList) <- mpTable$simLabel

  # Calculate total loss for variable/period
  totLossList  <- lapply( X = lossList,
                          FUN = calcTotalLossPeriod,
                          var = lossVar, period = period )
  names(totLossList) <- names(lossList)

  # Now pull dimensions from the blob
  nT  <- lossList[[1]]$nT
  nS  <- lossList[[1]]$nS
  nP  <- lossList[[1]]$nP
  pT  <- dim(lossList[[1]]$retroSB_itspt)[2]
  tMP <- lossList[[1]]$tMP

  speciesNames <- lossList[[1]]$speciesNames[1:3]
  stockNames   <- lossList[[1]]$stockNames[1:3]


  # Get number of MPs
  nSims <- length(lossList)

  AMs       <- unique(mpTable$AM)
  Fsources  <- unique(mpTable$Fsrce)
  eqbm      <- unique(mpTable$eqbm)
  MPs       <- unique(info.df$mp)
  scenarios <- unique(info.df$scenario)

  nAM       <- length(AMs)
  nSrce     <- length(Fsources)
  nEqbm     <- length(eqbm)
  nMP       <- length(MPs)
  nScen     <- length(scenarios)

  nReps_k <- numeric(length = nSims)
  goodRepsList <- vector(mode = "list", length = nSims)

  for( k in 1:nSims )
  {
    simID <- info.df$simLabel[k]
    simLoss <- totLossList[[simID]][[lossArrayName]] 
    nReps_k[k] <- dim(simLoss)[1]
    goodRepsList[[k]]$goodReps_isp <- lossList[[k]]$goodReps_isp
    goodRepsList[[k]]$nReps_sp     <- apply(X = lossList[[k]]$goodReps_isp, FUN = sum, MARGIN = c(2,3) )
    names(goodRepsList)[k] <- simID
  }

  nReps         <- max(nReps_k)
  speciesNames  <- lossList[[1]]$speciesNames
  stockNames    <- lossList[[1]]$stockNames

  # Make an array to hold loss function values
  totLossArray_SAisp <- array(NA,  dim = c(nScen, nAM, nReps, nS, nP ),
                                  dimnames = list(  scenario = scenarios,
                                                    AM = AMs,
                                                    rep = 1:nReps,
                                                    species = speciesNames[1:nS],
                                                    stock = stockNames[1:nP] ) )

  goodReps_SAisp <- array(NA,  dim = c(nScen, nAM, nReps, nS, nP ),
                                  dimnames = list(  scenario = scenarios,
                                                    AM = AMs,
                                                    rep = 1:nReps,
                                                    species = speciesNames[1:nS],
                                                    stock = stockNames[1:nP] ) )

  for( k in 1:nSims )
  {
    simID <- mpTable$simLabel[k]
    simLoss <- totLossList[[simID]][[lossArrayName]] 

    scenID  <- mpTable$scenario[k]
    amID    <- mpTable$AM[k]


    totLossArray_SAisp[scenID,amID,1:dim(simLoss)[1],,] <- simLoss[1:dim(simLoss)[1],1:nS,1:nP]

    thisGoodReps_isp <- goodRepsList[[simID]]$goodReps_isp
    goodReps_SAisp[scenID, amID,,,] <- thisGoodReps_isp
  }

  rankArray_ASisp <- apply( X = totLossArray_SAisp, FUN = calcRankVec,
                            MARGIN = c(1,3,4,5) )

  rankArray_SAisp <- aperm(rankArray_ASisp,c(2,1,3,4,5))
  dimnames( rankArray_SAisp ) <- dimnames( totLossArray_SAisp )

  attr(rankArray_SAisp,"conv") <- goodReps_SAisp

  rankArray_SAisp
} # END calcPairedLossRank()

# a function that takes in a vector, and ranks
# the values in that vector
calcRankVec <- function( x )
{
  rankOrder <- order(x)
  nRanks <- length(x)
  rankVec <- numeric(length = nRanks)
  rankVec[rankOrder] <- 1:nRanks

  rankVec
} # END calcRankVec()


# calcLoss()
# Loads blob and calculates yearly catch and biomass 
# loss values compared to a nominated baseline sim.
# By default, calculates biomass and catch loss
# on relative and absolute scale
calcLoss <- function( sim         = 2,
                      baseline    = "sim_omni_totCat_200reps",
                      groupFolder = "DERTACS_reruns_Jan5",
                      lossVars    = c("C_ispt","SB_ispt"),
                      output      = TRUE )
{
  # First, load the sim that we want
  # to calculate loss for
  simFolder <- .loadSim( sim = sim, folder = groupFolder )

  # Save the blob
  lossSim <- blob

  # figure out which reps we want
  goodReps_isp    <- blob$goodReps_isp
  maxRep_sp       <- apply(X = goodReps_isp, FUN = maxWhich, MARGIN = c(2,3))
  totReps         <- dim(goodReps_isp)[1]



  # Now load the baseline (omniscient manager)
  .loadSim( sim = baseline, folder = groupFolder )

  baseSim <- blob

  blob <- NULL
  gc()

  # Model dimensions
  tMP <- lossSim$om$tMP
  nT  <- lossSim$om$nT
  nF  <- lossSim$om$nF
  nS  <- lossSim$om$nS
  nP  <- lossSim$om$nP
  pT  <- lossSim$om$pT

  speciesNames  <- c(lossSim$om$speciesNames,"dataPooled")
  stockNames    <- c(lossSim$om$stockNames,"coastWide")
  fleetNames    <- lossSim$om$fleetNames

  # We want to calculate loss of given OM states
  # Make a list to hold baseline
  baseLineStates        <- vector(  mode = "list", 
                                    length = length(lossVars))
  names(baseLineStates) <- lossVars
  # copy for simulation states and loss
  simStates       <- baseLineStates
  lossRaw         <- baseLineStates
  lossRel         <- baseLineStates


  # Pull control list settings for
  # coastwide and data pooled

  lossCtl         <- lossSim$ctlList


  # Loop over 
  for( varIdx in 1:length(lossVars) )
  {
    # Clean arrays for base and sim states
    baseState_ispt <- array(NA, dim = c(totReps,nS+1,nP+1,nT),
                                dimnames = list(  rep = 1:totReps,
                                                  species = speciesNames,
                                                  stock   = stockNames,
                                                  tdx     = 1:nT ) )
    simState_ispt  <- baseState_ispt

    # Get variable
    var <- lossVars[varIdx]

    # Copy base and sim for individual stocks/species
    simReps <- (1:totReps)
    baseState_ispt[,1:nS,1:nP,] <- baseSim$om[[var]][simReps,,,1:nT]
    simState_ispt[,1:nS,1:nP,]  <- lossSim$om[[var]][simReps,,,1:nT]

    # Now do the aggregates
    # totally aggregated
    baseState_ispt[,nS+1,nP+1,] <- apply(   X = baseSim$om[[var]][simReps,,,1:nT],
                                            FUN = sum,
                                            MARGIN = c(1,4), na.rm = TRUE )

    simState_ispt[,nS+1,nP+1,]  <- apply(   X = lossSim$om[[var]][simReps,,,1:nT],
                                            FUN = sum,
                                            MARGIN = c(1,4), na.rm = TRUE )

    # data-pooled
    baseState_ispt[,nS+1,1:nP,]     <- apply(   X = baseSim$om[[var]][simReps,,,1:nT],
                                            FUN = sum,
                                            MARGIN = c(1,3,4), na.rm = TRUE )

    simState_ispt[,nS+1,1:nP,]      <- apply(   X = lossSim$om[[var]][simReps,,,1:nT],
                                            FUN = sum,
                                            MARGIN = c(1,3,4), na.rm = TRUE )
    
    # coastwide
    baseState_ispt[,1:nS,nP+1,]     <- apply(   X = baseSim$om[[var]][simReps,,,1:nT],
                                            FUN = sum,
                                            MARGIN = c(1,2,4), na.rm = TRUE )

    simState_ispt[,1:nS,nP+1,]      <- apply(   X = lossSim$om[[var]][simReps,,,1:nT],
                                            FUN = sum,
                                            MARGIN = c(1,2,4), na.rm = TRUE )

    baseLineStates[[var]] <- baseState_ispt
    simStates[[var]]      <- simState_ispt

    lossRaw[[var]]        <- baseLineStates[[var]] - simStates[[var]]
    lossRel[[var]]        <- lossRaw[[var]] / baseLineStates[[var]]
   
  }

  outList <- list(  simID         = lossSim$folder,
                    goodReps_isp  = goodReps_isp,
                    speciesNames  = speciesNames,
                    stockNames    = stockNames,
                    fYear         = lossSim$ctlList$opMod$fYear,
                    simStates     = simStates,
                    baseStates    = baseLineStates,
                    lossRaw       = lossRaw,
                    lossRel       = lossRel,
                    retroSB_itspt = lossSim$mp$assess$retroSB_itspt,
                    tMP           = tMP,
                    nT            = nT, 
                    nF            = nF, 
                    nS            = nS, 
                    nP            = nP, 
                    pT            = pT,
                    obsCVmult     = lossSim$ctlList$opMod$projObsErrMult,
                    simFolder     = simFolder  )


  # Save loss to sim folder
  save( outList, file = file.path(simFolder,"loss.RData") )

  if(output)
    return(outList)
} # END calcLoss()



makeLossTable <- function(  sim = 2, 
                            groupFolder = "DLSurveys7_.5tau",
                            save = TRUE )
{
  # Load loss reports
  objLoss <- .loadLoss(sim = sim, groupFolder)

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
  obsCVmult       <- objLoss$obsCVmult




  # Here's the loss table
  lossTableColNames <- c( "sim",
                          "scenario",
                          "mp",
                          "species",
                          "stock",
                          "obsCVmult",
                          "relLossCat_short",
                          "relLossCat_med",
                          "relLossCat_long",
                          "absLossCat_short",
                          "absLossCat_med",
                          "absLossCat_long",
                          "relLossBio_short",
                          "relLossBio_med",
                          "relLossBio_long",
                          "absLossBio_short",
                          "absLossBio_med",
                          "absLossBio_long" )

  lossTable <- matrix(NA, nrow = (nS * 1) * (nP + 1), ncol = length(lossTableColNames) )
  colnames(lossTable) <- lossTableColNames
  lossTable <- as.data.frame(lossTable)

  lossTable$sim       <- objLoss$sim
  lossTable$scenario  <- objLoss$scenario
  lossTable$mp        <- objLoss$mp
  lossTable$obsCVmult <- obsCVmult

  
  lossCat_shortTerm  <- calcTotalLossPeriod(objLoss,"C_ispt", period = tMP:(tMP+10))
  lossCat_medTerm    <- calcTotalLossPeriod(objLoss,"C_ispt", period = tMP:(tMP+20))
  lossCat_longTerm   <- calcTotalLossPeriod(objLoss,"C_ispt", period = tMP:nT)
  lossBio_shortTerm  <- calcTotalLossPeriod(objLoss,"SB_ispt", period = tMP:(tMP+10))
  lossBio_medTerm    <- calcTotalLossPeriod(objLoss,"SB_ispt", period = tMP:(tMP+20))
  lossBio_longTerm   <- calcTotalLossPeriod(objLoss,"SB_ispt", period = tMP:nT)

  for( s in 1:(nS+1) )
    for( p in 1:(nP+1) )
    {
      tabRow <- (s - 1)*( nP + 1 ) + p

      lossTable[tabRow,"species"] <- speciesNames[s]
      lossTable[tabRow,"stock"]   <- stockNames[p]
      lossTable[tabRow,"relLossCat_short"]    <- median(lossCat_shortTerm$totRelLoss_isp[,s,p])
      lossTable[tabRow,"relLossCat_med"]      <- median(lossCat_medTerm$totRelLoss_isp[,s,p])
      lossTable[tabRow,"relLossCat_long"]     <- median(lossCat_longTerm$totRelLoss_isp[,s,p])
      lossTable[tabRow,"absLossCat_short"]    <- median(lossCat_shortTerm$totAbsLoss_isp[,s,p])
      lossTable[tabRow,"absLossCat_med"]      <- median(lossCat_medTerm$totAbsLoss_isp[,s,p])
      lossTable[tabRow,"absLossCat_long"]     <- median(lossCat_longTerm$totAbsLoss_isp[,s,p])
      lossTable[tabRow,"relLossBio_short"]    <- median(lossBio_shortTerm$totRelLoss_isp[,s,p])
      lossTable[tabRow,"relLossBio_med"]      <- median(lossBio_medTerm$totRelLoss_isp[,s,p])
      lossTable[tabRow,"relLossBio_long"]     <- median(lossBio_longTerm$totRelLoss_isp[,s,p])
      lossTable[tabRow,"absLossBio_short"]    <- median(lossBio_shortTerm$totAbsLoss_isp[,s,p])
      lossTable[tabRow,"absLossBio_med"]      <- median(lossBio_medTerm$totAbsLoss_isp[,s,p])
      lossTable[tabRow,"absLossBio_long"]     <- median(lossBio_longTerm$totAbsLoss_isp[,s,p])
    
    }

  write.csv( lossTable, file = file.path(simFolder,"lossTable.csv") )

  lossTable
} # END makeLossTable()

# calcTotalLossPeriod()
# Takes a loss object, and variable name,
# and calculates total abs and relative loss
# for a given time period
calcTotalLossPeriod <- function(  obj, 
                                  var = "C_ispt",
                                  period = tMP:nT )
{
  # Pull abs and rel loss
  relLoss_ispt  <- obj$lossRel[[var]][,,,period]
  rawLoss_ispt  <- obj$lossRaw[[var]][,,,period]

  totRelLoss_isp    <- apply( X = abs(relLoss_ispt), FUN = sum, MARGIN = c(1,2,3) )
  totAbsLoss_isp    <- apply( X = abs(rawLoss_ispt), FUN = sum, MARGIN = c(1,2,3) )

  out <- list(  totRelLoss_isp = totRelLoss_isp,
                totAbsLoss_isp = totAbsLoss_isp )

  return(out)
} # END calcTotalLossPeriod



# calcBatchTable()
# Function to calculate a table of stats for 
# a whole batch (usually a given subfolder of
# ./Outputs/)
calcBatchTable <- function( batchFolder = "" )
{
  # set output folder
  simFolder <- here::here("Outputs",batchFolder)
  statsFolder <- file.path(simFolder, "statistics")

  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path=simFolder,full.names = FALSE,
                        recursive=FALSE)
  # Restrict to sim_ folders, count sims
  simList <- dirList[grep(pattern="sim",x=dirList)]
  nSims <- length(simList)

  statsTables <- lapply(  X = 1:nSims, FUN = makeStatTable,
                          folder = batchFolder )

  batchStatTable <- do.call("rbind", statsTables)

  # Now save to batchFolder
  if(!dir.exists(file.path(simFolder,"statistics")))
    dir.create(file.path(simFolder,"statistics"))

  outFile <- file.path(simFolder,"statistics","fullStatTable.csv")

  write.csv( batchStatTable, file = outFile)

  message(" (calcBatchTable) Stats table calculated and saved to ", outFile, "\n")

}


# makeStatTable()
# Wrapper for 
makeStatTable <- function(  sim = 1, folder = "",
                            save = FALSE )
{
  # First, load blob
  source("tools.R")

  .loadSim(sim, folder = folder)

  statTable <- .simPerfStats( obj = blob )

  simFolder <- statTable$simLabel[1]

  if( save )
  {
    outFile <- here::here("Outputs",folder,simFolder,"simPerfStats.csv")

    write.csv(statTable, file = outFile)

    message(paste("Stat table saved to ", outFile,"\n", sep ="" ) )
  }

  statTable
} # END makeStatTable()

# .simPerfStats()
# Produces a statistics table for a simulation from
# the output produced by a runMS3() call
# inputs:   sim=int indicating which simulation to compute stats for
# outputs:  statTable=data.frame showing conservation and catch performance
# usage:    in lapply to produce stats for a group of simulations
.simPerfStats <- function( obj, period = NULL  )
{
  # Pull sublists
  om        <- obj$om
  opMod     <- obj$ctlList$opMod
  mp        <- obj$mp
  ctlList   <- obj$ctlList
  rp        <- obj$rp[[1]]
  
  # Control info
  nS      <- om$nS
  nP      <- om$nP
  nT      <- om$nT
  tMP     <- om$tMP
  pT      <- opMod$pT
  nReps   <- ctlList$ctl$nReps

  if( is.null(period) )
    period <- tMP:nT

  # Species and stock labels
  speciesNames  <- opMod$species
  stockNames    <- opMod$stocks
  fYear         <- opMod$fYear

  # get the replicate numbers for succesful fits (PD Hessians) in
  # MPs
  # Get max number of reps
  maxReps <- obj$nSims  
  allConvReps_isp <- obj$goodReps_isp
  if( is.null(obj$goodReps_isp))
  {
    allConvReps_isp <- array(FALSE, dim = c(maxReps,nS,nP) )
    for( s in 1:nS )
      for(p in 1:nP )
        allConvReps_isp[,s,p] <- obj$goodReps
  }
  
  pdHess_itsp     <- mp$assess$pdHess_itsp


  # Calculate probability of a good replicate (all PD hessians)
  nGoodReps_sp <- apply(X = allConvReps_isp, FUN = sum, MARGIN =c(2,3))
  pGoodReps_sp <- signif(nGoodReps_sp/maxReps,2)
  maxGoodReps <- max(nGoodReps_sp)

  # And the median over replicates of the probability of a 
  # PD Hessian over time
  probPDHess_isp <- apply( X = pdHess_itsp, FUN = mean, MARGIN = c(1,3,4), na.rm = T )
  medProbPDH_sp  <- apply( X = probPDHess_isp, FUN = median, MARGIN = c(2,3), na.rm = T )

  if(is.null(obj$simLabel))
  {
    simLabel <- stringr::str_split(obj$path,"/")[[1]]
    simLabel <- simLabel[length(simLabel)]
  } else simLabel <- obj$simLabel

  # First, create a data.frame of NAs with a row for each of MRE,MARE
  colLabels <- c( "simLabel",
                  "scenario","mp",
                  "species","stock",
                  "projObsErrMult",
                  "nGoodReps","pGoodReps", "medProbPDH_sp",
                  "pBtGt.4Bmsy", "PBtGt.8Bmsy",
                  "pBtGtBmsy",
                  "pCtGtMSY", "pFtGtFmsy",
                  "avgCatch","AAV", "avgTACu",
                  "medDiscRent", "medEff",
                  "VAR5pc",
                  "pHistLowCatch" )

  statTable <- matrix( NA,  ncol = length(colLabels),
                            nrow = nS * nP )
  colnames(statTable) <- colLabels

  statTable <- as.data.frame(statTable)

  # Fill in global values
  statTable[,"simLabel"]        <- simLabel
  statTable[,"scenario"]        <- ctlList$ctl$scenarioName
  statTable[,"mp"]              <- ctlList$ctl$mpName
  statTable[,"projObsErrMult"]  <- opMod$projObsErrMult

  # Need to start layering in performance metrics


  # Pull reference points
  B0_sp     <- opMod$histRpt$B0_sp
  Bmsy_sp   <- rp$FmsyRefPts$BeqFmsy_sp
  Fmsy_sp   <- rp$FmsyRefPts$Fmsy_sp
  Umsy_sp   <- rp$FmsyRefPts$Umsy_sp
  MSY_sp    <- rp$FmsyRefPts$YeqFmsy_sp

  if( any(nGoodReps_sp) > 0 )
  { 
    # Calculate depletion wrt Bmsy
    SB_ispt   <- om$SB_ispt[1:maxGoodReps,,,,drop = FALSE]
    C_ispt    <- om$C_ispt[1:maxGoodReps,,,,drop = FALSE]
    TAC_ispt  <- mp$hcr$TAC_ispt[1:maxGoodReps,,,,drop = FALSE]
    TACu_ispt <- C_ispt / TAC_ispt
    F_ispt    <- om$F_ispft[1:maxGoodReps,,,2,]
    
    # Calculate probability that Bt above LRP
    pBtGt.4Bmsy_sp <- .calcStatsProportion( TS_ispt = SB_ispt,
                                            ref_sp = Bmsy_sp,
                                            tdx = period,
                                            prop = .4,
                                            nS = nS,
                                            nP = nP )
    # In a healthy state (Bt > .8Bmsy)
    pBtGt.8Bmsy_sp <- .calcStatsProportion( TS_ispt = SB_ispt,
                                            ref_sp = Bmsy_sp,
                                            tdx = period,
                                            prop = .8,
                                            nS = nS,
                                            nP = nP )

    # Above Bmsy
    pBtGtBmsy_sp <- .calcStatsProportion( TS_ispt = SB_ispt,
                                          ref_sp = Bmsy_sp,
                                          tdx = period,
                                          prop = 1.0,
                                          nS = nS,
                                          nP = nP )  

    # Overfishing is occuring (Ft > Fmsy)
    pFtGtFmsy_sp <- .calcStatsProportion( TS_ispt = F_ispt,
                                          ref_sp = Fmsy_sp,
                                          tdx = period,
                                          prop = 1,
                                          nS = nS,
                                          nP = nP )

    # Overfishing is occuring (Ct > MSY)
    pCtGtMSY_sp <- .calcStatsProportion(  TS_ispt = C_ispt,
                                          ref_sp = MSY_sp,
                                          tdx = period,
                                          prop = 1,
                                          nS = nS,
                                          nP = nP )

    # Find minimum catch
    minHistCatch_sp <- apply( X = C_ispt[,,,1:(tMP-1)], FUN = min,
                              MARGIN = c(2,3) )
    pCtLtHistMinC_sp <- .calcStatsProportion( TS_ispt = C_ispt,
                                              ref_sp = minHistCatch_sp,
                                              tdx = period,
                                              prop = 1,
                                              nS = nS,
                                              nP = nP )


    # Average catch and TAC utilisation
    Cbar_sp     <- apply( X = C_ispt[,,,tMP:nT], FUN = mean, MARGIN = c(2,3), na.rm = T)
    TACubar_sp  <- apply( X = TACu_ispt[,,,tMP:nT], FUN = mean, MARGIN = c(2,3), na.rm = T)

    # Catch variability
    AAV_sp      <- .calcStatsAAV( C_ispt = C_ispt,
                                  tdx = period,
                                  qProbs = c(0.5),
                                  margin = c(2,3) )

    # Calculate profit
    rev_ispft      <- obj$om$Rev_ispft[1:maxGoodReps,,,,]
    rev_ipft       <- apply( X = rev_ispft, FUN = sum, MARGIN = c(1,3,4,5))
    effCost_ipft   <- obj$om$effCost_ipft[1:maxGoodReps,,,]
    crewShare      <- opMod$crewShare
    discRate       <- opMod$discountRate

    profit_ipft    <- array(0, dim = dim(rev_ipft))
    profit_ipft    <- rev_ipft - effCost_ipft
    
    discProfit_ipft  <- profit_ipft
    for( t in tMP:nT )
      discProfit_ipft[,,,t] <- profit_ipft[,,,t] * (1 + discRate)^(-(t - tMP))


    discTotalProfit_ip      <- apply( X = discProfit_ipft[,,2,tMP:nT], FUN = sum, MARGIN = c(1,2) )
    VAR5pc_p                <- apply( X = discTotalProfit_ip, FUN = quantile, probs = 0.05, MARGIN = c(2), na.rm = T )
    medDiscTotalProfit_p    <- apply( X = discTotalProfit_ip, FUN = median, MARGIN = c(2), na.rm = T )

    E_ipft    <- obj$om$E_ipft[1:maxGoodReps,,,]
    medEff_ip <- apply( X = E_ipft[,,2,], FUN = median, MARGIN = c(1,2),na.rm = T)
    medEff_p  <- apply( X = medEff_ip, FUN = median, MARGIN = c(2),na.rm = T)

    for( s in 1:nS )
      for(p in 1:nP )
      {
        

        rowIdx <- (s-1) * nP + p

        statTable[rowIdx, "species"]        <- speciesNames[s]
        statTable[rowIdx, "stock"]          <- stockNames[p]
        statTable[rowIdx, "medProbPDH_sp"]  <- medProbPDH_sp[s,p]

        statTable[rowIdx, "nGoodReps"]      <- nGoodReps_sp[s,p]
        statTable[rowIdx, "pGoodReps"]      <- pGoodReps_sp[s,p]

        # Conservation performance
        statTable[rowIdx,"pBtGt.4Bmsy"]     <- round(pBtGt.4Bmsy_sp[s,p],2)
        statTable[rowIdx,"PBtGt.8Bmsy"]     <- round(pBtGt.8Bmsy_sp[s,p],2)
        statTable[rowIdx,"pBtGtBmsy"]       <- round(pBtGtBmsy_sp[s,p],2)
        statTable[rowIdx,"pCtGtMSY"]        <- round(pCtGtMSY_sp[s,p],2)
        statTable[rowIdx,"pFtGtFmsy"]       <- round(pFtGtFmsy_sp[s,p],2)

        # Catch statistics
        statTable[rowIdx,"avgCatch"]        <- round(Cbar_sp[s,p],2)
        statTable[rowIdx,"AAV"]             <- round(AAV_sp[s,p],2)
        statTable[rowIdx,"avgTACu"]         <- round(TACubar_sp[s,p],2)
        statTable[rowIdx,"medDiscRent"]     <- round(medDiscTotalProfit_p[p],2)
        statTable[rowIdx,"VAR5pc"]          <- round(VAR5pc_p[p],2)
        statTable[rowIdx,"medEff"]          <- round(medEff_p[p],2)
        statTable[rowIdx,"pHistLowCatch"]   <- round(1 - pCtLtHistMinC_sp[s,p],2)

      }
  }

  statTable
} # END .simPerfStats()

# .calcStatsProportion()
# Calculates the probability of a time series
# being above a certain proportion of a reference
# level. Used for depletion and overfishing statistics.
.calcStatsProportion <- function( TS_ispt = SB_ispt,
                                  ref_sp = Bmsy_sp,
                                  tdx = tMP:nT,
                                  prop = .4,
                                  nS = 3,
                                  nP = 3 )
{

  # Reduce to time period
  TS_ispt <- TS_ispt[,,,tdx]

  Quotient_ispt <- array(NA, dim = dim(TS_ispt) )

  # First, take quotient by reference
  for(s in 1:nS)
    for(p in 1:nP )
      Quotient_ispt[,s,p,] <- TS_ispt[,s,p,] / ref_sp[s,p] 

  # Set an indicator array
  Ind_ispt <- Quotient_ispt > prop

  probGtDep_sp <- apply( X = Ind_ispt, FUN = mean, MARGIN = c(2,3), na.rm = T )
      
  return( probGtDep_sp )
} # END .calcStatsProportion

# .calcStatsAAV()
# Calculate average annual variation
# in catch, measure of the variability
# in landings for each stock/species.
.calcStatsAAV <- function(  C_ispt = C_ispt,
                            tdx = tMP:nT,
                            qProbs = c(0.025,0.5,0.975),
                            margin = c(2,3) )
{
  C_ispt <- C_ispt[,,,tdx]
  # Append margin
  marg <- c(1,margin)
  # Add over margins, in case we are looking
  # complex/species/stock aggregates
  sumC_ispt <- apply( X = C_ispt, FUN = sum, MARGIN = c(marg,4), na.rm = T )


  # Calculate diff
  diffC_ispt        <- aperm(apply( X = sumC_ispt, FUN = diff, MARGIN = marg, na.rm = T ),c(margin,4,1))
  absDiffC_ispt     <- abs(diffC_ispt)
  sumAbsDiffC_isp   <- apply( X = absDiffC_ispt, FUN = sum, MARGIN = marg, na.rm = T )
  sumCatch_isp      <- apply( X = C_ispt, FUN = sum, MARGIN = marg, na.rm = T )

  AAV_isp           <- sumAbsDiffC_isp / sumCatch_isp
  AAV_isp[!is.finite(AAV_isp)] <- 0

  AAV_qsp   <- apply(X = AAV_isp, FUN = quantile, probs = qProbs, MARGIN = marg[-1] )

  return(AAV_qsp)
} # END .calcStatsAAV


# .getCplxStats()
# A pared down perf stats function for complex level
# quantities (i.e. total catch across all species/stocks etc.)
# Useful for comparing omniscient manager obj fun
# weightings
.getCplxStats <- function( obj )
{
  # Model dims
  tMP <- obj$om$tMP
  nT  <- obj$om$nT
  nS  <- obj$om$nS
  nP  <- obj$om$nP

  # Names of stuff
  scenarioName  <- obj$ctlList$ctl$scenarioName
  mpName        <- obj$ctlList$ctl$mpName

  # Pull simulation label
  if(is.null(obj$simLabel))
  {
    simLabel <- stringr::str_split(obj$path,"/")[[1]]
    simLabel <- simLabel[length(simLabel)]
  } else simLabel <- obj$simLabel

  allConvReps_isp <- obj$goodReps_isp
  allConvReps <- apply( X = allConvReps_isp, FUN = prod, MARGIN = 1 )

  projYrs <- tMP:nT

  C_ispt <- obj$om$C_ispt[allConvReps,,,projYrs]
  B_ispt <- obj$om$SB_ispt[allConvReps,,,projYrs]

  # Pull Bmsy
  Bmsy_sp <- obj$rp[[1]]$FmsyRefPts$BeqFmsy_sp

  # Turn into complex catch
  C_it <- apply( X = C_ispt, FUN = sum, MARGIN = c(1,4))

  # Calculate AAV at the complex level
  diffC_ti      <- apply( X = C_it, FUN = diff, MARGIN = 1)
  diffC_it      <- abs(t(diffC_ti))
  sumAbsDiff_i  <- apply( X = diffC_it, FUN = sum, MARGIN = 1)
  sumCatch_i    <- apply( X = C_it, FUN = sum, MARGIN = 1)
  AAVC_i        <- sumAbsDiff_i / sumCatch_i
  AAVCplx       <- round(median( AAVC_i),2)
  
  # We want median (over reps) average and total complex catch
  Cbar    <- round(median(apply( X = C_it, FUN = mean, MARGIN = 1 ) ),2)
  totC    <- round(median(apply( X = C_it, FUN = sum,  MARGIN = 1 ) ),2)

  # Now loop over species and stocks, calculate
  # depletion probs (make the following values inputs either
  # as function args or in the ctl file)
  LRP <- 0.4
  USR <- 2

  depBmsy_ispt <- B_ispt
  for( s in 1:nS )
    for( p in 1:nP )
    {
      depBmsy_ispt[,s,p,] <- B_ispt[,s,p,]/Bmsy_sp[s,p]
    }

  depGtLRP_ispt <- depBmsy_ispt > LRP
  depLtUSR_ispt <- depBmsy_ispt < USR

  probDepGtLRP_isp <- apply(X = depGtLRP_ispt, FUN = mean, MARGIN = c(1,2,3))
  probDepLtUSR_isp <- apply(X = depLtUSR_ispt, FUN = mean, MARGIN = c(1,2,3))

  # Get distribution of depletion within desired region 
  # across species, stocks, and replicates
  distProbDepLRP_q <- quantile( probDepGtLRP_isp, probs = c(0.025, 0.5, 0.975) )
  distProbDepUSR_q <- quantile( probDepLtUSR_isp, probs = c(0.025, 0.5, 0.975) )

  # Make above vectors into a distribution as "med (CI limits)"
  distLRP_chr <- paste(distProbDepLRP_q[2], " (", distProbDepLRP_q[1],", ", distProbDepLRP_q[3], ")", sep = "")
  distUSR_chr <- paste(distProbDepUSR_q[2], " (", distProbDepUSR_q[1],", ", distProbDepUSR_q[3], ")", sep = "")

  # Make into a list
  out.df <- data.frame( simLabel    = simLabel,
                        scenario    = scenarioName,
                        mp          = mpName,
                        Cbar        = round(Cbar,2),
                        totC        = round(totC,2),
                        AAV         = round(AAVCplx,2),
                        PGtLRP      = distLRP_chr,
                        PLtUSR      = distUSR_chr )

  out.df
} # END .getCplxStats()


# .getOmniInfo()
# Pulls omniscient manager information from the blob and control
# list.
.getOmniInfo <- function(obj)
{
  # Model dims
  tMP <- obj$om$tMP
  nT  <- obj$om$nT
  nS  <- obj$om$nS
  nP  <- obj$om$nP

  # Names of stuff
  scenarioName  <- obj$ctlList$ctl$scenarioName
  mpName        <- obj$ctlList$ctl$mpName
  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames

  # Pull simulation label
  if(is.null(obj$simLabel))
  {
    simLabel <- stringr::str_split(obj$path,"/")[[1]]
    simLabel <- simLabel[length(simLabel)]
  } else simLabel <- obj$simLabel



  allConvReps_isp <- obj$goodReps_isp
  allConvReps <- apply( X = allConvReps_isp, FUN = prod, MARGIN = 1)

  projYrs <- tMP:nT

  # Get omni sub-ctlList
  omniList    <- obj$ctlList$mp$omni
  omniObjFun  <- obj$omniObjFun

  # Pull obj fun targets
  loDepBmsy     <- omniList$loDepBmsy
  hiDepBmsy     <- omniList$hiDepBmsy
  maxF          <- omniList$maxF
  maxE          <- omniList$maxE
  maxAAV        <- omniList$maxAAV
  maxRelEffDiff <- omniList$maxRelEffDiff
  maxRelCatDiff <- omniList$maxRelCatDiff
  # component weights
  avgCatWt      <- omniList$avgCatWt
  totCatWt      <- omniList$totCatWt
  sumCatWt      <- omniList$sumCatWt
  hiDepBmsyWt   <- omniList$hiDepBmsyWt
  loDepBmsyWt   <- omniList$loDepBmsyWt
  probDepWt     <- omniList$probDepWt
  closedWt      <- omniList$closedWt
  AAVWt         <- omniList$AAVWt
  catDiffWt     <- omniList$catDiffWt
  effDiffWt     <- omniList$effDiffWt
  initCatDiffWt <- omniList$initCatDiffWt
  initEffDiffWt <- omniList$initEffDiffWt

  # Objective function penalty controls
  penType       <- omniList$penType
  linBeta       <- omniList$linBeta

  # Now pull objective function values

  # Complex level
  objFun_i            <- omniObjFun$objFun_i[allConvReps]
  totCbar_i           <- omniObjFun$totCbar_i[allConvReps]
  Csum_i              <- omniObjFun$Csum_i[allConvReps]
  barEffDiff_ip       <- omniObjFun$barEffDiff_ip[allConvReps,]


  # Species/stock
  objFun_isp          <- omniObjFun$objFun_isp[allConvReps,,]
  Cbar_isp            <- omniObjFun$Cbar_isp[allConvReps,,]
  barLoDep_isp        <- omniObjFun$barLoDep_isp[allConvReps,,]
  barHiDep_isp        <- omniObjFun$barHiDep_isp[allConvReps,,]
  barProbLoDep_isp    <- omniObjFun$barProbLoDep_isp[allConvReps,,]
  barProbHiDep_isp    <- omniObjFun$barProbHiDep_isp[allConvReps,,]
  barInitCatDiff_isp  <- omniObjFun$barInitCatDiff_isp[allConvReps,,]
  closedCount_isp     <- omniObjFun$closedCount_isp[allConvReps,,]
  barAAV_isp          <- omniObjFun$barAAV_isp[allConvReps,,]
  barCatDiff_isp      <- omniObjFun$barCatDiff_isp[allConvReps,,]

  # Sum effort diff over stocks
  barEffDiff_i        <- omniObjFun$barEffDiff_i[allConvReps]
  barInitEffDiff_i    <- omniObjFun$barInitEffDiff_i[allConvReps]

  # Now, we want to make a couple of tables
  # make a table of weights as well, this will be good for
  # tuning the weights later...
  wtTableColNames <- c( 'simLabel',
                        'scenario',
                        'mp',
                        'wt_CbarCplx',
                        'wt_CbarStock',
                        'wt_totC',
                        'wt_hiDep',
                        'wt_loDep',
                        'wt_probDep',
                        'wt_closure',
                        'wt_AAV',
                        'wt_catDiff',
                        'wt_effDiff',
                        'wt_initCatDiff',
                        'wt_initEffDiff')

  wtTable <- matrix(NA, nrow = 1, ncol = length(wtTableColNames) )
  colnames(wtTable) <- wtTableColNames

  wtTable[,'simLabel']        <- simLabel
  wtTable[,'scenario']        <- scenarioName
  wtTable[,'mp']              <- mpName
  wtTable[,'wt_CbarCplx']     <- totCatWt
  wtTable[,'wt_CbarStock']    <- avgCatWt
  wtTable[,'wt_totC']         <- sumCatWt
  wtTable[,'wt_hiDep']        <- hiDepBmsyWt
  wtTable[,'wt_loDep']        <- loDepBmsyWt
  wtTable[,'wt_probDep']      <- probDepWt
  wtTable[,'wt_closure']      <- closedWt
  wtTable[,'wt_AAV']          <- AAVWt
  wtTable[,'wt_catDiff']      <- catDiffWt
  wtTable[,'wt_effDiff']      <- effDiffWt
  wtTable[,'wt_initCatDiff']  <- initCatDiffWt
  wtTable[,'wt_initEffDiff']  <- initCatDiffWt


  # First, the complex level quantities
  cplxTableColNames <- c( 'simLabel',
                          'scenario',
                          'mp',
                          'objFun',
                          'wt_CbarCplx',
                          'of_CbarCplx',
                          'wt_totC',
                          'of_totC',
                          'wt_EffDiff',
                          'of_EffDiff',
                          'wt_initEffDiff',
                          'of_initEffDiff' )
  cplxOmniObjFunTable <- matrix(NA, nrow = 1, ncol = length(cplxTableColNames) )
  colnames(cplxOmniObjFunTable) <- cplxTableColNames
  cplxOmniObjFunTable <- as.data.frame(cplxOmniObjFunTable)

  cplxOmniObjFunTable[,'simLabel']    <- simLabel
  cplxOmniObjFunTable[,'scenario']    <- scenarioName
  cplxOmniObjFunTable[,'mp']          <- mpName
  cplxOmniObjFunTable[,'objFun']      <- round(mean(objFun_i),2)
  cplxOmniObjFunTable[,'of_CbarCplx'] <- round(mean(totCbar_i),2)
  cplxOmniObjFunTable[,'wt_CbarCplx'] <- totCatWt
  cplxOmniObjFunTable[,'of_totC']     <- round(mean(Csum_i),2)
  cplxOmniObjFunTable[,'wt_totC']     <- sumCatWt
  cplxOmniObjFunTable[,'of_EffDiff']  <- round(mean(barEffDiff_i),2)
  cplxOmniObjFunTable[,'wt_EffDiff']  <- effDiffWt
  cplxOmniObjFunTable[,'of_EffDiff']  <- round(mean(barInitEffDiff_i),2)
  cplxOmniObjFunTable[,'wt_EffDiff']  <- initEffDiffWt


  # Now a stock/species objective function table
  stockTableColNames <- c(  'simLabel',
                            'scenario',
                            'mp',
                            'species',
                            'stock',
                            'objFun',
                            'wt_CbarStock',
                            'of_CbarStock',
                            'wt_hiDep',
                            'of_hiDep',
                            'wt_loDep',
                            'of_loDep',
                            'wt_probDep',
                            'of_probHiDep',
                            'of_probLoDep',
                            'wt_closure',
                            'of_closure',
                            'wt_AAV',
                            'of_AAV',
                            'wt_catDiff',
                            'of_catDiff',
                            'wt_initCatDiff',
                            'of_initCatDiff' )

  stockObjFunTable <- matrix( NA, ncol = length(stockTableColNames), nrow = nS * nP )
  colnames(stockObjFunTable) <- stockTableColNames
  stockObjFunTable <- as.data.frame(stockObjFunTable)

  stockObjFunTable[,'simLabel']  <- simLabel
  stockObjFunTable[,'scenario']  <- scenarioName
  stockObjFunTable[,'mp']        <- mpName


  for( s in 1:nS )
    for( p in 1:nP )
    {
      rowIdx <- (s - 1) * nP + p
      stockObjFunTable[rowIdx,'species']        <- speciesNames[p]
      stockObjFunTable[rowIdx,'stock']          <- stockNames[p]
      stockObjFunTable[rowIdx,'objFun']         <- round(mean( objFun_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'wt_CbarStock']   <- avgCatWt
      stockObjFunTable[rowIdx,'of_CbarStock']   <- round(mean( Cbar_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'wt_hiDep']       <- hiDepBmsyWt
      stockObjFunTable[rowIdx,'of_hiDep']       <- round(mean( barHiDep_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'wt_loDep']       <- loDepBmsyWt
      stockObjFunTable[rowIdx,'of_loDep']       <- round(mean( barLoDep_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'wt_probDep']     <- probDepWt
      stockObjFunTable[rowIdx,'of_probHiDep']   <- round(mean( barProbHiDep_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'of_probLoDep']   <- round(mean( barProbLoDep_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'wt_closure']     <- closedWt
      stockObjFunTable[rowIdx,'of_closure']     <- round(mean( closedCount_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'wt_AAV']         <- AAVWt
      stockObjFunTable[rowIdx,'of_AAV']         <- round(mean( barAAV_isp[,s,p] ),2)
      stockObjFunTable[rowIdx,'wt_catDiff']     <- catDiffWt
      stockObjFunTable[rowIdx,'of_catDiff']     <- round(mean( barCatDiff_isp[,s,p]),2)
      stockObjFunTable[rowIdx,'wt_initCatDiff'] <- initCatDiffWt
      stockObjFunTable[rowIdx,'of_initCatDiff'] <- round(mean( barInitCatDiff_isp[,s,p]),2)
    }


  outList <- list(  wtTable     = wtTable,
                    cplxTable   = cplxOmniObjFunTable,
                    stockTable  = stockObjFunTable )


  outList
} # END .getOmniInfo()



# Envelopes of simulated assessment errors
calcREdists_AM <- function( simNum = 1,
                            obj = NULL,
                            groupFolder = "DLSurveys_.3tau_Long",
                            clearBadReps = FALSE,
                            abs = FALSE )
{

  if( is.null(obj))
  {
    .loadSim( simNum, groupFolder )
    obj <- blob
  }

  goodReps_isp  <- obj$goodReps_isp
  nReps         <- dim(goodReps_isp)[1]

  # Get biomass arrays
  SB_ispt        <- obj$om$SB_ispt
  VB_ispt        <- obj$om$vB_ispft[,,,2,]
  totB_ispt      <- obj$om$B_ispt
  retroSB_itspt  <- obj$mp$assess$retroSB_itspt
  retroUmsy_itsp <- blob$mp$assess$retroUmsy_itsp
  retroBmsy_itsp <- blob$mp$assess$retroBmsy_itsp
  omBmsy_sp      <- blob$rp[[1]]$FmsyRefPts$BeqFmsy_sp
  omYeqFmsy_sp   <- blob$rp[[1]]$FmsyRefPts$YeqFmsy_sp
  omUmsy_sp      <- omYeqFmsy_sp / omBmsy_sp



  ctlList <- obj$ctlList

  retroSB_itspt[retroSB_itspt < 0] <- NA
  # nReps     <- sum(goodReps)

  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT


  # Aggregate OM biomasses to match AM biomass
  if( ctlList$mp$data$spatialPooling )
  {
    newSB_ispt        <- apply( X = SB_ispt, FUN = sum, MARGIN = c(1,2,4), na.rm = T )

    omBmsy_sp_sum     <- apply( X = omBmsy_sp, FUN = sum, MARGIN = c(1) )
    omYeqFmsy_sp_sum  <- apply( X = omYeqFmsy_sp, FUN = sum, MARGIN = c(1) )
    omUmsy_sp_sum     <- omYeqFmsy_sp_sum / omBmsy_sp_sum

    omBmsy_sp <- omBmsy_sp[,1,drop = FALSE]
    omUmsy_sp <- omUmsy_sp[,1,drop = FALSE]
    omYeqFmsy_sp <- omYeqFmsy_sp[,1,drop = FALSE]

    # browser()

    omBmsy_sp[,1]   <- omBmsy_sp_sum
    omUmsy_sp[,1]   <- omUmsy_sp_sum
    omYeqFmsy_sp[,1]<- omYeqFmsy_sp_sum    

    SB_ispt    <- SB_ispt[,,1,,drop = FALSE]
    SB_ispt[,,1,] <- newSB_ispt

  }

  if( ctlList$mp$data$speciesPooling )
  {
    newSB_ispt        <- apply( X = SB_ispt, FUN = sum, MARGIN = c(1,3,4), na.rm = T )
  
    omBmsy_sp_sum     <- apply( X = omBmsy_sp, FUN = sum, MARGIN = c(2) )
    omYeqFmsy_sp_sum  <- apply( X = omYeqFmsy_sp, FUN = sum, MARGIN = c(2) )
    omUmsy_sp_sum     <- omYeqFmsy_sp_sum / omBmsy_sp_sum

    omBmsy_sp <- omBmsy_sp[1,,drop = FALSE]
    omUmsy_sp <- omUmsy_sp[1,,drop = FALSE]
    omYeqFmsy_sp <- omYeqFmsy_sp[1,,drop = FALSE]
    
    omBmsy_sp[1,]   <- omBmsy_sp_sum
    omUmsy_sp[1,]   <- omUmsy_sp_sum
    omYeqFmsy_sp[1,]<- omYeqFmsy_sp_sum



    SB_ispt    <- SB_ispt[,1,,,drop = FALSE]
    SB_ispt[,1,,] <- newSB_ispt
    
  }

  SB_ispt[SB_ispt == 0]     <- NA
  

  nSS     <- dim( SB_ispt)[2]
  nPP     <- dim( SB_ispt)[3]

  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fYear         <- obj$ctlList$opMod$fYear
  pT            <- obj$ctlList$opMod$pT
  scenario      <- obj$ctlList$ctl$scenarioName
  mp            <- obj$ctlList$ctl$mpName


  if( nPP == 1 )
    stockNames <- "Spatial Pooled"

  if( nSS == 1 )
    speciesNames <- "Species Pooled"

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")

  if( clearBadReps )
    for( s in 1:nSS )
      for( p in 1:nPP )
      {
        badIdx <- which(!goodReps_isp[,s,p])
        retroSB_itspt[badIdx,,s,p,] <- NA
        SB_ispt[badIdx,s,p,] <- NA
        retroBmsy_itsp[badIdx,,s,p] <- NA
        retroUmsy_itsp[badIdx,,s,p] <- NA
      }


  assErr_itspt <- array(NA, dim = c(nReps,pT,nSS,nPP,nT) )
  BmsyErr_itsp <- array(NA, dim = c(nReps,pT,nSS,nPP) )
  UmsyErr_itsp <- array(NA, dim = c(nReps,pT,nSS,nPP) )
  for( s in 1:nSS )
    for( p in 1:nPP )
    {
      for( tt in 1:pT )
      {
        # Now loop over projection years
        assErr_itspt[,tt,s,p,] <- (retroSB_itspt[1:nReps,tt,s,p,] - SB_ispt[1:nReps,s,p,])/SB_ispt[1:nReps,s,p,]
        BmsyErr_itsp[,tt,s,p]  <- (retroBmsy_itsp[1:nReps,tt,s,p] - omBmsy_sp[s,p])/omBmsy_sp[s,p]
        UmsyErr_itsp[,tt,s,p]  <- (retroUmsy_itsp[1:nReps,tt,s,p] - omUmsy_sp[s,p])/omUmsy_sp[s,p]
      }
    }


  absAssErr_itspt <- abs(assErr_itspt)
  MAREassErr_sp   <- apply( X = absAssErr_itspt, FUN = median, MARGIN = c(3,4), na.rm = T)
  MREassErr_sp    <- apply( X = assErr_itspt, FUN = median, MARGIN = c(3,4), na.rm = T)

  assErr_qspt <- apply( X = assErr_itspt,
                        FUN = quantile,
                        probs = c(0.025, 0.5, 0.975),
                        MARGIN = c(3,4,5),
                        na.rm = TRUE )

  absBmsyErr_itsp <- abs(BmsyErr_itsp)
  MAREBmsyErr_sp   <- apply( X = absBmsyErr_itsp, FUN = median, MARGIN = c(3,4), na.rm = T)
  MREBmsyErr_sp    <- apply( X = BmsyErr_itsp, FUN = median, MARGIN = c(3,4), na.rm = T)

  BmsyErr_qsp <- apply( X = BmsyErr_itsp,
                        FUN = quantile,
                        probs = c(0.025, 0.5, 0.975),
                        MARGIN = c(3,4),
                        na.rm = TRUE ) 

  absUmsyErr_itsp <- abs(UmsyErr_itsp)
  MAREUmsyErr_sp   <- apply( X = absUmsyErr_itsp, FUN = median, MARGIN = c(3,4), na.rm = T)
  MREUmsyErr_sp    <- apply( X = UmsyErr_itsp, FUN = median, MARGIN = c(3,4), na.rm = T)

  UmsyErr_qsp <- apply( X = UmsyErr_itsp,
                        FUN = quantile,
                        probs = c(0.025, 0.5, 0.975),
                        MARGIN = c(3,4),
                        na.rm = TRUE )  
  
  # Make an output table
  tabRows <- nSS * nPP
  colNames <- c(  "simID",
                  "Scenario",
                  "AM",
                  "Species",
                  "Stock",
                  "MARE_Bt",
                  "MRE_Bt",
                  "MARE_Bmsy",
                  "MRE_Bmsy",
                  "MARE_Umsy",
                  "MRE_Umsy" )

  mareTable <- matrix(NA, nrow = tabRows, ncol = length(colNames))
  colnames(mareTable) <- colNames

  mareTable <- as.data.frame(mareTable)

  rowIdx <- 1
  for( s in 1:nSS )
    for( p in 1:nPP )
    {
      mareTable[rowIdx,"simID"]     <- obj$simLabel
      mareTable[rowIdx,"Scenario"]  <- scenario
      mareTable[rowIdx,"AM"]        <- mp
      mareTable[rowIdx,"Species"]   <- speciesNames[s]
      mareTable[rowIdx,"Stock"]     <- stockNames[p]
      mareTable[rowIdx,"MARE_Bt"]   <- MAREassErr_sp[s,p]
      mareTable[rowIdx,"MRE_Bt"]    <- MREassErr_sp[s,p]
      mareTable[rowIdx,"MARE_Bmsy"] <- MAREBmsyErr_sp[s,p]
      mareTable[rowIdx,"MRE_Bmsy"]  <- MREBmsyErr_sp[s,p]
      mareTable[rowIdx,"MARE_Umsy"] <- MAREUmsyErr_sp[s,p]
      mareTable[rowIdx,"MRE_Umsy"]  <- MREUmsyErr_sp[s,p]

      rowIdx <- rowIdx + 1
    }

  return(mareTable)

} # END plotTulipAssError









