# # makeResultPlots.R
source("ms3R.R")

lapply(X = 2:7, FUN = calcLoss, groupFolder = "sens_hierSD")
lapply(X = 2:31, FUN = calcLoss, groupFolder = "sens_MSYCV")
lapply(X = 2:31, FUN = calcLoss, groupFolder = "sens_projObsErr")
lapply(X = 2:31, FUN = calcLoss, groupFolder = "sens_UmsySD")



makeResultPlots <- function(groupFolder, prefix = "parBat")
{

  batchDir <- file.path("./Outputs",groupFolder)
  outputDirList <- list.dirs(batchDir,recursive = FALSE, full.names = FALSE)
  outputDirList <- outputDirList[grepl("sim_",outputDirList)]
  outputDirList <- outputDirList[order(outputDirList)]

  # ourSimIndices <- which(grepl(prefix,outputDirList))

  # # abs and rel loss for SB
  lapply( X = ourSimIndices, FUN = plotLossTulip,
          groupFolder = groupFolder,
          lossType    = "rel",
          var         = "SB_ispt",
          save        = TRUE )

  lapply( X = ourSimIndices, FUN = plotLossTulip,
          groupFolder = groupFolder,
          lossType    = "abs",
          var         = "SB_ispt",
          save        = TRUE )

  # abs and rel loss for Catch
  lapply( X = ourSimIndices, FUN = plotLossTulip,
          groupFolder = groupFolder,
          lossType    = "rel",
          var         = "C_ispt",
          save        = TRUE )

  lapply( X = ourSimIndices, FUN = plotLossTulip,
          groupFolder = groupFolder,
          lossType    = "abs",
          var         = "C_ispt",
          save        = TRUE )

  saveBioPlots <- function( sim = 1,
                            groupFolder = groupFolder )
  {
    .loadSim(sim, folder=groupFolder)


    # First biomass tulips
    saveFile <- here("Outputs",groupFolder,blob$simLabel,"tulipBt.pdf")

    pdf( file = saveFile, width = 11, height = 8 )

    plotTulipBt(  obj = blob, nTrace = 3,
                  dep = TRUE,
                  ref = "B0",
                  Ct  = TRUE,
                  leg = FALSE )
    dev.off()

    # Then retrospective biomass
    goodReps_isp <- blob$goodReps_isp
    maxRep_sp <- apply( X = goodReps_isp, FUN = maxWhich, MARGIN =c(2,3))
    maxRep <- max(maxRep_sp)

    singleRepFolder <- here("Outputs",groupFolder,blob$simLabel,"singleReps")
    if(!dir.exists(singleRepFolder))
      dir.create(singleRepFolder)

    retroBioFolder <- file.path(singleRepFolder,"retroBt")
    if(!dir.exists(retroBioFolder))
      dir.create(retroBioFolder)

    for( k in 1:maxRep )
    {
      
      saveFile <- file.path(retroBioFolder,paste("retroBio",k,".pdf",sep = ""))
      pdf( file = saveFile, width = 11, height = 8 )

      plotRetroSBagg(  obj = blob, iRep = k )
      
      dev.off()
    }

  }

  # Now output tulip plots and retro SB
  lapply( X = ourSimIndices, FUN = saveBioPlots, groupFolder= groupFolder)

  # And do biomas and catch overlay plots
  lapply( X = ourSimIndices, FUN = plotTulipBtCtBaseSim, groupFolder = groupFolder, 
          save = TRUE, var = "SB_ispt")
  lapply( X = ourSimIndices, FUN = plotTulipBtCtBaseSim, save = TRUE, 
          groupFolder = groupFolder,
          var = "C_ispt")
  lapply( X = ourSimIndices, FUN = plotTulipAssError, 
          save = TRUE, groupFolder = groupFolder, obj = NULL )

  speciesNames <- c("Dover","English","Rock")
  stockNames <- c("HSHG","QCS","WCVI")
  scenNames  <- c("DERfit_HcMcAsSsIdx","DERfit_McAsSsIdx","DERfit_McSsIdx","DERfit_AsSsIdx")

  for( sIdx in 1:3)
    for( pIdx in 1:3 )
      for( scenIdx in 1:4)
      {
        saveFile <- paste("RetroBioSample_",scenNames[scenIdx],"_",speciesNames[sIdx],"_",stockNames[pIdx],".pdf",sep = "")
        savePath <- here::here("Outputs",groupFolder,"batchPlots")
        if(!dir.exists(savePath))
          dir.create(savePath)
        pdf( file = file.path(savePath,saveFile),
             width = 11, height = 8.5 )
        plotRetroBio_Scenario(  groupFolder = groupFolder,
                                iRep = 1,
                                scenName = scenNames[scenIdx],
                                prefix = "parBat",
                                species = speciesNames[sIdx],
                                stock = stockNames[pIdx] )
        stamp <- paste(scenNames[scenIdx],speciesNames[sIdx],stockNames[pIdx],sep =":")
        mtext( side = 1, text = stamp, adj = .65, cex = .8, col = "grey70", outer = TRUE )
        dev.off()
      }
}

makeResultPlots("sens_hierSD")
makeResultPlots("sens_MSYCV")
makeResultPlots("sens_projObsErr")
makeResultPlots("sens_UmsySD")


# speciesNames <- c("Dover","English","Rock")
# stockNames <- c("HSHG","QCS","WCVI")
# scenNames  <- c("DERfit_HcMcAsSsIdx","DERfit_McAsSsIdx","DERfit_McSsIdx","DERfit_AsSsIdx")

# for( sIdx in 1:3)
#   for( pIdx in 1:3 )
#     for( scenIdx in 1:4)
#     {
#       saveFile <- paste("RetroBioSample_",scenNames[scenIdx],"_",speciesNames[sIdx],"_",stockNames[pIdx],".pdf",sep = "")
#       savePath <- here::here("Outputs",groupFolder,"batchPlots")
#       if(!dir.exists(savePath))
#         dir.create(savePath)
#       pdf( file = file.path(savePath,saveFile),
#            width = 11, height = 8.5 )
#       plotRetroBio_Scenario(  groupFolder = groupFolder,
#                               iRep = 1,
#                               scenName = scenNames[scenIdx],
#                               prefix = "parBat",
#                               species = speciesNames[sIdx],
#                               stock = stockNames[pIdx] )
#       stamp <- paste(scenNames[scenIdx],speciesNames[sIdx],stockNames[pIdx],sep =":")
#       mtext( side = 1, text = stamp, adj = .65, cex = .8, col = "grey70", outer = TRUE )
#       dev.off()
#     }
