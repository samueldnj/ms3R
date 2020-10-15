# # makeResultPlots.R
source("ms3R.R")

lapply(X = 2:21, FUN = calcLoss, groupFolder = "DERTACs_reruns_Oct10")
lapply(X = 2:7, FUN = calcLoss, groupFolder = "sens_hierSD")
lapply(X = 2:31, FUN = calcLoss, groupFolder = "sens_MSYCV")
lapply(X = 2:31, FUN = calcLoss, groupFolder = "sens_projObsErr")
lapply(X = 2:31, FUN = calcLoss, groupFolder = "sens_UmsySD")



makeResultPlots <- function(groupFolder, prefix = "parBat", retroBioBatchPlot = TRUE)
{

  batchDir <- file.path("./Outputs",groupFolder)
  outputDirList <- list.dirs(batchDir,recursive = FALSE, full.names = FALSE)
  outputDirList <- outputDirList[grepl("sim_",outputDirList)]
  outputDirList <- outputDirList[order(outputDirList)]

  ourSimIndices <- which(grepl(prefix,outputDirList))

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

  if( retroBioBatchPlot)
  {

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
}


makeResultPlots("DERTACs_reruns_Oct10",retroBioBatchPlot = TRUE)
makeResultPlots("sens_hierSD",retroBioBatchPlot = FALSE)
makeResultPlots("sens_MSYCV", retroBioBatchPlot = FALSE)
makeResultPlots("sens_projObsErr",retroBioBatchPlot = FALSE)
makeResultPlots("sens_UmsySD", retroBioBatchPlot = FALSE)


# plotRankDists_sp( groupFolder = "sens_projObsErr",
#                   prefix = "parBat",
#                   lossType = "abs",
#                   var = "C_ispt",
#                   period = 73:82,
#                   lossList = NULL,
#                   dim1 = 1:3,   # species (D,E,R, DP)
#                   dim2 = 1:3,   # stocks (H,Q,W, CW)
#                   qProbs = c(.025,.5,.975),
#                   refPts = "MSrefPts",
#                   AMlabs = rev(c( SS = "singleStock",
#                                   HMS = "hierMultiStock",
#                                   SpecPool = "speciesPooling",
#                                   SpatPool = "spatialPooling",
#                                   TotPool = "totalAgg" ) ),
#                   scenLabs = c( Rich_obsErr0.1  = "Rich_obsErr0.1",
#                                 Poor_obsErr0.1  = "Poor_obsErr0.1",
#                                 Rich_obsErr0.5  = "Rich_obsErr0.5",
#                                 Poor_obsErr0.5  = "Poor_obsErr0.5",
#                                 Rich_obsErr1.0  = "Rich_obsErr1.0",
#                                 Poor_obsErr1.0  = "Poor_obsErr1.0"),
#                   clearBadReps = TRUE,
#                   minSampSize = 10,
#                   rotxlabs = 0,
#                   xlab = TRUE,
#                   vertLines = TRUE  )

# plotRankDists(  groupFolder = "sens_projObsErr",
#                 prefix = "parBat",
#                 lossType = "abs",
#                 var = "C_ispt",
#                 period = 73:82,
#                 lossList = NULL,
#                 dim1 = 1:3,   # species (D,E,R, DP)
#                 dim2 = 1:3,   # stocks (H,Q,W, CW)
#                 qProbs = c(.025,.5,.975),
#                 refPts = "MSrefPts",
#                 AMlabs = rev(c( SS = "singleStock",
#                                 HMS = "hierMultiStock",
#                                 SpecPool = "speciesPooling",
#                                 SpatPool = "spatialPooling",
#                                 TotPool = "totalAgg" ) ),
#                 scenLabs = c( Rich_obsErr0.1  = "Rich_obsErr0.1",
#                               Poor_obsErr0.1  = "Poor_obsErr0.1",
#                               Rich_obsErr0.5  = "Rich_obsErr0.5",
#                               Poor_obsErr0.5  = "Poor_obsErr0.5",
#                               Rich_obsErr1.0  = "Rich_obsErr1.0",
#                               Poor_obsErr1.0  = "Poor_obsErr1.0"),
#                 clearBadReps = FALSE,
#                 minSampSize = 49,
#                 rotxlabs = 0,
#                 xlab = TRUE,
#                 vertLines = TRUE,
#                 nGoodReps = 25 )

# plotBatchLossDists_Scenario(  groupFolder = "sens_projObsErr",
#                               prefix = "parBat",
#                               lossType = "abs",
#                               var = "C_ispt",
#                               period = 73:82,
#                               lossList = NULL,
#                               dim1 = 1:3,   # species (D,E,R, DP)
#                               dim2 = 1:3,   # stocks (H,Q,W, CW)
#                               qProbs = c(.025,.5,.975),
#                               refPts = "MSrefPts",
#                               AMlabs = rev(c( SS = "singleStock",
#                                               HMS = "hierMultiStock",
#                                               SpecPool = "speciesPooling",
#                                               SpatPool = "spatialPooling",
#                                               TotPool = "totalAgg" ) ),
#                               scenLabs = c( Rich_obsErr0.1  = "Rich_obsErr0.1",
#                                             Poor_obsErr0.1  = "Poor_obsErr0.1",
#                                             Rich_obsErr0.5  = "Rich_obsErr0.5",
#                                             Poor_obsErr0.5  = "Poor_obsErr0.5",
#                                             Rich_obsErr1.0  = "Rich_obsErr1.0",
#                                             Poor_obsErr1.0  = "Poor_obsErr1.0"),
#                               clearBadReps = TRUE,
#                               minSampSize = 40,
#                               rotxlabs = 45,
#                               xlab = TRUE,
#                               vertLines = TRUE  )

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
