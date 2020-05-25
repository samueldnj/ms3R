# multiBatchScript.R

# A wrapper script to run multiple parallel batches
# with limited user input - likely worth
# creating two versions of this for the two 
# Mac Pros

library(parallel)
source("ms3R.r")
source("makeResultPlots.R")

batchControlFiles <- c( "sensRuns_BmsyCV.bch",
                        "sensRuns_UmsyCV.bch",
                        "sensRuns_hierSD.bch",
                        "sensRuns_switchqF.bch" )
                       
nBatchJobs <- length( batchControlFiles )


baseControlFiles  <- c( "simCtlFileBase.txt",
                        "simCtlFileBase.txt",
                        "simCtlFileBase.txt",
                        "simCtlFileBase.txt")
                        

prefixes       <- c(  "sensRuns_BmsyCV",
                      "sensRuns_UmsyCV",
                      "sensRuns_hierSD",
                      "sensRuns_switchqF", )

saveDirName       <- prefixes

baseLine <- c( "sim_OmniRun_Long" )
                        
nCores  <- rep(detectCores()-1, nBatchJobs)


# Create a directory to hold completed mseR batch
# jobs
for( s in saveDirName)
  if(!dir.exists( here::here("Outputs",s) ))
    dir.create(here::here("Outputs",s))


for( bIdx in 1:nBatchJobs )
{
  message(  "Starting batch for ",
            batchControlFiles[bIdx], ".\n", sep = "")
  # Make the batch
  makeBatch(  baseCtlFile = baseControlFiles[bIdx],
              batchCtlFile = batchControlFiles[bIdx])


  # Run the batch in parallel - new 
  # format should stop errors from
  # crashing runParBatch
  tryPar <- try(.runBatchJob( par = TRUE,
                              nCores = nCores[bIdx],
                              prefix = prefixes[bIdx] ) )


  # Calculate loss functions - need to detect our
  # batch sims here
  outputDirList <- list.dirs("./Outputs",recursive = FALSE, full.names = FALSE)
  outputDirList <- outputDirList[grepl("sim_",outputDirList)]
  outputDirList <- outputDirList[order(outputDirList)]

  ourSimIndices <- which(grepl(prefixes[bIdx],outputDirList))
  
  lapply( X = ourSimIndices, FUN = calcLoss, baseline = baseLine, groupFolder = "" )
  gc()

  # Now copy fits to a subfolder location so
  # it doesn't grow the size of the copied WDs, and then clear
  # the sims
  destFolder <- file.path(".","Outputs",saveDirName[bIdx])

  projFolderContents <- list.dirs( "./Outputs", 
                                    recursive = FALSE,
                                    full.names = FALSE )

  # Restrict to the fits with this prefix
  projFolderContents <- projFolderContents[grepl("sim_",projFolderContents)]
  projFolderContents <- projFolderContents[grepl(prefixes[bIdx],projFolderContents)]

  oldPath <- here::here("Outputs/",projFolderContents)
  newPath <- file.path(destFolder,projFolderContents)

  for( k in 1:length(oldPath))
    fs::dir_copy(oldPath[k], newPath[k], overwrite = TRUE)

  fs::dir_delete(oldPath)

  # Run makeResultsPlots
  makeResultPlots(saveDirName[bIdx])

  
  message(  "Parallel batch complete for ",
            batchControlFiles[bIdx],
            " and fits copied and tidied up for next batch.\n",sep = "")
}

  message(  "All batch jobs complete!\n" )
