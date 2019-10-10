# loadFit()
# Loads the nominated fit reports object into memory, 
# so that plot functions can be called
# inputs:   fit=ordinal indicator of sim in project folder,
#               or character giving the name of the folder
# ouputs:   NULL
# usage:    Prior to plotting simulation outputs
.loadFit <- function( fit = 1, groupFolder = "." )
{
  fitFolder <- here::here("Outputs","fits",groupFolder)

  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path=fitFolder,full.names = FALSE,
                        recursive=FALSE)
  # Restrict to fit_ folders, pick the nominated simulation
  fitList <- dirList[grep(pattern="fit",x=dirList)]
  folder <- fitList[fit]

  # Load the nominated blob
  reportsFileName <- paste(folder,".RData",sep="")
  reportsPath <- file.path(fitFolder,folder,reportsFileName)
  load ( file = reportsPath )

  reports$repOpt <- renameReportArrays(reports$repOpt,reports$data)

  # Update the path object - not sure if this is necessary
  reports$path <- file.path(fitFolder, folder) 

  # Assign to global environment
  assign( "reports",reports,pos=1 )

  message("(.loadFit) Reports in ", folder, " loaded\n", sep="" )

  return( file.path(fitFolder, folder) )
} # END .loadFit()