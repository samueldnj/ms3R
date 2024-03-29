# <><><><><><><><><><><><><><><><><><><><><><><><><><><>
# tools.R
#
# Tools for the ms3R closed loop simulation platform.
#
# Author: SDN Johnson
# Date: October 10, 2019
#
# Last Update: October 10, 2019
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><>

# redoEconYield()
# recalculates economic yield parameters, good for
# appending coastwide curves post-simulation, avoiding
# optim calls between sims
redoEconYield <- function( sim, groupFolder = "",
                           iRep = 1, cwMEY = TRUE )
{
  source("ms3R.R")
  path <- .loadSim(sim, groupFolder = groupFolder)

  
  obj <- blob
  # Pull things we need to do it
  repObj  <- obj$ctlList$opMod$histRpt
  opMod   <- obj$ctlList$opMod

  opMod$cwMEY <- cwMEY

  rp      <- obj$rp[[iRep]]

  repObj$opMod    <- opMod
  repObj$om       <- list()
  repObj$om$C_spt <- obj$om$C_ispt[iRep,,,]
  repObj$om$E_pft <- obj$om$E_ipft[iRep,,,]

  econYieldCurves <- .calcEconomicYieldCurves(  repObj, 
                                                rp$refCurves$EffCurves,
                                                rp$FmsyRefPts )

  obj$rp[[iRep]]$EmeyRefPts <- econYieldCurves

  blob <- obj
  filename <- paste0(basename(path),".Rdata")
  save(blob, file = file.path(path,filename))
} # END redoEconYield


# getSplineVal()
# Takes in 2 vectors of the same length, fits
# a spline, and returns a given point value for
# the dependent variable at x = p
getSplineVal <- function( x, y, p)
{
  xySplineFun <- splinefun( x=x, y=y )  

  y.p <- xySplineFun(p)
  y.p
}

updateGoodReps <- function( groupFolder = "DLSurveys7_.5tau", 
                            prefix = "MPgrid" )
{
  # Load batch info
  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
                filter( grepl( prefix, simLabel ) )
  
  # Loop and
  nSims <- nrow(info.df)
  for( i in 1:nSims )
  {
    # 1. Load simulations
    blobFileName <- info.df$simLabel[i]
    simFolderPath <- here::here("Outputs",groupFolder,blobFileName)

    blobPath <- file.path(simFolderPath,paste(blobFileName,".RData",sep = ""))

    load(blobPath)

    # 2. Update goodReps structure
    repDim  <- dim(blob$om$SB_ispt)[1]
    totReps <- blob$nSims
    nS      <- blob$om$nS
    nP      <- blob$om$nP
    nT      <- blob$om$nT
    pT      <- blob$om$pT

    probPosSD_isp  <- apply( X = blob$mp$assess$posSDs_itsp, FUN = mean, MARGIN = c(1,3,4),na.rm = TRUE)
    probPDHess_isp <- apply( X = blob$mp$assess$pdHess_itsp, FUN = mean, MARGIN = c(1,3,4),na.rm = TRUE)

    probPosSD_isp[is.nan(probPosSD_isp)] <- 0
    probPDHess_isp[is.nan(probPosSD_isp)] <- 0

    goodReps_isp <- array(FALSE, dim = c(totReps,nS,nP))

    # Need the loop to preserve the array structure
    for( j in 1:repDim)
      for( s in 1:nS )
        for( p in 1:nP )
        {
          if( ( probPosSD_isp[j,s,p] >= .95 & probPDHess_isp[j,s,p] >= .95 ) )
          {
            goodReps_isp[j,s,p] <- TRUE 
          }
        }


    blob$goodReps_isp <- goodReps_isp
    # 3. Save over old blob

    save(blob, file = blobPath)
  }

  message("All blobs updated with new goodReps definition.")
}

# .saveBlob()
# Saves the blob to a sim folder
.saveBlob <- function( blob, ctlTable, outFolder )
{
  # save output to project folder
  # First, if a folder name isn't nominated, create a default sim folder
  if ( is.null(outFolder) )
  {
    stamp <- paste( format(Sys.time(),format="%d%m%Y%H%M%S" ),sep="" )
    folder <- paste ( "sim_",stamp, sep = "" )
  } else folder <- paste ("sim_", outFolder, sep = "")
  
  # Now paste together the path to the folder and create it
  path <- here::here("Outputs",folder)
  dir.create ( path )
  message( "\ (.saveFit) Created simulation folder ",folder,"in ./Outputs/.\n" )
  
  # Save path so we can access it later
  blob$path     <- path
  blob$simLabel <- folder

  # Save blob
  save(blob,file = file.path(path,paste(folder,".RData",sep="")))

  # Create a quick to read info file for the sim folder
  .makeInfoFile(blob)

  graphics.off()


  # Calculate and save stats table
  perfStats <- .simPerfStats(obj = blob)
  write.csv( perfStats, file = file.path(path,"simPerfStats.csv"))


  # Copy AM file
  if(blob$ctlList$mp$assess$method == "hierProd" )
    file.copy( "hierProd.cpp", file.path(path, "hierProd.cpp") )

  # Copy control file to sim folder for posterity
  cat(  "# simCtlFile.txt, written to ", folder, "on ", Sys.time(),"\n", sep = "", 
        file = file.path(path,"simCtlFile.txt"))
  write.table(  ctlTable,
                file = file.path(path,"simCtlFile.txt"),
                row.names = FALSE,
                quote = FALSE, qmethod = "double",
                append = TRUE )
  cat(  "# <End File>", sep = "", 
        file = file.path(path,"simCtlFile.txt"),
        append = TRUE)


  # Make html sim report
  .makeSimReport( simNum = folder, groupFolder = "" )

  message("Simulation complete!")
  
} # END .saveBlob

# .makeSimReport()
# Makes a .html report of the simulation
# showing basic performance plots
# and tables
.makeSimReport <- function( simNum, groupFolder = "",
                            tidy = TRUE )
{
  simFolder <- here::here("Outputs",groupFolder)

  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path=simFolder,full.names = FALSE,
                        recursive=FALSE)
  # Restrict to fit_ folders, pick the nominated simulation
  simList <- dirList[grep(pattern="sim",x=dirList)]

  if( is.character(simNum) )
    folder <- simList[grepl(x = simList, pattern = simNum) ]
  else
    folder <- simList[simNum]

  # Load the nominated blob
  simFileName <- paste(folder,".RData",sep="")
  simPath <- file.path(simFolder,folder,simFileName)
  simFolderPath <- here::here("Outputs",groupFolder,folder)

  # Create parameter list for rendering the document
  params <- list( rootDir= simFolderPath,
                  RdataFile = simFileName)
  # Make an output file name
  outFile <- paste( "simReport.html", sep = "")

  # Render
  rmarkdown::render(  input = here::here("docs/reports","simReportTemplate.Rmd"), 
                      output_file = outFile,
                      output_dir = simFolderPath,
                      params = params,
                      clean = tidy,
                      envir = new.env(),
                      output_format = "bookdown::html_document2" )

  # remove temporary files
  dataReportFiles <- "simReport_files"
  if( tidy )
    unlink(file.path(simFolderPath,dataReportFiles), recursive = TRUE)
} # END .makeSimReport()

# .makeBatchReport()
# Makes a .html report of the simulation
# showing basic performance plots
# and tables
.makeBatchReport <- function( batchFolder = "",
                              prefix = "" )
{
  simFolder <- here::here("Outputs",batchFolder)

  # Create parameter list for rendering the document
  params <- list( batchDir = simFolder,
                  groupFolder = batchFolder,
                  prefix = prefix  )


  # Make an output file name
  outFile <- paste( "batchReport.html", sep = "")
  outFolder <- file.path(simFolder,"batchReports")

  if(!dir.exists(outFolder))
    dir.create(outFolder)

  # Render
  rmarkdown::render(  input = here::here("docs/reports","batchReportTemplate.Rmd"), 
                      output_file = outFile,
                      output_dir = outFolder,
                      params = params,
                      envir = new.env(),
                      output_format = "bookdown::html_document2" )

  # remove temporary files
  dataReportFiles <- "batchReport_files"
  unlink(file.path(outFolder,dataReportFiles), recursive = TRUE)
} # END .makeBatchReport()

# .makeBatchReport()
# Makes a .html report of the simulation
# showing basic performance plots
# and tables
.makeBatchLossReport <- function( batchFolder = "",
                                  prefix = "" )
{
  simFolder <- here::here("Outputs",batchFolder)

  # Create parameter list for rendering the document
  params <- list( batchDir = simFolder,
                  groupFolder = batchFolder,
                  prefix = prefix )

  # Make an output file name
  outFile <- paste( "batchLossReport.html", sep = "")
  outFolder <- file.path(simFolder,"batchReports")

  if(!dir.exists(outFolder))
    dir.create(outFolder)

  # Render
  rmarkdown::render(  input = here::here("docs/reports","batchLossReportTemplate.Rmd"), 
                      output_file = outFile,
                      output_dir = outFolder,
                      params = params,
                      envir = new.env(),
                      output_format = "bookdown::html_document2" )

  # remove temporary files
  dataReportFiles <- "batchLossReport_files"
  unlink(file.path(outFolder,dataReportFiles), recursive = TRUE)
} # END .makeBatchLossReport()

# readBatchInfo()
# Loads a data.frame of info files for
# a given groupFolder
readBatchInfo <- function(batchDir = here("Outputs") )
{
  
  batchFitDirs <- list.dirs(batchDir,recursive = FALSE)
  batchFitDirs <- batchFitDirs[grepl(x = batchFitDirs, pattern = "sim_")]
  nFits <- length(batchFitDirs)

  # Read in the info files
  infoFiles <- file.path(batchFitDirs,"infoFile.txt")
  info.df <- lapply(X = infoFiles, FUN = lisread )
  info.df <- as.data.frame(do.call(rbind, info.df)) %>%
              mutate_all(as.character)


  info.df
}

# .makeInfoFile()
# Creates an info list for the
# fit object in the argument
.makeInfoFile <- function( obj )
{
  infoList <- list()

  # Now start populating
  infoList$scenario <- obj$ctlList$ctl$scenarioName
  infoList$mp       <- obj$ctlList$ctl$mpName
  infoList$path     <- obj$path
  infoList$simLabel <- obj$simLabel
  infoList$time     <- obj$time

  outFile <- file.path(obj$path, "infoFile.txt")


  # Create an output file
  cat('## Info file for ms3R simulation\n', file = outFile, append = FALSE, sep = "")
  cat('## Written ', Sys.time(), "\n", file = outFile, append = TRUE, sep = "")
  cat("", file = outFile, append = TRUE, sep = "")

  # Now start writing info list
  for( lIdx in 1:length(infoList) )
  {
    cat("# ", names(infoList)[lIdx], "\n", file = outFile, append = TRUE, sep = "")
    cat( infoList[[lIdx]], "\n", file = outFile, append = TRUE, sep = "")
    cat("", file = outFile, append = TRUE, sep = "")
  }

  cat("## End File\n", file = outFile, append = TRUE, sep = "")

  message("Info file created at ", outFile, ".\n", sep = "")
} # END .makeInfoFile()

# loadSim()
# Loads the nominated sim reports object into memory, 
# so that plot functions can be called
# inputs:   fit=ordinal indicator of sim in project folder,
#               or character giving the name of the folder
# ouputs:   NULL
# usage:    Prior to plotting simulation outputs
.loadSim <- function( sim = 1, 
                      baseDir = "Outputs",
                      groupFolder = "",
                      folder = "" )
{
  simFolder <- here::here(file.path(baseDir,groupFolder,folder))

  # List directories in project folder, remove "." from list
  if(is.numeric(sim))
  {
    dirList <- list.dirs (path=simFolder,full.names = FALSE,
                          recursive=FALSE)
    # Restrict to sim_ folders, pick the nominated simulation
    simList <- dirList[grep(pattern="sim_",x=dirList)]
    # Pick of folder number
    folder <- simList[sim]
  }
  else folder <- sim

  # Load the nominated blob
  simFileName <- paste(folder,".RData",sep="")
  simPath <- file.path(simFolder,folder,simFileName)
  load ( file = simPath )


  # Assign to global environment
  assign( "blob",blob,pos=1 )

  message(" (.loadSim) Simulation in ", folder, " loaded and placed in global environment as blob.\n", sep="" )

  return( file.path(simFolder,folder) )
} # END .loadSim()

# loadLoss()
# Loads the nominated loss reports object into memory, 
# so that plot functions can be called
# inputs:   sim=ordinal indicator of sim in project folder,
#               or character giving the name of the folder
# ouputs:   NULL
# usage:    Prior to plotting simulation outputs
.loadLoss <- function( sim = 1, folder = "" )
{
  simFolder <- here::here(file.path("Outputs",folder))


  # List directories in project folder, remove "." from list
  if(is.numeric(sim))
  {
    dirList <- list.dirs (path=simFolder,full.names = FALSE,
                          recursive=FALSE)
    # Restrict to sim_ folders, pick the nominated simulation
    simList <- dirList[grep(pattern="sim_",x=dirList)]
    # Pick of folder number
    folder <- simList[sim]
  }
  else folder <- sim

  # Load the nominated blob
  lossFileName <- paste("loss.RData",sep="")
  lossPath <- file.path(simFolder,folder,lossFileName)
  load ( file = lossPath )


  message(" (.loadLoss) Simulation loss calculations in ", folder, " loaded.\n", sep="" )

  return(outList)
} # END .loadLoss()


## -------------------------------------------- ##
# ---  Some helper functions from hierSCAL   --- #
## -------------------------------------------- ##

# loadFit()
# Loads the nominated fit reports object into memory, 
# so that plot functions can be called
# inputs:   fit=ordinal indicator of sim in project folder,
#               or character giving the name of the folder
# ouputs:   NULL
# usage:    Prior to plotting simulation outputs
.loadFit <- function( fit = 1,
                      baseDir = "history" )
{
  fitFolder <- here::here(baseDir)

  # List directories in project folder, remove "." from list
  if(is.numeric(fit))
  {
    dirList <- list.dirs (path=fitFolder,full.names = FALSE,
                          recursive=FALSE)
    # Restrict to fit_ folders, pick the nominated simulation
    fitList <- dirList[grep(pattern="fit_",x=dirList)]
    # Pick of folder number
    folder <- fitList[fit]
  }
  else folder <- fit

  # Load the nominated blob
  reportsFileName <- paste(folder,".RData",sep="")
  reportsPath <- file.path(fitFolder,folder,reportsFileName)
  load ( file = reportsPath )

  message(" (.loadFit) Reports in ", folder, " loaded\n", sep="" )

  return( reports )
} # END .loadFit()

# renameReportArrays()
# Updates the dimension names of the arrays in the 
# report lists, as a way of making plot code more efficient later.
renameReportArrays <- function( repObj = repInit, datObj = data )
{
  # Just go down the list, but first do the objects with the same names
  repNames <- names(repObj)
  datNames <- names(datObj)

  bothNames <- repNames[ repNames %in% datNames ]

  for( itemName in bothNames )
  {
    dimnames(repObj[[itemName]]) <- dimnames(datObj[[itemName]])
  }

  # Recover names
  specNames   <- dimnames(datObj$I_spft)[[1]]
  stockNames  <- dimnames(datObj$I_spft)[[2]]
  gearNames   <- dimnames(datObj$I_spft)[[3]]
  yearNames   <- dimnames(datObj$I_spft)[[4]]
  lenNames    <- dimnames(datObj$len_lxspft)[[1]]
  ageNames    <- dimnames(datObj$age_axspft)[[1]]
  sexNames    <- dimnames(datObj$age_axspft)[[2]]
  fltGrpNames <- c("commTrawl","HSAss","Synoptic")

  nS <- length(specNames)
  nP <- length(stockNames)

  # Create species/stock names vector
  specStock <- character(length = nS*nP)
  for( s in 1:nS )
    for( p in 1:nP )
      specStock[(s-1) * nP + p] <- paste( specNames[s],"_", stockNames[p],sep = "")


  # Ok, that's the data taken care of. There are still all the
  # new arrays that we created
  # Predicted data
  dimnames(repObj$I_spft_hat)       <- dimnames(datObj$I_spft)
  dimnames(repObj$aDist_axspft_hat) <- dimnames(datObj$age_axspft)
  dimnames(repObj$lDist_lxspft_hat) <- dimnames(datObj$len_lxspft)
  # State arrays
  dimnames(repObj$B_axspt)          <- dimnames(datObj$age_axspft)[c(1:4,6)]
  dimnames(repObj$N_axspt)          <- dimnames(datObj$age_axspft)[c(1:4,6)]
  dimnames(repObj$B_spt)            <- dimnames(datObj$age_axspft)[c(3,4,6)]
  dimnames(repObj$R_spt)            <- dimnames(datObj$age_axspft)[c(3,4,6)]
  dimnames(repObj$SB_spt)           <- dimnames(datObj$age_axspft)[c(3,4,6)]
  dimnames(repObj$vB_spft)          <- list(  species = specNames,
                                              stock = stockNames,
                                              fleet = gearNames,
                                              year = yearNames )
  dimnames(repObj$predC_spft)       <- dimnames(datObj$age_axspft)[c(3:6)]
  dimnames(repObj$predCw_spft)      <- dimnames(datObj$age_axspft)[c(3:6)]
  dimnames(repObj$C_axspft)         <- dimnames(datObj$age_axspft)
  dimnames(repObj$Cw_axspft)        <- dimnames(datObj$age_axspft)
  dimnames(repObj$F_spft)           <- dimnames(datObj$age_axspft)[c(3:6)]
  dimnames(repObj$Z_axspt)          <- dimnames(datObj$age_axspft)[c(1:4,6)]
  # Biological parameters
  dimnames(repObj$R0_sp)          <- dimnames(datObj$age_axspft)[c(3,4)]
  dimnames(repObj$B0_sp)          <- dimnames(datObj$age_axspft)[c(3,4)]
  dimnames(repObj$h_sp)           <- dimnames(datObj$age_axspft)[c(3,4)]
  dimnames(repObj$phi_sp)         <- dimnames(datObj$age_axspft)[c(3,4)]
  dimnames(repObj$reca_sp)        <- dimnames(datObj$age_axspft)[c(3,4)]
  dimnames(repObj$recb_sp)        <- dimnames(datObj$age_axspft)[c(3,4)]
  dimnames(repObj$sigmaR_sp)      <- dimnames(datObj$age_axspft)[c(3,4)]
  dimnames(repObj$omegaRmat_spt)  <- list(  specStock = specStock, 
                                            year = yearNames )
  dimnames(repObj$recCorrMat_sp)  <- list(  specStock = specStock,
                                            specStock = specStock )


  # Observation models
  dimnames(repObj$q_spf)        <- dimnames(datObj$age_axspft)[c(3:5)]  
  dimnames(repObj$q_spft)       <- dimnames(datObj$age_axspft)[c(3:6)]  
  dimnames(repObj$tau2Obs_spf)  <- list(  species = specNames,
                                          stock = stockNames,
                                          fleet = gearNames )  
  dimnames(repObj$tauObs_spf)   <- list(  species = specNames,
                                          stock = stockNames,
                                          fleet = gearNames )
  dimnames(repObj$sel_lspft)    <- list(  len = lenNames, 
                                          species = specNames,
                                          stock = stockNames,
                                          fleet = gearNames,
                                          year = yearNames )
  dimnames(repObj$sel_axspft)     <- list(  age = ageNames, 
                                            sex = sexNames,
                                            species = specNames,
                                            stock = stockNames,
                                            fleet = gearNames,
                                            year = yearNames )


  dimnames(repObj$ageRes_axspft)   <- dimnames(datObj$age_axspft)
  dimnames(repObj$tau2Age_spf)   <- list( species = specNames,
                                          stock = stockNames,
                                          fleet = gearNames )  
  dimnames(repObj$lenRes_lxspft)   <- dimnames(datObj$len_lxspft)
  dimnames(repObj$tau2Len_spf)   <- list( species = specNames,
                                          stock = stockNames,
                                          fleet = gearNames )  

  dimnames(repObj$probLenAge_laxsp)      <- list( len = lenNames,
                                                  age = ageNames,
                                                  sex = sexNames,
                                                  species = specNames,
                                                  stock = stockNames )

  dimnames(repObj$probAgeLen_laxspft)      <- list( len = lenNames,
                                                    age = ageNames,
                                                    sex = sexNames,
                                                    species = specNames,
                                                    stock = stockNames,
                                                    fleet = gearNames,
                                                    year = yearNames )

  dimnames(repObj$ageAtLenResids_laxspft)      <- list( len = lenNames,
                                                        age = ageNames,
                                                        sex = sexNames,
                                                        species = specNames,
                                                        stock = stockNames,
                                                        fleet = gearNames,
                                                        year = yearNames )

  dimnames(repObj$lenAge_axsp) <- list( age = ageNames,
                                        sex = sexNames,
                                        species = specNames,
                                        stock = stockNames )

  # Growth model parameters
  # vectors
  names(repObj$A1_s)      <- specNames
  names(repObj$A2_s)      <- specNames
  names(repObj$L1_s)      <- specNames
  
  names(repObj$sigmaLa_s) <- specNames
  names(repObj$sigmaLb_s) <- specNames
  # arrays
  dimnames(repObj$L2_sp)    <- list(  species = specNames,
                                      stock = stockNames )
  dimnames(repObj$L2_xs)    <- list(  sex = sexNames,
                                      species = specNames )
  dimnames(repObj$vonK_sp)  <- list(  species = specNames,
                                      stock = stockNames  )

  dimnames(repObj$vonK_xs)  <- list(  sex = sexNames,
                                      species = specNames )

  dimnames(repObj$L1_xsp)   <- list(  sex = sexNames,
                                      species = specNames,
                                      stock = stockNames )
  dimnames(repObj$L2_xsp)   <- list(  sex = sexNames,
                                      species = specNames,
                                      stock = stockNames )
  dimnames(repObj$vonK_xsp) <- list(  sex = sexNames,
                                      species = specNames,
                                      stock = stockNames )
  dimnames(repObj$M_xsp)    <- list(  sex = sexNames,
                                      species = specNames,
                                      stock = stockNames )




  return(repObj)
} # END renameReportArrays()


#-----------------------------------------------------------------------------##
#-- Helper Functions from mseR (some HIDDEN, e.g., .foo)                    --##
#-----------------------------------------------------------------------------##

lisread <- function( fname,quiet=TRUE )
{
  # lisread: Function to read a list of data objects from a file.
  # The initial characters "##" denote a comment line (ignored).
  # The initial characters "# " denote a variable name.
  # All other lines must contain scalars or vectors of numbers.
  # Furthermore, all rows in a matrix must contain the same number of
  # columns. Row and column vectors are not converted to matrices.
  #
  # fname  : File name.
  # quiet  : If true, shut up about reporting progress.
  # result : List object with components named in the file.

  # Original functions courtesy of Jon Schnute.
  # Modifications by A.R. Kronlund.

  lis2var <- function( x )
  {
    # lis2var: Makes global variables from the components of a list
    # x      : list object with named components.
    # result : global variables with names and contents extracted from x.

    namx <- names( x )
    nx <- length( namx )
    if (nx > 0) for (i in 1:nx)
    {
      if (namx[i] != "")
      {
        cmd <- paste( namx[i],"<<- x[[i]]" )
        eval( parse(text=cmd) )
      }
    }
    namx[namx != ""]
  }

  # numvecX functions:
  #
  # Function to convert a single string with white spaces into a numeric
  # vector. The string is parsed into separate components and converted
  # to char and numeric. A direct conversion to numeric fails.

  numvec <- function( x )
  {
    # Deprecated.
    xp <- parse( text=x,white=T )
    xc <- as.character( xp )
    as.numeric( xc )
  }

  numvec2 <- function( x )
  {
    # Patch required for S6.0 where white=T option is defunct in parse.
    # Deprecated:  xp <- parse( text=x,white=T )
    # ARK 30-Oct-03: R text connections get used up, must open/close.
    tc <- textConnection( x )
    xp <- scan( tc )
    close( tc )
    xc <- as.character( xp )
    as.numeric( xc )
  }

  numvec3 <- function( x,quiet )
  {
    # ARK 12-Jan-06: Patch to read characters because Rashmi asked nicely.
    # This is a largely untested hack, no expressed or implied warrantee.
    
    tmpwarn <- options( "warn" )
    options( warn=-1 )
    tc <- textConnection( x )
    xp <- scan( tc, what="character",quiet=quiet )
    close( tc )
    xc <- as.character( xp )
    if ( !all(is.na(as.numeric(xc))) )
      xc <- as.numeric( xc )
      
    options( tmpwarn )
    xc
  }
  
  #------------------------------------------------------------------#

  file <- scan( fname, what=character(), sep="\n", quiet=quiet )

  f2 <- file[ regexpr("##",file)==-1 ]           # remove comments
  nf2 <- length( f2 )                            # number of lines
  llab <- regexpr( "#",f2 )==1                   # identifies label lines
  vlab <- substring( f2[llab],3 )                # variable labels

  # ARK 30-Oct-03 R does not coerce logical to character for grep.
  ilab <- grep( "TRUE",as.character(llab) )      # label indices

  nvar <- length( vlab )                         # number of variables
  
  # ARK 19-Jan-10 When there is only one varaible in a file, the original
  # code does not work, namely:
  #    nrow <- c( ilab[2:nvar],nf2+1) - ilab - 1
  # returns an NA because the ilab vector is of length 1.
  #
  # Calculate the number of line for each variable.
  if ( nvar == 1 )
    nrow <- (nf2+1) - ilab - 1
  else  
    nrow <- c( ilab[2:nvar],nf2+1 ) - ilab - 1
    
  zout <- list( NULL )

  for ( i in 1:nvar )
  {
    i1 <- ilab[i] + 1                            # line of first var element
    i2 <- i1 + nrow[i] - 1                       # line of last  var element
    zstr <- paste(f2[i1:i2],collapse=" ")
#    zvec <- numvec2(zstr)                        # numeric vector
    zvec <- numvec3(zstr,quiet)                  # numeric or character vector

    nz <- length(zvec)
    zrow <- nrow[i]
    zcol <- nz / zrow                            # dimensions
    if ( (zrow>1) & (zcol>1) )                   # a true matrix
    {
      zvec <- matrix( zvec,nrow=zrow,ncol=zcol,byrow=T )
#      print( vlab[i] )
#      print( zvec )
#      scan()
    }
    
    zout[[i]] <- zvec
    if ( !quiet )
      cat( "vlab = ", vlab[i], "\n" )
  }
  names(zout) <- vlab
  zout
}

# panLab      (Place text labels in plot region)
# Purpose:    Place a text label in the plot region defined by (0,1), (0,1).
#             The ... notation allows all parameters available to "text" to be
#             passed.
# Parameters: x, y are the coordinates of the label
#             txt is the text
# Returns:    NULL (invisibly)
# Source:     A.R. Kronlund
# Revised:    K.Holt; 13-Jan-10 to accomodate axes on log scale
panLab <- function( x, y, txt, ... )
{
  # Allows text to be placed in plot panel at 0<x<1, 0<y<1.
  usr <- par( "usr" )
  
  yLog <- par("ylog")
  xLog <- par("xlog")
  
  # Check for log-transformed axes and adjust usr commands as needed
  # note: when a log scale is in use, 
  #           usr gives limits in the form 10 ^ par("usr")

  # Case 1: neither axis is on the log scale
  if (yLog==FALSE & xLog==FALSE)
  {
    par( usr=c(0,1,0,1) )
  }
  # Case 2: only the y-axis is on log scale
  if (yLog==TRUE & xLog==FALSE) 
  {
    usr[3:4]<-10 ^ par("usr")[3:4]
    par( usr=c(0,1,0,1), ylog=FALSE )
  } 
  # Case 3: only the x-axis is on log scale
  if (yLog==FALSE & yLog==TRUE) 
  {
    usr[1:2]<-10 ^ par("usr")[1:2]
    par( usr=c(0,1,0,1), xlog=FALSE )
  } 
  # Case 4: both axes are on the log scale
  if (yLog==TRUE & xLog==TRUE) 
  {
    usr[1:4]<-10 ^ par("usr")[1:4]
    par( usr=c(0,1,0,1), xlog=FALSE, ylog=FALSE )
  } 
  text( x, y, txt, ... )
  par( usr=usr )
  return( NULL )
}

#panLab <- function( x, y, txt, ... )
#{
#  # Allows text to be placed in plot panel at 0<x<1, 0<y<1.
#  usr <- par( "usr" )
#  par( usr=c(0,1,0,1) )
#  text( x, y, txt, ... )
#  par( usr=usr )
#  return( NULL )
#}

# panLegend   (Place legend in plot region)
# Purpose:    Place a legend in the plot region defined by (0,1), (0,1).
#             The ... notation allows all parameters available to "legend" to be
#             passed.
# Parameters: x, y are the coordinates of the legend
#             legTxt is the text associated with the legend
# Returns:    NULL (invisibly)
# Source:     A.R. Kronlund
# Revised:    K.Holt; 13-Jan-10 to accomodate axes on log scale
panLegend <- function( x, y, legTxt, ... )
{
  # Allows legend to be placed at 0<x<1, 0<y<1.
  usr <- par( "usr" )
  yLog<-par("ylog")
  xLog<-par("xlog")
  # Check for log-transformed axes and adjust usr commands as needed
    # note: when a log scale is in use, 
    #           usr gives limits in the form 10 ^ par("usr")
  # Case 1: neither axis is on the log scale
  if (yLog==FALSE & xLog==FALSE) {
    par( usr=c(0,1,0,1) )
  }
  # Case 2: only the y-axis is on log scale
  if (yLog==TRUE & xLog==FALSE) 
  {
     usr[3:4]<-10 ^ par("usr")[3:4]
     par( usr=c(0,1,0,1), ylog=FALSE )
   } 
  # Case 3: only the x-axis is on log scale
  if (yLog==FALSE & yLog==TRUE) 
  {
     usr[1:2]<-10 ^ par("usr")[1:2]
     par( usr=c(0,1,0,1), xlog=FALSE )
   } 
  # Case 4: both axes are on the log scale
  if (yLog==TRUE & xLog==TRUE) 
  {
     usr[1:4]<-10 ^ par("usr")[1:4]
     par( usr=c(0,1,0,1), xlog=FALSE, ylog=FALSE )
   } 
  legend( x, y, legend=legTxt, ... )
  par( usr=usr )
  return( NULL )
}

.addQuotes <- function( str ) 
{
  # Adds double-quotes to a string.
  return( paste("\"", str, "\"", sep = "") )
}

.convertSlashes <- function( expr, os = .Platform$OS.type ) 
{
  if ( os == "windows" ) 
    expr = gsub("/", "\\\\", expr)
  else expr = gsub("\\\\", "/", expr)
  return(expr)
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
    
    # Evaluate the parse string.
    eval( parse( text=listText ) )
  }
  result
}

# .excelTable (Creates and saves a dataframe to Microsoft Excel table)
# Purpose:    For a given connection and input data, create a dataframe and
#             save the dataframe as a worksheet in Excel.
#             If an Excel .xls file with the same name already exists then
#             delete the file and create a new file, else create a new file.
# Parameters: channel is a RODBC connection string with the .xls file name
#             dat       : a list, vector, or matrix containing worksheet data.
#             tablename : character variable containing the string to appear on
#                         Excel tabs for each worksheet.
#             colnam    : a vector of column names with length ncol(dat).
#             rownam    : a vector of row names with length nrow(dat).
# Returns:    NULL
# Source:     Modified from T.K. Deering (PopSim.r).
.excelTable <- function( channel, dat, tablename, colnam, rownam )
{
  dframe             <- as.data.frame( dat )
  names( dframe )    <- colnam
  rownames( dframe ) <- rownam
  sqlSave( channel, dframe, tablename=tablename )
  
  return()
}

.findFileName <- function( suffix ) 
{
  # Returns all file names with extension suffix.
  # Modified from PBSadmb .win.findTpl
  
  spat = gsub("\\.", "\\\\\\.", suffix)
  suff = list.files( pattern=paste( spat,"$",sep=""), ignore.case = TRUE )
  pref = substring(suff, 1, nchar(suff) - 4)
  return( pref )
}

# .closeActWin (close the active window)
# Purpose:     Closes the active window, say when the "exit" button is pressed.
# Parameters:  None
# Returns:     NULL (invisibly)
# Source:      A.R. Kronlund
.closeActWin <- function()
{
  closeWin( .getWinName() )
}

.getStamp <- function()
{
  stamp <- paste( format(Sys.time(),format="%d%m%Y%H%M%S" ),sep="" )
  return( stamp )
}

# .getWinName  (get the current winName)
# Purpose:     Determine which GUI is active (guiSim, guiView, guiPerf, etc.)
# Parameters:  None
# Returns:     A character containing the name of the current GUI window
# Source:      A.R. Kronlund, modified from PBSref (helper_funs.r)
.getWinName <- function()
{
  win <- .PBSmod$.activeWin
  
  # This is only required if PBSask is used, leave it for now.
  if(win == "PBSask")
  {
    win <- getWinVal("win", winName="PBSask")[[1]]   # Hidden field in PBSask
    win <- gsub("\n", "", win)                       # Remove the linefeed \n
  }
  return(win)
}

.intVal <- function( x )
{
  # Is the value of x an integer?  There must be a better way...
  result <- ifelse( (trunc(x)-x) == 0.0,TRUE,FALSE )
  result
}

.posVal <- function( x )
{
  # Sets all values of x < 0 to NA.
  x[ x < 0 ] <- NA
  x
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

# .unEvalList (converts possibly nested list to non-nested list)
# Purpose:    Converts a possibly nested list to a non-nested list.
# Parameters: obj is the possibly nested list to convert.
# Returns:    result, the non-nested list.
# Source:     A.R. Kronlund
.unEvalList <- function( obj ) 
{
  # Loop thru a (possibly) nested list, plucking out the name at the lowest
  # level and corresponding value.
  
  result   <- list()
  val      <- unlist( obj )
  valNames <-  names( val )
  
  for ( i in 1:length(val) )
  {
     tokenPos <- max(which(strsplit(valNames[i],'')[[1]]=='.')) + 1
     
     guiName <- substring( valNames[i], tokenPos,nchar(valNames[i]) )
     
     # Check for character or numeric.
     if ( is.character( val[i] ) )
       listText <- paste( "result$",guiName,"=\"",val[i],"\"",sep="" )
     else
       listText <- paste( "result$",guiName,"=",val[i],sep="" )
     eval( parse( text=listText ) )
  }
  result
}      

.updateGUI <- function()
{
   parentList <- ls( name=parent.frame(n=1) )
   
   win     <- .getWinName()                       # Get the current window name
   guiList <- getWinVal( scope="L", winName=win ) # GUI information local scope
  
   # Check for parent environment variables that match the GUI list.
   isMatch <- is.element( parentList,names(guiList) )
   parentList <- parentList[isMatch]
  
   # Now evaluate the variables into a list.
   nVals <- length( parentList )
   vals  <- as.list( 1:nVals )
   names( vals ) <- parentList
   
   for ( i in 1:length(vals) )
     vals[[i]] <- get( parentList[i], parent.frame(n=1) )
   
   setWinVal( vals )  
}
  
# .viewFile   (view a file saved in the mseRtemp directory)
# Purpose:    View a file that is stored in the mseR library directory in the
#             folder named "mseRtemp". This is the folder where copies of the R
#             code, the GUI description, the initial database, the ADMB
#             executable, and the documentation are kept.
# Parameters: fname is a character containing the name of the file to view
#             (default is based on the last action performed by the current
#             GUI window)
# Returns:    NULL
# Source:     PBSref (gui_funs.r")
.viewFile <- function(fname)
{
  # These two will be used when mseR is a proper R library.
  pckg  <- .PACKAGE                    # The name of this package
  dname <- paste( pckg,.FTEMP,sep="" ) # R directory where the file is located

  if( missing(fname) )
  {
    fname <- getWinAct(.getWinName())[1] # Name of the file to open
  }

  # This will be used when mseR is a proper R library.
  #rdir <- system.file(package = pckg)   # path to the R directory

  # Reference working directory.
  wkDir <- getwd()
  fname <- paste(wkDir, dname, fname, sep = "/")

  openFile(fname)

  return()
}

# .viewHelp   (view a help file or document)
# Purpose:    View a file that is stored in the mseR library directory in the
#             folder named "mseRtemp". This is the folder where copies of the R
#             code, the GUI description, the initial database, the ADMB
#             executable, and the documentation are kept.
# Parameters: fname is a character containing the name of the file to view
#             (default is based on the last action performed by the current
#             GUI window)
# Returns:    NULL
# Source:     PBSref (gui_funs.r")
.viewHelp <- function(fname)
{
  pckg  <- .PACKAGE                      # The name of this package
  dname <- paste( pckg,.FHELP,sep="" )   # R directory where the file is located

  if( missing(fname) )
  {
    fname <- getWinAct(.getWinName())[1] # Name of the file to open
  }

  # This will be used when mseR is a proper R library.
  #rdir <- system.file(package = pckg)   # path to the R directory

  # Reference working directory.
  wkDir <- getwd()                       # Path to R working directory
  fnam <- paste(wkDir, dname, fname, sep = "/")

  openFile(fnam)

  return()
}
