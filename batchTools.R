# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
# batchTools.R
#
# Tools for creating and running batch experiments with
# ms3R
# 
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# makeBatch()
# Takes a batch control file and produces all the necessary structure
# to run the batch of sim-est experiments.
# inputs:     batchCtl = character naming the batch control file
# ouputs:     batchDesign = data.frame containing the batch design
# usage:      to run batch jobs for multiple scenario/mp combinations
# author:     S.D.N. Johnson
makeBatch <- function ( batchCtlFile = "batchControlFile.bch", 
                        prjFld = ".",
                        batchFld = "Batch", 
                        baseCtlFile = "simCtlFileBase.txt")
{
  .subChar <<- "__"
  batchFilePath <- file.path(prjFld,batchFld,batchCtlFile)
  baseFilePath  <- file.path(prjFld,batchFld,baseCtlFile)
  # First, load the batch control file
  batchCtl <- .readParFile ( batchFilePath )
  # Now load the base control file
  baseCtl  <- .readParFile ( baseFilePath )

  # Set globals 
  #(THIS SHOULD BE MOVED TO ANOTHER FILE, OR DIFFERENT APPROACH FOUND)
  .PRJFLD <<- prjFld
  .DEFBATFLD <<- batchFld

  # Create Batch Design 
  batchDesign <- .createBatchDesign(  ctlPar = batchCtl, basePars = baseCtl )
  .FBATDES  <- file.path(getwd(),.PRJFLD,.DEFBATFLD,"batchDesign.txt")
  .CTLDES   <- file.path(getwd(),.PRJFLD,.DEFBATFLD,baseCtlFile)
  .BCTLDES  <- file.path(getwd(),.PRJFLD,.DEFBATFLD,batchCtlFile)
  .writeDesignFile (obj=batchDesign,out=.FBATDES)
  file.copy(baseCtlFile, .CTLDES)
  file.copy(batchCtlFile, .BCTLDES)

  return(batchDesign)
}

# Calls runSimEst in parallel inside a separate copy of the
# working directory, allowing for parallel system calls
doBatchRun <- function( arg )
{
  require(tools)
  cat("Running batchjob:", arg[1],"\n")
  # source control script to load DLL
  source("ms3R.R")
  
  # runMSE with the batch file
  # add random delay to offset simFolder names
  runMS3(ctlFile=arg[1], outFolder=arg[2])
  return(NULL)
}



# .runBatchJob  (Runs the simulations specified in a design dataframe)
# Purpose:      Loops through the rows of the design dataframe, each of which
#               specifies an input parameter file and labels for a simulation.
#               The mseR function runMSE is called for each row which generates
#               the simulation results (e.g., blob) and simulation folders.
# Parameters:   batchDesign is a batch design dataframe created by the function
#                 .createBatchDesign.
#               prefix is the Simulation File Prefix/
# Returns:      NULL (invisibly)
# Side Effects: A simulation folder containing *.info and *.Rdata file (blob) and
#               for each row of the design dataframe, i.e., for each simulation.
# Source:       A.R. Kronlund
.runBatchJob <- function( batchDesign=NULL, 
                          par=FALSE, 
                          prefix=NULL, 
                          initPar = 1,
                          subset = NULL, 
                          nCores = detectCores()-1)
{
  # Runs simulations from the design data.frame specified in batchDesign object.
  # 1. Does the mseR input parameter file exist? If YES then read the file.
  # 2. Run the simulation using runMSE().
  # 3. Update the design data.frame somehow...
  if (is.null(batchDesign))
  {
    browser()

    desPath <- file.path(getwd(),.PRJFLD,.DEFBATFLD,"batchDesign.txt")
    batchDesign <- read.csv(desPath, header=TRUE, skip=1, stringsAsFactors=FALSE)
  } 

  # Vector of simCtlFile names, like simCtlFile1.txt, simCtlFile2.txt, etc.
  batchParFile <- file.path( getwd(),.PRJFLD, .DEFBATFLD, basename( batchDesign$parFile ) )
  blobName     <- batchDesign$blobName  # Vector of blobNames.
  nJobs        <- length(blobName)      # Number of blobs (i.e., simulations).

  result <- data.frame( sim=character(nJobs), scenarioLabel=character(nJobs),
                        mpLabel=character(nJobs), elapsed=numeric(nJobs),
                        stringsAsFactors=FALSE )

  # This is part where shell script PPSS could be used to parallelize the batch job
  # Need to create a folder of inputParameters.par files to be processed via runMSE().
  # Might need to have the .par file be an argument to runMSE()...check what effects
  # that would have elsewhere.

  # if a parallel flag is set, run in parallel
  if( par )
  { 
    # First, load parallel package
    library(parallel)
    # turn off squawking
    options(warn=-1)
    # Get number of batch runs
    nBatchFiles   <- length(batchParFile)

    # Create folder names for batch running
    batchFolderNames <- paste("parBat",prefix,1:nBatchFiles,sep="")
    
    # combine folder and control file names
    if( !is.null( subset ) )
    {
      parBatchArgList <- vector(mode = "list", length = length(subset) )
      for( i in 1:length(subset) )
      {
        parBatchArgList[[i]] <- c( batchParFile[ subset[i] ],batchFolderNames[ subset[i] ] )   
      }
    } else {
      parBatchArgList <- vector(mode = "list", length = length(batchParFile) - initPar + 1)
      for(i in initPar:nBatchFiles)
      {
        parBatchArgList[[i - initPar + 1]] <- c( batchParFile[i],batchFolderNames[i] )
      }
    }

    nJobs <- length(parBatchArgList)

    # Now set # of cores and make a cluster
    nCores  <- min(length(parBatchArgList),nCores)
    cl      <- makeCluster(nCores, outFile = "./parBatchMsg.txt")
    # Run parallel batch
    cat ("Fitting ", nJobs, " scenario/hypothesis combinations in parallel on ",
          nCores, " cores.\n", sep = "" )
    tBegin    <- proc.time()
    startDate <- date()
    tmp       <- clusterApplyLB(cl, x=parBatchArgList, fun=doBatchRun)
    # tmp <-lapply(X=parBatchArgList, FUN=doBatchRun)
    stopCluster(cl)

    # Now copy the contents of each batchFolderName to the project folder
    elapsed <- (proc.time() - tBegin)[ "elapsed" ]
    cat( "\nMSG (.runBatchJob): Elapsed time for parallel batch = ",
      round(elapsed/60.0,digits=2)," minutes.\n" )

  } else for ( i in initPar:nJobs ) {
    
    if ( file.exists( batchParFile[i] ) )
    {
      fileName <- strsplit( batchParFile[i],"\\." )[[1]][1]
      cat( "\nRunning batch job: ",batchParFile[i],"...\n" )

      tBegin    <- proc.time()
      startDate <- date()

      cat( "\nMSG (.runBatchJob) Processing batchParFile = ", batchParFile[i], "\n" )
      
      folderName <- paste(prefix,"bat",i,sep ="")
     
      runMS3( batchParFile[i], outFolder = folderName )

      elapsed <- (proc.time() - tBegin)[ "elapsed" ]
      cat( "\nMSG (.runBatchJob): Elapsed time for simulation = ",
        round(elapsed/60.0,digits=2)," minutes.\n" )
    }
    
    result[ i,"sim" ]           <- batchDesign[ i,"parFile" ]
    result[ i,"scenarioLabel" ] <- batchDesign[ i,"scenarioLabel" ]
    result[ i,"mpLabel" ]       <- batchDesign[ i,"mpLabel" ]
    result[ i,"elapsed" ]       <- round( elapsed/60.0, digits=2)
    
    cat( "\nMSG (.runBatchJob) Progress as of ", date()," Fit ",i," of ",nJobs,"\n\n" )
    print( result[1:i,] )
  }

  # Play sound to signal that work is complete
  beep(sound = "complete")

  cat( "\n Work Complete! \n" )
  invisible()
}     # END function .runBatchJob


# .writeDesignFile (Writes the batch job design dataframe to a file)
# Purpose:         Given a dataframe, assumed to contain a batch job design,
#                  created by .createBatchDesign, write to text file.
# Parameters:      obj is a dataframe.
#                  outFile is the desired filename.
# Returns:         NULL (invisibly)
# Side Effects:    Call to this function may generate a warning, not sure why.
# Source:          A.R. Kronlund
.writeDesignFile <- function( obj, outFile=.FBATDES )
{
  nRow <- nrow( obj )
  if ( is.null(nRow) )
    nRow <- 1

  cat( file=outFile, "# mseR batch design file, ",date(),"\n" )
  
  # This generates a warning, not sure why.
  write.table( obj, file=outFile,
               col.names=TRUE, row.names=FALSE, quote=TRUE, sep=",", append=TRUE )
  cat( "\nMSG (.writeDesignFile): Batch design file written to ",outFile,"\n" )
  invisible()
}     # END function .writeDesignFile


# .createBatchDesign  (Creates design dataframe and writes input par files)
# Purpose:      Given a batch control list, base parameter file, and simulation
#               file prefix, this function creates the batch job design by
#               crossing each management procedure node with each scenario node.
#               If the number of scenarios is nScenario and the number of
#               management procedures is nMP, then the result is a design
#               dataframe with nScenario*nMP rows and a *.par file corresponding
#               to each row (combination of scenario and management procedure).
#               The output design dataframe is passed to the .runBatchJob to
#               control completing the simulations using runMSE.
# Parameters:   ctlList is a list created by .createKeyList from the Batch File
#               basePars is a dataframe containing the base parameters that are
#                 modified by the ctlList values;
#               prefix is the Simulation File Prefix string concatenated to the
#               design filename, output parameter files, and simulation results
#               in .Rdata files (e.g., blobs).
# Returns:      result is the batch job design dataframe.
# Side Effects: The *.par files that specify each simulation are written.
# Source:       A.R. Kronlund, and probably under questionable conditions.
.createBatchDesign <- function( ctlPar, basePars, prefix="job" )
{
  # Input a control list object, usually the output from .createList.
  ctlList <- .createKeyList(ctlPar)
  
  scenario  <- ctlList$scenario              # List of scenarios.
  nScenario <- length( scenario )            # Number of scenarios.

  mp  <- ctlList$mp                          # List of management procedures.
  nMP <- length( mp )                        # Number of management procedures.

  nBatch <- nScenario * nMP                  # Number of mseR simulations.

  # Design dataframe - each row identifies the simulation and adds colors,
  # line types, line widths and symbols to be used for plotting outside of mseR.
  # This file is needed as input to runBatchJob to control the job, and for use
  # in plotting and performance calculations are to be done outside of mseR.
  
  parFile      <- character( nBatch )        # Name of each input parameter file.
  blobName     <- character( nBatch )        # Name of each blob.

  scenarioLabel <- character( nBatch )       # Name of each scenario.
  mpLabel       <- character( nBatch )       # Name of each management procedure.
  
  dataName     <- rep( "data",     nBatch )  # Name of data method.
  methodName   <- rep( "method",   nBatch )  # Name of assessment method.
  ruleName     <- rep( "rule",     nBatch )  # Name of HCR.
  scCol        <- rep( "black",    nBatch )  # Color for scenario.
  scLty        <- rep( 1,          nBatch )  # Line type for scenario.
  scLwd        <- rep( 1,          nBatch )  # Line width for scenario.
  scSym        <- rep( 1,          nBatch )  # Symbol for scenario.
  mpCol        <- rep( "black",    nBatch )  # Color for procedure.
  mpLty        <- rep( 1,          nBatch )  # Line type procedure.
  mpLwd        <- rep( 1,          nBatch )  # Line width for procedure.
  mpSym        <- rep( 1,          nBatch )  # Symbol for procedure.
  join         <- rep( 0,          nBatch )  # Join... ???

  iBatch <- 1

  # Loop over the scenarios.
  for ( i in 1:nScenario )
  {
    # Loop over the management procedures.
    for ( j in 1:nMP )
    {
      parFile[iBatch]   <- paste( prefix,iBatch,".txt",sep="" )
      
      # Create a unique blobname.
      # NOTE: This step is deferred until the simulation is launched to obtain
      #       a unique date-time stamp at time of execution.  For now simply
      #       provide a name indexed to the simulation number.
      blobName[iBatch] <- paste( "sim",prefix,iBatch,sep="" )
      
      # Set default scenarioLabel and mpLabel values in case not supplied.
      scenarioLabel[iBatch] <- paste( "scenario",i,sep="" )
      mpLabel[iBatch]       <- paste( "mp",j,sep="" )

      # Use the scenarioLabel if provided.
      #if ( !is.null(scenario[[i]]$scenarioLabel) )
      #  scenarioLabel[iBatch] <- scenario[[i]]$scenarioLabel

      # Scenario colour.
      if ( !is.null(scenario[[i]]$scCol ) )
        scCol[iBatch] <- scenario[[i]]$scCol

      # Scenario line type.
      if ( !is.null(scenario[[i]]$scLty ) )
        scLty[iBatch] <- scenario[[i]]$scLty

      # Scenario line width.
      if ( !is.null(scenario[[i]]$scLwd ) )
        scLwd[iBatch] <- scenario[[i]]$scLwd

      # Scenario symbol.
      if ( !is.null(scenario[[i]]$scSym ) )
        scSym[iBatch] <- scenario[[i]]$scSym

      # Use the mpLabel if provided.
      #if ( !is.null(mp[[j]]$mpLabel) )
      #  mpLabel[iBatch] <- mp[[j]]$mpLabel

      # Use the dataName if provided.
      if ( !is.null(mp[[j]]$dataName) )
        dataName[iBatch] <- mp[[j]]$dataName

      # Use the methodName if provided.
      if ( !is.null(mp[[j]]$methodName ) )
        methodName[iBatch] <- mp[[j]]$methodName

      # Use the ruleName if provided.
      if ( !is.null(mp[[j]]$ruleName) )
        ruleName[iBatch] <- mp[[j]]$ruleName

      # Management procedure colour.
      if ( !is.null(mp[[j]]$mpCol ) )
        mpCol[iBatch] <- mp[[j]]$mpCol

      # Management procedure line type.
      if ( !is.null(mp[[j]]$mpLty ) )
        mpLty[iBatch] <- mp[[j]]$mpLty

      # Management procedure line width.
      if ( !is.null(mp[[j]]$mpLwd ) )
        mpLwd[iBatch] <- mp[[j]]$mpLwd

      # Management procedure symbol.
      if ( !is.null(mp[[j]]$mpSym ) )
        mpSym[iBatch] <- mp[[j]]$mpSym

      # Join the procedures (groups share common integer value).
      if ( !is.null(mp[[j]]$join ) )
        join[iBatch] <- mp[[j]]$join

      # Replace values for any keywords that match those in the mseR input
      # parameters file:
      # 1. Find keywords shared between the batch job control list and the mseR
      #    parameter file.
      # 2. Replace the values in the mseR parameter file with those from the
      #    batch control list.
      
      newPars <- basePars
      
      # This step compares the names in the i-th scenario to the names in the
      # first column of newPars to find the common names using intersect.
      # browser()

      # Change all .subChar symbols to "$".
      names(scenario[[i]]) <- gsub( .subChar,"$",names(scenario[[i]]), fixed=TRUE )
       
      sharedNames <- intersect( names(scenario[[i]]),newPars[,1] )

      # Create a batch root, and copy new entries from the batchPar list
      scenarioRoot <- paste("scenario$scenario",i,"$",sep="")
      sharedNamesPar <- paste(scenarioRoot,sharedNames,sep="")
 
      if ( length(sharedNames) > 0 )
        for ( k in 1:length(sharedNames) )
        {
          val <- ctlPar[,2][ ctlPar[,1] == sharedNamesPar[k] ]
          # if ( is.character(val) )
          #   val <- paste("\"",val,"\"",sep="")
          newPars[,2][ newPars[,1]==sharedNames[k] ] <- val
        }

      # Change all .subChar symbols to "$".
      names( mp[[j]] ) <- gsub( .subChar,"$",names(mp[[j]]), fixed=TRUE )
       
      sharedNames <- intersect( names(mp[[j]]),newPars[,1] )
      # Create an MP root for copying from batchPars table
      mpRoot <- paste ("mp$mp",j,"$",sep="")
      sharedNamesPar <- paste(mpRoot,sharedNames,sep="")
      
      if ( length(sharedNames) > 0 )
        for ( k in 1:length(sharedNames) )
        {
          val <- ctlPar[,2][ ctlPar[,1] == sharedNamesPar[k] ]
          # if ( is.character(val) )
          #   val <- dQuote( val )

          newPars[,2][ newPars[,1]==sharedNames[k] ] <- val
        }

      # Check to see if scenarioLabel and mpLabel are updated.
      scenarioLabel[iBatch] <- newPars[,2][ newPars[,1]=="ctl$scenarioName" ]
      mpLabel[iBatch] <- newPars[,2][ newPars[,1]=="ctl$mpName" ]
      
      # Remove leading white space: gsub('^[[:space:]]+', '', " a test ")
      newPars[,2] <- gsub( '^[[:space:]]+','',newPars[,2] )
      
      fName <- parFile[iBatch]          # mseR input parameters file.
      
      fName <- file.path( .PRJFLD, .DEFBATFLD, fName )
      batchDate <- date()               # Current date and time.      
      
      # Open a new parameter file and write the title and date/time stamp.
      cat( file=fName,
        "# ",parFile[iBatch],": mseR parameter file written ",batchDate,".\n", sep="" )

      # NOTE: write.table wants a matrix or data.frame, not a vector.
 
      # Write the header field names for the parameter file.
      
      colNames <- names( basePars )
      #write.table( file=fName, matrix( colNames, nrow=1 ), quote=FALSE,
      #             col.names=FALSE, row.names=FALSE,
      #             sep=" ", append=TRUE )
                   
      write.table( file=fName, newPars, quote=FALSE,
                   col.names=TRUE, row.names=FALSE,
                   sep=" ", append=TRUE )          #RF changed append=FALSE to append=TRUE
   
      #options( warn=-1 )                        # Turn off whining.
      #for ( k in 1:nrow(newPars) )
      #{
      #  isNumericVal <- !is.na( as.numeric( newPars[k,2] ) )  # Coerce non-numeric to NA.
      #  if ( isNumericVal )
      #    cat( file=fName, newPars[k,1]," ",newPars[k,2],"\n", append=TRUE, sep="" )
      #  else
      #    cat( file=fName, newPars[k,1]," ",dQuote(newPars[k,2]),"\n", append=TRUE, sep="" )        
      #}
      #options( warn=0 )                 # Turn on whining.
      
      cat( "\nMSG(.createBatchDesign) mseR parameter file ",fName," written...\n" )

      iBatch <- iBatch + 1           # Increment the batch job counter.
    }
  }
  
  # Bind the vectors that make up the design dataframe.
  result <- data.frame( parFile, blobName, scenarioLabel,
              mpLabel, dataName, methodName, ruleName,
              scCol, scLty, scLwd, scSym,
              mpCol, mpLty, mpLwd, mpSym, join,
              stringsAsFactors=FALSE )
  result
}    # END function .createBatchDesign


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

# .createKeyList (Creates a list from the Batch File for .createBatchDesign):
# Purpose:      Function to convert the Batch File dataframe (loaded by the
#               function .readParFile) into a list that structures parameters
#               into "scenario" and "mp" nodes, each of each contains the
#               parameters for unique scenarios and management procedures,
#               respectively.
# Parameters:   obj is a dataframe read by .readParFile with columns "parameter"
#                 and "value".
# Returns:      result, a list with the batch control structure required to form
#               the desired cross-combinations of scenarios and procedures.
# Source:       A.R. Kronlund
.createKeyList <- function( obj )
{
  # Input  a data frame with columns "parameter" and "value".
  # Output a list with elements named as first level of parameter and the
  # balance as the key, with values in "value".

  result <- list()

  options( warn=-1 )                        # Turn off whining.
  numericVal <- as.numeric( obj[,"value"] ) # Coerce non-numeric to NA.
  options( warn=0 )                         # Turn on whining.

  for ( i in 1:nrow(obj) )
  {
    # Value is numeric, build the parse string.
    if ( !is.na(numericVal[i]) )
    {
      parName <- obj[i,"parameter"]
      # Replace all the "$" with "&" in parameter name.
      parName <- gsub( "$",.subChar,parName, fixed=TRUE )
      # Replace ONLY the first "&" with "$" in parameter name, TWICE.
      parName <- sub( .subChar,"$",parName, fixed=TRUE )
      parName <- sub( .subChar,"$",parName, fixed=TRUE )
      
      listText <- paste( "result$",parName,"=",
                    obj[i,"value"],sep="" )
    }
    # Value is character, build the parse string.
    else
    {
      parName <- obj[i,"parameter"]
      # Replace all the "$" with "&" in parameter name.
      parName <- gsub( "$",.subChar,parName, fixed=TRUE )
      # Replace ONLY the first "&" with "$" in parameter name.
      parName <- sub( .subChar,"$",parName, fixed=TRUE )
      parName <- sub( .subChar,"$",parName, fixed=TRUE )
      
      listText <- paste( "result$",parName,"=",
                    obj[i,"value"],sep="" )
    }

    # ARK: At one point I had this code, presumably to handle a different
    #      input format convention, perhaps assuming "value" was all character.
    #                   sQuote(obj[i,"value"]),sep="" )
    # Evaluate the parse string.
    eval( parse( text=listText ) )
  }
  result
}    # END function .createKeyList