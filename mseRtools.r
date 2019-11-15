#-----------------------------------------------------------------------------##
#-- Helper Functions (some HIDDEN, e.g., .foo)                              --##
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
