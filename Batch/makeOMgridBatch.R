age1M   <- c("mean","1.25","1.64")
srSteep <- c(0.67, 0.5, 0.8)
areaM   <- c("identMDevs","diffMDevs")
rw      <- c(15,5)
PondM   <- c("lo","hi")

histFolderNames <- paste("fit_parBatHGherringGrid",1:16,sep = "")

outFile <- "HGherringGrid.bch"

gridList <- list( #age1M = age1M, srSteep = srSteep, areaM = areaM, 
                  histFolder = histFolderNames,
                  rw = rw, PondM = PondM )
omGrid <- expand.grid(gridList)
omGrid$histFolder <- as.character(omGrid$histFolder)
# omGrid$age1M <- as.character(omGrid$age1M)
# omGrid$areaM <- as.character(omGrid$areaM)


cat(  "# Batch Control File, created ", date(), " by makeBatchDesign() \n", 
        file = outFile, append = F, sep = "" )
cat( "parameter value\n", sep = "", append = T, file = outFile)
cat( "#\n", file = outFile, append = T )
cat( "# OM Scenarios \n", file = outFile, append = T )
cat( "#\n", file = outFile, append = T )

for( k in 1:nrow(omGrid))
{
  # age1Mlevel  <- omGrid$age1M[k]
  # steepVal    <- omGrid$srSteep[k]
  # areaMlevel  <- omGrid$areaM[k]
  histOMLabel <- omGrid$histFolder[k]

  # Pull the OM model hyp label
  omHypLabel <- read.csv(here::here("history",histOMLabel,"parEstTable.csv"))$modelHyp[1]

  # Split
  omHypLabs <- stringr::str_split(omHypLabel,"_")

  age1Mlevel  <- stringr::str_split(omHypLabs[[1]][1],"juveM")[[1]][2]
  steepVal    <- stringr::str_split(omHypLabs[[1]][2],"steep")[[1]][2]
  areaMlevel  <- omHypLabs[[1]][3]

  rwLevel     <- omGrid$rw[k]
  pondMlevel  <- omGrid$PondM[k]

  if(pondMlevel == "lo")
    pondMmu_f <- "c(0,0,0,0,0,0.315,0)"
  if(pondMlevel == "hi")
    pondMmu_f <- "c(0,0,0,0,0,1.05,0)"

  if(areaMlevel == "diffMDevs")
  {
    areaMlab <- "diffM"
    projMtype <- "correlated"
  }

  if(areaMlevel == "identMDevs")
  {
    areaMlab <- "identM"
    projMtype <- "identical"
  }

  scenLabel <- paste("a1M",age1Mlevel,"_sr",steepVal,"_",areaMlab,"_rw", rwLevel, "_",pondMlevel,"PondM",sep = "")

  cat( "# Scenario ", k, " :", scenLabel ,"\n", file = outFile, append = T, sep = "" )
  cat( "#\n", file = outFile, append = T )
  cat(  "scenario$scenario", k, "$ctl$scenarioName '", scenLabel, "'\n", 
        sep = "", append = T, file = outFile )
  cat(  "scenario$scenario", k, "$opMod$histFile '", histOMLabel, "'\n", 
        sep = "", append = T, file = outFile )
  cat(  "scenario$scenario", k, "$opMod$projMtrend ", rwLevel, "\n", 
        sep = "", append = T, file = outFile )
  cat(  "scenario$scenario", k, "$opMod$projMtype '", projMtype, "'\n", 
        sep = "", append = T, file = outFile )
  cat(  "scenario$scenario", k, "$opMod$pondMmu_f ", pondMmu_f, "\n", 
        sep = "", append = T, file = outFile )
  cat( "#\n", file = outFile, append = T )
}

cat( "#\n", file = outFile, append = T )
cat( "#\n", file = outFile, append = T )

cat( "# Model Hypotheses \n", file = outFile, append = T )
cat( "#\n", file = outFile, append = T )


cat( "#\n", file = outFile, append = T )
cat( "#\n", file = outFile, append = T )
cat( "# File Ends <not run>.\n", file = outFile, append = T )


