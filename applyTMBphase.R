# <><><><><><><><><><><><><><><><><><><><><><><><><><><>
# applyTMBphase.R
#
# Function for applying TMB in phases
#
# Author: SDN Johnson
# Date: Nov 13, 2019
#
# Last Update: Nov 13, 2019
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><>


# .applyTMBphase()
# Wrapper function that
# takes control lists and runs a TMB
# fit procedure with parameter phasing.
# inputs:   tmbLists
#           dllName
#           ...
# outputs:  amObj = list of fit info and perf
# usage:    within a model based MP
# Source:   SDN Johnson, modfied from TMB wiki
.applyTMBphase <- function( tmbLists,
                            dllName = "hierProd",
                            optimizer = "nlminb",
                            silent = FALSE,
                            calcSD = FALSE,
                            maxPhase = NULL,
                            base_map = list(),
                            maxEval = 1e3,
                            maxIter = 1e3,
                            savePhases = TRUE,
                            phaseMsg = TRUE )
{
  # Pull all th lists
  data    <- tmbLists$data
  pars    <- tmbLists$pars
  phases  <- tmbLists$phases
  random  <- tmbLists$random

  # function to fill list component with a factor
  # of NAs
  fill_vals <- function(x,vals)
  { 
    factor( rep( vals, length(x) ) ) 
  }

  #loop over phases
  if(!is.null(maxPhase))
    maxPhase <- min( maxPhase, max(unlist(phases) ) )
  else maxPhase <- max(unlist(phases))

  # Make a list of reports for each phase
  phaseReports <- vector(mode = "list", length = maxPhase)

  # generate a list of outputs to return
  # to runHierSCAL, initialise 
  # a success flag at TRUE
  outList <- list( success = TRUE )

  for( phase_cur in 1:maxPhase ) 
  {
    gc()
    # Start timing
    tBegin <- proc.time()[3]

    # work out the map for this phase
    # if the phase for a parameter is greater than the current phase 
    # or a negative value, then map will contain a factor filled with NAs
    map_use <- base_map
    j <- length(map_use)
    for( i in 1:length(pars) ) 
    {
      parName <- names(pars)[i]

      if( parName %in% names(phases) )
      {
        if( (phases[[parName]] > phase_cur) | phases[[parName]] < 0 ) 
        { 
          # Check if parName is included in the base_map
          if(parName %in% names(map_use))
            map_use[[parName]] <- fill_vals(pars[[i]],NA)
          else
          {
            j <- j + 1
            map_use[[j]] <- fill_vals(pars[[i]],NA)
            names(map_use)[j] <- parName
          }

        }
      } else {
        j <- j + 1
        map_use[[j]] <- fill_vals(pars[[i]],NA)
        names(map_use)[j] <- parName
      }

    }

    #remove the random effects if they are not estimated
    random_use <- random[ !random %in% names(map_use) ]
  
    # initialize the pars at values in previous phase
    params_use <- pars
    if( phase_cur > 1 ) 
      params_use <- obj$env$parList( opt$par )


    mapNames <- names(map_use)
    parNames <- names(params_use)

    # Check names in map correspond to par names
    if( any( ! names(map_use) %in% names(params_use) ) )
    {
      badNames <- names(map_use)[ !names(map_use) %in% names(params_use)]
      cat( badNames )
      browser()

    }

    # if(phase_cur == maxPhase)
      # browser()

    # Fit the model using fixed effects only
    obj <- TMB::MakeADFun(  data = data,
                            parameters = params_use,
                            random = NULL,
                            DLL= dllName,
                            map= map_use,
                            silent = silent )

    # Create a control list for the assessment model
    tmbCtrl <- list(  eval.max = maxEval, 
                      iter.max = maxIter  )

    repInit <- obj$report()

    startPar <- obj$par

    if(phaseMsg)
      cat("\nStarting optimisation for phase ", phase_cur, "\n\n")

    # Try the optimisation
    opt <- try( nlminb (  start     = obj$par,
                          objective = obj$fn,
                          gradient  = obj$gr,
                          control   = tmbCtrl ) )


    # break if there is an issue
    if( class(opt) == "try-error" )
    {

      if(phaseMsg)
        cat("\nOptimisation halted due to error\n")

      outList$success                   <- FALSE
      outList$maxPhaseComplete          <- phase_cur - 1

      if( savePhases )
      {
        phaseReports[[phase_cur]]$opt     <- opt
        phaseReports[[phase_cur]]$success <- FALSE
      }

      browser()

      break
    }

    if(phase_cur == maxPhase & !grepl("relative", opt$message) )
    {
      message("Non-convergent optimisation, attempting to refit.\n")
      for( kRefit in 1:4 )
      {
        lastPar <- opt$par
        if( grepl("false", opt$message) )
        {
          lastPar <- lastPar + rnorm(length(lastPar), mean = 0, sd = 0.2)
        }


        opt <- try( nlminb (  start     = lastPar,
                              objective = obj$fn,
                              gradient  = obj$gr,
                              control   = tmbCtrl ) )


        if( grepl("relative",opt$message))
          break
      }
    }

    if(phase_cur == maxPhase & !is.null(random) )
    {
      if(phaseMsg)
        message("Running optimisation with random effects.\n")
      # Run fit process with random effects
      # First, remove the random effects from
      # the map if they are in there
      for( randIdx in 1:length(random))
      {
        parName <- random[randIdx]
        if( parName %in% names(map_use))
          map_use[[parName]] <- NULL
      }

      params_use <- obj$env$parList( opt$par )

      # Then make the ADFun object with new stuff
      # Fit the model
      objRE <- TMB::MakeADFun(  data = data,
                                parameters = params_use,
                                random = random,
                                DLL= dllName,
                                map= map_use,
                                silent = TRUE )

      objFE <- TMB::MakeADFun(  data = data,
                                parameters = params_use,
                                random = NULL,
                                DLL= dllName,
                                map= map_use,
                                silent = TRUE )

      if( is.nan(objRE$fn() ) )
        randObj <- objFE
      else randObj <- objRE

      # Create a control list for the assessment model
      tmbCtrl <- list(  eval.max = maxEval, 
                        iter.max = maxIter  )

      startPar <- randObj$par


      optRE <- try( nlminb (  start     = startPar,
                              objective = randObj$fn,
                              gradient  = randObj$gr,
                              control   = tmbCtrl ) )

      if( class(optRE) == "try-error" )
      {
        browser()
        message("NaN evaluation of gradient function.\n")
        for( kRefit in 1:4 )
        {
          newPar <- startPar + rnorm(length(startPar), mean = 0, sd = 0.3)

          optRE <- try( nlminb (  start     = newPar,
                                  objective = randObj$fn,
                                  gradient  = randObj$gr,
                                  control   = tmbCtrl ) )

          if( class(optRE) != "try-error")
            if( grepl("relative",opt$message))
            {
              obj <- randObj
              opt <- optRE
              break
            }

        }
      } else {
        obj <- randObj
        opt <- optRE
      }
    }

    # Save sdreport object
    if(phase_cur == maxPhase & calcSD )
    {
      sdrep     <- try(TMB::sdreport(obj))
      if( class(sdrep) != "try-error")
      {
        pdHess    <- sdrep$pdHess
        SDs       <- sqrt(diag(sdrep$cov.fixed))
        posSDs    <- all(is.finite(SDs))

        if( !pdHess & !posSDs )
        {

          message("AM had non PD Hessian, attempting to refit. \n")
          for( kRefit in 1:4 )
          {
            lastPar <- opt$par +  rnorm(length(opt$par), mean = 0, sd = 2)
          

            opt <- try( nlminb (  start     = lastPar,
                                  objective = obj$fn,
                                  gradient  = obj$gr,
                                  control   = tmbCtrl ) )

            sdrepRetry   <- try(TMB::sdreport(obj))
            
            if( class(sdrepRetry) != "try-error")
            {
              pdHess  <- sdrepRetry$pdHess 
              SDs     <- sqrt(diag(sdrepRetry$cov.fixed))
              posSDs  <- all(is.finite(SDs))
              sdRep   <- sdrepRetry
            }

            if(pdHess | posSDs)
              break
          }

        }

        browser()

        gradTable <- data.frame(  par = names(sdrep$par.fixed),
                                  est = sdrep$par.fixed,
                                  sd  = sqrt(diag(sdrep$cov.fixed)),
                                  grad = as.numeric(sdrep$gradient.fixed) )

        outList$gradReport  <- gradTable
        outList$sdrep       <- sdrep
        outList$pdHess      <- pdHess
        outList$posSDs      <- posSDs
      }

      if( class( sdrep ) == "try-error" )
      {
        message("sdrep failed\n")
        browser()
      }
    }


    # Save max phase complete
    outList$maxPhaseComplete          <- phase_cur

    # Save reports and optimisation
    # output
    if(savePhases)
    {
      phaseReports[[phase_cur]]$pars    <- obj$env$parList(opt$par)
      phaseReports[[phase_cur]]$opt     <- opt
      phaseReports[[phase_cur]]$success <- TRUE
      phaseReports[[phase_cur]]$map     <- map_use
    }

    # browser()

    if(phaseMsg)
      cat(  "\nPhase ", phase_cur, " completed with code ",
          opt$convergence, " and following message:\n", sep = "" )
    if(phaseMsg)
      cat("\n", opt$message, "\n\n", sep = "" )
    
  } # close phase loop

  # if(!outList$pdHess)
  #   browser()



  # Save phase reports
  if( savePhases )
    outList$phaseReports      <- phaseReports

  # Now save report object
  if(outList$success)
  {
    outList$repOpt            <- obj$report()
    outList$optPar            <- obj$env$parList( opt$par )
  } else {
    outList$optPar         <- params_use
    outList$repOpt         <- repInit
  }


  # And the remainder of the details
  outList$objfun            <- obj$fn()
  outList$optOutput         <- opt
  outList$map               <- map_use
  if(!is.null(outList$sdrep))
    outList$maxGrad           <- max(abs(outList$sdrep$gradient.fixed),na.rm =T)
  else
    outList$maxGrad           <- max(abs(obj$gr(opt$par)),na.rm =T)
  # outList$totTime           <- sum(fitReport$time,na.rm = TRUE)
  
  return( outList )  

} # END .applyTMBphase

