// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>
// hierProd.cpp
// 
// A multi-species, multi-stock surplus production (Schaefer) state-space 
// stock assessment model with  joint prior distributions applied to 
// catchability (q) and productivity (U_msy), and a process error 
// year-effect component (eps_t) shared between stocks/species
// 
// Author: Samuel Johnson
// Date: 25 October, 2019
//
// Purpose: Assessment model for ms3R closed loop simulation package.
// 
// 
// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>

#include <TMB.hpp>       // Links in the TMB libraries
#include <iostream>

// posfun
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
}

// invLogit
template<class Type>
Type invLogit(Type x, Type scale, Type trans){
  return scale/(Type(1.0) + exp(-Type(1.0)*x)) - trans;
}

// invLogit
template<class Type>
Type square(Type x)
{
  return pow(x,2);
}

// dinvgamma()
// R-style inverse gamma probability distribution density
// function.
// inputs:    x = real number to calculate density of
//            alpha = shape parameter
//            beta = scale parameter
//            logscale = int determining whether function returns log (1) or 
//                        natural (0) scale density
// outputs:   dens = density on natural or log scale
// Usage:     usually for variance priors
// Source:    Wikipedia,  by S. D. N. Johnson
template<class Type>
Type dinvgamma( Type x,
                Type alpha,
                Type beta,
                int logscale )
{

  Type dens = -1 * ( (alpha + 1) * log(x) + beta / x );

  if( logscale == 0 )
    dens = exp( dens );

  return(dens);
}



// objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Call namespaces //
  using namespace density;

  /* ========== data section ========== */
  // Data Structures
  DATA_ARRAY(I_spft);           // Indices (CPUE/biomass)
  DATA_ARRAY(C_spt);            // Catch data
  
  // Model dimensions
  int nS = I_spft.dim(0);           // No. of species
  int nP = I_spft.dim(1);           // No. of species
  int nF = I_spft.dim(2);           // No. of surveys
  int nT = I_spft.dim(3);           // No of time steps

  // Model switches
  DATA_INTEGER(SigmaPriorCode); // 0 => IG on diagonal element, 1 => IW on cov matrix
  DATA_IVECTOR(condMLEq_f);     // 0 => q estimated, 1 => q derived from obs model (mean resid)
  DATA_IVECTOR(shrinkq_f);      // 0 => q a leading par, 1 => q shrank to complex mean
  DATA_IVECTOR(tvq_f);          // 0 => q constant, 1 => q time-varying
  DATA_IVECTOR(fleetq_f);       // 0 => no q estimated, 1 => q estimated for fleet f
  DATA_VECTOR(fleetWeights_f);  // scalars for fleet components of the likelihood function
  DATA_IVECTOR(condMLEobsErr_f);// 0 => obs error sd estimated as a free par, 1 => conditional MLE
  DATA_INTEGER(lnqPriorCode);   // 0 => hyperprior, 1 => multilevel
  DATA_INTEGER(lnUPriorCode);   // 0 => hyperprior, 1 => multilevel 
  DATA_INTEGER(BPriorCode);     // 0 => normal, 1 => Jeffreys 
  DATA_IARRAY(initT_sp);        // first year of assessment
  DATA_IARRAY(initProcErr_sp);  // first year of process error devs
  DATA_IARRAY(initBioCode_sp);  // initial biomass at 0 => unfished, 1=> fished
  DATA_IARRAY(calcIndex_spf);   // Indicator for calculating obs model variance for an index
  DATA_IARRAY(stockq_spf);      // Indicator for calculating catchability for a stock/species/fleet combo
  DATA_IARRAY(speciesq_sf);     // Indicator for calculating catchability for a species/fleet combo
  DATA_SCALAR(posPenFactor);    // Positive-penalty multiplication factor
  DATA_IARRAY(solveInitBio_sp); // Indicator for solving initial biomass from first data point and catchability. 0 => off, fIdx > 0 => fleet number


  /* ========== parameter section ==========*/
  // Leading Parameters
  PARAMETER_ARRAY(lnBmsy_sp);           // Biomass at MSY
  PARAMETER(lnUmsy);                    // Optimal complex exploitation rate
  PARAMETER_VECTOR(lnqFreespf_vec);     // species/stock/fleet catchability for freely estimated fleets (usually commercial)
  PARAMETER_VECTOR(lnqShrinksf_vec);    // fleet/species mean catchability for fleets that use a shrinkage prior
  PARAMETER_VECTOR(lntauspf_vec);       // survey obs error var
  PARAMETER_VECTOR(lnBinit_vec);        // Non-equilibrium initial biomass - vector so not all stocks need it
  // Priors
  PARAMETER_VECTOR(deltalnqspf_vec);    // deviation for stock catchability w/in a species/survey
  PARAMETER_VECTOR(lntauq_s);           // survey catchability sd among stocks w/in a species
  PARAMETER_VECTOR(tvlnqDevs_vec);      // Time-varying catchability devs
  PARAMETER(lntautvqDev);               // TV q SD
  PARAMETER(mlnq);                      // hyperprior mean catchability across surveys (tuning par)
  PARAMETER(sdlnq);                     // hyperprior sd in mean catchability (tuning par)
  PARAMETER_VECTOR(epslnUmsy_s);        // deviation in Umsy from complex mean to species level
  PARAMETER(lnsigUmsy);                 // complex level Umsy sd
  PARAMETER_ARRAY(epslnUmsy_sp);        // deviation in Umsy from species to stock
  PARAMETER_VECTOR(lnsigUmsy_s);        // complex level Umsy sd
  PARAMETER(mlnUmsy);                   // hyperprior mean Umsy (tuning par)
  PARAMETER(sdlnUmsy);                  // hyperprior Umsy sd (tuning par)
  PARAMETER_ARRAY(mBmsy_sp);            // prior mean eqbm biomass (tuning par)
  PARAMETER_ARRAY(sdBmsy_sp);           // prior eqbm biomass SD
  PARAMETER_ARRAY(mBinit_sp);           // Prior mean init biomass
  PARAMETER_ARRAY(sdBinit_sp);          // Prior SD for init biomass
  PARAMETER_VECTOR(tau2IGa_f);          // Inverse Gamma Prior parameters for tau2 prior
  PARAMETER_VECTOR(tau2IGb_f);          // Inverse Gamma Prior parameters for tau2 prior
  PARAMETER_VECTOR(sigma2IG);           // Inverse Gamma Prior parameters for Sigma2 prior
  PARAMETER(deltat);                    // Fractional time step used in pop dynamics to reduce chaotic behaviour
  
  // Random Effects
  PARAMETER(lnsigmaProc);               // Species effect cov matrix diag
  PARAMETER_VECTOR(zetaspt_vec);        // species-stock effect - in a vector for missing years
  PARAMETER_ARRAY(sigmaProcMult_sp);    // process error scalar mults - might be deprec later
  PARAMETER(logit_gammaYr);             // AR1 auto-corr on year effect (eps) - tuning?
  

  // State variables
  array<Type>       B_spt(nS,nP,nT+1);
  array<Type>       lnB_spt(nS,nP,nT+1);
  array<Type>       zeta_spt(nS,nP,nT+1);
  // Leading parameters
  array<Type>       Bmsy_sp(nS,nP);
  array<Type>       Binit_sp(nS,nP);
  array<Type>       tau_spf(nS,nP,nF);  
  array<Type>       tau2_spf(nS,nP,nF); 
  array<Type>       lnqhat_spf(nS,nP,nF);
  array<Type>       tau2hat_spf(nS,nP,nF);

  // Transform arrays
  Bmsy_sp = exp(lnBmsy_sp);

  // Leading scalars
  Type Umsy = exp(lnUmsy);
  
  // Prior hyperpars
  // Catchability
  array<Type>       q_spf(nS,nP,nF);
  array<Type>       lnqdev_spft(nS,nP,nF,nT);
  array<Type>       qhat_spf(nS,nP,nF);
  array<Type>       qhat_spft(nS,nP,nF,nT);
  array<Type>       deltalnq_spf(nS,nP,nF);
  array<Type>       lnq_sf(nS,nF);
  array<Type>       lnq_spf(nS,nP,nF);
  vector<Type>      tauq_s(nS);
  vector<Type>      tau2q_s(nS);
  vector<Type>      lnq_f(nF);

  lnq_f.setZero();
  q_spf.setZero();
  lnqdev_spft.setZero();
  qhat_spf.setZero();
  qhat_spft.fill(1.0);
  deltalnq_spf.setZero();
  lnq_sf.setZero();
  lnq_spf.setZero();
  tauq_s.setZero();
  tau2q_s.setZero();
  tau_spf.fill(0);

  // Umsy
  vector<Type>      Umsy_s(nS);
  array<Type>       Umsy_sp(nS,nP);
  vector<Type>      lnUmsy_s(nS);
  array<Type>       lnUmsy_sp(nS,nP);
  vector<Type>      sigUmsy_s(nS);
  // vector<Type>      sig2Umsy_s(nS);

  // Now transform and build 
  tauq_s    = exp(lntauq_s);
  tau2q_s    = exp(2. * lntauq_s);
  sigUmsy_s = exp(lnsigUmsy_s);


  Type sigUmsy = exp(lnsigUmsy);
  Type tautvqDev = exp(lntautvqDev);
  
  int tauVecIdx = 0;
  int deltalnqVecIdx = 0;
  int qDevVecIdx = 0;
  int BinitvecIdx = 0;
  int qShrinkIdx_f = 0;
  int qShrinkIdx_sf = 0;
  int qFreeIdx = 0;
  

  
  for( int s = 0; s < nS; s++ )
  {
    // First, do species level productivity
    lnUmsy_s(s) = lnUmsy + sigUmsy * epslnUmsy_s(s);

    // Loops to spread vector parameters over species/stock/fleet arrays
    for( int f = 0; f < nF; f++ )
    {
      // Now fill in speces/fleet level catchability pars
      if( shrinkq_f(f) == 1 & speciesq_sf(s,f) == 1 )
      {
        lnq_sf(s,f) = lnqShrinksf_vec(qShrinkIdx_sf);
        qShrinkIdx_sf++;
      }
    }



    for( int p = 0; p < nP; p++ )
    {

      // Umsy shrinkage decomp
      lnUmsy_sp(s,p)  = lnUmsy_s(s)  + sigUmsy_s(s) * epslnUmsy_sp(s,p);

      // Fill Binit if estimated
      if( initBioCode_sp(s,p) == 1 & solveInitBio_sp(s,p) == 0)
      {
        Binit_sp(s,p) = exp(lnBinit_vec(BinitvecIdx));
        BinitvecIdx++;
      }

      for( int f = 0; f < nF; f++ )
      {
        // Fill obs error variance
        int indicateFirstObs = 0;
        if( calcIndex_spf(s,p,f) == 1 & condMLEobsErr_f(f) == 0 )
        {
          tau_spf(s,p,f) = exp(lntauspf_vec(tauVecIdx));
          tauVecIdx++;
        }
        // Fill stock specific q devs
        if( stockq_spf(s,p,f) == 1 & condMLEq_f(f) == 0 )
        {
          deltalnq_spf(s,p,f) = deltalnqspf_vec(deltalnqVecIdx);
          deltalnqVecIdx++;
        }
        // Make stock specific survey qs using devs
        if( shrinkq_f(f) == 1 & condMLEq_f(f) == 0 )
        {
          lnq_spf(s,p,f)  = lnq_sf(s,f) + tauq_s(s) * deltalnq_spf(s,p,f);
        }
        // Now use free catchability if wanted
        if( shrinkq_f(f) == 0 & condMLEq_f(f) == 0 )
        {
          lnq_spf(s,p,f)  = lnqFreespf_vec(qFreeIdx);
          qFreeIdx++; 
        }
      }

    }

  }


  // Exponentiate for later
  q_spf     = exp(lnq_spf);
  Umsy_s    = exp(lnUmsy_s);
  Umsy_sp   = exp(lnUmsy_sp);

  // Square obs error SD
  tau2_spf = tau_spf * tau_spf;
  
  // Random Effects
  array<Type>       sigmaProc_sp(nS,nP);
  Type              gammaYr;
  
  gammaYr       = Type(1.96) / (Type(1.) + exp(Type(-2.)*logit_gammaYr) ) - Type(0.98);
  sigmaProc_sp  = exp(lnsigmaProc) * sigmaProcMult_sp;


  // Scalars
  Type              objFun  = 0.0;   // objective function (neg log likelihood)
  Type              nlpRE   = 0.0;   // process error likelihood
  Type              pospen  = 0.0;   // posfun penalty
  // Derived variables - not sure I need all this shit
  array<Type>       MSY_sp(nS,nP);
  array<Type>       lnMSY_sp(nS,nP);
  array<Type>       U_spt(nS,nP,nT);
  array<Type>       DnT_sp(nS,nP);
  array<Type>       lnDnT_sp(nS,nP);
  array<Type>       BnT_sp(nS,nP);
  array<Type>       lnBnT_sp(nS,nP);
  array<Type>       U_Umsy_spt(nS,nP,nT+1);
  array<Type>       lnU_Umsy_spt(nS,nP,nT+1);
  array<Type>       lnU_UmsyT_sp(nS,nP);

  // Fill arrays
  MSY_sp    = Bmsy_sp * Umsy_sp;
  lnMSY_sp  = log(MSY_sp);

  // ========== Procedure Section ========== //
  // Generate process errors
  
  // initialise first year effect at 0, then loop to fill in
  zeta_spt.fill(0);
  int zetaVecIdx = 0;
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
    {
      for( int t = initProcErr_sp(s,p); t < nT - 1; t++ )
      {
        Type invLogitZeta = -5 + 10 /( 1 + exp(-zetaspt_vec(zetaVecIdx)));

        zeta_spt(s,p,t) = gammaYr * zeta_spt(s,p,t-1) + sqrt(1 - square(gammaYr)) * sigmaProc_sp(s,p) * invLogitZeta;
        zetaVecIdx++;
      }
      zeta_spt(s,p,nT) = gammaYr * zeta_spt(s,p,nT-1) + 0;
    }
  
  // Initialise biomass and log-biomass
  B_spt.fill(-1);
  lnB_spt.fill(0.0);
  Type nSteps = 1 / deltat;
  // Now loop over species, reconstruct history from initT
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
    {
      // initialise population, if initBioCode(s)==0 this will be eqbm
      if( initBioCode_sp(s,p) == 0 ) B_spt(s,p,initT_sp(s,p)) = Type(2) * Bmsy_sp(s,p);
      if( initBioCode_sp(s,p) == 1 ) B_spt(s,p,initT_sp(s,p)) = Binit_sp(s,p);
      lnB_spt(s,p,initT_sp(s,p)) = log(B_spt(s,p,initT_sp(s,p)));
      for( int t = initT_sp(s,p) + 1; t < nT + 1; t++ )
      {
        Type tmpB_spt = B_spt(s,p,t-1);
        for( int dt = 0; dt < nSteps; dt ++ )
        {
          // Compute the deltat step of biomass (assuming catch is evenly distributed across the year)
          Type tmpBdt = 0;
          tmpBdt =  tmpB_spt + deltat * Umsy_sp(s,p) * tmpB_spt * (Type(2.0) - tmpB_spt / Bmsy_sp(s,p) ) - deltat*C_spt(s,p,t-1);
          tmpBdt *= exp( deltat *  zeta_spt(s,p,t) );
          // Update temp biomass
          tmpB_spt = posfun(tmpBdt, Type(1e-1), pospen);
        }
        if( t < nT)
          tmpB_spt = posfun( tmpB_spt, C_spt(s,p,t), pospen);
        
        B_spt(s,p,t) = tmpB_spt;

        lnB_spt(s,p,t) = log(B_spt(s,p,t));
      }

    }

  // Add (possibly correlated) process error devs
  nlpRE -= dnorm( zetaspt_vec, Type(0), Type(1), true).sum();
  
  // add REs to joint nll
  nlpRE += posPenFactor*pospen;
  objFun += nlpRE;

  // Concentrate species specific obs error likelihood?
  // Initialise arrays for concentrating conditional MLE qhat
  array<Type>   validObs_spf(nS,nP,nF);
  array<Type>   totObs_pf(nP,nF);
  array<Type>   z_spft(nS,nP,nF,nT);
  array<Type>   zSum_spf(nS,nP,nF);
  array<Type>   SS_spf(nS,nP,nF);
  array<Type>   totSS_pf(nP,nF);
  // Fill with 0s
  Type nllObs = 0.0;
  validObs_spf.fill(1e-6);
  totObs_pf.fill(0);
  zSum_spf.fill(0.0);
  z_spft.fill(0.0);
  SS_spf.fill(0.);
  totSS_pf.fill(0.);


  // Compute observation likelihood
  // Loop over surveys
  for( int f = 0; f < nF; f++ )
  {
    // Stock
    for( int p = 0; p < nP; p++)
    {  
      // species
      for( int s = 0; s < nS; s++ )
      {
        // time
        for( int t = initT_sp(s,p); t < nT; t++ )
        {
          // only add a contribution if the data exists (I_spft < 0 is missing)
          if( ( I_spft(s,p,f,t) > 0. ) & (calcIndex_spf(s,p,f) == 1) ) 
          {
            validObs_spf(s,p,f) += int(1);
            z_spft(s,p,f,t) = log( I_spft( s,p,f,t ) ) - log( B_spt( s,p,t ) );
            zSum_spf(s,p,f) += z_spft(s,p,f,t);
          }       
        }
        if( (condMLEq_f(f) == 1) & (calcIndex_spf(s,p,f) == 1))
        {
          // compute conditional MLE q from observation
          // Single stock model
          if( shrinkq_f(f) == 0 | nP == 1 ) 
            lnqhat_spf(s,p,f) = zSum_spf(s,p,f) / validObs_spf(s,p,f);

          // Multi-stock model
          if( nP > 1 & shrinkq_f(f) == 1 ) 
            lnqhat_spf(s,p,f) = ( zSum_spf(s,p,f) / tau2_spf(s,p,f) + lnq_sf(s,f)/tau2q_s(s) ) / ( validObs_spf(s,p,f) / tau2_spf(s,p,f) + 1 / tau2q_s(s) );  
        } else
          lnqhat_spf(s,p,f) = lnq_spf(s,p,f);

        // Exponentiate
        qhat_spf(s,p,f) = exp(lnqhat_spf(s,p,f));

        // Subtract lnq_sf from the resids
        for( int t = initT_sp(s,p); t < nT; t++ )
        {
          qhat_spft(s,p,f,t) = exp( lnqhat_spf(s,p,f) + lnqdev_spft(s,p,f,t));
          if( (I_spft(s,p,f,t) > 0.0) & (calcIndex_spf(s,p,f) == 1) )
          {
            z_spft(s,p,f,t) -= lnqhat_spf(s,p,f) - lnqdev_spft(s,p,f,t);
            SS_spf(s,p,f)   += square(z_spft(s,p,f,t));
          }
        }

        // Add valid Obs and SS to a total for each survey
        tau2hat_spf(s,p,f) = SS_spf(s,p,f) / validObs_spf(s,p,f);

        if( condMLEobsErr_f(f) == 1)
          tau2_spf(s,p,f) = tau2hat_spf(s,p,f);

        // Add to likelihood
        if( calcIndex_spf(s,p,f) == 1 )
          nllObs += fleetWeights_f(f) * 0.5 * ( validObs_spf(s,p,f)*log(tau2_spf(s,p,f)) + SS_spf(s,p,f)/tau2_spf(s,p,f));

        
      }
      
    }
  }
  objFun += nllObs;

  // Add priors
  Type nlpB = 0.0;
  Type nlpq = 0.0;
  Type nlpU = 0.0;

  // eqbm biomass prior - used for tuning MPs
  for( int p = 0; p < nP; p++)
    for( int s = 0; s < nS; s++ )
    {
      // Normal prior
      if(BPriorCode == 0)
      {
        // Bmsy
        nlpB -= dnorm(Bmsy_sp(s,p),mBmsy_sp(s,p), sdBmsy_sp(s,p), 1); 

        // Initial biomass
        if(initBioCode_sp(s,p) == 1 & solveInitBio_sp(s,p) == 0) 
          nlpB -= dnorm( Binit_sp(s,p), mBinit_sp(s,p), sdBinit_sp(s,p), 1);  

      }
      
      
    } // End biomass prior

  // Jeffreys prior
  if( BPriorCode == 1 )
  {
    nlpB += lnBmsy_sp.sum() + lnBinit_vec.sum();
  }

  // multispecies and multstock shared priors
  // I think there's no "most efficient" way to do this,
  // so, lots of conditionals!
  // First, let's do catchability
  if( lnqPriorCode == 1 )
  {
    if(nP > 1)
      nlpq -= dnorm( deltalnqspf_vec, Type(0), Type(1), true).sum();
    
  }

  if( tvq_f.sum() > 0 )
    nlpq -= dnorm( tvlnqDevs_vec, Type(0), Type(1), true).sum();

  // Then productivity
  // If number of species is greater than 1,
  // add species deviation
  if( lnUPriorCode == 1 )
  {
    if( nS > 1 )
      nlpU -= dnorm( epslnUmsy_s, Type(0), Type(1), true ).sum();

    // If number of stocks is > 1, add stock dev.
    if( nP > 1 )
      for( int p = 0; p < nP; p++ )
      {
        vector<Type> devVec = epslnUmsy_sp.col(p);
        nlpU -= dnorm( devVec, Type(0), Type(1), true ).sum();
      }
  }

  // Hyperpriors
  // Catchability
  nlpq -= dnorm( lnqShrinksf_vec, mlnq, sdlnq, true).sum();

  // Productivity
  nlpU -= dnorm( lnUmsy, mlnUmsy, sdlnUmsy, true);
  
  // Add all priors
  objFun += nlpB +  nlpq + nlpU;
  
  // Variance IG priors
  // Obs error var
  Type nllObsVarPrior = 0.;
  Type nllProcVarPrior = 0.;
  
  for( int f = 0; f < nF; f++ )
    for( int s = 0; s < nS; s++ )
      for( int p = 0; p < nP; p++)
        if( calcIndex_spf(s,p,f) == 1 & condMLEobsErr_f(f) == 0)
          nllObsVarPrior -= dinvgamma(tau2_spf(s,p,f),tau2IGa_f(f), tau2IGb_f(f), true);


  // Apply Sigma Prior
  if( SigmaPriorCode == 1 ) // Apply IG to estimated process error var element
    for( int s = 0; s < nS; s++ )
      for( int p = 0; p < nP; p++ )
        nllProcVarPrior -= dinvgamma( square(sigmaProc_sp(s,p)), sigma2IG(0), sigma2IG(1), true); 

  objFun += nllObsVarPrior + nllProcVarPrior;

  // Derive some output variables
  DnT_sp    = B_spt.col(nT-1)/Bmsy_sp/2;
  lnDnT_sp  = log(DnT_sp);
  lnBnT_sp  = lnB_spt.col(nT-1);
  for( int t = 0; t < nT; t ++ )
  {
    U_spt.col(t)      = C_spt.col(t) / B_spt.col(t);
    U_Umsy_spt.col(t) = U_spt.col(t) / Umsy_sp;
    lnU_Umsy_spt.col(t) = log(U_Umsy_spt.col(t));
  }
  lnU_UmsyT_sp = lnU_Umsy_spt.col(nT-1);



  // =========================================================== //
  // ==================== Reporting Section ==================== //
  // =========================================================== //
  
  // Variables we want SEs for
  ADREPORT(lnB_spt);
  ADREPORT(lnqhat_spf);
  ADREPORT(lnMSY_sp);
  ADREPORT(lnBinit_vec);
  ADREPORT(lnDnT_sp);
  ADREPORT(lnBnT_sp);
  ADREPORT(lnU_UmsyT_sp);

  // Model dims
  REPORT(nT);
  REPORT(nF);
  REPORT(nS);
  REPORT(nP);

  // Data
  REPORT(I_spft);
  REPORT(C_spt);

  // Switches
  REPORT( shrinkq_f );
  REPORT( stockq_spf );
  REPORT( speciesq_sf );

  // Leading Pars
  REPORT(Binit_sp);
  REPORT(Umsy);
  REPORT(Bmsy_sp);
  REPORT(tau2_spf);
  REPORT(tau_spf);
  REPORT(tau2hat_spf);
  REPORT(sigmaProc_sp);

  // State variables
  REPORT(B_spt);

  // Priors
  REPORT(mBmsy_sp);
  REPORT(sdBmsy_sp);
  REPORT(mlnUmsy);
  REPORT(sdlnUmsy);

  // Random effects
  REPORT(zeta_spt);
  REPORT(deltalnq_spf);
  REPORT(epslnUmsy_s);
  REPORT(epslnUmsy_sp);
  REPORT(gammaYr);
  REPORT(lnqdev_spft);

  // Derived variables
  REPORT(U_spt);
  REPORT(lnq_spf);
  REPORT(lnqhat_spf);
  REPORT(lnq_sf);
  REPORT(lnq_f);
  REPORT(MSY_sp);
  REPORT(qhat_spf);
  REPORT(qhat_spft);
  REPORT(q_spf);
  REPORT(Umsy_sp);
  REPORT(Umsy_s);
  REPORT(DnT_sp);
  REPORT(U_Umsy_spt);
  REPORT(qShrinkIdx_sf);
  REPORT(deltalnqspf_vec);
  REPORT(tauq_s);

  // Optimisation quantities
  REPORT(zSum_spf);
  REPORT(z_spft);
  REPORT(SS_spf);
  REPORT(validObs_spf);
  REPORT(nlpRE);
  REPORT(nllObs);
  REPORT(nlpq);
  REPORT(nlpU);
  REPORT(nlpB);
  // REPORT(nllProcVarPrior);
  // REPORT(nllObsVarPrior);
  REPORT(pospen);

  REPORT(objFun);

  // Everything else //
  
  // Return objective function value
  return objFun;
}
