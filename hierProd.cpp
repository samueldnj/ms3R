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



// objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Call namespaces //
  using namespace density;

  /*data section*/
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
  DATA_INTEGER(tauqPriorCode);  // 0 => IG on tauq2, 1 => normal
  DATA_INTEGER(sigUPriorCode);  // 0 => IG on sigU2, 1 => normal
  DATA_INTEGER(condMLEq);       // 0 => q leading par, 1 => q concentrated
  DATA_INTEGER(lnqPriorCode);   // 0 => hyperprior, 1 => multilevel
  DATA_INTEGER(lnUPriorCode);   // 0 => hyperprior, 1 => multilevel 
  DATA_INTEGER(BPriorCode);     // 0 => normal, 1 => Jeffreys 
  DATA_IARRAY(initT_sp);        // first year of assessment
  DATA_IARRAY(initBioCode_sp);  // initial biomass at 0 => unfished, 1=> fished
  DATA_SCALAR(posPenFactor);    // Positive-penalty multiplication factor


  /*parameter section*/
  // Leading Parameters
  PARAMETER_ARRAY(lnBmsy_sp);           // Biomass at MSY
  PARAMETER_ARRAY(lnUmsy_sp);           // Optimal exploitation rate            
  PARAMETER_ARRAY(lnq_spf);             // Survey-species-stock catchability
  PARAMETER_VECTOR(lntau_spf);          // survey obs error var
  PARAMETER_ARRAY(lnBinit_sp);          // Non-equilibrium initial biomass
  // Priors
  PARAMETER_ARRAY(lnqbar_sf);           // survey mean catchability
  PARAMETER_ARRAY(lntauq_sf);           // survey catchability sd 
  PARAMETER(mq);                        // hyperprior mean catchability
  PARAMETER(sq);                        // hyperprior catchability sd
  PARAMETER_VECTOR(lnUmsybar_s);        // shared prior mean Umsy
  PARAMETER_VECTOR(lnsigUmsy_s);        // shared prior Umsy sd
  PARAMETER(mUmsy);                     // hyperprior mean Umsy
  PARAMETER(sUmsy);                     // hyperprior Umsy sd
  PARAMETER_VECTOR(mBmsy);              // prior mean eqbm biomass
  PARAMETER_VECTOR(sBmsy);              // prior eqbm biomass var
  PARAMETER_VECTOR(tau2IGa);            // Inverse Gamma Prior parameters for tau2 prior
  PARAMETER_VECTOR(tau2IGb);            // Inverse Gamma Prior parameters for tau2 prior
  PARAMETER_VECTOR(tauq2Prior);         // Hyperparameters for tauq2 prior - (IGa,IGb) or (mean,var)
  PARAMETER_VECTOR(sigU2Prior);         // Hyperparameters for sigU2 prior - (IGa,IGb) or (mean,var)
  PARAMETER_VECTOR(kappa2IG);           // Inverse Gamma Prior parameters for kappa2 prior
  PARAMETER_VECTOR(Sigma2IG);           // Inverse Gamma Prior parameters for Sigma2 prior
  PARAMETER_MATRIX(wishScale);          // IW scale matrix for Sigma prior
  PARAMETER(nu);                        // IW degrees of freedom for Sigma prior    
  PARAMETER(deltat);                    // Fractional time step used in pop dynamics to reduce chaotic behaviour
  
  // Random Effects
  PARAMETER(lnsigmaProc);               // Species effect cov matrix diag
  PARAMETER_VECTOR(zetaspt_vec);        // species-stock effect - in a vector for missing years
  PARAMETER_ARRAY(sigmaProcMult_sp);    // process error scalar mults - might be deprec later
  PARAMETER(logit_gammaYr);             // AR1 auto-corr on year effect (eps) - tuning?
  

  // State variables
  array<Type>       B_spt(nS,nP,nT);
  array<Type>       lnBt_spt(nS,nP,nT);
  array<Type>       zeta_spt(nS,nP,nT);
  // Leading parameters
  array<Type>       Bmsy_sp(nS,nP);
  array<Type>       Umsy_sp(nS,nP);
  array<Type>       Binit_sp(nS,nP);
  array<Type>       tau_spf(nS,nP,nF);  
  array<Type>       tau2_spf(nS,nP,nF); 
  array<Type>       lnqhat_spf(nS,nP,nF);
  array<Type>       tau2hat_pf(nP,nF);

  Bmsy_sp = exp(lnBmsy_sp);
  Umsy_sp = exp(lnUmsy_sp);
  Binit_sp = exp(lnBinit_sp);
  tau_spf = exp(lntau_spf);
  tau2_spf = exp(2.*lntau_spf);
  
  // Prior hyperpars
  array<Type>       qbar_sf(nS,nF);
  array<Type>       tauq_sf(nS,nF);
  array<Type>       tauq2_sf(nS,nF);
  Type              Umsybar_s   = exp(lnUmsybar_s);
  Type              sigUmsy2    = exp(Type(2) * lnsigUmsy_s);
  
  // Random Effects
  array<Type>       sigmaProc_sp(nS,nP);
  Type              gammaYr;
  
  gammaYr       = Type(1.96) / (Type(1.) + exp(Type(-2.)*logit_gammaYr) ) - Type(0.98);
  qbar_sf       = exp(lnqbar_sf);
  tauq_sf       = exp(lntauq_sf);
  tauq2_sf      = exp(Type(2) * lntauq_sf);
  sigmaProc_sp  = exp(lnsigmaProc) * sigmaProcMult_sp;


  // Scalars
  Type              nll     = 0.0;   // objective function (neg log likelihood)
  Type              nllRE   = 0.0;   // process error likelihood
  Type              pospen  = 0.0;   // posfun penalty
  // Derived variables - not sure I need all this shit
  array<Type>       MSY_sp(nS,nP);
  array<Type>       lnMSY_sp(nS,nP);
  array<Type>       U_spt(nS,nP,nT);
  array<Type>       DnT_sp(nS,nP);
  array<Type>       lnDnT_sp(nS,nP);
  array<Type>       BnT_sp(nS,nP);
  array<Type>       lnBnT_sp(nS,nP);
  array<Type>       U_Umsy_spt(nS,nP,nT);
  array<Type>       lnU_Umsy_spt(nS,nP,nT);

  // Fill arrays
  MSY_sp    = Bmsy_sp * Umsy_sp;
  lnMSY_sp  = log(MSY_sp);

  // Procedure Section //
  // First create the correlation matrix for the species effects
  UNSTRUCTURED_CORR_t<Type> specEffCorr(SigmaOffDiag);
  SigmaCorr = specEffCorr.cov();
  // Standard deviation component
  SigmaD.fill(0.0);
  SigmaD.diagonal() = sqrt(SigmaDiag);
  // Now combine
  Sigma = (SigmaD * SigmaCorr);
  Sigma *= SigmaD;

  // Generate Population Dynamics
  // initialise first year effect at 0, then loop to fill in
  zeta_spt.fill(0);
  int zetaVecIdx = 0;
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
    {
      int corrMatEntry = s * (nP) + p;
      for( int t = initT_sp(s,p) + 1; t < nT; t++ )
      {
        zeta_spt(s,p,t) = gammaYr * zeta_spt(s,p,t) + (1 - gammaYr) * SigmaDiag(corrMatEntry) * zetaspt_vec(zetaVecIdx);
        zetaVecIdx++;
      }
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
      if( initBioCode_sp(s,p) == 0 ) Bt(s,p,initT_sp(s,p)) = Type(2) * Bmsy_sp(s,p);
      if( initBioCode_sp(s,p) == 1 ) Bt(s,p,initT_sp(s,p)) = Binit_sp(s,p);
      lnB_spt(s,p,initT_sp(s,p)) = log(B_spt(s,p,initT_sp(s,p)));
      for( int t = initT_sp(s,p)+1; t < nT; t++ )
      {
        Type tmpB_spt = B_spt(s,p,t-1);
        for( int dt = 0; dt < nSteps; dt ++ )
        {
          // Compute the deltat step of biomass (assuming catch is evenly distributed across the year)
          Type tmpBdt = 0;
          tmpBdt =  tmpB_spt + deltat * Umsy_sp(s,p) * tmpB_spt * (Type(2.0) - tmpB_spt / Bmsy(s) ) - deltat*Ct(s,p,t-1);
          tmpBdt *= exp( deltat *  zeta_spt(s,p,t-1) );
          // Now update tmpB_spt
          tmpB_spt = posfun(tmpBdt, Type(1e-3), pospen);
        }
        B_spt(s,p,t) = tmpB_spt;
        lnB_spt(s,p,t) = log(B_spt(s,p,t));
      }

    }

  // Add (possibly correlated) process error devs
  nllRE -= dnorm( zetaspt_vec, Type(0), Type(1), true).sum();
  
  // add REs to joint nll
  nllRE += posPenFactor*pospen;
  nll += nllRE;

  // Concentrate species specific obs error likelihood?
  // Initialise arrays for concentrating conditional MLE qhat
  array<Type>   validObs_spf(nS,nP,nF);
  vector<Type>  totObs_pf(nP,nF);
  array<Type>   qhat_sf(nS,nP,nF);
  array<Type>   z_spft(nS,nP,nF,nT);
  array<Type>   zSum_spf(nS,nP,nF);
  array<Type>   SS_spf(nS,nP,nF);
  vector<Type>  totSS_pf(nP,nF);
  // Fill with 0s
  Type nllObs = 0.0;
  validObs_spf.fill(1e-6);
  totObs_spf.fill(0);
  zSum_spf.fill(0.0);
  z_spft.fill(0.0);
  qhat_spf.fill(-1.0);
  SS_spf.fill(0.);
  totSS_spf.fill(0.);


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
          // only add a contribution if the data exists (Iost < 0 is missing)
          if( ( I_spft(s,p,f,t) > 0. ) ) 
          {
            validObs_spf(s,p,f) += int(1);
            z_spft(s,p,f,t) = log( I_spft( s,p,f,t ) ) - log( B_spt( s,p,t ) );
            zSum_spf(s,p,f) += z_spft(s,p,f,t);
          }       
        }
        if(condMLEq == 1)
        {
          // compute conditional MLE q from observation
          // Single stock model
          if( nP == 1 | lnqPriorCode == 0) 
            lnqhat_spf(s,p,f) = zSum_spf(s,p,f) / validObs_spf(s,p,f);

          // Multi-stock model
          if( nP > 1 & lnqPriorCode == 1 ) 
            lnqhat_spf(s,p,f) = ( zSum_spf(s,p,f) / tau2_spf(s,p,f) + lnqbar_sf(s,f)/tauq2_sf(s,f) ) / ( validObs_spf(s,p,f) / tau2_spf(s,p,f) + 1 / tauq2_sf(s,f) );  
        } else
          lnqhat_spf(s,p,f) = lnq_spf(s,p,f);

        // Exponentiate
        qhat_spf(s,p,f) = exp(lnqhat_spf(s,p,f));

        // Subtract lnq_sf from the resids
        for( int t = initT_sp(s,p); t < nT; t++ )
          if( (I_spft(s,p,f,t) > 0.0) )
          {
            z_spft(s,p,f,t) -= lnqhat_sf(o,s);
            SS_spf(s,p,f)   += square(z_spft(s,p,f,t));
          }

        // Add to likelihood
        nllObs += 0.5 * ( validObs_spf(s,p,f)*lntau2_spf(s,p,f) + SS_spf(s,p,f)/tau2_spf(s,p,f));

        // Add valid Obs and SS to a total for each survey
        totObs_pf(p,f) += validObs_spf(s,p,f);
        totSS_pf(p,f) += SS_spf(s,p,f);
      }
      tau2hat_pf(p,f) = totSS_pf(p,f) / totObs_pf(p,f);
    }
  }
  nll += nllObs;

  // Add priors
  Type nllBprior = 0.0;
  Type nllqPrior = 0.0;
  Type nllUprior = 0.0;

  // eqbm biomass prior - used for tuning MPs
  for( int p = 0; p < nP; p++)
    for( int s = 0; s < nS; s++ )
    {
      // Normal prior
      if(BPriorCode == 0)
      {
        // Bmsy
        nllBprior -= dnorm(Bmsy_sp(s,p),mBmsy_sp(s,p), sBmsy_sp(s,p), 1); 

        // Initial biomass
        if(initBioCode_sp(s,p) == 1) 
          nllBprior -= dnorm( Binit_sp(s,p), mBmsy_sp(s,p)/2, sBmsy_sp(s,p)/2, 1);  

      }
      // Jeffreys prior
      if( BPriorCode == 1 )
        nllBprior += lnBmsy_sp(s,p) + lnBinit_sp(s,p);
      
    } // End biomass prior

  // multispecies shared priors
  // If not using shared priors, use the hyperpriors for all
  // species specific parameters
  if( nS > 1 | nP > 1 )
  {
    // 1st level priors
    for( int s = 0; s < nS; s++ )
    {
      // In here, loop over surveys
      for( int o = 0; f < nF; o++ )
      {
        // catchability
        // Shared Prior
        if( lnqPriorCode == 1 & condMLEq == 0 )
          nllqPrior -= dnorm( lnqhat_sf(o,s), lnqbar_f(o), tauq_f(o), 1);

        // No shared prior (uses the same prior as the SS model)
        if( lnqPriorCode == 0 )
          nllqPrior += 0.5 * pow((qhat_sf(0,s) - mq)/sq,2);

      }
      // productivity
      // Shared Prior
      if( lnUPriorCode == 1 )
        nllUprior -= dnorm( lnUmsy(s), lnUmsybar, sqrt(sigUmsy2), 1);

      // No shared prior (uses SS model prior)
      if( lnUPriorCode == 0 )
        nllUprior += pow( (Umsy(s) - mUmsy)/ sUmsy, 2);

    }  
    // Hyperpriors
    // catchability
    if( lnqPriorCode == 1 )
      for( int o = 0; f < nF; o++ ) 
        nllqPrior -= dnorm( qbar_f(o), mq, sq, 1); 

    // productivity
    if( lnUPriorCode == 1 )
      nllUprior -= dnorm(Umsybar, mUmsy, sUmsy, 1);
    
  } // End multispecies shared priors    
  
  // Now for single species model
  if( nS == 1 ) 
  { 
    // catchability
    for (int f = 0; o<nF; o++)
      nllqPrior -= dnorm( qhat_sf(o,0), mq, sq, 1); 
    
    // productivity
    nllUprior -= dnorm( Umsy(0), mUmsy, sUmsy, 1);
  }
  // Add all priors
  nll += nllBprior +  nllqPrior + nllUprior;
  
  // Variance IG priors
  // Obs error var
  Type nllVarPrior = 0.;
  Type nllSigPrior = 0.;
  for( int o = 0; f < nF; o++ )
    nllVarPrior += (tau2IGa(o)+Type(1))*lntau2_f(o)+tau2IGb(o)/tau2_f(o);  

  // year effect deviations var
  if( kappaPriorCode == 1 )  
    nllVarPrior += (kappa2IG(0)+Type(1))*lnkappa2 + kappa2IG(1)/kappa2;

  // Now multispecies priors
  if (nS > 1)
  {
    // shared q prior variance
    if( tauqPriorCode == 0 )
    {
      for( int o = 0; f < nF; o++)
      {
        nllVarPrior += (tauq2Prior(0)+Type(1))*Type(2)*lntauq_f(o)+tauq2Prior(1)/tauq2_f(o);
      }
    }
    if( tauqPriorCode == 1 )
    {
      for( int o = 0; f < nF; o++)
      {
        nllVarPrior += Type(0.5) * pow( tauq2_f(o) - tauq2Prior(0), 2) / tauq2Prior(1);
      }
    }
    // shared U prior variance
    // IG
    if( sigUPriorCode == 0 )
    {
      nllVarPrior += (sigU2Prior(0)+Type(1))*Type(2.)*lnsigUmsy+sigU2Prior(1)/sigUmsy2;
    }
    // Normal
    if( sigUPriorCode == 1 )
    {
      nllVarPrior += Type(0.5) * pow( sigUmsy2 - sigU2Prior(0), 2) / sigU2Prior(1);
    }
    
    if( SigmaPriorCode == 1 ) // Apply IW prior to Sigma matrix
    {
      matrix<Type> traceMat = wishScale * Sigma.inverse();
      Type trace = 0.0;
      for (int s=0;s<nS;s++) trace += traceMat(s,s);
      nllSigPrior += Type(0.5) *( (nu + nS + 1) * atomic::logdet(Sigma) + trace);
    }
  }

  // Apply Sigma Prior
  if( SigmaPriorCode == 0 ) // Apply IG to estimated SigmaDiag element
    nllSigPrior += (Sigma2IG(0)+Type(1))*lnSigmaDiag+Sigma2IG(1)/exp(lnSigmaDiag);

  nll += nllVarPrior + nllSigPrior;

  // Derive some output variables
  Ut      = Ct / Bt;
  DnT     = Bt.col(nT-1)/Bmsy/2;
  lnDnT   = log(DnT);
  lnBnT   = log(Bt.col(nT-1));
  for( int t = 0; t < nT; t ++ )
  {
    U_Umsy.col(t) = Ut.col(t) / Umsy;
    lnU_Umsy.col(t) = log(U_Umsy.col(t));
  }
  lnU_UmsyT = lnU_Umsy.col(nT-1);

  // Reporting Section //
  // Variables we want SEs for
  ADREPORT(lnBt);
  ADREPORT(lnqhat_sf);
  ADREPORT(lnMSY);
  ADREPORT(lnBinit);
  ADREPORT(lntau2_f);
  ADREPORT(lnkappa2);
  ADREPORT(lnDnT);
  ADREPORT(lnBnT);
  ADREPORT(lnU_UmsyT);

  
  // Everything else //
  REPORT(Bt);
  REPORT(It);
  REPORT(Ut);
  REPORT(Ct);
  REPORT(eps_t);
  REPORT(omegat);
  REPORT(Binit);
  REPORT(DnT);
  REPORT(U_Umsy);
  REPORT(qhat_sf);
  REPORT(MSY);
  REPORT(Bmsy);
  REPORT(Umsy);
  REPORT(tau2_f);
  REPORT(kappa2);
  REPORT(nT);
  REPORT(nF);
  REPORT(nS);
  REPORT(zSum_sf);
  REPORT(z_sft);
  REPORT(validObs);
  REPORT(zeta_st);
  REPORT(Sigma);
  REPORT(SigmaDiag);
  if (nS > 1)
  {
    REPORT(SigmaCorr);
    REPORT(Umsybar);
    REPORT(sigUmsy2);
    REPORT(qbar_f);
    REPORT(tauq2_f);
  }
  REPORT(gammaYr);
  REPORT(nll);
  REPORT(nllRE);
  REPORT(nllObs);
  REPORT(nllqPrior);
  REPORT(nllUprior);
  REPORT(nllBprior);
  REPORT(nllSigPrior);
  REPORT(nllVarPrior);
  REPORT(pospen);
  
  return nll;
}
