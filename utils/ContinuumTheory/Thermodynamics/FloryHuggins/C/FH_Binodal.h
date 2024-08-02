int findBinaryChemicalPotential(FHPAR *FHPar, double tol, double mu_target, double phiL, double phiU)
{
  double mu, muBE1, muSE1, err=(phiU-phiL)/phiL;
  while(err>tol)
  {
    FHPar->phiA = 0.5*( phiL + phiU );
    if(!FH_ExchangeChemicalPotential(FHPar, &mu, &muBE1, &muSE1 ))
      {printf("Error: Failed to execute FH_ExchangeChemicalPotential()!\n"); return(0);}
    if(mu>mu_target)
      phiU=FHPar->phiA;
    else if(mu<mu_target)
      phiL=FHPar->phiA;
    else
      {phiL=(phiU=FHPar->phiA);}
    err=(phiU-phiL)/phiL;
  }

  FHPar->phiA = 0.5*( phiL + phiU ); // Result
  return(1);
}

int binaryCostFunction(FHPAR *FHPar, double tol, double phi_sp1, double phi_sp2, double muTarget, double *cost)
{
  double phi_1,phi_2, f_1, f_2;

  findBinaryChemicalPotential(FHPar, tol, muTarget,       0, phi_sp1);
  phi_1=FHPar->phiA;
  if(!FH_LocalFreeEnergy(FHPar, &f_1))
    {printf("Error: Failed to execute FH_LocalFreeEnergy()!\n"); return(0);}

  findBinaryChemicalPotential(FHPar, tol, muTarget, phi_sp2,       1);
  phi_2=FHPar->phiA;
  if(!FH_LocalFreeEnergy(FHPar, &f_2))
    {printf("Error: Failed to execute FH_LocalFreeEnergy()!\n"); return(0);}

  *cost = f_2 - f_1 - (phi_2-phi_1)*muTarget;

  return(1);
}

int FH_BinaryBinodal(FHPAR *FHPar, double tol, int *N_bn, double *phi_bn1, double *phi_bn2)
{
  int N_sp;
  double chiABcr, phiAcr, f;
  double phi_sp1, phi_sp2, f_sp1, f_sp2;

  double err, cost;
  double muBE1, muSE1, mu, mu_L, mu_U;


  // Numerical parameters
//  double   delta    =0.1,       // Initial step size
//          tolerance =1e-4,   // Tolerance level: Controls the step size
//              tolY  = 1e-7;

  //-----------------------------------------------------------------------------------------
  // Calculate critical point
  FH_BinaryCriticalChi(FHPar, &chiABcr);

  //-----------------------------------------------------------------------------------------
  // Below critical point
  if(FHPar->chiAB < chiABcr) 
  {
    printf("Error: chiAB=%f smaller than chiABcr=%f: No phase separation!\n", FHPar->chiAB, chiABcr);
    return(0);
  }
  FH_BinaryCriticalConcentration(FHPar, &phiAcr);

  //-----------------------------------------------------------------------------------------
  // At critical point
  if(FHPar->chiAB == chiABcr)
  {
    *phi_bn1=phiAcr; 
    *phi_bn2=phiAcr;
    return(1);
  }

  //-----------------------------------------------------------------------------------------
  // Above critical point: Common tangent construction
  // 
  // Method: determine chemical potential from interval [mu_L, mu_U] using the bisection method
  // the cost function is f2 - f1 - (phi2-phi1)*mu; indices 1 and 2 label the coexisting phases


  // Determine upper and lower values mu_L and mu_U using the spinodal concentrations
  FH_SpinodalSolventFraction(FHPar, 0, 0, &N_sp, &phi_sp1, &phi_sp2);
  if(N_sp!=2)
  {
    printf("Error: FH_SpinodalSolventFraction() predicts %d spinodal concentrations",N_sp );
  }
  FHPar->phiA=phi_sp2;
  if(!FH_ExchangeChemicalPotential(FHPar, &mu_L, &muBE1, &muSE1 ))
      {printf("Error: Failed to execute FH_ExchangeChemicalPotential()!\n"); return(0);}
  FHPar->phiA=phi_sp1;
  if(!FH_ExchangeChemicalPotential(FHPar, &mu_U, &muBE1, &muSE1 ))
      {printf("Error: Failed to execute FH_ExchangeChemicalPotential()!\n"); return(0);}

  // Bisection method
  err=(mu_U-mu_L)/mu_L;err=(err>0?err:-err);
  while (err>tol)
  {
    mu=0.5*( mu_L + mu_U);
    binaryCostFunction(FHPar, tol, phi_sp1, phi_sp2, mu, &cost);
    if(cost>0)
      mu_L=mu;
    else if (cost<0)
      mu_U=mu;
    else
      {mu_L=(mu_U=mu);}
    err=(mu_U-mu_L)/mu_L;err=(err>0?err:-err);
  }

  // Get concentrations at common chemical potential
  findBinaryChemicalPotential(FHPar, tol, mu,       0, phi_sp1);
  *phi_bn1=FHPar->phiA;
  findBinaryChemicalPotential(FHPar, tol, mu, phi_sp2, 1      );
  *phi_bn2=FHPar->phiA;


  return(1);
}




