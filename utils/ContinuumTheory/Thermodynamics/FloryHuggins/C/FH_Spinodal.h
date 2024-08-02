int FH_SpinodalSolventFraction(FHPAR *, double, double, int *, double *, double *);
void FH_BinarySpinodal(FHPAR *, double *, double *);


// eta:   fixed mixing ratio of A and B: eta = phiA/(phiA+phiB)
// phiG0: only used for four-component mixture -> uses approximation!
int FH_SpinodalSolventFraction(FHPAR *FHPar, double phiG0, double eta, int *N_sp, double *phi_sp1, double *phi_sp2)
{
  double  chiAC = FHPar->chiAC;
  double  chiBC = FHPar->chiBC;
  double  chiAB = FHPar->chiAB; 
  double  NA    = FHPar->NA; 
  double  NB    = FHPar->NB;
  double  tmpf1;

  int i, Nsolutions;
  double CHI, coeff[4], phiSp[3];
  double phispA,phispB;

  if(FHPar->Ncomp==2)
  {
    FH_BinarySpinodal(FHPar, &phispA , &phispB); 
    Nsolutions = (0<=phispA&phispA<=1) + (0<=phispB&phispB<=1);
  }
  else if(FHPar->Ncomp==3)
  {
    CHI = 2*(chiAC*chiBC + chiAB*(chiAC+chiBC)) 
        - chiAC*chiAC - chiBC*chiBC - chiAB*chiAB;
    coeff[3] = CHI;
    coeff[2] = -2*(CHI+chiAB - chiAC/((1.-eta)*NB) - chiBC/(eta*NA));
    coeff[1] = CHI + 4.*chiAB 
           - (2.*chiAC+1.)/((1.-eta)*NB)
            - (2.*chiBC+1.)/(eta*NA)
          + 1.0/(eta*(1-eta)*NA*NB);
    coeff[0] = 1.0/(eta*NA) + 1.0/((1.-eta)*NB) - 2.*chiAB;

    // Get volume fractions
    cubic(coeff, phiSp, &Nsolutions);
    
    Nsolutions = 0;
    for(i=0; i<3; i++)
    {
      if ( 0<= phiSp[i] & phiSp[i] <= 1)
      {
        if(Nsolutions==0)
          phispA = 1-phiSp[i];
        else
          phispB = 1-phiSp[i];
        Nsolutions++;
      }
    }  
    if(Nsolutions>1) // Two solutions: Get phispA such that this is the one with the minimum polymer fraction.
    { 
      tmpf1   = (phispA<phispB ? phispA : phispB);
      phispB = (phispA>phispB ? phispA : phispB);
      phispA = tmpf1;
    }
  }
  else if(FHPar->Ncomp==4)
  {
    chiAC *= (1.0-phiG0);
    chiBC *= (1.0-phiG0);
    chiAB *= (1.0-phiG0);
    printf("WARNING: spinodal concentration calculated based on perfectly mixed C and D components!\n");
    CHI = 2*(chiAC*chiBC + chiAB*(chiAC+chiBC)) 
        - chiAC*chiAC - chiBC*chiBC - chiAB*chiAB;
    coeff[3] = CHI;
    coeff[2] = -2*(CHI+chiAB - chiAC/((1.-eta)*NB) - chiBC/(eta*NA));
    coeff[1] = CHI + 4.*chiAB 
           -  (2.*chiAC+1.)/((1.-eta)*NB)
            - (2.*chiBC+1.)/(eta*NA)
          + 1.0/(eta*(1-eta)*NA*NB);
    coeff[0] = 1.0/(eta*NA) + 1.0/((1.-eta)*NB) - 2.*chiAB;

    // Get volume fractions
    cubic(coeff, phiSp, &Nsolutions);
    
    Nsolutions = 0;
    for(i=0; i<3; i++)
    {
      if ( 0<= phiSp[i] & phiSp[i] <= 1)
      {
        if(Nsolutions==0)
          phispA = 1-phiSp[i];
        else
          phispB = 1-phiSp[i];
        Nsolutions++;
      }
    }  

    if(Nsolutions>1) // Two solutions: Get phispA such that this is the one with the minimum polymer fraction.
    { 
      tmpf1   = (phispA<phispB ? phispA : phispB);
      phispB = (phispA>phispB ? phispA : phispB);
      phispA = tmpf1;
    }
  }
  else  
    {printf("ERROR: Ncomponents should habe value 2,3, or 4!\n"); return(0);}


  (*N_sp)=Nsolutions;
  (*phi_sp1) = phispA;
  (*phi_sp2) = phispB;


  return(1);
}

void FH_BinarySpinodal(FHPAR * FHPar, double *phi_sp1, double *phi_sp2)
{
  double A, B, C, sqrtD;

  A = 2.0*FHPar->chiAB;
  B = 1.0/FHPar->NB -1.0/FHPar->NA - A;
  C = 1.0/FHPar->NA;
  sqrtD = sqrt(B*B - 4.0*A*C);
  *phi_sp1 = (-B - sqrtD)/(2.0*A);
  *phi_sp2 = (-B + sqrtD)/(2.0*A);
}
