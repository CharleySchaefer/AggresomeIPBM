/*

  Flory-Huggins model for two, three, and four components, A, B, C, D.

  The functions in this header file calculate the
  
    - local free energy           (f)
    - exchange chemical potential (df/dphi) 
    - bare chemical potential     (df/dN)
    - Hessian of the free energy  (d2f/dphi^2)

  The parameters include the (1) volume fractions, (2) molecular volumes, and (3) interaction parameters.
  

  (1) Volume fraction of each component:
      phiA, phiB, phiC, phiD.
 
      FH theory assumes incompressibility, hence, phiA + phiB + phiC + phiD = 1 must hold.

  (2) Molecular volumes:
      NA, NB, NC, ND

  (3) interaction parameters:
      chiAB, chiAC, etc.
*/
int  FH_LocalFreeEnergy(FHPAR *,double *);

int  FH_ExchangeChemicalPotential(FHPAR *, double *, double *, double *);
int  FH_BareChemicalPotential(    FHPAR *, double *, double *, double *, double *);

void FH_FreeEnergyHessianTernary( FHPAR *,  double *, double *, double *);
int  FH_FreeEnergyHessianDeterminant(FHPAR *, double *);


int FH_LocalFreeEnergy(FHPAR *FHPar, double *f)
{
  double NA, NB, NC, ND;
  double phiA, phiB, phiC, phiD;
  double chiAB,chiAC,chiAD,chiBC, chiBD, chiCD;
  
  NA=FHPar->NA; 
  NB=FHPar->NB;
  NC=FHPar->NC;
  ND=FHPar->ND;
  phiA=FHPar->phiA;
  phiB=FHPar->phiB;
  phiC=FHPar->phiC;
  phiD=FHPar->phiD;
  chiAB=FHPar->chiAB;
  chiAC=FHPar->chiAC;
  chiAD=FHPar->chiAD;
  chiBC=FHPar->chiBC;
  chiBD=FHPar->chiBD;
  chiCD=FHPar->chiCD;

  if(FHPar->Ncomp==2)
  {
    phiB=1.0-phiA;
    if(phiA==0 | phiA ==1)
      *f=0.0;
    else
      *f = phiA*log(phiA)/NA + phiB*log(phiB)/NB +phiA*phiB*chiAB;
  }
  else if(FHPar->Ncomp==3)
  {
    phiC=1.0-phiA-phiB;
    *f = (phiA==0 ? 0 : phiA*log(phiA)/NA) + (phiB==0 ? 0 : phiB*log(phiB)/NB) + (phiC==0 ? 0 : phiC*log(phiC)/NC)
       + phiA*phiB*chiAB + phiA*phiC*chiAC + phiB*phiC*chiBC;  
  }
  else if(FHPar->Ncomp==4)
  {
    phiD=1.0-phiA-phiB-phiC;
    *f = (phiA==0 ? 0 : phiA*log(phiA)/NA) + (phiB==0 ? 0 : phiB*log(phiB)/NB) + (phiC==0 ? 0 : phiC*log(phiC)/NC) + (phiD==0 ? 0 : phiD*log(phiD)/ND)
       + phiA*phiB*chiAB   + phiA*phiC*chiAC   + phiB*phiC*chiBC
       + phiA*phiD*chiAD   + phiB*phiD*chiBD   + phiC*phiD*chiCD;
  }
  else
    {printf("ERROR: Ncomponents should have value 2,3, or 4!\n"); return(0);}

  return(1);
}



int FH_ExchangeChemicalPotential(FHPAR *FHPar, double *muA, double *muB, double *muC )
{
  double NA, NB, NC, ND;
  double phiA, phiB, phiC, phiD;
  double chiAB,chiAC,chiAD,chiBC, chiBD, chiCD;
  
  NA=FHPar->NA; 
  NB=FHPar->NB;
  NC=FHPar->NC;
  ND=FHPar->ND;
  phiA=FHPar->phiA;
  phiB=FHPar->phiB;
  phiC=FHPar->phiC;
  phiD=FHPar->phiD;
  chiAB=FHPar->chiAB;
  chiAC=FHPar->chiAC;
  chiAD=FHPar->chiAD;
  chiBC=FHPar->chiBC;
  chiBD=FHPar->chiBD;
  chiCD=FHPar->chiCD;

 // double phiApowm2, phiBpowm2, phiSpowm2,phiGpowm2;
  
  if(FHPar->Ncomp==2) 
  {
    phiB=1.0-phiA;
    *muA = (1.0 + log(phiA))/NA - (1.0+log(phiB) )/NB + (phiB - phiA)*chiAB ;
  }  
  else if(FHPar->Ncomp==3)  
  {
    phiC=1.0-phiA-phiB;
    *muA = ((1.0 + log(phiA))/NA - (1.0+log(phiC) )
                 + (phiC - phiA)*chiAC + phiB*(chiAB-chiBC));
    *muB = ((1.0 + log(phiB))/NB - (1.0+log(phiC) )
           + (phiC - phiB)*chiBC + phiA*(chiAB-chiAC));
  }
  else if(FHPar->Ncomp==4)
  {
    phiD=1.0-phiA-phiB-phiC;
    *muA = (1.0 + log(phiA))/NA - (1.0+log(phiD))
             + phiB*(chiAB - chiBD)
             + phiC*(chiAC - chiCD)
             + (phiD- phiA)*chiAD;
    *muB = (1.0 + log(phiB))/NB - (1.0+log(phiD))
             + phiA*(chiAB - chiAD)
             + phiC*(chiBC - chiCD)
             + (phiD- phiB)*chiBD;  
    *muC = (1.0 + log(phiC))    - (1.0+log(phiD))
             + phiA*(chiAC - chiAD)
             + phiB*(chiBC - chiBD)
             + (phiD- phiC)*chiCD;
             
    //------------------------
/*    if  (phiD < Num->threshphi) 
    {
      phiApowm2 = 1.0/(phiA*phiA)-Num->threshphi2i;
      phiBpowm2 = 1.0/(phiB*phiB)-Num->threshphi2i;
      phiSpowm2 = 1.0/(phiSlocal*phiSlocal) - Num->threshphi2i;
      phiGpowm2 = 1.0/(phiD*phiD) - Num->threshphi2i;
      *muA += Num->threshsoft*(phiGpowm2 - phiApowm2);
      *muB += Num->threshsoft*(phiGpowm2 - phiBpowm2);
      *muS += Num->threshsoft*(phiGpowm2 - phiSpowm2);
    }
    else 
    {
      if (phiA < Num->threshphi)
      {
        phiApowm2 = 1.0/(phiA*phiA)-Num->threshphi2i;
        phiGpowm2 = 1.0/(phiD*phiD) - Num->threshphi2i;
        *muA += Num->threshsoft*(phiGpowm2 - phiApowm2);
      }
      if (phiB < Num->threshphi)
      {
        phiBpowm2 = 1.0/(phiB*phiB)-Num->threshphi2i;
        phiGpowm2 = 1.0/(phiD*phiD) - Num->threshphi2i;
        *muB += Num->threshsoft*(phiGpowm2 - phiBpowm2);
      }
      if (phiC < Num->threshphi)
      {
        phiGpowm2 = 1.0/(phiD*phiD) - Num->threshphi2i;
        phiSpowm2 = 1.0/(phiC*phiC) - Num->threshphi2i;
        *muS += Num->threshsoft*(phiGpowm2 - phiSpowm2);
      }
    }*/
  }    
  else
    {printf("ERROR: Ncomponents should have value 2,3, or 4!\n"); return(0);}
  return(1);
}


int FH_BareChemicalPotential(FHPAR *FHPar, double *muA, double *muB, double *muC, double *muN )
{
  double f;
  double muAex,muBex,muCex; // exchange chemical potential
  
  FH_ExchangeChemicalPotential(FHPar, &muAex, &muBex, &muCex );
  FH_LocalFreeEnergy(FHPar, &f);
  
  if (FHPar->Ncomp==2)
  {
    *muN = f - FHPar->phiA*muAex;
    *muA = *muN + muAex;
  }
  else if (FHPar->Ncomp==3)
  {
    *muN = f - FHPar->phiA*muAex- FHPar->phiB*muBex;
    *muA = *muN + muAex;
    *muB = *muN + muBex;
  }
  else if (FHPar->Ncomp==4)
  {
    *muN = f - FHPar->phiA*muAex- FHPar->phiB*muBex- FHPar->phiC*muCex;
    *muA = *muN + muAex;
    *muB = *muN + muBex;
    *muC = *muN + muCex;
  }
  return(1);
}

// TODO: correct for NC and ND values other than 1
void FH_FreeEnergyHessianTernary(FHPAR *FHPar, double *FAA, double *FBB, double *FAB)
{
  double NA, NB, NC, ND;
  double phiA, phiB, phiC, phiD;
  double chiAB,chiAC,chiAD,chiBC, chiBD, chiCD;
  double firstterm;
  
  NA=FHPar->NA; 
  NB=FHPar->NB;
  NC=FHPar->NC;
  //ND=FHPar->ND;
  phiA=FHPar->phiA;
  phiB=FHPar->phiB;
  phiC=FHPar->phiC;
  //phiD=FHPar->phiD;
  chiAB=FHPar->chiAB;
  chiAC=FHPar->chiAC;
  //chiAD=FHPar->chiAD;
  chiBC=FHPar->chiBC;
  //chiBD=FHPar->chiBD;
  //chiCD=FHPar->chiCD;

  phiC=1.0-phiA-phiB;
  firstterm=1.0/(phiC);
  *FAA = firstterm + 1.0/(phiA*NA) - 2.0*chiAC;
  *FBB = firstterm + 1.0/(phiB*NB) - 2.0*chiBC;
  *FAB = firstterm + chiAB - chiAC - chiBC;
}

int FH_FreeEnergyHessianDeterminant(FHPAR *FHPar, double *determinant)
{
  double NA, NB, NC, ND;
  double phiA, phiB, phiC, phiD;
  double chiAB,chiAC,chiAD,chiBC, chiBD, chiCD;
  double firstterm;
  double faa, fbb, fss, fab, fas, fbs;
  double tmpf, tmpf2 ,tmpf3;
  
  NA=FHPar->NA; 
  NB=FHPar->NB;
  NC=FHPar->NC;
  ND=FHPar->ND;
  phiA=FHPar->phiA;
  phiB=FHPar->phiB;
  phiC=FHPar->phiC;
  phiD=FHPar->phiD;
  chiAB=FHPar->chiAB;
  chiAC=FHPar->chiAC;
  chiAD=FHPar->chiAD;
  chiBC=FHPar->chiBC;
  chiBD=FHPar->chiBD;
  chiCD=FHPar->chiCD;
  
  if(FHPar->Ncomp==2)
  {
    phiB = 1.0-phiA;
    faa = 1.0/(NA*phiA) + 1.0/(NB*phiB) - 2.0*chiAB;
    *determinant = faa;
  }
  else if(FHPar->Ncomp==3)
  {
    FH_FreeEnergyHessianTernary(FHPar, &faa, &fbb, &fab);
    *determinant = faa*fbb - fab*fab;
  }
  else if(FHPar->Ncomp==4)
  {
    phiD=1.0-phiA-phiB-phiC; tmpf = 1.0/(phiD);
    faa = tmpf + 1.0/(phiA*NA) - 2*chiAD;
    fbb = tmpf + 1.0/(phiB*NB) - 2*chiBD;
    fss = tmpf + 1.0/(phiC           ) - 2*chiCD;
    fab = tmpf +  2.0*chiAB - chiAD - chiBD;
    fas = tmpf +  2.0*chiAC - chiAD - chiCD;
    fbs = tmpf +  2.0*chiBC - chiBD - chiCD;
    tmpf  = faa*fbb - fab*fab;
    tmpf2 = faa*fbs - fab*fas;
    tmpf3 = fab*fbs - fas*fbb;
    *determinant = fss*tmpf - fbs*tmpf2 + fbb*tmpf3;
  }
  else
    {printf("ERROR: Ncomponents shoud have value 2,3, or 4!\n"); return(0);}

  return(1);
}
