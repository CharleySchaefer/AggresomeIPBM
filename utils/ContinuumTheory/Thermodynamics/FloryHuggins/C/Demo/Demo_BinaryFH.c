#include <stdio.h>
#include <stdlib.h>

#include "../FloryHugginsLIB.h"
//#include "FH_Binodal.h"

int main(int c, char *argv[])
{
  int i;
  int N=100, N_sp, N_bn;
  double chiABcr, phiAcr, f;
  double phi_sp1, phi_sp2, f_sp1, f_sp2;
  double phi_bn1, phi_bn2, f_bn1, f_bn2;
  double tol=1e-4;
 
  //--------------------------------------------------------
  // Define parameters
  FHPAR *FHPar=(FHPAR*)malloc(sizeof(FHPAR));
  FHPar->Ncomp=2; // Binary mixture
  FHPar->NA=1;
  FHPar->NB=1;
  FHPar->chiAB=0.2;
  FHPar->chiAB=(double)atof(argv[1]);
  FHPar->NA=(double)atof(argv[2]);
  printf("TEST: %e %e\n", FHPar->chiAB, FHPar->NB);
  //for (int i = 0; i < argc; ++i){
  //  FHPar->chiAB=atof(argv[i]);
  //}

  //--------------------------------------------------------
  // Determine Spinodal
  FH_SpinodalSolventFraction(FHPar, 0, 0, &N_sp, &phi_sp1, &phi_sp2);
  FHPar->phiA=phi_sp1;
  FH_LocalFreeEnergy(FHPar, &f_sp1);
  FHPar->phiA=phi_sp2;
  FH_LocalFreeEnergy(FHPar, &f_sp2);

  printf("#Spinodal:\n"); 
  printf("#phi_sp1: %f %f\n", phi_sp1, f_sp1); 
  printf("#phi_sp2: %f %f\n", phi_sp2, f_sp2); 


  //--------------------------------------------------------
  // Determine Binodal
  FH_BinaryBinodal(FHPar, tol, &N_bn, &phi_bn1, &phi_bn2);
  FHPar->phiA=phi_bn1;
  FH_LocalFreeEnergy(FHPar, &f_bn1);
  FHPar->phiA=phi_bn2;
  FH_LocalFreeEnergy(FHPar, &f_bn2);

  printf("#Binodal:\n"); 
  printf("#phi_bn1: %f %f\n", phi_bn1, f_bn1); 
  printf("#phi_bn2: %f %f\n", phi_bn2, f_bn2); 

  //--------------------------------------------------------
  // Determine critical point
  FH_BinaryCriticalConcentration( FHPar, &phiAcr  );
  FH_BinaryCriticalChi(           FHPar, &chiABcr   );

  printf("#Critical point:\n"); 
  printf("#phi_cr: %f\n", phiAcr); 
  printf("#chi_cr: %f\n", chiABcr); 
 

  //--------------------------------------------------------
  for (i=0; i<N; i++)
  {
    FHPar->phiA=(double)i/(N-1);
    FH_LocalFreeEnergy(FHPar, &f);
    printf("%8f %8f\n", FHPar->phiA, f);
  }

  return(0);
}
