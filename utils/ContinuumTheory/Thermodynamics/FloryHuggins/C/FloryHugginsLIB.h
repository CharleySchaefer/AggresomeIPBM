#ifndef  FLORY_HUGGINS_H
  #define FLORY_HUGGINS_H 1
#include <math.h>


typedef struct {
  int    Ncomp;          // Number of components
  double NA, NB, NC, ND; // Molecular sizes
  double chiAB, chiAC, chiAD, chiBC, chiBD, chiCD;
  double phiA, phiB, phiC, phiD;
} FHPAR;

int is_incompressible(FHPAR *FHPar)
{
  int i;
  double phitot;
  switch (FHPar->Ncomp)
  {
    case 2:
        phitot=FHPar->phiA+FHPar->phiB;
      break;
    case 3:
        phitot=FHPar->phiA+FHPar->phiB+FHPar->phiC;
      break;
    case 4:
        phitot=FHPar->phiA+FHPar->phiB+FHPar->phiC+FHPar->phiD;
      break;
    default:
      printf("ERROR: Ncomp=%d should have value 2, 3, or 4.\n", FHPar->Ncomp); return(-1);
      break;
  }

  if(phitot==1)
    return(1);
  else
    return(0);
}

#include "lib/Polynomial/Polynomial.h" 
#include "FH_LocalFreeEnergy.h"
#include "FH_CriticalPoint.h"
#include "FH_Spinodal.h"
#include "FH_Binodal.h"

#endif
