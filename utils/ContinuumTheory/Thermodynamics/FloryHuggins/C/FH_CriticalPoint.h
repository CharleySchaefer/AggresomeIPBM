

void FH_BinaryCriticalConcentration(FHPAR *FHPar, double *phiAcr)
{
  *phiAcr=1.0/( sqrt( (double)FHPar->NA/FHPar->NB ) + 1 );
}

void FH_BinaryCriticalChi(FHPAR *FHPar, double *chiABcr)
{
  double tmp=(1.0/sqrt(FHPar->NA) + 1.0/sqrt(FHPar->NB));
  (*chiABcr)=0.5*tmp*tmp;
}
