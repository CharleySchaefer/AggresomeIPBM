#include "GenerateEnergyLandscape.h" /* lattice structure*/

/*
  Gaussian Disorder Model:

    P(E) = exp(  -(E-Emean)^2/(2*sigma^2)  )/(2*pi*sigma^2)

  Cummulative P(E):

    p(E) = integral_{-inf}^{E}  P(e) de = ...

  Inverse:

    E(p) = Emean - sqrt(2)*sigma*erfinv(1-2*p)

*/

int GenerateEnergyLandscape(LATTICE_CUBE *Lattice,
 double Emean0, double sigma0, 
 double TrapFrac0, double Etrapmax0, double Etrapwidth0,
 double Emean1, double sigma1,
 double TrapFrac1, double Etrapmax1, double Etrapwidth1,
 double *EnergyLandscape)
{
  int i,j,k, siteID;
  double p;
  for(i=0; i<Lattice->Nx; i++) 
    for(j=0; j<Lattice->Ny; j++) 
      for(k=0; k<Lattice->Nz; k++) {

        siteID=coor2ind_cube(Lattice, i,j,k);
        /*CYTOSOL*/
        if(Lattice->site[siteID]<=-2 || Lattice->site[siteID]==0){
                   // Exponential Distribution
          if (( 1+ (double)rand() ) / ( (double)RAND_MAX+1.0)<TrapFrac0 ){
            p=( 1+ (double)rand() ) / ( (double)RAND_MAX+1.0);
            EnergyLandscape[siteID] = 
                Etrapmax0 + Etrapwidth0*log(p);
          } else { // Gaussian Distribution
            p=( 1+ (double)rand() ) / ( (double)RAND_MAX+1.0);
            EnergyLandscape[siteID] = 
                Emean0 - sqrt(2)*sigma0*(double)erfinv(1.0-2*p);
          }
        /* AGGRESOME */
        } else if(Lattice->site[siteID]>=2 || Lattice->site[siteID]==1) {    
                   // Exponential Distribution
          if (( 1+ (double)rand() ) / ( (double)RAND_MAX+1.0)<TrapFrac1 ){
            p=( 1+ (double)rand() ) / ( (double)RAND_MAX+1.0);
            EnergyLandscape[siteID] = 
                Etrapmax1 + Etrapwidth1*log(p);
          } else { // Gaussian Distribution
            p=( 1+ rand() ) / ((double)RAND_MAX+1.0); // rand (0,1]
            EnergyLandscape[siteID] =
                Emean1 - sqrt(2)*sigma1*(double)erfinv((1.0-2*p));
          }
        }
      }
  return(1);
}

int export_EnergyLandscape(char *outdir, char *fcnf,LATTICE_CUBE *Lattice, double *EnergyLandscape, char label) {
  int i,j,k,siteID;
  FILE *fp;
  sprintf(fcnf, "%s/EnergyLandscape%c.out", outdir, label);
  fp=fopen(fcnf, "w");


  for(k=0; k<Lattice->Nz; k++){
        for(i=0; i<Lattice->Nx; i++){
          for(j=0; j<Lattice->Ny; j++){
            siteID=coor2ind_cube(Lattice,i,j,k);
            fprintf(fp, " %9.2e", EnergyLandscape[siteID] );
          }
          fprintf(fp, "\n");
        }
  }

  fclose(fp);
  return(1);
}


