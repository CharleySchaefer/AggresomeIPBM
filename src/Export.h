#include "../lib/ZiltoidLIB/ZiltoidLIB.h"
/*
  CONTENTS:
  > export_settings
  > export_time_progress
  > export_particle_tracking
  > export_cnf

*/



int export_settings(char *fname, char *progname, double kmc_time, unsigned long long NTIME, 
int NX, int NY, int NZ, int rseed, 
int nA, double nu_A0, double nu_A1, double epsAA, double epsA1, 
int nB, double nu_B0, double nu_B1, double epsBB, double epsB1, double epsAB, 
int EnergyLandscapeModuleA, double EgaussmeanA0, double EgaussmeanA1, double EgausssigmaA0, double EgausssigmaA1, 
int EnergyLandscapeModuleB, double EgaussmeanB0, double EgaussmeanB1, double EgausssigmaB0, double EgausssigmaB1, 
double tbleach)
{
  FILE *fp=fopen(fname, "w");
  fprintf(fp, "Generated using %s %s\n",  progname, VERSION);
  /*------------------------------------------------*/
  /* INITILISATION                                  */
  fprintf(fp, "Simulation initialisation:\n");
  fprintf(fp, "  rseed   %d\n",   rseed);
  // TODO: print if continue-previous is selected
  fprintf(fp, "Simulation time:\n");
  fprintf(fp, "  time0   %e\n",   kmc_time);
  fprintf(fp, "  NTIME   %lld\n", NTIME);
  /*------------------------------------------------*/
  /* CELL TOPOLOGY                                  */
  fprintf(fp, "Simulation box:\n");
  fprintf(fp, "  NX      %d\n",   NX);
  fprintf(fp, "  NY      %d\n",   NY);
  fprintf(fp, "  NZ      %d\n",   NZ);
  // TODO: print name topology file
  /*------------------------------------------------*/
  /* PROTEIN PROPERTIES                             */
  fprintf(fp, "Properties protein A:\n");
  fprintf(fp, "  nA      %d\n",   nA);
  if(nA>0){
    fprintf(fp, "  nu_A0   %f\n",   nu_A0);
    fprintf(fp, "  nu_A1   %f\n",   nu_A1);
    fprintf(fp, "  epsAA   %f\n",   epsAA);
    /* interaction with regions of attraction */
    if(EnergyLandscapeModuleA){
      fprintf(fp, "  EgaussmeanA0   %f\n",   EgaussmeanA0);
      fprintf(fp, "  EgaussmeanA1   %f\n",   EgaussmeanA1);
      fprintf(fp, "  EgausssigmaA0  %f\n",   EgausssigmaA0);
      fprintf(fp, "  EgausssigmaA1  %f\n",   EgausssigmaA1);
// TODO: export info on A traps
    } else {
     fprintf(fp, "epsA1   %f\n",   epsA1);
    }
  }
  fprintf(fp, "Properties protein B:\n");
  fprintf(fp, "  nB      %d\n",   nB);
  if(nB>0){
    fprintf(fp, "  nu_B0   %f\n",   nu_B0);
    fprintf(fp, "  nu_B1   %f\n",   nu_B1);
    fprintf(fp, "  epsBB   %f\n",   epsBB);
    /* interaction with regions of attraction */
    if(EnergyLandscapeModuleB){
      fprintf(fp, "  EgaussmeanB0   %f\n",   EgaussmeanB0);
      fprintf(fp, "  EgaussmeanB1   %f\n",   EgaussmeanB1);
      fprintf(fp, "  EgausssigmaB0  %f\n",   EgausssigmaB0);
      fprintf(fp, "  EgausssigmaB1  %f\n",   EgausssigmaB1);
// TODO: export info on B traps
    } else {
     fprintf(fp, "epsB1   %f\n",   epsB1);
    }
  }
  if(nA>0 && nB>0)
    fprintf(fp, "  epsAB   %f\n",   epsAB);

  /*------------------------------------------------*/
  /* PHOTOBLEACHING                                 */
  fprintf(fp, "Photobleaching:\n");
  fprintf(fp, "  tbleach %e\n",   tbleach);
  fclose(fp);
  return(1);
}


/*
  input - mode=0: initialise file; mode=1: append data to file

*/
int export_time_progress(char *fname, unsigned long long iter, double kmc_time, unsigned long long count_accepted, int mode)
{
  FILE *fp=NULL;
  if(mode==0) {        // new file
    fp=fopen(fname, "w");
    fprintf(fp, "%12s %12s %12s\n", "kmc_steps", "time", "accepted(%)");
  } else if(mode==1) { // append data to file
    fp=fopen(fname, "a");
    fprintf(fp, "%12e %12e %12e\n", (double)iter, kmc_time, (iter==0? 0.0 : 100.0*count_accepted/iter));
  }
  if(fp!=NULL)
    fclose(fp);
  return(1);
}

/*
  Returns 0 if it failed to fopen file.
*/
int export_particle_tracking(char *fname, double kmc_time, int *x, int *y, int *z, int NA, int mode)
{
  int pid; /*particle ID*/
  FILE *fp=NULL;
  if (mode==0){ /* write file header*/
    fp = fopen(fname, "w");
    if (fp == NULL) {
      return(0);
    }

    // Write header  
//  fprintf(fp, "#NX\tNY\tNZ\n");
//  fprintf(fp, "%d\t%d\t%d\n#\n", NX, NY, NZ),
//  fprintf(fp, "#interactDeltaH: %3.1f\tParticle Number: %d\tNiter: %llu\n", interact_DeltaH,NA,NTIME);
    fprintf(fp, "#Time\t");
    for( pid = 0; pid < NA; pid++) {
      fprintf(fp, "x_part_%d\t", pid);
    }
    for( pid = 0; pid < NA; pid++) {
      fprintf(fp, "y_part_%d\t", pid);
    }
    for( pid = 0; pid < NA; pid++) {
      if(pid != NA - 1) {
        fprintf(fp, "z_part_%d\t", pid);
      } else {
        fprintf(fp, "z_part_%d\n", pid); // end of line
      }
    } 

  } else if(mode==1) { /* append data to file*/
    fp=fopen(fname, "a"); // append to file
    if (fp == NULL) {
      return(0);
    }

    // print time column
    fprintf(fp, "%f", kmc_time);
      
    // columns with x positions  
    for( pid = 0; pid<NA; pid++) {  
      fprintf(fp, "\t%d", x[pid]);
    }
    // columns with y positions  
    for( pid = 0; pid<NA; pid++) {
      fprintf(fp, "\t%d", y[pid]);
    }
    // columns with z positions  
    for( pid = 0; pid<NA; pid++) {
      fprintf(fp, "\t%d", z[pid]);
    }
    // end of line
    fprintf(fp, "\n");
  }
  if(fp!=NULL)
    fclose(fp);
  return(1);
}


/*
  label: 'A' or 'B'
*/
int export_cnf(char *outdir, char *fcnf, int count_cnf,LATTICE_CUBE *Lattice, int *bleached, char label) {
  int i,j,k,pid;
  FILE *fp;
  sprintf(fcnf, "%s/cnf/cnf%c%05d.out", outdir, label, count_cnf);

  fp=fopen(fcnf, "w");

  if( bleached==NULL ) {
    for(k=0; k<Lattice->Nz; k++){
        for(i=0; i<Lattice->Nx; i++){
          for(j=0; j<Lattice->Ny; j++){
            pid=get_particle_id(Lattice, coor2ind_cube(Lattice,i,j,k) );
            fprintf(fp, "%2d", (pid>=0?1:0) );
          }
          fprintf(fp, "\n");
        }
      }
  } else {
    for(k=0; k<Lattice->Nz; k++){
        for(i=0; i<Lattice->Nx; i++){
          for(j=0; j<Lattice->Ny; j++){
            pid=get_particle_id(Lattice, coor2ind_cube(Lattice,i,j,k) );
            fprintf(fp, "%2d", (pid>=0 ? 1+bleached[pid] : 0) );
          }
          fprintf(fp, "\n");
        }
    }
  }
  

  fclose(fp);
  return(1);
}


int export_bleach_profile(char *outdir, char *fcnf,LATTICE_CUBE *Lattice, int *bleach_profile) {
  int i,j,k,siteID;
  FILE *fp;
  sprintf(fcnf, "%s/BleachProfile.out", outdir);
  fp=fopen(fcnf, "w");

  for(k=0; k<Lattice->Nz; k++){
        for(i=0; i<Lattice->Nx; i++){
          for(j=0; j<Lattice->Ny; j++){
            siteID=coor2ind_cube(Lattice,i,j,k);
            fprintf(fp, "%2d", bleach_profile[siteID] );
          }
          fprintf(fp, "\n");
        }
  }

  fclose(fp);
  return(1);
}

