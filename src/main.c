/*
AggresomeIBM - Individual Based Model for Protein Dynamics

Copyright (C) 2020 Charley Schaefer <The University of York, UK>

This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License
version 3 as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/
#include "main.h"
/* see info on program input & program functions in main.h */
/*
  kmc=kinetic monte carlo

  Cuboid geometry:    NX x NY x NZ
  Boundary condition: no-flux/non-permeable wall, no interaction with the wall

  Site values:
  -1: excluded volume (e.g., location of nucleoid)
   0: vacant, cytosol
   1: vacant, region of attraction (e.g., aggresome / stress granule)
 <-2: abs(value)-2 = particle ID of protein in cytosol         
  >2: value-2      = particle ID of protein in region of attraction


  ALGORITHM:
  > Initialisation configuration:
       Proteins are placed at random positions in the cell
  > For every lattice site, the number of occupied neighbours is determined
  > The enabled processes are listed
      definition process: the movement of a protein towards on of its six nearest neighbours.
      The maximum number of enabled processes (Nenabled) for proteins A is 6*NA, with NA the number of A proteins.
      The rate of each process may be different. The maximum rate by which A moves is rA_max
  > Time loop - at every time step / iteration:
         (for detailed explanation on the method, see, e.g., 
           APJ Jansen, An Introduction to Kinetic Monte Carlo Simulations of Surface Reactions, section 3.4)
         1. time is updated with step delta_t = -log(u)/(Nenabled*rA_max), with u a uniform random value on the interval (0,1]
         2. an enabled process is selected, and accepted/rejected with a certain probability.
            (time is updated regardless if the process is accepted or rejected)
         3. if the process is accepted, the protein is moved
         4. the number of neighbours per site is updated
         5. the list of enabled processes is updated
*/

int main(int argc, char *argv[])
{

  int i;
  char *arg;
  char attempt_who='A';
  /*==============================================*/
  /* USER SETTINGS                                */
  int continue_previous=0;                     /* see --continue-previous                                  */
  int track_A=0,  track_B=0;                    /* see --track-A --track-B                                  */
  int rseed=1; //time(0);                      /* see --rseed <value>                                      */
  /* Simulation volume  */
  int NX=50, NY=100, NZ=1;                  /* see --Nx <value>, --Ny ...                               */
  int xperiodic=0, yperiodic=0, zperiodic=0;   /* non-periodic boundaries (currently a fixed setting)      */
  /* Simulation time  */
  double kmc_time=0.0;                         /* initial time (overridden if --continue-previous is used) */
  unsigned long long NTIME      =5000,         /* see --Niter <value>                                      */
      NTIME_PRINT     =  50, iter0=0, iter;    /* see --Niter-print <value>                                */
  int NITER_PRINT_LOG=0;                       /* see --Niter-print-log                                    */

  /* Energy landscape for protein A */
  int    EnergyLandscapeModuleA=0;
  double EgaussmeanA0=0.0, EgaussmeanA1=0.0, EgausssigmaA0=0.0, EgausssigmaA1=0.0;
  double TrapFracA0=0; double EtrapmaxA0=0; double EtrapwidthA0=0;
  double TrapFracA1=0; double EtrapmaxA1=0; double EtrapwidthA1=0;
  /* Energy landscape for protein B */
  int    EnergyLandscapeModuleB=0;
  double EgaussmeanB0=0.0, EgaussmeanB1=0.0, EgausssigmaB0=0.0, EgausssigmaB1=0.0;
  double TrapFracB0=0; double EtrapmaxB0=0; double EtrapwidthB0=0;
  double TrapFracB1=0; double EtrapmaxB1=0; double EtrapwidthB1=0;

  /* Protein properties*/
  int    nA    = 1, nB=0;                     /* --NA ; --NB    - Number of A,B proteins                                */
  double epsAA = 0.0,epsBB = 0.0,epsAB = 0.0;  /* --epsAA ; etc. - protein-protein interaction energies in units of kT   */
  double epsA1 = 0.0, epsB1=0.0;                   /* --epsA1 ; etc. - Binding energies [kT] to regions of attraction        */
  double nu_A0=1.0, nu_B0=1.0;                     /* --nuA0 ; etc.  - mobilities in cytoplasm                               */
  double nu_A1=1.0, nu_B1=1.0;                     /* --nuA1 ; etc.  - mobilities in regions of attraction                   */
  double nu_AB=1.0, nu_BA=1.0;                     /* --nuAB , --nuBA */
  double MAB=1, MBA=1;                             /* Relative mobilities in mixed phase */

  /* Photobleaching */
  double tbleach=-1;                           /* time at which photobleaching is performed */
  char    bleach_who='A'; //default
  /* Input/output directory and files  */
  char fcnf[MAX_STR_L], ftimeprogress[MAX_STR_L], /* input and output files */
       ftrackA[MAX_STR_L], ftrackB[MAX_STR_L],    ftopology[MAX_STR_L],  
       fsettings[MAX_STR_L], fbleach[MAX_STR_L],
       outdir[MAX_STR_L-100];                         /* output director (default working directory '.')*/

  /* site and coordinates for the selected particle before and after the move*/
  int pID; /*particle identifier*/
  int siteID, x_old, y_old, z_old;
  int sitef, xf, yf, zf;
  
  double interact_probA[5], interact_probB[5];
  int move_direction;  // 0:x+1; 1:x-1; 2:y+1; 3:y-1; 4:z+1; 5:z-1 
  int move_tmp = -1;
  int NenabledA = 0,NenabledB=0;   /* Number of enabled processes */
  int num_processesA, num_processesB;
int bleach_iter=-1, bleach_iter_target=2, bleach_capture=0;


  /* Dynamics */
  double kmc_dt;
  double nu_Amax;//=(nu_A0>nu_A1?nu_A0:nu_A1); // Acceptance probability will be weighted using the maximum rate
  double nu_Bmax;//=(nu_B0>nu_B1?nu_B0:nu_B1); // Acceptance probability will be weighted using the maximum rate
  double rsm_kA, rsm_kB, rsm_k,Delta_t; /* sum of all rates sets step size according to Random Selection Method */
  long long unsigned count_accepted=0;
  int output_counter = 0;

  /* sites & coordinates for all particles */
  int *XA0 = NULL, *YA0 = NULL,*ZA0 = NULL;           /* initial position holders */
  int *XB0 = NULL, *YB0 = NULL,*ZB0 = NULL;           
  int *siteA0=NULL, *xA=NULL, *yA=NULL, *zA=NULL;     /* current position holders */
  int *siteB0=NULL, *xB=NULL, *yB=NULL, *zB=NULL;
  int *xA_cross=NULL, *yA_cross=NULL, *zA_cross=NULL; /* track if protein crossed the periodic boundaries (if activated) */
  int *xB_cross=NULL, *yB_cross=NULL, *zB_cross=NULL; /* track if protein crossed the periodic boundaries (if activated) */
  int x_cr=0, y_cr=0, z_cr=0; /* integers for possible cross in current move. -1: downwards, 0: none, 1: upwards */
  int *enabledA=NULL,*enabledB=NULL, *process_regA=NULL,  *process_regB=NULL;
  int Nexcluded=0;

  LATTICE_CUBE *LatticeA = NULL, *LatticeB = NULL;  
  int    *ibuff=NULL;
  double *dbuff=NULL;
  double *EnergyLandscapeA=NULL,*EnergyLandscapeB=NULL;
  int    *bleached=NULL, *bleach_profile=NULL;

  sprintf(ftopology, "-1");  /* optional input file */
  sprintf(fbleach,   "-1") ; /* optional input file */
  sprintf(outdir, ".");      /* default output directory */
  /* DECLARATIONS DONE */
  /*==============================================*/


  /*==============================================*/
  /* READ INPUT */
  if(argc==1) {
    printf("No arguments given - exiting.\n");
    printf("For more information, run \'%s --help\'\n", argv[0]);
  }
  i=1;
  while( i < argc ){
    arg = argv[i];
    if ( strcmp( arg, "--help" ) == 0 ){       /* print help */
      printHelp(argv[0], rseed, NTIME, NTIME_PRINT, NX, NY, NZ);
      exit(1);
    } else if ( strcmp( arg, "--version" ) == 0 ){ /* Random number seed */
      printf("%s %s\n", argv[0]+2, VERSION);  exit(1);
      i += 1;
    /* SIMULATION: INITIALISATION, TIME STEPS, AND OUTPUT */
    } else if ( strcmp( arg, "--rseed" ) == 0 ){ /* Random number seed */
      rseed=atoi(argv[i+1]);  
      i += 2;
    } else if ( strcmp( arg, "--continue-previous" ) == 0 ) { 
      continue_previous=1; i++;
    } else if ( strcmp( arg, "--Niter" ) == 0 ){
      NTIME=atoll(argv[i+1]);  
      i += 2;
    } else if ( strcmp( arg, "--Niter-print" ) == 0 ){
      NTIME_PRINT=atoll(argv[i+1]);  
      i += 2;
    } else if ( strcmp( arg, "--Niter-print-log" ) == 0 ){
      NITER_PRINT_LOG=1; NTIME_PRINT=1; /* double NTIME_PRINT at every print step*/
      i += 1;
    } else if ( strcmp( arg, "--outdir" ) == 0 ){ /* dir */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); exit(1);
      } else {
        sprintf(outdir, "%s", argv[i+1]); i+=2; 
        printf("  Output directory:          %s\n", outdir);
      }
    } else if ( strcmp( arg, "--track-A" ) == 0 ) { 
      track_A=1; i++;
    } else if ( strcmp( arg, "--track-B" ) == 0 ) { 
      track_B=1; i++;
    /* SIMULATION VOLUME & SET CELL TOPOLOGY */
    } else if ( strcmp( arg, "--Nx" ) == 0 ){
      NX=atoi(argv[i+1]);  
      i += 2;
    } else if ( strcmp( arg, "--Ny" ) == 0 ){
      NY=atoi(argv[i+1]);  
      i += 2;
    } else if ( strcmp( arg, "--Nz" ) == 0 ){
      NZ=atoi(argv[i+1]);  
      i += 2;
    } else if ( strcmp( arg, "--periodic-boundary" ) == 0 ){
      xperiodic=1; yperiodic=1; zperiodic=1;  
      i += 1;
    } else if ( strcmp( arg, "--E-gauss-mean-A-0" ) == 0 ){  
      EgaussmeanA0=atof(argv[i+1]);  EnergyLandscapeModuleA=1;
      i += 2;
    } else if ( strcmp( arg, "--E-gauss-mean-A-1" ) == 0 ){ 
      EgaussmeanA1=atof(argv[i+1]);   EnergyLandscapeModuleA=1;
      i += 2;
    } else if ( strcmp( arg, "--E-gauss-sigma-A-0" ) == 0 ){ 
      EgausssigmaA0=atof(argv[i+1]);   EnergyLandscapeModuleA=1;
      i += 2;
    } else if ( strcmp( arg, "--E-gauss-sigma-A-1" ) == 0 ){ 
      EgausssigmaA1=atof(argv[i+1]);   EnergyLandscapeModuleA=1;
      i += 2;
    } else if ( strcmp( arg, "--E-trap-frac-A-0" ) == 0 ){ 
      TrapFracA0=atof(argv[i+1]);   EnergyLandscapeModuleA=1;
      i += 2;
    } else if ( strcmp( arg, "--E-trap-max-A-0" ) == 0 ){ 
      EtrapmaxA0=atof(argv[i+1]);   EnergyLandscapeModuleA=1;
      i += 2;
    } else if ( strcmp( arg, "--E-trap-width-A-0" ) == 0 ){ 
      EtrapwidthA0=atof(argv[i+1]);   EnergyLandscapeModuleA=1;
      i += 2;
    } else if ( strcmp( arg, "--E-trap-frac-A-1" ) == 0 ){ 
      TrapFracA1=atof(argv[i+1]);   EnergyLandscapeModuleA=1;
      i += 2;
    } else if ( strcmp( arg, "--E-trap-max-A-1" ) == 0 ){ 
      EtrapmaxA1=atof(argv[i+1]);   EnergyLandscapeModuleA=1;
      i += 2;
    } else if ( strcmp( arg, "--E-trap-width-A-1" ) == 0 ){ 
      EtrapwidthA1=atof(argv[i+1]);   EnergyLandscapeModuleA=1;
      i += 2;
    } else if ( strcmp( arg, "--E-gauss-mean-B-0" ) == 0 ){  
      EgaussmeanB0=atof(argv[i+1]);  EnergyLandscapeModuleB=1;
      i += 2;
    } else if ( strcmp( arg, "--E-gauss-mean-B-1" ) == 0 ){ 
      EgaussmeanB1=atof(argv[i+1]);   EnergyLandscapeModuleB=1;
      i += 2;
    } else if ( strcmp( arg, "--E-gauss-sigma-B-0" ) == 0 ){ 
      EgausssigmaB0=atof(argv[i+1]);   EnergyLandscapeModuleB=1;
      i += 2;
    } else if ( strcmp( arg, "--E-gauss-sigma-B-1" ) == 0 ){ 
      EgausssigmaB1=atof(argv[i+1]);   EnergyLandscapeModuleB=1;
      i += 2;
    } else if ( strcmp( arg, "--E-trap-frac-B-0" ) == 0 ){ 
      TrapFracB0=atof(argv[i+1]);   EnergyLandscapeModuleB=1;
      i += 2;
    } else if ( strcmp( arg, "--E-trap-max-B-0" ) == 0 ){ 
      EtrapmaxB0=atof(argv[i+1]);   EnergyLandscapeModuleB=1;
      i += 2;
    } else if ( strcmp( arg, "--E-trap-width-B-0" ) == 0 ){ 
      EtrapwidthB0=atof(argv[i+1]);   EnergyLandscapeModuleB=1;
      i += 2;
    } else if ( strcmp( arg, "--E-trap-frac-B-1" ) == 0 ){ 
      TrapFracB1=atof(argv[i+1]);   EnergyLandscapeModuleB=1;
      i += 2;
    } else if ( strcmp( arg, "--E-trap-max-B-1" ) == 0 ){ 
      EtrapmaxB1=atof(argv[i+1]);   EnergyLandscapeModuleB=1;
      i += 2;
    } else if ( strcmp( arg, "--E-trap-width-B-1" ) == 0 ){ 
      EtrapwidthB1=atof(argv[i+1]);   EnergyLandscapeModuleB=1;
      i += 2;
    } else if ( strcmp( arg, "--import-topology" ) == 0 ){ /* dir */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); exit(1);
      } else {
        sprintf(ftopology, "%s", argv[i+1]); i+=2; 
        printf("  Importing topology file:          %s\n", ftopology);
      }
    /* PHOTOBLEACHING */
    } else if ( strcmp( arg, "--bleach-time" ) == 0 ){
      tbleach = atof(argv[i+1]);
      i+=2;
    } else if ( strcmp( arg, "--bleach-profile" ) == 0 ){ /* dir */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); exit(1);
      } else {
        sprintf(fbleach, "%s", argv[i+1]); i+=2; 
        printf("  Importing bleaching profile:          %s\n", fbleach);
      }
    /* PROTEINS AND THEIR PHYSICAL PROPERTIES */
         /* PROTEIN A */
    } else if ( strcmp( arg, "--NA" ) == 0 ){
      nA = atoi(argv[i+1]); 
      i+=2;
    } else if ( strcmp( arg, "--epsAA" ) == 0 ){
      epsAA = atof(argv[i+1]); 
      i+=2;
    } else if ( strcmp( arg, "--epsA1" ) == 0 ){
      epsA1 = atof(argv[i+1]); 
      i+=2;
    } else if ( strcmp( arg, "--nuA0" ) == 0 ){
      nu_A0 = atof(argv[i+1]); 
      i+=2;
    } else if ( strcmp( arg, "--nuA1" ) == 0 ){
      nu_A1 = atof(argv[i+1]); 
      i+=2;
    } else if ( strcmp( arg, "--nuAB" ) == 0 ){
      nu_AB = atof(argv[i+1]); 
      i+=2;
         /* PROTEIN B */
    } else if ( strcmp( arg, "--NB" ) == 0 ){
      nB = atoi(argv[i+1]); 
      i+=2;
    } else if ( strcmp( arg, "--epsAB" ) == 0 ){
      epsAB = atof(argv[i+1]); 
      i+=2;
    } else if ( strcmp( arg, "--epsBB" ) == 0 ){
      epsBB = atof(argv[i+1]); 
      i+=2;
    } else if ( strcmp( arg, "--epsB1" ) == 0 ){
      epsB1 = atof(argv[i+1]); 
      i+=2;
    } else if ( strcmp( arg, "--nuB0" ) == 0 ){
      nu_B0 = atof(argv[i+1]); 
      i+=2;
    } else if ( strcmp( arg, "--nuB1" ) == 0 ){
      nu_B1 = atof(argv[i+1]); 
      i+=2;
    } else if ( strcmp( arg, "--nuBA" ) == 0 ){
      nu_BA = atof(argv[i+1]); 
      i+=2;
    /* ERROR CAPTURE */
    } else {
      printf( "Error: Argument \"%s\" not recognized!\n", arg );
      printf( "       For help, run \"%s --help\".\n", argv[0] );
      
      exit(1);
    }
  }

   nu_Amax=(nu_A0>nu_A1?nu_A0:nu_A1); // Acceptance probability will be weighted using the maximum rate
   nu_Bmax=(nu_B0>nu_B1?nu_B0:nu_B1); // Acceptance probability will be weighted using the maximum rate

  if(nA>0 && nB>0){
    MAB=(nu_AB/nu_A0);
    if(MAB>1){
      printf("Error: in the current implementation, the mobility nu_AB should be smaller than or equal to nu_A0\n");
    }
    MBA=(nu_BA/nu_B0);
    if(MBA>1){
      printf("Error: in the current implementation, the mobility nu_BA should be smaller than or equal to nu_B0\n");
    }
  }

  /* output directory 'outdir' should exist before running the program */
  DIR* dir = opendir(outdir);
  if (dir) {  /* Directory exists. */
    closedir(dir);
  } else if (ENOENT == errno) {
    printf("Error: directory %s does not exist - trying to create using 'mkdir'.\n", outdir); exit(1);
    sprintf(fcnf, "mkdir %s", outdir );
    if(system(fcnf)==-1) {
      printf("Failed!\n"); goto TERMINATE_PROGRAM;
    }
    closedir(dir);
  } else {
    printf("Error: opendir(%s) failed.\n", outdir); exit(1);
  }
  /* output directory 'outdir' should exist before running the program */
  sprintf(fcnf, "%s/cnf", outdir);
  dir = opendir(fcnf);
  if (dir) {  /* Directory exists. */
    closedir(dir);
  } else if (ENOENT == errno) {
    printf("Warning: directory %s/cnf does not exist - trying to create using 'mkdir'.\n", outdir); 
    sprintf(fcnf, "mkdir %s/cnf", outdir );
    if(system(fcnf)==-1) {
      printf("Failed!\n"); goto TERMINATE_PROGRAM;
    }
  } else {
    printf("Error: opendir(%s) failed.\n", outdir); exit(1);
  }
  
  printf("randmax: %d\n", RAND_MAX); ;//exit(1);  
  
   /*will only be used if --track-A / --track-B is activated */
  sprintf(ftrackA,        "%s/ParticleTrackingA.out", outdir);
  sprintf(ftrackB,        "%s/ParticleTrackingB.out", outdir);


  sprintf(ftimeprogress, "%s/timeprogress.out", outdir);
  sprintf(fsettings,     "%s/settings.out", outdir);
  export_settings(fsettings, argv[0], kmc_time, NTIME, NX, NY, NZ, rseed, 
    nA, nu_A0, nu_A1, epsAA, epsA1, nB, nu_B0, nu_B1, epsBB, epsB1, epsAB, 
    EnergyLandscapeModuleA, EgaussmeanA0, EgaussmeanA1, EgausssigmaA0, EgausssigmaA1, 
    EnergyLandscapeModuleB, EgaussmeanB0, EgaussmeanB1, EgausssigmaB0, EgausssigmaB1, 
    tbleach);
/*(fsettings, argv[0], kmc_time, NTIME, NX, NY, NZ, rseed, 
                   nA, epsAA, epsA1, nu_A0, nu_A1, 
                   nB, epsBB, epsB1, nu_B0, nu_B1, epsAB, 
                   EnergyLandscapeModule, EgaussmeanA0, EgaussmeanA1, EgausssigmaA0, EgausssigmaA1, 
                   tbleach);*/

  /* END USER SETTINGS                            */
  /*==============================================*/


  /*==============================================*/
  /* INITIALISE */



  if(NY>0) {ibuff=(int*)malloc(NY*sizeof(int)); }
   /* memory buffer; should be at least length NY to read CellTopology and BleachProfile files*/
 
  if(MAX_LINE_WIDTH>0){dbuff=(double*)malloc(MAX_LINE_WIDTH*sizeof(double));} /* memory buffer; should be at least length NY to read CellTopology and BleachProfile files*/

  iter0=0;
  kmc_time=0.0;
  if(continue_previous) { /* continue_previous: cnf file will be read later */
    readTimeProgress(ftimeprogress, &iter0, &kmc_time, &count_accepted, &output_counter, dbuff);
  }

  // Seed random number generator 
  sfmt_init_gen_rand(&sfmt, rseed);
  srand(rseed);

  /* initialise lattice */
 // if(nA>0){
    LatticeA=make_lattice_cube(NX, NY, NZ);
    LatticeA->is_periodic_x=xperiodic;
    LatticeA->is_periodic_y=yperiodic;
    LatticeA->is_periodic_z=zperiodic;
 // }
 // if(nB>0){
    LatticeB=make_lattice_cube(NX, NY, NZ);
    LatticeB->is_periodic_x=xperiodic;
    LatticeB->is_periodic_y=yperiodic;
    LatticeB->is_periodic_z=zperiodic;
 // }
  
  /* import cell topology (optional) - regions of attraction (aggresome/stress granule/microtubuli) and excluded volume (nucleoid/cell porosity) */
  if(strcmp(ftopology, "-1")) {
    printf("Importing topology file\n");
    // TODO: use assert to check file existence
    if( nA>0 && !initialise_cell_topology(ftopology, LatticeA, ibuff)){
      printf("Error: initialise_cell_topology() failed.\n" ); exit(1);
    }
    if( nB>0 &&!initialise_cell_topology(ftopology, LatticeB, ibuff)){
      printf("Error: initialise_cell_topology() failed.\n" ); exit(1);
    }
  }
  Nexcluded=0; /* number of excluded volume sites*/
  for (siteID=0; siteID<LatticeA->Nxyz; siteID++){
    if(nA>0 && LatticeA->site[siteID]==-1)
      Nexcluded++;
  }
  if(nA>0 && nA>LatticeA->Nxyz-Nexcluded) {
    printf("Error: number of A proteins exceeds number of available lattice sites.\n"); goto TERMINATE_PROGRAM;
  }
  if(nB>0 && nB>LatticeB->Nxyz-Nexcluded) {
    printf("Error: number of B proteins exceeds number of available lattice sites.\n"); goto TERMINATE_PROGRAM;
  }

  /* Energy Landscape */
  if(EnergyLandscapeModuleA) {
      printf("Creating Disordered Energy Landscape A\n");
     EnergyLandscapeA=(double *)malloc(LatticeA->Nxyz*sizeof(double));
    GenerateEnergyLandscape(LatticeA, EgaussmeanA0, EgausssigmaA0,
TrapFracA0, EtrapmaxA0, EtrapwidthA0, EgaussmeanA1, EgausssigmaA1,
TrapFracA1, EtrapmaxA1, EtrapwidthA1, EnergyLandscapeA);

    export_EnergyLandscape(outdir, fcnf,LatticeA, EnergyLandscapeA, 'A');
  }
  if(EnergyLandscapeModuleB) {
      printf("Creating Disordered Energy Landscape B\n");
     EnergyLandscapeB=(double *)malloc(LatticeA->Nxyz*sizeof(double));
    GenerateEnergyLandscape(LatticeA, EgaussmeanB0, EgausssigmaB0,
TrapFracB0, EtrapmaxB0, EtrapwidthB0, EgaussmeanB1, EgausssigmaB1,
TrapFracB1, EtrapmaxB1, EtrapwidthB1, EnergyLandscapeB);

    export_EnergyLandscape(outdir, fcnf,LatticeA, EnergyLandscapeB, 'B');
  }




  if (nA>0){
    XA0 = (int*)malloc(sizeof(*XA0) * nA);         
    YA0 = (int*)malloc(sizeof(*XA0) * nA);
    ZA0 = (int*)malloc(sizeof(*XA0) * nA);
    siteA0 = (int*)malloc(sizeof(*siteA0) * nA);
    xA = (int*)malloc(sizeof(*xA) * nA);
    yA = (int*)malloc(sizeof(*yA) * nA);
    zA = (int*)malloc(sizeof(*zA) * nA);
    xA_cross = (int*)calloc(nA, sizeof(*xA_cross));
    yA_cross = (int*)calloc(nA, sizeof(*yA_cross));
    zA_cross = (int*)calloc(nA, sizeof(*zA_cross));
  }
  if(nB>0){
    XB0 = (int*)malloc(sizeof(*XB0) * nB);         
    YB0 = (int*)malloc(sizeof(*XB0) * nB);
    ZB0 = (int*)malloc(sizeof(*XB0) * nB);
    siteB0 = (int*)malloc(sizeof(*siteB0) * nB);
    xB = (int*)malloc(sizeof(*xB) * nB);
    yB = (int*)malloc(sizeof(*yB) * nB);
    zB = (int*)malloc(sizeof(*zB) * nB);
    xB_cross = (int*)calloc(nB, sizeof(*xB_cross));
    yB_cross = (int*)calloc(nB, sizeof(*yB_cross));
    zB_cross = (int*)calloc(nB, sizeof(*zB_cross));
  }  
         



  /* Photobleaching */ // TODO: chose if A or B is bleached
  if(tbleach>=0){
    bleach_profile=(int*)calloc(LatticeA->Nxyz, sizeof(int));
    bleached=(int*)calloc((nA>nB?nA:nB), sizeof(int)); /* initialise with zero values */
    if(strcmp(fbleach, "-1")) {
      printf("Importing bleaching profile input file\n");
      /* LatticeA: only Nx/Ny/Nz used*/
      if(!initialise_bleach_profile(fbleach, LatticeA, ibuff, bleach_profile)){ 
        printf("Error: initialise_bleach_profile() failed.\n" ); exit(1);
      }
      // export
      /* LatticeA: only Nx/Ny/Nz used*/
      export_bleach_profile(outdir, fcnf,LatticeA, bleach_profile);
    } else {
      printf("Error: no bleaching profile imported; use --bleach-profile <fname> argument.\n");
      exit(1);
    }
  }


  /* Initialise particle positions */
  if(continue_previous){
    if( (nA>0) && !import_previous_cnf(LatticeA, ftimeprogress, fcnf, outdir,  siteA0, XA0, YA0, ZA0, ibuff, 'A')){
       printf("Error: import_previous_cnf() failed.\n"); exit(1);
    }
    if( (nB>0) && !import_previous_cnf(LatticeB, ftimeprogress, fcnf, outdir,  siteB0, XB0, YB0, ZB0, ibuff, 'B')){
       printf("Error: import_previous_cnf() failed.\n"); exit(1);
    }
  } else {

  /* Initialise particle positions*/
  for(int p = 0; p < nA; p++) // A-A EXCLUDED VOLUME INTERACTION -> ONLY ONE A PER SITE
  {
      /* pick random site outside nucleus and not on another particle */
      siteA0[p]=(int)(LatticeA->Nxyz-1)*((double)rand()/(double)RAND_MAX);
      while( LatticeA->site[siteA0[p]] <=-1 || LatticeA->site[siteA0[p]] >= 2) {
        siteA0[p]=(int)(LatticeA->Nxyz-1)*((double)rand()/(double)RAND_MAX); 
      }

      /* update LatticeA -> There's now the particle with ID p in the previously empty site */
      if(LatticeA->site[siteA0[p]] == 0) {LatticeA->site[siteA0[p]] = -2 - p;} // cytoplasm
      else if(LatticeA->site[siteA0[p]] == 1) {LatticeA->site[siteA0[p]] = 2 + p;} // aggresome
      else {printf("Error with assigning the particle ID.");}
      /* define initial positions. X0+p is pointer arithmetic! */
      ind2coor_cube(LatticeA, siteA0[p], XA0+p, YA0+p, ZA0+p);
    
  }
  // TODO: Add option to have an A-B excluded volume interaction
  for(int p = 0; p < nB; p++)  // B-B EXCLUDED VOLUME INTERACTION -> ONLY ONE B PER SITE
  {                            // HOWEVER, A AND B CAN OCCUPY THE SAME SITE!
      /* pick random site outside nucleus and not on another particle */
      siteB0[p]=(int)(LatticeB->Nxyz-1)*((double)rand()/(double)RAND_MAX);
      while( LatticeB->site[siteB0[p]] <=-1 || LatticeB->site[siteB0[p]] >= 2) {
        siteB0[p]=(int)(LatticeB->Nxyz-1)*((double)rand()/(double)RAND_MAX); 
      }

      /* update LatticeB -> There's now the particle with ID p in the previously empty site */
      if(LatticeB->site[siteB0[p]] == 0) {LatticeB->site[siteB0[p]] = -2 - p;} // cytoplasm
      else if(LatticeB->site[siteB0[p]] == 1) {LatticeB->site[siteB0[p]] = 2 + p;} // aggresome
      else {printf("Error with assigning the particle ID.");}
      /* define initial positions. X0+p is pointer arithmetic! */
      ind2coor_cube(LatticeB, siteB0[p], XB0+p, YB0+p, ZB0+p);
    
  } // end for loop
  } // end initialise particles




  // Calculate initial number of nearest neighbours at each lattice point
  if (nA>0) {
    calc_occupied_neighbours(LatticeA);
    /* bare acceptance probability just affected by A-A interactions */
    for(int neigh_d = 1; neigh_d <= 5; neigh_d++) {
      interact_probA[neigh_d - 1] = exp(-neigh_d*epsAA);
    }
  }
  if (nB>0) {
    calc_occupied_neighbours(LatticeB);
    /* bare acceptance probability just affected by B-B interactions */
    for(int neigh_d = 1; neigh_d <= 5; neigh_d++) {
      interact_probB[neigh_d - 1] = exp(-neigh_d*epsBB);
    }
  }

  /* Current Position = Starting Position. Memcpy faster than for loop as I understood */
  if(nA>0){
    memcpy(xA, XA0, sizeof(*XA0) * nA);
    memcpy(yA, YA0, sizeof(*YA0) * nA);
    memcpy(zA, ZA0, sizeof(*ZA0) * nA);
  }
  if(nB>0){
    memcpy(xB, XB0, sizeof(*XB0) * nB);
    memcpy(yB, YB0, sizeof(*YB0) * nB);
    memcpy(zB, ZB0, sizeof(*ZB0) * nB);
  }
  

  /* Tracking Enabled Processes */

  num_processesA = HOP_DIRECTIONS*nA;
  num_processesB = HOP_DIRECTIONS*nB;

    // memory allocation
  if(nA>0) {
    enabledA = (int*)malloc(sizeof(*enabledA) * num_processesA);           /* holds all enabled moves*/
    process_regA = (int*)malloc(sizeof(*process_regA) * num_processesA);   /* index: move x ; process_reg[index] : position of move x in "enabled" */
    for(int l = 0; l<num_processesA; l++) {
      enabledA[l]=-1;
      process_regA[l]=-1;
    }
  }
  if(nB>0){
    enabledB = (int*)malloc(sizeof(*enabledB) * num_processesB);           /* holds all enabled moves*/
    process_regB = (int*)malloc(sizeof(*process_regB) * num_processesB);   /* index: move x ; process_reg[index] : position of move x in "enabled" */
    for(int l = 0; l<num_processesB; l++) {
      enabledB[l]=-1;
      process_regB[l]=-1;
    }
  }
           
    // initialise enabled processes                                                                   /* number of possible moves */
    /* PROTEIN A */
  for (int process_ID_tmp = 0, temp_part_id = 0, allow = 0; 
        process_ID_tmp < num_processesA; process_ID_tmp++)                                         /* evaluate possible initial moves */
  {
    temp_part_id = (process_ID_tmp) / HOP_DIRECTIONS;
    move_tmp =  process_ID_tmp % HOP_DIRECTIONS;
    allow = possible_y1n0(LatticeA, xA[temp_part_id], yA[temp_part_id], zA[temp_part_id], move_tmp);
    if(allow == 1) {
      if(add_enabled(enabledA, &NenabledA, process_regA, process_ID_tmp) == 0)                       /* Add process_ID_tmp to "enabled array"; Also updates Nenabled */
      {
        printf("Error when initialising the enabled processes"
        "of the just enabled particles. Process was already enabled."
        "That cannot be possible as this is the first try to enable this process.");
      }
    }
    else if(allow==-1) {
      printf("Exiting. Error with possible_y1n0 function.\n"); exit(EXIT_FAILURE);  
    }
  }
    /* PROTEIN B */
  for (int process_ID_tmp = 0, temp_part_id = 0, allow = 0; 
        process_ID_tmp < num_processesB; process_ID_tmp++)                                         /* evaluate possible initial moves */
  {
    temp_part_id = (process_ID_tmp) / HOP_DIRECTIONS;
    move_tmp =  process_ID_tmp % HOP_DIRECTIONS;
    allow = possible_y1n0(LatticeB, xB[temp_part_id], yB[temp_part_id], zB[temp_part_id], move_tmp);
    if(allow == 1) {
      if(add_enabled(enabledB, &NenabledB, process_regB, process_ID_tmp) == 0)                       /* Add process_ID_tmp to "enabled array"; Also updates Nenabled */
      {
        printf("Error when initialising the enabled processes"
        "of the just enabled particles. Process was already enabled."
        "That cannot be possible as this is the first try to enable this process.");
      }
    }
    else if(allow==-1) {
      printf("Exiting. Error with possible_y1n0 function.\n"); exit(EXIT_FAILURE);  
    }
  }



  /* WRITE OUTPUT FILES */
  if(track_A) { /* write header to file */
    export_particle_tracking(ftrackA, kmc_time, xA, yA, zA, nA, 0);
  }
  if(track_B) { /* write header to file */
    export_particle_tracking(ftrackB, kmc_time, xB, yB, zB, nB, 0);
  }
  if(!continue_previous)
    export_time_progress(ftimeprogress, iter0, kmc_time, count_accepted, 0); /*initialise timeprogress file*/
  /* INITIALISATION & PREPARATIONS DONE */
  /*==============================================*/



  /*==============================================*/
  // START MONTE CARLO TIME STEPPER
  //###########################################################################################
  // ONLY PROTEIN A
  if(  (nA>0)&&(nB==0)  ) {
  for(iter=iter0; iter<=NTIME+iter0; iter++) 
  {
    /*--------------------------------------------------*/
    /* PHOTOBLEACH                                      */
    if(tbleach>=0 && bleach_iter==-1 && kmc_time>=tbleach ) { 
      for(siteID=0; siteID<LatticeA->Nxyz; siteID++) {
        pID=get_particle_id(LatticeA, siteID);
        /*bleach_profile[] gives the probability that the protein is bleached */
        if( pID>=0 && ( 1.0 + rand() ) / ((double)RAND_MAX+1.0) < bleach_profile[siteID] ) { 
          bleached[pID]=1;  // Protein label is set from 0 to 1
        }
      }
      bleach_iter=iter;     // Remember iteration at which was bleached
    }
    /*--------------------------------------------------*/

    /*--------------------------------------------------*/
    /* EXPORT DATA */
    bleach_capture=( bleach_iter!=-1 && ((iter-bleach_iter)%bleach_iter_target)==0);
    if( iter%NTIME_PRINT == 0 || iter==NTIME+iter0 || bleach_capture)
    {
      /* logarithmic time interval */
      if(NITER_PRINT_LOG)
        NTIME_PRINT*=2;
      if(bleach_capture)
        bleach_iter_target*=2;


      output_counter++;
      /* Export time progress */
      printf("iter: %12lld (%5.2f%%); time: %12e; accepted %5.2f%%\n", iter, 100.0*(iter-iter0)/NTIME, kmc_time, (iter==0?0:100.0*(double)count_accepted/iter));
      export_time_progress(ftimeprogress, iter, kmc_time, count_accepted, 1);
      /* Export configuration */
      export_cnf(outdir, fcnf, output_counter, LatticeA, bleached, 'A');

      /* Export tracked particles */
      if(track_A)
        export_particle_tracking(ftrackA, kmc_time, xA, yA, zA, nA, 1);
      if(track_B)
        export_particle_tracking(ftrackB, kmc_time, xB, yB, zB, nB, 1);
    }
    /*--------------------------------------------------*/


    /*---------------*/
    /* UPDATE TIME: Random selection -> time updated regardless if process is accepted or rejected */
    rsm_k=1.0/(NenabledA*nu_Amax); /* maximum time step size */
    kmc_dt = -rsm_k*log( ( 1.0 + rand() ) / ((double)RAND_MAX+1.0) );
    kmc_time+=kmc_dt;

    /*---------------*/
    /* SELECT EVENT  */
    SelectEnabledProcess(LatticeA, NenabledA, enabledA, process_regA, siteA0, xA,yA,zA, &pID, &move_direction, &siteID, &x_old,&y_old,&z_old, &sitef, &xf,&yf,&zf);

    /*---------------*/
    // ACCEPT/REJECT MOVE
    if(epsAA==0) { /* No A-A interactions*/
      if(EnergyLandscapeModuleA) { /* Disordered energy landscape */
        if(0==AcceptRejectS0NN_EL(LatticeA, siteID, sitef, EnergyLandscapeA, nu_Amax, nu_A0, nu_A1))
          continue; /* REJECTED --> End of current time iteraction */
      } else {
        if(0==AcceptRejectS0NN(LatticeA, siteID, sitef, epsA1,           nu_Amax, nu_A0, nu_A1))
          continue; /* REJECTED --> End of current time iteraction */ 
      }
    } else {
      if(EnergyLandscapeModuleA) {
         printf("Error: EnergyLandscapeModuleA and nearest neighbour interactions not yet implemented.\n");
      }else {
        if(0==AcceptRejectS6NN(LatticeA, siteID, sitef, epsAA, epsA1, nu_Amax, nu_A0, nu_A1,interact_probA))
        continue; /* REJECTED --> End of current time iteraction */
      }
    }
    // else: ACCEPTED
    count_accepted++;

    //---------------------
    /* EXECUTE MOVE AND UPDATE LATTICE & PARTICLE POS */
    // update configuration
    update_cnf(LatticeA, siteA0, xA, yA, zA, pID, siteID, sitef, xf, yf, zf);

    /* update enabled processes and update the number of occupied neighbours */
    update_enabled_and_neighbours(LatticeA, &NenabledA, enabledA, process_regA, pID, siteID, x_old, y_old, z_old,sitef, xf,yf,zf, move_direction);
    
    // update counter (how many times did protein cross a boundary in case of periodic boundaries)
    if(x_cr!=0)      {xA_cross[pID]+=x_cr;} 
    else if(y_cr!=0) {yA_cross[pID]+=y_cr;} 
    else if(z_cr!=0) {zA_cross[pID]+=z_cr;}
    x_cr=0;y_cr=0;z_cr=0; // reset
    //---------------------
  } // last time step
  //###########################################################################################
  // ONLY PROTEIN B
  } else if( (nA==0) && nB>0) {
  for(iter=iter0; iter<=NTIME+iter0; iter++) 
  {
    /*--------------------------------------------------*/
    /* PHOTOBLEACH                                      */
    if(tbleach>=0 && bleach_iter==-1 && kmc_time>=tbleach ) { 
      for(siteID=0; siteID<LatticeB->Nxyz; siteID++) {
        pID=get_particle_id(LatticeB, siteID);
        /*bleach_profile[] gives the probability that the protein is bleached */
        if( pID>=0 && ( 1.0 + rand() ) / ((double)RAND_MAX+1.0) < bleach_profile[siteID] ) { 
          bleached[pID]=1;  // Protein label is set from 0 to 1
        }
      }
      bleach_iter=iter;     // Remember iteration at which was bleached
    }
    /*--------------------------------------------------*/

    /*--------------------------------------------------*/
    /* EXPORT DATA */
    bleach_capture=( bleach_iter!=-1 && ((iter-bleach_iter)%bleach_iter_target)==0);
    if( iter%NTIME_PRINT == 0 || iter==NTIME+iter0 || bleach_capture)
    {
      /* logarithmic time interval */
      if(NITER_PRINT_LOG)
        NTIME_PRINT*=2;
      if(bleach_capture)
        bleach_iter_target*=2;


      output_counter++;
      /* Export time progress */
      printf("iter: %12lld (%5.2f%%); time: %12e; accepted %5.2f%%\n", iter, 100.0*(iter-iter0)/NTIME, kmc_time, (iter==0?0:100.0*(double)count_accepted/iter));
      export_time_progress(ftimeprogress, iter, kmc_time, count_accepted, 1);
      /* Export configuration */
      export_cnf(outdir, fcnf, output_counter, LatticeB, bleached, 'B');

      /* Export tracked particles */
      if(track_A)
        export_particle_tracking(ftrackA, kmc_time, xA, yA, zA, nA, 1);
      if(track_B)
        export_particle_tracking(ftrackB, kmc_time, xB, yB, zB, nB, 1);
    }
    /*--------------------------------------------------*/


    /*---------------*/
    /* UPDATE TIME: Random selection -> time updated regardless if process is accepted or rejected */
    rsm_k=1.0/(NenabledB*nu_Bmax); /* maximum time step size */
    kmc_dt = -rsm_k*log( ( 1.0 + rand() ) / ((double)RAND_MAX+1.0) );
    kmc_time+=kmc_dt;

    /*---------------*/
    /* SELECT EVENT  */
    SelectEnabledProcess(LatticeB, NenabledB, enabledB, process_regB, siteB0, xB,yB,zB, &pID, &move_direction, &siteID, &x_old,&y_old,&z_old, &sitef, &xf,&yf,&zf);

    /*---------------*/
    // ACCEPT/REJECT MOVE
    if(epsBB==0) { /* No A-A interactions*/
      if(EnergyLandscapeModuleB) { /* Disordered energy landscape */
        if(0==AcceptRejectS0NN_EL(LatticeB, siteID, sitef, EnergyLandscapeB, nu_Bmax, nu_B0, nu_B1))
          continue; /* REJECTED --> End of current time iteraction */
      } else {
        if(0==AcceptRejectS0NN(LatticeB, siteID, sitef, epsB1,           nu_Bmax, nu_B0, nu_B1))
          continue; /* REJECTED --> End of current time iteraction */ 
      }
    } else {
      if(EnergyLandscapeModuleA) {
         printf("Error: EnergyLandscapeModuleA and nearest neighbour interactions not yet implemented.\n");
      }else {
        if(0==AcceptRejectS6NN(LatticeB, siteID, sitef, epsBB, epsB1, nu_Bmax, nu_B0, nu_B1,interact_probB))
        continue; /* REJECTED --> End of current time iteraction */
      }
    }
    // else: ACCEPTED
    count_accepted++;

    //---------------------
    /* EXECUTE MOVE AND UPDATE LATTICE & PARTICLE POS */
    // update configuration
    update_cnf(LatticeB, siteB0, xB, yB, zB, pID, siteID, sitef, xf, yf, zf);

    /* update enabled processes and update the number of occupied neighbours */
    update_enabled_and_neighbours(LatticeB, &NenabledB, enabledB, process_regB, pID, siteID, x_old, y_old, z_old,sitef, xf,yf,zf, move_direction);
    
    // update counter (how many times did protein cross a boundary in case of periodic boundaries)
    if(x_cr!=0)      {xB_cross[pID]+=x_cr;} 
    else if(y_cr!=0) {yB_cross[pID]+=y_cr;} 
    else if(z_cr!=0) {zB_cross[pID]+=z_cr;}
    x_cr=0;y_cr=0;z_cr=0; // reset
    //---------------------

  } // last time step
  //###########################################################################################
  // BOTH A AND B PROTEINS
  } else  {
  for(iter=iter0; iter<=NTIME+iter0; iter++) 
  {
    /*--------------------------------------------------*/
    /* PHOTOBLEACH                                      */
    if(tbleach>=0 && bleach_iter==-1 && kmc_time>=tbleach ) { 
      // BLEACH A
      if(bleach_who=='A') {
        for(siteID=0; siteID<LatticeA->Nxyz; siteID++) {
          pID=get_particle_id(LatticeA, siteID);
          /*bleach_profile[] gives the probability that the protein is bleached */
          if( pID>=0 && ( 1.0 + rand() ) / ((double)RAND_MAX+1.0) < bleach_profile[siteID] ) { 
           bleached[pID]=1;  // Protein label is set from 0 to 1
          }
        }
      // BLEACH B
      } else if (bleach_who=='B') {
        for(siteID=0; siteID<LatticeB->Nxyz; siteID++) {
          pID=get_particle_id(LatticeB, siteID);
          /*bleach_profile[] gives the probability that the protein is bleached */
          if( pID>=0 && ( 1.0 + rand() ) / ((double)RAND_MAX+1.0) < bleach_profile[siteID] ) { 
           bleached[pID]=1;  // Protein label is set from 0 to 1
          }
        }
      }
      bleach_iter=iter;     // Remember iteration at which was bleached
    }
    /*--------------------------------------------------*/

    /*--------------------------------------------------*/
    /* EXPORT DATA */
    bleach_capture=( bleach_iter!=-1 && ((iter-bleach_iter)%bleach_iter_target)==0);
    if( iter%NTIME_PRINT == 0 || iter==NTIME+iter0 || bleach_capture)
    {
      output_counter++;

      /* logarithmic time interval */
      if(NITER_PRINT_LOG)
        NTIME_PRINT*=2;
      if(bleach_capture)
        bleach_iter_target*=2;

      /* Export time progress */
      printf("iter: %12lld (%5.2f%%); time: %12e; accepted %5.2f%%\n", iter, 100.0*(iter-iter0)/NTIME, kmc_time, (iter==0?0:100.0*(double)count_accepted/iter));
      export_time_progress(ftimeprogress, iter, kmc_time, count_accepted, 1);

      /* Export configuration */
      if(nA>0){
        export_cnf(outdir, fcnf, output_counter, LatticeA, bleached, 'A');
        /* Export tracked particles */
        if(track_A)
          export_particle_tracking(ftrackA, kmc_time, xA, yA, zA, nA, 1);

      }
      if(nB>0){
        export_cnf(outdir, fcnf, output_counter, LatticeB, bleached, 'B');
        if(track_B)
          export_particle_tracking(ftrackB, kmc_time, xB, yB, zB, nB, 1);
      }
    }
    /*--------------------------------------------------*/


    /*---------------*/
    /* UPDATE TIME: Random selection -> time updated regardless if process is accepted or rejected */
    rsm_kA=(NenabledA>0?NenabledA*nu_Amax:0);
    rsm_kB=(NenabledB>0?NenabledB*nu_Bmax:0);
    Delta_t=1.0/(rsm_kA+rsm_kB); /* maximum time step size */
    kmc_dt = -Delta_t*log( ( 1.0 + rand() ) / ((double)RAND_MAX+1.0) );
    kmc_time+=kmc_dt;

    // probability to attempt an A move: 
    if(NenabledA==0)
      attempt_who='A';
    else if(NenabledB==0)
      attempt_who='B';
    // probability to attempt A: rsm_kA*Delta_t
    else if ( ( rand() ) / ((double)RAND_MAX+1.0) <= rsm_kA*Delta_t ) {
      attempt_who='A';
    } else {
      attempt_who='B';
    }
    // ATTEMPT TO MOVE A
    if(attempt_who=='A'){

/*---------------*/
    /* SELECT EVENT  */
    SelectEnabledProcess(LatticeA, NenabledA, enabledA, process_regA, siteA0, xA,yA,zA, &pID, &move_direction, &siteID, &x_old,&y_old,&z_old, &sitef, &xf,&yf,&zf);

    /*---------------*/
    // ACCEPT/REJECT MOVE
    if(epsAA==0) { /* No A-A interactions*/
      if(EnergyLandscapeModuleA) { /* Disordered energy landscape */
        printf("Error: Energy landscape not implemented for two component simulations.\n"); goto TERMINATE_PROGRAM;
        if(0==AcceptRejectS0NN_EL(LatticeA, siteID, sitef, EnergyLandscapeA, nu_Amax, nu_A0, nu_A1))
          continue; /* REJECTED --> End of current time iteraction */
      } else {
        //if(0==AcceptRejectS0NN(LatticeA, siteID, sitef, epsA1,           nu_Amax, nu_A0, nu_A1))
        if(0==AcceptRejectS0NN_TwoComponent(LatticeA, LatticeB, epsAB, siteID, sitef, MAB))
          continue; /* REJECTED --> End of current time iteraction */ 
      }
    } else {
      if(EnergyLandscapeModuleA) {
         printf("Error: EnergyLandscapeModuleA and nearest neighbour interactions not yet implemented.\n"); goto TERMINATE_PROGRAM;
      }else {
        if(0==AcceptRejectS6NN_TwoComponent(LatticeA, LatticeB, epsAA, epsAB, siteID, sitef, MAB))
        continue; /* REJECTED --> End of current time iteraction */
      }
    }
    // else: ACCEPTED
    count_accepted++;

    //---------------------
    /* EXECUTE MOVE AND UPDATE LATTICE & PARTICLE POS */
    // update configuration
    update_cnf(LatticeA, siteA0, xA, yA, zA, pID, siteID, sitef, xf, yf, zf);

    /* update enabled processes and update the number of occupied neighbours */
    update_enabled_and_neighbours(LatticeA, &NenabledA, enabledA, process_regA, pID, siteID, x_old, y_old, z_old,sitef, xf,yf,zf, move_direction);
    
    // update counter (how many times did protein cross a boundary in case of periodic boundaries)
    if(x_cr!=0)      {xA_cross[pID]+=x_cr;} 
    else if(y_cr!=0) {yA_cross[pID]+=y_cr;} 
    else if(z_cr!=0) {zA_cross[pID]+=z_cr;}
    x_cr=0;y_cr=0;z_cr=0; // reset

    // ATTEMPT TO MOVE B
    } else if(attempt_who=='B'){
    /*---------------*/
    /* SELECT EVENT  */
    SelectEnabledProcess(LatticeB, NenabledB, enabledB, process_regB, siteB0, xB,yB,zB, &pID, &move_direction, &siteID, &x_old,&y_old,&z_old, &sitef, &xf,&yf,&zf);

    /*---------------*/
    // ACCEPT/REJECT MOVE
    if(epsBB==0) { /* No A-A interactions*/
      if(EnergyLandscapeModuleB) { /* Disordered energy landscape */
        printf("Error: Energy landscape not implemented for two component simulations.\n"); goto TERMINATE_PROGRAM;
        if(0==AcceptRejectS0NN_EL(LatticeB, siteID, sitef, EnergyLandscapeB, nu_Bmax, nu_B0, nu_B1))
          continue; /* REJECTED --> End of current time iteraction */
      } else {
        if(0==AcceptRejectS0NN_TwoComponent(LatticeB, LatticeA, epsAB, siteID, sitef, MBA))
          continue; /* REJECTED --> End of current time iteraction */ 
      }
    } else {
      if(EnergyLandscapeModuleA) {
         printf("Error: EnergyLandscapeModuleA and nearest neighbour interactions not yet implemented.\n"); goto TERMINATE_PROGRAM;
      }else {
        if(0==AcceptRejectS6NN_TwoComponent(LatticeB, LatticeA, epsBB, epsAB, siteID, sitef, MBA))
        continue; /* REJECTED --> End of current time iteraction */
      }
    }
    // else: ACCEPTED
    count_accepted++;

    //---------------------
    /* EXECUTE MOVE AND UPDATE LATTICE & PARTICLE POS */
    // update configuration
    update_cnf(LatticeB, siteB0, xB, yB, zB, pID, siteID, sitef, xf, yf, zf);

    /* update enabled processes and update the number of occupied neighbours */
    update_enabled_and_neighbours(LatticeB, &NenabledB, enabledB, process_regB, pID, siteID, x_old, y_old, z_old,sitef, xf,yf,zf, move_direction);
    
    // update counter (how many times did protein cross a boundary in case of periodic boundaries)
    if(x_cr!=0)      {xB_cross[pID]+=x_cr;} 
    else if(y_cr!=0) {yB_cross[pID]+=y_cr;} 
    else if(z_cr!=0) {zB_cross[pID]+=z_cr;}
    x_cr=0;y_cr=0;z_cr=0; // reset
    }
    //---------------------

  } // last time step
  }

  /* END CORE                                     */
  /*==============================================*/


  /*==============================================*/
  /* FREE MEMORY */
  /*===========================================*/
  TERMINATE_PROGRAM: {        /* FREE MEMORY */
  free_lattice_cube(LatticeA);
  free_lattice_cube(LatticeB);

  if(XA0!=NULL){free(XA0);} if(YA0!=NULL){free(YA0);} if(ZA0!=NULL){free(ZA0);}
  if(XB0!=NULL){free(XB0);} if(YB0!=NULL){free(YB0);} if(ZB0!=NULL){free(ZB0);}
  if(siteA0!=NULL){free(siteA0);} if(xA!=NULL){free(xA);} if(yA!=NULL){free(yA);} if(zA!=NULL){free(zA);}
  if(siteB0!=NULL){free(siteB0);} if(xB!=NULL){free(xB);} if(yB!=NULL){free(yB);} if(zB!=NULL){free(zB);}
  if(xA_cross!=NULL){free(xA_cross);} if(yA_cross!=NULL){free(yA_cross);} if(zA_cross!=NULL){free(zA_cross);}
  if(xB_cross!=NULL){free(xB_cross);} if(yB_cross!=NULL){free(yB_cross);} if(zB_cross!=NULL){free(zB_cross);}

  if(enabledA!=NULL){free(enabledA);} if(process_regA!=NULL){free(process_regA);}
  if(enabledB!=NULL){free(enabledB);} if(process_regB!=NULL){free(process_regB);}

  if(ibuff!=NULL){free(ibuff);};
  if(dbuff!=NULL){free(dbuff);};

  if(EnergyLandscapeA!=NULL) {free(EnergyLandscapeA);};
  if(EnergyLandscapeB!=NULL) {free(EnergyLandscapeB);};
  if(bleached       !=NULL) {free(bleached       );};
  if(bleach_profile !=NULL) {free(bleach_profile );};
  }; 
  /*==============================================*/

  printf("program terminated.\n");
  return(0);
}
