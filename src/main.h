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
#define VERSION "0.200705"  // based on date YYMMDD
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 
#include <assert.h>
#include <dirent.h>
#include <errno.h>

#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#include <limits.h>
#include <stdint.h>

/*------------------------------------------------*/
/* Random number generator */
// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/
//Copyright (c) 2006,2007 Mutsuo Saito, Makoto Matsumoto and Hiroshima University.
//Copyright (c) 2012 Mutsuo Saito, Makoto Matsumoto, Hiroshima University and The University of Tokyo.
// See SFMT/LICENSE.txt
#include "../lib/SFMT/SFMT.h" 
//static uint32_t r_sfmt;
sfmt_t sfmt;
// initialised in main.c; used in SelectEnabledProcess.h
/*------------------------------------------------*/

/*------------------------------*/
/*
  CONTENTS
  > Program Settings

  > Random Number Generator

  > Static Properties
   -lattice dimensions
   -interaction energies

  > Dynamic Properties
   -hop directions
   -enabled Monte Carlo processes

  > Exporting Simulation Results
   - timeprogress file
   - configuration files
   - particle tracking
*/
/*------------------------------*/
/* PROGRAM SETTINGS             */
void printHelp(char *pname, int rseed, unsigned long long NTIME, unsigned long long NTIME_PRINT, int NX, int NY, int NZ) {
  printf("\n  USAGE: %s <arguments>\n", pname);
    printf("\n  Arguments:\n");
    printf("  INFO\n");
    printf("  --help\n");
    printf("  --version\n");
    printf("  INITIALISATION\n");
    printf("  --continue-previous ; [under construction: limited functionality] continue previous simulation \n");
    printf("  --rseed <integer>  (random number seed; default value is 1)\n");
    printf("  TIME STEPS / NUMBER OF ITERATIONS\n");
    printf("  --Niter <integer> ; default: %llu\n", NTIME);
    printf("  SIMULATION BOX:\n");
    printf("  --Nx <integer>           ; default: %d\n", NX);
    printf("  --Ny <integer>           ; default: %d\n", NY);
    printf("  --Nz <integer>           ; default: %d\n", NZ);
    printf("  --import-topology <input file> ; create excluded volume and regions of attraction\n");
    printf("  PROPERTIES OF PROTEINS:\n");
    printf("  (examples for protein A ; for protein A replace e.g. --NA by --NB)\n");
    printf("  --NA <integer>  ; Number of proteins A; default 1\n");
    printf("  --nuA0 <float>  ; hop rate in cytosol;      default nuA0=1\n");
    printf("  --nuA1 <float>  ; hop rate in aggresome;    default nuA1=1\n");
    printf("  --nuAB <float>  ; hop rate of A in B phase; default nuAB=nuA0\n");
    printf("  --nuBA <float>  ; hop rate of B in A phase; default nuBA=nuB0\n");
    printf("  --epsAA <float> ; A-A nearest-neighbour interaction energy [kT]\n");
    printf("  --epsAB <float> ; A-B shared-site interaction energy [kT]\n");
    printf("  --epsA1 <float> ; interaction energy [kT] of protein with\n");
    printf("                    aggresome; overridden by disordered DOS.\n");
    printf("      Define disordered density of states (DOS; overrides --epsA1)\n");
    printf("      DOS consists of a Gaussian exp( -(E-Emean)^2/(2*sigma^2  ) )\n");
    printf("      and a fraction of exponential traps: exp( -(E-Emax)/width  )\n");
    printf("      the DOS can be specified for cytosol (--E-xxx-0) and regions of\n");
    printf("      attraction (--E-xxx-1) independently.\n");
    printf("  --E-gauss-mean-A-0  <float>  ; mean of Gaussian [kT]\n");
    printf("  --E-gauss-sigma-A-0 <float>  ; std of Gaussian [kT]\n");
    printf("  --E-trap-frac-A-0   <float>  ; fraction of exponential traps\n");
    printf("  --E-trap-max-A-0    <float>  ; maximum of exponential traps [kT]\n");
    printf("  --E-trap-width-A-0  <float>  ; width of exponential traps [kT]\n");
    printf("  PHOTOBLEACHING:\n");
    printf("  --bleach-time    ; [kmc units];\n");
    printf("  --bleach-profile  <input file>; create position-dependent bleaching probability\n");
    printf("  EXPORT:\n");
    printf("  --Niter-print <integer> (interval between \'snapshots\') ;\n");
    printf("                           default: %llu\n", NTIME_PRINT);
    printf("  --Niter-print-log  - logarithmic time interval for \'snapshots\'\n");
    
  printf("\n  --track-A ; export file with protein A positions\n");
  printf("\n  --track-B ; export file with protein B positions\n");
} 

/*------------------------------*/
/* RANDOM NUMBER GENERATORS     */
//#include "power_random.h" /* functions for random numbers */
/* Cryptographically secure, unbiased, uniform random numbers and more */
/*#include "lib/libsodium-win64/include/sodium.h"
#define SODIUM_STATIC 1 */

/*------------------------------*/
/* STATIC PROPERTIES OF LATTICE */
#define MAX_STR_L 10000 /* overrides MAX_STR_L in ZiltoidLIB; 
                          also used to allocate memory for strings 
                          in main.c 
                       */
#include "../lib/ZiltoidLIB/ZiltoidLIB.h" /* Cuboid lattice --> Lattice->site */
/*
  Cuboid geometry:    NX x NY x NZ
  Boundary condition: no-flux/non-permeable wall, no interaction with the wall

  Site values:
  -1: excluded volume (e.g., location of nucleoid)
   0: vacant, cytosol
   1: vacant, region of attraction (e.g., aggresome / stress granule)
 <-2: abs(value)-2 = particle ID of protein in cytosol         
  >2: value-2      = particle ID of protein in region of attraction    
      
*/
int get_particle_id(LATTICE_CUBE *Lattice, int site_id){
  int pID = Lattice->site[site_id];
  pID =(pID>0?pID:-pID)-2; // site values -1, 0 and 1 indicate the 'phase' and don't refer to any particle ID
  return pID;
}


/*------------------------------*/
/* HOPPING OF PROTEINS          */
#define HOP_DIRECTIONS 6
/*
  0: x->x+1
  1: x->x-1
  2: y->y+1
  3: y->y-1
  4: z->z+1
  5: z->z-1
*/
#include "move_track_decide.h"



/*------------------------------*/
/* EXPORTING RESULTS            */
#include "Import.h"
#include "Initialise.h"
#include "GenerateEnergyLandscape.h"
#include "AcceptReject.h"
#include "SelectEnabledProcess.h"
#include "Export.h"






