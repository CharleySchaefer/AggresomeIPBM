#include "../lib/ZiltoidLIB/ZiltoidLIB.h" /* used for coor2ind_square */
/*----------------------------------------------------------------------------*/
int readTimeProgress(char *, long long unsigned *, double *,  long long unsigned *, int*, double *);

/*----------------------------------------------------------------------------*/
/*
  READ TIMEPROGRESS FILE

  input: 
    ftimeprogress - name of the timeprogress file
    dbuff:          [char]memory buffer, should be long enough to read a line from timeprogress
  output:
    iter:    number of iterations read from the last line
    kmctime: simulation time (Monte Carlo units) read from the last line
*/
int readTimeProgress(char *ftimeprogress, long long unsigned *iter, double *kmctime, long long unsigned *count_accepted, int *output_counter, double *dbuff) // TODO: column recognition
{
  int    Nlines;
  double *row=dbuff;

  if( access( ftimeprogress, F_OK ) == -1 ) // file does not exist
      {printf("Error: No time progress file %s found!\n", ftimeprogress); return(0);}


  /* Read last line */
  countLines(ftimeprogress, &Nlines);      /* Get number of lines */
  dreadRow(ftimeprogress, Nlines-1, row);  /* Read last line */

  (*iter)   =(long long unsigned)row[0];  /* first  column: Number of iterations     */
  (*kmctime)=row[1];  /* second column: Time (Monte Carlo units) */
  *count_accepted=(long long unsigned)((row[2])*(*iter)*0.01); /*0.01: convert percentage to fraction */ 
  *output_counter=Nlines-1; // -1 : subtract header line
  return(1);
}

/*----------------------------------------------------------------------------*/
/*
  IMPORT CONFIGURATION FILE

  read_iCNF_file - read cnf file with integers
  read_dCNF_file - read cnf file with doubles

  Purpose: should be able to read either 
           i.  particle positions, or
           ii. structures (fixed aggresomes/stress granules/nucleoid)

  input:
    fname :  File name
    Lattice: Only expected lattice size (Lattice->Nx, Lattice->Ny, Lattice->Nz) is used
    ibuff:   [int] buffer of length Ny

  output:
    cnf=

*/
int read_iCNF_file(LATTICE_CUBE *Lattice, char *fname, int *ibuff, int *cnf)
{
  int i,j;
  int verbose=0;
  int nx, siteID, NX,NY,Nlines, Ndata, Nheader;
  int *arr=ibuff;
  if (Lattice->Nz>1) {
    printf("Error: Reading in 3D cnf file not yet implemented.\n"); return(0); 
    /* for 3D implementation, coor2ind_square below should be replaced by coor2ind_cube (see lib/ZiltoidLIB/LatticeLIB) */
  }

  /* Get cnf file properties */
  if(!analyse_data_file_properties(fname, &Nlines, &NY, &Nheader, &Ndata, verbose)) {
    printf("Error: analyse_data_file_properties() failed.\n"); return(0);
  }
  if(Nheader!=0) {
    printf("Warning: unexpected header lines in %s\n", fname);
    if(Nlines!=Ndata){ 
      printf("Warning: Number of data lines unequal to total number of lines in %s.\n", fname);
    }
  }
  NX=Ndata;
  if(NX!=Lattice->Nx || NY!=Lattice->Ny) {
    printf("Error: unexpected matrix size."); return(0); 
  }
	
  /* Read data */
  for (i=0; i<NX; i++) {
    ireadRow(fname, i, arr, &nx); // length line with NY integers (nx is unused)
    for (j=0; j<NY; j++) {
      siteID=coor2ind_square(Lattice, i, j);
      cnf[siteID]=arr[j];
    }
  }

  return(1);
}

