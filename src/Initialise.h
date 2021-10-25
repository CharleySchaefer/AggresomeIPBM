#include "../lib/ZiltoidLIB/ZiltoidLIB.h" /* Lattice structure in  ZiltoidLIB/LatticeLIB */
/*
  INITIALISE CELL TOPOLOGY

  Lattice of size Nx*Ny*Nz, with site values:
  -1: excluded volume (e.g., location of nucleoid)
   0: vacant, cytosol
   1: vacant, region of attraction (e.g., aggresome / stress granule)

  input:
    fname: cnf file
    Lattice->Nx, Lattice->Ny, Lattice->Nz, Lattice->Nxyz (Lattice->Nxyz=Lattice->Nx*Lattice->Ny*Lattice->Nz)
    Lattice->site - memory should be allocated 
  output
    Lattice->site - values
*/
int initialise_cell_topology(char *fname, LATTICE_CUBE *Lattice, int *ibuff) 
{
  int i,value;

  /* Read */
  if(!read_iCNF_file(Lattice, fname,ibuff, Lattice->site)) {
    printf("Error: read_iCNF_file() failed.\n"); return(0);
  }

  /* Check values */
  for(i=0; i<Lattice->Nxyz; i++) {
    value=Lattice->site[i];
    if( value!=-1 && value!=0 && value!=1 ) {
      printf("Error: unexpected value %d read from %s", value, fname); return(0);
    }
  }

  return(1);
}
/**/
/*
  ibuff: memory buffer of length Ny
  TODO: non-sharp boundaries / bleach proteins with probability -> bleach_profile should be of type double*
                                                                -> read_dCNF_file variant in Import.h should be written. 

*/
int initialise_bleach_profile(char *fname, LATTICE_CUBE *Lattice, int *ibuff, int *bleach_profile)
{
  int i,value;
  

  /* Read */
  if(!read_iCNF_file(Lattice, fname,ibuff, bleach_profile)) {
    printf("Error: read_iCNF_file() failed.\n"); return(0);
  }
  /* Check values */
  for(i=0; i<Lattice->Nxyz; i++) {
    value=bleach_profile[i];
    if( value!=0 && value!=1 ) {
      printf("Error: unexpected value %d read from %s", value, fname); return(0);
    }
  }
  return(1);
}

/* determine number of output from timeprogress file; read last cnf file
  NOTE: protein ID will be relabelled (reading ParticleTracking not yet implemented).
*/
int import_previous_cnf(LATTICE_CUBE *Lattice, char *ftimeprogress, char *fcnf, char *outdir, int *site0, int *X0, int *Y0, int *Z0, int *ibuff, char label)
{
  int i, Nlines, Ncol, Nheader, Ndata, verbose=0;
  int *cnf=(int*)malloc(Lattice->Nxyz*sizeof(int));
  /* Read time progress */ // TODO: check exist using assert
  if(!analyse_data_file_properties(ftimeprogress, &Nlines, &Ncol, &Nheader, &Ndata, verbose)) {
    printf("Error: analyse_data_file_properties() failed.\n"); return(0);
  }

  /* Read last cnf file */ // TODO: check exist using assert
  sprintf(fcnf, "%s/cnf/cnf%c%05d.out", outdir, label, Ndata); // label: 'A' or 'B'
  if(!read_iCNF_file(Lattice, fcnf, ibuff, cnf)) {
    printf("Error: read_iCNF_file() failed.\n"); return(0);
  }
  int pID=0;
  for (i=0; i<Lattice->Nxyz; i++){
    if(cnf[i]!=0){
      site0[pID]=i;
      if(     Lattice->site[i] == 0) {Lattice->site[i] = -2 - pID;} // cytoplasm
      else if(Lattice->site[i] == 1) {Lattice->site[i] =  2 + pID;} // aggresome
      else {printf("Error with assigning the particle ID."); return(0);}
      /* define initial positions. X0+p is pointer arithmetic! */
      ind2coor_cube(Lattice, i, X0+pID, Y0+pID, Z0+pID);
      pID++;
    }
  }
  free(cnf);
  return(1);
}
