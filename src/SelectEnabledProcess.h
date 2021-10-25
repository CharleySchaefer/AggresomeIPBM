
/*

  input:
    Nenabled, enabled, process_reg;   - Register of enabled events
    x,y,z                             - Arrays with all protein positions
  output:
    pID                               - ID of selected particle
    move_direction                    - Direction the selected particle moves to
    siteID, x_old,y_old,z_old         - Current position of the moving particle
    sitef,  xf,yf,f                   - New     position of the moving particle
    
*/
int SelectEnabledProcess(LATTICE_CUBE *Lattice, int Nenabled, int *enabled, int *process_reg, int *site0, int *x,int *y,int *z, int *pID, int *move_direction, int *siteID, int *x_old,int *y_old,int *z_old, int *sitef, int *xf,int *yf,int *zf)
{
  int x_cr,y_cr,z_cr; /*currently not used*/
    /* Select random process */
    //int sel_proc_enabled_index=(int)(Nenabled-1)*((double)rand()/(double)RAND_MAX); 
     // Tried rand(), random() and lrand48() before ; this gives a bias
     // if there are many proteins/processes in the system
     // if some weird bias is found again, may try sfmt_genrand_uint64?
     //init_gen_rand or init_by_array must be called before this function
     int sel_proc_enabled_index= (int)(sfmt_genrand_uint32(&sfmt)%Nenabled);
 
    int sel_proc_id = enabled[sel_proc_enabled_index]; /* Identify the corresponding event */
    if(sel_proc_id == -1) {
      printf("Error with process id selection\n"); exit(-1);
    } else if(process_reg[sel_proc_id] != sel_proc_enabled_index) { /* check consistency */
      printf("Process register and 'enabled'-array don't fit\n"); return(0);  
    }

    /* Interpret Process  */
    (*pID) = sel_proc_id / HOP_DIRECTIONS;             /* ID of the moving particle */
    (*move_direction) = sel_proc_id % HOP_DIRECTIONS;  /* Direction of movement     */

    /* Initial protein position */
    (*siteID) = site0[(*pID)];
    (*x_old) = x[(*pID)]; (*y_old) = y[(*pID)]; (*z_old) = z[(*pID)];

    /* New protein position */
    if (-1 == next_position(Lattice->Nx, Lattice->Ny, Lattice->Nz, (*move_direction), (*x_old), (*y_old), (*z_old), xf, yf, zf,
                  &x_cr, &y_cr, &z_cr)){
      printf("Exiting. Error in function next_position"); return(0);  
    }
    *sitef = coor2ind_cube(Lattice, (*xf), (*yf), (*zf));
  return(1);
}
