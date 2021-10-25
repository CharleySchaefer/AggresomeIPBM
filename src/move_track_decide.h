/* written by Johannes Wesseler 2020 */
#ifndef MOVE_TRACKER_H
    #define MOVE_TRACKER_H 1
    #include <stdio.h>
    #include <stdlib.h>
    
    #include "../lib/ZiltoidLIB/ZiltoidLIB.h"

    // get new (x,y,z) coordinate, given the initial one and the move_direction
    int next_position(int, int, int, int, int, int, int, int *,        // NX, NY, NZ, move_direction, x y z currently
                        int *, int *, int *, int *, int *);             // x y z pointers for next position, pointers to indicators
                                                                        // if boundary crossed in x y z
    int get_backward_direction(int); 

    // Update enabled events
    int add_enabled(int *, int *, int *, int ); // pointer to enabled moves array, pointer to number of currently possible moves,
                                                // Pointer to process register array, process ID

    int rm_enabled(int *, int *, int *, int );  // pointer to enabled moves array, pointer to number of currently possible moves,
                                                // Pointer to process register array, process ID

    // Check if move is allowed: target site should be vacant & should not be excluded volume 
    int possible_y1n0(LATTICE_CUBE *, int, int, int, int);              // Lattice, current x, current y, current z, move ID

    // Count number of occupied nearest neighbours
    int calc_occupied_neighbours(LATTICE_CUBE *);

    int update_cnf(LATTICE_CUBE *Lattice, int *, int *, int *, int *,  int , int ,int , int ,int ,int );
    int update_enabled_and_neighbours(LATTICE_CUBE *, int *, int *, int *, int , int , int , int , int ,int , int ,int ,int , int );


    /*-------------------------------------------------------------------------------------------------------------*/
    /* GET BACKWARD MOVE */
    /*-------------------------------------------------------------------------------------------------------------*/
    	
/*
  get_backward_direction(): 
  Get inverse of move_direction
  move_direction   backwards_direction
  0: x->x+1        1: x->x-1
  1: x->x-1        0: x->x+1
  2: y->y+1        3: y->y-1
  3: y->y-1        2: y->y+1
  4: z->z+1        5: z->z-1
  5: z->z-1        4: z->z+1 
*/
int get_backward_direction(int move_direction) {
  // A direct equation for the backwards direction would be
  // int backwards_direction = move_direction + ((move_direction+1)%2) - (move_direction % 2); 
  // but a switch is more efficient I think. Especially modulo operations cost much time, as they are
  // 3 real operations (C = A % B is equivalent to C = A – B * (A / B))
  // int backwards_direction;
  switch(move_direction)
      {
        case 0:                        // x->x+1
          return 1; //backwards_direction = 1;      
          break;
        case 1:
          return 0; //backwards_direction = 0;
          break;
        case 2:
          return 3; //backwards_direction = 3;
          break;
        case 3:
          return 2; //backwards_direction = 2;
          break;
        case 4:
          return 5; //backwards_direction = 5;
          break;
        case 5:
          return 4; //backwards_direction = 4;
          break;
        default:
          return -1;
          break;
      }
  //return backwards_direction;
}

        
    /*-------------------------------------------------------------------------------------------------------------*/
    /* CALCULATES OF OCCUPIED NEAREST NEIGHBOUR CELLS */
    /*-------------------------------------------------------------------------------------------------------------*/
    int calc_occupied_neighbours(LATTICE_CUBE *Lattice)
    {
        int x, y, z;
        int occ_count;

        int site_temp;
        int site_value_temp;
        int side_var;

        for(int cell = 0; cell < Lattice->Nxyz; cell++)
        {
            occ_count = 0;
            ind2coor_cube(Lattice, cell, &x, &y, &z);
            

            /* In x direction */
            for(int diff_x = -1; diff_x <=1; diff_x+=2)                                     /* Site below or above the considered site  */
            {
                side_var = x+diff_x;
                if(side_var == Lattice->Nx || side_var == -1) {continue;}                   /* The old spot was next to boundary        */
                site_temp = coor2ind_cube(Lattice, side_var, y, z);
                site_value_temp = Lattice->site[site_temp];

                if(site_value_temp <= -2 || site_value_temp >=2) {occ_count++;}             /* Check if neighbours occupied             */
            }

            /* In y direction */
            for(int diff_y = -1; diff_y <=1; diff_y+=2)                                     /* Site below or above the considered site  */
            {
                side_var = y+diff_y;
                if(side_var == Lattice->Ny || side_var == -1) {continue;}                   /* The old spot was next to boundary*/
                site_temp = coor2ind_cube(Lattice, x, side_var, z);
                site_value_temp = Lattice->site[site_temp];

                if(site_value_temp <= -2 || site_value_temp >=2) {occ_count++;}             /* Check if neighbours occupied               */
            }

            /* In z direction */
            for(int diff_z = -1; diff_z <=1; diff_z+=2)                                     /* Site below or above the considered site  */
            {
                side_var = z+diff_z;
                if(side_var == Lattice->Nz || side_var == -1) {continue;}                   /* The old spot was next to boundary*/
                site_temp = coor2ind_cube(Lattice, x, y, side_var);
                site_value_temp = Lattice->site[site_temp];

                if(site_value_temp <= -2 || site_value_temp >=2) {occ_count++;}             /* Check if neighbours occupied               */
            }

            // Add occ_count to Lattice variable
            Lattice->num_neighbours[cell] = occ_count;

        }
      return(1);
    }







    /*
    Selecting processes is rendered efficient by keeping track of events that are
    possible (a particle cannot move to a lattice site that is already occupied).
    While selecting a process, it is randomly selected from an 'enabled' array:

    eventID=enabled[ rand()%Nenabled ]

    A process is added to or removed from the enabled list using the functions
    add_enabled and rm_enabled below. 
    The arguments to these functions are: 

    processID:   the ID of a process (a process is, e.g., particle i moves right)
    Nenabled:    number of enabled processes (updated every time step)
    enabled:     array with the ID's of enabled processes
                its allocated length should equal the maximum number of processes
                should be initialised with value -1 on each entry
    process_reg: registers the position of an enabled process in the 'enabled' array
                should be initialised with value -1 on each entry
    */


    /*-------------------------------------------------------------------------------------------------------------*/
    /* ADD */
    /*-------------------------------------------------------------------------------------------------------------*/

    int add_enabled(int *enabled, int *Nenabled, int *process_reg, int processID) 
    {

        if(process_reg[processID]!=-1) {return(0);} /* processID is already enabled */
        else
        {
            enabled[(*Nenabled)]=processID;
            /* register place of processID in enabled list */
            process_reg[processID]=(*Nenabled); 
            (*Nenabled)++;
            return(1);
        }
    }



    /*-------------------------------------------------------------------------------------------------------------*/
    /* REMOVE */
    /*-------------------------------------------------------------------------------------------------------------*/
    int rm_enabled(int *enabled, int *Nenabled, int *process_reg, int processID) 
    {
        int reg_id = process_reg[processID];
        if(reg_id==-1) {return(0);}                     /* processID is already disabled */
            
        else if(reg_id != (*Nenabled)-1) 
        { /* move last entry */
            enabled[reg_id]= enabled[(*Nenabled)-1];    /* replace the removed entry */
            process_reg[enabled[reg_id]] = reg_id;
        }
        process_reg[processID]=-1; /* delete the removed entry */
        enabled[(*Nenabled)-1]=-1; /* remove last entry*/
        (*Nenabled)--;
        return(1);
    }








    /*-------------------------------------------------------------------------------------------------------------*/
    /* DECIDE IF PARTICLE MOVE ALLOWED, CONSIDERING NEIGHBOURS & LATTICE STRUCTURE */
    /*-------------------------------------------------------------------------------------------------------------*/

    int possible_y1n0(LATTICE_CUBE *Lattice, int xc, int yc, int zc, int move)
    {
        int tmp_h;                                                                              /* Site index of particle after move */
        int xfl = xc, yfl = yc, zfl = zc;                                                          /* Coor       of particle after move */

        switch(move) 
        {
            case 0:
                xfl++;
                if( xfl == Lattice->Nx) {return(0);}                                            /* Move outside forbidden if not periodic */  
                else {
                    tmp_h = Lattice->site[coor2ind_cube(Lattice, xfl, yfl, zfl)];
                    if(tmp_h == 0 || tmp_h == 1 )      {return(1);}                            /* Move allowed if free & not nucleus */
                    else                                {return(0);}
                }
                break;
            
            case 1:
                xfl--;
                if( xfl == -1) {return(0);}
                else {
                    tmp_h = Lattice->site[coor2ind_cube(Lattice, xfl, yfl, zfl)];
                    if(tmp_h == 0 || tmp_h == 1 )      {return(1);}                            /* Move allowed if free & not nucleus */
                    else                                {return(0);}
                }
                break;

            case 2:
                yfl++;
                if( yfl == Lattice->Ny) {return(0);}
                else {
                    tmp_h = Lattice->site[coor2ind_cube(Lattice, xfl, yfl, zfl)];
                    if(tmp_h == 0 || tmp_h == 1 )      {return(1);}                            /* Move allowed if free & not nucleus */
                    else                                {return(0);}
                }
                break;
            
            case 3:
                yfl--;
                if( yfl == -1) {return(0);}
                else {
                    tmp_h = Lattice->site[coor2ind_cube(Lattice, xfl, yfl, zfl)];
                    if(tmp_h == 0 || tmp_h == 1 )      {return(1);}                            /* Move allowed if free & not nucleus */
                    else                                {return(0);}
                }
                break;

            case 4:
                zfl++;
                if( zfl == Lattice->Nz) {return(0);}
                else {
                    tmp_h = Lattice->site[coor2ind_cube(Lattice, xfl, yfl, zfl)];
                    if(tmp_h == 0 || tmp_h == 1 )      {return(1);}                            /* Move allowed if free & not nucleus */
                    else                                {return(0);}
                }
                break;
            
            case 5:
                zfl--;
                if( zfl == -1) {return(0);}
                else {
                    tmp_h = Lattice->site[coor2ind_cube(Lattice, xfl, yfl, zfl)];
                    if(tmp_h == 0 || tmp_h == 1 )      {return(1);}                            /* Move allowed if free & not nucleus */
                    else                                {return(0);}
                }
                break;    
        }

        /* If none of the cases is true or no condition was met, then something went wrong */
        printf("Error in the function that decides which move is possible\n");
        return(-1);


    }






    /*-------------------------------------------------------------------------------------------------------------*/
    /* CALCULATE THE FOLLOWING SITE FOR PARTICLE MOVE WITH THE ASSUMPTION THAT THE MOVE IS POSSIBLE */
    /*-------------------------------------------------------------------------------------------------------------*/
    int next_position(int NX, int NY, int NZ, int move_direction, int xcur, int ycur, int zcur, int *xfol, int *yfol,
                        int *zfol, int *x_cross, int *y_cross, int *z_cross)
    {
        switch(move_direction)
        {
        case 0:
            *xfol = xcur+1; *yfol=ycur; *zfol=zcur;
            return(0);
            break;
            
        case 1:
            *xfol = xcur-1; *yfol=ycur; *zfol=zcur;
            return(0);
            break;

        case 2:
            *yfol = ycur+1; *xfol=xcur; *zfol=zcur;
            return(0);
            break;
            
        case 3:
            *yfol = ycur-1; *xfol=xcur; *zfol=zcur;
            return(0);
            break;

        case 4:
            *zfol = zcur+1; *xfol=xcur; *yfol=ycur;
            return(0);
            break;
            
        case 5:
            *zfol = zcur-1; *xfol=xcur; *yfol=ycur;
            return(0);
            break;
        }

        printf("Error in function that calculates next position\n");
        return(-1);
    }



/*
  update_cnf()
  // input
  pID:              particle identifier
  siteID:           initial position of protein
  sitef; xf,yf,zf:  target position of protein

  // input/output
  site0:   list with all particle positions (grid index)
  x,y,z:   lists with all particle positions (x,y,z coordinate)

  // output
  0: failed    execution of update_cnf()
  1: succesful execution of update_cnf()
  
*/
int update_cnf(LATTICE_CUBE *Lattice, int *site0, int *x, int *y, int *z,  int pID, int siteID,int sitef, int xf,int yf,int zf)
{
      /* UPDATE OCCUPATION OF OLD LATTICE SITE*/
      if(Lattice->site[siteID] <= -2) {Lattice->site[siteID] = 0;}                        /* Previously occupied cytoplasm          */
      else if(Lattice->site[siteID] >=2 ) {Lattice->site[siteID] = 1;}                    /* Previously occupied aggresome          */
      else  {
        printf("Error during updating the old Lattice site. Site was empty and "              /* Previously unoccupied !!!              */
                "program tries to move a particle from this empty cell to somewhere " 
                "else\n");
                return(0);
      }

      /* UPDATE OCCUPATION OF NEW LATTICE SITE*/
      if(Lattice->site[sitef] == 0) {Lattice->site[sitef] = -2 - pID ;}                   /* Previously occupied cytoplasm          */
      else if(Lattice->site[sitef] == 1 ) {Lattice->site[sitef] = 2 + pID ;}              /* Previously occupied aggresome          */
      else
      {
        printf("Error during updating the new Lattice site. Site isn't empty "                /* Forbidden cell                         */
                "or is nucleus and program tries to move a particle in this " 
                "forbidden cell\n");
                return(0);
      }

      /* UPDATE PARTICLE POSITION */
      site0[pID] = sitef;
      x[pID] = xf; y[pID] = yf; z[pID] = zf;


  return(1);
}

/*
  pID: particle identifier of the protein that moves
  
*/
int update_enabled_and_neighbours(LATTICE_CUBE *Lattice, int *Nenabled_pointer, int *enabled, int *process_reg, int pID, int site_old, int x_old, int y_old, int z_old,int sitef, int xf,int yf,int zf, int move_direction)
{
      int site_temp;
      int side_var; 
  int pID_nn; /* particle ID of a nearest neighbour */
  int NX=Lattice->Nx, NY=Lattice->Ny,NZ=Lattice->Nz;
  int Nenabled=(*Nenabled_pointer);


    // Update number of occupied nearest neighbours
    Lattice->num_neighbours[site_old]++;                                                    /* Particle is now neighbour of old cell  */
    Lattice->num_neighbours[sitef]--;                                                     /* Particle isn't neighbour of new cell   */
    


      /*--------------------------------------------*/
      /* UPDATE ALLOWED MOVES */
      int possible_var;
      int process_ID;


      /* Update moved particle itself */
      int backwards_direction=get_backward_direction(move_direction) ; /* inverse of move_direction */


      add_enabled(enabled, &Nenabled, process_reg, pID*HOP_DIRECTIONS + backwards_direction);              /* Add backwards move                 */

      for(int direction = 0; direction < HOP_DIRECTIONS; direction++)                                          
      {
        if (direction == backwards_direction) {continue;}                                         /* Skip backwards move here           */

        process_ID = pID*HOP_DIRECTIONS + direction;                                                       /* Test if move possible              */
        possible_var = possible_y1n0(Lattice, xf, yf, zf, direction);
        if(possible_var == -1)
        {
          printf("Error with possible_y1n0 function\n");
          exit(EXIT_FAILURE);
        } 
        

        if(possible_var == 1) {add_enabled(enabled, &Nenabled, process_reg, process_ID);}         /* Add if possible                    */
        else {rm_enabled(enabled, &Nenabled, process_reg, process_ID);}                           /* Remove if not                      */
      }


      /* Update all particles around the old position */
      /* In x direction */
      for(int diff_x = -1; diff_x <=1; diff_x+=2)                                                 /* Site below or above the old particle site  */
      {
        side_var = x_old+diff_x;
        if(side_var == NX || side_var == -1) {continue;}                                          /* The old spot was next to boundary*/
        site_temp = coor2ind_cube(Lattice, side_var, y_old, z_old);
        pID_nn=get_particle_id(Lattice, site_temp);
        if (pID_nn>=0)
        {
          if(pID_nn == pID) {continue;}                                                 /* Backwards already enabled & Neigh_num updated */
          add_enabled(enabled, &Nenabled, process_reg, HOP_DIRECTIONS*pID_nn + (diff_x+1)/2 );
        }
        Lattice->num_neighbours[site_temp]--;                                                     /* One neighbour less as old site now empty */
      }

      /* In y direction */
      for(int diff_y = -1; diff_y <=1; diff_y+=2)                                                 /* Site below or above the old particle site  */
      {
        side_var = y_old+diff_y;
        if(side_var == NY || side_var == -1) {continue;}
        site_temp = coor2ind_cube(Lattice, x_old, side_var, z_old);
        pID_nn=get_particle_id(Lattice, site_temp);
        if (pID_nn>=0) {
          if(pID_nn == pID) {continue;}                                               /* Backwards already enabled & Neigh_num updated */
          add_enabled(enabled, &Nenabled, process_reg, HOP_DIRECTIONS*pID_nn + (diff_y+1)/2 + 2);
        }
        Lattice->num_neighbours[site_temp]--;                                                     /* One neighbour less as old site now empty */
      }

      /* In z direction */
      for(int diff_z = -1; diff_z <=1; diff_z+=2)                                                 /* Site below or above the old particle site  */
      {
        side_var = z_old+diff_z;
        if(side_var == NZ || side_var == -1) {continue;}
        site_temp = coor2ind_cube(Lattice, x_old, y_old, side_var);
        pID_nn=get_particle_id(Lattice, site_temp);
        if (pID_nn>=0) {
          if(pID_nn == pID) {continue;}  /* Backwards already enabled & Neigh_num updated */
          add_enabled(enabled, &Nenabled, process_reg, HOP_DIRECTIONS*pID_nn + (diff_z+1)/2 + 4);
        }
        Lattice->num_neighbours[site_temp]--;                                                     /* One neighbour less as old site now empty */
      }


      /* Update all particles around the new position */

      /* In x direction */
      for(int diff_x = -1; diff_x <=1; diff_x+=2)                                                 /* Site below or above the new particle site  */
      {
        side_var = xf+diff_x;
        if(side_var == NX || side_var == -1) {continue;}                                          /* The old spot was next to boundary          */
        if(side_var == x_old) {continue;}                                                         /* No particle here                           */
        
        site_temp = coor2ind_cube(Lattice, side_var, yf, zf);
        pID_nn=get_particle_id(Lattice, site_temp);

        if(pID_nn>=0) {                                                                           /* Check if neighbours occupied         */                          
          rm_enabled(enabled, &Nenabled, process_reg, HOP_DIRECTIONS*pID_nn + (diff_x+1)/2 );     /* Stop neighbours from hopping on particle   */
        }

        Lattice->num_neighbours[site_temp]++;                                                     /* One neighbour more as new site now occupied */
      }

      /* In y direction */
      for(int diff_y = -1; diff_y <=1; diff_y+=2)                                                 /* Site below or above the new particle site  */
      {
        side_var = yf+diff_y;
        if(side_var == NY || side_var == -1) {continue;}                                          /* The old spot was next to boundary          */
        if(side_var == y_old) {continue;}                                                         /* No particle here                           */

        site_temp = coor2ind_cube(Lattice, xf, side_var, zf);
        pID_nn=get_particle_id(Lattice, site_temp);

        if(pID_nn>=0)   {                                                                         /* Check if neighbours occupied               */
          rm_enabled(enabled, &Nenabled, process_reg, HOP_DIRECTIONS*pID_nn + (diff_y+1)/2 + 2);
        }

        Lattice->num_neighbours[site_temp]++;                                                     /* One neighbour more as new site now occupied */
      }

      /* In z direction */
      for(int diff_z = -1; diff_z <=1; diff_z+=2)                                                 /* Site below or above the new particle site  */
      {
        side_var = zf+diff_z;
        if(side_var == NZ || side_var == -1) {continue;}                                          /* The old spot was next to boundary          */
        if(side_var == z_old) {continue;}                                                         /* No particle here                           */

        site_temp = coor2ind_cube(Lattice, xf, yf, side_var);
        pID_nn=get_particle_id(Lattice, site_temp);

        if(pID_nn>=0) {                                                                           /* Check if neighbours occupied               */
          rm_enabled(enabled, &Nenabled, process_reg, HOP_DIRECTIONS*pID_nn + (diff_z+1)/2 + 4);
        }
        Lattice->num_neighbours[site_temp]++;                                                     /* One neighbour more as new site now occupied */
      }

  (*Nenabled_pointer)=Nenabled;
  return(1);
}



#endif
