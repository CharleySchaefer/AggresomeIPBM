/*
  AcceptRejectS0NN - Single-component; binding energy to aggresome; NO nearest-neigbour interactions
  AcceptRejectS6NN - Single-component; binding energy to aggresome; With nearest-neigbour interactions

*/
/*
  ACCEPT/REJECT
  Single component with
  	> NO nearest-neighbour self-interactions
	> Binding interaction to aggresome

  Input:
    Lattice
    siteID:        site before hopping
    sitef:         site after hopping
    epsAA:         A-A interaction energy
    epsA1:         A binding energy to aggresome
*/
int AcceptRejectS0NN(LATTICE_CUBE *Lattice, int siteID, int sitef, double epsA1, double nu_Amax, double nu_A0, double nu_A1)
{
    double nuA, DeltaE;
    /* Protein is initially in cytosol */
    if (Lattice->site[siteID] <= -2) {   
      nuA=nu_A0;       
      /*   protein stays in cytosol  */
      if(       Lattice->site[sitef] == 0 ) {            
        
          if( nuA == nu_Amax ||( ((double)rand()) / ((double)RAND_MAX) <= nuA/nu_Amax ) ) {                                               // accept
            return(1);
          } else {
            return(0); // REJECT
          }
              
      /*   protein moves to aggresome - epsA1 has to be taken into account */
      }  else if(       Lattice->site[sitef] == 1 ) {
         DeltaE=-epsA1;  /* epsA1: binding energy */
         if(DeltaE<0) { /* decrease in energy */
           if(nuA>=nu_Amax || ( (double)rand() / ((double)RAND_MAX) <= nuA/nu_Amax ) ) { // accept
              return(1);
           } else {
             return(0); // REJECT
           }
         } else { /* increase in energy */
           if(nuA>=nu_Amax) { // accept 
             if( ((double)rand()) / ((double)RAND_MAX) <= exp(-DeltaE) ) {
               return(1);
             } else      {
               return(0); // REJECT
             }
           } else if( ((double)rand()) / ((double)RAND_MAX) <= nuA/nu_Amax*exp(-DeltaE) )   {
               return(1);
           } else  {
               return(0); // REJECT
           }
         }
      /*   protein tries to move neither to cytosol nor to aggresome  */
      }else {
        printf("Error: Trying to move particle to forbidden site\n");
        printf("Value of old site: %d, Value of attempted new site: %d", Lattice->site[siteID], Lattice->site[sitef]);
        return(-1);
      }
    /* Protein is initially in aggresome */  
    } else if (Lattice->site[siteID] >= 2) {            
      /*   protein stays in aggresome  */
      if(       Lattice->site[sitef] == 1 ){            
        
          if( nu_A1 == nu_Amax ) {                                               // accept
            return(1);
          } else if ( ((double)rand()) / ((double)RAND_MAX) <= nu_A1/nu_Amax ) { // accept with probability (weighted by rate)
            return(1);
          } else {
            return(0); // REJECT
          }        
      /*   protein moves to cytosol  */
      }  else if(       Lattice->site[sitef] == 0 ) {
         DeltaE=epsA1;  /* epsA1: binding energy */
         nuA=nu_A0;
         if(DeltaE<0) { /* decrease in energy */
           if ( (nuA>=nu_Amax) || ( (double)rand() /((double)RAND_MAX ) <= nuA/nu_Amax ) ) { // accept
              return(1);
           } else {
             return(0); // REJECT
           } 
         } else { /* increase in energy */
           if(nuA>=nu_Amax) { // accept 
             if(  (double)rand() /( (double)RAND_MAX ) <= exp(-DeltaE) ) {
              return(1);
             } else      {
               return(0); // REJECT
             } 
           } else if( ((double)rand()) / ((double)RAND_MAX) <= nuA/nu_Amax*exp(-DeltaE) )   {
             return(1);
           } else      {
             return(0); // REJECT
           }
         }
      /*   protein tries to move neither to cytosol nor to aggresome  */
      }else {
        printf("Error: Trying to move particle to forbidden site\n");
        printf("Value of old site: %d, Value of attempted new site: %d", Lattice->site[siteID], Lattice->site[sitef]);
        return(-1);
      }
    }
  return(1);
}



int AcceptRejectS0NN_EL(LATTICE_CUBE *Lattice, int siteID, int sitef, double *EnergyLandscape, double nu_Amax, double nu_A0, double nu_A1)
{
    double nuA, DeltaE, Ei;
    Ei=EnergyLandscape[siteID];
    /* Protein is initially in cytosol */
    if (Lattice->site[siteID] <= -2) {   
      nuA=nu_A0;    


   
      /*   protein stays in cytosol  */
      if((       Lattice->site[sitef] == 0 )||(       Lattice->site[sitef] == 1 )) {       
         DeltaE=EnergyLandscape[sitef ]-Ei;  
         if(DeltaE<=0) { /* decrease in energy */
           if(nuA>=nu_Amax || ( (double)rand() / ((double)RAND_MAX) <= nuA/nu_Amax ) ) { // accept
              return(1);
           } else {
             return(0); // REJECT
           }
         } else { /* increase in energy */
           if(nuA>=nu_Amax) { // accept 
             if( ((double)rand()) / ((double)RAND_MAX) <= exp(-DeltaE) ) {
               return(1);
             } else      {
               return(0); // REJECT
             }
           } else if( ((double)rand()) / ((double)RAND_MAX) <= nuA/nu_Amax*exp(-DeltaE) )   {
               return(1);
           } else  {
               return(0); // REJECT
           }
         }
      /*   protein tries to move neither to cytosol nor to aggresome  */
      }else {
        printf("Error: Trying to move particle to forbidden site\n");
        printf("Value of old site: %d, Value of attempted new site: %d", Lattice->site[siteID], Lattice->site[sitef]);
        return(-1);
      }
    /* Protein is initially in aggresome */  
    } else if (Lattice->site[siteID] >= 2) {            
      /*   protein stays in aggresome  */
      if(       Lattice->site[sitef] == 1 ){    
          DeltaE=EnergyLandscape[sitef ]-Ei;            
          if(DeltaE<=0) {
          if( nu_A1 == nu_Amax ) {                                               // accept
            return(1);
          } else if ( ((double)rand()) / ((double)RAND_MAX) <= nu_A1/nu_Amax ) { // accept with probability (weighted by rate)
            return(1);
          } else {
            return(0); // REJECT
          }} else{
          if( nu_A1 == nu_Amax ) {                                               // accept
            return(1);
          } else if ( ((double)rand()) / ((double)RAND_MAX) <= nu_A1/nu_Amax*exp(-DeltaE) ) { // accept with probability (weighted by rate)
            return(1);
          } else {
            return(0); // REJECT
          }    
}
      /*   protein moves to cytosol  */
      }  else if(       Lattice->site[sitef] == 0 ) { 
          DeltaE=EnergyLandscape[sitef ]-Ei;   
         nuA=nu_A0;
         if(DeltaE<=0) { /* decrease in energy */
           if ( (nuA>=nu_Amax) || ( (double)rand() /((double)RAND_MAX ) <= nuA/nu_Amax ) ) { // accept
              return(1);
           } else {
             return(0); // REJECT
           } 
         } else { /* increase in energy */
           if(nuA>=nu_Amax) { // accept 
             if(  (double)rand() /( (double)RAND_MAX ) <= exp(-DeltaE) ) {
              return(1);
             } else      {
               return(0); // REJECT
             } 
           } else if( ((double)rand()) / ((double)RAND_MAX) <= nuA/nu_Amax*exp(-DeltaE) )   {
             return(1);
           } else      {
             return(0); // REJECT
           }
         }
      /*   protein tries to move neither to cytosol nor to aggresome  */
      }else {
        printf("Error: Trying to move particle to forbidden site\n");
        printf("Value of old site: %d, Value of attempted new site: %d", Lattice->site[siteID], Lattice->site[sitef]);
        return(-1);
      }
    }
  return(1);
}


/*
  ACCEPT/REJECT
  Single component with
  	> Nearest-neighbour self-interactions
	> Binding interaction to aggresome

  Input:
    Lattice
    siteID:        site before hopping
    sitef:         site after hopping
    epsAA:         A-A interaction energy
    epsA1:         A binding energy to aggresome
    interact_prob: bare probability of moving 
*/
int AcceptRejectS6NN(LATTICE_CUBE *Lattice, int siteID, int sitef, double epsAA, double epsA1, double nu_Amax, double nu_A0, double nu_A1, double *interact_prob)
{
    double nuA, DeltaE;
    int neigh_dif = Lattice->num_neighbours[siteID] - (Lattice->num_neighbours[sitef] - 1); /* difference in number of neighbours after moving */

    /* Protein is initially in cytosol */
    if (Lattice->site[siteID] <= -2) {   
      nuA=nu_A0;       
      /*   protein stays in cytosol  */
      if(       Lattice->site[sitef] == 0 ) {            
        
        /* decrease in energy --*/
        if ( (neigh_dif <= 0) ){ 
          if( nuA == nu_Amax ||( ((double)rand()) / ((double)RAND_MAX) <= nuA/nu_Amax ) ) {                                               // accept
            return(1);
          } else {
            return(0); // REJECT
          }
        /* increase in energy --*/
        } else {                 
          if( nuA == nu_Amax ) {                                               // accept with probability
            if  ( (double)rand()  / ((double)RAND_MAX) <= interact_prob[neigh_dif - 1]  ) {
              return(1);
            }  else {
              return(0); // REJECT
            }                                                                // accept with probability (weighted by rate)
          } else if ( (double)rand() / ((double)RAND_MAX) <= interact_prob[neigh_dif - 1]*nuA/nu_Amax ) {
              return(1);
            }  else {
              return(0); // REJECT
            } 
          
        }         
      /*   protein moves to aggresome - epsA1 has to be taken into account */
      }  else if(       Lattice->site[sitef] == 1 ) {
         DeltaE=neigh_dif*epsAA-epsA1;  /* epsA1: binding energy */
         if(DeltaE<0) { /* decrease in energy */
           if(nuA>=nu_Amax || ( (double)rand() / ((double)RAND_MAX) <= nuA/nu_Amax ) ) { // accept
              return(1);
           } else {
             return(0); // REJECT
           }
         } else { /* increase in energy */
           if(nuA>=nu_Amax) { // accept 
             if( ((double)rand()) / ((double)RAND_MAX) <= exp(-DeltaE) ) {
               return(1);
             } else      {
               return(0); // REJECT
             }
           } else if( ((double)rand()) / ((double)RAND_MAX) <= nuA/nu_Amax*exp(-DeltaE) )   {
               return(1);
           } else  {
               return(0); // REJECT
           }
         }
      /*   protein tries to move neither to cytosol nor to aggresome  */
      }else {
        printf("Error: Trying to move particle to forbidden site\n");
        printf("Value of old site: %d, Value of attempted new site: %d", Lattice->site[siteID], Lattice->site[sitef]);
        return(-1);
      }
    /* Protein is initially in aggresome */  
    } else if (Lattice->site[siteID] >= 2) {            
      /*   protein stays in aggresome  */
      if(       Lattice->site[sitef] == 1 ){            
        
        /* decrease in energy --*/
        if ( (neigh_dif <= 0) ){ 
          if( nu_A1 == nu_Amax ) {                                               // accept
            return(1);
          } else if ( ((double)rand()) / ((double)RAND_MAX) <= nu_A1/nu_Amax ) { // accept with probability (weighted by rate)
            return(1);
          } else {
            return(0); // REJECT
          }
        /* increase in energy --*/
         } else {                 
        if( nu_A1 == nu_Amax ) {                                               // accept with probability
             if  ( (double)rand()  / ( (double)RAND_MAX ) <= interact_prob[neigh_dif - 1]  ) {
              return(1);
            }  else {
              return(0); // REJECT
            }                                                                 // accept with probability (weighted by rate)
          } else if ( ((double)rand()) / ((double)RAND_MAX) <= interact_prob[neigh_dif - 1]*nu_A1/nu_Amax ) {
              return(1);
          } else {
              return(0); // REJECT
          } 
        }         
      /*   protein moves to cytosol  */
      }  else if(       Lattice->site[sitef] == 0 ) {
         DeltaE=neigh_dif*epsAA+epsA1;  /* epsA1: binding energy */
         nuA=nu_A0;
         if(DeltaE<0) { /* decrease in energy */
           if ( (nuA>=nu_Amax) || ( (double)rand() /((double)RAND_MAX ) <= nuA/nu_Amax ) ) { // accept
              return(1);
           } else {
             return(0); // REJECT
           } 
         } else { /* increase in energy */
           if(nuA>=nu_Amax) { // accept 
             if(  (double)rand() /( (double)RAND_MAX ) <= exp(-DeltaE) ) {
              return(1);
             } else      {
               return(0); // REJECT
             } 
           } else if( ((double)rand()) / ((double)RAND_MAX) <= nuA/nu_Amax*exp(-DeltaE) )   {
             return(1);
           } else      {
             return(0); // REJECT
           }
         }
      /*   protein tries to move neither to cytosol nor to aggresome  */
      }else {
        printf("Error: Trying to move particle to forbidden site\n");
        printf("Value of old site: %d, Value of attempted new site: %d", Lattice->site[siteID], Lattice->site[sitef]);
        return(-1);
      }
    }
  return(1);
}


//###########################################################################################
// TWO COMPONENT
// TODO: upgrade f regions of attraction
// M12 should be <1 ; represents the relative mobility, e.g., M12=nuAB/nuA0
int AcceptRejectS0NN_TwoComponent(LATTICE_CUBE *Lattice1, LATTICE_CUBE *Lattice2, double eps12, int siteID, int sitef, double M12)
{
    double  DeltaE,   M;
    int vali=(Lattice2->site[siteID]!=0), 
        valf=(Lattice2->site[sitef]!=0);

    // Mobility prefactor
    if (vali) {
      if(valf)      /* protein '1' is in phase '2' and remains in phase '2' */
        M=M12;
      else
        M=M12*0.5;  /* protein '1' is in phase '2' and moves to phase '1'   */
    } else {
      if(valf) 
        M=M12*0.5;  /* protein '1' is in phase '1' and moves to phase '2'   */
      else
        M=1.0;      /* protein '1' is in phase '1' and remains in phase '1' */
    }

    /*Change in energy after hop*/
    DeltaE= -( (valf-1)  -  vali )*eps12; // eps12>0: attractive interaction 
    if(M==1) {
      if(DeltaE<=0 || ( ((double)rand()) / ((double)RAND_MAX+1) <= exp(-DeltaE) ) ){
        return(1); // accept
      } else {
        return(0); // reject
      }
    } else {
      if( (DeltaE<=0 && ((double)rand()) / ((double)RAND_MAX+1) <=M) || ( ((double)rand()) / ((double)RAND_MAX+1) <= M*exp(-DeltaE) ) ){
        return(1); // accept
      } else {
        return(0); // reject
      }
    }

  return(1);
}

int AcceptRejectS6NN_TwoComponent(LATTICE_CUBE *Lattice1, LATTICE_CUBE *Lattice2, double eps11, double eps12, int siteID, int sitef, double M12)
{
    double DeltaE,   M;
    int vali=(Lattice2->site[siteID]!=0), 
        valf=(Lattice2->site[sitef]!=0);

    // Mobility prefactor
    if (vali) {
      if(valf) 
        M=M12;      /* protein '1' is in phase '2' and remains in phase '2' */
      else
        M=M12*0.5;  /* protein '1' is in phase '2' and moves to phase '1'   */
    } else {
      if(valf) 
        M=M12*0.5;  /* protein '1' is in phase '1' and moves to phase '2'   */
      else
        M=1.0;      /* protein '1' is in phase '1' and remains in phase '1' */
    }

    /*Change in energy after hop*/
    DeltaE=  -(       (Lattice2->site[sitef]!=0) -  (Lattice2->site[siteID]!=0) )*eps12         // eps12>0: attractive interaction
             +( Lattice1->num_neighbours[siteID] - (Lattice1->num_neighbours[sitef] - 1))*eps11; 
    if(M==1) {
      if(DeltaE<=0 || ( ((double)rand()) / ((double)RAND_MAX+1) <= exp(-DeltaE) ) ){
        return(1); // accept
      } else {
        return(0); // reject
      }
    } else {
      if( (DeltaE<=0 && ((double)rand()) / ((double)RAND_MAX+1) <=M) || ( ((double)rand()) / ((double)RAND_MAX+1) <= M*exp(-DeltaE) ) ){
        return(1); // accept
      } else {
        return(0); // reject
      }
    }

  return(1);
}


