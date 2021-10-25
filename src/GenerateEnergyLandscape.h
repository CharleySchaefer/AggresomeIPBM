#ifndef GEN_ENERGY_LANDSCAPE_H
  #define GEN_ENERGY_LANDSCAPE_H 1
  #include "../lib/ZiltoidLIB/ZiltoidLIB.h" /* lattice structure*/

/*
  Gaussian Disorder Model:

    P(E) = exp(  -(E-Emean)^2/(2*sigma^2)  )/(2*pi*sigma^2)

  Cummulative P(E):

    p(E) = integral_{-inf}^{E}  P(e) de = ...

  Inverse:

    E(p) = Emean - sqrt(2)*sigma*erfinv(1-2*p)

*/
/*
int GenerateEnergyLandscape(LATTICE_CUBE *, 
 double, double , double, double,
 double *);
int export_EnergyLandscape(char *, char *,LATTICE_CUBE *, double *);
*/
  int GenerateEnergyLandscape(LATTICE_CUBE *,
 double , double , 
 double , double , double,
 double , double ,
 double , double , double,
 double *);

  int export_EnergyLandscape(char *, char *,LATTICE_CUBE *, double *, char);

#endif
