#!/bin/bash

#=====================================================
# USER SETTINGS
# To use the demo:
#     copy demo*.sh to the build directory, 
#     change location to the build directory, and run 
#     ./demo*.sh
do_compile=1
do_run=1
do_postprocess=1
do_debug=0

# Simulation input
outdir="demo3D_FreeDiffusion"
NX=30          #
NY=30         # spatial resolution in X and Y direction should be the same
NZ=30
NA=50        # number of A proteins
EPSAA=0.0      # A-A nearest-neighbour interaction energy
NITER=1000000 # Number of iterations

# Postprocessing input
dx=20           # resolution [nm]
DA=0.3 # Diffusivity in [micron^2/s]
# END USER SETTINGS
#=====================================================

#COMPILE
if [ $do_compile -eq 1 ]; then
  echo "Compile"
  if [ $do_debug -eq 1 ] ; then
    ./compile.sh -g
  else
  ./compile.sh
  fi
fi

#RUN
if [ $do_run -eq 1 ]; then
  mkdir -p $outdir
  mkdir -p $outdir/cnf ; rm $outdir/cnf/cnf*.out
  echo "Run simulation"
  if [ $do_debug -eq 1 ] ; then
valgrind -v  --track-origins=yes --leak-check=full --show-leak-kinds=all ./IBM --Nx $NX --Ny $NY --Nz $NZ \
        --NA $NA --Niter $NITER --Niter-print-log \
        --epsAA $EPSAA \
        --track-A \
        --outdir $outdir
  else
     ./IBM --Nx $NX --Ny $NY --Nz $NZ \
           --NA $NA --Niter $NITER --Niter-print-log \
           --epsAA $EPSAA \
           --track-A \
           --outdir $outdir
  fi
fi

#POSTPROCESS
# 1. Create .png configuration images
# 2. Analyse structure factor and radial correlation function
# 3. Extract characteristic length scales from correlation function
if [ $do_postprocess -eq 1 ]; then
  echo "Postprocess"

  echo "Get RMSD from ParticleTracking file: generates rmsd.out"
  pushd $outdir
    printf "addpath('../../utils'); getRMSDFromParticleTracking('.', $dx, $DA, 'A')" | octave --no-gui
  popd

  # Plot RMSD (uses gnuplot)
  plot="set terminal pngcairo enhanced ; set output '$outdir/RMSD.png'
  set log xy; set xlabel 'time [s]' ; set xtics format '10^{%%T}';
  set ylabel 'RMSD [{/Symbol m}m]' ; set ytics format '10^{%%T}';
  set grid 
  set key left top reverse Left
  plot '$outdir/rmsd.out' u 1:2 w p ps 1.4 lw 3 pt 6 lc rgb '#FF4444' title 'simulation', \\
  sqrt(6*$DA*x) w l lw 2.5 lc rgb 'black' title '(6Dt)^{1/2} (3D diffusion)'
  "
  printf "$plot" | gnuplot
fi

