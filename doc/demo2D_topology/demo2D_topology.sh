#!/bin/bash

#=====================================================
# USER SETTINGS
# To use the demo:
#     copy demo*.sh to the build directory, 
#     change location to the build directory, and run 
#     ./demo*.sh
do_compile=0
do_run=1
do_postprocess=1
do_debug=0

# Simulation input
seed1=1
seed2=3
outdir0="demo2D_topology"
topology_input_file="../doc/demo2D_topology/arrowTopology.txt" #../doc/demo2D_topology/TopologyTwoAggresomesNX50NY150.dat"
bleach_input_file="../doc/demo2D_topology/BleachProfile.dat"
NX=50          #
NY=150         # spatial resolution in X and Y direction should be the same
NZ=1
NA=1200        # number of A proteins
EPSAA=2.0      # A-A nearest-neighbour interaction energy
EPSA1=2.0      # Binding energy of A to aggresome
NITER=100000000 # Number of iterations
NPRINT=100000  # Number of iterations after which output is generated
tbleach=2e4    # Bleaching time in monte carlo units

# Postprocessing input
dx=20           # resolution [nm]
DA=0.2 # Diffusivity in [micron^2/s]
# END USER SETTINGS
#=====================================================

# CONVERT MONTE CARLO UNITS TO REAL UNITS
# Lattice size in micrometer (dx is given in nanometer)
LX=`awk -vdx=$dx -vNX=$NX 'BEGIN {print dx*NX/1000}'`
LY=`awk -vdx=$dx -vNY=$NY 'BEGIN {print dx*NY/1000}'`

# time [s] per unit of kmc time (DA=[micron^2/s]; dx=[nm])
timefac=`awk -vdx=$dx -vDA=$DA 'BEGIN {print dx*dx/DA/1000000}'`

#COMPILE
if [ $do_compile -eq 1 ]; then
  echo "Compile"
  if [ $do_debug -eq 1 ] ; then
  ./compile.sh -g
  else
  ./compile.sh
  fi
fi



# Get absolute path to utils
pushd .. > /dev/null
  utils=`pwd`"/utils"
popd > /dev/null

# LOOP OVER RANDOM NUMBER SEEDS
for rseed in `seq $seed1 $seed2` ; do
  outdir=$outdir0/seed$rseed

  #RUN
  if [ $do_run -eq 1 ]; then

    echo "Create cnf and img folders"
    mkdir -p $outdir0
    mkdir -p $outdir
    mkdir -p $outdir/cnf
    mkdir -p $outdir/img
 

    echo "Run simulation"
    if [ $do_debug -eq 1 ] ; then
      valgrind --leak-check=full ./IBM --outdir $outdir --Nx $NX --Ny $NY --Nz $NZ \
        --NA $NA --Niter $NITER --Niter-print $NPRINT \
        --rseed $rseed \
        --epsAA $EPSAA \
        --epsA1 $EPSA1 \
        --import-topology  $topology_input_file \
        --bleach-time $tbleach \
        --bleach-profile $bleach_input_file
    else
      ./IBM --outdir $outdir --Nx $NX --Ny $NY --Nz $NZ \
        --NA $NA --Niter $NITER --Niter-print $NPRINT \
        --rseed $rseed \
        --epsAA $EPSAA \
        --epsA1 $EPSA1 \
        --import-topology  $topology_input_file \
        --bleach-time $tbleach \
        --bleach-profile $bleach_input_file
    fi
  fi

#POSTPROCESS
# 1. Create .png configuration images
# 2. Analyse structure factor and radial correlation function
# 3. Extract characteristic length scales from correlation function
if [ $do_postprocess -eq 1 ]; then
  echo "Postprocess"

  echo "  Create images in img folder"
  pushd $outdir 
    cp $utils/cnf2png.m .
    echo "cnf2png" | octave --no-gui
  popd

  echo "  Calculate Structure Factor in SF folder"
  mkdir -p $outdir/SF ; rm $outdir/SF/SF*

  i=1 ; check=1
  > $outdir/AnalysedLengthScales.out
  while [ $check -eq 1 ] ; do
    numstr=`printf "%05d" $i`
    fin=$outdir/cnf/cnfA$numstr.out
    if [ -f $fin ] ; then

      # read time in units [s] from timeprogress.out
      kmctime=`awk -vtimefac=$timefac -vData=$i <"$outdir/timeprogress.out" 'BEGIN{Nheader=1} NR==(Data+Nheader) {print timefac*$2}'`




      echo "Calculating structure factor of $fin" 
      fout=$outdir/SF/SF$numstr.out # output file
      $utils/StructureFactorAndCorrelationFunction.o --file $fin --Lx $LX --Ly $LY --stretch-matrix > $fout  

      #echo "Analyse length scales"
      Nheader=1 #number of header lines in output file
      Rcol=3   # column with x values
      CRcol=4  # column with y values
      echo $kmctime " " `$utils/get_roots_from_datafile.sh $fout $Rcol $CRcol $Nheader` >> $outdir/AnalysedLengthScales.out
    else
      check=0;
    fi
    i=$((i * 2))
  done

  # Plot SF (uses gnuplot)
  pushd $outdir
    $utils/plotSF.sh
    $utils/plotCR.sh
  popd
  plot="set terminal pngcairo enhanced ; set output \"$outdir/LengthScales.png\"
  set log xy; set xlabel \"time [s]\" ; set ylabel \"R [{/Symbol m}m]\"
  set grid 
  set key right bottom
  plot \"$outdir/AnalysedLengthScales.out\" u 1:2 w p ps 1.2 lw 2 pt 6 title \"first zero\", \\
  \"\" u 1:4  w p ps 1.2 lw 2 pt 10 title \"first minimum\", \\
  \"\" u 1:3  w p ps 1.2 lw 2 pt 8  title \"second zero\", \\
  0.4*($DA*x)**(1.0/3) w l lw 2 lc rgb \"black\" title \"t^{1/3}\"
  "
  printf "$plot" | gnuplot
  echo "Calling octave to analyse FRAP transient"
  echo "addpath('$utils'); analyseFRAP('$outdir', 'A')" | octave --no-gui
  echo "Calling gnuplot to plot frap transients."
  $utils/plotfrap.sh $outdir $timefac
  $utils/plotfrap_loglog.sh $outdir $timefac
fi

done # FOR LOOP RANDOM NUMBER SEEDS


# output written to AverageLengthScale.out

if [ $do_postprocess -eq 1 ] ; then

# AVERAGE COARSENING STATISTICS
echo "Averaging Coarsening Data"
echo "addpath('$utils') ; average_lengthscale('$outdir0')" | octave --no-gui
plot="set terminal pngcairo enhanced
set output '$outdir0/AverageLengthScale.png'
set log xy; set xlabel \"time [s]\" ; set ylabel 'R [{/Symbol m}m]'
  set grid 
unset key
  plot '$outdir0/AverageLengthScale.out' u 1:(\$2-0.07*$dx/20):3:4 w xye ps 1.2 lw 2 pt 6 lc rgb 'black' notitle, \\
  0.3*($DA*x)**(1.0/3) w l lw 2 lc rgb \"black\" title 't^{1/3}'
"
echo "$plot" | gnuplot



# AVERAGE FRAP
echo "Averaging FRAP Data"
echo "addpath('$utils') ; average_frap('$outdir0')" | octave --no-gui
plot="set terminal pngcairo enhanced
set output '$outdir0/frap.png'
set log xy; set xlabel \"time [s]\" ; set ylabel 'frap [-]'
  set grid 
unset key
  plot '$outdir0/averagefrap.out' u 1:2:3:4 w xye ps 1.2 lw 2 pt 6 lc rgb 'black' notitle
"
echo "$plot" | gnuplot
fi

