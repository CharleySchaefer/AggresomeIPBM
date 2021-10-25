#!/bin/bash

# DEMO DESCRIPTION:
#
#  THE BIOPHYSICS: 
#    Protein B phase separates to form B-enriched droplets due to an attractive B-B interaction
#    Protein A diffuses into the B-enriched droplets due to an attractive A-B interaction
# THE EXPERIMENT:
#    At some point of time, A proteins inside a laser focus region are photobleached
#    After this time, we measure the fluorescent recovery in that region due to the diffusive influx
#    of unbleached A proteins.
#    As a bonus, also the coarsening dynamics of the B-enriched droplets are measured.
# THIS SCRIPT:
#    First, the cell is modelled (box size) and the number of proteins, their interactions,
#    and dynamical properties are set.
#    Also, the time of the photobleaching and the bleaching profile are specified.
#    Then, the random number seeds are specified (multiple simulations are averaged to collect
#    statistics).
#    For every of these simulations, postprocessing is performed (measuring coarsening dynamics
#    and the FRAP transient).

#-----------------------------------
# USER SETTINGS
do_debug=0
do_run=1
do_postprocess=1

outdir0='demo2D_2CFRAP'
NITER=10000000
NPRINT=20000
NX=25
NY=75
seed1=1
seed2=5

# Protein B: driving force for phase separation
NB=550    # Number of B proteins
epsBB=2.1  # B-B interaction energy [kT]

# Protein A: diffuses and binds to B
NA=125     # Number of A proteins
epsAB=2.2  # A-B interaction energy [kT]
nuAB=0.25  # mobility of A in B phase 

dx=40.0           # resolution [nm]
DA=0.001 # Diffusivity in [micron^2/s]

# PHOTOBLEACHING
xbleach="0.5" #micron
ybleach="0.5" #micron
Rbleach="0.4" #micron
fbleach="$outdir0/BleachProfile.out"
tbleach=5000

# Postprocessing input
#=====================================================

# CONVERT MONTE CARLO UNITS TO REAL UNITS
# Lattice size in micrometer (dx is given in nanometer)
LX=`awk -vdx=$dx -vNX=$NX 'BEGIN {print dx*NX/1000}'`
LY=`awk -vdx=$dx -vNY=$NY 'BEGIN {print dx*NY/1000}'`

# time [s] per unit of kmc time (DA=[micron^2/s]; dx=[nm])
timefac=`awk -vdx=$dx -vDA=$DA 'BEGIN {print dx*dx/DA/1000000}'`

# Get absolute path to utils
pushd .. > /dev/null
  utils=`pwd`"/utils"
popd > /dev/null

mkdir -p $outdir0
echo "addpath('$utils') ; createBleachProfile2D($NX,$NY,$LX,$LY,$xbleach, $ybleach, $Rbleach, '$fbleach');" | octave --no-gui



# LOOP OVER RANDOM NUMBER SEEDS
for rseed in `seq $seed1 $seed2` ; do

# MAKE OUTPUT DIRECTORIES
outdir=$outdir0/seed$rseed
mkdir -p $outdir
mkdir -p $outdir/img

#SIMULATION
if [ $do_run -eq 1 ] ; then
if [ $do_debug -eq 1 ] ; then
valgrind  --leak-check=full --track-origins=yes ./IBM --outdir $outdir --Nx $NX --Ny $NY --Nz 1 \
      --NA $NA --NB $NB --Niter $NITER --Niter-print $NPRINT
else
./IBM --outdir $outdir --Nx $NX --Ny $NY --Nz 1 \
      --rseed $rseed \
      --NA $NA --NB $NB --epsAB $epsAB --epsBB $epsBB --nuAB $nuAB \
      --Niter $NITER --Niter-print $NPRINT \
      --bleach-profile $fbleach --bleach-time $tbleach \
#      --continue-previous
fi
fi

#POSTPROCESSING
if [ $do_postprocess -eq 1 ] ; then
  # CREATE PNG CONFIGURATION IMAGES
  pushd $outdir > /dev/null
    cp $utils/cnf2png_frap.m .
    echo "cnf2png_frap" | octave --no-gui
  popd > /dev/null

  # ANALYSE FRAP TRANSIENTS
  echo "Calling octave to analyse FRAP transient"
  echo "addpath('$utils'); analyseFRAP('$outdir', 'A')" | octave --no-gui
  echo "Calling gnuplot to plot frap transients."
  $utils/plotfrap.sh $outdir $timefac
  $utils/plotfrap_loglog.sh $outdir $timefac

  # ANALYSE STRUCTURE FACTOR
  echo "  Calculate Structure Factor in SF folder"
  mkdir -p $outdir/SF ; rm $outdir/SF/SF*

  i=1 ; check=1
  > $outdir/AnalysedLengthScales.out
  while [ $check -eq 1 ] ; do
    numstr=`printf "%05d" $i`
    fin=$outdir/cnf/cnfB$numstr.out
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
  pushd $outdir > /dev/null
    $utils/plotSF.sh
    $utils/plotCR.sh
  popd > /dev/null 
  plot="set terminal pngcairo enhanced ; set output '$outdir/LengthScales.png'
  set log xy; set xlabel \"time [s]\" ; set ylabel 'R [{/Symbol m}m]'
  set grid 
  set key right bottom
  plot '$outdir/AnalysedLengthScales.out' u 1:2 w p ps 1.2 lw 2 pt 6 title 'first zero', \\
  '' u 1:4  w p ps 1.2 lw 2 pt 10 title 'first minimum', \\
  '' u 1:3  w p ps 1.2 lw 2 pt 8  title 'second zero', \\
  0.6*($DA*x)**(1.0/3) w l lw 2 lc rgb 'black' title 't^{1/3}'
  "
  printf "$plot" | gnuplot

fi # end postprocessing
done # end loop over random number seeds

# AVERAGE COARSENING STATISTICS
echo "addpath('$utils') ; average_lengthscale('$outdir0')" | octave --no-gui
# output written to AverageLengthScale.out

if [ $do_postprocess -eq 1 ] ; then
plot="set terminal pngcairo enhanced
set output '$outdir0/AverageLengthScale.png'
set log xy; set xlabel \"time [s]\" ; set ylabel 'R [{/Symbol m}m]'
  set grid 
unset key
  plot '$outdir0/AverageLengthScale.out' u 1:2:3:4 w xye ps 1.2 lw 2 pt 6 lc rgb 'black' notitle, \\
  0.6*($DA*x)**(1.0/3) w l lw 2 lc rgb \"black\" title 't^{1/3}'
"
echo "$plot" | gnuplot
fi
