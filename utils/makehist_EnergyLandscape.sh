#!/bin/bash
outdir=$1

#-------------------------------------------
# GET HISTOGRAM DATA
script="addpath('../utils/DensityOfStates'); 
Emat=importdata('$outdir/EnergyLandscape.out');

Nbins=500;
[Erow, prow]=get_E_hist(Emat, Nbins);
fp=fopen('$outdir/HistEnergyLandscape.out', 'w');
for i=1:length(Erow)
  fprintf(fp, '%12e %12e\n', Erow(i), prow(i));
end
fclose(fp);"

echo "$script" | octave --no-gui

#-------------------------------------------
# PLOT
script="set terminal pngcairo enhanced
set output '$outdir/HistEnergyLandscape.png'

set xlabel 'E/kT'
set ylabel 'P(E)'
set log y

plot '$outdir/HistEnergyLandscape.out' u 1:2 w l lw 2 lc rgb 'black'
"
echo "$script" | gnuplot



