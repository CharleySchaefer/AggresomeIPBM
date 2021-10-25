#!/bin/bash
outdir=$1
timefac=$2  # convert kmctime to real time
plot="set terminal pngcairo enhanced
set output \'$outdir/frap_loglog.png\'

set log xy
set xlabel \"time [s]\"
set xtics format \"10^{%%T}\"
set ylabel \"FRAP [%%]\"

set grid lt 1 lw 0.5 lc rgb \"grey\"

unset key

plot \'$outdir/frap.out\' u ($timefac*\$1):2 w p ps 2 pt 6 lc rgb \"black\""


printf "$plot" | gnuplot
