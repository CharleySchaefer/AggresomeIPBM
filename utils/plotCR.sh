#!/bin/bash

  # Plot Correlation Function (uses gnuplot)
  plot="set terminal pngcairo enhanced; set output \"CR.png\"
  # AXES
  unset log x; set xlabel \"r\";    # set xtics format \"10^{%%T}\"
  unset log y; set ylabel \"C(r)\"; # set ytics format \"10^{%%T}\"
  unset key
  plot 1.0/0"
  for file in `find SF/SF*` ; do
    plot="$plot,\\
\"$file\" u 3:4 w l lt 1 lw 2 lc rgb \"black\""
  done
  plot="$plot
  "
  printf "$plot" | gnuplot
