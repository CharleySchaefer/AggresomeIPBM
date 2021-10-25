#!/bin/bash
  plot="set terminal pngcairo enhanced; set output \"SF.png\"
  # AXES
  set log x; set xlabel \"q\";    set xtics format \"10^{%%T}\"
  set log y; set ylabel \"S(q)\"; set ytics format \"10^{%%T}\"
  set grid
  unset key
  plot 1.0/0"
  for file in `find SF/SF*` ; do
    plot="$plot,\\
\"$file\" u 1:2 w l lt 1 lw 2 lc rgb \"black\""
  done
  plot="$plot
  "
  printf "$plot" | gnuplot
