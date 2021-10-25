set terminal pngcairo enhanced
set output 'whole_frap_loglog.png'

set log xy
set xlabel "time [s]"
set xrange [1:1000]
set xtics format "10^{%T}"
set ylabel "frap [%]"

set grid lt 1 lw 0.5 lc rgb "grey"

set key left top
TIME0=600
plot \
'WholeFRAP.csv' u 1:2 w p ps 1.2 lw 2 pt 4 lc rgb "#444444" title 'experiment',\
'frap.out' u 1:2 w p ps 1.4 pt 8 lw 2 lc rgb "#FF4444" title 'simulation',\
x**0.25/25 w l lw 2 lc rgb "black" title '0.04*t^{1/4}'
