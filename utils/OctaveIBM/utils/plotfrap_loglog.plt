set terminal pngcairo enhanced
set output 'frap_loglog.png'

set log xy
set xlabel "time [s]"
set xtics format "10^{%T}"
set ylabel "frap [%]"

set grid lt 1 lw 0.5 lc rgb "grey"

unset key

plot 'frap.out' u 1:2 w p ps 2 pt 6,\
x**0.5,\
x**0.25
