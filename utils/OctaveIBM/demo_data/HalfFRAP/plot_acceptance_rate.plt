set terminal pngcairo enhanced
set output 'acceptance_rate.png'

set log x
set xlabel "time [s]"
set xtics format "10^{%T}"
set ylabel "accepted processes [%]"
set yrange [0:100]

set grid lt 1 lw 0.5 lc rgb "grey"

unset key

plot 'timeprogress.dat' u 1:5 w p ps 2 pt 6
