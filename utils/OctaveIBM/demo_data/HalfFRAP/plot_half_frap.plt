set terminal pngcairo enhanced
set output 'halffrap.png'

#set log xy
set xlabel "time [s]"
#set xtics format "10^{%T}"
set xrange [0:600]
set ylabel "frap [%]"
set yrange [0:1]

set grid lt 1 lw 0.5 lc rgb "grey"



plot \
'HalfFRAP.csv' u 1:2 w p ps 1.2 lw 2 pt 4 lc rgb "#444444" title 'experiment',\
'frap.out' u 1:2 w p ps 1.4 pt 8 lw 2 lc rgb "#FF4444" title 'simulation',\
x**0.5/40 w l lw 2 lc rgb "black" title '0.025*t^{1/2}'
