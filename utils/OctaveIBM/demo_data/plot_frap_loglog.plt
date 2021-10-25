set terminal pngcairo enhanced
set output 'frap_loglog.png'

set log xy
set xlabel "time [s]"
set xrange [1:1000]
set xtics format "10^{%T}"
set ylabel "frap [%]"

set grid lt 1 lw 0.5 lc rgb "grey"

set key left top

plot \
'HalfFRAP/HalfFRAP.csv' u 1:2 w p ps 1.2 lw 2 pt 4 lc rgb "#884444" title 'Half FRAP experiment',\
'WholeFRAP/WholeFRAP.csv' u 1:2 w p ps 1.2 lw 2 pt 6 lc rgb "#444488" title 'Whole FRAP experiment',\
'HalfFRAP/frap.out' u (($1>7)&&($1<300)?$1:1/0):2 w p ps 1.2 pt 9 lw 2 lc rgb "#FF4444" title 'Half FRAP simulation',\
'WholeFRAP/frap.out' u ($1>7?$1:1/0):2 w p ps 1.2 pt 11 lw 2 lc rgb "#4444FF" title 'Whole FRAP simulation',\
x**0.25/25 w l lw 2 lc rgb "blue" title '\~t^{1/4}',\
(x<300?x**0.5/40:1/0) w l lw 2 lc rgb "red" title '\~t^{1/2}'
