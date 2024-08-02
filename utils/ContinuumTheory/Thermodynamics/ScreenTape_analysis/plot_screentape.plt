set terminal pngcairo enhanced
set output 'Fig2c'

set xlabel 'RNA length in nucleotides'
set log x

set log y

set multiplot
set size 1,0.5
set origin 0,0.5
set ylabel 'Intensity'
plot \
"Figure_2c_screentape.csv" u 2:6:7 w e title 'Ars-CR',\
"" u 2:11:12 w e title 'control', \
"" u 2:16:17 w e title 'Ars-AR'

set origin 0,0
set ylabel 'Intensity ratio Ars-AR/Ars-CR'
unset key
set log y
set yrange [0:100]
plot \
"Figure_2c_screentape.csv" u 2:($16/$6):( ($16/$6)*sqrt( (($12+40)/$11)**2 + (($17+40)/$16)**2  ) ) w e,\
exp(x*0.01)/8
unset multiplot
