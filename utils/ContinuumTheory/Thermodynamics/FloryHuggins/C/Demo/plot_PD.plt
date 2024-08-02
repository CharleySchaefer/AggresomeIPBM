set terminal pngcairo enhanced
set output 'PD.png'

unset key

set xlabel 'Protein concentration, {/Symbol f}'
set ylabel 'Interaction Parameter, {/Symbol c}'

plot \
'PhaseDiagram_1.0.out' u 2:1 w l dt 2 lw 2 lc 'green',\
"" u 4:1 w l dt 2 lw 2 lc 'green',\
"" u 6:1 w l lw 2 lc 'blue',\
"" u 8:1 w l lw 2 lc 'blue',\
"<(sed -n '1,1p' PhaseDiagram_1.0.out)" u 2:1 w p pt 6 ps 2.2 lw 4 lc rgb 'red',\
'PhaseDiagram_2.0.out' u 2:1 w l dt 2 lw 2 lc 'green',\
"" u 4:1 w l dt 2 lw 2 lc 'green',\
"" u 6:1 w l lw 2 lc 'blue',\
"" u 8:1 w l lw 2 lc 'blue',\
"<(sed -n '1,1p' PhaseDiagram_2.0.out)" u 2:1 w p pt 6 ps 2.2 lw 4 lc rgb 'red',\
'PhaseDiagram_5.0.out' u 2:1 w l dt 2 lw 2 lc 'green',\
"" u 4:1 w l dt 2 lw 2 lc 'green',\
"" u 6:1 w l lw 2 lc 'blue',\
"" u 8:1 w l lw 2 lc 'blue',\
"<(sed -n '1,1p' PhaseDiagram_5.0.out)" u 2:1 w p pt 6 ps 2.2 lw 4 lc rgb 'red',\
'PhaseDiagram_10.0.out' u 2:1 w l dt 2 lw 2 lc 'green',\
"" u 4:1 w l dt 2 lw 2 lc 'green',\
"" u 6:1 w l lw 2 lc 'blue',\
"" u 8:1 w l lw 2 lc 'blue',\
"<(sed -n '1,1p' PhaseDiagram_10.0.out)" u 2:1 w p pt 6 ps 2.2 lw 4 lc rgb 'red',\
'PhaseDiagram_20.0.out' u 2:1 w l dt 2 lw 2 lc 'green',\
"" u 4:1 w l dt 2 lw 2 lc 'green',\
"" u 6:1 w l lw 2 lc 'blue',\
"" u 8:1 w l lw 2 lc 'blue',\
"<(sed -n '1,1p' PhaseDiagram_20.0.out)" u 2:1 w p pt 6 ps 2.2 lw 4 lc rgb 'red',\
'PhaseDiagram_50.0.out' u 2:1 w l dt 2  lw 2 lc 'green',\
"" u 4:1 w l dt 2 lw 2 lc 'green',\
"" u 6:1 w l lw 2 lc 'blue',\
"" u 8:1 w l lw 2 lc 'blue',\
"<(sed -n '1,1p' PhaseDiagram_50.0.out)" u 2:1 w p pt 6 ps 2.2 lw 4 lc rgb 'red'
