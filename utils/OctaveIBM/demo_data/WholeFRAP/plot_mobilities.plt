set terminal pngcairo enhanced
set output "mobilities.png"

set log x
set xlabel "D[{/Symbol m}m^2/s]"
set xtics format '10^{%T}'
set log y

plot "mobilities.dat" u (10**($1)):2,\
"" u (10**($3)):4,\
