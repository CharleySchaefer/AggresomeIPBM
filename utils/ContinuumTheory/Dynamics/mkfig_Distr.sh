#!/bin/bash

pltscript="
#==================================
# TERMINAL SETTINGS
#set terminal pngcairo enhanced dashed
#set output 'Fig_LvsR.png'

set terminal epslatex standalone #size 10 cm,  4 cm
set output 'tmp_figure.tex'
#==================================

set lmargin 9
set bmargin 4.5

#set multiplot

#set origin 0,0.55
#set size 1,0.45

#set xlabel '\huge time \$\\\\times K_\mathrm{deg}\$' offset 0,-0.5
#set ylabel '\huge Mean RNA length' offset -0.5,0

set log x
#set log y
set xrange [6:1000]
#set yrange [0.5e-1:1e2]
set xtics format '\\Large \$10^{%%T}\$' 
#set ytics format '\\Large \$10^{%%T}\$'

#set xrange [0:1]
#set yrange [0:2.5]
#set ytics 0,0.5,4
#set xtics 0,0.2,1

#set label '\Large \$N_\mathrm{P}=50\$' at 0.08,2.1
#set label '\Large \\\\textcolor{blue}{\$\phi_\mathrm{P,A}(\chi_{PS})\$}' at 0.62,1.05 rotate by 35
#set label '\Large \\\\textcolor{blue}{\$\phi_\mathrm{P,C}(\chi_{PS})\$}' at 0.06,1.05 
#set xtics format ''
#set arrow from 0.03, 0.68 to 0.07,0.95 lc rgb 'blue'nohead front


set yrange [0:5000]
set ytics 0,1000,5000

#dir='../FreeDiffusion_LatticeCube_S5_lK1.00/Nbound3'
set key right top #Left reverse spacing 1.5
#set grid
set multiplot
set size 1,0.5

set label '\$K_\mathrm{diff}/K_\mathrm{deg}=161\$' at screen 0.2,0.9

unset key
set origin 0,0.4
set ylabel '\\\\Large \$I_\mathrm{C}\$'
plot \
'Distr161.txt' u (6*\$0):(\$3*1e8) w l dt 1 lw 3 lc rgb '#000066' title 'early',\
'Distr161.txt' u (6*\$0):(\$17*1e8) w l dt 1 lw 3 lc rgb '#000077' title 'mid',\
'Distr161.txt' u (6*\$0):(\$21*1e8) w l dt 1 lw 3 lc rgb '#000088' title 'mid',\
'Distr161.txt' u (6*\$0):(\$25*1e8) w l dt 1 lw 3 lc rgb '#000099' title 'mid',\
'Distr161.txt' u (6*\$0):(\$29*1e8) w l dt 1 lw 3 lc rgb '#0000BB' title 'mid',\
'Distr161.txt' u (6*\$0):(\$33*1e8) w l dt 1 lw 3 lc rgb '#0000DD' title 'late',\
'Distr161.txt' u (6*\$0):(\$35*1e8) w l dt 1 lw 3 lc rgb '#0000FF' title 'late'


set origin 0,0.0
set xlabel '\\\\Large Nucleotides, \$N_\mathrm{R}\$'
set ylabel '\\\\Large \$I_\mathrm{A}\$'
plot \
'Distr161.txt' u (6*\$0):(\$4*1e8) w l dt 1 lw 3 lc rgb '#660000' title 'early',\
'Distr161.txt' u (6*\$0):(\$18*1e8) w l dt 1 lw 3 lc rgb '#990000' title 'mid',\
'Distr161.txt' u (6*\$0):(\$22*1e8) w l dt 1 lw 3 lc rgb '#990000' title 'mid',\
'Distr161.txt' u (6*\$0):(\$26*1e8) w l dt 1 lw 3 lc rgb '#990000' title 'mid',\
'Distr161.txt' u (6*\$0):(\$30*1e8) w l dt 1 lw 3 lc rgb '#990000' title 'mid',\
'Distr161.txt' u (6*\$0):(\$34*1e8) w l dt 1 lw 3 lc rgb '#FF0000' title 'late',\
'Distr161.txt' u (6*\$0):(\$36*1e8) w l dt 1 lw 3 lc rgb '#FF0000' title 'late',\


unset multiplot


"



#pltscript="$pltscript
#unset multiplot
#"

echo "% \documentclass[a5paper, landscape, 10 pt]{article}
\documentclass[a4paper, 10 pt]{article}
\usepackage[pdftex]{graphicx}
\usepackage{epstopdf}			% Use eps figures

\usepackage[english]{babel}
\usepackage[latin1]{inputenc}
\usepackage{a4wide}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{amsthm}
\usepackage{ctable}
\usepackage{wrapfig}
\usepackage{graphicx}
\usepackage{amsfonts}
\usepackage[small,bf]{caption}
\usepackage{cite}
\usepackage[T1]{fontenc}
\usepackage{ae,aecompl}
%\usepackage{eurosym}	% Euro symbol
\usepackage{url} 		% URL formatting package
\usepackage{hyperref}
\usepackage{float}
\usepackage{varioref}
\usepackage{ifthen}
\usepackage{color}
%\usepackage{times}

%\setlength{\parindent}{0cm}

%\bibliographystyle{unsrt} % Style for the References Section.
\bibliographystyle{unsrturl}
%\bibliographystyle{nprlrtn}


\begin{document}

\begin{figure}[h]
    \includegraphics*{Fig_Distr.eps}
\end{figure}


\end{document}
" > tmp_create_standalone.tex

rm -f tmp_create_standalone.pdf
printf "$pltscript"
printf "$pltscript" | gnuplot

latex tmp_figure.tex
dvips tmp_figure.dvi -o tmp_figure.ps
mv tmp_figure.ps "Fig_Distr.eps"

pdflatex tmp_create_standalone.tex #&>/dev/null
mv tmp_create_standalone.pdf Fig_RNA_length.pdf
rm tmp_*  Fig_RNA_length.pdf



