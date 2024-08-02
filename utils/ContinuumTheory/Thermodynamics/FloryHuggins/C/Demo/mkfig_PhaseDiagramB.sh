#!/bin/bash

pltscript="
#==================================
# TERMINAL SETTINGS
#set terminal pngcairo enhanced dashed
#set output 'Fig_LvsR.png'

set terminal epslatex standalone #size 10 cm,  4 cm
set output 'tmp_figure.tex'
#==================================

set lmargin 8

set multiplot

#set origin 0,0.55
#set size 1,0.45

set xlabel '\huge Total Protein Concentration'
set ylabel '\huge \$\chi\$'

set xrange [0:1]
set yrange [0:2.5]
set ytics 0,0.5,4
set xtics 0,0.2,1

#set label '\Large \$N_\mathrm{P}=50\$' at 0.08,2.1
set label '\Large \\\\textcolor{blue}{\$\phi_\mathrm{P,A}(\chi_{PS})\$}' at 0.62,1.05 rotate by 35
set label '\Large \\\\textcolor{blue}{\$\phi_\mathrm{P,C}(\chi_{PS})\$}' at 0.06,1.05 
#set xtics format ''
set arrow from 0.03, 0.68 to 0.07,0.95 lc rgb 'blue'nohead front


unset key
#set key left bottom Left reverse
plot 'PhaseDiagram_50.0.out' u 2:1 w l dt 1  lw 12 lc 'green',\
'' u 4:1 w l dt 1 lw 12 lc 'green',\
'' u 6:1 w l lw 12 lc 'blue',\
'' u 8:1 w l lw 12 lc 'blue'#,\
#\"<(sed -n '1,1p' PhaseDiagram_50.0.out)\" u 2:1 w p pt 6 ps 3.2 lw 8 lc rgb 'red'

"



pltscript="$pltscript
unset multiplot
"

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
    \includegraphics*{Fig_cropped.eps}
\end{figure}


\end{document}
" > tmp_create_standalone.tex

rm -f tmp_create_standalone.pdf
printf "$pltscript"
printf "$pltscript" | gnuplot

latex tmp_figure.tex
dvips tmp_figure.dvi -o tmp_figure.ps
mv tmp_figure.ps "Fig_PhaseDiagramB.eps"

#pdflatex tmp_create_standalone.tex &>/dev/null
#mv tmp_create_standalone.pdf Fig_a4paper.pdf
rm tmp_* 



