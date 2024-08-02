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




set xlabel '\Large \$\chi_\mathrm{PS}\$'
set ylabel '\Large fraction of RNA in aggresome'
set yrange [0:5]


VA_by_VAmin(phiA,phiC, phi0)=((1/phi0)*( phi0-phiC )/(phiA-phiC))

set key right bottom

set label '\Large \$N_\mathrm{P}=50\$' at screen 0.75,0.4
set label '\Large \$\phi_\mathrm{P}=0.2\$' at screen 0.75,0.45
set yrange [0.5:1]
plot 'PhaseDiagram_50.0.out' u 1:(1.0/(1+exp(-(\$8-\$6)*0.1*10))) w l lw 4 title '\$ N_\mathrm{R}(\chi_{RS}-\chi_{RP})=1\$',\
'' u 1:(1.0/(1+exp(-(\$8-\$6)*0.1*20))) w l lw 4 title '\$ N_\mathrm{R}(\chi_{RS}-\chi_{RP})=2\$',\
'' u 1:(1.0/(1+exp(-(\$8-\$6)*0.1*50))) w l lw 4 title '\$ N_\mathrm{R}(\chi_{RS}-\chi_{RP})=5\$',\
'' u 1:(1.0/(1+exp(-(\$8-\$6)*0.1*100))) w l lw 4 title '\$ N_\mathrm{R}(\chi_{RS}-\chi_{RP})=10\$',\
'' u 1:(1.0/(1+exp(-(\$8-\$6)*0.1*200))) w l lw 4 title '\$ N_\mathrm{R}(\chi_{RS}-\chi_{RP})=20\$'
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
mv tmp_figure.ps "Fig_RNA_partitioning.eps"

#pdflatex tmp_create_standalone.tex &>/dev/null
#mv tmp_create_standalone.pdf Fig_a4paper.pdf
rm tmp_* 



