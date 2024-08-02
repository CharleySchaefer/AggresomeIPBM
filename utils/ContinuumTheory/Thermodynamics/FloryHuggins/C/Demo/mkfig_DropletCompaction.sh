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


set xlabel '\Large Protein-Solvent interaction, \$\chi_\mathrm{PS}\$'
set ylabel '\Large \$V_\mathrm{A}/V_\mathrm{A,min}\$'
set yrange [0:5]


VA_by_VAmin(phiA,phiC, phi0)=((1/phi0)*( phi0-phiC )/(phiA-phiC))



set label '\Large \$N_\mathrm{P}=50\$' at screen 0.7,0.9
#set xtics format ''

set key at 3.8,4.2
#set key left bottom Left reverse
plot \
'PhaseDiagram_50.0.out' u 1:(VA_by_VAmin(\$8, \$6, 0.01)) w l dt 1  lw 6 title '\$\phi_\mathrm{P}=0.01\$',\
'PhaseDiagram_50.0.out' u 1:(VA_by_VAmin(\$8, \$6, 0.02)) w l dt 1  lw 6   title '\$\phi_\mathrm{P}=0.02\$',\
'PhaseDiagram_50.0.out' u 1:(VA_by_VAmin(\$8, \$6, 0.05)) w l dt 1  lw 6  title '\$\phi_\mathrm{P}=0.05\$',\
'PhaseDiagram_50.0.out' u 1:(VA_by_VAmin(\$8, \$6, 0.10)) w l dt 1  lw 6 title '\$\phi_\mathrm{P}=0.10\$',\
'PhaseDiagram_50.0.out' u 1:(VA_by_VAmin(\$8, \$6, 0.20)) w l dt 1  lw 6  title '\$\phi_\mathrm{P}=0.20\$'

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
mv tmp_figure.ps "Fig_DropletCompaction.eps"

#pdflatex tmp_create_standalone.tex &>/dev/null
#mv tmp_create_standalone.pdf Fig_a4paper.pdf
rm tmp_* 



