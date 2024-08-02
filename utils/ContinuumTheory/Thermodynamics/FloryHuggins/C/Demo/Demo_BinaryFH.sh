#!/bin/bash
LIB=../../lib
executeable=Demo_BinaryFH.o
tag="_test"
do_compile=1;

CHI=1.47
N=50.0
if [[ do_compile -eq 1 ]] ; then
echo "  Compiling"
if gcc -o  $executeable Demo_BinaryFH.c ../lib/Polynomial/Polynomial.c -lm ;
then
  echo "  Compiled $executeable"
else
  echo "Error: Failed to compiled $executeable"
  exit 1
fi
fi

CHICR=`awk -vN=$N  'BEGIN {print 0.5*(1+1.0/sqrt(N))**2}'`
PHICR=`awk -vN=$N  'BEGIN {print 1.0/(1+sqrt(N))}'`

> PhaseDiagram_$N.out
i=0
for chi in `seq $CHICR 0.01 4` ; do 
  i=$(( i + 1 ))
if [[ $i -eq 1 ]] ; then
echo "TEST"
  echo $CHICR" "$PHICR" "$PHICR" "$PHICR" "$PHICR" "$PHICR" "$PHICR" "$PHICR" "$PHICR   >> PhaseDiagram_$N.out
else

echo $chi
  ./$executeable $chi $N  > Demo_BinaryFH.out
  awk -vChi=$chi <"Demo_BinaryFH.out" '$1=="#phi_sp1:" {printf Chi" "$2"  "$3" "}' >> PhaseDiagram_$N.out
  awk <"Demo_BinaryFH.out" '$1=="#phi_sp2:" {printf $2"  "$3" "}' >> PhaseDiagram_$N.out
  awk <"Demo_BinaryFH.out" '$1=="#phi_bn1:" {printf $2"  "$3" "}' >> PhaseDiagram_$N.out
  awk <"Demo_BinaryFH.out" '$1=="#phi_bn2:" {printf $2"  "$3"\n"}' >> PhaseDiagram_$N.out
fi

done




if [[ 0 -eq 1 ]] ; then
echo "  Running"
./$executeable $CHI $N  > Demo_BinaryFH.out

> PhaseDiagram_$N.out
awk -vChi=$CHI <"Demo_BinaryFH.out" '$1=="#phi_sp1:" {print Chi" "$2"  "$3" "}' >> Demo_BinarySpinodal.out
awk <"Demo_BinaryFH.out" '$1=="#phi_sp2:" {printf $2"  "$3" "}' >> PhaseDiagram_$N.out
awk <"Demo_BinaryFH.out" '$1=="#phi_bn1:" {printf $2"  "$3" "}' >> PhaseDiagram_$N.out
awk <"Demo_BinaryFH.out" '$1=="#phi_bn2:" {printf $2"  "$3"\n"}' >> PhaseDiagram_$N.out

echo "  Extracting data from output"
> Demo_BinarySpinodal.out
awk <"Demo_BinaryFH.out" '$1=="#phi_sp1:" {print $2"  "$3" "}' >> Demo_BinarySpinodal.out
awk <"Demo_BinaryFH.out" '$1=="#phi_sp2:" {print $2"  "$3}' >> Demo_BinarySpinodal.out

> Demo_BinaryBinodal.out
awk <"Demo_BinaryFH.out" '$1=="#phi_bn1:" {print $2"  "$3" "}' >> Demo_BinaryBinodal.out
awk <"Demo_BinaryFH.out" '$1=="#phi_bn2:" {print $2"  "$3}' >> Demo_BinaryBinodal.out


echo "  Plotting"
plot="set terminal pngcairo enhanced dashed
set output \"Demo_BinaryFH.png\"

set xlabel \"phi\"
set xrange [0:1]

set ylabel \"f\"

set key right bottom

plot \"Demo_BinaryFH.out\" u 1:2 w l lt 1 lw 2 lc rgb \"#AA4444\" title \"free energy\",\
\"Demo_BinaryBinodal.out\" u 1:2 w l lt 2 lw 1 lc rgb \"#4444AA\" notitle,\
\"Demo_BinarySpinodal.out\" u 1:2 w p pt 6 ps 2 lw 2  lc rgb \"#44AA44\"  title \"spinodal\",\
\"Demo_BinaryBinodal.out\" u 1:2 w p pt 6 ps 2 lw 2  lc rgb \"#4444AA\"  title \"binodal\"
"
printf "$plot" | gnuplot

fi
