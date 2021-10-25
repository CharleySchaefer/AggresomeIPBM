#!/bin/bash
executable=IBM
src="../src"
lib="../lib"
utils="../utils"

CC=gcc # gcc or g++

echo "BUILD LIBRARY"
pushd $lib/ZiltoidLIB >/dev/null
  ./makeZiltoidLIB.sh --$CC
popd >/dev/null

echo "COMPILE $executable"
#gcc -o $executable $src/main.c $src/GenerateEnergyLandscape.c -L$lib/ZiltoidLIB -lZiltoidLIB -lm $1 
${CC}  -O3 -msse2 -DHAVE_SSE2 -DSFMT_MEXP=607 -Wall -o $executable  $lib/SFMT/SFMT.c $src/main.c $src/GenerateEnergyLandscape.c -L$lib/ZiltoidLIB -lZiltoidLIB -lm $1  


echo "COMPILE POSTPROCESSING TOOLS: StructureFactor"
SFdir=$lib/ZiltoidLIB/Applications/StructureFactor/Demo/
pushd $SFdir >/dev/null
  ./compile_StructureFactor.sh --$CC
popd >/dev/null
cp $SFdir/StructureFactorAndCorrelationFunction.o $utils

echo "USAGE"
echo "  To view program arguments, enter './$executable --help'"
echo "  For examples, copy demos to the build directory: 'cp ../doc/demo*/demo*.sh .'"
