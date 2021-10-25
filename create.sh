#!/bin/bash

mkdir -p build
pushd build  > /dev/null
  cp ../src/compile.sh .
  echo "COMPILE"
  ./compile.sh
  cp ../doc/demo*/demo*.sh .
  echo "DEMOS ARE COPIED TO BUILD"
popd > /dev/null
