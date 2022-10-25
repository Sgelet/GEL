#!/bin/bash

# Variations to run
declare -a variations=(
  "-DMULTISCALE=0"
  ""
  "-DTHICC_SEP=1"
  "-DRECALC=2"
  "-DTHICC_SEP=1 -DRECALC=2"
)

# Build default
#/usr/bin/cmake --build ../cmake-build-default --target skeltal -j 6


for i in "${variations[@]}"
do
  # Build with compile time option
  cmake -G Ninja -DCMAKE_MAKE_PROGRAM=ninja -DCMAKE_CXX_FLAGS="$i -O2 -DCORE_TEST=16 -DCORE_TEST_SEC=2" -S .. -B ../cmake-build-default -DCMAKE_BUILD_TYPE=Release
  cmake --build ../cmake-build-default

  # Run
  ls -1 /work3/etoga/3DMeshes/[^0-9]* | ./runtime_test.sh 3 $i
done
