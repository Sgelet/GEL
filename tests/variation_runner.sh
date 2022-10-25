#!/bin/bash

# Variations to run
declare -a variations=(
  ""
  "-DTHICC_SEP=1"
  "-DRECALC=2"
  "-DTHICC_SEP=1 -DRECALC=2"
)

variant=0

mkdir -p out_raw
mkdir -p skeletons

for i in "${variations[@]}"
do
  # Build with compile time option
  cmake -G Ninja -DCMAKE_MAKE_PROGRAM=ninja -DCMAKE_CXX_FLAGS="$i -O2 -DCORE_TEST=16 -DCORE_TEST_SEC=2" -S .. -B ../cmake-build-default -DCMAKE_BUILD_TYPE=Release
  cmake --build ../cmake-build-default

  # Make directories
  mkdir -p "skeletons/var$((++variant))"

  # Run
  ls -1 /work3/etoga/3DMeshes/[^0-9]* | ./runtime_test.sh 3 "var$variant"
done
