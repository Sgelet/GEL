#!/bin/bash

# Variations to run
declare -a variations=(
  "-DCORE_TEST=8"
)

variant=0

mkdir -p out_raw
mkdir -p skeletons

for i in "${variations[@]}"
do
  # Build with compile time option
  cmake -G Ninja -DCMAKE_MAKE_PROGRAM=ninja -DCMAKE_CXX_FLAGS="$i -O2 -DCORE_TEST_SEC=2" -S .. -B ../cmake-build-default -DCMAKE_BUILD_TYPE=Release
  cmake --build ../cmake-build-default

  # Make directories
  mkdir -p "skeletons/var$((++variant))"

  # Run
  ls -1 data/3DMeshes/asian* | ./runtime_test.sh 3 "var$variant"
done
