#!/bin/bash

# Build with compile time option
cmake -G Ninja -DCMAKE_MAKE_PROGRAM=ninja -DCMAKE_CXX_FLAGS="$i -O2 -DCORE_TEST_SEC=2" -S .. -B ../cmake-build-default -DCMAKE_BUILD_TYPE=Release
cmake --build ../cmake-build-default

#Gather files
mapfile -t files

length=${#files[@]}
for f in "${!files[@]}";
do
  full="${files[$f]##*/}"
  filename="${full%.*}"
  echo -n "${filename}"
  ../cmake-build-default/skeltal "${files[$f]}" "../results/skeletons/baseline/${filename}.skel" 1
done
