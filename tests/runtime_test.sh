#!/bin/bash

#Gather files
mapfile -t files

length=${#files[@]}
for f in "${!files[@]}";
do
  full="${files[$f]##*/}"
  filename="${full%.*}"
  for (( i=0 ; i<$1 ; i++ ));
  do
    echo "Running ${files[$f]} ${2:+skeletons/$2/$filename.skel} $(($i + 1)) of $1"
    #../cmake-build-default/skeltal "${files[$f]}" "${2:+skeletons/$2/$filename.skel}" >> "out_raw/$filename.log"
    #./progresser.sh $((f*$1+i+1)) $((length*$1)) 50
  done
done
echo
