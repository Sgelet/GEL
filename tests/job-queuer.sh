#!/bin/bash

#Gather files
mapfile -t files

length=${#files[@]}
for f in "${!files[@]}";
do
  full="${files[$f]##*/}"
  filename="${full%.*}"
  sed "s:FILENAME:"${files[$f]}":g" template.sh | bsub
done
echo
