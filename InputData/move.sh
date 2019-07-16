#!/bin/bash 
mv *GeVTracks_Mixed_300_Correlation_Lambda.root $1
mv *_SE_L0_Correlation.root $1
cd ../Event_Mixing/
./divide_mix.sh ../InputData/$1/

git pull
git add -f $1
git commit -m "${1}"
git push
