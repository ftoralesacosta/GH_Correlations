#!/bin/bash 
mv *GeVTracks_Mixed_300_Correlation_Lambda.root $1
mv *_SE_L0_Correlation.root $1
cd ../Event_Mixing/
./divide_mix.sh ../InputData/$1/
mkdir ../Analysis_Notebooks/pics/LO/$1/
touch ../Analysis_Notebooks/pics/LO/$1/.gitignore

git pull
git add -f ../InputData/$1
git add ../Corr_config.yaml
git add ../Analysis_Notebooks/pics/LO/$1/
git commit -m "${1}"
git push
