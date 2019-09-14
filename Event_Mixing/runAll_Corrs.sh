#!/bin/bash

date

module load root
module load hdf5

n_mix=300

for dataset in {13d,13e,13f,17q}

do
    
    for p in {0,4}

    do
	name=${dataset}
	echo $name
	
	if [ "17q" = "$name" ]; then #account for 17q-17p mixing
	    echo "./Parallel_Mix_Correlations ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/17p_MB_${p}GeV.hdf5 0 $n_mix $p"
	    ./Parallel_Mix_Correlations ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/17p_MB_${p}GeV.hdf5 0 $n_mix $p
	else
	    echo "./Parallel_Mix_Correlations ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/${name}_MB_${p}GeV.hdf5 0 $n_mix $p"
	    ./Parallel_Mix_Correlations ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/${name}_MB_${p}GeV.hdf5 0 $n_mix $p
	fi
	echo "$nmix Mixed Events DATASET =  $name"
	echo
    done
    
done

#mail -s "Mixing Done" ftoralesacosta@lbl.gov < ./
