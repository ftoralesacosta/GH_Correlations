#!/bin/bash

date

module load root

n_mix=300

for dataset in {13f,13d,13e,17q} #13f first since troublesome (only uses half block size)

do
    
    for p in {0,4}

    do
	name=${dataset}
	echo $name
	
	if [ "17q" = "$name" ]; then #account for 17q-17p mixing
	    echo "./Parallel_Mix_Correlations_pWeight ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/17p_MB_${p}GeV.hdf5 0 $n_mix $p pp"
	    ./Parallel_Mix_Correlations_pWeight ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/17p_MB_${p}GeV.hdf5 0 $n_mix $p pp

	else
	    echo "./Parallel_Mix_Correlations_pWeight ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/${name}_MB_${p}GeV.hdf5 0 $n_mix $p pPb"
	    ./Parallel_Mix_Correlations_pWeight ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/${name}_MB_${p}GeV.hdf5 0 $n_mix $p pPb
	fi
	echo "$nmix Mixed Events DATASET =  $name"
	echo
    done
    
done

#mail -s "Mixing Done" ftoralesacosta@lbl.gov < ./
