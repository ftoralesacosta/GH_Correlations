#!/bin/bash

date

module load root
export HDF5_USE_FILE_LOCKING=FALSE
n_mix=300

#for dataset in {13f_new_skimClusterMinE12,13}
#dataset="13f_new_skimClusterMinE12"
	       #for dataset in {13f_new_skimClusterMinE12,17q}
for p in {0,}
		 
do
    name=${1}
    echo $name
    
    if [ "17q" = "$name" ]; then #account for 17q-17p mixing
	echo "./parallel_pair_gale_shapley_hdf5 ../InputData/${name}.root ../InputData/17p_minbias_${p}GeVtracks.hdf5 0 $n_mix $p"
	./parallel_pair_gale_shapley_hdf5 ../InputData/${name}.root ../InputData/17p_minbias_${p}GeVtracks.hdf5 0 $n_mix $p
	
    elif [ "13fnew_skimClusterMinE12" = "$name" ];
    then
	echo "./parallel_pair_gale_shapley_hdf5 ../InputData/${name}.root ../InputData/13f_minbias_${p}GeVtracks.hdf5 0 $n_mix $p"
	./parallel_pair_gale_shapley_hdf5 ../InputData/${name}.root ../InputData/13f_minbias_${p}GeVtracks.hdf5 0 $n_mix $p
	
    else
	echo "./parallel_pair_gale_shapley_hdf5 ../InputData/${name}.root ../InputData/${name}_minbias_${p}GeVtracks.hdf5 0 $n_mix $p"
	./parallel_pair_gale_shapley_hdf5 ../InputData/${name}.root ../InputData/${name}_minbias_${p}GeVtracks.hdf5 0 $n_mix $p
    fi
    echo "$nmix Mixed Events DATASET =  $name"
    echo
done

mail -s "Pairing ${name} Done" ftoralesacosta@lbl.gov < ./
