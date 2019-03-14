#!/bin/bash

if [[ $# -eq 1 ]] ; then
    echo 'please give [{13d, 13e, 17q,...}.root] and [{___}_MB_{}GeV.hdf5] [test/full] as arguments'
    exit 0
fi

date

module add ROOT/6.08.00

for p in {0,4}
#for p in 0
do
    echo $1
    name=${1%.*}
    #name=$(basename ../$1 .root)
    echo $name
    if [ ! -f ${name}_${p}GeVTrack_paired.root ]; then
	echo "calling paired Injector with file $name.root and track Gev $p"
	./../pair_gale_shapley/paired_injector $name.root $p
    fi
    
    #for i in {0..4..2} #Mix 2 test events
    for i in {0..150..5} #Mix 300 events
    do
	mix_min=$i
	mix_max="$((i + 19))"
	if [[ $3 == full ]]; then
	    if [ "../InputData/17q" = "$name" ]; then #account for 17q-17p mixing
		sbatch -p shared-chos -t 20:00:00 HH_runCorr.sh ${name}_${p}GeVTrack_paired.root ${name::$((${#name} - 1))}p_MB_${p}GeV.hdf5 $mix_min $mix_max $p
	    else
		sbatch -p shared-chos -t 20:00:00 HH_runCorr.sh ${name}_${p}GeVTrack_paired.root ${name}_MB_${p}GeV.hdf5 $mix_min $mix_max $p	  
	    fi
	else
	    if [ "../InputData/17q" = "$name" ]; then
		./HH_runCorr.sh ${name}_${p}GeVTrack_paired.root ${name::$((${#name} - 1))}p_MB_${p}GeV.hdf5 $mix_min $mix_max $p
	    else
		./HH_runCorr.sh ${name}_${p}GeVTrack_paired.root ${name}_MB_${p}GeV.hdf5 $mix_min $mix_max $p
	    fi
	fi
    echo "$mix_min $mix_max $name"
    done
done
