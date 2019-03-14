#!/bin/bash

for m in {0..280..20}
do
    mplus="$((m+19))"
    echo "${mplus}"
    mv *_*GeVTrack_paired_*GeVTracks_Correlation_DNN_${m}_to_${mplus}.root "./inclusive/"
done