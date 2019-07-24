#!/bin/bash

source /Users/fernando/.bash_profile
make Draw2D

#file=../InputData/zT_Rebin_8_006zT06zT/pPb_SE_L0_Correlation.root
#file=../InputData/zT_Rebin_8_006zT06zT/pPb_0GeVTracks_ME.root
file=../InputData/zT_Rebin_8_006zT06zT/pPb_SE_L0_Correlation_GMB_Ratio.root

./Draw2D $file $1
