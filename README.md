# GH_Correlations


Azimutha correlation of isolated, high-pT trigger photon. This analysis uses NTuples created by https://github.com/yslai/ntuple-gj . The full analysis procedure is:

1. Run same event correlations on a gamma-triggered NTupple
2. Run Event Mixing
   1. Run pairing between gamma triggered 
   2. Convert minimum bias NTupple to hdf5 format
   3. Run mixed event correlatinos between paired,gamma-triggered NTupple (output of step 2.1) and min-bias hdf5 file (output of step 2.2)
   4. Repeat for track skimming for higher track-pt correlations with high statistics
   5. Divide Same event correlations by normalized mixed event correlations
3. Repeat steps 1-2 for all datasets individually
4. Combine correlations fo same collision system (pp, pPb, PbPb)
5. Analysis of correlations and fragmentation with Jupyter Notebooks