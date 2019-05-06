# GH_Correlations


Azimuthal correlation of isolated, high-pT trigger photon. This analysis uses NTuples created by https://github.com/yslai/ntuple-gj . The full analysis procedure is:

1. Run same event correlations on a gamma-triggered NTupple
2. Run Event Mixing
   1. Run pairing between gamma triggered and min-bias NTupples
   2. Convert minimum bias NTupple to hdf5 format
   3. Run mixed event correlations between paired, gamma-triggered NTupple (output of step 2.i) and min-bias hdf5 file (output of step 2.ii)
   4. Repeat for track skimming for higher track-pt correlations with high statistics
   5. Divide Same event correlations by normalized mixed event correlations
3. Repeat steps 1-2 for all datasets individually
4. Combine correlations fo same collision system (pp, pPb, PbPb)
5. Analysis of correlations and fragmentation with Jupyter Notebooks
   1. Open pythonic_correlatinos.ipynb
   2. edit defaulty_values.py for desired settings
   3. Run notebook

The notebook goes throgh the following steps:
1. Uncorrelated background subtraction
2. Correlated Background subtraction (using purity based on https://github.berkeley.edu/alwina/photon-correlations)
3. Away side integration of correltation for fragmentation
4. Comparasion between collision systems