#import sys
import root_numpy
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

import pandas as pd
import root_pandas as rpd
import ROOT
from default_values import *

def Get_NTriggers_2D(filename, ipt, Signal_DNN=True): 
    file = ROOT.TFile(filename)
    
    if (Signal_DNN == "NoIso_NoShape"):
        ntrig_histo = file.Get("N_Inclusive_Triggers_pT%1.0f_%1.0f" %(pTbins[ipt],pTbins[ipt+1]))
    
    elif (Signal_DNN == "Inclusive"):
        ntrig_histo = file.Get('N_Triggers_pT%1.0f_%1.0f' %(pTbins[ipt],pTbins[ipt+1]))
    
    else:
        DNN_Rgn = int(Signal_DNN) + 2*(1-int(Signal_DNN))
        ntrig_histo = file.Get('N_DNN%i_Triggers_pT%1.0f_%1.0f' %(DNN_Rgn,pTbins[ipt],pTbins[ipt+1]))
    
    NTriggers = 1
    if not(ntrig_histo == None):
        NTriggers = ntrig_histo.GetEntries()
    file.Close()
    return NTriggers

def Get_2D_Correlation(filename, ipt, izt,ntriggers,Signal_DNN=True):     #returns TH2D

    file = ROOT.TFile(filename)    
    
    if (Signal_DNN == "NoIso_NoShape"):
        histo2D = file.Get('Inclusive_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                  %(pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]))
                               
    elif (Signal_DNN == "Inclusive"):
        histo2D = file.Get('Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                  %(pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]))

    else:
        DNN_Rgn = int(Signal_DNN) + 2*(1-int(Signal_DNN)) #convert bool to DNN_1 (Sgn) or DNN_2 (Bkgd)
        histo2D = file.Get('DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                           %(DNN_Rgn,pTbins[ipt], pTbins[ipt+1], 100*zTbins[izt],100*zTbins[izt+1]))

    if not(ntriggers == None):
        histo2D.Scale(1.0/ntriggers)
    return histo2D


def GetPhiProj(histo2D,ntriggers):
            
    Eta_Axis = histo2D.GetYaxis()
    PhiProjection = histo2D.ProjectionX('DNN%i_PhiProjection__pT_%1.0f_%1.0f__zt_%1.0f_%1.0f' 
                            %(DNN_Rgn,pT_Bins[ipt],pT_Bins[ipt+1],100*zT_Bins[izt],
                            100*zTbins[izt+1]),Eta_Axis.FindBin(-eta_max),Eta_Axis.FindBin(eta_max))
                                                
    PhiProjection.SetDirectory(0)
    PhiProjection.Rebin(2)
    PhiProjection.Scale(1.0/(2*eta_max))
    
    #per trigger yield
    if not(ntriggers == None):
        PhiProjection.Scale(1.0/ntriggers)
    
    file.Close()
    return PhiProjection

def Get_Phi_Arrays(PhiProjection):
  
    Phi_Array = np.zeros(len(delta_phi_centers))
    Phi_Error_Array = np.zeros(len(delta_phi_centers))
    for bin in range(2,9):
        Phi_Array[bin-2] = PhiProjection.GetBinContent(bin)
        Phi_Error_Array[bin-2] = PhiProjection.GetBinError(bin)
     
    #Phi_Array = root_numpy.hist2array(PhiProjection) 
    return Phi_Array, Phi_Error_Array


def Get_N_2D_Pairs(histo2D):
    N_Pairs = histo2D.GetEntries()
    return N_Pairs

def Get_N_Phi_Pairs(PhiProjection):
    N_phi_pairs = PhiProjection.GetEntries()
    return N_phi_pairs

#Define function that makes array of N_phi pairs by looping over ipt and izt
# def Plot_pT_zT_bins(function):
#     N_Trig_Array = np.zeros(NpT_Bins)
#     N_Trig_Array_Scales = np.zeros(NpT_Bins)
#     for ipt in range (NpT_Bins):
#         N_Phi_Pair_Array = np.zeros((NpT_Bins,NzT))
#         N_2D_Pair_Array = np.zeros((NpT_Bins,NzT))

#         Signal_Triggers = Get_NTriggers(filename,ipt,True)
#         Background_Triggers = Get_NTriggers(filename,ipt,False)

#         for izt in range (NzT):
#             Sig2D = Get_2D_Correlation(filename, ipt, izt, True)
#             Bgk2D = Get_2D_Correlation(filename, ipt, izt, False)

#Define more generally root2numpy py file, where the functions loop over zT and pT. This reduces I/O substantially, and helps run 
