import numpy as np
import math

        #DEFAULTS:


Shower = "LO"

pPb_File = 'InputData/pPb_SE_L0_Correlation_GMB_Ratio_Track.root'
pp_File = 'InputData/pp_SE_L0_Correlation_GMB_Ratio_Track.root'

Systems = ["pp","p-Pb"]
Files = [pp_File,pPb_File]

purity = np.asarray([0.208095, 0.339478, 0.483944, 0.509])
purity_Uncertainty = np.asarray([0.0273652,0.0305499,0.0358471,0.036])

Uncorr_Estimate = "ZYAM"
Average_UE = False

Corrections = np.asarray([1.007,0.982,0.957,0.926,0.894,0.853,0.817,0.757,0.681,0.673,0.619,0.469,0.342,0.301])
Fake_Rate = np.asarray([0.0179,0.0197,0.0249,0.0355,0.0475,0.0722,0.0902,0.1298,0.1407,0.2130,0.2175,0.2376,0.2611,0.2611])
oneminFake = np.ones(len(Fake_Rate))-Fake_Rate


        # RANGE & BINNING:

#pT
pTbins = [12,15,19,26,40]
N_pT_Bins = len(pTbins)-1

#zT
zTbins = np.asarray([0.05, 0.07670637, 0.11767734, 0.18053204, 0.27695915, 0.42489062, 0.65183634]) #0.65-1 skipped for now
zT_centers = (zTbins[1:] + zTbins[:-1]) / 2
zT_widths = [(j-i)/2 for i, j in zip(zTbins[:-1], zTbins[1:])]
zT_offset = 0
NzT = len(zTbins)-1
zt_box = np.ones(NzT) * 0.03 #plotting Uncert. Boxes

#deta
eta_max = 1.2 #Range of Signal Correlations

#dPhi
N_dPhi_Bins = 8
dPhi_Bins = [i*math.pi/N_dPhi_Bins for i in range(0,8)]
delta_phi_centers = [i*math.pi/N_dPhi_Bins+math.pi/N_dPhi_Bins/2 for i in range(1,N_dPhi_Bins)] #skip first dPhi bin to avoid Isolation

#dPhi_Integration
N_Phi_Integrate = 3 #Number of dPhi Bins for away-side integration. 3 Corresponds to dphi > 2.1
Integration_Width = math.pi/(len(delta_phi_centers)+1) * N_Phi_Integrate
phi_width = math.pi/(N_dPhi_Bins)/2

#UE
ue_error_bar = [dPhi_Bins[1],dPhi_Bins[2]] #Horiz. width of UE at first plotted dphi point


        #CASE SWITCHING

#Shower = "NN"
Shower = "LO"
Use_Weights = True
CorrectedP = True  # FALSE FOR HARDPROBES
Use_MC = False
pT_Rebin = True

if (Shower == "LO"):
    pPb_File = 'InputData/pPb_SE_L0_Correlation_GMB_Ratio_Track.root'
    pp_File = 'InputData/pp_SE_L0_Correlation_GMB_Ratio_Track.root'

    if (CorrectedP):
        purity = np.asarray([0.208095, 0.339478, 0.483944, 0.509])
        purity_Uncertainty = np.asarray([0.0273652,0.0305499,0.0358471,0.036])
    else:
        purity = [0.35]
        
    if (pT_Rebin):
        pPb_File = 'InputData/pPb_SE_L0_Correlation_GMB_Ratio_pTRebin.root'
        pp_File = 'InputData/pp_SE_L0_Correlation_GMB_Ratio_pTRebin.root'
        
        purity = np.asarray([0.260799,0.511580])
        purity_Uncertainty = np.asarray([0.0287628,0.0362961])
        pTbins = [12,21,40]
        N_pT_Bins = len(pTbins)-1
        
        Files = [pp_File,pPb_File]

if (Shower == "NN"):
        pPb_File = 'InputData/pPb_SE_NN_Correlation_GMB_Ratio_Track.root'
        pp_File = 'InputData/pp_SE_NN_Correlation_GMB_Ratio_Track.root'
                
        if (CorrectedP):
            purity = np.asarray([0.238477, 0.341009, 0.479701, 0.532013])
            purity_Uncertainty = np.asarray([0.09,0.1,0.11,0.11])
        else:
            purity = [0.352546]

MC_File = 'InputData/18b10a_pthat_1_2_SE_NN_Correlation_GMB_Ratio.root'

if(Use_MC):
    Systems = ["pp","p-Pb","MC"]
    Files = [pp_File,pPb_File,MC_File]
    
print Files
