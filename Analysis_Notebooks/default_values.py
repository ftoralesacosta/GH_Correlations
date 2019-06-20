import numpy as np
import math
import ROOT
import scipy.stats

from variations import *


                 #####CASE SWITCHING#####

#Shower = "NN"
Shower = "LO"
Use_Weights = True #BKR pT Weight
Use_P_Weights = True
CorrectedP = True  # FALSE FOR HARDPROBES
Use_MC = False
pT_Rebin = False
N_dPhi_Bins = 16
Ped_Sub_First = False #Important after weight implementation
Average_UE = False
Show_Fits = True
Uncorr_Estimate = "ZYAM"

#rad_start = 1.6
#Phi_String = "\pi/2"

rad_start = 1.9
Phi_String = "5\pi/8"

#rad_start = 2.7
#Phi_String = "7\pi/8"

        #DEFAULTS:

Shower = "LO"

#description_string="05zT"
#description_string="05zT_Remake"
#description_string="05zT_AliUSA"
#description_string="05zT_3bins"
#description_string="05zT_2bins"

#description_string="05zT_working_old"
#description_string="1zT"
#description_string="15zT"
#description_string="2zT"

#description_string="pT_Rebin_1"
#description_string="pT_Rebin_1_15pT"
#description_string="pT_Rebin_1_20pT"
#description_string = "pT_Rebin_1_pDevPlus"
#description_string = "pT_Rebin_1_pDevMinus"
#description_string = "pT_Rebin_1_pDevNONE"
#description_string = "pT_Rebin_1_90p"
#description_string= "pT_Rebin_1_8dPhi"
#description_string= "pT_Rebin_1_16dPhi"

#description_string = "pT_Rebin_2"


#description_string="pT_Rebin_3"

#description_string="pT_Rebin_4"
#description_string="pT_Rebin_4_Lambda"
#description_string="pT_Rebin_4_Weights"
#description_string="pT_Rebin_4_05zT"
#description_string="pT_Rebin_4_Cut" #15-40eV
#description_string= "pT_Rebin_4_ErrWeights"
#description_string= "pT_Rebin_4_ErrWeights_16dPhi"


#description_string = "zT_Rebin_6"
#description_string = "zT_Rebin_8"
#description_string = "zT_Rebin_9"
description_string = "zT_Rebin_14"



pPb_File = '../InputData/%s/pPb_SE_L0_Correlation_GMB_Ratio.root'%(description_string)
pp_File = '../InputData/%s/pp_SE_L0_Correlation_GMB_Ratio.root'%(description_string)
print pp_File
Systems = ["pp","p-Pb"]
Files = [pp_File,pPb_File]

purity = np.asarray([0.208095, 0.339478, 0.483944, 0.509])
purity_Uncertainty = np.asarray([0.0376517,0.0476880,0.0489686,0.0480000])

#print("purities:")
#print(purity)
#print(purity_Uncertainty)



        # RANGE & BINNING:

#pT
#pTbins = [12,15,19,26,40]
pTbins = [12,40]
#N_pT_Bins = len(pTbins)-1

#zT
zTbins = np.asarray([0.05, 0.07670637, 0.11767734, 0.18053204, 0.27695915, 0.42489062, 0.65183634,1])
zT_offset = 0
if (len(zTbins)==8):
    ZT_OFF_PLOT = 0 #Offset for FF plotting

#deta
eta_max = 1.2 #Range of Signal Correlations

#dPhi
if ("8dPhi" in description_string):
    N_dPhi_Bins = 8
    
dPhi_Bins = [i*math.pi/N_dPhi_Bins for i in range(0,N_dPhi_Bins+1)]
dPhi_Width = dPhi_Bins[1]-dPhi_Bins[0]
delta_phi_centers = [i*math.pi/N_dPhi_Bins+math.pi/N_dPhi_Bins/2 for i in range(0,N_dPhi_Bins)] #skip first dPhi bin to avoid Isolation

#UE

ue_error_bar = [] #Horiz. width of UE at first plotted dphi point

#for i,dphi in enumerate(dPhi_Bins):
#    if (dphi > 0.39):
#        ue_error_bar.append(dPhi_Bins[i])
#        break;
        
#for i,dphi in enumerate(dPhi_Bins):
#    #if (dphi > 0.78):
#    if(dphi > 1.5):
#        ue_error_bar.append(dPhi_Bins[i+1])
#        break;

ZYAM_Min_i = 0
for i,dphi in enumerate(dPhi_Bins):
    #if (dphi > 1.17):
    if (dphi > 0.39):
        ZYAM_Min_i = i
        ue_error_bar.append(dPhi_Bins[ZYAM_Min_i])
        break

ZYAM_Max_i = 0       
for i,dphi in enumerate(dPhi_Bins):
    if (dphi > 1.7):
        ZYAM_Max_i = i
        ue_error_bar.append(dPhi_Bins[ZYAM_Max_i])
        break

dphi_start_integral = 0
for i,dphi in enumerate(dPhi_Bins):
    if (dphi >= rad_start):
        dphi_start_integral = i
        break
    
#if (Uncorr_Estimate == "ZYAM"):
#    print("ZYAM Index Range = %i to %i"%(ZYAM_Min_i,ZYAM_Max_i))

N_Phi_Integrate = len(delta_phi_centers)-dphi_start_integral
#N_Phi_Integrate = len(dPhi_Bins)-dphi_start_integral
#print(N_Phi_Integrate)
Integration_Width = math.pi/(len(delta_phi_centers)+1) * N_Phi_Integrate
phi_width = math.pi/(N_dPhi_Bins)/2

                 #####CASE SWITCHING#####
if (description_string == "zT_Rebin_5"):
    zTbins = np.asarray([0.05, 0.09, 0.17, 0.30, 0.55, 1.00])

if (description_string == "zT_Rebin_6"):
    zTbins = np.asarray([0.05, 0.08, 0.14, 0.22, 0.37, 0.61, 1.00])
    
if (description_string == "zT_Rebin_8"):
    zTbins = np.asarray([0.05, 0.07, 0.11, 0.15, 0.22, 0.33, 0.47, 0.69, 1.00])

if (description_string == "zT_Rebin_9"):
    zTbins = np.asarray([0.05, 0.07, 0.10, 0.14, 0.19, 0.26, 0.37, 0.51, 0.72, 1.00])
    
if (description_string == "zT_Rebin_14"):
    zTbins = np.asarray([0.05, 0.06, 0.08, 0.10, 0.12, 0.15, 0.18, 0.22, 0.28, 0.34, 0.42, 0.53, 0.65, 0.81, 1.00])

if ("pT_Rebin_1" in description_string):
    pTbins = [12.0,40.0]
    
if ("pT_Rebin_2" in description_string):
    pTbins = [12.0,22.0,40.0]
    
if ("pT_Rebin_3" in description_string):
    pTbins = [ 12.0, 17.93, 26.78, 40.0]
    
if ("pT_Rebin_4" in description_string):
    pTbins = [12.0, 15.0, 19.0, 26.0, 40.0]
    
if ("pT_Rebin_5" in description_string):
    pTbins = [12.00, 15.27, 19.42, 24.71, 31.44, 40.00]

N_pT_Bins = len(pTbins)-1

if (description_string == "dPhi_Rebin_8"):
    N_dPhi_Bins = 8

def Get_Purity(filename):
    file = ROOT.TFile(filename)
    purities = np.zeros(N_pT_Bins)
    p_uncertainties = np.zeros(N_pT_Bins)
    
    for ipt in range(N_pT_Bins):
        purity_histo = file.Get('H_Purities_pT%1.0f_%1.0f' %(pTbins[ipt],pTbins[ipt+1]))
        #print('H_Purities_pT%1.0f_%1.0f' %(pTbins[ipt],pTbins[ipt+1]))
        purity_uncertainty_histo = file.Get('H_Purity_Uncertanty_pT%1.0f_%1.0f' %(pTbins[ipt],pTbins[ipt+1]))
        purities[ipt] = purity_histo.GetMean()
        p_uncertainties[ipt] = purity_uncertainty_histo.GetMean()
    
    return purities, p_uncertainties

if not(description_string == "05zT_working_old"):
    purity,purity_Uncertainty = Get_Purity(Files[1])
#print(purity)
    
#print("purities:")
#print(purity)
#print(purity_Uncertainty)
    
zT_widths = [(j-i)/2 for i, j in zip(zTbins[zT_offset:-1], zTbins[zT_offset+1:])]
zT_widths = np.asarray(zT_widths)
zT_centers = (zTbins[1+zT_offset:] + zTbins[zT_offset:-1]) / 2
NzT = len(zTbins)-zT_offset-1
zt_box = np.ones(NzT) * 0.03 #plotting Uncert. Boxes

if (False):
    pPb_File = '../InputData/pPb_SE_L0_Correlation_GMB_Ratio_Track.root'
    pp_File = '../InputData/pp_SE_L0_Correlation_GMB_Ratio_Track.root'

    if (CorrectedP):
        purity = np.asarray([0.208095, 0.339478, 0.483944, 0.509])
        purity_Uncertainty = np.asarray([0.0376517,0.0476880,0.0489686,0.0480000])
        #purity_Uncertainty = np.asarray([0.0273652,0.0305499,0.0358471,0.036])
    else:
        purity = [0.35]
        
        Files = [pp_File,pPb_File]

if (Shower == "NN"):
        pPb_File = '../InputData/pPb_SE_NN_Correlation_GMB_Ratio_Track.root'
        pp_File = '../InputData/pp_SE_NN_Correlation_GMB_Ratio_Track.root'
                
        if (CorrectedP):
            purity = np.asarray([0.238477, 0.341009, 0.479701, 0.532013])
            purity_Uncertainty = np.asarray([0.09,0.1,0.11,0.11])
        else:
            purity = [0.352546]

MC_File = '../InputData/18b10a_pthat_1_2_SE_NN_Correlation_GMB_Ratio.root'

if(Use_MC):
    Systems = ["pp","p-Pb","MC"]
    Files = [pp_File,pPb_File,MC_File]
    
for file,SYS in zip(Files,Systems):
    print("%s File: %s"%(SYS,file))
    
    
def Get_pp_pPb_List_Chi2(array1,array1_E,SYS_1,array2,array2_E,SYS_2):
        
    Stat_1_array = np.diag(array1_E)
    Stat_2_array = np.diag(array2_E)

    Sys_Error_Matrix_1 = np.outer(SYS_1,SYS_1) #ABSOLUTE Errors
    Sys_Error_Matrix_2 = np.outer(SYS_2,SYS_2)
    
    Cov = Stat_1_array**2 + Stat_2_array**2 + Sys_Error_Matrix_1 + Sys_Error_Matrix_2
    Cov_Inverse = np.linalg.inv(Cov)

    Delta = np.absolute(array1-array2)
    Delta_T = np.transpose(Delta)
    
    Mult1 = np.matmul(Cov_Inverse,Delta)
    
    Chi2 = np.matmul(Delta_T,Mult1)
    NDF = len(array1) - 1
    #Pval = chisqprob(Chi2,NDF)
    Pval = scipy.stats.chi2.sf(Chi2,NDF)
    return Chi2,NDF,Pval

#Each dictionary value is a numpy array of dimension 2 or 3. Want to print uncertainties.        
def print_from_Dict(Dict):
    
    for key in Dict:
        if ("CBR" in key):
            continue
        if ("pip" in key):
            continue
        if ("Uncorr" in key):
            continue
        #if ("CBR" in key):
        #    continue
        
        nparr = Dict[key]
        print("%s:"%(key))
        if (len(nparr.shape) < 2):
            for i in nparr:
                print("%1.4f,"%(i)),
            print("")
        else:
            for sublist in nparr:
                if (len(sublist.shape)) < 2:
                    for i in sublist:
                        print("%1.4f,"%(i)),
                    print("")
                else:
                    for subsublist in sublist:
                        for i in subsublist:
                            print("%1.4f,"%(i)),
                        print("")
        print("")
                    