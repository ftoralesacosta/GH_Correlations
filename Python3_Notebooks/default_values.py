import numpy as np
import math
import ROOT
import scipy.stats
from hepdata_lib import Submission
submission = Submission()
                 #####CASE SWITCHING#####

#Shower = "NN"
Shower = "LO"
Use_Weights = True #BKR pT Weight
Use_P_Weights = True
CorrectedP = True  # FALSE FOR HARDPROBES
Use_MC = True
pT_Rebin = False
N_dPhi_Bins = 8
Ped_Sub_First = False #Important after weight implementation
Average_UE = False
Show_Fits = True
Uncorr_Estimate = "ZYAM"

#rad_start = 1.6
#Phi_String = "\pi/2"

#rad_start = 1.9
#Phi_String = "5\pi/8"

#rad_start = 2.3
#Phi_String = "3\pi/4"

#rad_start = 2.15
#Phi_String = "11\pi/16"

rad_start = 2.7
Phi_String = "7\pi/8"

#rad_start = 2.9
#Phi_String = "15\pi/16"


ZYAM_Start = 0.39
#ZYAM_Start = 0.41

ZYAM_End = 1.5
#ZYAM_End = 1.3

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
#description_string = "zT_Rebin_10"
#description_string = "zT_Rebin_12"
#description_string = "zT_Rebin_12pT150"
#description_string = "zT_Rebin_12_06zT"
#description_string = "zT_Rebin_8_06zT" 
#description_string = "zT_Rebin_9_06zT"
#description_string = "zT_Rebin_9_004zT06zT" 
#description_string = "zT_Rebin_9_003zT06zT"
#description_string = "zT_Rebin_9_004zT06zTcheck"

#description_string = "zT_Rebin_10_003zT06zT"
#description_string = "zT_Rebin_14"


#description_string = "zT_Rebin_7_006zT06zTpT2"

#description_string = "zT_Rebin_6_006zT06zT"
#description_string = "zT_Rebin_7_006zT06zT"

#description_string = "zT_Rebin_8_006zT06zT" 
#description_string = "zT_Rebin_8_006zT06zTOldBinNewPurity"
#description_string = "zT_Rebin_8_006zT06zTOldBinNewNewPurity"

#----------^^^ OLD ^^^ --------


#default_string = "zT_Rebin_8_006zT06zTOldBinNewNewPurity8dPhi"
#default_string = "zT_Rebin_8_006zT06zT"
default_string = "zT_Rebin_8_006zT06zT13fnew"
description_string = default_string
#description_string = "IRCProp"
#description_string="pT_Rebin_4"
#description_string = "zT_Rebin_8_006zT06zTpPb"
#description_string = "zT_Rebin_8_006zT06zTPbp"
#description_string = "zT_Rebin_8_006zT06zT16dPhi"
#description_string = "zT_Rebin_8_006zT06zTWeightTest"

#description_string = "zT_Rebin_8_006zT06zTwSDD"
#description_string = "zT_Rebin_8_006zT06zTnewCuts"
#description_string = "zT_Rebin_8_006zT06zTnewCutsOldIso"
#description_string = "zT_Rebin_8_006zT06zTnoEventCuts"
#description_string = "zT_Rebin_8_006zT06zTnoEeventnoTOF"
#description_string = "zT_Rebin_8_006zT06zTREPRODUCE"
#description_string = "zT_Rebin_8_006zT06zTnoNLM"
#description_string = "zT_Rebin_8_006zT06zTwNLMwAllElse"
#description_string = "zT_Rebin_8_006zT06zTREPRODUCE2"
#description_string = "zT_Rebin_8_006zT06zTREPRODUCE3"
#description_string = "zT_Rebin_8_006zT06zTjustTOF"
#description_string = "zT_Rebin_8_006zT06zTTOFLambdaMin"
#description_string = "zT_Rebin_8_006zT06zTEventSel"
#description_string = "zT_Rebin_8_006zT06zTPileUp"
#description_string = "zT_Rebin_8_006zT06zTZvTX"
#description_string = "zT_Rebin_8_006zT06zTREPRODUCE4"
#description_string = "zT_Rebin_8_006zT06zTZvTX"
#description_string = "zT_Rebin_8_006zT06zTTPCIso"
#description_string = "zT_Rebin_8_006zT06zTITSSub16dPhi"
#description_string = "zT_Rebin_8_006zT06zTpileCut"
#description_string = "zT_Rebin_8_006zT06zTTPCTRACKS"
#description_string = "zT_Rebin_8_006zT06zTTPCTRACKS16dPhi"
#description_string = "zT_Rebin_8_006zT06zTITSSub" #Current Default
#description_string = "zT_Rebin_8_006zT06zTITSSubpPb"
#description_string = "zT_Rebin_8_006zT06zTITSSubPbp"
#description_string = "zT_Rebin_6_006zT06zTITSSub"
#description_string = "zT_Rebin_6_006zT06zTpileCut"
#description_string = "zT_Rebin_6_006zT06zT8GeV"
#description_string = "zT_Rebin_8_006zT06zT6GeV"
#description_string = "zT_Rebin_8_006zT06zT8GeV"
#description_string = "zT_Rebin_8_006zT06zT13fnew"
#description_string = "zT_Rebin_8_006zT06zT13fnewONLY"
#description_string = "zT_Rebin_8_006zT06zT13fnewChi2"
#description_string = "zT_Rebin_8_006zT06zT13fnewChi36noWeights"
#description_string = "zT_Rebin_8_006zT06zT13fnewChi2noWeights"
#escription_string = "zT_Rebin_8_006zT06zT13fnewChi3GJ"
#description_string = "zT_Rebin_8_006zT06zT13fnewTPCChi3GJ"
#description_string = "zT_Rebin_8_006zT06zT13fnewTPCChi36GJ"

 
#description_string = "zT_Rebin_7_006zT06zT"
#description_string = "zT_Rebin_9_006zT06zT"
#description_string = "pT_Rebin_2_006zT06zTOldBinNewNewPurity"
#description_string = "zT_Rebin_8_006zT06zTOldBinNewNewPurity_Small_Zyam"
#description_string= "zT_Rebin_8_006zT06zTpT2"


#description_string = "zT_Rebin_8_006zT06zTminpT15"
#description_string = "zT_Rebin_8_006zT06zT_Small_Zyam_Avg"

pPb_File = '../InputData/%s/pPb_SE_L0_Correlation_GMB_Ratio.root'%(description_string)
pp_File = '../InputData/%s/pp_SE_L0_Correlation_GMB_Ratio.root'%(description_string)
print(pp_File)
Systems = ["pp","p-Pb"]
Files = [pp_File,pPb_File]


#Apply Eta correction:
Eta_Correction = 1.05 #5% Correction
Eta_Correction_Uncertainty = 0.025
Apply_Eta_Correction = True
        # RANGE & BINNING:

#pT
#pTbins = [12,15,19,26,40]
pTbins = [12,40]
#N_pT_Bins = len(pTbins)-1

#zT
zTbins = np.asarray([0.060, 0.080, 0.107, 0.142, 0.190, 0.253, 0.337, 0.450, 0.600])
zT_offset = 0
ZT_OFF_PLOT = 0 #Offset for FF plotting

#deta
eta_max = 1.19 #Range of Signal Correlations

#dPhi
if ("8dPhi" in description_string):
    N_dPhi_Bins = 8
    
if ("16dPhi" in description_string):
    N_dPhi_Bins = 16
    
dPhi_Bins = [i*math.pi/N_dPhi_Bins for i in range(0,N_dPhi_Bins+1)]
dPhi_Width = dPhi_Bins[1]-dPhi_Bins[0]
delta_phi_centers = [i*math.pi/N_dPhi_Bins+math.pi/N_dPhi_Bins/2 for i in range(0,N_dPhi_Bins)] #skip first dPhi bin to avoid Isolation
delta_phi_edges = []
for p in range(0,N_dPhi_Bins):
    delta_phi_edges.append(tuple([(p*dPhi_Width),((p+1)*dPhi_Width)]))
#UE

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
    
if ("zT_Rebin_7_006zT06zT" in description_string):
    zTbins = np.asarray([0.060, 0.083, 0.116, 0.161, 0.224, 0.311, 0.432, 0.600])

if ("zT_Rebin_6_006zT06zT" in description_string):
    zTbins = np.asarray([0.060, 0.088, 0.129, 0.190, 0.278, 0.409, 0.600])

if ("zT_Rebin_8_006zT06zT" in description_string or "zT_Rebin_8" in description_string):
    zTbins = np.asarray([0.060, 0.080, 0.107, 0.142, 0.190, 0.253, 0.337, 0.450, 0.600])
    
if ("zT_Rebin_9_006zT06zT" in description_string):
    zTbins = np.asarray([0.060, 0.077, 0.100, 0.129, 0.167, 0.216, 0.278, 0.360, 0.465, 0.600])
    
if (description_string == "zT_Rebin_8"):
    zTbins = np.asarray([0.05, 0.07, 0.11, 0.15, 0.22, 0.33, 0.47, 0.69, 1.00])
    
if (description_string == "zT_Rebin_8_06zT"):
    zTbins = np.asarray([0.050, 0.068, 0.093, 0.127, 0.173, 0.236, 0.322, 0.440, 0.600])

if (description_string == "zT_Rebin_9"):
    zTbins = np.asarray([0.05, 0.07, 0.10, 0.14, 0.19, 0.26, 0.37, 0.51, 0.72, 1.00])
    
if(description_string == "zT_Rebin_9_06zT"):
    zTbins = np.asarray([0.050, 0.066, 0.087, 0.114, 0.151, 0.199, 0.262, 0.345, 0.455, 0.600])
    
if("zT_Rebin_9_004zT06zT" in description_string):
    zTbins = np.asarray([0.040, 0.054, 0.073, 0.099, 0.133, 0.180, 0.243, 0.329, 0.444, 0.600])
        
if (description_string == "zT_Rebin_10"):
    zTbins = np.asarray([0.05, 0.0675, 0.0910, 0.1228, 0.1657, 0.2236, 0.3017, 0.4071, 0.5493, 0.7411, 1.00])
    
if ("zT_Rebin_12" in description_string):
    zTbins = np.asarray([0.05, 0.0642, 0.0824, 0.1057, 0.1357,0.1742, 0.2236 , 0.287, 0.368, 0.473 ,0.607, 0.78, 1.0])
    
if ("zT_Rebin_12_06zT" in description_string):
    zTbins = np.asarray([0.050, 0.062, 0.076, 0.093, 0.114, 0.141, 0.173, 0.213, 0.262, 0.322, 0.397, 0.488, 0.600])
                        
if (description_string == "zT_Rebin_10_003zT06zT"):
    zTbins = np.asarray([0.033, 0.044, 0.059, 0.079, 0.105, 0.141, 0.188, 0.251, 0.336, 0.449, 0.600])
    
if (description_string == "zT_Rebin_9_003zT06zT"):
    zTbins = np.asarray([0.033, 0.046, 0.063, 0.087, 0.120, 0.165, 0.228, 0.315, 0.435, 0.600])
        
if (description_string == "zT_Rebin_14"):
    zTbins = np.asarray([0.05, 0.06, 0.08, 0.10, 0.12, 0.15, 0.18, 0.22, 0.28, 0.34, 0.42, 0.53, 0.65, 0.81, 1.00])
    
if ("pT_Rebin_1" in description_string):
    pTbins = [12.0,40.0]
    
if ("pT_Rebin_2" in description_string or "pT2" in description_string):
    pTbins = [12.0,22.0,40.0]
    
if ("pT_Rebin_3" in description_string):
    pTbins = [ 12.0, 17.93, 26.78, 40.0]
    
if ("pT_Rebin_4" in description_string):
    pTbins = [12.0, 15.0, 19.0, 26.0, 40.0]
    
if ("pT_Rebin_5" in description_string):
    pTbins = [12.00, 15.27, 19.42, 24.71, 31.44, 40.00]

if ("minpT13" in description_string):
    pTbins = [13.0,40.0]
    
N_pT_Bins = len(pTbins)-1

if (description_string == "dPhi_Rebin_8"):
    N_dPhi_Bins = 8

if ("Small_Zyam_Avg" in description_string):
    ZYAM_Start = 0.41

    ZYAM_End = 1.3
    
ue_error_bar = [] #Horiz. width of UE at first plotted dphi point

ZYAM_Min_i = 0
for i,dphi in enumerate(dPhi_Bins):
    #if (dphi > 1.17):
    if (dphi > ZYAM_Start):
        ZYAM_Min_i = i
        ue_error_bar.append(dPhi_Bins[ZYAM_Min_i])
        break

ZYAM_Max_i = 0       
for i,dphi in enumerate(dPhi_Bins):
    if (dphi > ZYAM_End):
        ZYAM_Max_i = i
        ue_error_bar.append(dPhi_Bins[ZYAM_Max_i])
        break
        
        
Max_Hadron_pT = 10.0
Min_Hadron_pT = 0.5    

def Get_Purity(filename):

    p_uncertainties = []
    purities = []
    for iter_file in Files:
        file = ROOT.TFile(iter_file)
        purities_pTs = np.zeros(N_pT_Bins)
        p_uncertainties_pTs = np.zeros(N_pT_Bins)
    
        for ipt in range(N_pT_Bins):
            purity_histo = file.Get('H_Purities_pT%1.0f_%1.0f' %(pTbins[ipt],pTbins[ipt+1]))
            purities_pTs[ipt] = purity_histo.GetMean()
            
            purity_uncertainty_histo = file.Get('H_Purity_Uncertanty_pT%1.0f_%1.0f' %(pTbins[ipt],pTbins[ipt+1]))
            p_uncertainties_pTs[ipt] = purity_uncertainty_histo.GetMean()
        
            purities.append(purities_pTs)
            p_uncertainties.append(p_uncertainties_pTs) #np array of pt bin purity for each system

    return purities, p_uncertainties

p,p_uncert = Get_Purity(Files)
purity = dict(zip(Systems,p))
purity_Uncertainty = dict(zip(Systems,p_uncert))
Rel_pUncert = dict((k, float(purity_Uncertainty[k]) / purity[k]) for k in purity_Uncertainty)

print(purity.keys())
print(purity.values())
temp_p_uncert = np.array(purity_Uncertainty.values())
temp_p = np.array(purity.values())
#print(temp_p_uncert/temp_p)                

zT_widths = [(j-i)/2. for i, j in zip(zTbins[zT_offset:-1], zTbins[zT_offset+1:])]
zT_widths = np.asarray(zT_widths)
zT_centers = (zTbins[1+zT_offset:] + zTbins[zT_offset:-1]) / 2
NzT = len(zTbins)-zT_offset-1
zt_box = np.ones(NzT) * 0.03 #plotting Uncert. Boxes

zT_edges = []
for i in range(0,NzT):
    zT_edges.append(tuple([zTbins[i],zTbins[i+1]]))
zT_edges

import itertools
#marker = itertools.cycle((',', '+', '.', 'o', '*'))
marker = itertools.cycle(('.', 'o', '^', '+', '*')) 



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
                    
            

pythia = [[-0.0016102345883117469,-0.0006754869947770241,-0.0014505300420320222,8.697118151178296e-07,-0.0005723781659106187,0.0002319217485895423,-3.3887573769554244e-05,-0.00048770845917885976,0.0004877084591788615,0.0028804160910520567,0.00310653763232498,0.006895373958291436,0.010373820086637762,0.014915717863947752,0.020746763135237066,0.022318536875564062],[-0.0016619206977794435,0.0003999247589527525,0.0006099119560241256,0.000532312518306462,0.0003748486655376932,0.00028788191910829666,0.0006461576827871216,-0.0008634114551029004,0.0008634114551028995,0.0028475420767304083,0.0030981573898599245,0.0044974008038530765,0.008126045997189439,0.013639483322021572,0.02225766732404616,0.024367476213950565],[-0.001502694024509732,-0.00026276796689440104,0.0018814775230559825,0.0014229966655862798,-0.0005880993919460476,0.0009149027724786474,0.0003389107361687427,-0.0005994486172264512,0.0005994486172264512,0.0014391692568924384,0.0025320070071560823,0.0032653648866322335,0.008159579986473525,0.011989068107837463,0.01959014638635437,0.023081737202538546],[-0.0013410521225892657,-0.000727873753142429,-0.00010326932363808163,0.0008247925132202367,0.0012770344829802624,-3.8905464038635745e-05,-0.00024808297725446675,-5.6869413837063766e-05,5.68694138370642e-05,0.0008928336338385572,0.0013334021205060644,0.003758603135659619,0.005370728028139511,0.01038052886186025,0.014856978437016458,0.01899010139708396],[-0.0010117870381265801,-0.00012827690790699012,0.0003766736699152679,-0.0003704915242387723,0.00010433849611320677,-0.00014007123231873746,-0.001409254874172857,0.00012113141061082964,-0.00012113141061082921,-0.00015582803524515987,0.0011944719069993445,0.0013206773321110619,0.0037066626744888583,0.005610583133844517,0.010781752821106556,0.015794753090953975],[-0.0004721328056127951,-0.0008234400715871529,0.0010684391475931539,0.00019817906078323146,0.00016567700159373744,9.256798402917151e-06,1.8968817659410298e-05,-0.00022329948280979122,0.00022329948280979143,0.0007938144699302443,0.001343167748793859,0.0016333209957574016,0.004009616480960281,0.005736579569882582,0.008313509350745487,0.012934911233870531],[-0.0006920225696556129,-0.00042059509142724626,0.000499508078802313,0.0006138334718964252,0.0005123292550067748,-0.00024871255303665553,-0.00010697585257998202,-0.00015928259203781033,0.00015928259203781044,7.987565715372952e-05,0.00032109944186991006,0.001551944653311061,0.0021963107517722377,0.0030488312317554682,0.004583386171947362,0.008129098530744436],[-0.0005128646941286784,-0.000401646034352389,0.0004960424586265302,3.6097039871697866e-05,-0.0003868860182594831,-4.885715059975361e-05,-0.00047122858889796443,9.521529706052858e-05,-9.521529706052858e-05,0.000647960545344772,-2.042188578302111e-05,0.0003670025254467686,0.0007155020994642126,0.0018079721580679022,0.003605308774773443,0.004766081781293104]]


pythia_error = [[0.0012303622478521583,0.0012330561932697342,0.0011990254077779048,0.0012673968733941845,0.0012421248619426388,0.0012650007368390957,0.0012520774710275719,0.0012150068549236848,0.001240820642732318,0.0012824694649687248,0.001283607944660397,0.0013373801624987634,0.0013904918017342367,0.00145262714407455,0.0014984645357907418,0.0014993704589513165],[0.001037681392612962,0.0010938511258161972,0.001096958306183931,0.0010892270524774002,0.0010939336860799381,0.001095473419610736,0.0011038094720328762,0.0010274432519906255,0.0010845484791806527,0.0011363859357670855,0.0011319199921518447,0.001136943066681792,0.001190173689843855,0.0012872063033582543,0.0013897347821918062,0.0013730219208489527],[0.0008632886434643671,0.0009021730706101397,0.0009746687156638884,0.0009574391249139033,0.0008721165554090953,0.0009274558967668641,0.0008996815396510124,0.0008586148491423139,0.0009193978755658163,0.0009214910828329808,0.000961072734237045,0.000935463499180012,0.0010480763989039901,0.0011148287532904235,0.0012358768976291192,0.0012316243455911715
],[0.0008151727568052208,0.0008467606917622548,0.0008657878191986545,0.0009122114825061199,0.0009097483107055916,0.0008755878640608482,0.0008532948367654673,0.0008685804830049017,0.0008560484938469718,0.000885322117950491,0.0008820138636204805,0.0009458452487014607,0.0009748016708206666,0.0010800987269872437,0.00112059842824301,0.001140180871605562
],[0.0006689638085290101,0.0007168427603673966,0.0007364017125343849,0.0006905205582789034,0.0007236319524342524,0.0006794481214093358,0.0006307105872354927,0.0006937125153302192,0.0006903789873219305,0.0006720713199994423,0.000732178083627205,0.0007218260387124358,0.0007746984856783489,0.0008105704261508647,0.0008923380393544973,0.000989934457537705],[0.0005282666600671808,0.0005011093808726591,0.0006262572007234288,0.0005774013933479912,0.0005847539595135945,0.0005606691109252617,0.0005514258003046244,0.0005391255358387238,0.0005816625304021708,0.0005908115172107477,0.0006035638444179746,0.000607288797777647,0.0006894264571366028,0.0007136004327116789,0.0007851656226851943,0.000853179873321096],[0.00042810982076691546,0.000446048179723369,0.0005092021059053372,0.000522415781913585,0.0004937541358818619,0.000439733483942661,0.00046441649950578514,0.00045744250021813337,0.00048730020682122597,0.0004773327482437276,0.00046501168076932013,0.0005547959053835643,0.0005575512205329529,0.0005598657789041747,0.0005892270850505426,0.0007257270960205526],[0.0002565985410419714,0.00026377123892711084,0.00037741502175339355,0.0003137870658491026,0.00025251354221806336,0.00028042515894366736,0.0002485464856006564,0.0003259753360766937,0.0003170969493561224,0.0003678704135602845,0.00031186867908816887,0.0003186984460260079,0.0003382560222876157,0.000415601292238216,0.0004973582727263059,0.0005151827299237345]]

pythia_FF = [5.483244297963561, 4.397405669002265, 3.1046588362122227, 1.7956433549410293, 1.0742304218451844, 0.6441506988451846, 0.28647855283731016, 0.14211714000541883]

if (N_dPhi_Bins == 8):

    for z in range(len(pythia)): #loop over zT
        ipyth = pythia[z]
        ipyth_error = pythia_error[z]
        e = 0
        for i in range(N_dPhi_Bins):
            pythia[z][i] = ipyth[e] + ipyth[e+1]
            pythia_error[z][i] = (ipyth_error[e] + ipyth_error[e+1])
            e = e+2
        
        pythia[z] = pythia[z][:8]
        pythia_error[z] = pythia_error[z][:8]

pythia  = np.asarray(pythia)/dPhi_Width
pythia_error = np.asarray(pythia_error)/dPhi_Width

pythia_FF_Errors = np.zeros(N_dPhi_Bins)
for izt in range(0,NzT):
    pythia_FF_Errors[izt] = pythia_error[izt][-1:]/zT_widths[izt]


#Models
CNM_zt = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2], dtype='float64')
CNM_IpPb = np.array([1.064, 1.036, 1.056, 1.100, 1.141, 1.180, 1.157, 1.181, 1.202, 1.164, 1.206, 1.237], dtype='float64')

#With parton energy loss
QGP_zt = np.array([0.15, 0.35, 0.55, 0.75, 0.95, 1.15], dtype='float64')
QGP_IpPb_LowAll = np.array([0.974025, 0.909975, 0.95335, 0.943225, 0.954, 0.9837], dtype='float64')
QGP_IpPb_HighAll = np.array([1.0018, 0.963775, 0.9727, 0.9789, 0.995675, 1.00545], dtype='float64')
