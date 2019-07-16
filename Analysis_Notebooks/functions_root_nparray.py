### ROOT -> numpy ###

import ROOT
from default_values import *
import matplotlib.pyplot as plt
import pandas as pd

def Get_NTriggers(filename,ipt, Signal_DNN=True): 
    file = ROOT.TFile(filename)
    if (Signal_DNN == "Inclusive"):
        ntrig_histo = file.Get('N_Triggers_pT%1.0f_%1.0f' %(pTbins[ipt],pTbins[ipt+1]))
    else:
        DNN_Rgn = int(Signal_DNN) + 2*(1-int(Signal_DNN))
        ntrig_histo = file.Get('N_DNN%i_Triggers_pT%1.0f_%1.0f' %(DNN_Rgn,pTbins[ipt],pTbins[ipt+1]))
    NTriggers = 1
    if not(ntrig_histo == None):
        if (Signal_DNN):
            NTriggers = ntrig_histo.GetEntries()
        else:
            if (Use_Weights):
                NTriggers = ntrig_histo.GetBinContent(1)
            else:
                NTriggers = ntrig_histo.GetEntries()
    file.Close()
    return NTriggers
    
    
#UNCORRELATED BACKGROUND:

def ZYAM_Line(Phi_Array, Phi_Error_Array):
    
    z_temp = Phi_Array[ZYAM_Min_i:ZYAM_Max_i]
    
    z_temp_error = Phi_Error_Array[ZYAM_Min_i:ZYAM_Max_i]**2

    Z_Value = z_temp.mean()

    Z_Error = z_temp_error.sum()
    Z_Error = math.sqrt(Z_Error)
    
    return Z_Value,Z_Error

def GetLEProj(filename, ipt, izt, Signal_DNN=True,DoAverage=True):
    file = ROOT.TFile(filename)
    
    eta_min = 0.9
    eta_max = 1.1

    DNN_Rgn = int(Signal_DNN) + 2*(1-int(Signal_DNN)) #convert bool to DNN_1 (Sgn) or DNN_2 (Bkgd)


    histo2D_1 = file.Get('DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                                %(DNN_Rgn,pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]))

    histo2D_2 = file.Get('DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                                %(DNN_Rgn,pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]))

        
    #Project:
    Eta_Axis = histo2D_1.GetYaxis()
  
    ntriggers_DNN1 = Get_NTriggers(filename, ipt, Signal_DNN)   
        
        #Should be 2 lines, but root can only integrate continuous distributions...
    LE_Projection = histo2D_1.ProjectionX('NegEta_DNN1_PhiProjection__pT_%1.0f_%1.0f__zt_%1.0f_%1.0f' 
                                  %(pTbins[ipt],pTbins[ipt+1],10*zTbins[izt],
                                    10*zTbins[izt+1]),Eta_Axis.FindBin(-eta_max),Eta_Axis.FindBin(-eta_min))
        
    LE_Projection_pos = histo2D_1.ProjectionX('PosEta_DNN1_PhiProjection__pT_%1.0f_%1.0f__zt_%1.0f_%1.0f' 
                                  %(pTbins[ipt],pTbins[ipt+1],10*zTbins[izt],
                                    10*zTbins[izt+1]),Eta_Axis.FindBin(eta_min),Eta_Axis.FindBin(eta_max))    
        
    #print("%1.3f, %1.3f"%(Eta_Axis.FindBin(-1.1),Eta_Axis.FindBin(-0.7)))
    #print("%1.3f, %1.3f"%(Eta_Axis.FindBin(0.7),Eta_Axis.FindBin(1.1)))
    LE_Projection.Scale(1.0/ntriggers_DNN1)
    LE_Projection_pos.Scale(1.0/ntriggers_DNN1)
        
        
    
    #Add,scale :
    
    LE_Projection.SetDirectory(0)
    LE_Projection.Add(LE_Projection_pos,1)

    LE_Projection.Scale(1.0/((eta_max+0.1-(eta_min-0.1))*2)) #scale by eta region
    #print((eta_max+0.1-(eta_min-0.1))*2)
    
    file.Close()
    
    LE_Phi_Array = np.zeros(len(delta_phi_centers))
    LE_Error_Array = np.zeros(len(delta_phi_centers))
    for bin in range(0,N_dPhi_Bins):
                LE_Phi_Array[bin] = LE_Projection.GetBinContent(bin+1)
                LE_Error_Array[bin] = LE_Projection.GetBinError(bin+1)
    
    return LE_Phi_Array, LE_Error_Array


def GetLE_Val( LE_Phi_Array, LE_Error_Array):
    LE_temp = LE_Phi_Array[ZYAM_Min_i:ZYAM_Max_i]
    LE_value = LE_temp.mean()
    
    LE_temp_error = LE_Error_Array[ZYAM_Min_i:ZYAM_Max_i]**2
    LE_Error = LE_temp_error.sum()
    LE_Error = math.sqrt(LE_Error)
    LE_Error = LE_Error/len(LE_temp_error)
    
    return LE_value, LE_Error


#CORRELATIONS, PIH PROJECTIONS,& PLOTTING:


def GetPhiProj(filename,prfx,ipt, izt, Signal_DNN=True):
    
    file = ROOT.TFile(filename)
    
    if (Signal_DNN == "Inclusive"):
        histo2D = file.Get('Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                        %(pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]))
        
        Eta_Axis = histo2D.GetYaxis()
        PhiProjection = histo2D.ProjectionX('Inclusive_PhiProjection__pT_%1.0f_%1.0f__zt_%1.0f_%1.0f' 
                                      %(pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],
                                        100*zTbins[izt+1]),Eta_Axis.FindBin(-eta_max),Eta_Axis.FindBin(eta_max))
    
    else:
        DNN_Rgn = int(Signal_DNN) + 2*(1-int(Signal_DNN)) #convert bool to DNN_1 (Sgn) or DNN_2 (Bkgd)
        
        if (not(Signal_DNN) and not(Use_Weights) ):
        #if not(Signal_DNN or Use_Weights):
            histo2D = file.Get('Unweighted_DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                            %(DNN_Rgn,pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]))
            
        else:
            histo2D = file.Get('DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                            %(DNN_Rgn,pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]))
        
        Eta_Axis = histo2D.GetYaxis()
        PhiProjection = histo2D.ProjectionX('DNN%i_PhiProjection__pT_%1.0f_%1.0f__zt_%1.0f_%1.0f' 
                                      %(DNN_Rgn,pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],
                                        100*zTbins[izt+1]),Eta_Axis.FindBin(-eta_max),Eta_Axis.FindBin(eta_max))
                                            
    PhiProjection.SetDirectory(0)
#    if (N_dPhi_Bins == 8):
#        PhiProjection.Rebin(2)
    #PhiProjection.Scale(1.0/(2*eta_max))
    PhiProjection.Scale(1.0/(2*1.2))
    
    #per trigger yield
    ntriggers = Get_NTriggers(filename,ipt, Signal_DNN)
    if not(ntriggers == None or ntriggers == 0):
        PhiProjection.Scale(1.0/ntriggers)
    
    file.Close()
    
    Phi_Array = np.zeros(len(delta_phi_centers))
    Phi_Error_Array = np.zeros(len(delta_phi_centers))

    for bin in range(N_dPhi_Bins):
        Phi_Array[bin] = PhiProjection.GetBinContent(bin+1)
        Phi_Error_Array[bin] = PhiProjection.GetBinError(bin+1)
    
    return Phi_Array, Phi_Error_Array


def Plot_UB():
    for sys,ifile in zip(Systems,Files):
        print("%s: $z_T$ interval   & LE Signal Region & LE Background Region & ZYAM Signal Region & ZYAM Background Region"%(sys))
        for ipt in range (N_pT_Bins):
            fig = plt.figure(figsize=(34,28))
            #if (ipt > 0): continue
            #ipt = ipt+2
            for izt in range (0,NzT):

                ztb = izt-zT_offset

                Sig_LE_Phi_Array, Sig_LE_Error_Array = GetLEProj(ifile, ipt, izt, True,False)#2 Bools: Signal_Region, LE Region Averegaing
                Bkg_LE_Phi_Array, Bkg_LE_Error_Array = GetLEProj(ifile, ipt, izt, False,False)

                S_LE, S_LE_Error = GetLE_Val(Sig_LE_Phi_Array, Sig_LE_Error_Array)
                
                print("Large Eta Error per unit phi")
                print(S_LE/(phi_width))
                print(S_LE_Error/phi_width)
                
                Bkg_LE, Bkg_LE_Error = GetLE_Val(Bkg_LE_Phi_Array, Bkg_LE_Error_Array)


                Sig_Phi_Array, Sig_Phi_Error_Array = GetPhiProj(ifile,sys,ipt,izt,True)
                Bkg_Phi_Array, Bkg_Phi_Error_Array = GetPhiProj(ifile,sys,ipt,izt,False)


                Sig_Z_Value,Sig_Z_Error = ZYAM_Line(Sig_Phi_Array, Sig_Phi_Error_Array)
                Sig_Z_Error = Sig_Z_Error/(ZYAM_Max_i-ZYAM_Min_i)
                Bkg_Z_Value,Bkg_Z_Error = ZYAM_Line(Bkg_Phi_Array, Bkg_Phi_Error_Array)
                
                
                print("ZYAM per unit phi")
                print(Sig_Z_Value/phi_width)
                print(Sig_Z_Error/phi_width)

                print("%1.2f - %1.2f & %1.3f $\pm$ %1.3f & %1.3f $\pm$ %1.3f \\\\"
                      %(zTbins[izt],zTbins[izt+1],S_LE,S_LE_Error,Sig_Z_Value,Sig_Z_Error))


                                       #--------------plot--------------------#

                ax = fig.add_subplot(3,3,izt+1)
                if (NzT ==4):
                    ax = fig.add_subplot(2,2,izt+1)
                elif (NzT ==6):
                    ax = fig.add_subplot(2,3,izt+1)
                elif (NzT >=7 and NzT <=9):
                    ax = fig.add_subplot(3,3,izt+1)
                elif (NzT >9 and NzT <=12):
                    ax = fig.add_subplot(4,3,izt+1)
                elif (NzT >12):
                    ax = fig.add_subplot(5,3,izt+1)
                    

                fsize = 20

                #sig
                ax.plot(delta_phi_centers,Sig_Phi_Array,'bo',ms=10)
                s_plot = ax.errorbar(delta_phi_centers,Sig_Phi_Array,xerr=phi_width,yerr=Sig_Phi_Error_Array,fmt=None,ecolor='b',label='Signal Region (stat. error)')

                ax.plot(delta_phi_centers,Sig_LE_Phi_Array,'s',color="Grey",alpha=0.6,ms=10)
                s_le_plot = ax.errorbar(delta_phi_centers,Sig_LE_Phi_Array,xerr=phi_width,yerr=Sig_LE_Error_Array,fmt=None,ecolor='Grey',alpha=0.8,label="0.8 <|$\Delta\eta$| < %1.1f"%(eta_max))

                plt.xlabel(r'|$\Delta \varphi$|',fontsize=fsize+4)
                plt.xticks(fontsize=(fsize))
                plt.xlim((0.39269908169872414,3.14159))
                plt.ylabel(r'$1/N_{\mathrm{trig}} \: \: \mathrm{d}N/\mathrm{d}\Delta\eta$',fontsize=fsize+2)
                #plt.ylim((-0.001,1.2*max(Sig_LE_Phi_Array)))
                empt, = ax.plot([], []," ")
                empt2, = ax.plot([],[]," ")
                plt.yticks(fontsize=fsize-5)

                fill_x = [0,3.14159]
                s_z_line = ax.fill_between(ue_error_bar, Sig_Z_Value-Sig_Z_Error,Sig_Z_Value+Sig_Z_Error,interpolate=False,edgecolor='blue',linewidth=0.0, alpha=0.6,facecolor='blue')
                s_le_line = ax.fill_between(ue_error_bar, S_LE-S_LE_Error,S_LE+S_LE_Error,interpolate=False,edgecolor='green',linewidth=0.0, alpha=0.5,facecolor='green')

                leg = ax.legend([s_plot,s_le_plot,s_le_line,s_z_line,empt,empt2],['Shower Sig. Region (stat. error)','0.8 <|$\Delta\eta$| < %1.1f'%(eta_max),'Large $\Delta\eta$ Estimate',
                    'ZYAM',r'%1.2f < $z_\mathrm{T}$ < %1.2f'%(zTbins[izt],zTbins[izt+1]),r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1])],
                    loc='best',title = "Alice %s 5 TeV",fontsize=14,frameon=False,numpoints=1)
                if (sys == 'pp'):
                    leg.set_title("ALICE Work in Progress, $\sqrt{s}=$5 TeV %s"%(sys))
                else:
                    leg.set_title("ALICE Work in Progress, $\sqrt{s_{\mathrm{_{NN}}}}=$5 TeV %s"%(sys))                
                plt.setp(leg.get_title(),fontsize=14)

                continue
                
                
                #bkg
                if (NzT ==4):
                    ax = fig.add_subplot(2,4,(2*ztb+2))
                elif (NzT ==6):
                    ax = fig.add_subplot(3,4,(2*ztb+2))
                elif (NzT ==7):
                    ax = fig.add_subplot(4,4,(2*ztb+2))
                    
                #ax = fig.add_subplot(1,2,1)
                plt.xlabel(r'|$\Delta \varphi$|',fontsize=fsize+4)
                plt.xticks(fontsize=(fsize))
                plt.xlim((0.39269908169872414,3.14159))
                plt.ylabel(r'$1/N_{\mathrm{trig}} \: \: \mathrm{d}N/\mathrm{d}\Delta\eta$',fontsize=fsize+2)
                #plt.ylim((-0.01,1.2*max(Bkg_LE_Phi_Array)))
                plt.yticks(fontsize=fsize-5)

                fill_x = [0,3.14149]
                b_z_line = plt.fill_between(fill_x, Bkg_Z_Value-Bkg_Z_Error,Bkg_Z_Value+Bkg_Z_Error,interpolate=False,edgecolor='cyan',linewidth=0.0, alpha=0.3,facecolor='cyan')
                b_le_line = plt.fill_between(fill_x, Bkg_LE-Bkg_LE_Error,Bkg_LE+Bkg_LE_Error,interpolate=False,edgecolor='grey',linewidth=0.0, alpha=0.5,facecolor='grey')

                ax.plot(delta_phi_centers,Bkg_Phi_Array,'ro',ms=10)
                b_plot = ax.errorbar(delta_phi_centers,Bkg_Phi_Array,xerr=phi_width,yerr=Bkg_Phi_Error_Array,fmt=None,ecolor='r')

                plt.plot(delta_phi_centers,Bkg_LE_Phi_Array,'s',color="Grey",alpha=0.6,ms=10)

                b_le_plot =plt.errorbar(delta_phi_centers,Bkg_LE_Phi_Array,xerr=phi_width,yerr=Bkg_LE_Error_Array,fmt=None,ecolor='Grey',alpha=0.8)
                leg = ax.legend([b_plot,b_le_plot,b_le_line,b_z_line,empt,empt2],['Shower Bkg Region (stat. error)','0.8 <|$\Delta\eta$| < 1.4','UB Estimate',
                    'ZYAM',r'%1.2f < $z_\mathrm{T}$ < %1.2f'%(zTbins[izt],zTbins[izt+1]),r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1])],
                    loc='best',fontsize=14,frameon=False,numpoints=1)
                if (sys == 'pp'):
                    leg.set_title("ALICE Work in Progress, $\sqrt{s}=$5 TeV %s"%(sys))
                else:
                    leg.set_title("ALICE Work in Progress, $\sqrt{s_{\mathrm{_{NN}}}}=$5 TeV %s"%(sys))
                plt.setp(leg.get_title(),fontsize=15)
            fig.savefig('pics/%s/%s/UE_Plot_%s_pT_%i__zT_%i.pdf'%(Shower,description_string,sys,ipt,izt), bbox_inches='tight')        
            print("")
        print("")
        #return

                
def ROOT_to_nparray():
    Keys = []
    Corr_Arrays = []

    for SYS,ifile in zip(Systems,Files):

        Keys.append("%s_CSR"%(SYS))
        Keys.append("%s_CSR_Errors"%(SYS))
        Keys.append("%s_CBR"%(SYS))
        Keys.append("%s_CBR_Errors"%(SYS))
        Keys.append("%s_Uncorr_Error"%(SYS))
        Keys.append("%s_Uncorr_Estimate"%(SYS))

    Corr_Arrays = np.zeros((len(Keys),N_pT_Bins,NzT,N_dPhi_Bins))

    Dict = dict(zip(Keys, Corr_Arrays))


    for SYS,ifile in zip(Systems,Files):    
        for ipt in range (N_pT_Bins):

            for izt in range (0,NzT):

                SR_UB = 0
                SR_UB_Error = 0
                BR_UB = 0
                BR_UB_Error = 0

                UB_Error = 0 

                Sig_Phi_Array, Sig_Phi_Error_Array = GetPhiProj(ifile,SYS,ipt,izt,True)
                Bkg_Phi_Array, Bkg_Phi_Error_Array = GetPhiProj(ifile,SYS,ipt,izt,False)


                if (Uncorr_Estimate == "LargeEta"):

                    Sig_LE_Phi_Array, Sig_LE_Error_Array = GetLEProj(ifile, ipt, izt, True,True) # Using AVERAGE LE
                    Bkg_LE_Phi_Array, Bkg_LE_Error_Array = GetLEProj(ifile, ipt, izt, False,True) # Using AVERAGE LE

                    SR_UB, SR_UB_Error = GetLE_Val(Sig_LE_Phi_Array, Sig_LE_Error_Array)
                    BR_UB, BR_UB_Error = GetLE_Val(Bkg_LE_Phi_Array, Bkg_LE_Error_Array)

                elif (Uncorr_Estimate == "ZYAM"):
                    SR_UB,SR_UB_Error = ZYAM_Line(Sig_Phi_Array, Sig_Phi_Error_Array)
                    BR_UB,BR_UB_Error = ZYAM_Line(Bkg_Phi_Array, Bkg_Phi_Error_Array)                                                                


                if (Average_UE):
                    Avg_UB = (SR_UB + BR_UB)/2
                    UB_Error = np.sqrt(SR_UB_Error**2 + BR_UB_Error**2)/2
                    
                    if (Ped_Sub_First):
                        Sig_Phi_Array -= Avg_UB
                        Bkg_Phi_Array -= Avg_UB

                else:

                    if (Ped_Sub_First):
                        Sig_Phi_Array -= SR_UB
                        Bkg_Phi_Array -= BR_UB
                        
                    UB_Error = np.sqrt(SR_UB_Error**2 + BR_UB_Error**2)

                    Dict["%s_CSR"%(SYS)][ipt][izt] = Sig_Phi_Array
                    Dict["%s_CSR_Errors"%(SYS)][ipt][izt] = Sig_Phi_Error_Array
                    Dict["%s_CBR"%(SYS)][ipt][izt] = Bkg_Phi_Array
                    Dict["%s_CBR_Errors"%(SYS)][ipt][izt] = Bkg_Phi_Error_Array
                    Dict["%s_Uncorr_Error"%(SYS)][ipt][izt] = UB_Error #constant in for zT bins
                    Dict["%s_Uncorr_Estimate"%(SYS)][ipt][izt] = SR_UB

    return Dict    
        
