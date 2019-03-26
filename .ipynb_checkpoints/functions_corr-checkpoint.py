### ROOT -> numpy ###

import ROOT
from default_values import *
import matplotlib.pyplot as plt

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
    
    z_temp = Phi_Array[2:3]
    Z_Value = z_temp.mean()
    
    z_temp_error = Phi_Error_Array[2:3]**2
    Z_Error = z_temp_error.sum()
    Z_Error = math.sqrt(Z_Error)
    Z_Error = Z_Error/len(z_temp_error)
    
    return Z_Value,Z_Error

def GetLEProj(filename, ipt, izt, Signal_DNN=True,DoAverage=True):
    file = ROOT.TFile(filename)
    if (Signal_DNN == "Inclusive"):
        histo2D = file.Get('Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                        %(pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]))
        
        Eta_Axis = histo2D.Get
        Yaxis()
        
        LE_Projection = histo2D.ProjectionX('Inclusive_PhiProjection__pT_%1.0f_%1.0f__zt_%1.0f_%1.0f' 
                                      %(pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],
                                        100*zTbins[izt+1]),Eta_Axis.FindBin(-1.4),Eta_Axis.FindBin(-0.8))
                                        #10*zTbins[izt+1]),5,11)
        LE_Projection_pos = histo2D.ProjectionX('PosEta_inclusive_PhiProjection__pT_%1.0f_%1.0f__zt_%1.0f_%1.0f' 
                                  %(pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],
                                    100*zTbins[izt+1]),Eta_Axis.FindBin(0.8),Eta_Axis.FindBin(1.4))
  
        ntriggers = Get_NTriggers(filename, ipt, Signal_DNN)  
        if not(ntriggers == None):
            LE_Projection.Scale(1.0/ntriggers) #per trigger yield
            LE_Projection_pos.Scale(1.0/ntriggers)
    else:
        DNN_Rgn = int(Signal_DNN) + 2*(1-int(Signal_DNN)) #convert bool to DNN_1 (Sgn) or DNN_2 (Bkgd)

        if (DoAverage):
            histo2D_1 = file.Get('DNN1_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                                 %(pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]))

            if not(Use_Weights):
                histo2D_2 = file.Get('Unweighted_DNN2_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                                    %(pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]))
            else:    
                histo2D_2 = file.Get('DNN2_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                                    %(pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]))

        else:
            if (not(Use_Weights) and not(Signal_DNN)):
                histo2D_1 = file.Get('Unweighted_DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                                %(DNN_Rgn,pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]))

                histo2D_2 = file.Get('Unweighted_DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                                %(DNN_Rgn,pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]))
            else:
                histo2D_1 = file.Get('DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                                %(DNN_Rgn,pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]))

                histo2D_2 = file.Get('DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f' 
                                %(DNN_Rgn,pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]))

        
        #Project:
        Eta_Axis = histo2D_1.GetYaxis()
  
        ntriggers_DNN2 = Get_NTriggers(filename, ipt, False) 
        
        if (DoAverage):
            ntriggers_DNN1 = Get_NTriggers(filename, ipt, True) 
        else:
            ntriggers_DNN1 = Get_NTriggers(filename, ipt, Signal_DNN)   
        
        #Should be 2 lines, but root can only integrate continuous distributions...
        LE_Projection = histo2D_1.ProjectionX('NegEta_DNN1_PhiProjection__pT_%1.0f_%1.0f__zt_%1.0f_%1.0f' 
                                  %(pTbins[ipt],pTbins[ipt+1],10*zTbins[izt],
                                    10*zTbins[izt+1]),Eta_Axis.FindBin(-1.4),Eta_Axis.FindBin(-0.8))
        
        LE_Projection_pos = histo2D_1.ProjectionX('PosEta_DNN1_PhiProjection__pT_%1.0f_%1.0f__zt_%1.0f_%1.0f' 
                                  %(pTbins[ipt],pTbins[ipt+1],10*zTbins[izt],
                                    10*zTbins[izt+1]),Eta_Axis.FindBin(0.8),Eta_Axis.FindBin(1.4))    
        
        LE_Projection_DNN2 = histo2D_2.ProjectionX('NegEta2_DNN2_PhiProjection__pT_%1.0f_%1.0f__zt_%1.0f_%1.0f' 
                                  %(pTbins[ipt],pTbins[ipt+1],10*zTbins[izt],
                                    10*zTbins[izt+1]),Eta_Axis.FindBin(-1.4),Eta_Axis.FindBin(-0.8))

        LE_Projection_DNN2_pos = histo2D_2.ProjectionX('PosEta2_DNN2_PhiProjection__pT_%1.0f_%1.0f__zt_%1.0f_%1.0f' 
                                  %(pTbins[ipt],pTbins[ipt+1],10*zTbins[izt],
                                    10*zTbins[izt+1]),Eta_Axis.FindBin(0.8),Eta_Axis.FindBin(1.4))
        
        LE_Projection.Scale(1.0/ntriggers_DNN1)
        LE_Projection_pos.Scale(1.0/ntriggers_DNN1)
        
        LE_Projection_DNN2.Scale(1.0/ntriggers_DNN2)
        LE_Projection_DNN2_pos.Scale(1.0/ntriggers_DNN2)
        
        if (DoAverage):
            LE_Projection.Add(LE_Projection_DNN2,1.0)
            LE_Projection.Scale(0.5)
        
            LE_Projection_pos.Add(LE_Projection_DNN2_pos,1.0)
            LE_Projection_pos.Scale(0.5)
        
    
    #Add,scale :
    
    #LE_Projection.SetDirectory(0)
    LE_Projection.Add(LE_Projection_pos,1)
    LE_Projection.Rebin(2)
    LE_Projection.Scale(1.0/1.2) #scale by eta region
    
    file.Close()
    
    LE_Phi_Array = np.zeros(len(delta_phi_centers))
    LE_Error_Array = np.zeros(len(delta_phi_centers))
    for bin in range(2,9):
                LE_Phi_Array[bin-2] = LE_Projection.GetBinContent(bin)
                LE_Error_Array[bin-2] = LE_Projection.GetBinError(bin)
    
    return LE_Phi_Array, LE_Error_Array

def GetLE_Val( LE_Phi_Array, LE_Error_Array):
    LE_temp = LE_Phi_Array[:3]
    LE_value = LE_temp.mean()
    
    LE_temp_error = LE_Error_Array[:3]**2
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
    PhiProjection.Rebin(2)
    PhiProjection.Scale(1.0/(2*eta_max))
    
    #per trigger yield
    ntriggers = Get_NTriggers(filename,ipt, Signal_DNN)
    if not(ntriggers == None):
        PhiProjection.Scale(1.0/ntriggers)
    
    file.Close()
    
    Phi_Array = np.zeros(len(delta_phi_centers))
    Phi_Error_Array = np.zeros(len(delta_phi_centers))
    for bin in range(2,9):
        Phi_Array[bin-2] = PhiProjection.GetBinContent(bin)
        Phi_Error_Array[bin-2] = PhiProjection.GetBinError(bin)
    
    return Phi_Array, Phi_Error_Array


def Plot_UB():
    for sys,ifile in zip(Systems,Files):
        print("%s: $z_T$ interval   & LE Signal Region & LE Background Region & ZYAM Signal Region & ZYAM Background Region"%(sys))
        for ipt in range (N_pT_Bins):
            fig = plt.figure(figsize=(34,17))
            #if (ipt > 0): continue
            #ipt = ipt+2
            for izt in range (zT_offset,NzT+zT_offset):

                ztb = izt-zT_offset

                Sig_LE_Phi_Array, Sig_LE_Error_Array = GetLEProj(ifile, ipt, izt, True,False)#2 Bools: Signal_Region, LE Region Averegaing
                Bkg_LE_Phi_Array, Bkg_LE_Error_Array = GetLEProj(ifile, ipt, izt, False,False)

                S_LE, S_LE_Error = GetLE_Val(Sig_LE_Phi_Array, Sig_LE_Error_Array)
                Bkg_LE, Bkg_LE_Error = GetLE_Val(Bkg_LE_Phi_Array, Bkg_LE_Error_Array)


                Sig_Phi_Array, Sig_Phi_Error_Array = GetPhiProj(ifile,sys,ipt,izt,True)
                Bkg_Phi_Array, Bkg_Phi_Error_Array = GetPhiProj(ifile,sys,ipt,izt,False)


                Sig_Z_Value,Sig_Z_Error = ZYAM_Line(Sig_Phi_Array, Sig_Phi_Error_Array)
                Bkg_Z_Value,Bkg_Z_Error = ZYAM_Line(Bkg_Phi_Array, Bkg_Phi_Error_Array)

                print("%1.2f - %1.2f & %1.3f $\pm$ %1.3f & %1.3f $\pm$ %1.3f & %1.3f $\pm$ %1.3f & %1.3f $\pm$ %1.4f \\\\"
                      %(zTbins[izt],zTbins[izt+1],S_LE,S_LE_Error,Bkg_LE,Bkg_LE_Error,Sig_Z_Value,Sig_Z_Error,Bkg_Z_Value,Bkg_Z_Error))


                                            #--------------plot--------------------#

                if (NzT ==4):
                    ax = fig.add_subplot(2,4,(2*ztb+1))
                elif (NzT ==6):
                    ax = fig.add_subplot(3,4,(2*ztb+1))
                    

                fsize = 20

                #sig
                ax.plot(delta_phi_centers,Sig_Phi_Array,'bo',ms=10)
                s_plot = ax.errorbar(delta_phi_centers,Sig_Phi_Array,xerr=phi_width,yerr=Sig_Phi_Error_Array,fmt=None,ecolor='b',label='Signal Region (stat. error)')

                ax.plot(delta_phi_centers,Sig_LE_Phi_Array,'s',color="Grey",alpha=0.6,ms=10)
                s_le_plot = ax.errorbar(delta_phi_centers,Sig_LE_Phi_Array,xerr=phi_width,yerr=Sig_LE_Error_Array,fmt=None,ecolor='Grey',alpha=0.8,label="0.8 <|$\Delta\eta$| < 1.4")

                plt.xlabel(r'|$\Delta \varphi$|',fontsize=fsize+4)
                plt.xticks(fontsize=(fsize))
                plt.xlim((0.39269908169872414,3.14159))
                plt.ylabel(r'$1/N_{\mathrm{trig}} \: \: \mathrm{d}N/\mathrm{d}\Delta\eta$',fontsize=fsize+2)
                plt.ylim((-0.01,1.2*max(Sig_LE_Phi_Array)))
                empt, = ax.plot([], []," ")
                empt2, = ax.plot([],[]," ")
                plt.yticks(fontsize=fsize-5)

                fill_x = [0,3.14159]
                s_z_line = ax.fill_between(fill_x, Sig_Z_Value-Sig_Z_Error,Sig_Z_Value+Sig_Z_Error,interpolate=False,edgecolor='cyan',linewidth=0.0, alpha=0.3,facecolor='cyan')
                s_le_line = ax.fill_between(fill_x, S_LE-S_LE_Error,S_LE+S_LE_Error,interpolate=False,edgecolor='grey',linewidth=0.0, alpha=0.5,facecolor='grey')

                leg = ax.legend([s_plot,s_le_plot,s_le_line,s_z_line,empt,empt2],['Shower Sig. Region (stat. error)','0.8 <|$\Delta\eta$| < 1.4','UB Estimate',
                    'ZYAM',r'%1.2f < $z_\mathrm{T}$ < %1.2f'%(zTbins[izt],zTbins[izt+1]),r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1])],
                    loc='best',title = "Alice %s 5 TeV",fontsize=14,frameon=False,numpoints=1)
                if (sys == 'pp'):
                    leg.set_title("ALICE, $\sqrt{s}=$5 TeV %s"%(sys))
                else:
                    leg.set_title("ALICE, $\sqrt{s_{\mathrm{_{NN}}}}=$5 TeV %s"%(sys))                
                plt.setp(leg.get_title(),fontsize=14)


                #bkg
                if (NzT ==4):
                    ax = fig.add_subplot(2,4,(2*ztb+2))
                elif (NzT ==6):
                    ax = fig.add_subplot(3,4,(2*ztb+2))

                #ax = fig.add_subplot(1,2,1)
                plt.xlabel(r'|$\Delta \varphi$|',fontsize=fsize+4)
                plt.xticks(fontsize=(fsize))
                plt.xlim((0.39269908169872414,3.14159))
                plt.ylabel(r'$1/N_{\mathrm{trig}} \: \: \mathrm{d}N/\mathrm{d}\Delta\eta$',fontsize=fsize+2)
                plt.ylim((-0.01,1.2*max(Bkg_LE_Phi_Array)))
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
                    leg.set_title("ALICE, $\sqrt{s}=$5 TeV %s"%(sys))
                else:
                    leg.set_title("ALICE, $\sqrt{s_{\mathrm{_{NN}}}}=$5 TeV %s"%(sys))
                plt.setp(leg.get_title(),fontsize=15)
            fig.savefig('pics/%s/%s_%s_Gamma_hadron_UE_pT_%i__zT_%i.pdf'%(Shower,Shower,sys,ipt,izt), bbox_inches='tight')        
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

    Corr_Arrays = np.zeros((len(Keys),N_pT_Bins,NzT,N_dPhi_Bins-1))

    Dict = dict(zip(Keys, Corr_Arrays))


    for SYS,ifile in zip(Systems,Files):    
        for ipt in range (N_pT_Bins):

            for izt in range (zT_offset,NzT+zT_offset):
                ztb = izt-zT_offset

                SR_UB = 0
                SR_UB_Error = 0
                BR_UB = 0
                BR_UB_Error = 0

                UB_Error = 0 #Final Value

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

                    Sig_Phi_Array -= Avg_UB
                    Bkg_Phi_Array -= Avg_UB

                else:
                    Sig_Phi_Array -= SR_UB
                    Bkg_Phi_Array -= BR_UB
                    UB_Error = np.sqrt(SR_UB_Error**2 + BR_UB_Error**2)

                    Dict["%s_CSR"%(SYS)][ipt][izt] = Sig_Phi_Array
                    Dict["%s_CSR_Errors"%(SYS)][ipt][izt] = Sig_Phi_Error_Array
                    Dict["%s_CBR"%(SYS)][ipt][izt] = Bkg_Phi_Array
                    Dict["%s_CBR_Errors"%(SYS)][ipt][izt] = Bkg_Phi_Error_Array
                    Dict["%s_Uncorr_Error"%(SYS)][ipt][izt] = UB_Error #constant in for zT bins

    return Dict