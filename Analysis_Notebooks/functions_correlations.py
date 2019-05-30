### Corr -> Frag ###

import matplotlib.pyplot as plt
import matplotlib
import math
from default_values import *
from functions_root_nparray import ZYAM_Line

def Plot_UB_Subtraction_old(Dict):
    fsize = 20
    for SYS,ifile in zip(Systems,Files):
        if (SYS == "p-Pb"):
            print("\n \n                                       PROTON-LEAD")
        else:
            SYS=SYS
        if (SYS == "pp"):
            print("\n \n                                       PROTON-PROTON")


        for ipt in range (N_pT_Bins):
            fig = plt.figure(figsize=(34,17))
            #if (ipt > 0): continue
            #ipt = ipt+2
            for izt in range (zT_offset,NzT+zT_offset):
                ztb = izt-zT_offset

                #sig
                if (NzT ==4):
                    ax = fig.add_subplot(2,4,(2*ztb+1))
                elif (NzT ==6):
                    ax = fig.add_subplot(3,4,(2*ztb+1))
                #ax.plot(delta_phi_centers,Dict["%s_CSR"%(SYS)][ipt][ztb],'bo',ms=10)
                s_plot = ax.errorbar(delta_phi_centers,Dict["%s_CSR"%(SYS)][ipt][ztb],xerr=phi_width,
                    yerr=Dict["%s_CSR_Errors"%(SYS)][ipt][ztb],fmt='bo',ecolor='b',label='Signal Region (stat. error)')

                UE_Band = ax.fill_between(ue_error_bar,-Dict["%s_Uncorr_Error"%(SYS)][ipt][ztb][0],
                     Dict["%s_Uncorr_Error"%(SYS)][ipt][ztb][0],facecolor='purple',alpha=0.35) 

                plt.axhline(y=0,color='gray',linestyle='--',linewidth=1.3,alpha=0.8)


                plt.xlabel(r'|$\Delta \varphi$|',fontsize=fsize+4)
                plt.xticks(fontsize=(fsize))
                plt.xlim((0.39269908169872414,3.14159))
                plt.ylabel(r'$1/N_{\mathrm{trig}} \: \: \mathrm{d}N/\mathrm{d}\Delta\eta$',fontsize=fsize+2)
                empt, = ax.plot([], [], ' ')
                empt2, = ax.plot([],[],' ')
                plt.yticks(fontsize=fsize-5)

                leg = ax.legend([s_plot,UE_Band,empt,empt2],['Shower Sig. Region (stat. error)',"UB Error",r'%1.2f < $z_\mathrm{T}$ < %1.2f'
                    %(zTbins[izt],zTbins[izt+1]),r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1])],
                    loc='best',title = "Alice %s 5 TeV",fontsize=14,frameon=False,numpoints=1)

                if (SYS == 'pp'):
                    leg.set_title("ALICE Work in Progress, $\sqrt{s}=$5 TeV %s"%(SYS))
                else:
                    leg.set_title("ALICE Work in Progress, $\sqrt{s_{\mathrm{_{NN}}}}=$5 TeV %s"%(SYS))                
                plt.setp(leg.get_title(),fontsize=14)


                #bkg
                if (NzT ==4):
                    ax = fig.add_subplot(2,4,(2*ztb+2))
                elif (NzT ==6):
                    ax = fig.add_subplot(3,4,(2*ztb+2))

                plt.xlabel(r'|$\Delta \varphi$|',fontsize=fsize+4)
                plt.xticks(fontsize=(fsize))
                plt.xlim((0.39269908169872414,3.14159))
                plt.ylabel(r'$1/N_{\mathrm{trig}} \: \: \mathrm{d}N/\mathrm{d}\Delta\eta$',fontsize=fsize+2)
                plt.yticks(fontsize=fsize-5)

                ax.plot(delta_phi_centers,Dict["%s_CBR"%(SYS)][ipt][ztb],'ro',ms=10)
                b_plot = ax.errorbar(delta_phi_centers,Dict["%s_CBR"%(SYS)][ipt][ztb],xerr=phi_width,yerr=Dict["%s_CBR_Errors"%(SYS)][ipt][ztb],fmt=None,ecolor='r')
                UE_Band = ax.fill_between(ue_error_bar,-Dict["%s_Uncorr_Error"%(SYS)][ipt][ztb][0],Dict["%s_Uncorr_Error"%(SYS)][ipt][ztb][0],facecolor='purple',alpha=0.35) 
                plt.axhline(y=0,color='gray',linestyle='--',linewidth=1.3,alpha=0.8)


                leg = ax.legend([b_plot,UE_Band,empt,empt2],[' Shower Bkg Region (stat. error)',"UB Error",r'%1.2f < $z_\mathrm{T}$ < %1.2f'
                    %(zTbins[izt],zTbins[izt+1]),r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1])],
                    loc='best',fontsize=14,frameon=False,numpoints=1)
                if (SYS == 'pp'):
                    leg.set_title("ALICE Work in Progress, $\sqrt{s}=$5 TeV %s"%(SYS))
                else:
                    leg.set_title("ALICE Work in Progress, $\sqrt{s_{\mathrm{_{NN}}}}=$5 TeV %s"%(SYS))
                plt.setp(leg.get_title(),fontsize=15)
                fig.savefig('pics/%s/%s_%s_Gamma_hadron_UE_sub_zT_%i.pdf'%(Shower,Shower,SYS,izt), bbox_inches='tight')

def Plot_UB_Subtraction(Dict):
    fsize = 20
    for SYS,ifile in zip(Systems,Files):
        if (SYS == "p-Pb"):
            print("\n \n                                       PROTON-LEAD")
        else:
            SYS=SYS
        if (SYS == "pp"):
            print("\n \n                                       PROTON-PROTON")


        for ipt in range (N_pT_Bins):
            fig = plt.figure(figsize=(37,17))
            if not(Ped_Sub_First):
                fig = plt.figure(figsize=(18,12))
            #if (ipt > 0): continue
            #ipt = ipt+2
            for izt in range (0,NzT):
                ztb = izt-zT_offset

                UE_Band = "None"
                if not(Ped_Sub_First):
                    #ax = fig.add_subplot(2,3,(ztb+1))
                    if (NzT ==4):
                        ax = fig.add_subplot(2,3,(ztb+1))
                    elif (NzT ==6):
                        ax = fig.add_subplot(2,3,(ztb+1))
                    elif(NzT ==7):
                        ax = fig.add_subplot(3,3,(ztb+1))
                else:
                    #Sig
                    if (NzT ==4):
                        ax = fig.add_subplot(2,4,(2*ztb+1))
                    elif (NzT ==6):
                        ax = fig.add_subplot(3,4,(2*ztb+1))
                    elif(NzT ==7):
                        ax = fig.add_subplot(4,4,(2*ztb+1))
                        #ax.plot(delta_phi_centers,Dict["%s_CSR"%(SYS)][ipt][ztb],'bo',ms=10)
                    
                    UE_Band = ax.fill_between(ue_error_bar,-Dict["%s_Uncorr_Error"%(SYS)][ipt][ztb][0],
                                Dict["%s_Uncorr_Error"%(SYS)][ipt][ztb][0],facecolor='purple',alpha=0.35) 

                    plt.axhline(y=0,color='gray',linestyle='--',linewidth=1.3,alpha=0.8)

                s_plot = ax.errorbar(delta_phi_centers,Dict["%s_CSR"%(SYS)][ipt][ztb],xerr=phi_width,
                                     yerr=Dict["%s_CSR_Errors"%(SYS)][ipt][ztb],fmt='bo',ecolor='b',
                                     label='Signal Region (stat. error)')


                plt.xlabel(r'|$\Delta \varphi$|',fontsize=fsize+4)
                plt.xticks(fontsize=(fsize))
                plt.xlim((0.39269908169872414,3.14159))
                plt.ylabel(r'$1/N_{\mathrm{trig}} \: \: \mathrm{d}N/\mathrm{d}\Delta\eta$',fontsize=fsize+2)
                empt, = ax.plot([], [], ' ')
                empt2, = ax.plot([],[],' ')
                plt.yticks(fontsize=fsize-5)

                if not(Ped_Sub_First):
                    leg = ax.legend([s_plot,UE_Band,empt,empt2],['Shower Sig. Region (stat. error)',"UB Error",r'%1.2f < $z_\mathrm{T}$ < %1.2f'
                            %(zTbins[izt],zTbins[izt+1]),r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1])],
                            loc='best',title = "Alice %s 5 TeV",fontsize=14,frameon=False,numpoints=1)
                else:
                    leg = ax.legend([s_plot,empt,empt2],['Shower Sig. Region (stat. error)',r'%1.2f < $z_\mathrm{T}$ < %1.2f'
                            %(zTbins[izt],zTbins[izt+1]),r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1])],
                            loc='best',title = "Alice %s 5 TeV",fontsize=14,frameon=False,numpoints=1)
                    

                if (SYS == 'pp'):
                    leg.set_title("ALICE Work in Progress, $\sqrt{s}=$5 TeV %s"%(SYS))
                else:
                    leg.set_title("ALICE Work in Progress, $\sqrt{s_{\mathrm{_{NN}}}}=$5 TeV %s"%(SYS))                
                plt.setp(leg.get_title(),fontsize=14)
                if not(Ped_Sub_First):
                    plt.setp(leg.get_title(),fontsize=14)

                
                if (Ped_Sub_First):
                    #bkg
                    if (NzT ==4):
                        ax = fig.add_subplot(2,4,(2*ztb+2))
                    elif (NzT ==6):
                        ax = fig.add_subplot(3,4,(2*ztb+2))
                    elif(NzT ==7):
                        ax = fig.add_subplot(4,4,(2*ztb+1))

                    plt.xlabel(r'|$\Delta \varphi$|',fontsize=fsize+4)
                    plt.xticks(fontsize=(fsize))
                    plt.xlim((0.39269908169872414,3.14159))
                    plt.ylabel(r'$1/N_{\mathrm{trig}} \: \: \mathrm{d}N/\mathrm{d}\Delta\eta$',fontsize=fsize+2)
                    plt.yticks(fontsize=fsize-5)

                    ax.plot(delta_phi_centers,Dict["%s_CBR"%(SYS)][ipt][ztb],'ro',ms=10)
                    b_plot = ax.errorbar(delta_phi_centers,Dict["%s_CBR"%(SYS)][ipt][ztb],xerr=phi_width,yerr=Dict["%s_CBR_Errors"%(SYS)][ipt][ztb],fmt=None,ecolor='r')
                    UE_Band = ax.fill_between(ue_error_bar,-Dict["%s_Uncorr_Error"%(SYS)][ipt][ztb][0],Dict["%s_Uncorr_Error"%(SYS)][ipt][ztb][0],facecolor='purple',alpha=0.35) 
                    plt.axhline(y=0,color='gray',linestyle='--',linewidth=1.3,alpha=0.8)


                    leg = ax.legend([b_plot,UE_Band,empt,empt2],[' Shower Bkg Region (stat. error)',"UB Error",r'%1.2f < $z_\mathrm{T}$ < %1.2f'
                            %(zTbins[izt],zTbins[izt+1]),r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1])],
                            loc='best',fontsize=14,frameon=False,numpoints=1)
                    if (SYS == 'pp'):
                        leg.set_title("ALICE Work in Progress, $\sqrt{s}=$5 TeV %s"%(SYS))
                    else:
                        leg.set_title("ALICE Work in Progress, $\sqrt{s_{\mathrm{_{NN}}}}=$5 TeV %s"%(SYS))
                    plt.setp(leg.get_title(),fontsize=15)
                fig.savefig('pics/%s/%s_%s_Gamma_hadron_UE_sub_zT_%i.pdf'%(Shower,Shower,SYS,izt), bbox_inches='tight')

                
def Plot_Sub_UB_Overlay(Dict):
    #NOTE: UB Errors already propagated!
    fsize = 20
    for SYS,ifile in zip(Systems,Files):
        
        sys_color = 'red' #pp
        if (SYS == "p-Pb"):
            sys_color = 'blue'
        elif(SYS == "MC"):
             sys_color = 'green'

        for ipt in range (N_pT_Bins):
            fig = plt.figure(figsize=(50,8))
            #if (ipt > 0): continue
            #ipt = ipt+2
            for izt in range (zT_offset,NzT+zT_offset):
                ztb = izt-zT_offset

                if (NzT == 4):
                    ax = fig.add_subplot(1,4,(ztb+1))
                elif (NzT == 6):
                    ax = fig.add_subplot(1,6,(ztb+1))
                elif (NzT == 7):
                    ax = fig.add_subplot(1,7,(ztb+1))

                ax.plot(delta_phi_centers,Dict["%s_CSR"%(SYS)][ipt][ztb],'bo',color="blue",ms=10)
                s_plot = ax.errorbar(delta_phi_centers,Dict["%s_CSR"%(SYS)][ipt][ztb],xerr=phi_width,
                    yerr=Dict["%s_CSR_Errors"%(SYS)][ipt][ztb],fmt='bo',ecolor="blue",label='Signal Region (stat. error)')

                ax.plot(delta_phi_centers,Dict["%s_CBR"%(SYS)][ipt][ztb]*(1-purity[ipt]),'ro',color="red",ms=10) #Scale UE Error by purity!
                b_plot = ax.errorbar(delta_phi_centers,Dict["%s_CBR"%(SYS)][ipt][ztb]*(1-purity[ipt]),xerr=phi_width,
                    yerr=Dict["%s_CBR_Errors"%(SYS)][ipt][ztb]*(1-purity[ipt]),fmt='ro',ecolor="red",label='Background Region (stat. error)')

                UE_Band = ax.fill_between(ue_error_bar,-Dict["%s_Uncorr_Error"%(SYS)][ipt][ztb][0],Dict["%s_Uncorr_Error"%(SYS)][ipt][ztb][0],facecolor="purple",alpha=0.35) 
                plt.axhline(y=0,color='gray',linestyle='--',linewidth=1.3,alpha=0.8)


                plt.xlabel(r'|$\Delta \varphi$|',fontsize=fsize+4)
                plt.xticks(fontsize=(fsize))
                plt.xlim((0.39269908169872414,3.14159))
                plt.ylabel(r'$1/N_{\mathrm{trig}} \: \: \mathrm{d}N/\mathrm{d}\Delta\eta$',fontsize=fsize+2)
                #plt.ylim((0,1.2*max(Sig_LE_Phi_Array)))
                empt, = ax.plot([], [], ' ')
                empt2, = ax.plot([],[],' ')
                plt.yticks(fontsize=fsize-5)

                leg = ax.legend([s_plot,b_plot,UE_Band,empt,empt2],['Shower Sig Region (stat. error)','Shower Bkg Region (scaled) (stat. error)',
                    "UB Error",r'%1.2f < $z_\mathrm{T}$ < %1.2f'%(zTbins[izt],zTbins[izt+1]),r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1])],
                    loc='best',title = "Alice %s 5 TeV",fontsize=14,frameon=False,numpoints=1)
                if (SYS == 'pp'):
                    leg.set_title("ALICE Work in Progress, $\sqrt{s}=$5 TeV %s"%(SYS))
                else:
                    leg.set_title("ALICE Work in Progress, $\sqrt{s_{\mathrm{_{NN}}}}=$5 TeV %s"%(SYS))                
                plt.setp(leg.get_title(),fontsize=14)
            fig.savefig('pics/%s/%s_%s_Region_Overlays_UE_sub_pT_%i.pdf'%(Shower,Shower,SYS,ipt), bbox_inches='tight')

def Correlated_Subtraction(Dict):   
    for SYS in Systems:
        for ipt in range (N_pT_Bins):
            #if (ipt > 0): continue
            
            Dict["%s_CSR"%(SYS)][ipt] = (Dict["%s_CSR"%(SYS)][ipt] - (1-purity[ipt])*Dict["%s_CBR"%(SYS)][ipt])/purity[ipt]
            Dict["%s_CSR_Errors"%(SYS)][ipt] = Dict["%s_CSR_Errors"%(SYS)][ipt]/purity[ipt]
            Dict["%s_Uncorr_Error"%(SYS)][ipt] = Dict["%s_Uncorr_Error"%(SYS)][ipt]/purity[ipt]
            
            
            #if not(Ped_Sub_First):
            #    Ped_Sub_After_Cs(Dict)
            
def Correlated_Subtraction_Weights(Dict):   
    for SYS in Systems:
        for ipt in range (N_pT_Bins):
            #if (ipt > 0): continue
            
            Dict["%s_CSR"%(SYS)][ipt] = (Dict["%s_CSR"%(SYS)][ipt] - Dict["%s_CBR"%(SYS)][ipt])
            Dict["%s_CSR_Errors"%(SYS)][ipt] = np.sqrt((Dict["%s_CSR_Errors"%(SYS)][ipt])**2 + (Dict["%s_CBR_Errors"%(SYS)][ipt])**2)
            Dict["%s_Uncorr_Error"%(SYS)][ipt] = Dict["%s_Uncorr_Error"%(SYS)][ipt]
            
            
            #if not(Ped_Sub_First):
            #    Ped_Sub_After_Cs(Dict)



def Ped_Sub_After_Cs(Dict):
#This function calculates ZYAM after correlated subtraction (Cs)
    for SYS,ifile in zip(Systems,Files):
        for ipt in range (N_pT_Bins):
            #if (ipt > 0): continue

            for izt in range (NzT):
                ZYAM_Cs,ZYAM_Cs_Error = ZYAM_Line(Dict["%s_CSR"%(SYS)][ipt][izt],
                                                  Dict["%s_CSR_Errors"%(SYS)][ipt][izt])

                Dict["%s_CSR"%(SYS)][ipt][izt] = Dict["%s_CSR"%(SYS)][ipt][izt] - ZYAM_Cs
                Dict["%s_Uncorr_Error"%(SYS)][ipt][izt] = ZYAM_Cs_Error
                                                  

def Get_pp_pPb_List_Chi2(array1,array1_E,array2,array2_E):
    
    hist1 = ROOT.TH1F("hist1","histo1",N_dPhi_Bins, 0, 1)
    for i in range(len(array1)-1):
        hist1.SetBinContent(i+1,array1[i])
        hist1.SetBinError(i+1,array1_E[i])
        
    hist2 = ROOT.TH1F("hist2","histo2",N_dPhi_Bins, 0, 1)
    for i in range(len(array2)-1):
        hist2.SetBinContent(i+1,array2[i])
        hist2.SetBinError(i+1,array2_E[i])

    chi2 = ROOT.Double(0.0)
    ndf = ROOT.Long(0.0)
    igood = ROOT.Long(0.0)
    pval = hist1.Chi2TestX(hist2,chi2,ndf,igood,"WW")

    return pval,chi2,ndf,igood

def getfitvals(Dict):
    for ipt in range(N_pT_Bins):
        for ztb in range(NzT):
            pval,chi2,ndf,igood = Get_pp_pPb_List_Chi2(Dict["p-Pb_CSR"][ipt][ztb],Dict["p-Pb_CSR_Errors"][ipt][ztb],Dict["pp_CSR"][ipt][ztb],Dict["pp_CSR_Errors"][ipt][ztb])
            print("zT %i: pval = %f, chi2 = %f, ndf = %f"%(ztb+1,pval,chi2,ndf))
                    
                    
def Plot_pp_pPb_Cs(Dict):
    
    for ipt in range (N_pT_Bins):
        #if (ipt > 0): continue
        #ipt = ipt+2
        #plt.figure(figsize=(10,7))
        fig = plt.figure(figsize=(24,12))
        if (NzT==7):
            fig = plt.figure(figsize=(24,18))
        for izt in range (zT_offset,NzT+zT_offset):
            ztb = izt-zT_offset

            if (NzT ==4):
                ax = fig.add_subplot(2,2,ztb+1)
            elif (NzT ==6):
                ax = fig.add_subplot(2,3,ztb+1)
            elif (NzT ==7):
                ax = fig.add_subplot(3,3,ztb+1)

            pPb = plt.errorbar(delta_phi_centers,Dict["p-Pb_CSR"][ipt][ztb],xerr=phi_width,yerr=Dict["p-Pb_CSR_Errors"][ipt][ztb],fmt='bo',capsize=4,markersize=11)
            pp = plt.errorbar(delta_phi_centers,Dict["pp_CSR"][ipt][ztb],xerr=phi_width,yerr=Dict["pp_CSR_Errors"][ipt][ztb],fmt='ro',capsize=4,markersize=11)

            if(Use_MC):
                MC = plt.errorbar(delta_phi_centers,["MC_CSR"][ipt][ztb],xerr=phi_width,yerr=["MC_CSR_Errors"][ipt][ztb],fmt='go',capsize=4,markersize=11)



            plt.xlabel(r'|$\Delta \varphi$|',fontsize=28)
            plt.xticks(fontsize=18)
            plt.xlim((0.39269908169872414,3.14159))
            plt.ylabel(r'$1/N_{\gamma} \: \: \mathrm{d}N/\mathrm{d}\Delta \eta$',fontsize=28)
            plt.yticks(fontsize=18)
            plt.axhline(y=0,color='gray',linestyle='--',linewidth=1.3,alpha=0.8)        

            pp_UE = ax.fill_between(ue_error_bar,-Dict["pp_Uncorr_Error"][ipt][ztb][0],Dict['pp_Uncorr_Error'][ipt][ztb][0],facecolor='red',alpha=0.35) #Other for pp
            pPb_UE = ax.fill_between(ue_error_bar,-Dict["p-Pb_Uncorr_Error"][ipt][ztb][0],Dict['p-Pb_Uncorr_Error'][ipt][ztb][0],facecolor='blue',alpha=0.35)#One for p-Pb
            #MC_UE = ax.fill_between(ue_error_bar,-Dict["MC_Uncorr_Error"][ipt][ztb][0],Dict["MC_Uncorr_Error"][ipt][ztb][0],facecolor='green',alpha=0.35)#One for p-Pb

            empt, = ax.plot([], [],' ')
            empt2, = ax.plot([],[],' ')
            
            empt3, = ax.plot([],[],' ')
            pval,chi2,ndf,igood = Get_pp_pPb_List_Chi2(Dict["p-Pb_CSR"][ipt][ztb],Dict["p-Pb_CSR_Errors"][ipt][ztb],Dict["pp_CSR"][ipt][ztb],Dict["pp_CSR_Errors"][ipt][ztb])
            fit_string = "pval = %1.2f, chi2 = %1.1f, ndf = %i"%(pval,chi2,ndf)

            if(Use_MC):
                leg = plt.legend([pp,pPb,MC,pp_UE,pPb_UE,MC_UE,empt,empt2,empt3],['pp $\sqrt{s}= 5$ TeV (stat. error)',
                    'p-Pb $\sqrt{s_{\mathrm{_{NN}}}}=5$ TeV (stat. error)','pythia GJ $\sqrt{s}=5$ TeV (stat. error)', 
                    'pp UE Error', 'p-Pb UE Error','pythia UE Error', r'%1.2f < $z_\mathrm{T}$ < %1.2f'%(zTbins[izt],zTbins[izt+1]),
                    r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1]),fit_string],loc = "upper left",fontsize=16,frameon=False,numpoints=1)
            else:    
                leg = plt.legend([pp,pPb,pp_UE,pPb_UE,empt,empt2,empt3],['pp $\sqrt{s}= 5$ TeV (stat. error)',
                    'p-Pb $\sqrt{s_{\mathrm{_{NN}}}}=5$ TeV (stat. error)', 'pp UB Error', 'p-Pb UB Error',
                    r'%1.2f < $z_\mathrm{T}$ < %1.2f'%(zTbins[izt],zTbins[izt+1]),
                    r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1]),fit_string],loc = "upper left",fontsize=16,frameon=False,numpoints=1)

            
            
            leg.set_title("ALICE Work in Progress")
            leg._legend_box.align = "left"
            plt.setp(leg.get_title(),fontsize=22)
            fig.tight_layout()

        plt.show()
        #fig.savefig('pics/Gamma_hadron_corr_zT_%i.pdf'%(ztb))
        fig.savefig('pics/%s/%s_Gamma_hadron_corr_%i.pdf'%(Shower,Shower,ipt))
        
        
        #Above can be adapted for pp & PbPb comparisons
        
def Integrate_Away_Side(Phi_array,Phi_Errors,LE_Error):
    
    Use_Uncorr_Error = True
    FF_zt = np.zeros((N_pT_Bins, NzT))
    FF_zt_Errors = np.zeros((N_pT_Bins, NzT))
    if Use_Uncorr_Error:
        LE_Error = LE_Error*Integration_Width
    for ipt in range(N_pT_Bins):
        for izt in range(0, NzT):
            ztb = izt
            zT_width = zTbins[izt+1]-zTbins[izt]

            temp_phi = Phi_array[ipt][ztb][(len(Phi_Errors[ipt][ztb])-N_Phi_Integrate):]
            FF_zt[ipt][ztb] = temp_phi.sum()/zT_width
            temp_error = (Phi_Errors[ipt][ztb][(len(Phi_Errors[ipt][ztb])-N_Phi_Integrate):])**2
            if (Use_Uncorr_Error):
                FF_zt_Errors[ipt][ztb] = (math.sqrt(temp_error.sum() + (LE_Error[ipt][izt][0])**2))/zT_width
            else:
                FF_zt_Errors[ipt][ztb] = math.sqrt(temp_error.sum())/zT_width
            
    return FF_zt, FF_zt_Errors

def Get_Fragmentation(Dict):
    
    Keys = []
    
    for SYS in Systems:
        Keys.append("%s_FF"%(SYS))
        Keys.append("%s_FF_Errors"%(SYS))
        Keys.append("%s_purity_FF_Errors"%(SYS))
    
    FF_Vals = np.zeros((len(Keys),N_pT_Bins,NzT))
    FF_Dict = dict(zip(Keys,FF_Vals))

    for SYS in Systems:
        
        temp_FF, temp_FF_Errors = Integrate_Away_Side(Dict["%s_CSR"%(SYS)],Dict["%s_CSR_Errors"%(SYS)],Dict["%s_Uncorr_Error"%(SYS)])
        temp_purity_Errors = []
        for ipt in range(N_pT_Bins):
            temp_purity_Errors.append(temp_FF[ipt]*(purity_Uncertainty[ipt]/purity[ipt])*Integration_Width)  # abs. FF purity uncertainty
        
        FF_Dict["%s_FF"%(SYS)], FF_Dict["%s_FF_Errors"%(SYS)] = temp_FF, temp_FF_Errors
        FF_Dict["%s_purity_FF_Errors"%(SYS)] = np.asarray(temp_purity_Errors)
    
    for SYS in Systems:
        if (Use_Weights):
            np.save("npy_files/%s_%s_Fragmentation_Functions.npy"%(Shower,SYS),FF_Dict["%s_FF"%(SYS)])
            np.save("npy_files/%s_%s_Fragmentation_Function_Errors.npy"%(Shower,SYS),FF_Dict["%s_FF_Errors"%(SYS)])
            np.save("npy_files/%s_%s_FF_purity_Uncertainty.npy"%(Shower,SYS),FF_Dict["%s_purity_FF_Errors"%(SYS)])
        else:
            np.save("npy_files/%s_%s_Fragmentation_Functions_Unweight.npy"%(Shower,SYS),FF_Dict["%s_FF"%(SYS)])
            np.save("npy_files/%s_%s_Fragmentation_Function_Errors_Unweight.npy"%(Shower,SYS),FF_Dict["%s_FF_Errors"%(SYS)])
            np.save("npy_files/%s_%s_FF_purity_Uncertainty_Unweight.npy"%(Shower,SYS),FF_Dict["%s_purity_FF_Errors"%(SYS)])
    
    return FF_Dict

def Plot_FF(FF_Dict):
    fig = plt.figure(figsize=(17,13))
    for ipt in range(N_pT_Bins):

        zt_box = np.ones(NzT) * 0.03
        #pPb_bar = plt.bar(zT_centers, pPb_sys[ipt]+pPb_sys[ipt], bottom=(pPb_FF[ipt])-pPb_sys[ipt], width=zt_box, align='center',edgecolor="blue",color='white',)
        #pp_bar = plt.bar(zT_centers, pp_sys[ipt]+pp_sys[ipt], bottom=pp_FF[ipt]-pp_sys[ipt], width=zt_box, align='center',edgecolor="red",color='white',)

        zT_max = 0
        for izt in range(zT_offset, NzT+zT_offset):
            if (zTbins[izt]*pTbins[ipt] > 15.0):
                zT_max = zTbins[izt]
                break
        
        if (N_pT_Bins < 5):
            ax = fig.add_subplot(2,2,ipt+1)
        elif (N_pT_Bins >=5):
            ax = fig.add_subplot(3,2,ipt+1)
            
        
        #fig = plt.figure(figsize=(10,7))
        
        pPb_plot = plt.errorbar(zT_centers, FF_Dict["p-Pb_FF"][ipt],xerr=zT_widths,yerr=FF_Dict["p-Pb_FF_Errors"][ipt],linewidth=1, fmt='bo',capsize=1,label='p-Pb')
        pp_plot = plt.errorbar(zT_centers, FF_Dict["pp_FF"][ipt],xerr=zT_widths,yerr=FF_Dict["pp_FF_Errors"][ipt],linewidth=1,fmt='ro',capsize=1,label='pp')

        if(Use_MC):
            plt.errorbar(zT_centers, FF_Dict["MC_FF"][ipt],xerr=zT_widths,yerr=FF_Dict["MC_FF_Errors"][ipt],linewidth=1, fmt='go',capsize=1,label='MC')
            MC_plot = plt.fill_between(zT_centers, FF_Dict["MC_FF"][ipt]-FF_Dict["MC_FF_Errors"][ipt], FF_Dict["MC_FF"][ipt]+FF_Dict["MC_FF_Errors"][ipt],color='green',alpha=0.6)

        empt, = plt.plot([], [],' ')

        plt.yscale('log')                                                                                                                                                                                                                                                                                                                                          
        plt.ylabel(r"$\frac{1}{N_{\mathrm{\gamma}}}\frac{\mathrm{d}N}{\mathrm{d}z_{\mathrm{T}} \mathrm{d}\Delta\eta}$",fontsize=20)
        plt.xlabel("${z_\mathrm{T}} = p_\mathrm{T}^\mathrm{h}/p_\mathrm{T}^\mathrm{\gamma}$",fontsize=20)
        #plt.xlim(xmin = 0.1,xmax=0.7)
        plt.ylim(ymin = 0.001,ymax=10)

        if(Use_MC):
            leg = plt.legend([pp_plot,pPb_plot,MC_plot,pp_bar,pPb_bar,empt],["pp $\sqrt{s} = 5$ TeV","p-Pb $\sqrt{s_{\mathrm{_{NN}}}}=5$ TeV","Pythia GJ $\sqrt{s} = 5$ TeV",
                "pp Systematic","p-Pb Systematic","Normalization $\pm 25\%$"],frameon=False,numpoints=1,title=' ',loc="upper right",prop={'size':14})
        else:
            #leg = plt.legend([pp_plot,pPb_plot,pp_bar,pPb_bar,empt],["pp $\sqrt{s} = 5$ TeV","p-Pb $\sqrt{s_{\mathrm{_{NN}}}}=5$ TeV",
            #    "pp Systematic","p-Pb Systematic","Normalization $\pm 25\%$"],frameon=False,numpoints=1,title=' ',loc="upper right",prop={'size':14})
            leg = plt.legend([pp_plot,pPb_plot,empt],["pp $\sqrt{s} = 5$ TeV","p-Pb $\sqrt{s_{\mathrm{_{NN}}}}=5$ TeV",
                r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1])],frameon=False,numpoints=1,title=' ',loc="upper right",prop={'size':14})

        leg.set_title("ALICE Work in Progress")
        plt.setp(leg.get_title(),fontsize=20)

        Title = plt.title(r'Integrated $\mathrm{\gamma}$-Hadron Correlation: $2\pi/3 < \Delta\varphi < \pi, |\Delta\eta| < %1.1f$ '%(eta_max),fontsize=15)
        plt.gcf()
        plt.savefig("pics/%s_Systems_FFunction_%i.pdf"%(description_string,ipt), bbox_inches='tight')
        #plt.show()

    plt.gcf()
    plt.savefig("pics/Systems_FFunction_All.pdf", bbox_inches='tight')
    plt.show()
