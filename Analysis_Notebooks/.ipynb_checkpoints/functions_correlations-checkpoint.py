### Corr -> Frag ###

import matplotlib.pyplot as plt
import scipy.stats
#from scipy.stats import chisqprob
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
            for izt in range (0,NzT-ZT_OFF_PLOT):
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
            Dict["%s_CSR_Errors"%(SYS)][ipt] = np.sqrt((Dict["%s_CSR_Errors"%(SYS)][ipt]/purity[ipt])**2 + (1-purity[ipt]*Dict["%s_CBR"%(SYS)][ipt]/purity[ipt]**2))
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
                #Recalculate ZYAM for decay-subtracted Correlations
                ZYAM_Cs,ZYAM_Cs_Error = ZYAM_Line(Dict["%s_CSR"%(SYS)][ipt][izt],
                                                  Dict["%s_CSR_Errors"%(SYS)][ipt][izt])

                #print("%s ZYAM %i %1.4f"%(SYS,izt,ZYAM_Cs))
                Dict["%s_CSR"%(SYS)][ipt][izt] = Dict["%s_CSR"%(SYS)][ipt][izt] - ZYAM_Cs
                Dict["%s_Uncorr_Error"%(SYS)][ipt][izt] = ZYAM_Cs_Error
                
        np.save("npy_files/%s_%s_%s_Cs"%(Shower,description_string,SYS),Dict["%s_CSR"%(SYS)])
        np.save("npy_files/%s_%s_%s_Cs_Errors"%(Shower,description_string,SYS),Dict["%s_CSR_Errors"%(SYS)])
        np.save("npy_files/%s_%s_%s_Cs_Uncorr_Error"%(Shower,description_string,SYS),Dict["%s_Uncorr_Error"%(SYS)])                                                                
    
                    
def Plot_pp_pPb_Cs(Dict):
    
    Quad_UE = True
    
    for ipt in range (len(Dict["p-Pb_CSR"])):
        #if (ipt > 0): continue
        #ipt = ipt+2
        #plt.figure(figsize=(10,7))
        fig = plt.figure(figsize=(24,12))
        if (NzT==7):
            fig = plt.figure(figsize=(22,18))
        for izt in range (NzT-ZT_OFF_PLOT):

            if (NzT ==4):
                ax = fig.add_subplot(2,2,izt+1)
            elif (NzT ==6):
                ax = fig.add_subplot(2,3,izt+1)
            elif (NzT >=7):
                ax = fig.add_subplot(3,3,izt+1)
            elif (NzT >=12):
                ax = fig.add_subplot(3,5,izt+1)

            pPb = plt.errorbar(delta_phi_centers,Dict["p-Pb_CSR"][ipt][izt],xerr=phi_width,yerr=Dict["p-Pb_CSR_Errors"][ipt][izt],fmt='bo',capsize=4,markersize=11)
            pp = plt.errorbar(delta_phi_centers,Dict["pp_CSR"][ipt][izt],xerr=phi_width,yerr=Dict["pp_CSR_Errors"][ipt][izt],fmt='ro',capsize=4,markersize=11)

            if(Use_MC):
                MC = plt.errorbar(delta_phi_centers,["MC_CSR"][ipt][izt],xerr=phi_width,yerr=["MC_CSR_Errors"][ipt][ztb],fmt='go',capsize=4,markersize=11)


            if (izt>2):
                plt.xlabel(r'|$\Delta \varphi$|',fontsize=28)
            plt.xticks(fontsize=18)
            plt.xlim((0.39269908169872414,3.14159))
            if (izt%3 == 0):
                plt.ylabel(r'$1/N_{\gamma} \: \: \mathrm{d}N/\mathrm{d}\Delta \eta$',fontsize=28)
            plt.yticks(fontsize=18)
            plt.axhline(y=0,color='gray',linestyle='--',linewidth=1.3,alpha=0.8)        

            if not(Quad_UE):
                pp_UE = ax.fill_between(ue_error_bar,-Dict["pp_Uncorr_Error"][ipt][izt][0],
                Dict['pp_Uncorr_Error'][ipt][izt][0],facecolor='red',alpha=0.35) #Other for pp
                
                pPb_UE = ax.fill_between(ue_error_bar,-Dict["p-Pb_Uncorr_Error"][ipt][izt][0],
                Dict['p-Pb_Uncorr_Error'][ipt][izt][0],facecolor='blue',alpha=0.35)#One for p-Pb
            
            else:
                pp_UE = Dict["pp_Uncorr_Error"][ipt][izt][0]/(ZYAM_Max_i-ZYAM_Min_i)
                pPb_UE = Dict["p-Pb_Uncorr_Error"][ipt][izt][0]/(ZYAM_Max_i-ZYAM_Min_i)
                UE_Val =math.sqrt(pp_UE**2 + pPb_UE**2)
                Combined_UE = ax.fill_between(ue_error_bar,-UE_Val,UE_Val,facecolor='purple',alpha=0.35)
                
            Int_Window = ax.axvline(x=dPhi_Bins[-N_Phi_Integrate-1],linestyle='--',color="lightgreen",alpha=1.0)
            
            #MC_UE = ax.fill_between(ue_error_bar,-Dict["MC_Uncorr_Error"][ipt][ztb][0],Dict["MC_Uncorr_Error"][ipt][ztb][0],facecolor='green',alpha=0.35)#One for p-Pb
            
            
            Chi2,NDF,Pval = Get_pp_pPb_List_Chi2(Dict["p-Pb_CSR"][ipt][izt],Dict["p-Pb_CSR_Errors"][ipt][izt],Dict["p-Pb_Uncorr_Error"][ipt][izt],
                                        Dict["pp_CSR"][ipt][izt],Dict["pp_CSR_Errors"][ipt][izt],Dict["pp_Uncorr_Error"][ipt][izt])
                        
            plt.annotate("$\chi^2$ = %1.1f, ndf = %i, p = %1.2f"%(Chi2,NDF,Pval), xy=(0.99, 0.06), xycoords='axes fraction', ha='right', va='top', fontsize=16)

            if(Use_MC):
                leg = plt.legend([pp,pPb,MC,pp_UE,pPb_UE,MC_UE],['pp $\sqrt{s}= 5$ TeV (stat. error)',
                    'p-Pb $\sqrt{s_{\mathrm{_{NN}}}}=5$ TeV (stat. error)','pythia GJ $\sqrt{s}=5$ TeV (stat. error)', 
                    'pp UE Error', 'p-Pb UE Error','pythia UE Error'],
                    loc = "upper left",fontsize=16,frameon=False,numpoints=1)
            else:    
                if not(Quad_UE):
                    leg = plt.legend([pp,pPb,pp_UE,pPb_UE],['pp $\sqrt{s}= 5$ TeV (stat. error)',
                    'p-Pb $\sqrt{s_{\mathrm{_{NN}}}}=5$ TeV (stat. error)', 'pp UB Error', 'p-Pb UB Error'],
                    loc = "upper left",fontsize=16,frameon=False,numpoints=1)
                else:
                    leg = plt.legend([pp,pPb,Combined_UE],['pp $\sqrt{s}= 5$ TeV (stat. error)',
                    'p-Pb $\sqrt{s_{\mathrm{_{NN}}}}=5$ TeV (stat. error)', 'UB Error'],
                    loc = "upper left",fontsize=16,frameon=False,numpoints=1)

            plt.annotate(r'%1.2f < $z_\mathrm{T}$ < %1.2f'%(zTbins[izt],zTbins[izt+1]), xy=(0.05, 0.68), xycoords='axes fraction', ha='left', va='top', fontsize=16)
            
            if (len(Dict["p-Pb_CSR"]) > 2):
                plt.annotate(r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1]), xy=(0.05, 0.6), xycoords='axes fraction', ha='left', va='top', fontsize=16)
            else:
                plt.annotate(r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[0],pTbins[N_pT_Bins]), xy=(0.05, 0.6), xycoords='axes fraction', ha='left', va='top', fontsize=16)
            
            leg.set_title("ALICE Work in Progress")
            leg._legend_box.align = "left"
            plt.setp(leg.get_title(),fontsize=22)
            fig.tight_layout()

        plt.show()
        #fig.savefig('pics/Gamma_hadron_corr_zT_%i.pdf'%(ztb))
        fig.savefig('pics/%s/%s_Gamma_hadron_corr_%i.pdf'%(Shower,Shower,ipt))
        
        
        #Above can be adapted for pp & PbPb comparisons
        
def Cs_Weighted_Average(Dict):
    
    #Change this to empty dictionary with keys taking SYS. Then just loop over systems
    
    Averaged_pPb = np.zeros((1,NzT,N_dPhi_Bins))
    Averaged_pPb_Error = np.zeros((1,NzT,N_dPhi_Bins))
    Averaged_UB_pPb_Error = np.zeros((1,NzT,N_dPhi_Bins))
    
    Averaged_pp = np.zeros((1,NzT,N_dPhi_Bins))
    Averaged_pp_Error = np.zeros((1,NzT,N_dPhi_Bins))
    Averaged_UB_pp_Error = np.zeros((1,NzT,N_dPhi_Bins))
    #Dimension of size 1 for pT loop compatibility in plotters...
    
    for izt in range (NzT):
        for dphi in range(N_dPhi_Bins):
            for ipt in range (N_pT_Bins):

                pPb_Cs = Dict["p-Pb_CSR"][ipt][izt][dphi]
                pPb_Cs_Error = Dict["p-Pb_CSR_Errors"][ipt][izt][dphi]
                pp_Cs = Dict["pp_CSR"][ipt][izt][dphi]
                pp_Cs_Error = Dict["pp_CSR_Errors"][ipt][izt][dphi]
    
                if (pp_Cs_Error == 0 or pPb_Cs_Error == 0):
                    continue
                pPb_Cs_Weight = 1/(pPb_Cs_Error**2) #check if zero
                pp_Cs_Weight = 1/(pp_Cs_Error**2)
                
                Averaged_pPb[0][izt][dphi] += pPb_Cs_Weight*pPb_Cs
                Averaged_pp[0][izt][dphi] += pp_Cs_Weight*pp_Cs
                
                Averaged_pPb_Error[0][izt][dphi] += pPb_Cs_Weight
                Averaged_pp_Error[0][izt][dphi] += pp_Cs_Weight
                
                Averaged_UB_pPb_Error[0][izt][dphi] += (Dict["p-Pb_Uncorr_Error"][ipt][izt][0])**2
                Averaged_UB_pp_Error[0][izt][dphi] += (Dict["pp_Uncorr_Error"][ipt][izt][0])**2
                
        Averaged_pPb[0][izt] = Averaged_pPb[0][izt]/Averaged_pPb_Error[0][izt]
        Averaged_pp[0][izt] = Averaged_pp[0][izt]/Averaged_pp_Error[0][izt]
        
        Averaged_pPb_Error[0][izt] = np.sqrt(1/Averaged_pPb_Error[0][izt])
        Averaged_pp_Error[0][izt] = np.sqrt(1/Averaged_pp_Error[0][izt])
        

        Averaged_UB_pPb_Error[0][izt] = np.sqrt(Averaged_UB_pPb_Error[0][izt])/N_pT_Bins
        Averaged_UB_pp_Error[0][izt] = np.sqrt(Averaged_UB_pp_Error[0][izt])/N_pT_Bins
        
    
    #for izt in range(NzT):
    #    for ipt in range(N_pT_Bins):
    #        print Dict["p-Pb_Uncorr_Error"][ipt][izt][0]
    #        #Averaged_UB_Error[izt] += np.sqrt(Dict["p-Pb_Uncorr_Error"][ipt][izt][0]**2 + Dict["pp_Uncorr_Error"][ipt][izt][0]**2)
    #        Averaged_UB_pPb_Error[izt] += (Dict["p-Pb_Uncorr_Error"][ipt][izt][0])**2
    #        Averaged_UB_pp_Error[izt] += (Dict["pp_Uncorr_Error"][ipt][izt][0])**2
    #    
    #    #Averaged_UB_Error[izt] = Averaged_UB_Error[izt]/N_pT_Bins
    #    Averaged_UB_pPb_Error[izt] = np.sqrt(Averaged_UB_pPb_Error[izt])/N_pT_Bins
    #    Averaged_UB_pp_Error[izt] = np.sqrt(Averaged_UB_pp_Error[izt])/N_pT_Bins

        
    Keys = []
    Corr_Arrays = []
    
    #Keys.append("Combined_p-Pb_Cs")
    #Keys.append("Combined_p-Pb_Cs_Errors")
    #Keys.append("Combined_p-Pb_Cs_Uncorr_Error")
    
    Keys.append("p-Pb_CSR")
    Keys.append("p-Pb_CSR_Errors")
    Keys.append("p-Pb_Uncorr_Error")

    #Keys.append("Combined_pp_Cs")
    #Keys.append("Combined_pp_Cs_Errors")
    #Keys.append("Combined_pp_Cs_Uncorr_Error")
    
    Keys.append("pp_CSR")
    Keys.append("pp_CSR_Errors")
    Keys.append("pp_Uncorr_Error")

    Corr_Arrays.append(Averaged_pPb)
    Corr_Arrays.append(Averaged_pPb_Error)
    Corr_Arrays.append(Averaged_UB_pPb_Error)
    
    Corr_Arrays.append(Averaged_pp)
    Corr_Arrays.append(Averaged_pp_Error)
    Corr_Arrays.append(Averaged_UB_pp_Error)
    
    Weighted_Corr = dict(zip(Keys,Corr_Arrays))
    
    Save_Avg_Cs_npy(Weighted_Corr)
    
    return Weighted_Corr

        
def Save_Avg_Cs_npy(Corr):
    for SYS in Systems:
        np.save("npy_files/%s_%s_Combined_%s_Cs"%(Shower,description_string,SYS),Corr["%s_CSR"%(SYS)])
        np.save("npy_files/%s_%s_Combined_%s_Cs_Errors"%(Shower,description_string,SYS),Corr["%s_CSR_Errors"%(SYS)])
        np.save("npy_files/%s_%s_Combined_%s_Cs_Uncorr_Error"%(Shower,description_string,SYS),Corr["%s_Uncorr_Error"%(SYS)])
                
                
def Compare_Cs_Averages(strings,string_descrp_list,colors):
    
    #shapes = ["o","x","s"]
    for SYS in Systems:
        fig = plt.figure(figsize=(22,18))
        
        for (string,string_descr,colr) in zip(strings,string_descrp_list,colors):
            print(string)
            
            CS_Avg = np.load("npy_files/%s_%s_Combined_%s_Cs.npy"%(Shower,string,SYS))[0]
            CS_Avg_Err = np.load("npy_files/%s_%s_Combined_%s_Cs_Errors.npy"%(Shower,string,SYS))[0]
            CS_Avg_Uncorr_Err = np.load("npy_files/%s_%s_Combined_%s_Cs_Uncorr_Error.npy"%(Shower,string,SYS))[0]
            
            for izt in range(NzT-ZT_OFF_PLOT):
                
                if (NzT ==4):
                    ax = fig.add_subplot(2,2,izt+1)
                elif (NzT ==6):
                    ax = fig.add_subplot(2,3,izt+1)
                elif (NzT ==7):
                    ax = fig.add_subplot(3,3,izt+1)
   
                N_dPhi = len(CS_Avg)+1
                dPhi_Centers = [i*math.pi/N_dPhi+math.pi/N_dPhi/2 for i in range(0,N_dPhi)] #skip first dPhi bin to avoid Isolation
                    
                plt.errorbar(dPhi_Centers,CS_Avg[izt],xerr=phi_width,yerr=CS_Avg_Err[izt],fmt='o',color = colr,capsize=4,markersize=11,label = "average %s"%(string_descr))
                
                plt.xlim((0.39269908169872414,3.14159))
                
                #Labels
                if (izt>2):
                    plt.xlabel(r'|$\Delta \varphi$|',fontsize=28)
                if (izt%3 == 0):
                    plt.ylabel(r'$1/N_{\gamma} \: \: \mathrm{d}N/\mathrm{d}\Delta \eta$',fontsize=28)
                
                #Chi2,NDF,Pval = Get_pp_pPb_List_Chi2(Dict["Combined_p-Pb_Cs"][izt],Dict["Combined_p-Pb_Cs_Errors"][izt],Dict["Combined_p-Pb_Cs_Uncorr_Error"][izt],
                #                        Dict["Combined_pp_Cs"][izt],Dict["Combined_pp_Cs_Errors"][izt],Dict["Combined_pp_Cs_Uncorr_Error"][izt])
            
                
                leg = plt.legend(numpoints=1,frameon=False,loc="best")
                leg.set_title("%s :%1.2f < $z_\mathrm{T}$ < %1.2f"%(SYS,zTbins[izt],zTbins[izt+1]))
                plt.setp(leg.get_title(),fontsize=18)
                    
            
def Compare_Cs_pTBins():
    
    #shapes = ["o","x","s"]
    for SYS in Systems:
        fig = plt.figure(figsize=(22,18))
        
        CS = np.load("npy_files/%s_%s_%s_Cs.npy"%(Shower,description_string,SYS))
        CS_ERR = np.load("npy_files/%s_%s_%s_Cs_Errors.npy"%(Shower,description_string,SYS))
        CS_Uncorr_ERR = np.load("npy_files/%s_%s_%s_Cs_Uncorr_Error.npy"%(Shower,description_string,SYS))    
        
        CS_Avg = np.load("npy_files/%s_%s_Combined_%s_Cs.npy"%(Shower,description_string,SYS))[0]
        CS_Avg_Err = np.load("npy_files/%s_%s_Combined_%s_Cs_Errors.npy"%(Shower,description_string,SYS))[0]
        CS_Avg_Uncorr_Err = np.load("npy_files/%s_%s_Combined_%s_Cs_Uncorr_Error.npy"%(Shower,description_string,SYS))[0]
            
        for izt in range(NzT-ZT_OFF_PLOT):
                
            if (NzT ==4):
                ax = fig.add_subplot(2,2,izt+1)
            elif (NzT ==6):
                ax = fig.add_subplot(2,3,izt+1)
            elif (NzT ==7):
                ax = fig.add_subplot(3,3,izt+1)
                
            plt.errorbar(delta_phi_centers,CS_Avg[izt],xerr=phi_width,yerr=CS_Avg_Err[izt],fmt='ko',capsize=4,markersize=11,
            label = '%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[0],pTbins[N_pT_Bins]))
                
            for ipt in range(N_pT_Bins):
                plt.errorbar(delta_phi_centers,CS[ipt][izt],xerr=phi_width,yerr=CS_ERR[ipt][izt],fmt='o',capsize=4,markersize=11,
                label = '%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'
                %(pTbins[ipt],pTbins[ipt+1]))
                
                plt.xlim((0.39269908169872414,3.14159))
                
                #Labels
                if (izt>2):
                    plt.xlabel(r'|$\Delta \varphi$|',fontsize=28)
                if (izt%3 == 0):
                    plt.ylabel(r'$1/N_{\gamma} \: \: \mathrm{d}N/\mathrm{d}\Delta \eta$',fontsize=28)                
                
                leg = plt.legend(numpoints=1,frameon=False,loc="best")
                leg.set_title("%s :%1.2f < $z_\mathrm{T}$ < %1.2f"%(SYS,zTbins[izt],zTbins[izt+1]))
                plt.setp(leg.get_title(),fontsize=18)           
            
                
        
def Integrate_Away_Side(Phi_array,Phi_Errors,LE_Error):
    
    Use_Uncorr_Error = True
    FF_zt = np.zeros((len(Phi_array), NzT))
    FF_zt_Errors = np.zeros((len(Phi_array), NzT))
    
    if Use_Uncorr_Error:
        LE_Error = LE_Error/(dPhi_Width*(ZYAM_Max_i-ZYAM_Min_i))
    
    for ipt in range(len(Phi_array)):
        
        for izt in range(0, NzT):
            
            zT_width = zTbins[izt+1]-zTbins[izt]
            #zT_width = 1
            
            #temp_phi = Phi_array[ipt][izt][(len(Phi_Errors[ipt][izt])-N_Phi_Integrate):]/(dPhi_Width*N_Phi_Integrate)
            temp_phi = Phi_array[ipt][izt][-N_Phi_Integrate:]/(dPhi_Width*N_Phi_Integrate)
            #print(temp_phi)
            
            FF_zt[ipt][izt] = temp_phi.sum()/zT_width
            temp_error = (Phi_Errors[ipt][izt][-N_Phi_Integrate:]/(dPhi_Width*N_Phi_Integrate))**2
            
            if (Use_Uncorr_Error):
                FF_zt_Errors[ipt][izt] = (math.sqrt(temp_error.sum() + (LE_Error[ipt][izt][0])**2))/zT_width
            else:
                FF_zt_Errors[ipt][izt] = math.sqrt(temp_error.sum())/zT_width
            
    return FF_zt, FF_zt_Errors


def Get_Fragmentation(Dict,Use_Avg_Cs=False):
    
    Keys = []
    
    for SYS in Systems:
        
        if not Use_Avg_Cs:
            Keys.append("%s_FF"%(SYS))
            Keys.append("%s_FF_Errors"%(SYS))
            Keys.append("%s_purity_FF_Errors"%(SYS))
        
        if Use_Avg_Cs:
            Keys.append("%s_Combined_FF"%(SYS))
            Keys.append("%s_Combined_FF_Errors"%(SYS))
            Keys.append("%s_purity_Uncertainty"%(SYS))
            
    
    #FF_Vals = np.zeros((len(Keys),N_pT_Bins,NzT))
    FF_Vals = []

    for index,SYS in enumerate(Systems):
        
        temp_FF, temp_FF_Errors = Integrate_Away_Side(Dict["%s_CSR"%(SYS)],Dict["%s_CSR_Errors"%(SYS)],Dict["%s_Uncorr_Error"%(SYS)])
        temp_purity_Errors = []
        
        for ipt in range(len(Dict["%s_CSR"%(SYS)])): 
            temp_purity_Errors.append(temp_FF[ipt]*(purity_Uncertainty[ipt]/purity[ipt]))  # abs. FF purity uncertainty
        
        if (Use_Avg_Cs):
            temp_FF,temp_FF_Errors,temp_purity_Errors = temp_FF[0],temp_FF_Errors[0],temp_purity_Errors[0]
            
        FF_Vals.append(temp_FF)
        FF_Vals.append(temp_FF_Errors)
        FF_Vals.append(np.asarray(temp_purity_Errors))
        
        #FF_Dict["%s_FF"%(SYS)], FF_Dict["%s_FF_Errors"%(SYS)] = temp_FF, temp_FF_Errors
        #FF_Dict["%s_purity_FF_Errors"%(SYS)] = np.asarray(temp_purity_Errors)
    
    FF_Dict = dict(zip(Keys,FF_Vals))
    
    for SYS in Systems:
        
        if (Use_Avg_Cs):
            
            if (Use_Weights):
                np.save("npy_files/%s_%s_Fragmentation_Functions_Combined_Cs.npy"%(Shower,SYS),FF_Dict["%s_Combined_FF"%(SYS)])
                np.save("npy_files/%s_%s_Fragmentation_Function_Errors_Combined_Cs.npy"%(Shower,SYS),FF_Dict["%s_Combined_FF_Errors"%(SYS)])
                np.save("npy_files/%s_%s_FF_purity_Uncertainty_Combined_Cs.npy"%(Shower,SYS),FF_Dict["%s_purity_Uncertainty"%(SYS)])
            
            else:
                np.save("npy_files/%s_%s_Fragmentation_Functions_Unweight_Combined_Cs.npy"%(Shower,SYS),FF_Dict["%s_Combined_FF"%(SYS)])
                np.save("npy_files/%s_%s_Fragmentation_Function_Errors_Unweight_Combined_Cs.npy"%(Shower,SYS),FF_Dict["%s_Combined_FF_Errors"%(SYS)])
                np.save("npy_files/%s_%s_FF_purity_Uncertainty_Unweight_Combined_Cs.npy"%(Shower,SYS),FF_Dict["%s_purity_Uncertainty"%(SYS)])
        
        else:
            
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
    
    
    for ipt in range(len(FF_Dict["p-Pb_FF"])):

        zt_box = np.ones(NzT) * 0.03
        #pPb_bar = plt.bar(zT_centers, pPb_sys[ipt]+pPb_sys[ipt], bottom=(pPb_FF[ipt])-pPb_sys[ipt], width=zt_box, align='center',edgecolor="blue",color='white',)
        #pp_bar = plt.bar(zT_centers, pp_sys[ipt]+pp_sys[ipt], bottom=pp_FF[ipt]-pp_sys[ipt], width=zt_box, align='center',edgecolor="red",color='white',)

        zT_max = 0
        for izt in range(0, NzT-ZT_OFF_PLOT):
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
