from default_values import *
import matplotlib.pyplot as plt
import matplotlib
from ROOT import TGraphErrors

def FF_Ratio(FF_Dict):
    
    fig = plt.figure(figsize=(17,13/(4-N_pT_Bins+1)))
    
    for ipt in range (N_pT_Bins):
        #if (ipt > 0): continue

        Ratio = FF_Dict["p-Pb_FF"][ipt]/FF_Dict["pp_FF"][ipt]

        #Stat. Error in ratio
        Ratio_Error = np.sqrt((FF_Dict["p-Pb_FF_Errors"][ipt]/FF_Dict["p-Pb_FF"][ipt])**2 + (FF_Dict["pp_FF_Errors"][ipt]/FF_Dict["pp_FF"][ipt])**2)*Ratio
        
        ax = fig.add_subplot((N_pT_Bins/2),2,ipt+1)

        #Sys_Plot = plt.bar(zT_centers, Ratio_Systematic+Ratio_Systematic, bottom=Ratio-Ratio_Systematic, width=zt_box, align='center',edgecolor="black",color='white',)

        empt4, = plt.plot([], [],' ') #Labels pT range

        Ratio_Plot = plt.errorbar(zT_centers, Ratio, yerr=Ratio_Error,xerr=zT_widths, fmt='ko',capsize=3, ms=6,lw=1)
        plt.xlabel("${z_\mathrm{T}} = p_\mathrm{T}^{\mathrm{h}}/p_\mathrm{T}^\gamma$",fontsize=20)
        plt.ylabel(r"$\frac{\mathrm{p-Pb}}{\mathrm{pp}}$",fontsize=20)
        
        leg = plt.legend([Ratio_Plot,empt4],["Statistical Error",r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1])],frameon=False,numpoints=1,title=' ',prop={'size':18})
        leg.set_title("ALICE Work in Progress\n ")
        plt.setp(leg.get_title(),fontsize=20)
        
        plt.xlim(xmin = 0.0,xmax=0.7)
        plt.axhline(y=1, color='r', linestyle='--')

        plt.gcf()
        plt.savefig("pics/%s_pp_FFunction_%i.pdf"%(Shower,ipt), bbox_inches='tight')

    plt.gcf()
    plt.savefig("pics/%s_pp_FFunction.pdf"%(Shower), bbox_inches='tight')
    plt.show()
    

def Overlay_pT_FF(FF_Dict):
    
    for SYS in Systems:
        
        plt.figure(figsize=(10,7))
        
        for ipt in range (N_pT_Bins):
            if(ipt>2):
                #plt.errorbar(zT_centers[:5], pPb_FF[ipt][:5],xerr=zT_widths[:5],yerr=pPb_FF_Errors[ipt][:5],linewidth=1, fmt='o',capsize=1,label=r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1]))
                plt.errorbar(zT_centers[1:5], FF_Dict["%s_FF"%(SYS)][ipt][1:5],xerr=zT_widths[1:5],yerr=FF_Dict["%s_FF_Errors"%(SYS)][ipt][1:5],linewidth=1, fmt='o',capsize=1,label=r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1]))

            #elif(ipt<2):
            #    plt.errorbar(zT_centers[1:], pPb_FF[ipt][1:],xerr=zT_widths[1:],yerr=pPb_FF_Errors[ipt][1:],linewidth=1, fmt='o',capsize=1,label=r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1]))
            else:
                plt.errorbar(zT_centers[1:], FF_Dict["%s_FF"%(SYS)][ipt][1:],xerr=zT_widths[1:],yerr=FF_Dict["%s_FF_Errors"%(SYS)][ipt][1:],linewidth=1, fmt='o',capsize=1,label=r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1]))
                #plt.errorbar(zT_centers, pPb_FF[ipt],xerr=zT_widths,yerr=pPb_FF_Errors[ipt],linewidth=1, fmt='o',capsize=1,label=r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1]))
            #plt.errorbar(zT_centers, pPb_FF[ipt],xerr=zT_widths,yerr=pPb_FF_Errors[ipt],linewidth=1, fmt='o',capsize=1,label=r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1]))

            plt.yscale('log')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
            plt.ylabel(r"$\frac{1}{N_{\mathrm{\gamma}}}\frac{\mathrm{d}N}{\mathrm{d}z_{\mathrm{T}} \mathrm{d}\Delta\eta}$",fontsize=20)
            plt.xlabel("${z_\mathrm{T}} = p_\mathrm{T}^\mathrm{h}/p_\mathrm{T}^\mathrm{\gamma}$",fontsize=20)
            #plt.xlim(xmin = 0.1,xmax=0.7)
            plt.ylim(ymin = 0.001,ymax=20)

        leg = plt.legend(numpoints=1,frameon=False)
        leg.set_title("ALICE Work in Progress\n  %s $\sqrt{s_{\mathrm{_{NN}}}} = $ 5 TeV"%(SYS))
        plt.setp(leg.get_title(),fontsize=18)

        plt.title(r'Integrated $\mathrm{\gamma}$-Hadron Correlation: $2\pi/3 < \Delta\varphi < \pi, |\Delta\eta| < %1.1f$ '%(eta_max),fontdict = {'fontsize' : 20})
        plt.gcf()
        plt.savefig("pics/All_pT_FFunction_pPb_%s.pdf"%(Shower), bbox_inches='tight')
        plt.show()

        if (SYS == "p-Pb"):
            print("                              PROTON-LEAD:")
            
        elif (SYS == "pp"):
            print("                             PROTON-PROTON:")
        
        print("Central Values")
        print(FF_Dict["%s_FF"%(SYS)])
        print("Statistical Errors (Relative)")
        print(FF_Dict["%s_FF_Errors"%(SYS)]/FF_Dict["%s_FF"%(SYS)])
        print("Relative Uncertainty from Purity")
        print(FF_Dict["%s_purity_FF_Errors"%(SYS)])
        print("\n")
        
def Weighted_Average(FF,FF_Errors,purity_FF_Errors):

    Combined = np.zeros(len(FF[0]))
    Combined_Errors = np.zeros(len(FF_Errors[0]))
    purity_Combined_Errors = np.zeros(len(purity_FF_Errors[0]))

    Rel_Stat_Erorr = FF_Errors/FF
    Rel_Purity_Error = purity_FF_Errors/FF
    
    for izt in range (NzT):
        weight_sum = 0

        for ipt in range(0,N_pT_Bins):

            if (ipt>2 and izt>4): continue #Tracking efficiency limit
            #if (izt<1):
            #    continue
            #if(ipt<2 and izt<1): continue
            
            Combined[izt] += 1/((Rel_Stat_Erorr[ipt][izt]**2) + (Rel_Purity_Error[ipt][izt]**2)) * FF[ipt][izt] 
            weight_sum += 1/((Rel_Stat_Erorr[ipt][izt]**2) + (Rel_Purity_Error[ipt][izt]**2))
            
            #OLD METHOD: INVERSE VARIANCE WEIGHTING with ABSOLUTE UNCERTAINTIES
            #Combined[izt]+= 1/((FF_Errors[ipt][izt])**2 + (purity_FF_Errors[ipt][izt]**2)) * FF[ipt][izt] #Weighting by stat. error**2 = sqrt(N)**2
            #weight_sum += 1/((FF_Errors[ipt][izt]**2) + (purity_FF_Errors[ipt][izt]**2))

        Combined[izt] = Combined[izt]/weight_sum

    for ipt in range (N_pT_Bins):
        Combined_Errors += FF_Errors[ipt]**2
        purity_Combined_Errors += purity_FF_Errors[ipt]**2

    Combined_Errors = np.sqrt(Combined_Errors)/N_pT_Bins
    purity_Combined_Errors = np.sqrt(purity_Combined_Errors)/N_pT_Bins
    
    return(Combined,Combined_Errors,purity_Combined_Errors)

def Average_FF(FF_Dict):
    
    Keys = []
    Combined_FF = []
    
    for SYS in Systems:
        Keys.append("%s_Combined_FF"%(SYS))
        Keys.append("%s_Combined_FF_Errors"%(SYS))
        Keys.append("%s_purity_Uncertainty"%(SYS))
        
        temp_combined,temp_combined_Errors,temp_purity_Uncert = Weighted_Average(FF_Dict["%s_FF"%(SYS)],FF_Dict["%s_FF_Errors"%(SYS)],FF_Dict["%s_purity_FF_Errors"%(SYS)])
        Combined_FF.append(temp_combined)
        Combined_FF.append(temp_combined_Errors)
        Combined_FF.append(temp_purity_Uncert)
        
    return dict(zip(Keys,Combined_FF))


def Plot_pp_pPb_Avg_FF(Comb_Dict):
    
    Colors = ["red","blue","black"]
    plt.figure(figsize=(10,7))
    
    for SYS,sys_col in zip(Systems,Colors):
        
        Efficiency_Uncertainty = 0.05*Comb_Dict["%s_Combined_FF"%(SYS)]
        Sys_Uncertainty = np.sqrt(Efficiency_Uncertainty**2 + Comb_Dict["%s_purity_Uncertainty"%(SYS)]**2)


        #---------------- Plot ------------------------------#

        #plt.errorbar(zT_centers, Comb_Dict["%s_Combined_FF"%(SYS)],xerr=zT_widths,yerr=Comb_Dict["%s_Combined_FF_Errors"%(SYS)],
            #linewidth=1, fmt='o',color=sys_col,capsize=1,label=r' %s %1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(SYS,pTbins[0],pTbins[N_pT_Bins]))

        #Sys_Plot_pp = plt.bar(zT_centers, Sys_Uncertainty+Sys_Uncertainty, bottom=Comb_Dict["%s_Combined_FF"%(SYS)]-Sys_Uncertainty, 
        #                      width=zt_box, align='center',edgecolor="black",color='white',)
        
        plt.errorbar(zT_centers[ZT_OFF_PLOT:], Comb_Dict["%s_Combined_FF"%(SYS)][ZT_OFF_PLOT:],xerr=zT_widths[ZT_OFF_PLOT:],
            yerr=Comb_Dict["%s_Combined_FF_Errors"%(SYS)][ZT_OFF_PLOT:],linewidth=1, fmt='o',color=sys_col,capsize=1,
            label=r' %s %1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(SYS,pTbins[0],pTbins[N_pT_Bins]))

        Sys_Plot_pp = plt.bar(zT_centers[ZT_OFF_PLOT:], Sys_Uncertainty[ZT_OFF_PLOT:]+Sys_Uncertainty[ZT_OFF_PLOT:], 
            bottom=Comb_Dict["%s_Combined_FF"%(SYS)][ZT_OFF_PLOT:]-Sys_Uncertainty[ZT_OFF_PLOT:],width=zt_box[ZT_OFF_PLOT:], align='center',color='white',)

        plt.yscale('log')                                                                                                                                                                                                                                                              
        plt.ylabel(r"$\frac{1}{N_{\mathrm{\gamma}}}\frac{\mathrm{d}N}{\mathrm{d}z_{\mathrm{T}} \mathrm{d}\Delta\eta}$",fontsize=20)
        plt.xlabel("${z_\mathrm{T}} = p_\mathrm{T}^\mathrm{h}/p_\mathrm{T}^\mathrm{\gamma}$",fontsize=20)
        #plt.xlim(xmin = 0.1,xmax=0.7)
        plt.ylim(ymin = 0.01,ymax=20)

    leg = plt.legend(numpoints=1,frameon=False)
    leg.set_title("ALICE Work in Progress\n  $\sqrt{s_{\mathrm{_{NN}}}} = $ 5 TeV")
    plt.setp(leg.get_title(),fontsize=18)

    plt.title(r'Integrated $\mathrm{\gamma}$-Hadron Correlation: $2\pi/3 < \Delta\varphi < \pi, |\Delta\eta| < %1.1f$ '%(eta_max),fontdict = {'fontsize' : 19})
    plt.gcf()
    plt.savefig("pics/%s/Averaged_pT_FFunction_%s.pdf"%(Shower,Shower), bbox_inches='tight')
    plt.show()
    
    
    for SYS in Systems:
        print("                    %s Central Values:"%(SYS))
        print(Comb_Dict["%s_Combined_FF"%(SYS)][ZT_OFF_PLOT:])
        print("")
        print("                    %s Stat. Uncertainty:"%(SYS))
        print(Comb_Dict["%s_Combined_FF_Errors"%(SYS)][ZT_OFF_PLOT:])
        print("")
        
    print("                        LaTeX Table")
    
    pp_stat_min = np.amin(Comb_Dict["pp_Combined_FF_Errors"]/Comb_Dict["pp_Combined_FF"])*100
    pp_stat_max = np.amax(Comb_Dict["pp_Combined_FF_Errors"]/Comb_Dict["pp_Combined_FF"])*100
    pPb_stat_min = np.amin(Comb_Dict["p-Pb_Combined_FF_Errors"]/Comb_Dict["p-Pb_Combined_FF"])*100
    pPb_stat_max = np.amax(Comb_Dict["p-Pb_Combined_FF_Errors"]/Comb_Dict["p-Pb_Combined_FF"])*100
    
    pp_purity_min = np.amin(Comb_Dict["pp_purity_Uncertainty"]/Comb_Dict["pp_Combined_FF"])*100
    pp_purity_max = np.amax(Comb_Dict["pp_purity_Uncertainty"]/Comb_Dict["pp_Combined_FF"])*100
    pPb_purity_min = np.amin(Comb_Dict["p-Pb_purity_Uncertainty"]/Comb_Dict["p-Pb_Combined_FF"])*100
    pPb_purity_max = np.amax(Comb_Dict["p-Pb_purity_Uncertainty"]/Comb_Dict["p-Pb_Combined_FF"])*100
    
    print("Source   &  pp data & \pPb~data  \\\\")
    print("Statistical Uncertainty & {0}\%-{1}\% & {2}\%-{3}\% \\\\".format(int(pp_stat_min+0.5),
                        int(pp_stat_max+0.5),int(pPb_stat_min+0.5),int(pPb_stat_max+0.5)) )
    print("\hline")
    
    print("Purity & {0}\%-{1}\% & {2}\%-{3}\% \\\\".format(int(pp_purity_min+0.5),
        int(pp_purity_max+0.5),int(pPb_purity_min+0.5),int(pPb_purity_max+0.5)) )
    
    print("Tracking Efficiency &  5\% & 5\%  \\\\ ")
    
#      Source   &  pp data & \pPb~data  \\
#  Statistical uncertainty &  8 - 30 \% & 11 - 35\%  \\
#  \hline 
#  Purity (sys.) &  25\% & 25\%  \\
#  Purity (stat.)  & 6\% & 6\%\\
#  Tracking efficiency &  5\% & 5\%  \\
#  UE Estimate &  3-10\% & 3-9\%  \\

def pp_pPB_Avg_Ratio(Comb_Dict,pT_Start):
    
    pPb_Combined = Comb_Dict["p-Pb_Combined_FF"]
    pPb_Combined_Errors = Comb_Dict["p-Pb_Combined_FF_Errors"]
    pPb_purity_Uncertainty = Comb_Dict["p-Pb_purity_Uncertainty"]
    
    pp_Combined = Comb_Dict["pp_Combined_FF"]
    pp_Combined_Errors = Comb_Dict["pp_Combined_FF_Errors"]
    pp_purity_Uncertainty = Comb_Dict["pp_purity_Uncertainty"]
    
    Ratio = pPb_Combined/pp_Combined
    
    #Stat. Error in ratio
    Ratio_Error = np.sqrt((pPb_Combined_Errors/pPb_Combined)**2 + (pp_Combined_Errors/pp_Combined)**2)*Ratio
    Purity_Uncertainty = np.sqrt((pp_purity_Uncertainty/pPb_Combined)**2 + (pPb_purity_Uncertainty/pPb_Combined)**2)
    Efficiency_Uncertainty = np.ones(len(pPb_Combined))*0.05*math.sqrt(2)*Ratio

    if (CorrectedP):
        Ratio_Systematic = np.sqrt(Purity_Uncertainty**2 + Efficiency_Uncertainty**2)

    plt.figure(figsize=(10,7)) 


    Sys_Plot = plt.bar(zT_centers[ZT_OFF_PLOT:], Ratio_Systematic[ZT_OFF_PLOT:]+Ratio_Systematic[ZT_OFF_PLOT:], 
            bottom=Ratio[ZT_OFF_PLOT:]-Ratio_Systematic[ZT_OFF_PLOT:], width=zt_box[ZT_OFF_PLOT:], align='center',edgecolor="black",color='white',)
    #Sys_Plot = plt.bar(zT_centers, Ratio_Systematic+Ratio_Systematic, bottom=Ratio-Ratio_Systematic, width=zt_box, align='center',edgecolor="black",color='white',)

    empt4, = plt.plot([], [],' ')

    Sys_Box = [0.65,0.69]
    xfill = [0.65,0.7]

    Ratio_Plot = plt.errorbar(zT_centers[ZT_OFF_PLOT:], Ratio[ZT_OFF_PLOT:], yerr=Ratio_Error[ZT_OFF_PLOT:],xerr=zT_widths[ZT_OFF_PLOT:], fmt='ko',capsize=3, ms=6,lw=1)
    #Ratio_Plot = plt.errorbar(zT_centers, Ratio, yerr=Ratio_Error,xerr=zT_widths, fmt='ko',capsize=3, ms=6,lw=1)

    plt.xlabel("${z_\mathrm{T}} = p_\mathrm{T}^{\mathrm{h}}/p_\mathrm{T}^\gamma$",fontsize=20)
    plt.ylabel(r"$\frac{\mathrm{p-Pb}}{\mathrm{pp}}$",fontsize=20)
    plt.ylim((0, 2))
    plt.yticks(np.arange(-0, 2, step=0.2))
    plt.xlim(xmin = 0.0,xmax=0.7)
    plt.axhline(y=1, color='r', linestyle='--')

    ### ROOT LINEAR and CONSTANT FITS ###
    Ratio_TGraph = TGraphErrors()
    for izt in range (1,len(Ratio)):
        Ratio_TGraph.SetPoint(izt,zT_centers[izt],Ratio[izt])
        Ratio_TGraph.SetPointError(izt,0,Ratio_Error[izt])

    Ratio_TGraph.Fit("pol0","S")
    f = Ratio_TGraph.GetFunction("pol0")
    chi2_red  = f.GetChisquare()/f.GetNDF()
    pval = f.GetProb()
    p0 = f.GetParameter(0)
    p0e = f.GetParError(0)
    p0col = "blue"
    plt.text(0.01,1.9,"Constant Fit",color=p0col,fontsize=20,alpha=.7)
    plt.text(0.01,1.81,r"$p0 = {0:.2f} \pm {1:.2f}$".format(p0,p0e),color=p0col,fontsize=18,alpha=.7)
    plt.text(0.01,1.71,r"$\chi^2_{red} = %1.2f$"%(chi2_red),color=p0col,fontsize=18,alpha=.7)
    plt.text(0.01,1.63,r"$p_{val} = %1.2f$"%(pval),color=p0col,fontsize=18,alpha=.7)
    plt.fill_between(np.arange(0,1,0.1), p0+p0e, p0-p0e,color=p0col,alpha=.3)

    Ratio_TGraph.Fit("pol1","S")
    f2 = Ratio_TGraph.GetFunction("pol1")
    chi2_red  = f2.GetChisquare()/f2.GetNDF()
    pval = f2.GetProb()
    p0 = f2.GetParameter(0)
    p0e = f2.GetParError(0)
    p1 = f2.GetParameter(1)
    p1e = f2.GetParError(1)
    p1col = "cyan"
    plt.text(0.01,0.32,"Linear Fit",color=p1col,fontsize=20,alpha=.7)
    plt.text(0.01,0.23,r"$p1 = {0:.2f} \pm {1:.2f}$".format(p1,p1e),color=p1col,fontsize=18,alpha=.7)
    plt.text(0.01,0.13,r"$\chi^2_{red} = %1.2f$"%(chi2_red),color=p1col,fontsize=18,alpha=.7)
    plt.text(0.01,0.03,r"$p_{val} = %1.2f$"%(pval),color=p1col,fontsize=18,alpha=.7)
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = p0 + p1 * x_vals
    plt.plot(x_vals, y_vals, '--',color=p1col,linewidth=2)

    ### ROOT DONE ###


    leg = plt.legend([Ratio_Plot,empt4],["Statistical Error",r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[pT_Start],pTbins[N_pT_Bins])],frameon=False,numpoints=1,loc="best",title=' ',prop={'size':18})


    #leg = plt.legend(numpoints=1)
    leg.set_title("ALICE Work in Progress\n  $\sqrt{s_{\mathrm{_{NN}}}} = $ 5 TeV")
    plt.setp(leg.get_title(),fontsize=20)
    #plt.figtext(0.39,0.85,"ALICE Work in Progress\n  $\sqrt{s_{\mathrm{_{NN}}}} = $ 5 TeV",color='Black', fontsize=20)

    plt.gcf()
    plt.savefig("pics/Averaged_pT_FFunction_ratio_%s.pdf"%(Shower), bbox_inches='tight')
    plt.show()

    print("                Central Values:")
    print(Ratio[ZT_OFF_PLOT:])


    print("\n                Ratio Uncertainty from Purity:")
    print(Purity_Uncertainty[ZT_OFF_PLOT:])

    print("\n                Ratio Uncertainty from Single Track Efficiency:")
    print(Efficiency_Uncertainty[ZT_OFF_PLOT:])

    print("\n                Full Systematic Uncertainty:")

    print(Ratio_Systematic[ZT_OFF_PLOT:])
