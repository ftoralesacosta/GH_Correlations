from default_values import *
from functions_root_nparray import *
from functions_correlations import *
import matplotlib.pyplot as plt
import matplotlib
from ROOT import TGraphErrors
import scipy
import iminuit
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def FF_Ratio(FF_Dict):
    
    fig = plt.figure(figsize=(17,13))
    
    for ipt in range (len(FF_Dict["p-Pb_FF"])):

        Ratio = FF_Dict["p-Pb_FF"][ipt]/FF_Dict["pp_FF"][ipt]
        print(Ratio)
        #Stat. Error in ratio
        Ratio_Error = np.sqrt((FF_Dict["p-Pb_FF_Errors"][ipt]/FF_Dict["p-Pb_FF"][ipt])**2 + (FF_Dict["pp_FF_Errors"][ipt]/FF_Dict["pp_FF"][ipt])**2)*Ratio
        
        if (N_pT_Bins < 5):
            ax = fig.add_subplot(2,2,ipt+1)
        elif (N_pT_Bins >=5):
            ax = fig.add_subplot(3,2,ipt+1)

        #Sys_Plot = plt.bar(zT_centers, Ratio_Systematic+Ratio_Systematic, bottom=Ratio-Ratio_Systematic, width=zt_box, align='center',edgecolor="black",color='white',)

        empt4, = plt.plot([], [],' ') #Labels pT range

        Ratio_Plot = plt.errorbar(zT_centers, Ratio, yerr=Ratio_Error,xerr=zT_widths, fmt='ko',capsize=3, ms=6,lw=1)
        plt.xlabel("${z_\mathrm{T}} = p_\mathrm{T}^{\mathrm{h}}/p_\mathrm{T}^\gamma$",fontsize=20)
        plt.ylabel(r"$\frac{\mathrm{p-Pb}}{\mathrm{pp}}$",fontsize=20)
        
        plt.xlim(xmin = 0.0,xmax=zT_centers[NzT-ZT_OFF_PLOT])
        plt.ylim(ymin = -0.5, ymax=2.5)
        
        leg = plt.legend([Ratio_Plot,empt4],["Statistical Error",r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1])],frameon=False,numpoints=1,title=' ',prop={'size':18})
        leg.set_title("ALICE Work in Progress\n ")
        plt.setp(leg.get_title(),fontsize=20)
        
        plt.axhline(y=1, color='r', linestyle='--')
    
    plt.gcf()
    plt.savefig("pics/%s_pp_FFunction.pdf"%(Shower), bbox_inches='tight')
    plt.show()
    
    if (len(FF_Dict["p-Pb_FF"]) < 2): #Because Cs avge and FF avg methods the same for 1 pT bin 
        Ratio = np.save("npy_files/LO_Cs_Avg_FF_Ratio_%s.npy"%(description_string),Ratio)
        Ratio_Error = np.save("npy_files/LO_Cs_Avg_FF_Ratio_Errors_%s.npy"%(description_string),Ratio_Error)
    
def Overlay_pT_FF(FF_Dict):
    
    for SYS in Systems:
        
        plt.figure(figsize=(10,7))
        
        for ipt in range (N_pT_Bins):
            if(ipt>2):
                #plt.errorbar(zT_centers[:5], pPb_FF[ipt][:5],xerr=zT_widths[:5],yerr=pPb_FF_Errors[ipt][:5],linewidth=1, fmt='o',capsize=1,label=r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1]))
                plt.errorbar(zT_centers[ZT_OFF_PLOT:5], FF_Dict["%s_FF"%(SYS)][ipt][ZT_OFF_PLOT:5],xerr=zT_widths[ZT_OFF_PLOT:5],yerr=FF_Dict["%s_FF_Errors"%(SYS)][ipt][ZT_OFF_PLOT:5],linewidth=1, fmt='o',capsize=1,label=r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1]))

            #elif(ipt<2):
            #    plt.errorbar(zT_centers[1:], pPb_FF[ipt][1:],xerr=zT_widths[1:],yerr=pPb_FF_Errors[ipt][1:],linewidth=1, fmt='o',capsize=1,label=r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1]))
            else:
                plt.errorbar(zT_centers[ZT_OFF_PLOT:], FF_Dict["%s_FF"%(SYS)][ipt][ZT_OFF_PLOT:],xerr=zT_widths[ZT_OFF_PLOT:],yerr=FF_Dict["%s_FF_Errors"%(SYS)][ipt][ZT_OFF_PLOT:],linewidth=1, fmt='o',capsize=1,label=r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1]))
                #plt.errorbar(zT_centers, pPb_FF[ipt],xerr=zT_widths,yerr=pPb_FF_Errors[ipt],linewidth=1, fmt='o',capsize=1,label=r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[ipt],pTbins[ipt+1 
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
    
    for izt in range (0,NzT):
        weight_sum = 0
        
        for ipt in range(0,N_pT_Bins):
                    
            if (description_string == "05zT_2bins"):
                if (ipt <2 and izt <1): continue
            if (description_string == "05zT_3bins"):
                if (ipt <1 and izt <1): continue
            
            if (math.isnan(Rel_Stat_Erorr[ipt][izt]) or math.isnan(Rel_Purity_Error[ipt][izt]) or math.isnan(FF[ipt][izt])): continue
            
            Combined[izt] += (1/(FF_Errors[ipt][izt]**2)) * FF[ipt][izt]
            weight_sum += 1/(FF_Errors[ipt][izt]**2)
            #Combined[izt] += 1/((Rel_Stat_Erorr[ipt][izt]**2) + (Rel_Purity_Error[ipt][izt]**2)) * FF[ipt][izt] 
            #weight_sum += 1/((Rel_Stat_Erorr[ipt][izt]**2) + (Rel_Purity_Error[ipt][izt]**2))            

    
        Combined[izt] = Combined[izt]/weight_sum

    for izt in range(0,NzT):
        for ipt in range (N_pT_Bins):
            
            if (math.isnan(FF_Errors[ipt][izt]) or math.isnan(purity_FF_Errors[ipt][izt])): continue
                
            #Combined_Errors[izt] += FF_Errors[ipt][izt]**2
            #purity_Combined_Errors[izt] += purity_FF_Errors[ipt][izt]**2
            Combined_Errors[izt] += 1/FF_Errors[ipt][izt]**2
            purity_Combined_Errors[izt] += 1/purity_FF_Errors[ipt][izt]**2

    #Combined_Errors = np.sqrt(Combined_Errors)/N_pT_Bins
    #purity_Combined_Errors = np.sqrt(purity_Combined_Errors)/N_pT_Bins
    Combined_Errors = np.sqrt(1/Combined_Errors)
    purity_Combined_Errors = np.sqrt(1/purity_Combined_Errors)
    
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
        
    Comb_Dict = dict(zip(Keys,Combined_FF))
    
    for SYS in Systems:
                
        np.save("npy_files/%s_%s_Averaged_Fragmentation_Functions_%s.npy"%(Shower,SYS,description_string),Comb_Dict["%s_Combined_FF"%(SYS)])
        np.save("npy_files/%s_%s_Averaged_Fragmentation_Functions_Errors_%s.npy"%(Shower,SYS,description_string),Comb_Dict["%s_Combined_FF_Errors"%(SYS)])
        #print("saved to npy_files/%s_%s_Averaged_Fragmentation_Functions_%s.npy"%(Shower,SYS,description_string))
        
        Efficiency_Uncertainty = 0.05*Comb_Dict["%s_Combined_FF"%(SYS)]
        Sys_Uncertainty = np.sqrt(Efficiency_Uncertainty**2 + Comb_Dict["%s_purity_Uncertainty"%(SYS)]**2)
        np.save("npy_files/%s_%s_Averaged_Fragmentation_Functions_Systematics_%s.npy"%(Shower,SYS,description_string),Sys_Uncertainty)
        
    return Comb_Dict
        


def singleparameterPowerlawFunction(p, xmin, xmax, norm):
    def function(x):
        A = scipy.integrate.quad(lambda x: np.power(x, -p), xmin, xmax)[0]
        return norm * np.power(x, -p) / A
    return function

# hist must be normalized within binCenters!
def getSingleparameterPowerlawParamsAndErrors(hist, histerr, binCenters, xmin, xmax,norm):
    def Chi2(p):
        model = map(singleparameterPowerlawFunction(p, xmin, xmax, norm), binCenters)

        return np.sum(np.power(np.divide(np.subtract(hist, model), histerr, where=histerr != 0), 2.0))

    mt = iminuit.Minuit(Chi2, p=4.0, error_p=0.4, errordef=1, print_level=0)
    mt.migrad()

    if not mt.migrad_ok():
        print 'Warning: single-parameter power law fit did not converge'

    fitParams = {}
    fitParams['p'] = mt.values['p']
    fitParams['xmin'] = xmin
    fitParams['xmax'] = xmax
    fitParams['norm'] = norm
    fitErrors = mt.errors
    chi2dof = Chi2(**mt.values) / (len(hist) - 1)

    return fitParams, fitErrors, chi2dof, mt.migrad_ok()
    
    
def Fit_FF_PowerLaw(FF_Dictionary,SYS):
    
    print("%s Integrating %s"%(description_string,Phi_String))

    #fig = plt.figure(figsize=(12,6))
    
    #for i,SYS in enumerate(Systems):
        
    #ax = fig.add_subplot(1,2,(i+1))
    hist = FF_Dictionary["%s_Combined_FF"%(SYS)][:NzT-ZT_OFF_PLOT]
    histerr = FF_Dictionary["%s_Combined_FF_Errors"%(SYS)][:NzT-ZT_OFF_PLOT]
    binwidths = zT_widths[:NzT-ZT_OFF_PLOT]
    bincenters = zT_centers[:NzT-ZT_OFF_PLOT]
    histnorm = sum(np.multiply(hist, binwidths*2))
    Params, fiterrors, chi2dof, fitisok = getSingleparameterPowerlawParamsAndErrors(hist, histerr, bincenters, bincenters[0] - binwidths[0], bincenters[-1] - binwidths[-1], histnorm)
        
    power = Params["p"]
    print(fiterrors)

    if (fitisok == False):
        print("WARNING: POWER LAW FIT DID NOT CONVERGE")

    print(r"%s: p = %1.2f, chi2/dof = %1.2f"%(SYS,power,chi2dof))
    #print("fit Error"),
    #print(fiterrors)
        
    model = map(singleparameterPowerlawFunction(**Params), zT_centers[:NzT-ZT_OFF_PLOT])
    return model,power,chi2dof
        
        #plt.plot(zT_centers[:NzT-ZT_OFF_PLOT], model, 'g:')
        #plt.yscale("log")
        #plt.plot(zT_centers[:NzT-ZT_OFF_PLOT],FF_Dictionary["%s_Combined_FF"%(SYS)][:NzT-ZT_OFF_PLOT])

    #pp_sys_Error = (FF_Dictionary["pp_Combined_FF"][:NzT-ZT_OFF_PLOT])*math.sqrt(0.15**2+0.05**2)
    #p_Pb_sys_Error = (FF_Dictionary["p-Pb_Combined_FF"][:NzT-ZT_OFF_PLOT])*math.sqrt(0.15**2+0.05**2)
    #Chi2,NDF,Pval = Get_pp_pPb_List_Chi2(FF_Dictionary["pp_Combined_FF"][:NzT-ZT_OFF_PLOT],
    #                                     FF_Dictionary["pp_Combined_FF_Errors"][:NzT-ZT_OFF_PLOT],
    #                                     pp_sys_Error,
    #                                     FF_Dictionary["p-Pb_Combined_FF"][:NzT-ZT_OFF_PLOT],
    #                                     FF_Dictionary["p-Pb_Combined_FF_Errors"][:NzT-ZT_OFF_PLOT],
    #                                     p_Pb_sys_Error)
    
    print("Chi/NDF = %1.2f, pvalue = %1.2f"%(Chi2/NDF,Pval))

def Plot_pp_pPb_Avg_FF(Comb_Dict):
    
    Colors = ["red","blue","black"]
    fig = plt.figure(figsize=(8,8))
    
    fig.add_axes((0.1,0.3,0.88,0.6))
    for SYS,sys_col in zip(Systems,Colors):

        #Systematics
        Efficiency_Uncertainty = 0.05*Comb_Dict["%s_Combined_FF"%(SYS)]
        Sys_Uncertainty = np.sqrt(Efficiency_Uncertainty**2 + Comb_Dict["%s_purity_Uncertainty"%(SYS)]**2)
        
        #Plots
        plt.errorbar(zT_centers[:NzT-ZT_OFF_PLOT], Comb_Dict["%s_Combined_FF"%(SYS)][:NzT-ZT_OFF_PLOT],xerr=zT_widths[:NzT-ZT_OFF_PLOT],
            yerr=Comb_Dict["%s_Combined_FF_Errors"%(SYS)][:NzT-ZT_OFF_PLOT],linewidth=1, fmt='o',color=sys_col,capsize=1,
            label=r' %s %1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(SYS,pTbins[0],pTbins[N_pT_Bins]))

        Sys_Plot_pp = plt.bar(zT_centers[:NzT-ZT_OFF_PLOT], Sys_Uncertainty[:NzT-ZT_OFF_PLOT]+Sys_Uncertainty[:NzT-ZT_OFF_PLOT], 
            bottom=Comb_Dict["%s_Combined_FF"%(SYS)][:NzT-ZT_OFF_PLOT]-Sys_Uncertainty[:NzT-ZT_OFF_PLOT],width=zt_box[:NzT-ZT_OFF_PLOT], align='center',color='white',edgecolor="black",label="Systematic Uncertainty")

    plt.yscale('log')                                                                                                                                                                                                                                                              
    pplt.ylabel(r"$\frac{1}{N_{\mathrm{\gamma}}}\frac{\mathrm{d}N}{\mathrm{d}z_{\mathrm{T}}\mathrm{d}\Delta\phi\mathrm{d}\Delta\eta}$",fontsize=24)
    plt.xlabel("${z_\mathrm{T}} = p_\mathrm{T}^\mathrm{h}/p_\mathrm{T}^\mathrm{\gamma}$",fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

        
    #Chi2 and Labels
    #purity_normalization = np.full(len(Comb_Dict["pp_Combined_FF"][:NzT-ZT_OFF_PLOT]),math.sqrt(0.25**2+0.05**2)) #estimate of constant normalization Error
    pp_sys_Error = (Comb_Dict["pp_Combined_FF"][:NzT-ZT_OFF_PLOT])*math.sqrt(0.15**2+0.05**2)
    p_Pb_sys_Error = (Comb_Dict["p-Pb_Combined_FF"][:NzT-ZT_OFF_PLOT])*math.sqrt(0.15**2+0.05**2)
    #pp_sys_Error = purity_normalization
    #p-Pb_sys_Error = purity_normalization
    
    Chi2,NDF,Pval = Get_pp_pPb_List_Chi2(Comb_Dict["pp_Combined_FF"][:NzT-ZT_OFF_PLOT],
                                         Comb_Dict["pp_Combined_FF_Errors"][:NzT-ZT_OFF_PLOT],
                                         pp_sys_Error,
                                         Comb_Dict["p-Pb_Combined_FF"][:NzT-ZT_OFF_PLOT],
                                         Comb_Dict["p-Pb_Combined_FF_Errors"][:NzT-ZT_OFF_PLOT],
                                         p_Pb_sys_Error)
    #Labels
    plt.annotate("$\chi^2$ = %1.1f, ndf = %i, p = %f"%(Chi2,NDF,Pval), xy=(0.99, 0.06), xycoords='axes fraction', ha='right', va='top', fontsize=16)
        
    leg = plt.legend(numpoints=1,frameon=True,edgecolor='white', framealpha=0.0, fontsize=16)
    leg.set_title("ALICE Work in Progress\n  $\sqrt{s_{\mathrm{_{NN}}}} = $ 5 TeV")
    plt.setp(leg.get_title(),fontsize=20)

    plt.title(r'Integrated $\mathrm{\gamma}$-Hadron Correlation: $\pi/2< \Delta\varphi < \pi, |\Delta\eta| < %1.1f$ '%(eta_max),fontdict = {'fontsize' : 19})
    plt.gcf()
    plt.savefig("pics/%s/Averaged_pT_FFunction_%s.pdf"%(Shower,Shower), bbox='tight')
    plt.show()
    
    
    

    
    
    Printing = True
    
    if (Printing):
        #Printing
        for SYS in Systems:
            print("                    %s Central Values:"%(SYS))
            print(Comb_Dict["%s_Combined_FF"%(SYS)][:NzT-ZT_OFF_PLOT])
            print("")
            print("                    %s Stat. Uncertainty:"%(SYS))
            print(Comb_Dict["%s_Combined_FF_Errors"%(SYS)][:NzT-ZT_OFF_PLOT])
            print("")

            Efficiency_Uncertainty = 0.05*Comb_Dict["%s_Combined_FF"%(SYS)][:NzT-ZT_OFF_PLOT]
            Sys_Uncertainty = np.sqrt(Efficiency_Uncertainty**2 + Comb_Dict["%s_purity_Uncertainty"%(SYS)][:NzT-ZT_OFF_PLOT]**2)

            print("              %s Systematic Uncertainty:"%(SYS))
            print(Sys_Uncertainty)
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
        
    Comb_Dict = dict(zip(Keys,Combined_FF))
    
    for SYS in Systems:
                
        np.save("npy_files/%s_%s_Averaged_Fragmentation_Functions_%s.npy"%(Shower,SYS,description_string),Comb_Dict["%s_Combined_FF"%(SYS)])
        np.save("npy_files/%s_%s_Averaged_Fragmentation_Functions_Errors_%s.npy"%(Shower,SYS,description_string),Comb_Dict["%s_Combined_FF_Errors"%(SYS)])
        #print("saved to npy_files/%s_%s_Averaged_Fragmentation_Functions_%s.npy"%(Shower,SYS,description_string))
        
        Efficiency_Uncertainty = 0.05*Comb_Dict["%s_Combined_FF"%(SYS)]
        Sys_Uncertainty = np.sqrt(Efficiency_Uncertainty**2 + Comb_Dict["%s_purity_Uncertainty"%(SYS)]**2)
        np.save("npy_files/%s_%s_Averaged_Fragmentation_Functions_Systematics_%s.npy"%(Shower,SYS,description_string),Sys_Uncertainty)
        
    return Comb_Dict
        

def Plot_pp_pPb_Avg_FF_and_Ratio(Comb_Dict):
    
    plot_power = False
    Colors = ["red","blue","black"]
    fig = plt.figure(figsize=(8,8))
    
    fig.add_axes((0.1,0.3,0.88,0.6))
    for SYS,sys_col in zip(Systems,Colors):

        #Systematics
        Efficiency_Uncertainty = 0.05*Comb_Dict["%s_Combined_FF"%(SYS)]
        Sys_Uncertainty = np.sqrt(Efficiency_Uncertainty**2 + Comb_Dict["%s_purity_Uncertainty"%(SYS)]**2)
        
        #Plots
        plt.errorbar(zT_centers[:NzT-ZT_OFF_PLOT], Comb_Dict["%s_Combined_FF"%(SYS)][:NzT-ZT_OFF_PLOT],xerr=zT_widths[:NzT-ZT_OFF_PLOT],
            yerr=Comb_Dict["%s_Combined_FF_Errors"%(SYS)][:NzT-ZT_OFF_PLOT],linewidth=1, fmt='o',color=sys_col,capsize=1,
            label=SYS)
            #label=r' %s %1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(SYS,pTbins[0],pTbins[N_pT_Bins]))
            
        Sys_Plot_pp = plt.bar(zT_centers[:NzT-ZT_OFF_PLOT], Sys_Uncertainty[:NzT-ZT_OFF_PLOT]+Sys_Uncertainty[:NzT-ZT_OFF_PLOT],
            bottom=Comb_Dict["%s_Combined_FF"%(SYS)][:NzT-ZT_OFF_PLOT]-Sys_Uncertainty[:NzT-ZT_OFF_PLOT],width=zT_widths[:NzT-ZT_OFF_PLOT], align='center',color=sys_col,alpha=0.3)
            #bottom=Comb_Dict["%s_Combined_FF"%(SYS)][:NzT-ZT_OFF_PLOT]-Sys_Uncertainty[:NzT-ZT_OFF_PLOT],width=zt_box[:NzT-ZT_OFF_PLOT], align='center',color=sys_col,alpha=0.3)
        
        
        if (plot_power):
            model,p,chi2dof = Fit_FF_PowerLaw(Comb_Dict,SYS)
            plt.plot(zT_centers[:NzT-ZT_OFF_PLOT], model, sys_col,label=r"%s $\alpha = %1.2f\pm 0.1 \chi^2 = %1.2f$"%(SYS,p,chi2dof))
    
    plt.errorbar(zT_centers[:NzT-ZT_OFF_PLOT],pythia_FF,xerr=zT_widths[:NzT-ZT_OFF_PLOT],fmt='-g',label="Pythia 8.2 Monash")   
    
    plt.yscale('log')                             
    plt.ylabel(r"$\frac{1}{N_{\mathrm{\gamma}}}\frac{\mathrm{d}N}{\mathrm{d}z_{\mathrm{T}}\mathrm{d}\Delta\phi\mathrm{d}\Delta\eta}$",fontsize=24)
    plt.ylim(0.037,15)
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=0)
    plt.xlim(0,0.65)

        
    #Chi2 and Labels
    #purity_normalization = np.full(len(Comb_Dict["pp_Combined_FF"][:NzT-ZT_OFF_PLOT]),math.sqrt(0.15**2+0.05**2)) #estimate of constant normalization Error
    pp_sys_Error = (Comb_Dict["pp_Combined_FF"][:NzT-ZT_OFF_PLOT])*math.sqrt(0.15**2+0.05**2)
    p_Pb_sys_Error = (Comb_Dict["p-Pb_Combined_FF"][:NzT-ZT_OFF_PLOT])*math.sqrt(0.15**2+0.05**2)
    #pp_sys_Error = purity_normalization
    #p-Pb_sys_Error = purity_normalization
    
    Chi2,NDF,Pval = Get_pp_pPb_List_Chi2(Comb_Dict["pp_Combined_FF"][:NzT-ZT_OFF_PLOT],
                                         Comb_Dict["pp_Combined_FF_Errors"][:NzT-ZT_OFF_PLOT],
                                         pp_sys_Error,
                                         Comb_Dict["p-Pb_Combined_FF"][:NzT-ZT_OFF_PLOT],
                                         Comb_Dict["p-Pb_Combined_FF_Errors"][:NzT-ZT_OFF_PLOT],
                                         p_Pb_sys_Error)

    plt.annotate("$\chi^2$ = %1.1f, ndf = %i, p = %1.2f"%(Chi2,NDF,Pval), xy=(0.01, 0.06), xycoords='axes fraction', ha='left', va='top', fontsize=16)
    #plt.annotate("$\chi^2$ = %1.1f, ndf = %i, p = %f"%(Chi2,NDF,Pval), xy=(0.99, 0.06), xycoords='axes fraction', ha='right', va='top', fontsize=16)
    
    #plt.annotate("%s"%(description_string),xy=(0.01,0.1),xycoords="axes fraction",ha="left",va="top",fontsize=12)
    
    leg = plt.legend(numpoints=1,frameon=True,edgecolor='white', framealpha=0.0, fontsize=16)
    leg.set_title("ALICE Work in Progress\n  $\sqrt{s_{\mathrm{_{NN}}}} = $ 5 TeV \n")
    plt.setp(leg.get_title(),fontsize=20)
    plt.annotate("%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$"%(pTbins[0],pTbins[N_pT_Bins]),xy=(0.53, 0.81), xycoords='axes fraction', ha='left', va='top', fontsize=16)
    
    crap_boxes = False
    if (crap_boxes):
        plt.text(0.441, 2.92, '__',
         {'color': 'black', 'alpha': 0.0, 'fontsize': 16, 'ha': 'left', 'va': 'top',
          'bbox': dict(boxstyle="square", fc='red',alpha=0.3, ec="None", pad=0.2)})
        plt.text(0.441, 1.95, '__',
         {'color': 'black','alpha':0.0, 'fontsize': 16, 'ha': 'left', 'va': 'top',
          'bbox': dict(boxstyle="square", fc="blue",alpha=0.3, ec="None", pad=0.2)})

    plt.title(r'Integrated $\mathrm{\gamma}$-Hadron Correlation: $%s < \Delta\varphi < \pi$ '%(Phi_String),fontdict = {'fontsize' : 19})
    
    fig.add_axes((0.1,0.1,0.88,0.2))

    pPb_Combined = Comb_Dict["p-Pb_Combined_FF"]
    pPb_Combined_Errors = Comb_Dict["p-Pb_Combined_FF_Errors"]
    pPb_purity_Uncertainty = Comb_Dict["p-Pb_purity_Uncertainty"]
    
    pp_Combined = Comb_Dict["pp_Combined_FF"]
    pp_Combined_Errors = Comb_Dict["pp_Combined_FF_Errors"]
    pp_purity_Uncertainty = Comb_Dict["pp_purity_Uncertainty"]
    
    Ratio = pPb_Combined/pp_Combined
    Ratio_Error = np.sqrt((pPb_Combined_Errors/pPb_Combined)**2 + (pp_Combined_Errors/pp_Combined)**2)*Ratio
    Ratio_Plot = plt.errorbar(zT_centers[:NzT-ZT_OFF_PLOT], Ratio[:NzT-ZT_OFF_PLOT], yerr=Ratio_Error[:NzT-ZT_OFF_PLOT],xerr=zT_widths[:NzT-ZT_OFF_PLOT], fmt='ko',capsize=3, ms=6,lw=1)
    
        #Save
    np.save("npy_files/%s_Averaged_FF_Ratio_%s.npy"%(Shower,description_string),Ratio)
    np.save("npy_files/%s_Averaged_FF_Ratio_Errors_%s.npy"%(Shower,description_string),Ratio_Error)
    
    Purity_Uncertainty = np.sqrt((pp_purity_Uncertainty/pp_Combined)**2 + (pPb_purity_Uncertainty/pPb_Combined)**2)*Ratio
    Efficiency_Uncertainty = np.ones(len(pPb_Combined))*0.05*math.sqrt(2)*Ratio 
    if (CorrectedP):
        Ratio_Systematic = np.sqrt(Purity_Uncertainty**2 + Efficiency_Uncertainty**2)
    
    Sys_Plot = plt.bar(zT_centers[:NzT-ZT_OFF_PLOT], Ratio_Systematic[:NzT-ZT_OFF_PLOT]+Ratio_Systematic[:NzT-ZT_OFF_PLOT],
            #bottom=Ratio[:NzT-ZT_OFF_PLOT]-Ratio_Systematic[:NzT-ZT_OFF_PLOT], width=zt_box[:NzT-ZT_OFF_PLOT], align='center',edgecolor="k",color='w')
            bottom=Ratio[:NzT-ZT_OFF_PLOT]-Ratio_Systematic[:NzT-ZT_OFF_PLOT], width=zT_widths[:NzT-ZT_OFF_PLOT], align='center',color='black',alpha=0.2)
    
    #Sys_Plot = plt.bar(zT_centers[:NzT-ZT_OFF_PLOT], 2*Ratio_Systematic[:NzT-ZT_OFF_PLOT], 
    #   bottom=1.0-Ratio_Systematic[:NzT-ZT_OFF_PLOT], width=2*zT_widths[:NzT-ZT_OFF_PLOT], align='center',color='black',alpha = 0.2)
    
    
    plt.axhline(y=1, color='k', linestyle='--')
    
    plt.xlabel("${z_\mathrm{T}} = p_\mathrm{T}^{\mathrm{h}}/p_\mathrm{T}^\gamma$",fontsize=20)
    plt.ylabel(r"$\frac{\mathrm{p-Pb}}{\mathrm{pp}}$",fontsize=24)
    #plt.ylim((-0.5, 2.5))
    plt.xlabel("${z_\mathrm{T}} = p_\mathrm{T}^\mathrm{h}/p_\mathrm{T}^\mathrm{\gamma}$",fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks([0.5,1.0,1.5,2.0,2.5],fontsize=16)
    plt.xlim(0,0.65)

    #model_p-Pb = Fit_FF_PowerLaw(Comb_Dict,"p-Pb")
    #plt.plot(zT_centers[:NzT-ZT_OFF_PLOT], model_pp, 'g:')
        #plt.yscale("log")
    
    plt.gcf()
    #plt.tight_layout()
    plt.savefig("pics/%s/%s/Final_FFunction_and_Ratio.pdf"%(Shower,description_string), bbox_inches = "tight")
    plt.show()
        
def pp_pPB_Avg_Ratio(Comb_Dict,pT_Start):
    
    pPb_Combined = Comb_Dict["p-Pb_Combined_FF"]
    pPb_Combined_Errors = Comb_Dict["p-Pb_Combined_FF_Errors"]
    pPb_purity_Uncertainty = Comb_Dict["p-Pb_purity_Uncertainty"]
    
    pp_Combined = Comb_Dict["pp_Combined_FF"]
    pp_Combined_Errors = Comb_Dict["pp_Combined_FF_Errors"]
    pp_purity_Uncertainty = Comb_Dict["pp_purity_Uncertainty"]
    
    #Ratio
    Ratio = pPb_Combined/pp_Combined
    
    #Stat. Error in ratio
    Ratio_Error = np.sqrt((pPb_Combined_Errors/pPb_Combined)**2 + (pp_Combined_Errors/pp_Combined)**2)*Ratio
    
    #Save
    np.save("npy_files/%s_Averaged_FF_Ratio_%s.npy"%(Shower,description_string),Ratio)
    np.save("npy_files/%s_Averaged_FF_Ratio_Errors_%s.npy"%(Shower,description_string),Ratio_Error)
    
    #Sys. Error in ratio
    Purity_Uncertainty = np.sqrt((pp_purity_Uncertainty/pPb_Combined)**2 + (pPb_purity_Uncertainty/pPb_Combined)**2)*Ratio
    Efficiency_Uncertainty = np.ones(len(pPb_Combined))*0.05*math.sqrt(2)*Ratio 
    if (CorrectedP):
        Ratio_Systematic = np.sqrt(Purity_Uncertainty**2 + Efficiency_Uncertainty**2)

    plt.figure(figsize=(10,7)) 


    #Sys_Plot = plt.bar(zT_centers[:NzT-ZT_OFF_PLOT], Ratio_Systematic[:NzT-ZT_OFF_PLOT]+Ratio_Systematic[:NzT-ZT_OFF_PLOT],
    #        bottom=Ratio[:NzT-ZT_OFF_PLOT]-Ratio_Systematic[:NzT-ZT_OFF_PLOT], width=zt_box[:NzT-ZT_OFF_PLOT], align='center',edgecolor="k",color='w')
            #bottom=Ratio[:NzT-ZT_OFF_PLOT]-Ratio_Systematic[:NzT-ZT_OFF_PLOT], width=zt_box[:NzT-ZT_OFF_PLOT], align='center',color='black',alpha=0.2)

    Sys_Plot = plt.bar(zT_centers[:NzT-ZT_OFF_PLOT], 2*Ratio_Systematic[:NzT-ZT_OFF_PLOT], 
       bottom=1.0-Ratio_Systematic[:NzT-ZT_OFF_PLOT], width=2*zT_widths[:NzT-ZT_OFF_PLOT], align='center',color='black',alpha = 0.2)
    
    empt4, = plt.plot([], [],' ')

    Ratio_Plot = plt.errorbar(zT_centers[:NzT-ZT_OFF_PLOT], Ratio[:NzT-ZT_OFF_PLOT], yerr=Ratio_Error[:NzT-ZT_OFF_PLOT],xerr=zT_widths[:NzT-ZT_OFF_PLOT], fmt='ko',capsize=3, ms=6,lw=1)

    plt.xlabel("${z_\mathrm{T}} = p_\mathrm{T}^{\mathrm{h}}/p_\mathrm{T}^\gamma$",fontsize=20)
    plt.ylabel(r"$\frac{\mathrm{p-Pb}}{\mathrm{pp}}$",fontsize=20)
    plt.ylim((-0.49, 2.9))
    #plt.yticks(np.arange(-0, 2, step=0.2))
    
    if(NzT == 6):
        plt.xlim(xmin = 0.0,xmax=0.7)
    elif(NzT==7):
        plt.xlim(xmin = 0.0,xmax=1.0)
    #plt.xlim(xmin = 0.0,xmax=zTbins[NzT-ZT_OFF_PLOT])
    plt.xlim(xmin = 0.0,xmax=0.67)
    plt.axhline(y=1, color='k', linestyle='--')

    ### ROOT LINEAR and CONSTANT FITS ###
    Ratio_TGraph = TGraphErrors()
    for izt in range (len(Ratio)-ZT_OFF_PLOT):
        Ratio_TGraph.SetPoint(izt,zT_centers[izt],Ratio[izt])
        Ratio_TGraph.SetPointError(izt,0,Ratio_Error[izt])

    Ratio_TGraph.Fit("pol0","S")
    f = Ratio_TGraph.GetFunction("pol0")
    chi2_red  = f.GetChisquare()/f.GetNDF()
    pval = f.GetProb()
    p0 = f.GetParameter(0)
    p0e = f.GetParError(0)
    p0col = "blue"
    if (Show_Fits):
        plt.annotate("Constant Fit", xy=(0.05, 0.99), xycoords='axes fraction', ha='left', va='top', color=p0col,fontsize=20,alpha=.7)
        plt.annotate(r"$p0 = {0:.2f} \pm {1:.2f}$".format(p0,p0e), xy=(0.05, 0.94), xycoords='axes fraction', ha='left', va='top', color=p0col,fontsize=18,alpha=.7)
        plt.annotate(r"$\chi^2_{red} = %1.2f$"%(chi2_red), xy=(0.05, 0.89), xycoords='axes fraction', ha='left', va='top', color=p0col,fontsize=18,alpha=.7)
        plt.annotate(r"$p_{val} = %1.2f$"%(pval), xy=(0.05, 0.84), xycoords='axes fraction', ha='left', va='top', color=p0col,fontsize=18,alpha=.7)

        plt.fill_between(np.arange(0,1.1,0.1), p0+p0e, p0-p0e,color=p0col,alpha=.2)

    
    
    Ratio_TGraph.Fit("pol1","S")
    #zT_Points = np.linspace(0.05,zTbins[NzT-1],20)
    zT_Points = np.linspace(0.0,1,20)
    
    Fit_Band = ROOT.TGraphErrors(len(zT_Points));
    for i in range(len(zT_Points)-ZT_OFF_PLOT):
        Fit_Band.SetPoint(i, zT_Points[i], 0)
    (ROOT.TVirtualFitter.GetFitter()).GetConfidenceIntervals(Fit_Band,0.68)
    
    band_errors = np.zeros(len(zT_Points))
    for i in range (len(zT_Points)):
        band_errors[i] = Fit_Band.GetErrorY(i)
    print(band_errors)


    f2 = Ratio_TGraph.GetFunction("pol1")
    chi2_red  = f2.GetChisquare()/f2.GetNDF()
    pval = f2.GetProb()
    p0 = f2.GetParameter(0)
    p0e = f2.GetParError(0)
    p1 = f2.GetParameter(1)
    p1e = f2.GetParError(1)
    print(p1e)
    p1col = "Green"
    if (Show_Fits):
        
        plt.annotate("Linear Fit", xy=(0.4, 0.99), xycoords='axes fraction', ha='left', va='top', color=p1col,fontsize=20,alpha=.7)
        plt.annotate(r"$p0 = %1.2f \pm %1.2f$"%(p0,p0e), xy=(0.4, 0.94), xycoords='axes fraction', ha='left', va='top', color=p1col,fontsize=18,alpha=.7)
        plt.annotate(r"$p1 = %1.2f \pm %1.2f$"%(p1,p1e), xy=(0.4, 0.89), xycoords='axes fraction', ha='left', va='top', color=p1col,fontsize=18,alpha=.7)
        plt.annotate(r"$\chi^2_{red} = %1.2f$"%(chi2_red), xy=(0.4, 0.84), xycoords='axes fraction', ha='left', va='top', color=p1col,fontsize=18,alpha=.7)
        plt.annotate(r"$p_{val} = %1.2f$"%(pval), xy=(0.4, 0.79), xycoords='axes fraction', ha='left', va='top', color=p1col,fontsize=18,alpha=.7)
        
        axes = plt.gca()
        x_vals = np.array(zT_Points)
        y_vals = p0 + p1 * zT_Points
        plt.plot(zT_Points, y_vals, '--',color=p1col,linewidth=2,alpha=0.5)
        plt.fill_between(zT_Points,y_vals+band_errors,y_vals-band_errors,color=p1col,alpha=0.2)

    ### ROOT DONE ###


    leg = plt.legend([Ratio_Plot,empt4],["Statistical Error",r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[pT_Start],pTbins[N_pT_Bins])],frameon=False,numpoints=1,loc="lower left",title=' ',prop={'size':18})

    leg.set_title("ALICE Work in Progress\n  $\sqrt{s_{\mathrm{_{NN}}}} = $ 5 TeV")
    plt.setp(leg.get_title(),fontsize=20)

    plt.gcf()
    plt.savefig("pics/%s/%s/Ratio_Fits.pdf"%(Shower,description_string), bbox='tight')
    plt.show()

    print("                Central Values:")
    print(Ratio[:NzT-ZT_OFF_PLOT])

    print("\n                Satistical Uncertainty Absolute:")
    print(Ratio_Error[:NzT-ZT_OFF_PLOT])
    
    print("\n               Relative Satistical Uncertainty:")
    print(Ratio_Error[:NzT-ZT_OFF_PLOT]/Ratio[:NzT-ZT_OFF_PLOT])

    print("\n                Ratio Uncertainty from Purity:")
    print(Purity_Uncertainty[:NzT-ZT_OFF_PLOT])

    print("\n                Ratio Uncertainty from Single Track Efficiency:")
    print(Efficiency_Uncertainty[:NzT-ZT_OFF_PLOT])

    print("\n                Full Systematic Uncertainty:")

    print(Ratio_Systematic[:NzT-ZT_OFF_PLOT])
    

    
def Compare_pp_pPB_Avg_Ratio_lists(strings,string_descrp_list,colors,Show_Fits = False,Avg_Cs = False):
    
    plt.figure(figsize=(8,8))
    
    for (string,string_descr,colr) in zip(strings,string_descrp_list,colors):
        
        if not(Avg_Cs):
            Ratio = np.load("npy_files/LO_Averaged_FF_Ratio_%s.npy"%(string))
            Ratio_Error = np.load("npy_files/LO_Averaged_FF_Ratio_Errors_%s.npy"%(string))
        
        else:
            Ratio = np.load("npy_files/LO_Cs_Avg_FF_Ratio_%s.npy"%(string))
            Ratio_Error = np.load("npy_files/LO_Cs_Avg_FF_Ratio_Errors_%s.npy"%(string))

        Zbins = np.geomspace(0.06, 0.6, num=len(Ratio)+1)
        zT_centers = (Zbins[1:] + Zbins[:-1]) / 2
        zT_widths = [(j-i)/2 for i, j in zip(Zbins[:-1], Zbins[1:])]
        
        
        ### ROOT LINEAR and CONSTANT FITS ###
        Ratio_TGraph = TGraphErrors()
        for izt in range (ZT_OFF_PLOT,len(Ratio)):
        #for izt in range (2,len(Ratio)-1):
            Ratio_TGraph.SetPoint(izt,zT_centers[izt],Ratio[izt])
            Ratio_TGraph.SetPointError(izt,0,Ratio_Error[izt])

        Ratio_TGraph.Fit("pol0","S")
        f = Ratio_TGraph.GetFunction("pol0")
        chi2_red  = f.GetChisquare()/f.GetNDF()
        pval = f.GetProb()
        p0 = f.GetParameter(0)
        p0e = f.GetParError(0)
        p0col = colr
        if (Show_Fits):
            plt.fill_between(np.arange(0,1.1,0.1), p0+p0e, p0-p0e,color=p0col,alpha=.3)
        
        plt.errorbar(zT_centers[:NzT-ZT_OFF_PLOT], Ratio[:NzT-ZT_OFF_PLOT], yerr=Ratio_Error[:NzT-ZT_OFF_PLOT],xerr=zT_widths[:NzT-ZT_OFF_PLOT],capsize=3, fmt ="o",color=colr,alpha=0.7,ms=6,lw=1,label=string_descr)

    empt4, = plt.plot([], [],' ',label=r'%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$'%(pTbins[0],pTbins[N_pT_Bins]))
    
    plt.xlabel("${z_\mathrm{T}} = p_\mathrm{T}^{\mathrm{h}}/p_\mathrm{T}^\gamma$",fontsize=20)
    plt.ylabel(r"$\frac{\mathrm{p-Pb}}{\mathrm{pp}}$",fontsize=20)
    plt.ylim((-0, 3))
    #plt.yticks(np.arange(-0, 2, step=0.2))
    
    #if(NzT == 6):
    #    plt.xlim(xmin = 0.0,xmax=0.7)
    #elif(NzT==7):
    plt.xlim(xmin = 0.0,xmax=1.0)
    plt.xlim(xmin = 0.0,xmax = zTbins[NzT-ZT_OFF_PLOT])
    plt.axhline(y=1, color='r', linestyle='--')

    leg = plt.legend(frameon=False,numpoints=1,loc="best",title=' ',prop={'size':18})
    leg.set_title("ALICE Work in Progress\n  $\sqrt{s_{\mathrm{_{NN}}}} = $ 5 TeV")
    plt.setp(leg.get_title(),fontsize=20)

    plt.gcf()
    #plt.savefig("pics/Comprison_Averaged_pT_FFunction_ratio_%s_%s_and_%s.pdf"%(Shower,string_descr,string_descr2), bbox_inches='tight')
    plt.show()

    print("                Central Values:")
    print(Ratio[ZT_OFF_PLOT:])
    
    
def Compare_pp_pPB_Avg_lists(strings,string_descrp_list,colors):
        
    plt.figure(figsize=(8,8)) 
    shapes = ["o","x","s"]
    
    for (string,string_descr,colr) in zip(strings,string_descrp_list,colors):
        print(string)
        for SYS,shape in zip(Systems,shapes):  
                
            FF = np.load("npy_files/%s_%s_Averaged_Fragmentation_Functions_%s.npy"%(Shower,SYS,string))
            FF_Errors = np.load("npy_files/%s_%s_Averaged_Fragmentation_Functions_Errors_%s.npy"%(Shower,SYS,string))
            print("loading npy_files/%s_%s_Averaged_Fragmentation_Functions_Errors_%s.npy"%(Shower,SYS,string))
        
            Zbins = np.geomspace(0.05, 1.0, num=len(FF)+1)
            zT_centers = (Zbins[1:] + Zbins[:-1]) / 2
            zT_widths = [(j-i)/2 for i, j in zip(Zbins[:-1], Zbins[1:])]
        
            plt.errorbar(zT_centers[:NzT-ZT_OFF_PLOT]+0.2, FF[:NzT-ZT_OFF_PLOT],xerr=zT_widths[:NzT-ZT_OFF_PLOT],
                yerr=FF_Errors[:NzT-ZT_OFF_PLOT],linewidth=1, fmt=shape,color=colr,capsize=1,label="%s (%s)"%(string_descr,SYS))

            plt.yscale('log')                                                                                                                                                                                                                                                              
            plt.ylabel(r"$\frac{1}{N_{\mathrm{\gamma}}}\frac{\mathrm{d}N}{\mathrm{d}z_{\mathrm{T}} \mathrm{d}\Delta\eta}$",fontsize=20)
            plt.xlabel("${z_\mathrm{T}} = p_\mathrm{T}^\mathrm{h}/p_\mathrm{T}^\mathrm{\gamma}$",fontsize=20)
            #plt.xlim(xmin = 0.1,xmax=0.7)
            #plt.ylim(ymin = 0.001,ymax=20)
            
    leg = plt.legend(numpoints=1,frameon=False)
    leg.set_title("ALICE Work in Progress\n  $\sqrt{s_{\mathrm{_{NN}}}} = $ 5 TeV")
    plt.setp(leg.get_title(),fontsize=18)

    plt.title(r'Integrated $\mathrm{\gamma}$-Hadron Correlation: $2\pi/3 < \Delta\varphi < \pi, |\Delta\eta| < %1.1f$ '%(eta_max),fontdict = {'fontsize' : 19})
    plt.gcf()
    plt.savefig("pics/%s/Averaged_pT_FFunction_%s.pdf"%(Shower,Shower), bbox_inches='tight')
    plt.show()

    print("                Central Values:")
    
    
def LaTeX_Table(FF_Dictionary):
    
    plot = True
    
    print(r"%s Intg. %s"%(description_string,Phi_String)),

    fig = plt.figure(figsize=(12,6))

    for i,SYS in enumerate(Systems):
        
        hist = FF_Dictionary["%s_Combined_FF"%(SYS)][:NzT-ZT_OFF_PLOT]
        histerr = FF_Dictionary["%s_Combined_FF_Errors"%(SYS)][:NzT-ZT_OFF_PLOT]
        binwidths = zT_widths[:NzT-ZT_OFF_PLOT]
        bincenters = zT_centers[:NzT-ZT_OFF_PLOT]
        histnorm = sum(np.multiply(hist, binwidths*2))
        Params, fiterrors, chi2dof, fitisok = getSingleparameterPowerlawParamsAndErrors(hist, histerr, bincenters, bincenters[0] - binwidths[0], bincenters[-1] - binwidths[-1], histnorm)

        p = Params["p"]
        p_error = fiterrors["p"]
        
        if (fitisok == False):
            print("WARNING: POWER LAW FIT DID NOT CONVERGE")
        if plot:
            ax = fig.add_subplot(1,2,(i+1))
            model = map(singleparameterPowerlawFunction(**Params), zT_centers[:NzT-ZT_OFF_PLOT])
        
        plt.plot(zT_centers[:NzT-ZT_OFF_PLOT], model, 'g:')
        plt.yscale("log")
        plt.plot(zT_centers[:NzT-ZT_OFF_PLOT],FF_Dictionary["%s_Combined_FF"%(SYS)][:NzT-ZT_OFF_PLOT])

        string = r" & $%1.2f \pm %1.2f$ & %1.2f"%(p,p_error,chi2dof)
        print(string),

    pp_sys_Error = (FF_Dictionary["pp_Combined_FF"][:NzT-ZT_OFF_PLOT])*math.sqrt(0.15**2+0.05**2)
    p_Pb_sys_Error = (FF_Dictionary["p-Pb_Combined_FF"][:NzT-ZT_OFF_PLOT])*math.sqrt(0.15**2+0.05**2)
    Chi2,NDF,Pval = Get_pp_pPb_List_Chi2(FF_Dictionary["pp_Combined_FF"][:NzT-ZT_OFF_PLOT],
                                             FF_Dictionary["pp_Combined_FF_Errors"][:NzT-ZT_OFF_PLOT],
                                             pp_sys_Error,
                                             FF_Dictionary["p-Pb_Combined_FF"][:NzT-ZT_OFF_PLOT],
                                             FF_Dictionary["p-Pb_Combined_FF_Errors"][:NzT-ZT_OFF_PLOT],
                                             p_Pb_sys_Error)
        
    print(" & $%1.2f/%i\ %1.2f $\\\\"%(Chi2,NDF,Pval))
    
def Compare_FF_Integration(ranges,strings):

    Corr = ROOT_to_nparray() # Make a dictionary from numpy arrays
    Correlated_Subtraction_Weights(Corr)
    Ped_Sub_After_Cs(Corr)
    
    colors = ["b","r","g","p"]
    markers = ["o","^","s","*"]
    
    for SYS in Systems:
        
        fig = plt.figure(figsize=(8,8))
        fig.add_axes((0.1,0.3,0.88,0.6))
        
        denominator_window = 2.7
        for i,dphi in enumerate(dPhi_Bins):
            if (dphi >= denominator_window):
                dphi_start = i
                break
        N_Phi_Bin_start = len(delta_phi_centers)-dphi_start
        
        Frag_denominator = Get_Fragmentation(Corr,N_Phi_Bin_start) # R=04 is the denominator
        denominator = Frag_denominator["%s_FF"%(SYS)][0]/np.sum(Frag_denominator["%s_FF"%(SYS)][0])
        denominator_error = Frag_denominator["%s_FF_Errors"%(SYS)][0]/np.sum(Frag_denominator["%s_FF"%(SYS)][0])
        
        for window,phi_string,markr in zip(ranges,strings,markers):
            
            for i,dphi in enumerate(dPhi_Bins):
                if (dphi >= window):
                    dphi_start = i
                    break
            N_Phi_Bin_start = len(delta_phi_centers)-dphi_start
            
            Frags = Get_Fragmentation(Corr,N_Phi_Bin_start)
            plt.errorbar(zT_centers,(Frags["%s_FF"%(SYS)][0]/np.sum(Frags["%s_FF"%(SYS)][0])),yerr=(Frags["%s_FF_Errors"%(SYS)][0]/np.sum(Frags["%s_FF"%(SYS)][0])),xerr=zT_widths[:NzT-ZT_OFF_PLOT],linewidth=1,fmt=".",marker=markr,capsize=1,label=r'Integral $\Delta\phi > %s$'%(phi_string))

            plt.yscale('log')   
            plt.yticks(fontsize=16)
            plt.ylabel(r"$\frac{1}{N_{\mathrm{\gamma}}}\frac{\mathrm{d}N}{\mathrm{d}z_{\mathrm{T}}\mathrm{d}\Delta\phi\mathrm{d}\Delta\eta}\ \mathrm{(normalized)}$",fontsize=24)
            plt.xlabel("${z_\mathrm{T}} = p_\mathrm{T}^\mathrm{h}/p_\mathrm{T}^\mathrm{\gamma}$",fontsize=20)
                    #plt.xlim(xmin = 0.1,xmax=0.7)
            plt.ylim(ymin = 0.001,ymax=0.6)

            #leg = plt.legend(numpoints=1,frameon=False)
            #leg.set_title("ALICE Work in Progress\n  $\sqrt{s_{\mathrm{_{NN}}}} = $ 5 TeV %s"%(SYS))
            #
            #plt.setp(leg.get_title(),fontsize=18)
            
            leg = plt.legend(numpoints=1,frameon=True,edgecolor='white', framealpha=0.0, fontsize=14)
            leg.set_title("ALICE Work in Progress\n$\sqrt{s_{\mathrm{_{NN}}}} = $ 5 TeV %s \n"%(SYS))
            plt.setp(leg.get_title(),fontsize=20)
            plt.annotate("%1.0f < $p_\mathrm{T}^{\mathrm{trig}}$ < %1.0f GeV/$c$"%(pTbins[0],pTbins[N_pT_Bins]),xy=(0.52, 0.81), xycoords='axes fraction', ha='left', va='top', fontsize=14,fontweight=10)

            plt.title(r'Integrated $\mathrm{\gamma}$-Hadron Correlation',fontdict = {'fontsize' : 19})
            
            
        fig.add_axes((0.1,0.1,0.88,0.2))
        
        for window,phi_string,colr,markr in zip(ranges,strings,colors,markers):
            
            for i,dphi in enumerate(dPhi_Bins):
                if (dphi >= window):
                    dphi_start = i
                    break
            N_Phi_Bin_start = len(delta_phi_centers)-dphi_start
            
            Frags = Get_Fragmentation(Corr,N_Phi_Bin_start)

            numerator = Frags["%s_FF"%(SYS)][0]/np.sum(Frags["%s_FF"%(SYS)][0])
            numerator_error = Frags["%s_FF_Errors"%(SYS)][0]/np.sum(Frags["%s_FF"%(SYS)][0])
            Ratio = numerator/denominator
            Ratio_error = np.sqrt((numerator_error/numerator)**2 + (denominator_error/denominator)**2)*Ratio

            plt.errorbar(zT_centers,Ratio,yerr=0,xerr=zT_widths[:NzT-ZT_OFF_PLOT],linewidth=1, fmt=".",marker=markr,capsize=1)

            plt.axhline(y=1, color='k', linestyle='--')
    
            plt.xlabel("${z_\mathrm{T}} = p_\mathrm{T}^{\mathrm{h}}/p_\mathrm{T}^\gamma$",fontsize=20)
            plt.ylabel(r"$\frac{R}{R=0.4}$",fontsize=24)
            plt.ylim((0.61,1.39))
            plt.xlabel("${z_\mathrm{T}} = p_\mathrm{T}^\mathrm{h}/p_\mathrm{T}^\mathrm{\gamma}$",fontsize=20)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=14)
            plt.xlim(0,0.65)
            plt.gcf()
            #plt.tight_layout()
            plt.savefig("pics/%s/%s/%s_Integrations_FFunction_Comparison.pdf"%(Shower,description_string,SYS),bbox_inches = "tight")

            

        #Plot_pp_pPb_Avg_FF_and_Ratio(Combined_Frags)
        #Combined_Frags = Average_FF(Frags)