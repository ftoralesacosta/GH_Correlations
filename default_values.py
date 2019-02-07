def Set_Defaults(Shower="NN", UnCorrected_P=False,No_Weights=False, Use_MC = False):
    Shower = "NN"
    #Shower = "LO"
    Use_Weights = False
    CorrectedP = True     #FALSE FOR HARDPROBES

    if (Shower == "NN"):
        if (Use_Weights):
            pPb_File = 'InputData/pPb_SE_NN_Correlation_GMB_Ratio.root'
            pp_File = 'InputData/pp_SE_NN_Correlation_GMB_Ratio.root'
        else:
            pPb_File = 'InputData/pPb_SE_NN_Correlation_GMB_Ratio_UnWeight.root'
            pp_File = 'InputData/pp_SE_NN_Correlation_GMB_Ratio_UnWeight.root'           
        
        if (CorrectedP):
            purity = [0.276899, 0.358741, 0.456807, 0.476192]
        else:
            purity = 0.352546
            
    if (Shower == "LO"):
        if (Use_Weights):
            pPb_File = 'InputData/pPb_SE_L0_Correlation_GMB_Ratio.root'
            pp_File = 'InputData/pp_SE_L0_Correlation_GMB_Ratio.root'
        else:
            pPb_File = 'InputData/pPb_SE_L0_Correlation_GMB_Ratio_UnWeight.root'
            pp_File = 'InputData/pp_SE_L0_Correlation_GMB_Ratio_UnWeight.root'
        if (CorrectedP):
            #purity = 0.277
            purity = [0.276791, 0.358627, 0.456378, 0.476337]
        else:
            purity = 0.35
            
        print purity
            
    MC_File = 'InputData/18b10a_pthat_1_2_SE_NN_Correlation_GMB_Ratio.root'

    Use_MC = False

    if(Use_MC):
        Systems = ["pp","p-Pb","MC"]
        Files = [pp_File,pPb_File,MC_File]
    
    else:
        Systems = ["pp","p-Pb"]
        Files = [pp_File,pPb_File]
        #Systems = ["pPb"]
        #Files = [pPb_File]


        #pPb_File = 'InputData/13def_EMax_SE_GMB_Ratio.root'
        #pp_File = 'InputData/17q_SE_EMax_Correlation_GMB_Ratio.root'
        print(pp_File)
        print(pPb_File)

        return Files,purity,Systems #define as regular arrays outside of this function

pT_Bins = [12,15,19,26,40]
zT_Bins = [0.05, 0.07670637, 0.11767734, 0.18053204, 0.27695915, 0.42489062, 0.65183634, 1]
