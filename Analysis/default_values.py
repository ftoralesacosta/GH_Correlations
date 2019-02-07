###Files###
NN_pPb_File = 'InputData/pPb_SE_NN_Correlation_GMB_Ratio.root'
NN_pp_File = 'InputData/pp_SE_NN_Correlation_GMB_Ratio.root'

NN_purity = [0.276899, 0.358741, 0.456807, 0.476192]

L0_pPb_File = 'InputData/pPb_SE_L0_Correlation_GMB_Ratio.root'
L0_pp_File = 'InputData/pp_SE_L0_Correlation_GMB_Ratio.root'

L0_purity = [0.276791, 0.358627, 0.456378, 0.476337]


###Plotting Parameters##
delta_phi_centers= np.array([0.5890486225480862, 0.9817477042468103, 1.3744467859455345, 1.7671458676442586, 2.1598449493429825, 2.552544031041707, 2.945243112740431])

zT_centers = np.zeros(NzT)
zT_widths = np.zeros(NzT)
for ztbin in range(zT_offset,NzT+zT_offset):
    zT_centers[ztbin-zT_offset] = (zTbins[ztbin]+ zTbins[ztbin+1])/2
    zT_widths[ztbin-zT_offset] = (zTbins[ztbin+1]-zTbins[ztbin])/2

phi_width = [0.39269908169872414/2]*len(delta_phi_centers)

ue_error_bar = np.array([0,0.05,0.1,0.15,0.2,0.25,0.3,0.39269908169872414,2*0.39269908169872414])

###Binning
zTbins = [0.05, 0.07670637, 0.11767734, 0.18053204, 0.27695915, 0.42489062, 0.65183634, 1]
pTbins = [12, 15, 19, 26, 40]

NzT = 4
zT_offset = 2


###Corrections###
Corrections = [1,1.007,0.982,0.970,0.942,0.830,0.640]
oneminFake = [1,0.982,0.980,0.978,0.970,0.915,0.812]
