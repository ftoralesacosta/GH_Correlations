#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <THStack.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>

#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_ISO_ITS_04_SUB, CLUSTER_ISO_TPC_04_SUB, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};


Float_t Get_Purity_ErrFunction(Float_t pT_GeV, std::string deviation,bool Is_pp, bool Use_TPC) {

  // fprintf(stderr,"\n USE TPC FLAG = ");
  // fputs(Use_TPC ? "true \n" : "false \n", stdout);

  // fprintf(stderr,"\n Proton-Proton FLAG = ");
  // fputs(Is_pp ? "true \n" : "false \n", stdout);
  
  Float_t purity_val = 0;

  Float_t par[3] = {0.549684905516,
		     8.44338685256,
		    13.3454091464};

  //purity_val = par[0]*TMath::Erf((pT_GeV-par[1])/par[2]);
  //fprintf(stderr,"\n %d: Purity ITS = %f",__LINE__,purity_val);
  
  if(Use_TPC){ //pPb
    par[0] = 0.569429959156;
    par[1] = 8.1906528936;
    par[2] = 11.8993694765;

    // purity_val = par[0]*TMath::Erf((pT_GeV-par[1])/par[2]);
    // fprintf(stderr,"\n %d: Purity TPC = %f",__LINE__,purity_val);
    // fprintf(stderr,"\n %d: Purity TPC = %f",__LINE__,purity_val);

  }
		    
  if (Is_pp){
    // fprintf(stderr,"\n");
    // fprintf(stderr,"\n PP SELECTED \n");
    par[0] = 0.500229283252;
    par[1] = 9.016920902665;
    par[2] = 11.373299838596;
  }
  //order of conditionals ensures pp always overwrights TPC purity
  
  if (strcmp(deviation.data(),"Plus")==0){
    par[0] = 0.60750016509;
    par[1] = 7.05184155403;
    par[2] = 13.6116163603;
  }

  if (strcmp(deviation.data(),"Minus")==0){
    par[0] = 0.479958593235;
    par[1] = 9.05392932723;
    par[2] = 10.2061359452;
  }

  purity_val = par[0]*TMath::Erf((pT_GeV-par[1])/par[2]);
  //fprintf(stderr,"\n %d: Cluster pT = %f, Purity = %f \n",__LINE__,pT_GeV,purity_val);
  return purity_val;
}


int main(int argc, char *argv[])
{
  if (argc < 3) {
    fprintf(stderr,"Format: [command] [root file] [pp or pPb] \n");
    exit(EXIT_FAILURE);
  }
  int dummyc = 1;
  char **dummyv = new char *[1];

  dummyv[0] = strdup("main");


  bool Is_pp = false;


  std::string coll_system = argv[2];
  if (strcmp(coll_system.c_str(), "pp") == 0)
    Is_pp = true;

  if (Is_pp)
      fprintf(stderr,"\n PROTON PROTON SELECTED \n \n");
  
  //Config File
  FILE* config = fopen("../Corr_config.yaml", "r");
  double DNN_min = 0;
  double DNN_max = 0;
  double DNN_Bkgd = 0;
  double Lambda0_cut = 0; 
  double Emax_min = 0;
  double Emax_max = 0;
  double pT_min = 0;
  double pT_max = 0;
  double Eta_max = 0;
  double Cluster_min = 0;
  float Cluster_DtoBad = 0;
  UChar_t Cluster_NLocal_Max = 0;
  double EcrossoverE_min = 0;
  double cluster_time = 20;

  bool do_pile = false;

  float track_pT_min = 0.0;
  float track_pT_max = 0.0;
  int Track_Cut_Bit = 0;
  int track_chi_max = 0;
  double iso_max = 0;
  double noniso_min = 0;
  double noniso_max = 0;
  //double deta_max = 0;
  isolationDet determiner = CLUSTER_ISO_ITS_04;
  int n_eta_bins = 0;
  int n_phi_bins = 0;  
  std::string shower_shape = "DNN";
  std::string purity_deviation = "None";

  bool TPC_Iso_Flag = false;
  
  // Zt bins
  //FIXME: Will have to likely set nztbins first, then initialize array
  int nztbins = 7;
  float* ztbins;
  ztbins = new float[nztbins+1];
  ztbins[0] = 0.0; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;
   
  int nptbins = 3;
  float* ptbins;
  ptbins = new float[nptbins+1];
  ptbins[0] = 10.0; ptbins[1] = 11; ptbins[2] = 12.5; ptbins[3] = 16;

  //FIXME: Obviously needs to be put in a header file.
  // Loop through config file
  char line[MAX_INPUT_LENGTH];
  while (fgets(line, MAX_INPUT_LENGTH, config) != NULL) {
      if (line[0] == '#') continue;
      
      char key[MAX_INPUT_LENGTH];
      char dummy[MAX_INPUT_LENGTH];
      char value[MAX_INPUT_LENGTH];
      
      // Cap off key[0] and value[0] with null characters and load the key, dummy-characters, and value of the line into their respective arrays
      key[0] = '\0';
      value[0] = '\0';
      sscanf(line, "%[^:]:%[ \t]%100[^\n]", key, dummy, value);
      
      //Read Config File: Detect Keys 
      if (strcmp(key, "DNN_min") == 0) {
          DNN_min = atof(value);
          std::cout << "DNN_min: " << DNN_min << std::endl; }

      else if (strcmp(key, "DNN_max") == 0) {
          DNN_max = atof(value);
          std::cout << "DNN_max: " << DNN_max << std::endl; }

      else if (strcmp(key, "DNN_BKGD") == 0) {
	DNN_Bkgd = atof(value);
	std::cout << "DNN_BKGD: " << DNN_Bkgd << std::endl; }

      else if (strcmp(key, "Lambda0_cut") == 0){
	Lambda0_cut = atof(value);
	std::cout << "Lambda0_cut: " << Lambda0_cut << std::endl;}

      if (strcmp(key, "EMax_EClus_min") == 0) {
	Emax_min = atof(value);
	std::cout << "EMax_EClus_min:" << Emax_min << std::endl; }

      else if (strcmp(key, "EMax_EClus_max") == 0) {
	Emax_max = atof(value);
	std::cout << "EMax_EClus_max: " << Emax_max << std::endl; }

      else if (strcmp(key, "pT_min") == 0) {
          pT_min = atof(value);
          std::cout << "pT_min: " << pT_min << std::endl; }

      else if (strcmp(key, "pT_max") == 0) {
          pT_max = atof(value);
          std::cout << "pT_max: " << pT_max << std::endl; }

      else if (strcmp(key, "Eta_max") == 0) {
          Eta_max = atof(value);
          std::cout << "Eta_max: " << Eta_max << std::endl;
      }
      else if (strcmp(key, "Cluster_min") == 0) {
          Cluster_min = atof(value);
          std::cout << "Cluster_min: " << Cluster_min << std::endl; }

      else if (strcmp(key, "Cluster_dist_to_bad_channel") == 0){
	Cluster_DtoBad = atof(value);
	std::cout << "Cluster_DtoBad: "<< Cluster_DtoBad << std::endl;}

      else if (strcmp(key, "Cluster_N_Local_Maxima") == 0){
        Cluster_NLocal_Max = atof(value);
	std::cout << "Cluster_NLocal_Max: "<< Cluster_NLocal_Max << std::endl;}

      else if (strcmp(key, "EcrossoverE_min") == 0) {
          EcrossoverE_min = atof(value);
          std::cout << "EcrossoverE_min; " << EcrossoverE_min << std::endl; }
      
      else if (strcmp(key, "Cluster_Time") == 0){
        cluster_time = atof(value);
        std::cout << "Cluster_Time: "<< cluster_time << std::endl;}
      
      else if (strcmp(key, "iso_max") == 0) {
          iso_max = atof(value);
          std::cout << "iso_max: " << iso_max << std::endl; }

      else if (strcmp(key, "noniso_min") == 0) {
          noniso_min = atof(value);
          std::cout << "noniso_min: " << noniso_min << std::endl; }

      else if (strcmp(key, "noniso_max") == 0) {
          noniso_max = atof(value);
          std::cout << "noniso_max: " << noniso_max << std::endl; }

      else if (strcmp(key, "do_pileup_cut") == 0) {
	if (strcmp(value,"true") == 0)
	  do_pile = true;
	std::cout << "do_pileup_cut: " << do_pile << std::endl; }
      
      // else if (strcmp(key, "deta_max") == 0) {
      //     deta_max = atof(value);
      //     std::cout << "deta_max: " << deta_max << std::endl; }

      else if (strcmp(key, "N_Phi_Bins") == 0) {
	n_phi_bins = atoi(value);
	std::cout << "Number of Phi Bins: " << n_phi_bins << std::endl; }

      else if (strcmp(key, "N_Eta_Bins") == 0) {
        n_eta_bins = atoi(value);
	std::cout << "Number of Eta Bins: " << n_eta_bins << std::endl; }

      else if (strcmp(key, "Track_pT_Min") == 0) {
          track_pT_min = atof(value);
          std::cout << "Track Min pT: " << track_pT_min << std::endl; }

      else if (strcmp(key, "Track_pT_Max") == 0) {
          track_pT_max = atof(value);
          std::cout << "Track Max pT: " << track_pT_max << std::endl; }

      else if (strcmp(key, "Track_Cut_Bit") == 0) {
          Track_Cut_Bit = atoi(value);
          std::cout << "Track Cut Bit: " << Track_Cut_Bit << std::endl; }

      else if (strcmp(key, "Track_Chi_Max") == 0) {
          track_chi_max = atoi(value);
          std::cout << "Track Chi^2 Max: " << track_chi_max << std::endl; }

      
      else if (strcmp(key, "Zt_bins") == 0) {
          nztbins = -1;
          for (const char *v = value; *v != ']';) {
              while (*v != ']' && !isdigit(*v)) v++;
	      nztbins++;
              while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }

          ztbins = new float[nztbins + 1];
          int i = 0;
          for (const char *v = value; *v != ']' ;) {
              while (*v != ']' && !isdigit(*v)) v++;
              ztbins[i] = atof(v);
              i++;              
              while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }

          std::cout << "Number of Zt bins: " << nztbins << std::endl << "Zt bins: {";
          for (int i = 0; i <= nztbins; i++)
	    std::cout << ztbins[i] << ", ";
          std::cout << "}\n"; 
      }

      else if (strcmp(key, "Pt_bins") == 0) {
          nptbins = -1;
          for (const char *v = value; *v != ']';) {
            while (*v != ']' && !isdigit(*v)) v++;
	    nptbins++;
	    while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
	
          ptbins = new float[nptbins + 1];
	  int i = 0;
	  for (const char *v = value; *v != ']' ;) {
	     while (*v != ']' && !isdigit(*v))  v++;
	     ptbins[i] = atof(v);
	     i++;
	     while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
	
	  std::cout << "Number of Pt bins: " << nptbins << std::endl << "Pt bins: {";
          for (int i = 0; i <= nptbins; i++)
	    std::cout << ptbins[i] << ", ";
	  std::cout << "}\n";
      }

      else if (strcmp(key, "Cluster_isolation_determinant") == 0) {
          if (strcmp(value, "cluster_iso_tpc_04") == 0){
              determiner = CLUSTER_ISO_TPC_04;
              std::cout << "Isolation Variable: cluster_iso_tpc_04" << std::endl; }

          else if (strcmp(value, "cluster_iso_its_04") == 0){
              determiner = CLUSTER_ISO_ITS_04;
              std::cout << "Isolation Variable: cluster_iso_its_04" << std::endl; }

          else if (strcmp(value, "cluster_iso_its_04_sub") == 0){
              determiner = CLUSTER_ISO_ITS_04_SUB;
              std::cout << "Isolation Variable: cluster_iso_its_04_sub" << std::endl; }

	  else if (strcmp(value, "cluster_iso_tpc_04_sub") == 0){
              determiner = CLUSTER_ISO_TPC_04_SUB;
	      TPC_Iso_Flag = true;
	      std::cout << "Isolation Variable: cluster_iso_tpc_04_sub" << std::endl; }
	  
          else if (strcmp(value, "cluster_frixione_tpc_04_02") == 0){
              determiner = CLUSTER_FRIXIONE_TPC_04_02;
              std::cout << "Isolation Variable: cluster_frixione_tpc_04_02" << std::endl; }

          else if (strcmp(value, "cluster_frixione_its_04_02") == 0){
              determiner = CLUSTER_FRIXIONE_ITS_04_02;
              std::cout << "Isolation Variable: cluster_frixione_its_04_02" << std::endl; }

          else {
              std::cout << "ERROR: Cluster_isolation_determinant in configuration file must be \"cluster_iso_tpc_04\", \"cluster_iso_its_04\", \"cluster_frixione_tpc_04_02\", or \"cluster_frixione_its_04_02\"" << std::endl << "Aborting the program" << std::endl;
              exit(EXIT_FAILURE); }
      }

      else if (strcmp(key, "Shower_Shape") == 0){
	  shower_shape = value;
	  std::cout<<"Shower Shape: "<<shower_shape.data()<<std::endl;
	  //if (strcmp(shower_shape.data(),"Lambda")== 0) std::cout<<"test worked"<<std::endl;
      }

      else if (strcmp(key, "Purity_Dev") == 0){
	purity_deviation = value;
	std::cout<<"Purity Deviation Change: "<<purity_deviation.data()<<std::endl;
      }

      else std::cout << "WARNING: Unrecognized keyvariable " << key << std::endl;
  
  }
  //end Config Loop

  fclose(config);
  
  for (int i = 0; i <= nztbins; i++)
    std::cout << "zt bound: " << ztbins[i] << std::endl;
  for (int i = 0; i <= nptbins; i++)
    std::cout << "pt bound: " << ptbins[i] << std::endl;


  //HISTOGRAMS
  TCanvas canvas("canvas", "");

  TH1D h_purity("h_purity","purity distribution",100,0,1);

  TH1D *h_cluster_phi = new TH1D("cluster_phi","#phi distribution of paired clusters",32,1,M_PI);  
  TH1D *h_cluster_eta = new TH1D("cluster_eta","#eta distribution of paired clusters",28,-0.8,0.8);  

  h_cluster_phi->Sumw2();
  h_cluster_eta->Sumw2();

  //Purity Handling

  //Following purities for pT range: 12.5,13.2,14.4,15.8
  int  N_pT_Ranges = 5;
  float pT_Ranges[5] = {12.0,15.0,20.0,25.0,40.0};
  float purities[4] = {0};
  float purity_Uncertainties[4] = {0};
  float Cluster_Purity = 0;
  float Cluster_Purity_Uncertainty = 0;

  // if (strcmp(shower_shape.data(), "DNN") == 0){
  //   purities[0] = 0.207;
  //   purities[1] = 0.255;
  //   purities[2] = 0.326; 
  //   purities[3] = 0.372;
  //   purities[4] = 0.447;
  //   purities[5] = 0.502;
  //   purities[6] = 0.533; //Extrapolating last bin
  //   purities[7] = 0.533;


  //   purity_Uncertainties[0] = 0.030;
  //   purity_Uncertainties[1] = 0.037;
  //   purity_Uncertainties[2] = 0.042;
  //   purity_Uncertainties[3] = 0.050;
  //   purity_Uncertainties[4] = 0.056;
  //   purity_Uncertainties[5] = 0.062;
  //   purity_Uncertainties[6] = 0.058;
  //   purity_Uncertainties[7] = 0.058;

  // }

  if (strcmp(shower_shape.data(), "Lambda") == 0){
    purities[0] = 0.206;
    purities[1] = 0.341;
    purities[2] = 0.471;
    purities[3] = 0.546;

    purity_Uncertainties[0]= 0.0301;
    purity_Uncertainties[1]= 0.0305;
    purity_Uncertainties[2]= 0.0503;
    purity_Uncertainties[3]= 0.0572;

    if (Is_pp){

      purities[0] = 0.198;
      purities[1] = 0.318;
      purities[2] = 0.470;
      purities[3] = 0.487;
      
      purity_Uncertainties[0]= 0.0490;
      purity_Uncertainties[1]= 0.0400;
      purity_Uncertainties[2]= 0.0712;
      purity_Uncertainties[3]= 0.1221;
      

    }
  }



  // std::array< float,38 > track_pT_Correction_bins = {0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
  // 						0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90,
  // 						0.95, 1.00, 1.10, 1.20, 1.40, 1.60, 1.80, 2.00,
  // 						2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.60, 4.00,
  // 						5.00, 6.00, 8.00, 10.0, 13.0, 20.0};

  
  // //p-Pb
  // float pPb_Efficiency[37] = {0.701584, 0.794149, 0.811056, 0.824281, 0.833784, 0.8411, 0.846235, 0.850521,
  // 			      0.853948, 0.855174, 0.857245, 0.858315, 0.858334, 0.859078, 0.860819, 0.860223,
  // 			      0.860368, 0.860849, 0.862532, 0.864867, 0.867081, 0.869168, 0.869016, 0.872149,
  // 			      0.872669, 0.874919, 0.876616, 0.876491, 0.87365, 0.877567, 0.8777, 0.882617,
  // 			      0.88169, 0.881348, 0.880515, 0.87605, 0.882653};

  // float pPb_FakeRate[37] = {0.0444643, 0.0329245, 0.0274773, 0.0241898, 0.0222288, 0.0207158, 0.0195863,
  // 			    0.0189526, 0.0184162, 0.0180536, 0.0175808, 0.0177315, 0.017413, 0.0172551,
  // 			    0.0175637, 0.017715, 0.0172888, 0.0178326, 0.0178794, 0.0179673, 0.0177509,
  // 			    0.0182945, 0.0181996, 0.0188647, 0.0190429, 0.0202744, 0.0209381, 0.0213728,
  // 			    0.0237848, 0.0251898, 0.0297774, 0.0350703, 0.0518793, 0.0804111, 0.137953,
  // 			    0.214123, 0.301459};

  // float pPb_Smearing_Correction[37] = {0.894311, 0.994658, 1.00122, 1.0087, 1.00914, 1.01009, 1.00939,
  // 				       1.01073, 1.00998, 1.00622, 1.00519, 1.00425, 1.00182, 0.99553,
  // 				       0.992115, 0.991085, 0.989869, 0.988516,0.987711, 0.989906, 0.990167,
  // 				       0.994517, 0.994247, 0.99273, 0.983717, 0.986158, 0.978648, 0.986939,
  // 				       0.963266, 0.959108, 0.945863, 0.923434, 0.896206, 0.86263, 0.71868,
  // 				       0.629909, 0.425061};
  
  std::array< float,52 > track_pT_Correction_bins = {0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14, 15};

    int N_Track_pT_Bins = track_pT_Correction_bins.size()-1;

    float pPb_Efficiency[51],pPb_FakeRate[51],pPb_Smearing_Correction[51];
    
  //Chi < 36:
   float ITS_Chi_36_pPb_Efficiency[51] = {0.706923, 0.801196, 0.818368, 0.832127, 0.840773, 0.847301, 0.851753, 0.855662, 0.858772, 0.859778, 0.861794, 0.862719, 0.862874, 0.863599, 0.865444, 0.86517, 0.865167, 0.865733, 0.867507, 0.869371, 0.869943, 0.871633, 0.872051, 0.873205, 0.873924, 0.873286, 0.873352, 0.876119, 0.876512, 0.878243, 0.879848, 0.879726, 0.876458, 0.879516, 0.882163, 0.882294, 0.878792, 0.88552, 0.883079, 0.878306, 0.893859, 0.879929, 0.89569, 0.881874, 0.886792, 0.874687, 0.86, 0.899281, 0.885057, 0.892857, 0.894737};

   float ITS_Chi_36_pPb_FakeRate[51] = {0.06586, 0.0500067, 0.0407897, 0.0352041, 0.0315345, 0.0288352, 0.0269184, 0.0257262, 0.0247625, 0.0240762, 0.0235506, 0.023464, 0.0232308, 0.0230312, 0.0235422, 0.0235193, 0.022949, 0.0234972, 0.023682, 0.0238825, 0.0233844, 0.023943, 0.0240052, 0.0246918, 0.0249658, 0.0244902, 0.0250503, 0.025728, 0.0270794, 0.0292155, 0.0310878, 0.0324949, 0.036317, 0.0385375, 0.0393392, 0.0471944, 0.0530674, 0.0587577, 0.0699657, 0.0886427, 0.10224, 0.143758, 0.15557, 0.18744, 0.236691, 0.266279, 0.323333, 0.396181, 0.346505, 0.368421, 0.427184};

   float ITS_Chi_36_pPb_Smearing_Correction[51] = {0.891805, 0.993346, 1.00097, 1.01, 1.01008, 1.01077, 1.00976, 1.01116, 1.01027, 1.00654, 1.00555, 1.00457, 1.0023, 0.995985, 0.99255, 0.991825, 0.990458, 0.989069, 0.988395, 0.990534, 0.990073, 0.987179, 0.994019, 0.990603, 0.998289, 0.993022, 0.99427, 0.991535, 0.981762, 0.983794, 0.974509, 0.981303, 0.956655, 0.958883, 0.943399, 0.943427, 0.922941, 0.908759, 0.891862, 0.873412, 0.850649, 0.811177, 0.812353, 0.79159, 0.678889, 0.57686, 0.555556, 0.50813, 0.385, 0.316456, 0.443478};

   if (track_chi_max == 36){
     fprintf(stderr,"%d: Using Chi_Max < 36 Corrections \n",__LINE__);
     memcpy(pPb_Efficiency,ITS_Chi_36_pPb_Efficiency,sizeof(ITS_Chi_36_pPb_Efficiency));
     memcpy(pPb_FakeRate,ITS_Chi_36_pPb_FakeRate,sizeof(ITS_Chi_36_pPb_FakeRate));
     memcpy(pPb_Smearing_Correction,ITS_Chi_36_pPb_Smearing_Correction,sizeof(ITS_Chi_36_pPb_Smearing_Correction));
   }

   //Chi < 2
   float ITS_Chi_2_pPb_Efficiency[51] = {0.679585, 0.76674, 0.78565, 0.800042, 0.812381, 0.822312, 0.829837, 0.835912, 0.840904, 0.843102, 0.846413, 0.847949, 0.848675, 0.849882, 0.851767, 0.851662, 0.851702, 0.853024, 0.854968, 0.857285, 0.858209, 0.860071, 0.860713, 0.862577, 0.863482, 0.862668, 0.862379, 0.866158, 0.867516, 0.869688, 0.871099, 0.870675, 0.867989, 0.871116, 0.873867, 0.873728, 0.870625, 0.878394, 0.876006, 0.87091, 0.887415, 0.873455, 0.885345, 0.872369, 0.87373, 0.869674, 0.852, 0.899281, 0.873563, 0.892857, 0.877193};

    float ITS_Chi_2_pPb_FakeRate[51] = {0.0288587, 0.0227023, 0.0201382, 0.0185892, 0.0175023, 0.0167976, 0.0160735, 0.0157722, 0.015473, 0.0152811, 0.0149019, 0.0148865, 0.014695, 0.0146859, 0.0149738, 0.0150937, 0.014766, 0.0150073, 0.0149588, 0.0150989, 0.0149009, 0.0144755, 0.0144847, 0.0148771, 0.014925, 0.0149947, 0.0149865, 0.0154721, 0.0154909, 0.015774, 0.0164298, 0.0160145, 0.0177727, 0.0184202, 0.0185393, 0.0206371, 0.0210202, 0.0203409, 0.0239161, 0.0273942, 0.0345207, 0.0436947, 0.0468119, 0.0464516, 0.0710901, 0.0819367, 0.0946372, 0.15566, 0.15, 0.158333, 0.179487};

    float ITS_Chi_2_pPb_Smearing_Correction[51] = {0.900361, 0.99487, 1.00041, 1.0064, 1.00725, 1.00898, 1.00901, 1.01101, 1.01011, 1.00651, 1.00588, 1.00471, 1.00274, 0.996152, 0.992602, 0.991865, 0.990361, 0.98932, 0.988507, 0.990986, 0.990221, 0.987836, 0.994203, 0.991718, 0.999723, 0.994458, 0.995187, 0.993518, 0.985992, 0.988016, 0.980657, 0.988024, 0.96586, 0.969533, 0.955084, 0.959389, 0.9422, 0.93353, 0.925015, 0.912635, 0.912671, 0.883333, 0.896161, 0.894224, 0.790026, 0.730526, 0.766187, 0.722543, 0.584615, 0.526316, 0.793651};
    
    if (track_chi_max == 2){
     fprintf(stderr,"%d: Using Chi_Max < 2 Corrections \n",__LINE__);
     memcpy(pPb_Efficiency,ITS_Chi_2_pPb_Efficiency,sizeof(ITS_Chi_2_pPb_Efficiency));
     memcpy(pPb_FakeRate,ITS_Chi_2_pPb_FakeRate,sizeof(ITS_Chi_2_pPb_FakeRate));
     memcpy(pPb_Smearing_Correction,ITS_Chi_2_pPb_Smearing_Correction,sizeof(ITS_Chi_2_pPb_Smearing_Correction));
   }


    
    fprintf(stderr,"NUMBER OF TRACK PT BINS IS %i \n \n",N_Track_pT_Bins);

    if (track_chi_max == 3) {

      fprintf(stderr, "Using track correction numbers for Track Chi < 3 \n");
    }
  
      float ITS_Chi_3_pPb_Efficiency[51] = {0.680807, 0.769416, 0.786667, 0.801124, 0.812683, 0.822154, 0.829054, 0.833275, 0.837709, 0.842318, 0.838253, 0.844133, 0.845795, 0.845452, 0.8437, 0.844988, 0.850204, 0.851091, 0.851471, 0.851059, 0.852722, 0.856077, 0.856633, 0.855644, 0.85509, 0.85838, 0.856906, 0.855471, 0.857095, 0.860265, 0.860686, 0.866243, 0.861802, 0.864242, 0.865098, 0.861285, 0.858213, 0.861279, 0.856347, 0.866282, 0.853067, 0.860997, 0.855398, 0.860901, 0.842331, 0.861237, 0.856387, 0.861898, 0.865817, 0.828281, 0.852313};

      float ITS_Chi_3_pPb_FakeRate[51] = {0.03276, 0.0245696, 0.0221121, 0.0196002, 0.0184712, 0.0173392, 0.0164418, 0.0164659, 0.0160136, 0.015766, 0.0150485, 0.015234, 0.014941, 0.0154199, 0.0148796, 0.0141062, 0.0148753, 0.0152038, 0.0143306, 0.0149377, 0.0135737, 0.0146271, 0.0152528, 0.0151404, 0.0147567, 0.0152333, 0.0177789, 0.0147575, 0.0168503, 0.0169324, 0.0164895, 0.0151945, 0.0167499, 0.0166498, 0.0141885, 0.0154582, 0.0214637, 0.0174561, 0.0215349, 0.0253443, 0.0208158, 0.0263207, 0.0179423, 0.0341796, 0.0339603, 0.0458528, 0.0268865, 0.0376692, 0.0603408, 0.0873913, 0.0447524};

      float ITS_Chi_3_pPb_Smearing_Correction[51] = {0.897708, 0.99395, 1.00205, 1.00427, 1.0062, 1.01035, 1.00878, 1.01121, 1.0075, 1.01121, 1.00012, 1.00833, 1.00139, 0.999896, 0.993365, 0.989383, 1.00257, 0.990837, 0.995561, 0.996833, 0.984558, 0.99655, 0.999187, 0.998964, 1.00884, 0.990964, 0.989008, 0.990163, 0.996841, 0.988733, 0.976975, 1.00046, 0.957124, 0.959472, 1.01762, 0.968777, 0.922066, 0.984912, 0.965775, 0.938275, 0.982422, 0.924246, 0.89762, 0.936461, 0.875103, 0.987666, 0.854273, 0.877012, 0.727192, 0.860683, 0.78266};

      if (track_chi_max == 3){
      fprintf(stderr,"%d: Using Chi_Max < 3 Corrections \n",__LINE__);
      memcpy(pPb_Efficiency,ITS_Chi_3_pPb_Efficiency,sizeof(ITS_Chi_3_pPb_Efficiency));
      memcpy(pPb_FakeRate,ITS_Chi_3_pPb_FakeRate,sizeof(ITS_Chi_3_pPb_FakeRate));
      memcpy(pPb_Smearing_Correction,ITS_Chi_3_pPb_Smearing_Correction,sizeof(ITS_Chi_3_pPb_Smearing_Correction));
    }

      
      std::array< float,52 > TPC_track_pT_Correction_bins = {0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
							     0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
							     0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
							     1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8,
							     3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5,
							     5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14,15};

    int N_TPC_Track_pT_Bins = TPC_track_pT_Correction_bins.size()-1;

    //TPC Chi<36

    float TPC_pPb_Efficiency[N_TPC_Track_pT_Bins] = {0.547368, 0.645805, 0.729414, 0.754735, 0.743207, 0.74378, 0.836156, 0.846914, 0.854331, 0.858756, 0.904079, 0.864668, 0.865609, 0.86637, 0.868146, 0.895326, 0.871075, 0.869368, 0.874902, 0.878735, 0.880171, 0.876702, 0.876428, 0.878704, 0.878008, 0.874098, 0.871206, 0.861439, 0.854897, 0.851077, 0.847729, 0.839184, 0.841624, 0.841437, 0.841226, 0.83761, 0.838738, 0.844055, 0.838054, 0.839841, 0.836119, 0.841824, 0.833175, 0.827837, 0.829056, 0.830772, 0.817571, 0.812508, 0.81811, 0.813479, 0.804271};

    float TPC_pPb_FakeRate[N_TPC_Track_pT_Bins] = {0.0385663, 0.0265005, 0.0224742, 0.0207308, 0.0199845, 0.0204949, 0.0194557, 0.0186314, 0.0177608, 0.0174342, 0.0161894, 0.015664, 0.0161282, 0.0160003, 0.0158429, 0.0149451, 0.0150958, 0.0150021, 0.0148423, 0.014283, 0.0140022, 0.0139122, 0.0134801, 0.0132932, 0.0138496, 0.013443, 0.0134683, 0.0133246, 0.0134612, 0.0122121, 0.0125056, 0.012649, 0.0124756, 0.0121123, 0.0102316, 0.0126523, 0.0115912, 0.0108272, 0.012054, 0.0101419, 0.0105359, 0.011024, 0.00824151, 0.0105761, 0.00958554, 0.010172, 0.0101555, 0.0105039, 0.0114032, 0.00975607, 0.00928096};

      
      float TPC_pPb_Smearing_Correction[N_TPC_Track_pT_Bins] = {0.964013, 0.994686, 0.998026, 0.996843, 0.997778, 0.999472, 1.00381, 1.00472, 1.00455, 1.00453, 1.00413, 1.00522, 1.00504, 1.00531, 1.00373, 1.00708, 1.00725, 1.00329, 1.00379, 1.00425, 1.00252, 1.00493, 1.00465, 1.00418, 1.00159, 1.00254, 1.00779, 1.00462, 1.00395, 1.00784, 1.00474, 1.00217, 1.0055, 1.00271, 0.989254, 0.997873, 0.98817, 0.988905, 0.990132, 0.986167, 0.987632, 0.979275, 0.982885, 0.985685, 0.973379, 0.986632, 0.982638, 0.991744, 0.965312, 0.995285, 1.02297};

      //TPC Chi < 3

      // float TPC_pPb_Efficiency[N_TPC_Track_pT_Bins] = {0.506012, 0.593831, 0.634368, 0.660664, 0.680452, 0.691147, 0.778763, 0.792615, 0.802296, 0.822707, 0.853774, 0.815112, 0.817322, 0.819871, 0.823506, 0.851564, 0.830688, 0.831277, 0.841348, 0.841784, 0.843824, 0.841167, 0.841117, 0.845067, 0.843809, 0.840738, 0.835867, 0.829102, 0.822877, 0.819697, 0.818317, 0.810206, 0.811689, 0.810477, 0.811292, 0.808944, 0.806744, 0.812477, 0.80801, 0.808969, 0.806741, 0.808339, 0.798921, 0.794134, 0.796974, 0.798189, 0.77824, 0.771855, 0.773063, 0.779684, 0.764668};

      // float TPC_pPb_FakeRate[N_TPC_Track_pT_Bins] = {0.0301071, 0.0221534, 0.0183075, 0.0173483, 0.0172073, 0.0165191, 0.0158807, 0.0154625, 0.0147511, 0.0145365, 0.0135284, 0.0133251, 0.0137556, 0.0136217, 0.013694, 0.012744, 0.0130261, 0.012958, 0.0129025, 0.0124018, 0.0120807, 0.0122043, 0.0119859, 0.0116575, 0.0121953, 0.0120978, 0.0120622, 0.0120504, 0.0121441, 0.0110372, 0.0113105, 0.0113682, 0.0114341, 0.0110452, 0.0097054, 0.0112374, 0.0105443, 0.0101021, 0.0110572, 0.00971418, 0.00945253, 0.0102266, 0.00774202, 0.0102232, 0.00941761, 0.00829773, 0.00927938, 0.00874754, 0.0106023, 0.00867711, 0.0080265};

      
      // float TPC_pPb_Smearing_Correction[N_TPC_Track_pT_Bins] = {0.964056, 0.994522, 0.99792, 0.996716, 0.997695, 0.999594, 1.00345, 1.00438, 1.00446, 1.00467, 1.00432, 1.005, 1.00521, 1.0055, 1.00377, 1.00654, 1.00764, 1.00348, 1.0038, 1.00423, 1.00268, 1.00515, 1.00384, 1.00475, 1.00147, 1.00281, 1.00755, 1.00483, 1.00399, 1.00759, 1.0055, 1.00163, 1.0052, 1.00224, 0.988787, 1.00121, 0.98637, 0.990902, 0.989828, 0.986665, 0.991921, 0.976119, 0.984256, 0.988111, 0.974778, 0.987447, 0.980335, 0.993181, 0.964971, 1.00398, 1.02836};

    
      //pp
    std::array< float,52 > pp_track_pT_Correction_bins {0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14, 15};

    float pp_Efficiency[51] = {0.68011, 0.76461, 0.778984, 0.796766, 0.807108, 0.819888, 0.82531, 0.831658, 0.833392, 0.835791, 0.83586, 0.846695, 0.848953, 0.846094, 0.842743, 0.841815, 0.847231, 0.848502, 0.85282, 0.856009, 0.854111, 0.857407, 0.861889, 0.863067, 0.85911, 0.866662, 0.871307, 0.865416, 0.872027, 0.878082, 0.880938, 0.882207, 0.880407, 0.888934, 0.888041, 0.890855, 0.894479, 0.881659, 0.894079, 0.894791, 0.893812, 0.899538, 0.912514, 0.909233, 0.890754, 0.904707, 0.901652, 0.908604, 0.905241, 0.916258, 0.941168};

    float pp_FakeRate[51] = {0.0401015, 0.0341742, 0.0293513, 0.0275427, 0.0258545, 0.0258758, 0.0264281, 0.0230541, 0.0242865, 0.0237406, 0.0246322, 0.021845, 0.0214785, 0.0219201, 0.0242008, 0.0225091, 0.0228031, 0.0217989, 0.0218427, 0.0220313, 0.0181306, 0.0209162, 0.0184084, 0.0219567, 0.0208506, 0.0206115, 0.0188644, 0.0209044, 0.0208765, 0.0169143, 0.0184842, 0.0161239, 0.0184484, 0.0196252, 0.0223822, 0.0209666, 0.019068, 0.0270674, 0.0239608, 0.0227336, 0.0299992, 0.0229086, 0.0264195, 0.0256361, 0.0237618, 0.0249181, 0.0207626, 0.0281898, 0.0214536, 0.0266104, 0.0283818};

    float pp_Smearing_Correction[51] = {0.901894, 0.994355, 0.996727, 1.01158, 0.997977, 0.999246, 1.01004, 1.01041, 0.993775, 1.01544, 0.994904, 1.00169, 1.02182, 1.01117, 0.982041, 0.989534, 1.00133, 1.00846, 0.979539, 1.03044, 0.988919, 1.00315, 0.987368, 0.996547, 1.01726, 1.00539, 1.00079, 1.00315, 0.992049, 0.99654, 0.978189, 0.988968, 1.0146, 0.99889, 1.00061, 1.03534, 0.955467, 1.0114, 1.0047, 0.96447, 0.971297, 1.08372, 1.05976, 1.01012, 0.905184, 1.00907, 1.00258, 1.12743, 0.889183, 0.942362, 1.03666};
    
    int N_Track_pT_Bins_pp = pp_track_pT_Correction_bins.size()-1;
    fprintf(stderr,"NUMBER OF TRACK PT BINS IS %i \n \n",N_Track_pT_Bins_pp);
    

  TH2D* Corr[nztbins*nptbins];
  TH2D* IsoCorr[nztbins*nptbins];
  TH2D* BKGD_IsoCorr[nztbins*nptbins];
  TH2D* BKGD_IsoCorr_UW[nztbins*nptbins];
  TH2D* Inclusive_Corr[nztbins*nptbins];

  TH1D* H_Signal_Triggers[nptbins];
  TH1D* H_BKGD_Triggers[nptbins];
  TH1D* H_Signal_Overlap_Triggers[nptbins];
  TH1D* H_BKGD_Overlap_Triggers[nptbins];
  TH1D* Triggers[nptbins];
  TH1D* Inclusive_Triggers[nptbins];

  TH1F* H_Purities[nptbins];
  TH1F* H_Purity_Uncertainties[nptbins];

  TH2D * h_track_phi_eta[nztbins*nptbins];
  TH1D * h_track_eta[nztbins*nptbins];
  TH1D * h_track_phi[nztbins*nptbins];

  float N_Signal_Triggers = 0;
  float N_BKGD_Triggers = 0;
  
  TH1F* Signal_pT_Dist = new TH1F("Signal_pT_Dist","Cluster Pt Spectrum For Isolation (its_04) bins 0.55 < DNN < 0.85",(pT_max-pT_min)*2,pT_min,pT_max);
  TH1F* Signal_pT_Dist_OnlypTWeight = new TH1F("Signal_pT_Dist_OnlypTWeight","Cluster Pt Spectrum For Isolation (its_04) bins 0.55 < DNN < 0.85",(pT_max-pT_min)*2,pT_min,pT_max);
  
  TH1F* BKGD_pT_Dist = new TH1F("BKGD_pT_Dist","Cluster Pt Spectrum For Isolation (its_04) bins 0.0 < DNN < 0.3",(pT_max-pT_min)*2,pT_min,pT_max);
    TH1F* BKGD_pT_Dist_OnlypTWeight = new TH1F("BKGD_pT_Dist_OnlypTWeight","Cluster Pt Spectrum For Isolation (its_04) bins 0.0 < DNN < 0.3",(pT_max-pT_min)*2,pT_min,pT_max);
  TH1F* BKGD_pT_Dist_Weighted = new TH1F("BKGD_pT_Dist_Weighted","Weighted Cluster Pt Spectrum For Isolation (its_04) bins 0.0 < DNN < 0.3",(pT_max-pT_min)*2,pT_min,pT_max);

  Signal_pT_Dist->Sumw2();
  BKGD_pT_Dist->Sumw2();
  BKGD_pT_Dist_Weighted->Sumw2();

  TH1F hBR("hBR", "Isolated cluster, bkg region", 40, 10.0, 50.0);
  TH1F hweight("hweight", "Isolated cluster, signal region", 40, 10.0, 50.0);
  TH1F* Weights_Sum = new TH1F("Weights_Sum", "Sum of Weights. XBin Centers = pt Bin", nptbins,0.5,4.5);
  float weight_sum[nptbins] = {0};
  Float_t purity_weight_sum[nptbins] = {0};
  Float_t BR_purity_weight_sum[nptbins] = {0};

  TH1F Track_Weight_Spectra("Weight_Track_Spectra", "Weighted Track Spectra",30,0.5,8);
  TH1F Track_Raw_Spectra("Raw_Track_Spectra", "Raw Track Spectra",30,0.5,8);
  TH1F Chi_Distro("Chi_Track_Distribution","Distribution of track #chi^2",100,0,50);
  TH1F* Track_Uncertainty[nztbins];
    //("Efficiency Uncertainty","Single Track Efficiency zT Bin Uncertainty",nztbins,0.)

  for (int izt = 0; izt < nztbins; izt++){
    Track_Uncertainty[izt] = new TH1F(Form("Track_Uncertainty_zT_%1.0f_%1.0f",100*ztbins[izt],100*ztbins[izt+1]),Form("Track Uncertainty zT %1.2f %1.2f",ztbins[izt],ztbins[izt+1]),100,0.0,0.1);
    Track_Uncertainty[izt]->Sumw2();
  }//Used to implement the 8% Uncertainty for tracks between with pT 0.5-0.85GeV/c
  
  Track_Weight_Spectra.Sumw2();
  Track_Raw_Spectra.Sumw2();
  Chi_Distro.Sumw2();
  //  TH1F track_above_ten("tracks_above_ten","Tracks above 10 GeV/c",)
  Long64_t tracks_tenGev = 0;
    

  hweight.Sumw2();
  hBR.Sumw2();
  
    for (int ipt = 0; ipt <nptbins; ipt++) {
      H_Signal_Triggers[ipt] = new TH1D(
      Form("N_DNN%i_Triggers_pT%1.0f_%1.0f",1,ptbins[ipt],ptbins[ipt+1]),
      "Number of Isolated Photon Triggers", 2, -0.5,1.0);

      H_BKGD_Triggers[ipt] = new TH1D(
      Form("N_DNN%i_Triggers_pT%1.0f_%1.0f",2,ptbins[ipt],ptbins[ipt+1]),
      "Number of Isolated Low DNN Photon Triggers", 1, -0.5,3.5);

      H_Signal_Overlap_Triggers[ipt] = new TH1D(
      Form("N_DNN%i_Overlap_Triggers_pT%1.0f_%1.0f",1,ptbins[ipt],ptbins[ipt+1]),
      "Number of Isolated Low DNN Photon Triggers", 1, -0.5,3.5);

      H_BKGD_Overlap_Triggers[ipt] = new TH1D(
      Form("N_DNN%i_Overlap_Triggers_pT%1.0f_%1.0f",2,ptbins[ipt],ptbins[ipt+1]),
      "Number of Isolated Low DNN Photon Triggers", 1, -0.5,3.5);

      Triggers[ipt] = new TH1D(
      Form("N_Triggers_pT%1.0f_%1.0f",ptbins[ipt],ptbins[ipt+1]),
      "Number of Isolated Low DNN Photon Triggers", 2, -0.5,1.0);

      Inclusive_Triggers[ipt] = new TH1D(
      Form("N_Inclusive_Triggers_pT%1.0f_%1.0f",ptbins[ipt],ptbins[ipt+1]),
      "Number of Isolated Low DNN Photon Triggers", 2, -0.5,1.0);

      H_Purities[ipt] = new TH1F(
      Form("H_Purities_pT%1.0f_%1.0f",ptbins[ipt],ptbins[ipt+1]),
      "yields weighted average purity for pT bin", 100, 0.0,1.0);

      H_Purity_Uncertainties[ipt] = new TH1F(
      Form("H_Purity_Uncertanty_pT%1.0f_%1.0f",ptbins[ipt],ptbins[ipt+1]),
      "yields weighted average purity uncertainty for pT bin", 100, 0.0,0.10);
      

      for (int izt = 0; izt<nztbins; izt++){

      Inclusive_Corr[izt+ipt*nztbins] = new TH2D(Form("Inclusive_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",ptbins[ipt],ptbins[ipt+1],
      100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H Isclusive Correlation", n_phi_bins,0,M_PI, n_eta_bins, -1.4, 1.4);


      Corr[izt+ipt*nztbins] = new TH2D(Form("Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",ptbins[ipt],ptbins[ipt+1],
      100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H Isolated Correlation", n_phi_bins,0,M_PI, n_eta_bins, -1.4, 1.4);

      Corr[izt+ipt*nztbins]->Sumw2();

      Corr[izt+ipt*nztbins]->SetMinimum(0.);

      IsoCorr[izt+ipt*nztbins] = new TH2D(Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",1,ptbins[ipt],ptbins[ipt+1],
      100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H Signal Region Correlation", n_phi_bins,0,M_PI,n_eta_bins, -1.4, 1.4);

      IsoCorr[izt+ipt*nztbins]->Sumw2();
      IsoCorr[izt+ipt*nztbins]->SetMinimum(0.);

      BKGD_IsoCorr[izt+ipt*nztbins] = new TH2D(Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",2,ptbins[ipt],ptbins[ipt+1],
      100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H Background Region Correlation", n_phi_bins,0,M_PI, n_eta_bins, -1.4, 1.4);

      BKGD_IsoCorr[izt+ipt*nztbins]->Sumw2();
      BKGD_IsoCorr[izt+ipt*nztbins]->SetMinimum(0.);


      BKGD_IsoCorr_UW[izt+ipt*nztbins] = new TH2D(Form("Unweighted_DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",2,ptbins[ipt],ptbins[ipt+1],
						    100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H [AntiIso] Correlation", n_phi_bins,0,M_PI, n_eta_bins, -1.4, 1.4);

      BKGD_IsoCorr_UW[izt+ipt*nztbins]->Sumw2();
      BKGD_IsoCorr_UW[izt+ipt*nztbins]->SetMinimum(0.);


      h_track_phi_eta[izt+ipt*nztbins] = new TH2D(Form("track_phi_eta__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",2,ptbins[ipt],ptbins[ipt+1],
      100*ztbins[izt],100*ztbins[izt+1]),"Paired Track #phi #eta distribution", n_phi_bins*2,0,M_PI, n_eta_bins*2, -1, 1);

      h_track_phi_eta[izt+ipt*nztbins]->Sumw2();

      h_track_eta[izt+ipt*nztbins] = new TH1D(Form("track_eta__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",2,ptbins[ipt],ptbins[ipt+1],
      100*ztbins[izt],100*ztbins[izt+1]),"Paired Track #eta distribution",n_eta_bins*2, -0.8, 0.8);

      h_track_phi[izt+ipt*nztbins] = new TH1D(Form("track_phi__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",2,ptbins[ipt],ptbins[ipt+1],
      100*ztbins[izt],100*ztbins[izt+1]),"Paired Track #phi distribution", n_phi_bins*2,0,M_PI);

    }//zt bins
  }//pt bins                           
  
    //for (int iarg = 1; iarg < argc; iarg++) {
    int iarg = 1;
    TString root_file = (TString)argv[iarg];
    std::cout << "Opening: " << (TString)argv[iarg] << std::endl;

    TFile *file = TFile::Open(root_file);

    if (file == NULL) {
      std::cout << " fail" << std::endl;
      exit(EXIT_FAILURE);
    }
    file->Print();
    
    TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));

    if (_tree_event == NULL) {
      _tree_event = dynamic_cast<TTree *>(file->Get("AliAnalysisTaskNTGJ/_tree_event"));
      if (_tree_event == NULL) {
	std::cout << " fail " << std::endl;
	exit(EXIT_FAILURE);
      }
    }  

    //Events
    Bool_t is_pileup_from_spd_5_08;
    Double_t primary_vertex[3];
    Float_t ue_estimate_its_const;
    Float_t ue_estimate_tpc_const;
    
    //Tracks
    UInt_t ntrack;
    Float_t track_e[NTRACK_MAX];
    Float_t track_pt[NTRACK_MAX];
    Float_t track_eta[NTRACK_MAX];
    Float_t track_phi[NTRACK_MAX];
    Float_t track_eta_emcal[NTRACK_MAX];
    Float_t track_phi_emcal[NTRACK_MAX];
    UChar_t track_quality[NTRACK_MAX];
    UChar_t track_its_ncluster[NTRACK_MAX];
    Float_t track_its_chi_square[NTRACK_MAX];
    Float_t track_dca_xy[NTRACK_MAX];
    Float_t track_dca_z[NTRACK_MAX];

    //Clusters
    UInt_t ncluster;
    Float_t cluster_e[NTRACK_MAX];
    Float_t cluster_e_max[NTRACK_MAX];
    Float_t cluster_e_cross[NTRACK_MAX];
    Float_t cluster_pt[NTRACK_MAX];
    Float_t cluster_eta[NTRACK_MAX];
    Float_t cluster_phi[NTRACK_MAX]; 
    Float_t cluster_iso_tpc_04[NTRACK_MAX];
    Float_t cluster_iso_its_04[NTRACK_MAX];
    Float_t cluster_frixione_tpc_04_02[NTRACK_MAX];
    Float_t cluster_frixione_its_04_02[NTRACK_MAX];
    Float_t cluster_s_nphoton[NTRACK_MAX][4];
    Int_t cluster_ncell[NTRACK_MAX];
    UShort_t  cluster_cell_id_max[NTRACK_MAX];
    Float_t cluster_lambda_square[NTRACK_MAX][2];   
    Float_t cell_e[17664];
    Float_t cluster_distance_to_bad_channel[NTRACK_MAX];
    UChar_t cluster_nlocal_maxima[NTRACK_MAX];

    Float_t cluster_tof[NTRACK_MAX];
    Float_t cluster_iso_its_04_ue[NTRACK_MAX];
    Float_t cluster_iso_tpc_04_ue[NTRACK_MAX];
    
    //MC
    unsigned int nmc_truth;
    Float_t mc_truth_pt[NTRACK_MAX];
    Float_t mc_truth_eta[NTRACK_MAX];  
    Float_t mc_truth_phi[NTRACK_MAX];
    unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
    unsigned short track_mc_truth_index[NTRACK_MAX];
    UChar_t mc_truth_status[NTRACK_MAX];
    
    //Not using as of Feb 2020
    Float_t mc_truth_first_parent_e[NTRACK_MAX];
    Float_t mc_truth_first_parent_pt[NTRACK_MAX];
    Float_t mc_truth_first_parent_eta[NTRACK_MAX];
    Float_t mc_truth_first_parent_phi[NTRACK_MAX];
    short mc_truth_pdg_code[NTRACK_MAX];
    short mc_truth_first_parent_pdg_code[NTRACK_MAX];
    char mc_truth_charge[NTRACK_MAX];


    //Float_t eg_cross_section;
    //Int_t   eg_ntrial;
     
    // Set the branch addresses of the branches in the TTrees
    _tree_event->SetBranchStatus("*mc*", 0);

    //event Addresses
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
    _tree_event->SetBranchAddress("is_pileup_from_spd_5_08", &is_pileup_from_spd_5_08);
    _tree_event->SetBranchAddress("ue_estimate_its_const", &ue_estimate_its_const);
    _tree_event->SetBranchAddress("ue_estimate_tpc_const", &ue_estimate_tpc_const);
    
    //track Addresses
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
    _tree_event->SetBranchAddress("ntrack", &ntrack);
    _tree_event->SetBranchAddress("track_e", track_e);
    _tree_event->SetBranchAddress("track_pt", track_pt);
    _tree_event->SetBranchAddress("track_eta", track_eta);
    _tree_event->SetBranchAddress("track_phi", track_phi);
    _tree_event->SetBranchAddress("track_eta_emcal", track_eta_emcal);
    _tree_event->SetBranchAddress("track_phi_emcal", track_phi_emcal);
    _tree_event->SetBranchAddress("track_quality", track_quality);
    _tree_event->SetBranchAddress("track_its_ncluster", &track_its_ncluster);
    _tree_event->SetBranchAddress("track_its_chi_square", &track_its_chi_square);
    _tree_event->SetBranchAddress("track_dca_xy", &track_dca_xy);
    _tree_event->SetBranchAddress("track_dca_z", &track_dca_z);

    //Cluster Addresses
    _tree_event->SetBranchAddress("ncluster", &ncluster);
    _tree_event->SetBranchAddress("cluster_e", cluster_e);
    _tree_event->SetBranchAddress("cluster_e_max", cluster_e_max);
    _tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
    _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
    _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
    _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
    _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton);
    _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
    _tree_event->SetBranchAddress("cluster_iso_tpc_04",cluster_iso_tpc_04);
    _tree_event->SetBranchAddress("cluster_iso_its_04",cluster_iso_its_04);
    _tree_event->SetBranchAddress("cluster_frixione_tpc_04_02",cluster_frixione_tpc_04_02);
    _tree_event->SetBranchAddress("cluster_frixione_its_04_02",cluster_frixione_its_04_02);
    _tree_event->SetBranchAddress("cluster_distance_to_bad_channel", cluster_distance_to_bad_channel);
    _tree_event->SetBranchAddress("cluster_nlocal_maxima", cluster_nlocal_maxima);

    _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
    _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
    _tree_event->SetBranchAddress("cell_e", cell_e);

    _tree_event->SetBranchAddress("cluster_tof", cluster_tof);
    _tree_event->SetBranchAddress("cluster_iso_its_04_ue",cluster_iso_its_04_ue);
    _tree_event->SetBranchAddress("cluster_iso_tpc_04_ue",cluster_iso_tpc_04_ue);
    
    //MONTE CARLO STUFF
    _tree_event->SetBranchAddress("mc_truth_pt", mc_truth_pt);
    _tree_event->SetBranchAddress("mc_truth_phi", mc_truth_phi);
    _tree_event->SetBranchAddress("mc_truth_eta", mc_truth_eta);
    _tree_event->SetBranchAddress("mc_truth_status", mc_truth_status);        
    _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
    _tree_event->SetBranchAddress("track_mc_truth_index", track_mc_truth_index);
    
    //_tree_event->SetBranchAddress("eg_cross_section",&eg_cross_section);
    //_tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);


    //IMPORTANT BOOLEAN VARIABLES
    Bool_t Signal = false;
    Bool_t Background = false;
    Bool_t Isolated = false;

    Long64_t nentries = _tree_event->GetEntries();         
    std::cout << " Total Number of entries in TTree: " << nentries << std::endl;

    //Cluster Cut Summary
    fprintf(stderr,"%d: CLUSTER CUT SUMMARY \n ",__LINE__);
    fprintf(stderr,"%d: pT_max =  %f \n ",__LINE__,pT_max);
    fprintf(stderr,"%d: eta max = %f \n ",__LINE__,Eta_max);
    fprintf(stderr,"%d: SR Lambda max = %f \n ",__LINE__,Lambda0_cut);
    fprintf(stderr,"%d: ncell min = %f \n ",__LINE__,Cluster_min);
    fprintf(stderr,"%d: Ecross/E = %f \n ",__LINE__,EcrossoverE_min);
    fprintf(stderr,"%d: Dist. bad channel = %f \n ",__LINE__,Cluster_DtoBad);
    fprintf(stderr,"%d: cluster tof = %f \n ",__LINE__,cluster_time);

    
    //WEIGHTING and CLUSTER SPECTRA LOOP

    fprintf(stderr,"Looping to determine weights and pT spectra \n");
    for(Long64_t ievent = 0; ievent < nentries ; ievent++){     
      //for(Long64_t ievent = 0; ievent < 1000 ; ievent++){
      _tree_event->GetEntry(ievent);
      fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);

      //Event selection
      if(TMath::Abs(primary_vertex[2])>10) continue;
      if(primary_vertex[2]==0.00) continue;

      if(do_pile && is_pileup_from_spd_5_08) continue;
      
      bool first_cluster = true;
      for (ULong64_t n = 0; n < ncluster; n++) {
	if( not(cluster_pt[n]>pT_min and cluster_pt[n]<pT_max)) continue;   //select pt of photons
	if( not(TMath::Abs(cluster_eta[n])<Eta_max)) continue;              //cut edges of detector
	if( not(cluster_ncell[n]>Cluster_min)) continue;                    //removes clusters with 1 or 2 cells
	if( not(cluster_e_cross[n]/cluster_e[n]>EcrossoverE_min)) continue; //removes "spiky" clusters
	if( not(cluster_distance_to_bad_channel[n]>=Cluster_DtoBad)) continue; //removes clusters near bad channels
	if( not(cluster_nlocal_maxima[n] < 3)) continue; //require to have at most 2 local maxima.
	if (not(abs(cluster_tof[n]) < cluster_time)) continue;
	
	float isolation;
	if (determiner == CLUSTER_ISO_TPC_04) isolation = cluster_iso_tpc_04[n];
	else if (determiner == CLUSTER_ISO_ITS_04) isolation = cluster_iso_its_04[n];
	else if (determiner == CLUSTER_ISO_ITS_04_SUB) 
	  isolation = cluster_iso_its_04[n] + cluster_iso_its_04_ue[n] - ue_estimate_its_const*3.1416*0.4*0.4;
	else if (determiner == CLUSTER_ISO_TPC_04_SUB) {
	  isolation = cluster_iso_tpc_04[n] + cluster_iso_tpc_04_ue[n] - ue_estimate_tpc_const*3.1416*0.4*0.4;	  
	  if (Is_pp)
	    isolation = cluster_iso_its_04[n] + cluster_iso_its_04_ue[n]- ue_estimate_its_const*3.1416*0.4*0.4;
	    //pp does not use TPC (17q)
	}
	else if (determiner == CLUSTER_FRIXIONE_TPC_04_02) isolation = cluster_frixione_tpc_04_02[n];
	else isolation = cluster_frixione_its_04_02[n];
	
	Isolated = (isolation<iso_max);

	h_cluster_phi->Fill(cluster_phi[n]);
	h_cluster_eta->Fill(cluster_eta[n]);
	
	if (strcmp(shower_shape.data(),"Lambda")== 0) {

	  Signal = ((cluster_lambda_square[n][0] > 0.1) && (cluster_lambda_square[n][0] < Lambda0_cut));	  
	  //Background =  (cluster_lambda_square[n][0] > Lambda0_cut);
	  //Background =  ((cluster_lambda_square[n][0] > 0.4) && (cluster_lambda_square[n][0] < 1.0)); //DOES NOT WORK!!!!
	  Background =  ((cluster_lambda_square[n][0] > 0.4));
	}
	
	else if (strcmp(shower_shape.data(),"DNN")==0){
	  Signal =  ( (cluster_s_nphoton[n][1] > DNN_min) && (cluster_s_nphoton[n][1]<DNN_max));
	  Background = (cluster_s_nphoton[n][1] > 0.0 && cluster_s_nphoton[n][1] < DNN_Bkgd);

	}

	else if (strcmp(shower_shape.data(),"EMax")==0){

          Signal = (cluster_e_max[n]/cluster_e[n] > Emax_max);
          Background = (cluster_e_max[n]/cluster_e[n] < Emax_min);

        }

	  //High DNN Trigger SIGNAL
	  if (Signal and Isolated){  	    
	    
	    N_Signal_Triggers += 1;
	    Signal_pT_Dist->Fill(cluster_pt[n]);
	    hweight.Fill(cluster_pt[n]);

	    //N_pT_Ranges corresponds to pT binning of purity, not pT binning of correlations
	    for (int i = 0; i < (N_pT_Ranges-1); i++ ){
	      if ((cluster_pt[n] >= pT_Ranges[i]) && (cluster_pt[n] < pT_Ranges[i+1])){
		//fprintf(stderr,"\n%d: purity = %f; pT_Cluster = %f",__LINE__,purities[ipt],cluster_pt[n]);
		h_purity.Fill(purities[i]);
		Cluster_Purity = purities[i];
		Cluster_Purity_Uncertainty = purity_Uncertainties[i];
	      }
	    }

	    for (int ipt = 0; ipt < nptbins; ipt++){
	      if (cluster_pt[n] >ptbins[ipt] && cluster_pt[n] <ptbins[ipt+1]){
		H_Signal_Triggers[ipt]->Fill(1);
		H_Purities[ipt]->Fill(Cluster_Purity); 
		H_Purity_Uncertainties[ipt]->Fill(Cluster_Purity_Uncertainty);
	      }
	    }

	    //fprintf(stderr,"\n %d: Signal Cluster pT = %f \n",__LINE__,cluster_pt[n]);

	  
	  }//Signal

	  //Low DNN Trigger BKGD
	  if (Background and Isolated){

	    N_BKGD_Triggers += 1;
	    hBR.Fill(cluster_pt[n]);
	    
	  //   for (int ipt = 0; ipt < nptbins; ipt++){
	  //     if (cluster_pt[n] >ptbins[ipt] && cluster_pt[n] <ptbins[ipt+1])
	  // 	//H_BKGD_Triggers[ipt]->Fill(1); 
	  //   }
	  } //Background
	
	  //no dnn
	  if (Isolated){
	    //fprintf(stderr,"\n %d: Isolation = %f \n",__LINE__,isolation);
	    for (int ipt = 0; ipt < nptbins; ipt++)
	      Triggers[ipt]->Fill(1);
	  }

	  for (int ipt = 0; ipt < nptbins; ipt++)
	    Inclusive_Triggers[ipt]->Fill(1);

	  //DNN and L0 Signal Overlap
	  if (Isolated and (( (cluster_s_nphoton[n][1] > DNN_min) && (cluster_s_nphoton[n][1]<DNN_max))) 
	      and (((cluster_lambda_square[n][0] > 0.05) && (cluster_lambda_square[n][0] < Lambda0_cut))) ){
	    for (int ipt = 0; ipt < nptbins; ipt++){
	      if (cluster_pt[n] >ptbins[ipt] && cluster_pt[n] <ptbins[ipt+1])
		H_Signal_Overlap_Triggers[ipt]->Fill(1);
	    }
	  }
	  //DNN and L0 Background Overlap
	  if (Isolated and (((cluster_lambda_square[n][0] > 0.4))) 
	      and ((cluster_s_nphoton[n][1] > 0.0 && cluster_s_nphoton[n][1] < DNN_Bkgd)) ){
	    for (int ipt = 0; ipt < nptbins; ipt++){
	      if (cluster_pt[n] >ptbins[ipt] && cluster_pt[n] <ptbins[ipt+1])
	      H_BKGD_Overlap_Triggers[ipt]->Fill(1);		
	    }
	  }
      }//Clusters
    } //Events

    hweight.Divide(&hBR);
    std::cout<<"Clusters Passed Iosalation and Shower Shape: "<<N_Signal_Triggers<<std::endl;
    N_Signal_Triggers = 0; //helps check 2 loops have same cluster criteria

    
    //MAIN CORRELATION LOOP

    fprintf(stderr,"\n Looping for main correlation functions \n");
    for(Long64_t ievent = 0; ievent < nentries ; ievent++){     
      //for(Long64_t ievent = 0; ievent < 10000 ; ievent++){     
      _tree_event->GetEntry(ievent);
      fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);

      Float_t purity_weight = 0;
      Float_t BR_purity_weight = 0;
      bool first_cluster = true;
      //if (not(first_cluster)) continue;

      //Event Selection
      if(TMath::Abs(primary_vertex[2])>10) continue;
      if(primary_vertex[2]==0.00) continue;
      if(do_pile && is_pileup_from_spd_5_08) continue;
    

      for (ULong64_t n = 0; n < ncluster; n++) {
	if( not(cluster_pt[n]>pT_min and cluster_pt[n]<pT_max)) continue;   //select pt of photons
	if( not(TMath::Abs(cluster_eta[n])<Eta_max)) continue;              //cut edges of detector
	if( not(cluster_ncell[n]>Cluster_min)) continue;                    //removes clusters with 1 or 2 cells
	if( not(cluster_e_cross[n]/cluster_e[n]>EcrossoverE_min)) continue; //removes "spiky" clusters
	if( not(cluster_distance_to_bad_channel[n]>=Cluster_DtoBad)) continue; //removes clusters near bad channels
	if( not(cluster_nlocal_maxima[n] < 3)) continue; //require to have at most 2 local maxima.
	if (not(abs(cluster_tof[n]) < cluster_time)) continue;
	
	float isolation;
	if (determiner == CLUSTER_ISO_TPC_04) isolation = cluster_iso_tpc_04[n];
	else if (determiner == CLUSTER_ISO_ITS_04) isolation = cluster_iso_its_04[n];
	else if (determiner == CLUSTER_ISO_ITS_04_SUB) 
	  isolation = cluster_iso_its_04[n] + cluster_iso_its_04_ue[n] - ue_estimate_its_const*3.1416*0.4*0.4;
	else if (determiner == CLUSTER_ISO_TPC_04_SUB) {
	  isolation = cluster_iso_tpc_04[n] + cluster_iso_tpc_04_ue[n] - ue_estimate_tpc_const*3.1416*0.4*0.4;	  
	  if (Is_pp)
	    isolation = cluster_iso_its_04[n] + cluster_iso_its_04_ue[n]- ue_estimate_its_const*3.1416*0.4*0.4;
	    //pp does not use TPC (17q)
	}
	else if (determiner == CLUSTER_FRIXIONE_TPC_04_02) isolation = cluster_frixione_tpc_04_02[n];
	else isolation = cluster_frixione_its_04_02[n];
	
	Isolated = (isolation<iso_max);
	if (strcmp(shower_shape.data(),"Lambda")== 0) {
	  Signal = ((cluster_lambda_square[n][0] > 0.1) and (cluster_lambda_square[n][0] < Lambda0_cut));
	  Background = (cluster_lambda_square[n][0] > Lambda0_cut);
	}
	
	else if (strcmp(shower_shape.data(),"DNN")==0){
	  Signal = ( (cluster_s_nphoton[n][1] > DNN_min) && (cluster_s_nphoton[n][1]<DNN_max));
	  Background = (cluster_s_nphoton[n][1] > 0.0 && cluster_s_nphoton[n][1] < DNN_Bkgd);
	}

	else if (strcmp(shower_shape.data(),"EMax")==0){
          Signal = (cluster_e_max[n]/cluster_e[n] > Emax_max);
          Background = (cluster_e_max[n]/cluster_e[n] < Emax_min);
        }
	
	Bool_t isTrue = false;
	Float_t cluster_truth_pt = 0;
	Float_t cluster_truth_eta = 0;
	Float_t cluster_truth_phi = 0;
	for(int counter = 0 ; counter<32; counter++){
	  unsigned short index = cluster_mc_truth_index[n][counter];                   
	  if(isTrue) break;
	  if(index==65535) continue;
	  if( not (mc_truth_status[index] >0)) continue;        
	  isTrue = true;
	  cluster_truth_pt     = mc_truth_pt[index];
	  cluster_truth_phi    = mc_truth_phi[index];
	  cluster_truth_eta    = mc_truth_eta[index];
	}//end loop over indices
	   
	if (not(isTrue)) continue;
	
	float bkg_weight = 1.0;
	float track_weight = 1.0; //Fake Rate, smearing, efficiency

	if(Background and Isolated){
	  bkg_weight = hweight.GetBinContent(hweight.FindBin(cluster_pt[n]));
	  if (strcmp(shower_shape.data(),"Lambda")!= 0)
	    fprintf(stderr,"%s %f: WARNING \n \n WARNING: Using purity for LAMBDA");

	  BR_purity_weight = (1.0/Get_Purity_ErrFunction(cluster_pt[n],purity_deviation,Is_pp,TPC_Iso_Flag) - 1); //(1-p)/p = 1/p - 1
	  for (int ipt = 0; ipt < nptbins; ipt++){
	    if (cluster_pt[n] >= ptbins[ipt] && cluster_pt[n] < ptbins[ipt+1]){
	      Weights_Sum->Fill((ipt+1),bkg_weight); //Integrate histo for pTbin for sum of weights
	      weight_sum[ipt] += bkg_weight;
	      H_BKGD_Triggers[ipt]->Fill(1,bkg_weight); //Fill Single Bin, Sum = Bin Contents
	      BR_purity_weight_sum[ipt] += BR_purity_weight;
	    }
	  }
	  //fprintf(stderr,"\n %d: weight = %f \n",__LINE__,bkg_weight);
	  BKGD_pT_Dist->Fill(cluster_pt[n]);
	  BKGD_pT_Dist_Weighted->Fill(cluster_pt[n],bkg_weight);
	  }

	if (Signal and Isolated){
	  N_Signal_Triggers += 1;
	  purity_weight = 1.0/Get_Purity_ErrFunction(cluster_pt[n],purity_deviation,Is_pp,TPC_Iso_Flag);
	  for (int ipt = 0; ipt < nptbins; ipt++){
	    if (cluster_pt[n] >= ptbins[ipt] && cluster_pt[n] < ptbins[ipt+1]){
	    purity_weight_sum[ipt] += purity_weight;
	    }
	  }
	  //fprintf(stderr,"\n%d: Purity From Fit = %f\n",__LINE__,Get_Purity_ErrFunction(cluster_pt[n],purity_deviation));
	}

	//Track Loop

	for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {            
 	  if(track_pt[itrack] < track_pT_min) continue; //500 MeV Tracks or 1GeV Tracks
 	  if(track_pt[itrack] > track_pT_max) continue;
 	  if((track_quality[itrack]&Track_Cut_Bit)==0) continue; //select only tracks that pass selection 
	  if(abs(track_eta[itrack]) > 0.8) continue;
	  if( not(track_its_ncluster[itrack]>4)) continue;
	  if( not(track_its_chi_square[itrack]/track_its_ncluster[itrack] < track_chi_max)) continue;
	  //if( not(track_its_chi_square[itrack]/track_its_ncluster[itrack] <36)) continue;
	  //if( not(TMath::Abs(track_dca_xy[itrack])<0.0231+0.0315/TMath::Power(track_pt[itrack],1.3 ))) continue;
	  if (not(TMath::Abs(track_dca_xy[itrack]<2.4))) continue;
	  if (not(TMath::Abs(track_dca_z[n]) < 3.2)) continue;

	  unsigned int track_mc_index = track_mc_truth_index[itrack];
	  if (track_mc_index < 65534) continue;

	  //Electron Veto for associated tracks outside of isolation cone
	  double dRmin = 0.02;
	  bool Track_HasMatch = false;
	  for (ULong64_t c = 0; c < ncluster; c++){
	    Float_t deta =  cluster_eta[c]-track_eta_emcal[itrack];
	    Float_t dphi =  TVector2::Phi_mpi_pi(cluster_phi[c]-track_phi_emcal[itrack])/TMath::Pi();
	    float dR = sqrt(dphi*dphi + deta*deta);
	    if (dR < dRmin) {
	      Track_HasMatch = true;
	      break;
	    }
	  }
 	  //if (Track_HasMatch) continue;

	  //Apply corrections as weights. 1/Eff, *FakeRate, *SmearingAffect

	  if (first_cluster){
		    Track_Weight_Spectra.Fill(track_pt[itrack],track_weight);
		    Track_Raw_Spectra.Fill(track_pt[itrack]);
		    Chi_Distro.Fill(track_its_chi_square[itrack]/track_its_ncluster[itrack]);
	  }
	  
	  //Observables:
	  Double_t zt = track_pt[itrack]/cluster_pt[n];
	  Float_t DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(cluster_truth_phi - mc_truth_phi[track_mc_index]));
	  Float_t DeltaEta = cluster_truth_eta - mc_truth_eta[track_mc_index];
	  if ((TMath::Abs(DeltaPhi) < 0.005) && (TMath::Abs(DeltaEta) < 0.005)) continue; //Match Mixing Cut

	  for (int ipt = 0; ipt < nptbins; ipt++){
	    if (cluster_pt[n] >= ptbins[ipt] && cluster_pt[n] < ptbins[ipt+1]){
	      for(int izt = 0; izt<nztbins ; izt++){
		if(zt>ztbins[izt] and  zt<ztbins[izt+1]){
		  if (first_cluster){
		    h_track_phi_eta[izt+ipt*nztbins]->Fill(mc_truth_phi[track_mc_index],mc_truth_eta[track_mc_index]);
		    if (track_pt[itrack] > 10)
		    tracks_tenGev += 1;
		  }
		  //2 DNN Regions

		  if (Signal and Isolated)
		    //IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta,track_weight*purity_weight);
		    IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta,purity_weight);

		  if (Background and Isolated){
		    BKGD_IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
		    //not weighted with pT distro
		    BKGD_IsoCorr_UW[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
		  }
		
	    	  
		  //no shower shape selection
		  if(Isolated)
		    Corr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);

		  Inclusive_Corr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);	       
		  
		  if(track_pt[itrack] < 0.85)
		    Track_Uncertainty[izt]->Fill(0.08);
		  else
		    Track_Uncertainty[izt]->Fill(0.05);
		  
		}//if in zt bin 
	      } // for zt bins
	    }//if in pt bin
	  }//for pt bins
	}//for itracks
	first_cluster = false;
      }//for nclusters
    } //for nevents  
    //}//end loop over samples

  // Write to fout
  
  //TFile* fout = new TFile(Form("fout_Corr_config%s.root", opened_files.c_str()),"RECREATE");

  size_t lastindex = std::string(root_file).find_last_of(".");
  std::string rawname = std::string(root_file).substr(0, lastindex);
  fprintf(stderr,"%s: %d: Creating new file %s_SE_L0_Correlation.root",__FILE__,__LINE__,rawname.data());

  TFile* fout;

  if (strcmp(shower_shape.data(),"Lambda")== 0) 
    fout = new TFile(Form("%s_SE_L0_Correlation.root",rawname.data()),"RECREATE");
  else if (strcmp(shower_shape.data(),"DNN")== 0)
    fout = new TFile(Form("%s_SE_NN_Correlation.root",rawname.data()),"RECREATE");
  else if (strcmp(shower_shape.data(),"EMax")== 0)
    fout = new TFile(Form("%s_SE_EMax_Correlation.root",rawname.data()),"RECREATE");
  else
    fout = new TFile(Form("%s_SE_Correlation.root",rawname.data()),"RECREATE");

  std::cout<<"Clusters Passed Iosalation and Shower Shape: "<<N_Signal_Triggers<<std::endl;
  fprintf(stderr,"Tracks Above 10 GeV/c = %llu \n",tracks_tenGev);
			  
  h_purity.Write("purities");

  h_cluster_phi->Write();
  h_cluster_eta->Write();

  hweight.Write();

  Track_Weight_Spectra.Write();
  Track_Raw_Spectra.Write();
  Chi_Distro.Write();
  Signal_pT_Dist->Write();
  BKGD_pT_Dist->Write();
  BKGD_pT_Dist_Weighted->Scale(1.0/(hweight.Integral(1,40))); //Divide by sum of weights
  BKGD_pT_Dist_Weighted->Write();

  for (int ipt = 0; ipt<nptbins; ipt++)
    H_Signal_Triggers[ipt]->Write();
  
  for (int ipt = 0; ipt<nptbins; ipt++)
    H_BKGD_Triggers[ipt]->Write(); //NTriggers is filled with weights. Later I divide by this -> divide by sum of weights
    
  for (int ipt = 0; ipt < nptbins; ipt++)
    Triggers[ipt]->Write();

  for (int ipt = 0; ipt < nptbins; ipt++)
    Inclusive_Triggers[ipt]->Write();

  for (int ipt = 0; ipt < nptbins; ipt ++)
    H_Signal_Overlap_Triggers[ipt]->Write();

  for (int ipt = 0; ipt < nptbins; ipt ++)
    H_BKGD_Overlap_Triggers[ipt]->Write();


  for (int ipt = 0; ipt < nptbins; ipt++)
    H_Purities[ipt]->Write();

  for (int ipt = 0; ipt < nptbins; ipt++)
    H_Purity_Uncertainties[ipt]->Write();


  for (int ipt = 0; ipt<nptbins; ipt++){
    //Weights_Sum->SetBinContent(ipt+1,weight_sum[ipt]);
    fprintf(stderr,"\n %d: Weight sum in pt bin %i: %f \n",__LINE__,ipt,weight_sum[ipt]);    
    fprintf(stderr,"\n %d: Weight sum in pt bin %i: %f \n",__LINE__,ipt,Weights_Sum->GetBinContent(ipt+1));    
    fprintf(stderr,"\n %d: Weight sum in pt bin %i: %f \n",__LINE__,ipt,H_BKGD_Triggers[ipt]->GetBinContent(1));    
    
    for (int izt = 0; izt<nztbins; izt++)
      h_track_phi_eta[izt+ipt*nztbins]->Write();

    for (int izt = 0; izt<nztbins; izt++){
      Corr[izt+ipt*nztbins]->Write();
    }
    for (int izt = 0; izt<nztbins; izt++){
      IsoCorr[izt+ipt*nztbins]->Write();
    }
    for (int izt = 0; izt<nztbins; izt++){
      BKGD_IsoCorr[izt+ipt*nztbins]->Write(); //Not divided by sum of weights here. py notebooks divide by # entries -> sum of weights. See Bkd trigger histo
    }
    for (int izt = 0; izt<nztbins; izt++){
      BKGD_IsoCorr_UW[izt+ipt*nztbins]->Write();
    }
    for (int izt = 0; izt<nztbins; izt++){
      Inclusive_Corr[izt+ipt*nztbins]->Write();
    }
  }
  
  for(int izt = 0; izt < nztbins; izt++){
    Track_Uncertainty[izt]->Write();
  }
  Weights_Sum->Write();

  //Seperate zt loops for easier file reading
  fout->Close();     
  file->Close();  
  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}
