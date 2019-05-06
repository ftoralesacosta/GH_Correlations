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

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};

int main(int argc, char *argv[])
{
  if (argc < 2) {
    exit(EXIT_FAILURE);
  }
  int dummyc = 1;
  char **dummyv = new char *[1];

  dummyv[0] = strdup("main");
  
  //Config File
  FILE* config = fopen("Corr_config.yaml", "r");
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
  int Track_Cut_Bit = 0;
  double iso_max = 0;
  double noniso_min = 0;
  double noniso_max = 0;
  double deta_max = 0;
  isolationDet determiner = CLUSTER_ISO_ITS_04;
  int n_eta_bins = 0;
  int n_phi_bins = 0;  
  std::string shower_shape = "DNN";

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

      else if (strcmp(key, "iso_max") == 0) {
          iso_max = atof(value);
          std::cout << "iso_max: " << iso_max << std::endl; }

      else if (strcmp(key, "noniso_min") == 0) {
          noniso_min = atof(value);
          std::cout << "noniso_min: " << noniso_min << std::endl; }

      else if (strcmp(key, "noniso_max") == 0) {
          noniso_max = atof(value);
          std::cout << "noniso_max: " << noniso_max << std::endl; }

      else if (strcmp(key, "deta_max") == 0) {
          deta_max = atof(value);
          std::cout << "deta_max: " << deta_max << std::endl; }

      else if (strcmp(key, "N_Phi_Bins") == 0) {
	n_phi_bins = atoi(value);
	std::cout << "Number of Phi Bins: " << n_phi_bins << std::endl; }

      else if (strcmp(key, "N_Eta_Bins") == 0) {
        n_eta_bins = atoi(value);
	std::cout << "Number of Eta Bins: " << n_eta_bins << std::endl; }

      else if (strcmp(key, "Track_Cut_Bit") == 0) {
          Track_Cut_Bit = atoi(value);
          std::cout << "Track Cut Bit: " << Track_Cut_Bit << std::endl; }

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
  int  N_pT_Ranges = 9;
  float pT_Ranges[9] = {12,13.46,15.09,16.92,18.97,21.28,23.86,30,40};
  float purities[8] = {0};
  float Cluster_Purity = 0;

  if (strcmp(shower_shape.data(), "DNN") == 0){
    purities[0] = 0.229;
    purities[1] = 0.255;
    purities[2] = 0.326; 
    purities[3] = 0.372;
    purities[4] = 0.447;
    purities[5] = 0.502;
    purities[6] = 0.533;
    purities[7] = 0.531;
  }

  else if (strcmp(shower_shape.data(), "Lambda") == 0){
    purities[0] = 0.207;
    purities[1] = 0.210;
    purities[2] = 0.328;
    purities[3] = 0.367;
    purities[4] = 0.447;
    purities[5] = 0.524;
    purities[6] = 0.509;
    purities[7] = 0.597;
  }

  else {
    purities[0]=0.27;purities[1]=0.32;purities[2]=0.37;}//Emax/Ecluster

  //Track Corrections (Calculated in 1GeV pT bins, used as weights when filling)

  float Smearing_Correction[15] = {1.007,1.007,0.982,0.957,0.926,0.894,0.853,0.817,0.757,0.681,0.673,0.619,0.469,0.342,0.301};
  
  float OneMinFakeRate[15] = {0.9821,0.9821,0.9803,0.9751,0.9645,0.9525,0.9278,0.9098,0.8702,0.8593,0.7870,0.7825,0.7624,0.7389,0.6710};

  const float Efficiency = 0.85;


  TH2D* Corr[nztbins*nptbins];
  TH2D* IsoCorr[nztbins*nptbins];
  TH2D* BKGD_IsoCorr[nztbins*nptbins];
  TH2D* BKGD_IsoCorr_UW[nztbins*nptbins];
  TH2D* Inclusive_Corr[nztbins*nptbins];

  //Dihadron Stuff                                                                                                                                            
  int n_track_ptbins = 14;
  float track_ptbins[15] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  TH2D* dihadron_Corr[n_track_ptbins*nztbins];
  TH2D * All_Dihadron_Corr = new TH2D("All_Dihadron_Corr","H-H All Correlation" ,n_phi_bins,0,M_PI,n_eta_bins, -1.4, 1.4);

  for (int ipt = 0 ; ipt < n_track_ptbins; ipt++){
    for (int izt = 0; izt < nztbins; izt++){
      dihadron_Corr[izt+ipt*nztbins] = new TH2D(
      Form("DiHadron_Correlation__pT%1.0f_%1.0f__zT_%1.0f_%1.0f",track_ptbins[ipt],track_ptbins[ipt+1],100*ztbins[izt],100*ztbins[izt+1]),
      Form("H-H Correlation Track pT %1.0f - %1.0f GeV, zT %1.0f - %1.0f",track_ptbins[ipt]*10,track_ptbins[ipt+1]*10,100*ztbins[izt],100*ztbins[izt+1]),
      n_phi_bins,0,M_PI, n_eta_bins, -1.4, 1.4);
    }
  }


  TH1D* H_Signal_Triggers[nptbins];
  TH1D* H_BKGD_Triggers[nptbins];
  TH1D* H_Signal_Overlap_Triggers[nptbins];
  TH1D* H_BKGD_Overlap_Triggers[nptbins];
  TH1D* Triggers[nptbins];

  TH1F* H_Purities[nptbins];

  TH2D * h_track_phi_eta[nztbins*nptbins];
  TH1D * h_track_eta[nztbins*nptbins];
  TH1D * h_track_phi[nztbins*nptbins];

  float N_Signal_Triggers = 0;
  float N_BKGD_Triggers = 0;
  
  TH1F* Signal_pT_Dist = new TH1F("Signal_pT_Dist","Cluster Pt Spectrum For Isolation (its_04) bins 0.55 < DNN < 0.85",(pT_max-pT_min)*2,pT_min,pT_max);
  TH1F* BKGD_pT_Dist = new TH1F("BKGD_pT_Dist","Cluster Pt Spectrum For Isolation (its_04) bins 0.0 < DNN < 0.3",(pT_max-pT_min)*2,pT_min,pT_max);
  TH1F* BKGD_pT_Dist_Weighted = new TH1F("BKGD_pT_Dist_Weighted","Weighted Cluster Pt Spectrum For Isolation (its_04) bins 0.0 < DNN < 0.3",(pT_max-pT_min)*2,pT_min,pT_max);

  Signal_pT_Dist->Sumw2();
  BKGD_pT_Dist->Sumw2();
  BKGD_pT_Dist_Weighted->Sumw2();

  TH1F hBR("hBR", "Isolated cluster, bkg region", 40, 10.0, 50.0);
  TH1F hweight("hweight", "Isolated cluster, signal region", 40, 10.0, 50.0);
  TH1F* Weights_Sum = new TH1F("Weights_Sum", "Sum of Weights. XBin Centers = pt Bin", nptbins,0.5,4.5);
  float weight_sum[nptbins];

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

      H_Purities[ipt] = new TH1F(
      Form("H_Purities_pT%1.0f_%1.0f",ptbins[ipt],ptbins[ipt+1]),
      "yields weighted average purity for pT bin", 100, 0.0,1.0);
      
      weight_sum[ipt] = 0;

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

    //Tracks
    Double_t primary_vertex[3];
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
    unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
    Int_t cluster_ncell[NTRACK_MAX];
    UShort_t  cluster_cell_id_max[NTRACK_MAX];
    Float_t cluster_lambda_square[NTRACK_MAX][2];   
    Float_t cell_e[17664];
    Float_t cluster_distance_to_bad_channel[NTRACK_MAX];
    UChar_t cluster_nlocal_maxima[NTRACK_MAX];

    //MC
    unsigned int nmc_truth;
    Float_t mc_truth_pt[NTRACK_MAX];
    Float_t mc_truth_eta[NTRACK_MAX];  
    Float_t mc_truth_phi[NTRACK_MAX];
    short mc_truth_pdg_code[NTRACK_MAX];
    short mc_truth_first_parent_pdg_code[NTRACK_MAX];
    char mc_truth_charge[NTRACK_MAX];

    Float_t mc_truth_first_parent_e[NTRACK_MAX];
    Float_t mc_truth_first_parent_pt[NTRACK_MAX];
    Float_t mc_truth_first_parent_eta[NTRACK_MAX];
    Float_t mc_truth_first_parent_phi[NTRACK_MAX];
    UChar_t mc_truth_status[NTRACK_MAX];
    //Float_t eg_cross_section;
    //Int_t   eg_ntrial;
     
    // Set the branch addresses of the branches in the TTrees
    _tree_event->SetBranchStatus("*mc*", 0);
  
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

    //Cluster Addresses
    _tree_event->SetBranchAddress("ncluster", &ncluster);
    _tree_event->SetBranchAddress("cluster_e", cluster_e);
    _tree_event->SetBranchAddress("cluster_e_max", cluster_e_max);
    _tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
    _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
    _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
    _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
    _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton);
    _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
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

    //_tree_event->SetBranchAddress("eg_cross_section",&eg_cross_section);
    //_tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);


    //IMPORTANT BOOLEAN VARIABLES
    Bool_t Signal = false;
    Bool_t Background = false;
    Bool_t Isolated = false;

    Long64_t nentries = _tree_event->GetEntries();         
    //Long64_t nentries = 200;
    std::cout << " Total Number of entries in TTree: " << nentries << std::endl;

    //MAIN CORRELATION LOOP

    fprintf(stderr,"\n Looping for main correlation functions \n");
    for(Long64_t ievent = 0; ievent < nentries ; ievent++){     
      //for(Long64_t ievent = 0; ievent < 1000 ; ievent++){
      _tree_event->GetEntry(ievent);
      fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);

      bool first_track = true;
      float track_weight = 1.0;

      for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {            
	if(track_pt[itrack] < 1.0) continue; //500 MeV Tracks
	if(track_pt[itrack] > 15) continue;
	if((track_quality[itrack]&Track_Cut_Bit)==0) continue; //select only tracks that pass selection 
	if(abs(track_eta[itrack]) > 0.8) continue;
	if( not(track_its_ncluster[itrack]>4)) continue;
	if( not(track_its_chi_square[itrack]/track_its_ncluster[itrack] <36)) continue;
	if( not(TMath::Abs(track_dca_xy[itrack])<0.0231+0.0315/TMath::Power(track_pt[itrack],1.3 ))) continue;
	
	//Track Loop
	for (ULong64_t jtrack = 0; jtrack < ntrack; jtrack++) {            
 	  if(track_pt[jtrack] < 0.5) continue; //500 MeV Tracks
 	  if(track_pt[jtrack] > 15) continue;
 	  if((track_quality[jtrack]&Track_Cut_Bit)==0) continue; //select only tracks that pass selection 
	  if(abs(track_eta[jtrack]) > 0.8) continue;
	  if( not(track_its_ncluster[jtrack]>4)) continue;
	  if( not(track_its_chi_square[jtrack]/track_its_ncluster[jtrack] <36)) continue;
	  if( not(TMath::Abs(track_dca_xy[jtrack])<0.0231+0.0315/TMath::Power(track_pt[jtrack],1.3 ))) continue;

	  //Electron Veto for associated tracks outside of isolation cone
	  double dRmin = 0.02;
	  bool Track_HasMatch = false;
	  for (ULong64_t c = 0; c < ncluster; c++){
	    Float_t deta =  track_eta[itrack]-track_eta_emcal[jtrack];
	    Float_t dphi =  TVector2::Phi_mpi_pi(track_phi[itrack]-track_phi_emcal[jtrack])/TMath::Pi();
	    float dR = sqrt(dphi*dphi + deta*deta);
	    if (dR < dRmin) {
	      Track_HasMatch = true;
	      break;
	    }
	  }
 	  //if (Track_HasMatch) continue;

	  //Apply corrections as weights. 1/Eff, *FakeRate, *SmearingAffect

	  for (int ipt = 0; ipt < 15; ipt++){
	    if ( (track_pt[jtrack] >= ipt) && (track_pt[jtrack] < (ipt+1)) ) {
	      track_weight = Smearing_Correction[ipt]*OneMinFakeRate[ipt]/Efficiency;
	    }
	  }
	
	  //Observables:
	  Double_t zt = track_pt[jtrack]/track_pt[itrack];
	  Float_t DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(track_phi[itrack] - track_phi[jtrack]));
	  Float_t DeltaEta = track_eta[itrack] - track_eta[jtrack];
	  if ((TMath::Abs(DeltaPhi) < 0.005) && (TMath::Abs(DeltaEta) < 0.005)) continue; //Match Mixing Cut

	  All_Dihadron_Corr->Fill(DeltaPhi,DeltaEta);
	  
	  //fprintf(stderr,"\n %s %d: %f \n",__FILE__,__LINE__,track_pt[itrack]);
	  for (int ipt = 0; ipt < n_track_ptbins; ipt++){
	    if (track_pt[itrack] > track_ptbins[ipt] && track_pt[itrack] < track_ptbins[ipt+1]){
	      for (int izt =0; izt < nztbins; izt++){
		if(zt>ztbins[izt] and  zt<ztbins[izt+1]){
		  dihadron_Corr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
		  // if (first_track)
		  //   h_track_phi_eta[izt+ipt*nztbins]->Fill(track_phi[jtrack],track_eta[jtrack]);	  
		}//if in zt bin 
	      } // for zt bins
	    }//if in pt bin
	  }//for pt bins
	}//for jtracks
	first_track = false;
      }//for nclusters
    } //for nevents  
//}//end loop over samples

  // Write to fout
  
  //TFile* fout = new TFile(Form("fout_Corr_config%s.root", opened_files.c_str()),"RECREATE");

  size_t lastindex = std::string(root_file).find_last_of(".");
  std::string rawname = std::string(root_file).substr(0, lastindex);
  fprintf(stderr,"%s: %d: Creating new file",__FILE__,__LINE__);

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
			  
  h_purity.Write("purities");

  h_cluster_phi->Write();
  h_cluster_eta->Write();

  hweight.Write();

  Signal_pT_Dist->Write();
  BKGD_pT_Dist->Scale(1.0/(hweight.Integral(1,40))); //Divide by sum of weights
  BKGD_pT_Dist->Write();
  BKGD_pT_Dist_Weighted->Write();

  for (int ipt = 0; ipt<nptbins; ipt++)
    H_Signal_Triggers[ipt]->Write();
  
  for (int ipt = 0; ipt<nptbins; ipt++)
    H_BKGD_Triggers[ipt]->Write();
    
  for (int ipt = 0; ipt < nptbins; ipt++)
    Triggers[ipt]->Write();

  for (int ipt = 0; ipt < nptbins; ipt ++)
    H_Signal_Overlap_Triggers[ipt]->Write();

  for (int ipt = 0; ipt < nptbins; ipt ++)
    H_BKGD_Overlap_Triggers[ipt]->Write();


  for (int ipt = 0; ipt < nptbins; ipt++)
    H_Purities[ipt]->Write();

  All_Dihadron_Corr->Write();
  for (int ipt = 0; ipt < n_track_ptbins-1; ipt++){
    for (int izt = 0; izt<nztbins; izt++){
      dihadron_Corr[izt+ipt*nztbins]->Write();
    }
  }

  for (int ipt = 0; ipt<nptbins; ipt++){
    //Weights_Sum->SetBinContent(ipt+1,weight_sum[ipt]);
    fprintf(stderr,"\n %d: Weight sum in pt bin %i: %f \n",__LINE__,ipt,weight_sum[ipt]);    
    fprintf(stderr,"\n %d: Weight sum in pt bin %i: %f \n",__LINE__,ipt,Weights_Sum->GetBinContent(ipt+1));    
    fprintf(stderr,"\n %d: Weight sum in pt bin %i: %f \n",__LINE__,ipt,H_BKGD_Triggers[ipt]->GetBinContent(1));    
    
  //   for (int izt = 0; izt<nztbins; izt++)
  //     h_track_phi_eta[izt+ipt*nztbins]->Write();

  //   for (int izt = 0; izt<nztbins; izt++){
  //     Corr[izt+ipt*nztbins]->Write();
  //   }
  //   for (int izt = 0; izt<nztbins; izt++){
  //     IsoCorr[izt+ipt*nztbins]->Write();
  //   }
  //   for (int izt = 0; izt<nztbins; izt++){
  //     BKGD_IsoCorr[izt+ipt*nztbins]->Write();
  //   }
  //   for (int izt = 0; izt<nztbins; izt++){
  //     BKGD_IsoCorr_UW[izt+ipt*nztbins]->Write();
  //   }
  //   for (int izt = 0; izt<nztbins; izt++){
  //     Inclusive_Corr[izt+ipt*nztbins]->Write();
  //   }
  }

  // Weights_Sum->Write();
  
  //Seperate zt loops for easier file reading
  fout->Close();     
  file->Close();  
  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}
