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


  // Loop through config file
  void Read_Config() {
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

  return;
  }
