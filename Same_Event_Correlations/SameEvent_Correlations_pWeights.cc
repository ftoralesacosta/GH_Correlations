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
  if (argc < 2) {
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



  std::array< float,38 > track_pT_Correction_bins = {0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
						0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90,
						0.95, 1.00, 1.10, 1.20, 1.40, 1.60, 1.80, 2.00,
						2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.60, 4.00,
						5.00, 6.00, 8.00, 10.0, 13.0, 20.0};

  int N_Track_pT_Bins = track_pT_Correction_bins.size()-1;
  fprintf(stderr,"NUMBER OF TRACK PT BINS IS %i \n \n",N_Track_pT_Bins);
    

  //p-Pb
  float pPb_Efficiency[37] = {0.701584, 0.794149, 0.811056, 0.824281, 0.833784, 0.8411, 0.846235, 0.850521,
			      0.853948, 0.855174, 0.857245, 0.858315, 0.858334, 0.859078, 0.860819, 0.860223,
			      0.860368, 0.860849, 0.862532, 0.864867, 0.867081, 0.869168, 0.869016, 0.872149,
			      0.872669, 0.874919, 0.876616, 0.876491, 0.87365, 0.877567, 0.8777, 0.882617,
			      0.88169, 0.881348, 0.880515, 0.87605, 0.882653};

  float pPb_FakeRate[37] = {0.0444643, 0.0329245, 0.0274773, 0.0241898, 0.0222288, 0.0207158, 0.0195863,
			    0.0189526, 0.0184162, 0.0180536, 0.0175808, 0.0177315, 0.017413, 0.0172551,
			    0.0175637, 0.017715, 0.0172888, 0.0178326, 0.0178794, 0.0179673, 0.0177509,
			    0.0182945, 0.0181996, 0.0188647, 0.0190429, 0.0202744, 0.0209381, 0.0213728,
			    0.0237848, 0.0251898, 0.0297774, 0.0350703, 0.0518793, 0.0804111, 0.137953,
			    0.214123, 0.301459};

  float pPb_Smearing_Correction[37] = {0.894311, 0.994658, 1.00122, 1.0087, 1.00914, 1.01009, 1.00939, 1.01073, 1.00998, 1.00622, 1.00519, 1.00425, 1.00182, 0.99553, 0.992115, 0.991085, 0.989869, 0.988516, 0.987711, 0.989906, 0.990167, 0.994517, 0.994247, 0.99273, 0.983717, 0.986158, 0.978648, 0.986939, 0.963266, 0.959108, 0.945863, 0.923434, 0.896206, 0.86263, 0.71868, 0.629909, 0.425061};
  


    std::array< float,51 > TPC_track_pT_Correction_bins = {0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
							   0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
							   0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
							   1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8,
							   3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5,
							   5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14};

    int N_TPC_Track_pT_Bins = TPC_track_pT_Correction_bins.size()-1;

    float TPC_pPb_Efficiency[N_TPC_Track_pT_Bins] = {0.5867693424, 0.6911254525,0.7336804867, 0.7650601268,
						     0.7925264239, 0.8127673268,0.8275696635, 0.8387191296,
						     0.8461661935, 0.8519890308,0.85583359, 0.8587060571,
						     0.8605732322, 0.8621085882,0.8631671667, 0.8642088175,
						     0.8647749424, 0.8643072844,0.8677230477, 0.8718497753,
						     0.8744822741, 0.8703474402,0.8715935946, 0.8735370636,
						     0.8726101518, 0.8705785275,0.8655129671, 0.8545784354,
						     0.8449353576, 0.8368939161,0.8319132924, 0.8291873336,
						     0.8325293064, 0.834860146,0.8377019167, 0.8426615,
						     0.8430748582, 0.8475763798,0.856606245, 0.8528164625,
						     0.8603565693, 0.8575622439,0.8744791746, 0.8765413761,
						     0.8653901219, 0.8669821024,0.8549746871, 0.8416666389,
						     0.8819277287, 0.8582089543};

    float TPC_pPb_FakeRate[N_TPC_Track_pT_Bins] = {0.02952835709, 0.02360016666,
					     0.02065094374, 0.01941830292,0.01874902286, 0.01761478558,
					     0.01669583283, 0.01596949063,0.01551650651, 0.01526045054,
					     0.01498838048, 0.01479246933,0.01460200734, 0.01449344773,
					     0.01432006247, 0.0139457956,0.01425391622, 0.01395579707,
					     0.0137891043, 0.01371407881,0.013224978, 0.01300536282,
					     0.01288809255, 0.01266333833,0.01300011016, 0.01260772068,
					     0.01240321435, 0.01243534964,0.01235494018, 0.01211020909,
					     0.01210790593, 0.01196851861,0.01157244854, 0.01151115727,
					     0.0119260475, 0.01144340821,0.01208098885, 0.0112716658,
					     0.01107564662, 0.01196976099,0.01059919689, 0.01005556993,
					     0.01106236875, 0.01042175572,0.009913654998, 0.01078320108,
					     0.01424501464, 0.01892744564,0.008196720853, 0.007968127728};

    float TPC_pPb_Smearing_Correction[N_TPC_Track_pT_Bins] =  {0.9633022736,1.000103796,0.9990538363,0.9980777293,
							       0.9992338076,1.00126975,1.002764827,1.003955855,
							       1.00396996,1.004819844,1.004234793,1.003327841,
							       1.004936572,1.004297549,1.004503464,1.004763271,
							       1.006574166,1.003742764,1.003480562,1.003139268,
							       1.00315121,1.003883008,1.003246503,1.002383757,
							       1.003722925,1.002248603,1.003694197,1.003056597,
							       1.004307807,1.003282579,1.005886125,1.0035825,
							       1.001216322,1.00245251,0.9995267124,0.9980988593,
							       1.000157233,0.9919125185,0.9941652735,0.9886572324,
							       0.9887442842,0.9957731115,0.9868730408,0.9894348482,
							       0.997693575,0.9665885111,0.9941176471,0.9983525535,
							       1.045714286,0.9426229508};

    
  //pp

  std::array< float,38 > pp_track_pT_Correction_bins = {0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
                                                0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90,
                                                0.95, 1.00, 1.10, 1.20, 1.40, 1.60, 1.80, 2.00,
                                                2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.60, 4.00,
                                                5.00, 6.00, 8.00, 10.0, 13.0, 20.0};

  float pp_Efficiency[37] = {0.686749, 0.774554, 0.791098, 0.806701, 0.815491, 0.824043, 0.83069, 0.833836,
			     0.836869, 0.839235, 0.841083, 0.843128, 0.841643, 0.844369, 0.843158, 0.842828,
			     0.844946, 0.846174, 0.849676, 0.85162, 0.851656, 0.855973, 0.859772, 0.856858,
			     0.85658, 0.866535, 0.862416, 0.864511, 0.863729, 0.871034, 0.870704, 0.872091,
			     0.86718, 0.848246, 0.86087, 0.887701, 0.863014};

  float pp_FakeRate[37] = {0.0500267, 0.0398189, 0.0350848, 0.0316228, 0.0292172, 0.0279222, 0.0269255,
			   0.0262388, 0.0250246, 0.0245797, 0.0242376, 0.0240682, 0.023281, 0.0246955,
			   0.0235495, 0.0230017, 0.0236568, 0.022912, 0.0237135, 0.0229813, 0.0223819,
			   0.0223993, 0.0238591, 0.0252843, 0.0253992, 0.0250655, 0.0254672, 0.0271102,
			   0.0272849, 0.0233524, 0.0321137, 0.040308, 0.0611087, 0.0730563, 0.13447, 0.148867, 0.292818};


  float pp_Smearing_Correction[37] = {0.899306, 0.990662, 0.999535, 1.0065, 1.00578, 1.01049, 1.0104, 1.00846, 1.00905, 1.00653, 1.00859, 1.00616, 0.997317, 1.00387, 0.993037, 1.00025, 0.994664, 0.99143, 1.00412, 0.994017, 0.989386, 0.995787, 0.987843, 0.990099, 0.972225, 0.997532, 0.975834, 0.971739, 0.939599, 0.951686, 0.971399, 0.904032, 0.915311, 0.897727, 0.670429, 0.666667, 0.512195};
    
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

    Float_t cluster_tof[NTRACK_MAX];
    Float_t cluster_iso_its_04_ue[NTRACK_MAX];
    Float_t cluster_iso_tpc_04_ue[NTRACK_MAX];
    
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

    _tree_event->SetBranchAddress("cluster_tof", cluster_tof);
    _tree_event->SetBranchAddress("cluster_iso_its_04_ue",cluster_iso_its_04_ue);
    _tree_event->SetBranchAddress("cluster_iso_tpc_04_ue",cluster_iso_tpc_04_ue);
    
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
	  if( not(track_its_chi_square[itrack]/track_its_ncluster[itrack] <36)) continue;
	  if( not(TMath::Abs(track_dca_xy[itrack])<0.0231+0.0315/TMath::Power(track_pt[itrack],1.3 ))) continue;

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

	    for (int ipt = 0; ipt < N_Track_pT_Bins; ipt++){
	      if( (track_pt[itrack] >= track_pT_Correction_bins[ipt]) && (track_pt[itrack] < track_pT_Correction_bins[ipt+1])){
		track_weight = pPb_Smearing_Correction[ipt]*(1.0-pPb_FakeRate[ipt])/pPb_Efficiency[ipt];
		//fprintf(stderr,"\n %d: Low Edge=%f, High Edge=%f, Smear=%f, FakeRake=%f, Efficiency=%f\n",__LINE__,track_pT_Correction[ipt],track_pT_Correction[ipt+1],pPb_Smearing_Correction[ipt],pPb_FakeRate[ipt],pPb_Efficiency[ipt]);
	      }
	    }
	
	    if (TPC_Iso_Flag){
	      for (int ipt = 0; ipt < N_TPC_Track_pT_Bins; ipt++){
		if( (track_pt[itrack] >= TPC_track_pT_Correction_bins[ipt]) && (track_pt[itrack] < TPC_track_pT_Correction_bins[ipt+1])){
		  track_weight = TPC_pPb_Smearing_Correction[ipt]*(1.0-TPC_pPb_FakeRate[ipt])/TPC_pPb_Efficiency[ipt];
		  //fprintf(stderr,"\n %d: Low Edge=%f, High Edge=%f, Smear=%f, FakeRake=%f, Efficiency=%f\n",__LINE__,TPC_track_pT_Correction_bins[ipt],TPC_track_pT_Correction_bins[ipt+1],TPC_pPb_Smearing_Correction[ipt],TPC_pPb_FakeRate[ipt],TPC_pPb_Efficiency[ipt]);
		}
	      }
	    }
	    
	  if (Is_pp){
	    
	    track_weight = 1.0;

	    for (int ipt = 0; ipt < N_Track_pT_Bins; ipt++){
	      if( (track_pt[itrack] >= pp_track_pT_Correction_bins[ipt]) && (track_pt[itrack] < pp_track_pT_Correction_bins[ipt+1])){
		track_weight = pp_Smearing_Correction[ipt]*(1.0-pp_FakeRate[ipt])/pp_Efficiency[ipt];
	      }
	    }
	  }//pp
	
	  //fprintf(stderr,"\n Track pT = %f, Track weight = %f\n",track_pt[itrack],track_weight);
	  
		  
	  //Observables:
	  Double_t zt = track_pt[itrack]/cluster_pt[n];
	  Float_t DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[n] - track_phi[itrack]));
	  Float_t DeltaEta = cluster_eta[n] - track_eta[itrack];
	  if ((TMath::Abs(DeltaPhi) < 0.005) && (TMath::Abs(DeltaEta) < 0.005)) continue; //Match Mixing Cut

	  for (int ipt = 0; ipt < nptbins; ipt++){
	    if (cluster_pt[n] >= ptbins[ipt] && cluster_pt[n] < ptbins[ipt+1]){
	      for(int izt = 0; izt<nztbins ; izt++){
		if(zt>ztbins[izt] and  zt<ztbins[izt+1]){
		  if (first_cluster){
		    h_track_phi_eta[izt+ipt*nztbins]->Fill(track_phi[itrack],track_eta[itrack]);
		    if (track_pt[itrack] > 10)
		    tracks_tenGev += 1;
		  }
		  //2 DNN Regions

		  if (Signal and Isolated)
		    IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta,track_weight*purity_weight);

		  if (Background and Isolated){
		    BKGD_IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta,bkg_weight*track_weight*BR_purity_weight);
		    //not weighted with pT distro
		    BKGD_IsoCorr_UW[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta,track_weight);
		  }
		
	    	  
		  //no shower shape selection
		  if(Isolated)
		    Corr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);

		  Inclusive_Corr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);	       

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
  fprintf(stderr,"Tracks Above 10 GeV/c = %llu \n",tracks_tenGev);
			  
  h_purity.Write("purities");

  h_cluster_phi->Write();
  h_cluster_eta->Write();

  hweight.Write();

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

  Weights_Sum->Write();

  //Seperate zt loops for easier file reading
  fout->Close();     
  file->Close();  
  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}
