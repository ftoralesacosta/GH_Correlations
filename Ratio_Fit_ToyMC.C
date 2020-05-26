#include <TCanvas.h>
#include <vector>
#include <math.h>

void Ratio_Fit_ToyMC(){

  gStyle->SetOptFit(0100);

  const int NzT = 8;
  float zT_centers[NzT] = {0.07 , 0.0935, 0.1245, 0.166, 0.2215, 0.295, 0.3935, 0.525};
  float Ratio[NzT] = {1.35399601, 0.98875884, 0.89971162, 1.89903617, 
		    0.81284172, 1.11789588, 0.50233703, 0.6904798};
  float Ratio_Stat[NzT] = {0.47663597, 0.24809698, 0.29891981, 0.82527748, 
			 0.25208921, 0.56261092, 0.22150608, 0.45576316};
  float Ratio_Sys[NzT] = {0.22787401, 0.22788458, 0.22796965, 0.22804687, 
			0.22871356, 0.23272182, 0.257535, 0.37089975};

  TH1F* h = new TH1F("Constant_Fits","Distribution of Constant Fit",200,0.4,1.2);
  TH1F* chi = new TH1F("chi","Distribution of Constant Fit",200,0,10);
  TH1F* p = new TH1F("p","Distribution of Constant Fit",200,0,1);
      
  TFile* fout = new TFile("ToyMC_RatioFit.root","RECREATE");

  for (int i = 0; i < 1000000; ++i){

    float varry_array[NzT];
    for (int izt = 0; izt < NzT; ++izt){
      float rando = gRandom->Uniform(-Ratio_Sys[izt],Ratio_Sys[izt]);
      varry_array[izt] = Ratio[izt]+rando;   
    }
        
    TGraphErrors Ratio_TGraph = TGraphErrors(NzT,zT_centers,varry_array,0,Ratio_Stat);
    Ratio_TGraph.Fit("pol0","S");
    TF1* f = Ratio_TGraph.GetFunction("pol0");
    float c = f->GetParameter(0);
    float chi2_red = f->GetChisquare()/f->GetNDF();
    float pval = f->GetProb();
    h->Fill(c);
    chi->Fill(chi2_red);
    p->Fill(pval);
  }

  h->SetXTitle("Constant");
  h->Fit("gaus");
  h->Draw();

  h->Write();
  chi->Write();
  p->Write();
}
