#include <TCanvas.h>
#include <vector>
#include <math.h>
#include "TVirtualFitter.h"

void Ratio_Fit_MinMax(){

  gStyle->SetOptFit(0100);

  const int NzT = 8;
  float zT_centers[NzT] = {0.07, 0.0935, 0.1245, 0.166, 0.2215, 0.295, 0.3935, 0.525};
  float zT_widths[NzT] = {0.01, 0.0135, 0.0175, 0.024, 0.0315, 0.042, 0.0565, 0.075};

  float Ratio[NzT] = {1.35399601, 0.98875884, 0.89971162, 1.89903617, 
		    0.81284172, 1.11789588, 0.50233703, 0.6904798};
  float Ratio_Stat[NzT] = {0.47663597, 0.24809698, 0.29891981, 0.82527748, 
			 0.25208921, 0.56261092, 0.22150608, 0.45576316};
  float Ratio_Sys[NzT] = {0.22787401, 0.22788458, 0.22796965, 0.22804687, 
			0.22871356, 0.23272182, 0.257535, 0.37089975};
      
  TFile* fout = new TFile("ToyMC_RatioFit.root","RECREATE");

  float max_array[NzT];
  float min_array[NzT];
  for (int izt = 0; izt < NzT; ++izt){
    float rando = Ratio_Sys[izt];
    max_array[izt] = Ratio[izt]+rando;
    min_array[izt] = Ratio[izt]-rando;   
  }
        
  TGraphErrors* Ratio_TGraph = new TGraphErrors(NzT,zT_centers,Ratio,zT_widths,Ratio_Stat);
  TGraphErrors* Max_Ratio_TGraph = new TGraphErrors(NzT,zT_centers,max_array,0,Ratio_Stat);
  TGraphErrors* Min_Ratio_TGraph = new TGraphErrors(NzT,zT_centers,min_array,0,Ratio_Stat);
  auto c1 = new TCanvas("c1","c1",200,10,600,400);
  Ratio_TGraph->Draw("AP");
  Ratio_TGraph->Fit("pol0","S0");
  TF1* fit = Ratio_TGraph->GetFunction("pol0");
  float C = fit->GetParameter(0);

  TGraphErrors *grint = new TGraphErrors(NzT);
  grint->SetTitle("Fitted line with .68 conf. band");
  for (int i=0; i<100; i++)
    grint->SetPoint(i, (0.66/100)*i, 0);
  /*Compute the confidence intervals at the x points of the created graph*/
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,0.68);
  //Now the "grint" graph contains function values as its y-coordinates
  //and confidence intervals as the errors on these coordinates
  //Draw the graph, the function and the confidence intervals

  grint->SetLineColorAlpha(3,1);
  grint->SetFillColorAlpha(3,0.2);
  grint->Draw("Same3");

  Ratio_TGraph->SetMarkerStyle(5);
  Ratio_TGraph->SetMarkerSize(0.7);
  

  Max_Ratio_TGraph->Fit("pol0","S");
  Min_Ratio_TGraph->Fit("pol0","S");

  TF1* max_fit = Max_Ratio_TGraph->GetFunction("pol0");
  TF1* min_fit = Min_Ratio_TGraph->GetFunction("pol0");

  float max_const = max_fit->GetParameter(0);
  float min_const = min_fit->GetParameter(0);
 
  std::cout<<max_const<<std::endl;
  std::cout<<min_const<<std::endl;

  TLine* Max_Line = new TLine(0.005,max_const,0.655,max_const);
  TLine* Min_Line = new TLine(0.005,min_const,0.655,min_const);
  Max_Line->SetLineColor(kBlue);
  Min_Line->SetLineColor(kRed);
  Max_Line->Draw("Same");
  Min_Line->Draw("Same");
  gPad->Modified(); gPad->Update();
}
