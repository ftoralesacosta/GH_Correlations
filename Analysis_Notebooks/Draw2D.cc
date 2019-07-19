#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TH2D.h>

int main(int argc, char *argv[]){

  TString Corr_Type = (TString)argv[2];
  TFile *file = TFile::Open((TString)argv[1]);
  
  TH2D *histo2D;
  TH1D *ntrig_histo;

  if (strcmp(Corr_Type.Data(),"Inclusive") == 0){
    histo2D = (TH2D*)file->Get(Form("Inclusive_Correlation__pT%1.0d_%1.0d__zT%1.0d_zT%1.0d",12,40,6,8));
			      //,pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]
    ntrig_histo = (TH1D*)file->Get(Form("N_Inclusive_Triggers_pT%1.0d_%1.0d",12,40));

    histo2D->Scale(1.0/ntrig_histo->GetEntries());
  }

  TCanvas *c = new TCanvas("c","c",200,200);
  histo2D->Draw("SURF2");
  gStyle->SetOptStat(0);
  c->Update();
  c->SaveAs("./2D_Test.pdf");
  //fprintf(stderr,"%s",Corr_Type.Data());


}
