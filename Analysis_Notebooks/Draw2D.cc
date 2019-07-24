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

  if (argc < 2){
    fprintf(stderr,"[File] [Inclusive/Isolated/SR/BR] \n");
    return 0;
  }

  TString Corr_Type = (TString)argv[2];
  TFile *file = TFile::Open((TString)argv[1]);
  
  TH2D *histo2D;
  TH1D *ntrig_histo;

  if (strcmp(Corr_Type.Data(),"Inclusive") == 0){
    histo2D = (TH2D*)file->Get(Form("Inclusive_Correlation__pT%1.0d_%1.0d__zT%1.0d_zT%1.0d",12,40,6,8));
			      //,pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]

    ntrig_histo = (TH1D*)file->Get(Form("N_Inclusive_Triggers_pT%1.0d_%1.0d",12,40));

    histo2D->SetTitle("#gamma-H Inclusive Correlation");
  }

  if (strcmp(Corr_Type.Data(),"Isolated") == 0){
    histo2D = (TH2D*)file->Get(Form("Correlation__pT%1.0d_%1.0d__zT%1.0d_zT%1.0d",12,40,6,8));
    //,pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]                                      
 
    ntrig_histo = (TH1D*)file->Get(Form("N_Triggers_pT%1.0d_%1.0d",12,40));
 
   histo2D->SetTitle("#gamma-H Isolated Correlation");
  }

  if (strcmp(Corr_Type.Data(),"SR") == 0){
    histo2D = (TH2D*)file->Get(Form("DNN1_Correlation__pT%1.0d_%1.0d__zT%1.0d_zT%1.0d",12,40,6,8));
    //,pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]                                      

    ntrig_histo = (TH1D*)file->Get(Form("N_DNN1_Triggers_pT%1.0d_%1.0d",12,40));

    histo2D->SetTitle("#gamma-H Signal Region Correlation");
  }

  if (strcmp(Corr_Type.Data(),"BR") == 0){
    histo2D = (TH2D*)file->Get(Form("DNN2_Correlation__pT%1.0d_%1.0d__zT%1.0d_zT%1.0d",12,40,6,8));
    //,pTbins[ipt],pTbins[ipt+1],100*zTbins[izt],100*zTbins[izt+1]                                                                
    ntrig_histo = (TH1D*)file->Get(Form("N_DNN2_Triggers_pT%1.0d_%1.0d",12,40));

    histo2D->SetTitle("#gamma-H Background Region Correlation");
  }

  if (ntrig_histo != NULL){
      histo2D->Scale(1.0/ntrig_histo->GetEntries());}
  //histo2D->SetZTitle("#frac{1}{N}");

  TCanvas *c = new TCanvas("c","c",200,200);
  histo2D->Draw("SURF2");

  TAxis* x = histo2D->GetXaxis();
  x->SetTitle("#Delta#varphi");
  x->SetTitleOffset(1.8);
  x->SetTitleSize(0.04);
  x->SetLabelSize(0.03);
  x->CenterTitle();

  TAxis* y = histo2D->GetYaxis();
  y->SetTitle("#Delta#eta");
  y->SetTitleOffset(1.8);
  y->SetTitleSize(0.04);
  y->SetLabelSize(0.03);
  y->CenterTitle();

  TAxis* z = histo2D->GetZaxis();
  //z->SetTitle("#frac{1}{#it{N}_{#gamma}} #frac{d^{2}}{d#Delta#etad#Delta#varphi}");
  z->SetTitle("#frac{1}{#it{N}_{#gamma}}");
  z->SetTitleOffset(1.8);
  z->SetTitleSize(0.04);
  z->SetLabelSize(0.03);
  z->CenterTitle();

  gStyle->SetOptStat(0);
  c->SetRightMargin(0.09);
  c->SetLeftMargin(0.17);
  c->SetBottomMargin(0.1);
  c->Update();
  c->SaveAs(Form("2D_pics/2D_%s.pdf",Corr_Type.Data()));
  //fprintf(stderr,"%s",Corr_Type.Data());
}
