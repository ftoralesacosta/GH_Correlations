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

void gamma_hadron_mixed_events() {

  TH2F* h = new TH2F("h", "delta eta, delta phi", 50, 0., TMath::Pi(), 50, -1., 1.);
  
  for (Int_t i_pair=0; i_pair< 100000000; ++i_pair) {

    // pick gamma eta and phi
    Float_t eta_gamma = gRandom->Uniform(-0.67, 0.67);
    Float_t phi_gamma = gRandom->Uniform(80. * TMath::DegToRad(), 187. * TMath::DegToRad()); 

    // pick charged hadron eta and phi
    Float_t eta_hadron = gRandom->Uniform(-0.8, 0.8);
    Float_t phi_hadron = gRandom->Uniform(0., TMath::TwoPi());

    // calculate delta phi and delta eta
    Float_t delta_eta = eta_gamma - eta_hadron;
    Float_t delta_phi = phi_gamma - phi_hadron;
    if (delta_phi < 0) delta_phi += TMath::TwoPi();
    if (delta_phi > TMath::Pi()) delta_phi = TMath::TwoPi() - delta_phi;
    Float_t delta_phia = TMath::Abs(TVector2::Phi_mpi_pi(phi_gamma - phi_hadron));
    h->Fill(delta_phi, delta_eta);

  }

  h->SetXTitle("#Delta #phi");
  h->SetYTitle("#Delta #eta");
  h->Draw("lego2");
}
