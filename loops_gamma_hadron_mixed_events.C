#include <TCanvas.h>
#include <vector>
#include <math.h>

void loops_gamma_hadron_mixed_events() {

  TH2F* h = new TH2F("h", "delta eta, delta phi", 16, 0., TMath::Pi(), 50, -1., 1.);
  
  for (Int_t i_gamma = 0; i_gamma < 1000; ++i_gamma){
    // pick gamma eta and phi
    Float_t eta_gamma = gRandom->Uniform(-0.67, 0.67);
    Float_t phi_gamma = gRandom->Uniform(0. * TMath::DegToRad(), 107. * TMath::DegToRad()); 
    
    // pick charged hadron eta and phi
    for (Int_t i_hadron = 0; i_hadron < 10000; ++i_hadron){
      Float_t eta_hadron = gRandom->Uniform(-0.8, 0.8);
      Float_t phi_hadron = gRandom->Uniform(0., TMath::TwoPi());

      //****Roughly mimmic holes in ITS (only relative to EMCal)****
       if ((phi_hadron>25*TMath::DegToRad()) && (phi_hadron<45*TMath::DegToRad())) continue;
       if ((phi_hadron>245*TMath::DegToRad()) && (phi_hadron<270*TMath::DegToRad())) continue;
               
      // calculate delta phi and delta eta
      Float_t delta_eta = eta_gamma - eta_hadron;
      Float_t delta_phi = phi_gamma - phi_hadron;
      if (delta_phi < 0) delta_phi += TMath::TwoPi();
      if (delta_phi > TMath::Pi()) delta_phi = TMath::TwoPi() - delta_phi;
      h->Fill(delta_phi, delta_eta);
    }
  }

  h->SetXTitle("#Delta #phi");
  h->SetYTitle("#Delta #eta");
  h->Draw("lego2");
}
