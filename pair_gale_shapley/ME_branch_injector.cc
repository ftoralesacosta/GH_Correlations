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
#include <sstream>

#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>

int main(int argc, char *argv[])
{
  if (argc < 2) {
    fprintf(stderr,"\nSyntax is [Command] [root file with ME] [Unmixed root file] ");
    exit(EXIT_FAILURE);
  }

  //open rootfile with mixed event branch
  std::cout << "Opening: " << (TString)argv[1] << std::endl;
  TFile *me_file = TFile::Open((TString)argv[1]);

  if (me_file == NULL) {
    std::cout << " fail" << std::endl;
    exit(EXIT_FAILURE);
  }
  me_file->Print();

  TTree *me_tree_event = dynamic_cast<TTree *> (me_file->Get("_tree_event"));
  if (me_tree_event == NULL) {
    std::cout << "Failed to grab tree, perhaps AliAnalysisTaskNTGJ does not exist, trying again" << std::endl;
    me_tree_event = dynamic_cast<TTree *> (dynamic_cast<TDirectoryFile *>   (me_file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
    if (me_tree_event == NULL) {
      std::cout << " fail " << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  //_tree_event->Print();                                                                                                                                                                                                                                               
  std::cout<<"TTree successfully acquired" << std::endl;
  std::cout << " Total Number of entries in TTree: " << me_tree_event->GetEntries() << std::endl;

  //Get Branch with mixed events in it

  Long64_t mix_events[300];
  me_tree_event->SetBranchAddress("mixed_events", mix_events);

  //Fill  2D array with mixed events
  const Long64_t n_paired_events = me_tree_event->GetEntries();
  std::cout<<"got entries ok"<<std::endl;
  
  std::vector<vector<Long64_t> > mixed_events_all;
  std::cout<<"made 2D array ok"<<std::endl;

  for (Long64_t ievent = 0; ievent < n_paired_events; ievent++){
    //fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, n_paired_events);
    vector<Long64_t> temp_vec;
    me_tree_event->GetEntry(ievent);
    for (int imix = 0; imix < 300; imix ++)
      temp_vec.push_back(mix_events[imix]);
    mixed_events_all.push_back(temp_vec);
    std::cout<<mixed_events_all[ievent][149]<<std::endl;
  }
}
