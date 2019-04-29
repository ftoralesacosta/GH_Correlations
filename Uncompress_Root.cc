#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <iostream>

#define HI_TREE "_tree_event"
#define HI_TREE_2 "AliAnalysisTaskNTGJ/_tree_event"

void write_root(const char *filename_0)
{

  TFile *root_file = new TFile(filename_0,"update");
  TTree *hi_tree = dynamic_cast<TTree *>(root_file->Get(HI_TREE));
  if (hi_tree == NULL) {
    hi_tree = dynamic_cast<TTree *>(root_file->Get(HI_TREE_2));
    if(hi_tree == NULL){
      fprintf(stderr, "%s:%d: TREE FAIL\n",__FILE__, __LINE__);
      return;
    }
  }
  

  ULong64_t nentries = hi_tree->GetEntries();    
  size_t lastindex = std::string(filename_0).find_last_of("."); 
  std::string rawname = std::string(filename_0).substr(0, lastindex);
    
  TFile *newfile = new TFile(Form("%s_%s.root",rawname.data(),"Uncompressed"),"recreate","uncomp",0);// 0 means no compression
  TTree *newtree = hi_tree->CloneTree(0);
    
  for (ULong64_t t = 0; t<nentries;t++){ //Event # is key used in map <Matches>
    fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, t, nentries);
    hi_tree->GetEntry(t);
                
    newtree->Fill();  
        
  }//End loop over entries
  newtree->Write();

  delete root_file;
  delete newfile;
  
  gSystem->Exit(0);
}

int main(int argc, char *argv[])
{
  if (argc < 1) {
    fprintf(stderr,"%s\n","Argument Syntax is [Command] [File]");
    return EXIT_FAILURE;
  }

  write_root(argv[1]);
}
