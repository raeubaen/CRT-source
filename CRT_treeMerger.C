#include <TTree.h>
#include <TFile.h>
#include <TList.h>


void CRT_treeMerger(){

  TString inPath =    "./step2/";
  TString outPath =   "./step2/";
  TString outName =   "run174to175_s2.root";
  TString treeName =  "CRT";

  TFile *files[2] = {  new TFile( inPath + "run174_s2.root"), 
                       new TFile( inPath + "run175_s2.root")};

  TFile *outf =        new TFile( outPath + outName, "CREATE");

  TList *lst = new TList;

  for(auto &f: files){
    TTree *tree = (TTree*)f->Get(treeName);
    lst->Add(tree);
  }

  TTree *newtree = TTree::MergeTrees(lst);
  newtree->SetName(treeName);
  outf->cd();
  newtree->Write();
}
