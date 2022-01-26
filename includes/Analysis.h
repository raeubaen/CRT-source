#ifndef Analysis_h
#define Analysis_h

#define ROOT_CLASS ana
#include "ana.h"
#include "ana.C"

#define TREE_NAME "CRT"



class Analysis: public ROOT_CLASS
{
public:
  Analysis(TString infileName, TFile *f, TString, TString, int, TTree *tree = 0);
  virtual void Loop() override; //if not needed, fill with ROOT_CLASS::Loop();
  static void Run(TString infile, TString outfile, TString runName, TString calName, int window_close_handle);
  void LoopOverEntries();
  void ProcessHistos();
  void ProcessPlots();
  TFile *outFile;
  TString inFileName;
  TString runName;
  TString calName;
};

#endif



TTree * GetTree(TString infileName) {
  TTree *tree;
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(infileName.Data());
      if (!f || !f->IsOpen()) {
         f = new TFile(infileName.Data());
      }
      f->GetObject(TREE_NAME, tree);
  return tree;
}



Analysis::Analysis(TString infileName, TFile *f, TString runN, TString calN, int window_close_handle, TTree *tree) : ROOT_CLASS(GetTree(infileName)) { 
  outFile = f;
  inFileName = infileName;
  runName = runN;
  calName = calN;

};



void Analysis::Run(TString infile, TString outfile, TString runName, TString calName, int window_close_handle){

  TApplication *myapp;

  if (window_close_handle != -1) myapp = new TApplication("myapp", 0, 0);

  TFile *f = new TFile(outfile, "RECREATE");

  Analysis *a = new Analysis(infile, f, runName, calName, window_close_handle);

  a->Loop();

  if (window_close_handle != -1)  myapp->Run(true);
}
