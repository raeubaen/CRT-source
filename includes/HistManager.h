#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <vector>
#include <iostream>
#include <TApplication.h>
#include <TH1.h>
#include <stdlib.h>
#include <TFrame.h>
#include <unordered_map>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TObjString.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TF1.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <vector>
#include <TString.h>
#include <TTimer.h>
#include <TLatex.h>
#include <TGClient.h>
#include <TRootCanvas.h>
#include <TDatime.h>
#include <time.h>

using namespace std;


//void processHistBox(TH1* histObj, int histN, TString& histTag, int& histSkipFlag) {}


void _processHistObj_doNothing(TH1* obj, int n, TString& tag, int& skip) {};


class HistBox{

  private:

    TString _label1, _label2, _histTitle, _histosName;
    TDirectory  *_histDir, *_histosDir;
    TFile *_outfile;
    TH1 *_histosObj;
    TH1 *_histObj;
    int _histType, _histN;

    TString histTag, histTagFormat = "[%d] ";
    int skipFlag = 0;

    void (* _processHistObj)(TH1*, int, TString&, int&) = &_processHistObj_doNothing;

  public:
   
    HistBox(
      TFile *f,
      TDirectory *histosDir,
      TString histosName,
      int histType,
      int histN, 
      TString histTitle, TString label1, TString label2,
      int xBins, double xMin, double xMax, 
      int yBins = 1, double yMin = 0., double yMax = 0.,
      void (*processFunction)(TH1*, int, TString&, int&) = &_processHistObj_doNothing
    ):  
      _histosDir(histosDir),
      _processHistObj(processFunction),
      _histTitle(histTitle),
      _histType(histType),
      _histN(histN),
      _histosName(histosName),
      _outfile(f)
    { 
      cout<<"   HistBox constructor: ["<<histType<<"] ["<<histN<<"] ["<<histosName<<"] ["<<histTitle<<"]"<<endl;

      if (histType == 1) {
        if (label2.CompareTo("")) {
          _label2 = Form("Entries / %.1f %s", ((xMax-xMin)/(double)xBins), label2.Data());
          _label1 = label1 + " [" + label2 + "]";
        } else {
         _label2 = Form("Entries");
         _label1 = label1;
        }
        _histosObj = new TH2F(histosName + "_", histTitle, histN, 0, histN, xBins, xMin, xMax);
      } else if (histType == 2) { 
        _label1 = label1;
        _label2 = label2;
        _histosObj = new TH3F(histosName + "_", histTitle, histN, 0, histN, xBins, xMin, xMax, yBins, yMin, yMax);       
      }
    }

    TH1* GetHistosObj(){ return _histosObj; }

    void ProcessBox(){
      
      _outfile->cd();
      _histDir = _outfile->mkdir(_histosName, "recreate");
      _histDir->cd();

      for (int k = 0; k < _histN; k++) { 
        
        if (_histType == 1) {
          TH2F *histosObj_tmp = (TH2F*)_histosObj; 
          TH1F *_histObj_tmp = (TH1F*)(histosObj_tmp->ProjectionY(_histosName, k+1, k+1));
          _histObj = new TH1F();  
          _histObj = (TH1F*)_histObj_tmp;
        } else if (_histType == 2) {
          TH3F *histosObj_tmp = (TH3F*)_histosObj;
          histosObj_tmp->GetXaxis()->SetRange(k+1,k+1);
          TH2D* _histObj_tmp = static_cast<TH2D*>(histosObj_tmp->Project3D("zy"));
          _histObj = new TH2F();  
          _histObj = (TH2F*)_histObj_tmp;
        } else {
          continue;
        }

        skipFlag = 0;
        histTag = Form(histTagFormat, k);

        _histObj->SetName(_histosName);
        _histObj->SetTitle(histTag + _histTitle);
        _histObj->GetXaxis()->SetTitle(_label1.Data());
        _histObj->GetYaxis()->SetTitle(_label2.Data());

        _outfile->cd();
        _histDir->cd();

        _processHistObj(_histObj, k, histTag, skipFlag);

        if (skipFlag) {continue;}

        _outfile->cd();
        _histDir->cd();
        _histObj->Write(histTag + _histosName);

      }

      _outfile->cd();
      _histosDir->cd();

      if(_histType == 1) {
      TCanvas *temp_c = new TCanvas(_histosName, _histosName); 
      temp_c->cd(); 
      _histosObj->Draw("zcol");
      temp_c->Write(_histosName);
      } else {
      _histosObj->Write(_histosName);
      }
      
    }
};


class HistManager{

  private:

    const TString _histosDirName = "HistBoxes";
    TFile *_outfile;

  public:

    HistManager(){};

    TDirectory *_histosDir;

    void SetOutFile(TFile *outfile){
      _outfile = outfile; 
      _outfile->cd();
      _histosDir = _outfile->mkdir(_histosDirName);
    };

    pair<string, HistBox*> AddHistBox(
      TString, int, int, TString, TString, TString,
      int, double, double, 
      int, double, double,
      void (*processFunction)(TH1*, int, TString&, int&)
    );

    unordered_map<string, HistBox *> HistBoxes;
    
    TH1  *GetHist(string name) {return (TH1*)HistBoxes[name]->GetHistosObj();}
    TH3F *GetHist2d(string name) {return (TH3F*)HistBoxes[name]->GetHistosObj();}
    TH2F *GetHist1d(string name) {return (TH2F*)HistBoxes[name]->GetHistosObj();}

    void Fill1d(string name, int n, double val1) {GetHist1d(name)->Fill(n, val1);}
    void Fill2d(string name, int n, double val1, double val2) {GetHist2d(name)->Fill(n, val1, val2);}

    void ProcessBoxes() {
      for (auto& HistBox: HistBoxes) {
        HistBox.second->ProcessBox();
      }
    };

};


pair<string, HistBox *> HistManager::AddHistBox(
  TString histosName,
  int histType,
  int histN, 
  TString histTitle, TString label1, TString label2,
  int xBins, double xMin, double xMax, 
  int yBins = 1, double yMin = 0., double yMax = 0.,
  void (*processFunction)(TH1*, int, TString&, int&) = &_processHistObj_doNothing
) {
  return {
    histosName.Data(), new HistBox(_outfile, _histosDir, histosName, histType, histN, 
    histTitle, label1, label2, xBins, xMin, xMax, yBins, yMin, yMax,
    processFunction)
  };
};









