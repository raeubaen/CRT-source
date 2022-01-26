#include <fstream>
#include <chrono>
#include <TLine.h>
#include <iostream>
#include <list>

#include "TApplication.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1F.h"
#include "TF1.h"

#include "includes/Analysis.h"
#include "includes/langaus.h"
#include "includes/AnaPars.h"
#include "includes/CsvHandler.h"
#include "includes/MipSelection.h"
#include "includes/templ2charge.h"
#include "includes/HistManager.h"


using namespace std;

CsvHandler CSV;
MipSelection Selection; 
HistManager HM;

Long64_t nentries, nbytes, nb, ientry, jentry, jTrig_out;
double **timeOffset, **zetaOffset, **chargeEqual, **chargeEqualErr;
double tDiff = 0, zeta=0;
double *Q, *V, *teQ, *teT, *teA, *teB, *teX2, *ped;
list<double ** > arrayList = {&Q, &V, &teQ, &teT, &teA, &teB, &teX2, &ped};
void InitVectors() { for(double** &arr: arrayList) { *arr = new double[2*scintNum](); } }
int iSc_out; double_t Z_out; double_t Q_out[2], X2_out[2], T_out[2];



void chargeMip_proc(TH1* histObj, int histN, TString& histTag, int& histSkipFlag) {   

  int iSd = (int)((histN+1)>scintNum), iSc = histN - (iSd==1)*scintNum; 
  histTag = Form("[%d:%d] ",  iSd, iSc);

  gStyle->SetOptFit(1); 

  //TSpectrum s(2);
  //s.Search(histObj, 2, "nodraw");
  //double *pks = s.GetPositionX();
  //double qpeak = *std::max_element(pks, pks + 2);

  double min_tmp = histObj->GetXaxis()->GetXmin();
  double max_tmp = histObj->GetXaxis()->GetXmax();
  histObj->GetXaxis()->SetRangeUser(390,700);
  double qpeak = histObj->GetBinCenter(histObj->GetMaximumBin());
  histObj->GetXaxis()->SetRangeUser(qFrom, qTo);
  double qmax = qpeak + 500, qmin = qpeak-20; float pk,sigma;

  TF1 l1 = TF1("l", "landau", qmin, qmax);                 l1.SetParameters(histObj->Integral()/2, 400, 100);                histObj->Fit(&l1, "R"); 
  pk = l1.GetMaximumX(); sigma = l1.GetParameter(2);
  TF1 l2 = TF1("l", "landau", pk-40, pk+3*sigma);          l2.SetParameters(l1.GetParameter(0), l1.GetParameter(1), sigma);  histObj->Fit(&l2, "R");
  pk = l2.GetParameter(1); sigma = l2.GetParameter(2); 
  TF1 l3 = TF1("l", "landau", pk-1*sigma, pk+4*sigma);     l3.SetParameters(l2.GetParameter(0), l2.GetParameter(1), sigma);  histObj->Fit(&l3, "R");
  pk = l3.GetParameter(1); sigma = l3.GetParameter(2);
  TF1 l4 = TF1("l", "landau", pk-0.8*sigma, pk+4*sigma);   l4.SetParameters(l3.GetParameter(0), l3.GetParameter(1), sigma);  histObj->Fit(&l4, "R");
  
  chargeEqual[iSd][iSc] = l4.GetParameter(1);
  chargeEqualErr[iSd][iSc] = l4.GetParError(1);

}

void timeMip_proc(TH1* histObj, int histN, TString& histTag, int& histSkipFlag) {   

  int iSd = (int)((histN+1)>scintNum), iSc = histN - (iSd==1)*scintNum; 
  histTag = Form("[%d:%d] ",  iSd, iSc);

  gStyle->SetOptFit(1);

  double tpeak = histObj->GetBinCenter(histObj->GetMaximumBin());
  double tmax = tpeak + 40, tmin = tpeak - 40;
  //TF1 timeFit = TF1("logn", "[0]*ROOT::Math::lognormal_pdf(x,log([1]),log([2]))", tmin, tmax);  
  TF1 timeFit = TF1("g", "gaus", tmin, tmax); timeFit.SetParameter(1, tpeak); timeFit.SetParameter(2, 2);

  histObj->Fit(&timeFit, "R");

  timeOffset[iSd][iSc] = timeFit.GetParameter(1);

}

void zetaMip_proc(TH1* histObj, int histN, TString& histTag, int& histSkipFlag) {   

  histSkipFlag = 1;

  int iSd = (int)((histN+1)>scintNum), iSc = histN - (iSd==1)*scintNum; 
  histTag = Form("[%d:%d] ",  iSd, iSc);

  gStyle->SetOptFit(1);

  TF1 zFit = TF1("l", "( TMath::TanH( (x-[0])/[1] ) - TMath::TanH( (x-[2]) /[1] ) ) * ([3]+[4]*x)", -1.2*scintL, 1.2*scintL);
  zFit.SetParLimits(0, -100, -60);
  zFit.SetParLimits(1, 0, 20);
  zFit.SetParLimits(2, 60, 100);
  zFit.SetParLimits(3, 1, 1000);
  zFit.SetParLimits(4, 0.01, 2);
  histObj->Fit(&zFit, "R");
  histObj->Fit(&zFit, "R");

  zetaOffset[0][histN] = (zFit.GetParameter(2) + zFit.GetParameter(0))/2;;

}

void setHistTag2_proc(TH1* histObj, int histN, TString& histTag, int& histSkipFlag) {   

  int iSd = (int)((histN+1)>scintNum), iSc = histN - (iSd==1)*scintNum; 
  histTag = Form("[%d:%d] ",  iSd, iSc);

}


void createHistBoxes() {
  HM.HistBoxes = {

    HM.AddHistBox("chargeRaw",   1, 2*scintNum,  "Raw charges",      "charge", "pC",    qBins, 20, qTo,         1, 0, 0,   &setHistTag2_proc),
    HM.AddHistBox("chargeMip",   1, 2*scintNum,  "MIP charges",      "charge", "pC",    qBins, qFrom, qTo,      1, 0, 0,   &chargeMip_proc),
    HM.AddHistBox("chargeTeMip", 1, 2*scintNum,  "MIP template q",   "charge", "pC",    qBins, qFrom, qTo,      1, 0, 0,   &setHistTag2_proc),
    HM.AddHistBox("voltPeak",    1, 2*scintNum,  "Wave peak",        "ampl", "V",       100, 0, 2000,           1, 0, 0,   &setHistTag2_proc),
    HM.AddHistBox("timeMip",     1, 2*scintNum,  "Mip times",        "time", "ns",      100, 100, 400,          1, 0, 0,   &timeMip_proc),
    HM.AddHistBox("zetaMip",     1,   scintNum,  "Mip zetas",        "zeta", "cm",      320, -scintL, scintL,   1, 0, 0,   &zetaMip_proc),

  };
}


void fill_raw(int hitN) {

  HM.Fill1d("chargeRaw", hitN, Q[hitN]); 
  HM.Fill1d("chargeRaw", hitN, Q[hitN]);
  HM.Fill1d("voltPeak", hitN, V[hitN]); 
  HM.Fill1d("voltPeak", hitN, V[hitN]);

}

void fill_mip(int iScHit) {

  HM.Fill1d("chargeMip", iScHit, Q[iScHit]);
  HM.Fill1d("chargeMip", iScHit+scintNum, Q[iScHit+scintNum]);
  HM.Fill1d("chargeTeMip", iScHit, teQ[iScHit]);
  HM.Fill1d("chargeTeMip", iScHit+scintNum, teQ[iScHit+scintNum]);
  HM.Fill1d("timeMip", iScHit, teT[iScHit]);
  HM.Fill1d("timeMip", iScHit+scintNum, teT[iScHit+scintNum]);
  HM.Fill1d("zetaMip", iScHit, zeta);

}



void Analysis::ProcessPlots() {

  TDirectory* calib_dir = outFile->mkdir("calibration");
  calib_dir->cd();  

  //ChEq
    TGraphErrors *chEqGraph = new TGraphErrors(2*scintNum); chEqGraph->SetTitle("Charge equal");
    double qMean = 0;
    for (int k = 0; k < scintNum; k++) { 
      chEqGraph->SetPoint(k, (float)k+1-0.07, chargeEqual[0][k]);
      chEqGraph->SetPoint(k+scintNum,  (float)k+1+0.07, chargeEqual[1][k]);
      chEqGraph->SetPointError(k, 0, chargeEqualErr[0][k]);
      chEqGraph->SetPointError(k+scintNum, 0, chargeEqualErr[0][k]);
      qMean += chargeEqual[0][k] + chargeEqual[1][k];
    }
    qMean = qMean/(2*scintNum);
    TLine *line = new TLine(0.5, qMean, scintNum+0.5, qMean); line->SetLineColor(kRed);
    TCanvas *eq_c = new TCanvas("chargeEqual", "chargeEqual"); eq_c->cd(); 
    chEqGraph->SetLineColor(kBlue); chEqGraph->SetMarkerColor(kBlue); chEqGraph->SetMarkerSize(1.4); chEqGraph->SetMarkerStyle(25); 
    chEqGraph->GetXaxis()->SetRangeUser(0, scintNum+1); chEqGraph->Draw("AP"); 
    line = new TLine(0.5, qMean, scintNum+0.5, qMean); line->SetLineColor(kRed); line->Draw("same");
    line = new TLine(0.5, 1.1*qMean, scintNum+0.5, 1.1*qMean); line->Draw("same");
    line = new TLine(0.5, 0.9*qMean, scintNum+0.5, 0.9*qMean); line->Draw("same");
    eq_c->Write("chargeEqual");
  //ChEq

  //TimeOff
    TGraphErrors *tim = new TGraphErrors(2*scintNum); tim->SetTitle("Time offset");
    for (int k = 0; k < scintNum; k++) { 
      tim->SetPoint(k, (float)k+1-0.07, timeOffset[0][k]);
      tim->SetPoint(k+scintNum,  (float)k+1+0.07, timeOffset[1][k]); 
    }
    TCanvas *tim_c = new TCanvas("timeOff", "timeOff"); tim_c->cd(); 
    tim->SetLineColor(kBlue); tim->SetMarkerColor(kBlue); tim->SetMarkerSize(1.4); tim->SetMarkerStyle(25); 
    tim->GetXaxis()->SetRangeUser(0, scintNum+1);
    tim->Draw("AP"); tim_c->Write("timeOff");
  //TimeOff

  //zetaOff
    // TGraphErrors *zOffGraph = new TGraphErrors(scintNum); zOffGraph->SetTitle("Zeta offset");
    // for (int k = 0; k < scintNum; k++) { 
    //   zOffGraph->SetPoint(k, k+1, zetaOffset[0][k]);
    // }
    // line = new TLine(0.5, 0, scintNum+0.5, 0); 
    // TCanvas *zet_c = new TCanvas("zetaOffset", "zetaOffset"); zet_c->cd(); 
    // zOffGraph->SetLineColor(kBlue); zOffGraph->SetMarkerColor(kBlue); zOffGraph->SetMarkerSize(1.4); zOffGraph->SetMarkerStyle(25); 
    // zOffGraph->GetXaxis()->SetRangeUser(0, scintNum+1);
    // line->SetLineColor(kRed);
    // zOffGraph->Draw("AP"); line->Draw("same"); zet_c->Write("zetaOffset");
  //zetaOff
 
}



void Analysis::LoopOverEntries() {

  nentries = fChain->GetEntriesFast(); 
  Long64_t etp = min(maxEvToProcess3p1, nentries);
  cout << "Number of events to process: " << etp << endl << endl;
  nbytes = 0, nb = 0;
  
  for (jentry = 0; jentry < etp; jentry++) {

    ientry = LoadTree(jentry);

    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;

    if (!(jentry%10000)) {cout << Form( "     processing evt %lld / %lld  ( %.0f%% )", jentry, etp, (float)(100*jentry/etp) ) << endl;}

    InitVectors();
    int skipFlag = 0, m = 0, iScHit = -1;
    
    for(int hit = 0; hit < nCry; hit++){

      int hitSide=iSide[hit], hitScint = iScint[hit], hitN = hitSide*scintNum + hitScint;
      Q[hitN] = Qval[hit];
      V[hitN] = Vmax[hit];
      ped[hitN] = pedL[hit];
      teA[hitN] = templFit[hit][0];
      teQ[hitN] = templ2charge(teA[hitN]);
      teB[hitN] = templFit[hit][2];
      teT[hitN] = templTime[hit];
      teX2[hitN] = templChi2[hit];                            

      if (Selection.isSaturated(V[hitN])) {skipFlag = 1; continue;}

      fill_raw(hitN);
    } 
    
    if (skipFlag) {continue;}  

    for(int isc = 0; isc < scintNum; isc++) {

      if ( Selection.hitPrecheck(isc, iScint, nCry) && Selection.isChargeGood(Q, isc) ) {

        if ( Selection.isX2Good(teX2, isc) && !Selection.isShared(Q, isc) ) { iScHit = isc; m++; }
      } 
    } 
    
    if ( m != 1 ) {continue;}

    fill_mip(iScHit);   
  }

  cout<<endl;
}



void Analysis::Loop(){

  cout<<endl<<endl<<"::::::::::::::::::::: CRT analysis step 3p1 :::::::::::::::::::::"<<endl<<endl;

  if (fChain == 0) return;

  cout<<"Creating histograms:"<<endl;
  HM.SetOutFile(outFile);
  createHistBoxes();
  cout<<"...done"<<endl<<endl;

  chargeEqual = CSV.InitMatrix(2, scintNum);
  zetaOffset  = CSV.InitMatrix(1, scintNum);
  timeOffset  = CSV.InitMatrix(2, scintNum);
  chargeEqual = CSV.InitMatrix(2, scintNum);
  chargeEqualErr = CSV.InitMatrix(2, scintNum);
  // cout<<"Retrieving calibration data from [" + lutPrefix3p + "] ..."<<endl;
  // CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutChEqName + "*"),      ',', chargeEqual, 2, scintNum);
  // CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutTimeOffsName + "*"),  ',', timeOffset,  2, scintNum);
  // CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutZetaOffName + "*"),   ',', zetaOffset,  1, scintNum);
  // cout<<"...done"<<endl<<endl;

  Analysis::LoopOverEntries();
  HM.ProcessBoxes();
  Analysis::ProcessPlots();

  cout<<"Writing calibration data to [" + lutPrefix3p + "] ..."<<endl;
  CSV.Write(lutPrefix3p + runName + lutChEqName         + ".csv", ',', chargeEqual, 2, scintNum, 4);
  CSV.Write(lutPrefix3p + runName + lutTimeOffsName     + ".csv", ',', timeOffset,  2, scintNum, 4);
  cout<<"...done"<<endl;

  outFile->Close();

  cout<<endl<<endl<<"::::::::::::::::::::: analysis done :::::::::::::::::::::"<<endl<<endl;

}



int main(int argc, char*argv[]) { 

  if (argc != 5) {
    printf("Usage: %s [infile_name] [outfile_name] [run_name] [calib_name]\n", argv[0]);
    exit(-1);
  }

  Analysis::Run(argv[1], argv[2], argv[3], argv[4], -1);
}
