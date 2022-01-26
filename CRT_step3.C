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
double **timeOffset, **timeOffset_out, **zetaOffset, **zetaOffset_out, **chargeEqual, **chargeEqualErr, **chargeEqual_out, **chargeEqualErr_out;
double tDiff = 0, zeta=0;
double *intQ, *pkV, *teQ, *teT, *teA, *teB, *teX2, *ped;
list<double ** > arrayList = {&intQ, &pkV, &teQ, &teT, &teA, &teB, &teX2, &ped};
void InitVectors() { for(double** &arr: arrayList) { *arr = new double[2*scintNum](); } }
int iSc_out; double_t Z_out; double_t Q_out[2], X2_out[2], T_out[2];
TTree *CRTs3;



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
  
  chargeEqual_out[iSd][iSc] = l4.GetParameter(1);
  chargeEqualErr_out[iSd][iSc] = l4.GetParError(1);

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

  timeOffset_out[iSd][iSc] = timeFit.GetParameter(1);

}

void pedMip_proc(TH1* histObj, int histN, TString& histTag, int& histSkipFlag) {   

  int iSd = (int)((histN+1)>scintNum), iSc = histN - (iSd==1)*scintNum; 
  histTag = Form("[%d:%d] ",  iSd, iSc);

  gStyle->SetOptFit(1);

  double qpeakPed = histObj->GetBinCenter(histObj->GetMaximumBin());
  TF1 l = TF1("l", "gaus", qpeakPed-3, qpeakPed+3);
  histObj->Fit(&l, "R");

}

void zetaMip_proc(TH1* histObj, int histN, TString& histTag, int& histSkipFlag) {   

  int iSd = (int)((histN+1)>scintNum), iSc = histN - (iSd==1)*scintNum; 
  histTag = Form("[%d:%d] ",  iSd, iSc);

  gStyle->SetOptFit(1);

  TF1 zFit;
  if (centerMode) {
    zFit = TF1("l", "gaus", -20, 20);
    histObj->Fit(&zFit, "R");
    histObj->Fit(&zFit, "R");
  } else {
    zFit = TF1("l", "( TMath::TanH( (x-[0])/[1] ) - TMath::TanH( (x-[2]) /[1] ) ) * ([3]+[4]*x)", -1.2*scintL, 1.2*scintL); // non fitta mai
    zFit.SetParLimits(0, -100, -60);
    zFit.SetParLimits(1, 0, 20);
    zFit.SetParLimits(2, 60, 100);
    zFit.SetParLimits(3, 1, 1000);
    zFit.SetParLimits(4, 0.01, 2);
    histObj->Fit(&zFit, "R");
    histObj->Fit(&zFit, "R");
  }

  zetaOffset[0][histN] = centerMode ? zFit.GetParameter(1) : (zFit.GetParameter(2) + zFit.GetParameter(0))/2 ;

}

void qSharing_proc(TH1* histObj, int histN, TString& histTag, int& histSkipFlag) {   

  if (histN == 0 || histN == scintNum - 1) {histSkipFlag=1;}
  int iSd = (int)((histN+1)>scintNum), iSc = histN - (iSd==1)*scintNum; 
  histTag = Form("[%d:%d] ",  iSd, iSc);

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
    HM.AddHistBox("pedMip",      1, 2*scintNum,  "Pedestal",         "charge", "pC",    100, -10, 10,           1, 0, 0,   &pedMip_proc),
    HM.AddHistBox("voltPeak",    1, 2*scintNum,  "Wave peak",        "ampl", "V",       100, 0, 2000,           1, 0, 0,   &setHistTag2_proc),
    HM.AddHistBox("timeMip",     1, 2*scintNum,  "MIP times",        "time", "ns",      100, -30, 30,           1, 0, 0,   &timeMip_proc),
    HM.AddHistBox("zetaMip",     1,   scintNum,  "MIP zetas",        "zeta", "cm",      320, -scintL, scintL,   1, 0, 0,   &zetaMip_proc),

    HM.AddHistBox("q_chi2",      2,   scintNum,  "MIP q vs chi2",    "charge [pC]", "chi2",           0.5*qBins, qFrom, qTo, 100, 0, 40),
    HM.AddHistBox("zeta_q",      2,   scintNum,  "MIP q vs Z",       "zeta [cm]", "charge [pC]",      160, -scintL, scintL, 0.5*qBins, qFrom, qTo),
    HM.AddHistBox("qSharing",    2,   scintNum,  "Sharing",          "Q_i [pC]", "Q_neighbours [pC]", 50, qFrom, qTo, 50, qFrom, qTo, &qSharing_proc)

  };
}


void fill_raw(int hitN) {

  HM.Fill1d("chargeRaw", hitN, intQ[hitN]); 
  HM.Fill1d("chargeRaw", hitN, intQ[hitN]);
  HM.Fill1d("voltPeak", hitN, pkV[hitN]); 
  HM.Fill1d("voltPeak", hitN, pkV[hitN]);

}

void fill_mip(int iScHit) {

  HM.Fill1d("chargeMip", iScHit, intQ[iScHit]);
  HM.Fill1d("chargeMip", iScHit+scintNum, intQ[iScHit+scintNum]);
  HM.Fill1d("chargeTeMip", iScHit, teQ[iScHit]);
  HM.Fill1d("chargeTeMip", iScHit+scintNum, teQ[iScHit+scintNum]);
  HM.Fill1d("timeMip", iScHit, teT[iScHit]);
  HM.Fill1d("timeMip", iScHit+scintNum, teT[iScHit+scintNum]);
  HM.Fill1d("pedMip", iScHit, ped[iScHit]);
  HM.Fill1d("pedMip", iScHit+scintNum, ped[iScHit+scintNum]);
  HM.Fill1d("zetaMip", iScHit, zeta);

  HM.Fill2d("zeta_q", iScHit, zeta, teQ[iScHit]);
  HM.Fill2d("q_chi2", iScHit,  teX2[iScHit], teQ[iScHit]);

  jTrig_out = jentry; iSc_out = iScHit;
  Z_out = zeta;
  Q_out[0] = intQ[iScHit]; Q_out[1] = intQ[iScHit+scintNum];
  T_out[0] = teT[iScHit]; T_out[1] = teT[iScHit+scintNum];
  X2_out[0] = teX2[iScHit]; X2_out[1] = teX2[iScHit+scintNum];
  CRTs3->Fill();
}



void Analysis::ProcessPlots() {

  TDirectory* calib_dir = outFile->mkdir("calibration");
  calib_dir->cd();  

  //ChEq
    TGraphErrors *chEqGraph = new TGraphErrors(2*scintNum); chEqGraph->SetTitle("Charge equal");
    double qMean = 0;
    for (int k = 0; k < scintNum; k++) { 
      chEqGraph->SetPoint(k, (float)k+1-0.07, chargeEqual_out[0][k]);
      chEqGraph->SetPoint(k+scintNum,  (float)k+1+0.07, chargeEqual_out[1][k]);
      chEqGraph->SetPointError(k, 0, chargeEqualErr_out[0][k]);
      chEqGraph->SetPointError(k+scintNum, 0, chargeEqualErr_out[0][k]);
      qMean += chargeEqual_out[0][k] + chargeEqual_out[1][k];
    }
    qMean = qMean/(2*scintNum);
    TCanvas *eq_c = new TCanvas("chargeEqual", "chargeEqual"); eq_c->cd(); 
    chEqGraph->SetLineColor(kBlue); chEqGraph->SetMarkerColor(kBlue); chEqGraph->SetMarkerSize(1.4); chEqGraph->SetMarkerStyle(25); 
    chEqGraph->GetXaxis()->SetRangeUser(0, scintNum+1); chEqGraph->Draw("AP"); 
    TLine *line = new TLine(0.5, qMean, scintNum+0.5, qMean); line->SetLineColor(kRed); line->Draw("same");
    line = new TLine(0.5, 1.1*qMean, scintNum+0.5, 1.1*qMean); line->Draw("same");
    line = new TLine(0.5, 0.9*qMean, scintNum+0.5, 0.9*qMean); line->Draw("same");
    eq_c->Write("chargeEqual");
  //ChEq

  //ChEqOld
    chEqGraph = new TGraphErrors(2*scintNum); chEqGraph->SetTitle("chargeEqual_s3p1 (loaded from calib)");
    qMean = 0;
    for (int k = 0; k < scintNum; k++) { 
      chEqGraph->SetPoint(k, (float)k+1-0.07, chargeEqual[0][k]);
      chEqGraph->SetPoint(k+scintNum,  (float)k+1+0.07, chargeEqual[1][k]);
      qMean += chargeEqual[0][k] + chargeEqual[1][k];
    }
    qMean = qMean/(2*scintNum);
    eq_c = new TCanvas("chargeEqual_s3p1", "chargeEqual_s3p1"); eq_c->cd(); 
    chEqGraph->SetLineColor(kBlue); chEqGraph->SetMarkerColor(kBlue); chEqGraph->SetMarkerSize(1.4); chEqGraph->SetMarkerStyle(25); 
    chEqGraph->GetXaxis()->SetRangeUser(0, scintNum+1); chEqGraph->Draw("AP");  
    line = new TLine(0.5, qMean, scintNum+0.5, qMean); line->SetLineColor(kRed); line->Draw("same");
    line = new TLine(0.5, 1.1*qMean, scintNum+0.5, 1.1*qMean); line->Draw("same");
    line = new TLine(0.5, 0.9*qMean, scintNum+0.5, 0.9*qMean); line->Draw("same");
    eq_c->Write("chargeEqual_s3p1");
  //ChEqOld

  //TimeOff
    TGraphErrors *tim = new TGraphErrors(2*scintNum);  tim->SetTitle("timeOffset");
    for (int k = 0; k < scintNum; k++) { 
      tim->SetPoint(k, (float)k+1-0.07, timeOffset_out[0][k]);
      tim->SetPoint(k+scintNum,  (float)k+1+0.07, timeOffset_out[1][k]); 
    }
    TCanvas *tim_c = new TCanvas("timeOff", "timeOff"); tim_c->cd(); 
    tim->SetLineColor(kBlue); tim->SetMarkerColor(kBlue); tim->SetMarkerSize(1.4); tim->SetMarkerStyle(25); 
    tim->GetXaxis()->SetRangeUser(0, scintNum+1);
    tim->Draw("AP"); tim_c->Write("timeOff");
  //TimeOff

  //TimeOffOld
    tim = new TGraphErrors(2*scintNum);  tim->SetTitle("timeOffset_s3p1 (loaded from calib)");
    for (int k = 0; k < scintNum; k++) { 
      tim->SetPoint(k, (float)k+1-0.07, timeOffset[0][k]);
      tim->SetPoint(k+scintNum,  (float)k+1+0.07, timeOffset[1][k]); 
    }
    tim_c = new TCanvas("timeOff", "timeOff"); tim_c->cd(); 
    tim->SetLineColor(kBlue); tim->SetMarkerColor(kBlue); tim->SetMarkerSize(1.4); tim->SetMarkerStyle(25); 
    tim->GetXaxis()->SetRangeUser(0, scintNum+1);
    tim->Draw("AP"); tim_c->Write("timeOff_s3p1");
  //TimeOffOld

  //zetaOff
    TGraphErrors *zOffGraph = new TGraphErrors(scintNum); zOffGraph->SetTitle("Zeta offset");
    for (int k = 0; k < scintNum; k++) { 
      zOffGraph->SetPoint(k, k+1, zetaOffset_out[0][k]);
    }
    line = new TLine(0.5, 0, scintNum+0.5, 0); 
    TCanvas *zet_c = new TCanvas("zetaOffset", "zetaOffset"); zet_c->cd(); 
    zOffGraph->SetLineColor(kBlue); zOffGraph->SetMarkerColor(kBlue); zOffGraph->SetMarkerSize(1.4); zOffGraph->SetMarkerStyle(25); 
    zOffGraph->GetXaxis()->SetRangeUser(0, scintNum+1);
    line->SetLineColor(kRed);
    zOffGraph->Draw("AP"); line->Draw("same"); zet_c->Write("zetaOffset");
  //zetaOff

  //zetaOffOld
    zOffGraph = new TGraphErrors(scintNum); zOffGraph->SetTitle("zetaOffset_s3p2 (loaded from calib)");
    for (int k = 0; k < scintNum; k++) { 
      zOffGraph->SetPoint(k, k+1, zetaOffset[0][k]);
    }
    line = new TLine(0.5, 0, scintNum+0.5, 0); 
    zet_c = new TCanvas("zetaOffset", "zetaOffset"); zet_c->cd(); 
    zOffGraph->SetLineColor(kBlue); zOffGraph->SetMarkerColor(kBlue); zOffGraph->SetMarkerSize(1.4); zOffGraph->SetMarkerStyle(25); 
    zOffGraph->GetXaxis()->SetRangeUser(0, scintNum+1);
    line->SetLineColor(kRed);
    zOffGraph->Draw("AP"); line->Draw("same"); zet_c->Write("zetaOffset_s3p2");
  //zetaOffOld

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
      double chCal = chEqReference/chargeEqual[hitSide][hitScint];
      //chCal = (chCal>0.8 && chCal<1.2)?chCal:1;

      intQ[hitN] = Qval[hit]*chCal;
      pkV[hitN] = Vmax[hit];
      ped[hitN] = pedL[hit];
      teA[hitN] = templFit[hit][0];
      teQ[hitN] = templ2charge(teA[hitN])*chCal;
      teB[hitN] = templFit[hit][2];
      teT[hitN] = templTime[hit] - timeOffset[hitSide][hitScint];
      teX2[hitN] = templChi2[hit]; 

      if (Selection.isSaturated(pkV[hitN])) {skipFlag = 1; continue;}

      fill_raw(hitN);
    } 
    
    if (skipFlag) {continue;}

    for(int isc = 0; isc < scintNum; isc++) {

      if ( Selection.hitPrecheck(isc, iScint, nCry) && Selection.isChargeGood(intQ, isc) ) {

        //Fill
          HM.Fill2d("q_chi2", isc, teQ[isc], teX2[isc]);
          HM.Fill2d("q_chi2", isc+scintNum, teQ[isc+scintNum],  teX2[isc+scintNum]); 
          if( (isc+1)%scintNum > 1 ) { HM.Fill2d("qSharing", isc, intQ[isc] + intQ[isc+scintNum],  intQ[isc-1] + intQ[isc+scintNum-1] + intQ[isc+1] + intQ[isc+scintNum+1]); }
        //Fill

        if ( Selection.isX2Good(teX2, isc) && !Selection.isShared(intQ, isc) ) { iScHit = isc; m++; }
      } 
    } 
    
    if ( m != 1 ) {continue;}

    if ( !Selection.isTimeGood(teT[iScHit])) {continue;}

    tDiff = teT[iScHit] - teT[scintNum+iScHit]; 
    zeta = tDiff*scintVp/2-zetaOffset[0][iScHit];

    if ( !Selection.isZetaGood(zeta) ) {continue;}

    fill_mip(iScHit);    
  }

  cout<<endl;
}


void Analysis::Loop(){

  cout<<endl<<endl<<"::::::::::::::::::::: CRT analysis step 3 :::::::::::::::::::::"<<endl<<endl;

  if (fChain == 0) return;

  cout<<"Creating histograms:"<<endl;
  HM.SetOutFile(outFile);
  createHistBoxes();
  cout<<"...done"<<endl<<endl;

  cout<<"Retrieving calibration data from [" + lutPrefix3p + "] ..."<<endl;
  timeOffset  = CSV.InitMatrix(2, scintNum); timeOffset_out = CSV.InitMatrix(2, scintNum);
  chargeEqual = CSV.InitMatrix(2, scintNum); chargeEqual_out = CSV.InitMatrix(2, scintNum);
  chargeEqualErr = CSV.InitMatrix(2, scintNum); chargeEqualErr_out = CSV.InitMatrix(2, scintNum);
  zetaOffset  = CSV.InitMatrix(1, scintNum); zetaOffset_out  = CSV.InitMatrix(1, scintNum);
  CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutChEqName + "*"),      ',', chargeEqual, 2, scintNum);
  CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutTimeOffsName + "*"),  ',', timeOffset,  2, scintNum);
  CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutZetaOffName + "*"),   ',', zetaOffset,  1, scintNum);
  cout<<"...done"<<endl<<endl;

  outFile->cd();
  CRTs3 = new TTree("CRT","CRT");          
  CRTs3->SetAutoSave(1000);
  CRTs3->Branch("crt_iTrig",   &jTrig_out,  "crt_jTrig/I"); 
  CRTs3->Branch("crt_iSc",     &iSc_out,    "crt_iSc/I");
  CRTs3->Branch("crt_Z",       &Z_out,      "crt_Z/D");
  CRTs3->Branch("crt_Q",       &Q_out,      "crt_Q[2]/D");
  CRTs3->Branch("crt_T",       &T_out,      "crt_T[2]/D");
  CRTs3->Branch("crt_X2",      &X2_out,     "crt_X2[2]/D");

  Analysis::LoopOverEntries();
  HM.ProcessBoxes();
  Analysis::ProcessPlots();

  cout<<endl<<"Writing data to [" + lutPrefix3 + "] ..."<<endl;
  CSV.Write(lutPrefix3 + runName + lutChEqName + ".csv", ',', chargeEqual_out, 2, scintNum, 4);
  cout<<"...done"<<endl;

  outFile->cd();
  CRTs3->Write();
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
