#ifndef AnaPars_H
#define AnaPars_H


    const int centerMode = 1;

    const Long64_t maxEvToProcess3p1 = 1e3;
    const Long64_t maxEvToProcess3p2 = 1e3;
    const Long64_t maxEvToProcess3 = 1e3;

    const TString lutPrefix3p = "./data/calibration/luts_s3p/";
    const TString lutPrefix3 = "./data/calibration/luts_s3/";
    const TString lutChEqName = "_chargEq";
    const TString lutTimeOffsName = "_timeOff";
    const TString lutZetaOffName = "_zetaOff";








    const int sideNum = 2;
    const int scintNum = 8;

    const float scintVp = 12.5;
    const float scintL = 160.0;

    const double minQCut = 50;
    const double maxQCut = 100000;
    const double maxVpeak = 1800;
    const double maxChi2 = 20;
    const double maxQSharing = 40;
    const double chEqReference = 450;

    const float qFrom = 50;
    const float qTo = 1050;
    const int qBins = 200;






    

    
#endif
