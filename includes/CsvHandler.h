#pragma once



#include <glob.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "TObjString.h"



using namespace std;



class CsvHandler
{
public:

  CsvHandler() {};



  TString GetFirstFile(TString filtStr) {

    glob_t globbuf;
    glob(filtStr , 0, NULL, &globbuf);
    return globbuf.gl_pathv[0];
  }



  int Read(TString filename, char sep, double** data, int nrow, int ncol) {

    //double** data = new double*[nrow];
    string line, val; 
    ifstream inf(filename);
    if ( !inf.is_open() ) { cout << "Error: could not open " << filename << endl; return 1;} //if (!inf) //exit(EXIT_FAILURE);
    cout<< "Opened " << filename << endl;
 
    for(int row=0; row < nrow; row++) {      //inf >> line //while (getline(inf, line)) { 
      
      getline(inf, line); if (1) { cout<<"        parsed line["<<row<<"]:  "<<line<<endl; } 
      istringstream iss(line);
      //data[row] = new double[ncol];
      
      for (int col = 0; col < ncol; ++col) {
      
        getline(iss, val, sep); 
        data[row][col] = strtod(val.c_str(), NULL); if (0) { cout<<"        parsed token ["<<col<<"]: "<<data[row][col]<<endl; }
      }
    }
    return 0;
  }



  int Write(TString filename, char sep, double **vals, int nrow, int ncol, int precision) {

    ofstream outf(filename);
    if ( !outf.is_open() ) { cout << "Error: could not open " << filename << endl; return 1;}
    for (int irow = 0; irow < nrow; irow++) {
      for (int icol = 0; icol < ncol; icol++) {
        outf << std::setprecision(precision) << vals[irow][icol];
        if (icol != ncol - 1)
          outf << sep;
        else
          outf << endl;
      }
    }
    outf.close();
    return 0;
  }



  double ** InitMatrix(int nrow, int ncol) {

    double** arr = new double*[nrow];
    for(int i = 0; i < nrow; i++) {  
      arr[i] = new double[ncol]();
    }
    return arr;
  }












};




