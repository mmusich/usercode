#ifndef SAMPLE_H
#define SAMPLE_H

#include <iostream> 

#include <TObject.h>
#include <TXMLNode.h>
#include <TXMLAttr.h>
#include <TFile.h>
#include <TH1F.h>
#include <TString.h>

using namespace std;

class Sample : public TObject {

public: 
  Sample() {} 
  Sample (Bool_t isMC, TString inputrootfile, 
	  TString histolabel, Int_t histocolor,
	  Float_t lumi, Float_t xsection) ;

  Bool_t isMC_;
  TString inputrootfile_;
  TString histolabel_;
  Int_t histocolor_;
  Float_t lumi_;
  Float_t xsection_;
  Float_t weight_;
  Float_t nevents_;
  void dump() { cout << "DUMP:" << endl;} 
  void clear() {
    isMC_=0;  inputrootfile_="";
    histolabel_=""; histocolor_=0;
    lumi_=0;  xsection_=0;
    weight_=0; nevents_=0;
  }
};

Sample::Sample(Bool_t isMC, TString inputrootfile, 
	       TString histolabel, Int_t histocolor,
	       Float_t lumi, Float_t xsection) : 
  isMC_(isMC), inputrootfile_(inputrootfile), histolabel_(histolabel), histocolor_(histocolor), lumi_(lumi), xsection_(xsection) {

  TFile f(inputrootfile_.Data(),"READ");
  nevents_ = ((TH1F*)f.Get("analyzePat/Selevents/SelectedEvts"))->GetBinContent(1);
  // now compute the weight
  if ( isMC ) {
    weight_ = lumi*xsection/nevents_;
  } else {
    weight_ = 1;
  }
  f.Close();
}

#endif
