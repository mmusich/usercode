#ifndef ZbbAnalysis_AnalysisStep_Zbbstruct4JEC_h
#define ZbbAnalysis_AnalysisStep_Zbbstruct4JEC_h

#include "TObject.h"
//Define the structures for JEC studies
const Int_t nMaxBjets = 6;

struct Event_info
{
  Int_t run_Event;
  Int_t lumisec_Event;
  Int_t num_Event;
  Int_t njet_Event;
  Int_t nbjet_Event;
  Int_t nvert_Event;
  Float_t weight_Event;
  Bool_t isMuon;
  void reset(){
    run_Event=0;lumisec_Event=0;num_Event=0;weight_Event=0;nvert_Event=0;
    njet_Event=0;nbjet_Event=0;
    isMuon=false;
  }
};

struct Z_info
{ 
  Float_t pT_Z;
  Float_t pX_Z;
  Float_t pY_Z;
  Float_t eta_Z;
  Float_t phi_Z;
  Float_t eta_d1;
  Float_t phi_d1;
  Float_t pt_d1;
  Float_t eta_d2;
  Float_t phi_d2;
  Float_t pt_d2;
  void reset(){
    pT_Z=0;pX_Z=0;pY_Z=0;eta_Z=0;phi_Z=0;eta_d1=0;phi_d1=0;pt_d1=0;eta_d2=0;phi_d2=0;pt_d2=0;
  }
};

struct bjet_info
{ 
  Float_t pT_bjet;
  Float_t eta_bjet;
  Float_t phi_bjet;
  Float_t dzDilept_bjet;
  Float_t frac1_bjet;
  Float_t frac2_bjet;
  Float_t frac3_bjet;
  Float_t frac2b_bjet;
  Float_t ssvhediscr_bjet; 
  Float_t ssvhpdiscr_bjet;
  Float_t csvdiscr_bjet;
  Float_t jpdiscr_bjet;
  Float_t tchediscr_bjet;
  Float_t tchpdiscr_bjet;
  Float_t svmass_bjet;
  Int_t   trueflavour_bjet;
  
  void reset(){
    pT_bjet=0;eta_bjet=0;phi_bjet=0;dzDilept_bjet=0;frac1_bjet=0;frac2_bjet=0;frac3_bjet=0;frac2b_bjet=0;
    ssvhediscr_bjet=-1; ssvhpdiscr_bjet=-1; csvdiscr_bjet=-1; jpdiscr_bjet=-1; tchediscr_bjet=-100; tchpdiscr_bjet=-100; svmass_bjet=-1;	   
    trueflavour_bjet=-9;
  }
};

struct btagjets_info
{
  Int_t   trueflavour_bjets[nMaxBjets];
  Float_t pT_bjets[nMaxBjets];
  Float_t eta_bjets[nMaxBjets];
  Float_t phi_bjets[nMaxBjets];
  Float_t ssvhediscr_bjets[nMaxBjets]; 
  Float_t ssvhpdiscr_bjets[nMaxBjets];
  Float_t csvdiscr_bjets[nMaxBjets];
  Float_t jpdiscr_bjets[nMaxBjets];   
  Float_t jpdiscr_fromtags_bjets[nMaxBjets];
  Float_t tchediscr_bjets[nMaxBjets];
  Float_t tchpdiscr_bjets[nMaxBjets];
  Float_t svmass_bjets[nMaxBjets];

  void reset(){
    for(Int_t i=0;i<nMaxBjets;i++){
      pT_bjets[i]=-10;eta_bjets[i]=-10;phi_bjets[i]=-10;
      ssvhediscr_bjets[i]=-1; ssvhpdiscr_bjets[i]=-1; csvdiscr_bjets[i]=-1; jpdiscr_bjets[i]=-1; jpdiscr_fromtags_bjets[i]=-1; svmass_bjets[i]=-1; tchediscr_bjets[i]=-100; tchpdiscr_bjets[i]=-100;	   	   
      trueflavour_bjets[i]=-10;
    }
  }  
};

struct jet2_info
{ 
  Float_t pT_jet2;
  Float_t eta_jet2;
  Float_t phi_jet2;
  Float_t dzDilept_jet2;
  Float_t frac1_jet2;
  Float_t frac2_jet2;
  Float_t frac3_jet2;
  Float_t frac2b_jet2;
  void reset(){
    pT_jet2=0;eta_jet2=0;phi_jet2=0;dzDilept_jet2=0;frac1_jet2=0;frac2_jet2=0;frac3_jet2=0;frac2b_jet2=0;
  }
};

struct MET_info
{
  Float_t pT_MET;
  Float_t px_MET;
  Float_t py_MET;
  void reset(){
    pT_MET=0;px_MET=0;py_MET=0;
  }
};

#endif
