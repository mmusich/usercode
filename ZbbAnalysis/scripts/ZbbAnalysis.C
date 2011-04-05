//////////////////////////////////////////////////////////
//
// Simple script to run ROOT analysis job:
//
// Usage:
//
// root -b
// ZbbAnalysis("DATA")     #for DATA
// ZbbAnalysis("MCZBB")    #for MC Zbb
// ZbbAnalysis("MCZCC")    #for MC Zcc
// ZbbAnalysis("MCDYJets") #for MC DY+Jets
// ZbbAnalysis("MCTTJets") #for MC TT+Jets
//
// Original Author: M. Musich INFN Torino
//
//////////////////////////////////////////////////////////

#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TCut.h"
#include "TTree.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
//#if !defined(__CINT__) && !defined(__MAKECINT__)                                                                                                          
#include <string>
#include <iostream>
//#endif   

void setGraphics(TH1F *histo,TString zedType);
Double_t deltaPhi(Double_t phi1, Double_t phi2);
Double_t deltaR(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2);
Int_t hasNoOverlap(Double_t etaJ, Double_t phiJ, Double_t phiLetpt1,Double_t etaLept1, Double_t phiLept2, Double_t etaLept2);
Int_t jetId(Double_t nhf,Double_t nef,Int_t nconstituents,Double_t chf,Int_t nch,Double_t cef);
Bool_t isBJet(Double_t theVariable,TString workingPoint,TString algo);
void ZbbAnalysis(TString _theSwitch);

void setGraphics(TH1F *histo,TString zedType){

  if(zedType=="ee"){
    histo->SetFillColor(kAzure+7);
    //histo->SetLineWidth(2);
    histo->SetLineColor(kBlue+1);
  } else if(zedType="mumu"){
    histo->SetFillColor(kOrange+7);
    //histo->SetLineWidth(2);
    histo->SetLineColor(kRed+1);
  }
}

// SOME helper functions

Double_t deltaPhi(Double_t phi1, Double_t phi2){
  // in CMSSW phi = [0,2pi], in TLorentzVector phi = [-pi,pi].
  // With the conversion below deltaPhi works ok despite the
  // 2*pi difference in phi definitions.
  Double_t PI = 3.1415;
  if(phi1 < 0) phi1 += 2*PI;
  if(phi2 < 0) phi2 += 2*PI;
  
  Double_t dphi = fabs(phi1-phi2);
  
  if(dphi > PI) dphi = 2*PI - dphi;
  return dphi;
}

Double_t deltaR(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2){
  // function calculates cone size
  Double_t dphi = deltaPhi(phi1,phi2);
  Double_t deta = fabs(eta1-eta2);
  return sqrt(dphi*dphi + deta*deta);
}

Int_t hasNoOverlap(Double_t etaJ, Double_t phiJ, Double_t phiLept1,Double_t etaLept1, Double_t phiLept2, Double_t etaLept2){
  
  Int_t theOverlap(0);

  Double_t dr1 = deltaR(etaJ,etaLept1,phiJ,phiLept1); 
  Double_t dr2 = deltaR(etaJ,etaLept2,phiJ,phiLept2); 
  if (dr1>0.5 and dr2>0.5){
    theOverlap=0;
  }
  else{
    theOverlap=1;
  }
  return theOverlap;
}

Int_t jetId(Double_t nhf,Double_t nEF,Int_t nconstituents,Double_t chf,Int_t nch,Double_t cef){

  Int_t theJetId = -9;
  //std::cout<<"nhf: "<<nhf<<" nEF: "<<nEF<<" nconstituents: "<<nconstituents<<" chf: "<<chf<<" nch: "<<nch<<" cef: "<<cef<<std::endl;

  if(nhf<0.90 && nEF<0.90 && nconstituents>1 && chf>0 && nch>0 && cef<0.99){
    //level=="tight"
    //std::cout<<"level==tight"<<endl;
    theJetId=2;
  } else if(nhf<0.95 && nEF<0.95 && nconstituents>1 && chf>0 && nch>0 && cef<0.99){
    //level=="medium"
    //std::cout<<"level==medium"<<endl;
    theJetId=1;
  } else if (nhf<0.99 && nEF<0.99 && nconstituents>1 && chf>0 && nch>0 && cef<0.99){
    //level=="loose"
    //std::cout<<"level==loose"<<endl;
    theJetId=0;
  } else {
    //std::cout<<"level==none"<<endl;
    theJetId=-1;
  }
  return theJetId; 
}

Bool_t isBJet(Double_t theVariable,TString workingPoint,TString algo){

  Bool_t theDiscr(false);

  if (algo=="SSV"){
    if( workingPoint=="HE" && theVariable>1.7){
      theDiscr=true;
    }
    else if(workingPoint=="HP" && theVariable>2.0){
      theDiscr=true;
    }
    else{
      std::cout<<"Error: unforeseen working point for b-tagging. Use HE or HP"<<std::endl;
    }
  } else if(algo=="TC"){
    if (workingPoint=="HE" && theVariable>3.3){
      theDiscr=true;
    }
    else if(workingPoint=="HP" && theVariable>3.41){
      theDiscr=true;
    }
    else{
      std::cout<<"Error: unforeseen working point for b-tagging. Use HE or HP"<<std::endl;
    }
  } else{
    std::cout<<"Error: unforeseen algo for b-tagging. Use SSV or TC"<<std::endl;
  }
  return theDiscr;
}

//  __  __      _         __              _   _          
// |  \/  |__ _(_)_ _    / _|_  _ _ _  __| |_(_)___ _ _  
// | |\/| / _` | | ' \  |  _| || | ' \/ _|  _| / _ \ ' \ 
// |_|  |_\__,_|_|_||_| |_|  \_,_|_||_\__|\__|_\___/_||_|
                                                      
void ZbbAnalysis(TString _theSwitch){

  TH1::StatOverflows(kTRUE);
  
  gStyle->SetOptStat();
  gROOT->SetStyle("Plain");
  using namespace std;

  TChain Events("Events"); 
  TFile * output_file(0); 

  if(_theSwitch.Contains("DATA")){
    Events.Add("rfio:/castor/cern.ch/user/m/musich/Zbb2010Ntuples/v4/ZbbAnalysisCandNtuples_DATA_Merged.root"); 
    output_file = TFile::Open("plots_DATA.root","RECREATE");
  } else if (_theSwitch.Contains("MCZBB") ){
    Events.Add("rfio:/castor/cern.ch/user/m/musich/Zbb2010Ntuples/v4/ZbbAnalysisCandNtuples_MC_ZbbToLL.root"); 
    output_file = TFile::Open("plots_MCZBB.root","RECREATE"); 
  } else if (_theSwitch.Contains("MCZCC")){
    Events.Add("rfio:/castor/cern.ch/user/m/musich/Zbb2010Ntuples/v4/ZbbAnalysisCandNtuples_MC_ZccToLL.root"); 
    output_file = TFile::Open("plots_MCZCC.root","RECREATE"); 
  } else if (_theSwitch.Contains("MCDYJets")){
    Events.Add("rfio:/castor/cern.ch/user/m/musich/Zbb2010Ntuples/v4/ZbbAnalysisCandNtuples_MC_DYJets.root"); 
    output_file = TFile::Open("plots_MCDYJets.root","RECREATE"); 
  } else if (_theSwitch.Contains("MCTTJets")){
    Events.Add("rfio:/castor/cern.ch/user/m/musich/Zbb2010Ntuples/v4/ZbbAnalysisCandNtuples_MC_ZccToLL.root"); 
    output_file = TFile::Open("plots_MCTTJets.root","RECREATE"); 
  }

  TDirectory * dir = output_file->mkdir("ZbbPlots");
  dir->cd();

  //*****************************************************************
  //  All these interesting plots are present in Louvain Analysis
  //  to be includeded someday!
  //*****************************************************************

//   TH1F *h_triggerSelection = new TH1F("triggerSelection","triggerSelection ",2,0,2);
//   TH1F *h_triggerBit = new TH1F("triggerBits","trigger bits",20,0,20);
//   TH1F *h_zmassMu = new TH1F("zmassMu","zmassMu",2000,0,200);
//   TH1F *h_massBestMu = new TH1F("bestzmassMu","bestzmassMu",2000,0,200);
//   TH1F *h_zmassEle = new TH1F("zmassEle","zmassEle",2000,0,200);
//   TH1F *h_massBestEle = new TH1F("bestzmassEle","bestzmassEle",2000,0,200);
//   TH1F *h_zptMu = new TH1F("zptMu","zptMu",500,0,500);
//   TH1F *h_ptBestMu = new TH1F("bestzptMu","bestzptMu",500,0,500);
//   TH1F *h_zptEle = new TH1F("zptEle","zptEle",500,0,500);
//   TH1F *h_ptBestEle = new TH1F("bestzptEle","bestzptEle",500,0,500);
//   TH1F *h_scaldptZbj1 = new TH1F("scaldptZbj1","scaldptZbj1",1000,-500,500);
//   TH1F *h_drZbj1 = new TH1F("drZbj1","distance between Z and leading jet",100,0,5);
//   TH1F *h_vecdptZbj1 = new TH1F("vecdptZbj1","vecdptZbj1",500,0,500);
//   TH1F *h_dphiZbj1 = new TH1F("dphiZbj1","dphiZbj1",40,0,4);
//   TH1F *h_dijetM = new TH1F("dijetM","b bbar invariant mass",1000,0,1000);
//   TH1F *h_dijetPt = new TH1F("dijetPt","b bbar Pt",500,0,500);
//   TH1F *h_ZbM = new TH1F("ZbM","Zb invariant mass",1000,0,1000);
//   TH1F *h_ZbPt = new TH1F("ZbPt","Zb Pt",500,0,500);
//   TH1F *h_ZbbM = new TH1F("ZbbM","Zbb invariant mass",1000,0,1000);
//   TH1F *h_ZbbPt = new TH1F("ZbbPt","Zbb Pt",500,0,500);
//   TH2F *h_ZbbM2D = new TH2F("ZbbM2D","Zbb mass vs bb mass",100,0,1000,100,0,1000);
//   TH1F *h_category = new TH1F("category","event category",10,0,10);  
//   TH1F *h_mu1pt = new TH1F("mu1pt","leading muon Pt",500,0,500);
//   TH1F *h_mu2pt = new TH1F("mu2pt","subleading muon Pt",500,0,500);
//   TH1F *h_mu1eta = new TH1F("mu1eta","leading muon Eta",25,0,2.5);
//   TH1F *h_mu2eta = new TH1F("mu2eta","subleading muon Eta",25,0,2.5);
//   TH1F *h_mu1etapm = new TH1F("mu1etapm","leading muon Eta",50,-2.5,2.5);
//   TH1F *h_mu2etapm = new TH1F("mu2etapm","subleading muon Eta",50,-2.5,2.5);
//   TH1F *h_el1pt = new TH1F("el1pt","leading electron Pt",500,0,500);
//   TH1F *h_el2pt = new TH1F("el2pt","subleading electron Pt",500,0,500);
//   TH1F *h_el1eta = new TH1F("el1eta","leading electron Eta",25,0,2.5);
//   TH1F *h_el2eta = new TH1F("el2eta","subleading electron Eta",25,0,2.5);
//   TH1F *h_el1etapm = new TH1F("el1etapm","leading electron Eta",50,-2.5,2.5);
//   TH1F *h_el2etapm = new TH1F("el2etapm","subleading electron Eta",50,-2.5,2.5);
//   TH1F *h_SSVHEdisc = new TH1F("SSVHEdisc","SSVHEdisc",200,-10,10);
//   TH1F *h_SSVHPdisc = new TH1F("SSVHPdisc","SSVHPdisc",200,-10,10);
//   TH1F *h_met = new TH1F("MET","MET",100,0,200);
//   TH1F *h_phimet = new TH1F("METphi","MET #phi",70,-3.5,3.5);
//   TH1F *h_jet1pt = new TH1F("jet1pt","leading jet Pt",500,15,515);
//   TH1F *h_jet1eta = new TH1F("jet1eta","leading jet Eta",25,0,2.5);
//   TH1F *h_jet1etapm = new TH1F("jet1etapm","leading jet Eta",50,-2.5,2.5);
//   TH1F *h_jet2pt = new TH1F("jet2pt","subleading jet Pt",500,15,515);
//   TH1F *h_jet2eta = new TH1F("jet2eta","subleading jet Eta",25,0,2.5);
//   TH1F *h_jet2etapm = new TH1F("jet2etapm","subleading jet Eta",50,-2.5,2.5);
//   TH1F *h_bjet1pt = new TH1F("bjet1pt","leading bjet Pt",500,15,515);
//   TH1F *h_bjet1eta = new TH1F("bjet1eta","leading bjet Eta",25,0,2.5);
//   TH1F *h_bjet1etapm = new TH1F("bjet1etapm","leading bjet Eta",50,-2.5,2.5);
//   TH1F *h_bjet2pt = new TH1F("bjet2pt","subleading bjet Pt",500,15,515);
//   TH1F *h_bjet2eta = new TH1F("bjet2eta","subleading bjet Eta",25,0,2.5);
//   TH1F *h_bjet2etapm = new TH1F("bjet2etapm","subleading bjet Eta",50,-2.5,2.5);
//   TH1F *h_nj = new TH1F("nj","jet count",15,0,15);
//   TH1F *h_nb = new TH1F("nb","b-jet count",5,0,5);
//   TH1F *h_nbP = new TH1F("nbP","pure b-jet count",5,0,5);
//   TH2I *h_njb = new TH2I("njb","number of bjets vs number of jets",15,0,15,5,0,5);
  
  //  __  __ ___ _____                _      _    _        
  // |  \/  | __|_   _| __ ____ _ _ _(_)__ _| |__| |___ ___
  // | |\/| | _|  | |   \ V / _` | '_| / _` | '_ \ / -_|_-<
  // |_|  |_|___| |_|    \_/\__,_|_| |_\__,_|_.__/_\___/__/
                                                      
  TH1F * h_MET = new TH1F("h_MET","MET;  #slash{E}_{T} [GeV]; events/GeV",100,0,100);
  TH1F * h_METPhi = new TH1F("h_METPhi","MET #phi; #phi of #slash{E}_{T}",80,-TMath::Pi(),TMath::Pi());
  TH1F * h_METEta = new TH1F("h_METEta","MET #eta; #eta of #slash{E}_{T}",100,-5,5);

  //     _     _    __   __        _      _    _        
  //  _ | |___| |_  \ \ / /_ _ _ _(_)__ _| |__| |___ ___
  // | || / -_)  _|  \ V / _` | '_| / _` | '_ \ / -_|_-<
  //  \__/\___|\__|   \_/\__,_|_| |_\__,_|_.__/_\___/__/
 

  // All jets variables 
  TH1F *h_jetpt = new TH1F("jetpt","Jet p_{T}; p_{T} of all jets [GeV];jets/GeV",100,15,215);
  TH1F *h_jeteta = new TH1F("jeteta","Jet #eta; #eta of all jets; jets/0.1",50,-2.5, 2.5);
  TH1F *h_jetphi = new TH1F("jetphi","Jet #phi; #phi of all jets [rad];",80,-TMath::Pi(),TMath::Pi());

  // Jet multiplicity 
  TH1F * JetMultiplicity = new TH1F("JetMultiplicity","Jet multiplicity; jet multiplicity; events",15,-0.5,14.5);
  
  // Leading Jet 
  TH1F *LeadingJetPt  = new TH1F("LeadingJetPt","Leading Jet p_{T};leading jet p_{T} [GeV];jets/GeV",100,15.,215.);
  TH1F *LeadingJetEta = new TH1F("LeadingJetEta","Leading Jet #eta;leading jet #eta; jets/0.1",50,-2.5, 2.5);
  TH1F *LeadingJetPhi = new TH1F("LeadingJetPhi","Leading Jet #phi;leading jet #phi",80,-TMath::Pi(),TMath::Pi());

  // Sub-Leading Jet   
  TH1F *SubLeadingJetPt  = new TH1F("SubLeadingJetPt","Sub-leading Jet p_{t};sub-leading jet p_{T} [GeV];jets/GeV",100,15.,215.);
  TH1F *SubLeadingJetEta = new TH1F("SubLeadingJetEta","Sub-leading Jet #eta;sub-leading jet #eta; jets/0.1",50,-2.5, 2.5);
  TH1F *SubLeadingJetPhi = new TH1F("SubLeadingJetPhi","Sub-leading Jet #phi;sub-leading jet #phi [rad]",80,-TMath::Pi(),TMath::Pi());

  // Jet b-tagging (SSVHP)
  TH1F * h_JetSSVHP = new TH1F("h_JetSSVHP","SSVHP discriminant; SSVHP discriminant; Jets",100,-5,5);
  // Jet b-tagging (SSVHE)
  TH1F * h_JetSSVHE = new TH1F("h_JetSSVHE","SSVHE discriminant; SSVHE discriminant; Jets",100,-5,5);

  // Jet ID variables
  TH1F *h_nhf = new TH1F("nhf","neutral hadron energy fraction;neutral hadron energy fraction",101,0,1.01);
  TH1F *h_nef = new TH1F("nef","neutral EmEnergy fraction;neutral EmEnergy fraction;",101,0,1.01);
  TH1F *h_nconstituents = new TH1F("npf","total multiplicity;total multiplicity",50,-0.5,49.5);
  TH1F *h_chf = new TH1F("chf","charged hadron energy fraction;charged hadron energy fraction",101,0,1.01);
  TH1F *h_nch = new TH1F("nch","charged multiplicity;charged multiplicity",50,-0.5,49.5);
  TH1F *h_cef = new TH1F("cef","charged EmEnergy fraction;charged EmEnergy fraction",101,0,1.01);
  TH1F *h_jetid = new TH1F("jetid","Jet Id level (none, loose, medium, tight);",4,-1.5,2.5);
  
  TString jetIdBinLabels[4] ={"none","loose","medium","tight"};

  for(UInt_t bin=1; bin<=4; bin++){
    h_jetid->GetXaxis()->SetBinLabel(bin,jetIdBinLabels[bin-1]);    
  }

  // Check for overlaps with leptons from Z
  TH1F *h_jetoverlapmu = new TH1F("jetoverlapmu","jets overlaps with muons; overlap bit",2,-0.5,1.5);
  TH1F *h_jetoverlapele = new TH1F("jetoverlapele","jets overlaps with electrons; overlap bit",2,-0.5,1.5);

  TString jetOverlapBinLabels[2] ={"false","true"};

  for(UInt_t bin=1; bin<=2; bin++){
    h_jetoverlapmu->GetXaxis()->SetBinLabel(bin,jetOverlapBinLabels[bin-1]);    
    h_jetoverlapele->GetXaxis()->SetBinLabel(bin,jetOverlapBinLabels[bin-1]);
  }
  
  //  ___  _ __ _  _   _ __ _  _  __ ____ _ _ _(_)__ _| |__| |___ ___
  // |_ / | '  \ || | | '  \ || | \ V / _` | '_| / _` | '_ \ / -_|_-<
  // /__| |_|_|_\_,_| |_|_|_\_,_|  \_/\__,_|_| |_\__,_|_.__/_\___/__/
  
  TH1F * ZMMCandMultiplicity = new TH1F("ZMMCandMultiplicity","Z #rightarrow #mu#mu Candidate multiplicity; Z #rightarrow #mu#mu candidate multiplicity; events",15,-0.5,14.5);

  TH1F * zmmMass = new TH1F("zmmMass","Z #rightarrow #mu#mu Mass; #mu^{+}#mu^{-} mass [GeV]; events/GeV",130,20.,150.);
  TH1F * zmmPt = new TH1F("zmmPt","Z #rightarrow #mu#mu p_{T}; #mu^{+}#mu^{-} p_{T} [GeV]; events/6 GeV",50,0,300);
  TH1F * zmmEta = new TH1F("zmmEta","Z #rightarrow #mu#mu #eta; #mu^{+}#mu^{-} #eta; events/0.14",100 ,-7,7);
  TH1F * zmmPhi = new TH1F("zmmPhi","Z #rightarrow #mu#mu #phi; #mu^{+}#mu^{-} #phi; events/0.08 rad",80,-TMath::Pi() ,TMath::Pi());
  TH1F * zmmY = new TH1F("zmmY","Z #rightarrow #mu#mu y; #mu^{+}#mu^{-} y; events/0.0625",80 ,-2.5,2.5);
  
  // daughters sorted by pt
  
  TH1F * zmmDau1Pt = new TH1F("zmmDau1Pt","leading #mu from Z p_{T};leading #mu p_{T} [GeV]; muons/5 GeV",50,0 ,250);
  TH1F * zmmDau2Pt = new TH1F("zmmDau2Pt","sub-leading #mu from Z p_{T};subleading #mu p_{T} [GeV];muons/5 GeV",50,0 ,250);
  TH1F * zmmDau1Eta= new TH1F("zmmDau1Eta","leading #mu from Z #eta;leading #mu #eta; muons/0.0625",80 ,-2.5, 2.5);
  TH1F * zmmDau1Phi= new TH1F("zmmDau1Phi","leading #mu from Z #phi;leading #mu #phi; muons/0.08 rad",80,-TMath::Pi(),TMath::Pi() );
  TH1F * zmmDau2Eta= new TH1F("zmmDau2Eta","sub-leading #mu from Z #eta;subleading #mu #eta;muons/0.0625;",80,-2.5, 2.5);
  TH1F * zmmDau2Phi= new TH1F("zmmDau2Phi","sub-leading #mu from Z #phi;leading #mu #phi; muons/0.08 rad",80,-TMath::Pi(),TMath::Pi() );
  TH1F * zmmDeltaR = new TH1F("zmmDeltaR","#Delta R of  #mu from Z;#DeltaR(#mu^{+}#mu^{-})",100,0,10); 

  //  ___  ___ ___  __ ____ _ _ _(_)__ _| |__| |___ ___
  // |_ / / -_) -_) \ V / _` | '_| / _` | '_ \ / -_|_-<
  // /__| \___\___|  \_/\__,_|_| |_\__,_|_.__/_\___/__/

  TH1F * ZEECandMultiplicity = new TH1F("ZEECandMultiplicity","Z #rightarrow ee Candidate multiplicity; Z #rightarrow ee candidate multiplicity; events",15,-0.5,14.5);

  TH1F * zeeMass = new TH1F("zeeMass","Z #rightarrow ee mass; e^{+}e^{-} mass [GeV];events/GeV",130,20.,150.);
  TH1F * zeePt = new TH1F("zeePt","Z #rightarrow ee p_{T}; e^{+}e^{-} p_{T} [GeV]; events/6 GeV",50,0 ,300);
  TH1F * zeeEta = new TH1F("zeeEta","Z #rightarrow ee #eta; e^{+}e^{-} #eta; events/0.14",100 ,-7,7);
  TH1F * zeePhi = new TH1F("zeePhi","Z #rightarrow ee #phi; e^{+}e^{-} #phi; events/0.08 rad",80,-TMath::Pi() ,TMath::Pi());
  TH1F * zeeY = new TH1F("zeeY","Z #rightarrow ee y; e^{+}e^{-} y; events/0.0625",80 ,-2.5,2.5);
  
  // daughters sorted by pt
  
  TH1F * zeeDau1Pt = new TH1F("zeeDau1Pt","leading e from Z p_{T};leading e p_{T} [GeV]; electrons/5 GeV",50,0 ,250);
  TH1F * zeeDau2Pt = new TH1F("zeeDau2Pt","sub-leading e from Z p_{T};subleading e p_{T} [GeV];electrons/5 GeV",50,0,250);
  TH1F * zeeDau1Eta= new TH1F("zeeDau1Eta","leading e from Z #eta;leading e #eta; electrons/0.00625",80 ,-2.5, 2.5);
  TH1F * zeeDau1Phi= new TH1F("zeeDau1Phi","leading e form Z #phi;leading e #phi; electrons/0.08 rad",80,-TMath::Pi(),TMath::Pi() );
  TH1F * zeeDau2Eta= new TH1F("zeeDau2Eta","sub-leading e from Z #eta;subleading e #eta;electrons/0.0625;",80,-2.5, 2.5);
  TH1F * zeeDau2Phi= new TH1F("zeeDau2Phi","sub-leading e from Z #phi;subleading e #phi;electrons/0.08 rad",80,-TMath::Pi(),TMath::Pi() );
  TH1F * zeeDeltaR = new TH1F("zeeDeltaR","#Delta R of e from Z;#DeltaR(e^{+}e^{-})",100,0,10); 

//  _ _    _          __           _      
// | (_)__| |_   ___ / _|  __ _  _| |_ ___
// | | (_-<  _| / _ \  _| / _| || |  _(_-<
// |_|_/__/\__| \___/_|   \__|\_,_|\__/__/
  
  TCut selection_jets("floats_JetNtuplizer_JetEta_TEST.@obj.size()>1 && JetPt>25 && TMath::Abs(JetEta)<2.1");

  TCut selection_mu("(TMath::Abs(ZmmLeptDau1Eta)<2.1 && TMath::Abs(ZmmLeptDau2Eta<2.1)) && (ZmmLeptDau1Pt > 20 &&  ZmmLeptDau1Pt > 20) && (ZmmLeptDau1CombRelIso<0.15 && ZmmLeptDau2CombRelIso<0.15)"+selection_jets);
  
  TCut selection_e("(TMath::Abs(ZeeLeptDau1Eta)<2.5 && TMath::Abs(ZeeLeptDau2Eta)<2.5) && (ZeeLeptDau1Pt > 25 &&  ZeeLeptDau1Pt > 25) && (ZeeLeptDau1CombRelIso<0.15 && ZeeLeptDau2CombRelIso<0.15)"+selection_jets);

//     _       _   _               _     _   _   _           
//  __| |___  | |_| |_  ___   _ __| |___| |_| |_(_)_ _  __ _ 
// / _` / _ \ |  _| ' \/ -_) | '_ \ / _ \  _|  _| | ' \/ _` |
// \__,_\___/  \__|_||_\___| | .__/_\___/\__|\__|_|_||_\__, |
//                           |_|                       |___/ 

  //Zmm
  
  Events.Draw("floats_ZMMNtuplizer_ZmmEta_TEST.@obj.size()>>ZMMCandMultiplicity",selection_mu);

  Events.Project("zmmMass","ZmmMass",selection_mu);
  Events.Project("zmmPt","ZmmPt",selection_mu); 
  Events.Project("zmmEta","ZmmEta",selection_mu); 
  Events.Project("zmmPhi","ZmmPhi",selection_mu); 
  Events.Project("zmmY","ZmmY",selection_mu); 
  
  Events.Project("zmmDau1Pt","ZmmLeptDau1Pt",selection_mu);
  Events.Project("zmmDau2Pt","ZmmLeptDau2Pt",selection_mu);
  Events.Project("zmmDau1Eta","ZmmLeptDau1Eta",selection_mu);
  Events.Project("zmmDau2Eta","ZmmLeptDau2Eta",selection_mu);
  Events.Project("zmmDau1Phi","ZmmLeptDau1Phi",selection_mu);
  Events.Project("zmmDau2Phi","ZmmLeptDau2Phi",selection_mu);
  Events.Project("zmmDeltaR","deltaR(ZmmLeptDau1Eta,ZmmLeptDau2Eta,ZmmLeptDau1Phi,ZmmLeptDau2Phi)",selection_mu);
  
  // Zee

  Events.Draw("floats_ZEENtuplizer_ZeeEta_TEST.@obj.size()>>ZEECandMultiplicity",selection_e);

  Events.Project("zeeMass","ZeeMass",selection_e);
  Events.Project("zeeMass","ZeeMass",selection_e);
  Events.Project("zeePt","ZeePt",selection_e); 
  Events.Project("zeeEta","ZeeEta",selection_e); 
  Events.Project("zeePhi","ZeePhi",selection_e); 
  Events.Project("zeeY","ZeeY",selection_e); 
  
  Events.Project("zeeDau1Pt","ZeeLeptDau1Pt",selection_e);
  Events.Project("zeeDau2Pt","ZeeLeptDau2Pt",selection_e);
  Events.Project("zeeDau1Eta","ZeeLeptDau1Eta",selection_e);
  Events.Project("zeeDau2Eta","ZeeLeptDau2Eta",selection_e);
  Events.Project("zeeDau1Phi","ZeeLeptDau1Phi",selection_e);
  Events.Project("zeeDau2Phi","ZeeLeptDau2Phi",selection_e);
  Events.Project("zeeDeltaR","deltaR(ZeeLeptDau1Eta,ZeeLeptDau2Eta,ZeeLeptDau1Phi,ZeeLeptDau2Phi)",selection_e);

  cout<<"Number of Z -> mm candidate : "<<zmmMass->GetEntries()<<endl;
  cout<<"Number of Z -> ee candidate : "<<zeeMass->GetEntries()<<endl;
  
  //    _     _                  _      _    _
  //   (_)___| |_  __ ____ _ _ _(_)__ _| |__| |___ ___
  //   | / -_)  _| \ V / _` | '_| / _` | '_ \ / -_|_-<
  //  _/ \___|\__|  \_/\__,_|_| |_\__,_|_.__/_\___/__/
  // |__/                                             
  
  Events.Project("jetpt","JetPt",selection_jets);
  Events.Project("jeteta","JetEta",selection_jets); 
  Events.Project("jetphi","JetPhi",selection_jets); 
  			 
  Events.Project("LeadingJetPt","JetPt[0]",selection_jets);
  Events.Project("LeadingJetEta","JetEta[0]",selection_jets); 
  Events.Project("LeadingJetPhi","JetPhi[0]",selection_jets); 
  
  Events.Project("SubLeadingJetPt","JetPt[1]",selection_jets);
  Events.Project("SubLeadingJetEta","JetEta[1]",selection_jets); 
  Events.Project("SubLeadingJetPhi","JetPhi[1]",selection_jets); 		 

  Events.Draw("floats_JetNtuplizer_JetEta_TEST.@obj.size()>>JetMultiplicity",selection_jets);
  
  Events.Project("h_JetSSVHP","JetJetSSVHP",selection_jets);   
  Events.Project("h_JetSSVHE","JetJetSSVHE",selection_jets); 

  Events.Project("nhf","JetJetNhf",selection_jets);  
  Events.Project("nef","JetJetNef",selection_jets);   
  Events.Project("npf","JetJetNconstituents",selection_jets);   
  Events.Project("chf","JetJetChf",selection_jets);   
  Events.Project("nch","JetJetNch",selection_jets);   
  Events.Project("cef","JetJetCef",selection_jets);   
  Events.Project("jetid","jetId(JetJetNhf,JetJetNef,JetJetNconstituents,JetJetChf,JetJetNch,JetJetCef)",selection_jets);  

  Events.Project("jetoverlapmu","hasNoOverlap(JetEta,JetPhi,ZmmLeptDau1Phi,ZmmLeptDau1Eta,ZmmLeptDau2Phi,ZmmLeptDau2Eta)",selection_jets);
  Events.Project("jetoverlapele","hasNoOverlap(JetEta,JetPhi,ZeeLeptDau1Phi,ZeeLeptDau1Eta,ZeeLeptDau2Phi,ZeeLeptDau2Eta)",selection_jets);

//  __  __ ___ _____                _      _    _        
// |  \/  | __|_   _| __ ____ _ _ _(_)__ _| |__| |___ ___
// | |\/| | _|  | |   \ V / _` | '_| / _` | '_ \ / -_|_-<
// |_|  |_|___| |_|    \_/\__,_|_| |_\__,_|_.__/_\___/__/
                                                      
  Events.Project("h_MET","pfMetPt",selection_jets); 
  Events.Project("h_METPhi","pfMetPhi",selection_jets); 
  Events.Project("h_METEta","pfMetEta",selection_jets); 

//     _                                _                    
//  __| |_ _ __ ___ __ __  __ _ _ _  __| |  ___ __ ___ _____ 
// / _` | '_/ _` \ V  V / / _` | ' \/ _` | (_-</ _` \ V / -_)
// \__,_|_| \__,_|\_/\_/  \__,_|_||_\__,_| /__/\__,_|\_/\___|
                                                          
  //Zmm

  ZMMCandMultiplicity->Draw();
  ZMMCandMultiplicity->Write();
  
  zmmMass->Draw(); 
  zmmMass->Write();
  
  zmmPt->Draw(); 
  zmmPt->Write(); 
  
  zmmEta->Draw(); 
  zmmEta->Write(); 
  
  zmmPhi->Draw(); 
  zmmPhi->Write(); 
  
  zmmY->Draw(); 
  zmmY->Write(); 
  
  zmmDeltaR->Draw(); 
  zmmDeltaR->Write(); 
 
  zmmDau1Pt->Draw(); 
  zmmDau1Pt->Write(); 
  
  zmmDau2Pt->Draw(); 
  zmmDau2Pt->Write(); 
  
  zmmDau1Eta->Draw(); 
  zmmDau1Eta->Write(); 
  
  zmmDau1Phi->Draw(); 
  zmmDau1Phi->Write(); 
  
  zmmDau2Eta->Draw(); 
  zmmDau2Eta->Write(); 
  
  zmmDau2Phi->Draw(); 
  zmmDau2Phi->Write();
   
  //Zee
  
  ZEECandMultiplicity->Draw();
  ZEECandMultiplicity->Write();
 
  zeeMass->Draw(); 
  zeeMass->Write();
  
  zeePt->Draw();
  zeePt->Write(); 
  
  zeeEta->Draw(); 
  zeeEta->Write(); 
  
  zeePhi->Draw(); 
  zeePhi->Write(); 
  
  zeeY->Draw(); 
  zeeY->Write(); 
  
  zeeDeltaR->Draw(); 
  zeeDeltaR->Write(); 

  zeeDau1Pt->Draw(); 
  zeeDau1Pt->Write(); 
  
  zeeDau2Pt->Draw(); 
  zeeDau2Pt->Write(); 
  
  zeeDau1Eta->Draw(); 
  zeeDau1Eta->Write();
  
  zeeDau1Phi->Draw(); 
  zeeDau1Phi->Write(); 
  
  zeeDau2Eta->Draw(); 
  zeeDau2Eta->Write(); 
  
  zeeDau2Phi->Draw(); 
  zeeDau2Phi->Write(); 

  // MET

  h_MET->Draw();
  h_MET->Write();

  h_METPhi->Draw();
  h_METPhi->Write();

  h_METEta->Draw();
  h_METEta->Write();

  // Leading Jet

  LeadingJetPt->Draw(); 
  LeadingJetPt->Write(); 
  
  LeadingJetEta->Draw();
  LeadingJetEta->Write();

  LeadingJetPhi->Draw();
  LeadingJetPhi->Write();

  // SubLeading Jet

  SubLeadingJetPt->Draw(); 
  SubLeadingJetPt->Write(); 

  SubLeadingJetEta->Draw();
  SubLeadingJetEta->Write();

  SubLeadingJetPhi->Draw();
  SubLeadingJetPhi->Write();
  
  JetMultiplicity->Draw(); 
  JetMultiplicity->Write(); 
        
  // all jets

  h_JetSSVHP->Draw();
  h_JetSSVHP->Write();

  h_JetSSVHE->Draw();
  h_JetSSVHE->Write(); 
  
  h_jetpt->Draw();
  h_jetpt->Write();

  h_jeteta->Draw();
  h_jeteta->Write();

  h_jetphi->Draw();
  h_jetphi->Write(); 
    				    			        				  
  h_nhf->Draw();
  h_nhf->Write();

  h_nef->Draw();
  h_nef->Write();

  h_nconstituents->Draw();
  h_nconstituents->Write();

  h_chf->Draw();
  h_chf->Write();
  
  h_nch->Draw();
  h_nch->Write();

  h_cef->Draw();
  h_cef->Write();

  h_jetid->Draw();
  h_jetid->Write();

  h_jetoverlapmu->Draw();
  h_jetoverlapmu->Write();

  h_jetoverlapele->Draw();
  h_jetoverlapele->Write();

  //close outputfile

  output_file->Close();

}
