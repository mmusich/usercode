//////////////////////////////////////////////////////////
//
// Simple script to run ROOT analysis job:
//
// Usage:
//
// root -b
// .x ZBB("ele") or .x ZBB("mu")
//
// Original Author: M. Musich INFN Torino
//
//////////////////////////////////////////////////////////

#include "TFile.h"
#include "TH1F.h"
#include "TCut.h"
#include "TTree.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#if !defined(__CINT__) && !defined(__MAKECINT__)                                                                                                          
#include <string>
#include <iostream>
#endif   


void setGraphics(TH1F *histo,TString zedType){

  if(zedType=="ee"){
    histo->SetFillColor(kAzure+7);
    histo->SetLineWidth(2);
    histo->SetLineColor(kBlue+1);
  } else if(zedType="mumu"){
    histo->SetFillColor(kOrange+7);
    histo->SetLineWidth(2);
    histo->SetLineColor(kRed+1);
  }
}


void ZBB(TString _theSwitch){

  TH1::StatOverflows(kTRUE);
  
  gStyle->SetOptStat();
  gROOT->SetStyle("Plain");
  using namespace std;

  TChain Events("Events"); 
  TFile * output_file; 

  //Events.Add("rfio:/castor/cern.ch/user/m/musich/Zbb2010Ntuples/ZJJCandNtuples_Mu.root", 65);
  //Events.Add("rfio:/castor/cern.ch/user/m/musich/Zbb2010Ntuples/ZJJCandNtuples_Ele.root", 65);

  if(_theSwitch.Contains("ele")){
    Events.Add("/tmp/musich/ZJJCandNtuples_Ele_CleanedJets.root"); 
    Events.Add("/tmp/musich/ZJJCandNtuples_Egamma_CleanedJets.root"); 
    output_file = TFile::Open("plots_datasetEle.root", "RECREATE");
  } else if (_theSwitch.Contains("mu") ){
    Events.Add("/tmp/musich/ZJJCandNtuples_Mu_CleanedJets.root", 65);
    output_file = TFile::Open("plots_datasetMu.root", "RECREATE"); 
  } else if (_theSwitch.Contains("egamma")){
    Events.Add("/tmp/musich/ZJJCandNtuples_Egamma_CleanedJets.root", 65); 
    output_file = TFile::Open("plots_datasetEgamma.root", "RECREATE"); 
  }

  TDirectory * dir = output_file->mkdir("ZPlots");
  dir->cd();

  //  ___  _ __ _  _   _ __ _  _  __ ____ _ _ _(_)__ _| |__| |___ ___
  // |_ / | '  \ || | | '  \ || | \ V / _` | '_| / _` | '_ \ / -_|_-<
  // /__| |_|_|_\_,_| |_|_|_\_,_|  \_/\__,_|_| |_\__,_|_.__/_\___/__/
  
  TH1F * zmmMass = new TH1F("zmmMass","Z #rightarrow #mu#mu Mass; #mu^{+}#mu^{-} mass [GeV]; events/GeV",130,20.,150.);
  TH1F * zmmPt = new TH1F("zmmPt","Z #rightarrow #mu#mu p_{T}; #mu^{+}#mu^{-} p_{T} [GeV]; events/3 GeV",100,0 ,300);
  TH1F * zmmEta = new TH1F("zmmEta","Z #rightarrow #mu#mu #eta; #mu^{+}#mu^{-} #eta; events/0.14",100 ,-7,7);
  TH1F * zmmPhi = new TH1F("zmmPhi","Z #rightarrow #mu#mu #phi; #mu^{+}#mu^{-} #phi; events/0.08 rad",80,-TMath::Pi() ,TMath::Pi());
  TH1F * zmmY = new TH1F("zmmY","Z #rightarrow #mu#mu y; #mu^{+}#mu^{-} y; events/0.0625",80 ,-2.5,2.5);
  
  // daughters sorted by pt
  
  TH1F * zmmDau1Pt = new TH1F("zmmDau1Pt","leading #mu from Z p_{T};leading #mu p_{T} [GeV]; muons/2.5 GeV",100 ,0 ,250);
  TH1F * zmmDau2Pt = new TH1F("zmmDau2Pt","sub-leading #mu from Z p_{T};subleading #mu p_{T} [GeV];muons/2.5 GeV",100 ,0 ,250);
  TH1F * zmmDau1Eta= new TH1F("zmmDau1Eta","leading #mu from Z #eta;leading #mu #eta; muons/0.0625",80 ,-2.5, 2.5);
  TH1F * zmmDau1Phi= new TH1F("zmmDau1Phi","leading #mu from Z #phi;leading #mu #phi; muons/0.08 rad",80,-TMath::Pi(),TMath::Pi() );
  TH1F * zmmDau2Eta= new TH1F("zmmDau2Eta","sub-leading #mu from Z #eta;subleading #mu #eta;muons/0.0625;",80,-2.5, 2.5);
  TH1F * zmmDau2Phi= new TH1F("zmmDau2Phi","sub-leading #mu from Z #phi;leading #mu #phi; muons/0.08 rad",80,-TMath::Pi(),TMath::Pi() );
  
  // daughters sorted by charge
  
  TH1F * zmmDauplusPt = new  TH1F("zmmDauplusPt",  "#mu^{+} p_{T}; #mu^{+} p_{T} [GeV]; muons/2.5 GeV",100 ,0 ,250);
  TH1F * zmmDauMinusPt = new TH1F("zmmDauMinusPt", "#mu^{-} p_{T};#mu^{-} p_{T} [GeV];muons/ 2.5 GeV ",100 ,0 ,250);
  TH1F * zmmDauplusEta= new  TH1F("zmmDauplusEta", "#mu^{+} #eta;#mu^{+} #eta; muons/0.05",100 ,-2.5, 2.5);
  TH1F * zmmDauplusPhi= new  TH1F("zmmDauplusPhi", "#mu^{-} #phi;#mu^{-} #phi; muons",100,-TMath::Pi(),TMath::Pi() );
  TH1F * zmmDauMinusEta= new TH1F("zmmDauMinusEta","#mu^{+} #eta;#mu^{+} #eta;muons/0.05;",100,-2.5, 2.5);
  TH1F * zmmDauMinusPhi= new TH1F("zmmDauMinusPhi","#mu^{-} #phi;#mu^{-} #phi; muons",100,-TMath::Pi(),TMath::Pi() );
  
  TH1F * zmmDauplusPtclone = new  TH1F("zmmDauplusPtclone",  ";#mu^{+} p_{T};#mu^{+} p_{T} [GeV]; muons/2.5 GeV",100 ,0 ,250);
  TH1F * zmmDauMinusPtclone = new TH1F("zmmDauMinusPtclone", ";#mu^{-} p_{T};#mu^{-} p_{T} [GeV];muons/ 2.5 GeV",100 ,0 ,250);
  TH1F * zmmDauplusEtaclone= new  TH1F("zmmDauplusEtaclone", ";#mu^{+} #eta;#mu^{+} #eta; muons/0.05",100 ,-2.5, 2.5);
  TH1F * zmmDauplusPhiclone= new  TH1F("zmmDauplusPhiclone", ";#mu^{-} #phi;#mu^{-} #phi; muons",100,-TMath::Pi(),TMath::Pi() );
  TH1F * zmmDauMinusEtaclone= new TH1F("zmmDauMinusEtaclone",";#mu^{+} #eta;#mu^{+} #eta;muons/0.05;",100,-2.5, 2.5);
  TH1F * zmmDauMinusPhiclone= new TH1F("zmmDauMinusPhiclone",";#mu^{-} #phi;#mu^{-} #phi; muons",100,-TMath::Pi(),TMath::Pi() );
 
  //  ___  ___ ___  __ ____ _ _ _(_)__ _| |__| |___ ___
  // |_ / / -_) -_) \ V / _` | '_| / _` | '_ \ / -_|_-<
  // /__| \___\___|  \_/\__,_|_| |_\__,_|_.__/_\___/__/
  
  TH1F * zeeMass = new TH1F("zeeMass","Z #rightarrow ee mass; e^{+}e^{-} mass [GeV];events/GeV",130,20.,150.);
  TH1F * zeePt = new TH1F("zeePt","Z #rightarrow ee p_{T}; e^{+}e^{-} p_{T} [GeV]; events/3 GeV",100 ,0 ,300);
  TH1F * zeeEta = new TH1F("zeeEta","Z #rightarrow ee #eta; e^{+}e^{-} #eta; events/0.14",100 ,-7,7);
  TH1F * zeePhi = new TH1F("zeePhi","Z #rightarrow ee #phi; e^{+}e^{-} #phi; events/0.08 rad",80,-TMath::Pi() ,TMath::Pi());
  TH1F * zeeY = new TH1F("zeeY","Z #rightarrow ee y; e^{+}e^{-} y; events/0.0625",80 ,-2.5,2.5);
  
  // daughters sorted by pt
  
  TH1F * zeeDau1Pt = new TH1F("zeeDau1Pt","leading e from Z p_{T};leading e p_{T} [GeV]; electrons/2.5 GeV",100 ,0 ,250);
  TH1F * zeeDau2Pt = new TH1F("zeeDau2Pt","sub-leading e from Z p_{T};subleading e p_{T} [GeV];electrons/2.5 GeV",100 ,0 ,250);
  TH1F * zeeDau1Eta= new TH1F("zeeDau1Eta","leading e from Z #eta;leading e #eta; electrons/0.00625",80 ,-2.5, 2.5);
  TH1F * zeeDau1Phi= new TH1F("zeeDau1Phi","leading e form Z #phi;leading e #phi; electrons/0.08 rad",80,-TMath::Pi(),TMath::Pi() );
  TH1F * zeeDau2Eta= new TH1F("zeeDau2Eta","sub-leading e from Z #eta;subleading e #eta;electrons/0.0625;",80,-2.5, 2.5);
  TH1F * zeeDau2Phi= new TH1F("zeeDau2Phi","sub-leading e from Z #phi;subleading e #phi;electrons/0.08 rad",80,-TMath::Pi(),TMath::Pi() );
  
  // daughters sorted by charge
  
  TH1F * zeeDauplusPt = new TH1F("zeeDauplusPt","e^{+} p_{T};e^{+} p_{T} [GeV]; electrons/2.5 GeV",100 ,0 ,250);
  TH1F * zeeDauMinusPt = new TH1F("zeeDauMinusPt","e^{-} p_{T};e^{-} p_{T} [GeV];electrons/ 2.5 GeV",100 ,0 ,250);
  TH1F * zeeDauplusEta= new TH1F("zeeDauplusEta","e^{+} #eta;e^{+} #eta; electrons/0.05",100 ,-2.5, 2.5);
  TH1F * zeeDauplusPhi= new TH1F("zeeDauplusPhi","e^{-} #phi;e^{-} #phi; electrons",100,-TMath::Pi(),TMath::Pi() );
  TH1F * zeeDauMinusEta= new TH1F("zeeDauMinusEta","e^{+} #eta;e^{+} #eta;electrons/0.05;",100,-2.5, 2.5);
  TH1F * zeeDauMinusPhi= new TH1F("zeeDauMinusPhi","e^{-} #phi;e^{-} #phi; electrons",100,-TMath::Pi(),TMath::Pi() );

  TH1F * zeeDauplusPtclone = new TH1F("zeeDauplusPtclone","e^{+} p_{T};e^{+} p_{T} [GeV]; electrons/2.5 GeV",100 ,0 ,250);
  TH1F * zeeDauMinusPtclone = new TH1F("zeeDauMinusPtclone","e^{-} p_{T};e^{-} p_{T} [GeV];electrons/ 2.5 GeV",100 ,0 ,250);
  TH1F * zeeDauplusEtaclone= new TH1F("zeeDauplusEtaclone","e^{+} #eta;e^{+} #eta; electrons/0.05",100 ,-2.5, 2.5);
  TH1F * zeeDauplusPhiclone= new TH1F("zeeDauplusPhiclone","e^{-} #phi;e^{-} #phi; electrons",100,-TMath::Pi(),TMath::Pi() );
  TH1F * zeeDauMinusEtaclone= new TH1F("zeeDauMinusEtaclone","e^{+} #eta;e^{+} #eta;electrons/0.05;",100,-2.5, 2.5);
  TH1F * zeeDauMinusPhiclone= new TH1F("zeeDauMinusPhiclone","e^{-} #phi;e^{-} #phi; electrons",100,-TMath::Pi(),TMath::Pi() );

  TCut goldenCutOnNumberOjJets("floats_JetNtuplizer_JetEta_TEST.@obj.size()>1");

  TCut selection_mu("ZmmLeptDau1Eta<5 && ZmmLeptDau2Eta<5"+goldenCutOnNumberOjJets);

  TCut selection_chargemuplus1("ZmmLeptDau1Q>0"+goldenCutOnNumberOjJets);
  TCut selection_chargemuplus2("ZmmLeptDau2Q>0"+goldenCutOnNumberOjJets);
  TCut selection_chargemuminus1("ZmmLeptDau1Q<0"+goldenCutOnNumberOjJets);
  TCut selection_chargemuminus2("ZmmLeptDau2Q<0"+goldenCutOnNumberOjJets);

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

  //histograms by charge

  Events.Project("zmmDauplusPt","ZmmLeptDau1Pt",selection_chargemuplus1);
  Events.Project("zmmDauplusPtclone","ZmmLeptDau2Pt",selection_chargemuplus2);

  zmmDauplusPt->Add(zmmDauplusPtclone);
 
  Events.Project("zmmDauMinusPt","ZmmLeptDau1Pt",selection_chargemuminus1);
  Events.Project("zmmDauMinusPtclone","ZmmLeptDau2Pt",selection_chargemuminus2);

  zmmDauMinusPt->Add(zmmDauMinusPtclone);

  Events.Project("zmmDauplusEta","ZmmLeptDau1Eta",selection_chargemuplus1);
  Events.Project("zmmDauplusEtaclone","ZmmLeptDau2Eta",selection_chargemuplus2);

  zmmDauplusEta->Add(zmmDauplusEtaclone);

  Events.Project("zmmDauMinusEta","ZmmLeptDau1Phi",selection_chargemuminus1);
  Events.Project("zmmDauMinusEtaclone","ZmmLeptDau2Phi",selection_chargemuminus2);  

  zmmDauMinusEta->Add(zmmDauMinusEtaclone);

  Events.Project("zmmDauplusPhi","ZmmLeptDau1Eta",selection_chargemuplus1);
  Events.Project("zmmDauplusPhiclone","ZmmLeptDau2Eta",selection_chargemuplus2);
  
  zmmDauplusPhi->Add(zmmDauplusPhiclone);

  Events.Project("zmmDauMinusPhi","ZmmLeptDau1Phi",selection_chargemuminus1);
  Events.Project("zmmDauMinusPhiclone","ZmmLeptDau2Phi",selection_chargemuminus2); 

  zmmDauMinusPhi->Add(zmmDauMinusPhiclone);
  
  TCut selection_e("ZeeLeptDau1Eta<5 && ZeeLeptDau2Eta<5"+goldenCutOnNumberOjJets);
  TCut selection_chargeeleplus1("ZeeLeptDau1Q>0"+goldenCutOnNumberOjJets);
  TCut selection_chargeeleplus2("ZeeLeptDau2Q>0"+goldenCutOnNumberOjJets);
  TCut selection_chargeeleminus1("ZeeLeptDau1Q<0"+goldenCutOnNumberOjJets);
  TCut selection_chargeeleminus2("ZeeLeptDau2Q<0"+goldenCutOnNumberOjJets);

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

  //histograms by charge

  Events.Project("zeeDauplusPt","ZeeLeptDau1Pt",selection_chargeeleplus1);
  Events.Project("zeeDauplusPtclone","ZeeLeptDau2Pt",selection_chargeeleplus2);

  zeeDauplusPt->Add(zeeDauplusPtclone);
 
  Events.Project("zeeDauMinusPt","ZeeLeptDau1Pt",selection_chargeeleminus1);
  Events.Project("zeeDauMinusPtclone","ZeeLeptDau2Pt",selection_chargeeleminus2);

  zeeDauMinusPt->Add(zeeDauMinusPtclone);

  Events.Project("zeeDauplusEta","ZeeLeptDau1Eta",selection_chargeeleplus1);
  Events.Project("zeeDauplusEtaclone","ZeeLeptDau2Eta",selection_chargeeleplus2);

  zeeDauplusEta->Add(zeeDauplusEtaclone);

  Events.Project("zeeDauMinusEta","ZeeLeptDau1Phi",selection_chargeeleminus1);
  Events.Project("zeeDauMinusEtaclone","ZeeLeptDau2Phi",selection_chargeeleminus2);  

  zeeDauMinusEta->Add(zeeDauMinusEtaclone);

  Events.Project("zeeDauplusPhi","ZeeLeptDau1Eta",selection_chargeeleplus1);
  Events.Project("zeeDauplusPhiclone","ZeeLeptDau2Eta",selection_chargeeleplus2);
  
  zeeDauplusPhi->Add(zeeDauplusPhiclone);

  Events.Project("zeeDauMinusPhi","ZeeLeptDau1Phi",selection_chargeeleminus1);
  Events.Project("zeeDauMinusPhiclone","ZeeLeptDau2Phi",selection_chargeeleminus2); 

  zeeDauMinusPhi->Add(zeeDauMinusPhiclone);

  cout<<"Number of Z -> mm candidate : "<<zmmMass->GetEntries()<<endl;
  cout<<"Number of Z -> ee candidate : "<<zeeMass->GetEntries()<<endl;
  
  //    _     _                  _      _    _
  //   (_)___| |_  __ ____ _ _ _(_)__ _| |__| |___ ___
  //   | / -_)  _| \ V / _` | '_| / _` | '_ \ / -_|_-<
  //  _/ \___|\__|  \_/\__,_|_| |_\__,_|_.__/_\___/__/
  // |__/                                             
  
  TH1F * LeadingJetPt = new TH1F("LeadingJetPt","Leading Jet p_{T};leading jet p_{T} [GeV];jets/6GeV",50,0.,300.);
  Events.Project("LeadingJetPt","JetPt[0]",goldenCutOnNumberOjJets);
  
  TH1F * SubLeadingJetPt = new TH1F("SubLeadingJetPt","Sub-leading Jet p_{t};sub-leading jet p_{T} [GeV];jets/6GeV",50,0.,300.);
  Events.Project("SubLeadingJetPt","JetPt[1]",goldenCutOnNumberOjJets);
  
  // Check on dijet invariant mass
  TH1F * dijetMass = new TH1F("dijetMass","jj Mass; jj mass [GeV];events/10GeV",100,0.,1000.);
  Events.Project("dijetMass","TMath::Sqrt((Jetenergy[0]+Jetenergy[1])*(Jetenergy[0]+Jetenergy[1]) - (Jetpx[0]+Jetpx[1])*(Jetpx[0]+Jetpx[1]) -  (Jetpy[0]+Jetpy[1])*(Jetpy[0]+Jetpy[1]) -  (Jetpz[0]+Jetpz[1])*(Jetpz[0]+Jetpz[1]))","floats_JetNtuplizer_JetEta_TEST.@obj.size()==2");
  
  TH1F * diElectronMass = new TH1F("diElectronMass","ee Mass; ee mass [GeV];events/GeV",130,20.,150.);
  Events.Project("diElectronMass","TMath::Sqrt((ZeeLeptDau1energy+ZeeLeptDau2energy)*(ZeeLeptDau1energy+ZeeLeptDau2energy) - (ZeeLeptDau1px+ZeeLeptDau2px)*(ZeeLeptDau1px+ZeeLeptDau2px) -  (ZeeLeptDau1py+ZeeLeptDau2py)*(ZeeLeptDau1py+ZeeLeptDau2py) -  (ZeeLeptDau1pz+ZeeLeptDau2pz)*(ZeeLeptDau1pz+ZeeLeptDau2pz))","");
  
  // Check on leading lepton-jet delta phi
  
  TH1F * DeltaPhiLeptDau1Jet= new TH1F("DeltaPhiLeptDau1Jet","#Delta#phi(l-jet) leading (di-jet);#Delta #phi (l-jet) leading [rad]; lepton-jet pair",100,0,TMath::Pi() );
  TH1F * h2 = new TH1F("h2", "h2", 100, 0,TMath::Pi());
  TH1F * h3 = new TH1F("h3", "h3", 100, 0,TMath::Pi());
  
  Events.Project("DeltaPhiLeptDau1Jet","abs(ZeeLeptDau1Phi - JetPhi[0])"," -TMath::Pi() < (ZeeLeptDau1Phi - JetPhi[0]) < TMath::Pi() && floats_JetNtuplizer_JetEta_TEST.@obj.size()==2"); 
  Events.Project("h2", "abs(ZeeLeptDau1Phi - JetPhi[0] - 2 * TMath::Pi())", "(ZeeLeptDau1Phi - JetPhi[0]) > TMath::Pi() && floats_JetNtuplizer_JetEta_TEST.@obj.size()==2" );
  Events.Project("h3", "abs(ZeeLeptDau1Phi - JetPhi[0] + 2 * TMath::Pi())", "(ZeeLeptDau1Phi - JetPhi[0]) <= -TMath::Pi() && floats_JetNtuplizer_JetEta_TEST.@obj.size()==2" );
  
  DeltaPhiLeptDau1Jet->Add(h2,h3);
  DeltaPhiLeptDau1Jet->Draw();
  DeltaPhiLeptDau1Jet->Write();

  // Check on sub-leading lepton-jet delta phi
  
  TH1F * DeltaPhiLeptDau2Jet= new TH1F("DeltaPhiLeptDau2Jet","#Delta#phi(l-jet) sub-leading (di-jet);#Delta #phi (l-jet) sub-leading [rad]; lepton-jet pair",100,0,TMath::Pi() );
  TH1F * h4 = new TH1F("h4", "h4", 100, 0,TMath::Pi());
  TH1F * h5 = new TH1F("h5", "h5", 100, 0,TMath::Pi());
  
  Events.Project("DeltaPhiLeptDau2Jet","abs(ZeeLeptDau2Phi - JetPhi[1])"," -TMath::Pi() < (ZeeLeptDau2Phi - JetPhi[1]) < TMath::Pi() && floats_JetNtuplizer_JetEta_TEST.@obj.size()==2"); 
  Events.Project("h4", "abs(ZeeLeptDau2Phi - JetPhi[1] - 2 * TMath::Pi())", "(ZeeLeptDau2Phi - JetPhi[1]) > TMath::Pi() && floats_JetNtuplizer_JetEta_TEST.@obj.size()==2" );
  Events.Project("h5", "abs(ZeeLeptDau2Phi - JetPhi[1] + 2 * TMath::Pi())", "(ZeeLeptDau2Phi - JetPhi[1]) <= -TMath::Pi() && floats_JetNtuplizer_JetEta_TEST.@obj.size()==2" );
  
  DeltaPhiLeptDau2Jet->Add(h4,h5);
  DeltaPhiLeptDau2Jet->Draw();
  DeltaPhiLeptDau2Jet->Write();

  // Check on leading lepton-jet delta eta
  
  TH1F * DeltaEtaLeptDau1Jet= new TH1F("DeltaEtaLeptDau1Jet","#Delta#eta(l-jet) leading (di-jet);#Delta #eta (l-jet) leading; lepton-jet pair",100,-5,5);
  Events.Project("DeltaEtaLeptDau1Jet","ZeeLeptDau1Eta - JetEta[0]","floats_JetNtuplizer_JetEta_TEST.@obj.size()==2"); 
  
  // Check on sub-leading lepton-jet delta eta
 
  TH1F * DeltaEtaLeptDau2Jet= new TH1F("DeltaEtaLeptDau2Jet","#Delta#eta(l-jet) sub-leading (di-jet);#Delta #eta (l-jet) sub-leading; lepton-jet pair",100,-5,5);
  Events.Project("DeltaEtaLeptDau2Jet","ZeeLeptDau2Eta - JetEta[1]","floats_JetNtuplizer_JetEta_TEST.@obj.size()==2"); 
  
  // Check on leading lepton-jet delta pt
  
  TH1F * DeltaptLeptDau1Jet= new TH1F("DeltaptLeptDau1Jet","#Delta p_{T}(l-jet) leading (di-jet);#Delta p_{T} (l-jet) leading [GeV]; lepton-jet pair",100,0,100);
  Events.Project("DeltaptLeptDau1Jet","abs(ZeeLeptDau1Pt - JetPt[0])","floats_JetNtuplizer_JetPt_TEST.@obj.size()==2"); 
  
 // Check on sub-leading lepton-jet delta pt
  
  TH1F * DeltaptLeptDau2Jet= new TH1F("DeltaptLeptDau2Jet","#Delta p_{T}(l-jet) sub-leading (di-jet);#Delta p_{T} (l-jet) sub-leading [GeV]; lepton-jet pair",100,0,100);
  Events.Project("DeltaptLeptDau2Jet","abs(ZeeLeptDau2Pt - JetPt[1])","floats_JetNtuplizer_JetPt_TEST.@obj.size()==2"); 
  
  // Jet multiplicity (suboptimal)
  TH1F * JetMultiplicity = new TH1F("JetMultiplicity","Jet multiplicity; jet multiplicity; events",15,-0.5,14.5);
  Events.Draw("floats_JetNtuplizer_JetEta_TEST.@obj.size()>>JetMultiplicity",goldenCutOnNumberOjJets);
  
  // Jet b-tagging (SSVHP)
  TH1F * h_JetSSVHP = new TH1F("h_JetSSVHP","SSVHP discriminant; SSVHP discriminant; Jets",100,-5,5);
  Events.Project("h_JetSSVHP","JetJetSSVHP","floats_JetNtuplizer_JetPt_TEST.@obj.size()>1");   
  
  // Jet b-tagging (SSVHE)
  TH1F * h_JetSSVHE = new TH1F("h_JetSSVHE","SSVHE discriminant; SSVHE discriminant; Jets",100,-5,5);
  Events.Project("h_JetSSVHE","JetJetSSVHE","floats_JetNtuplizer_JetPt_TEST.@obj.size()>1"); 

//  __  __ ___ _____                _      _    _        
// |  \/  | __|_   _| __ ____ _ _ _(_)__ _| |__| |___ ___
// | |\/| | _|  | |   \ V / _` | '_| / _` | '_ \ / -_|_-<
// |_|  |_|___| |_|    \_/\__,_|_| |_\__,_|_.__/_\___/__/
                                                      
  TH1F * h_MET = new TH1F("h_MET","MET; MET [GeV]; events/GeV",100,0,100);
  Events.Project("h_MET","pfMetPt","floats_JetNtuplizer_JetPt_TEST.@obj.size()>1"); 
  
//                       _    _       
//  __ _ _ _ __ _ _ __| |_ (_)__ ___
// / _` | '_/ _` | '_ \ ' \| / _(_-<
// \__, |_| \__,_| .__/_||_|_\__/__/
// |___/         |_|                

  //  setGraphics(zmmMass,"mumu");
  //  setGraphics(zmmPt,"mumu"); 
  //  setGraphics(zmmEta,"mumu"); 
  //  setGraphics(zmmPhi,"mumu"); 
  //  setGraphics(zmmY,"mumu"); 
  //  setGraphics(zmmDau1Pt,"mumu"); 
  //  setGraphics(zmmDau2Pt,"mumu"); 
  //  setGraphics(zmmDau1Eta,"mumu");
  //  setGraphics(zmmDau1Phi,"mumu");
  //  setGraphics(zmmDau2Eta,"mumu");
  //  setGraphics(zmmDau2Phi,"mumu");
  
  //  setGraphics(zeeMass,"ee");
  //  setGraphics(zeePt,"ee"); 
  //  setGraphics(zeeEta,"ee"); 
  //  setGraphics(zeePhi,"ee"); 
  //  setGraphics(zeeY,"ee"); 
  //  setGraphics(zeeDau1Pt,"ee"); 
  //  setGraphics(zeeDau2Pt,"ee"); 
  //  setGraphics(zeeDau1Eta,"ee");
  //  setGraphics(zeeDau1Phi,"ee");
  //  setGraphics(zeeDau2Eta,"ee");
  //  setGraphics(zeeDau2Phi,"ee");
  
  //Zmm
  
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

  // Daughters by charge

  zmmDauplusPt->Draw();  
  zmmDauplusPt->Write();  

  zmmDauMinusPt->Draw();  
  zmmDauMinusPt->Write();  

  zmmDauplusEta->Draw(); 
  zmmDauplusEta->Write();

  zmmDauplusPhi->Draw(); 
  zmmDauplusPhi->Write(); 

  zmmDauMinusEta->Draw();
  zmmDauMinusEta->Write();

  zmmDauMinusPhi->Draw(); 
  zmmDauMinusPhi->Write(); 

  zeeDauplusPt->Draw();  
  zeeDauplusPt->Write();  

  zeeDauMinusPt->Draw();  
  zeeDauMinusPt->Write();  

  zeeDauplusEta->Draw(); 
  zeeDauplusEta->Write();

  zeeDauplusPhi->Draw(); 
  zeeDauplusPhi->Write(); 

  zeeDauMinusEta->Draw();
  zeeDauMinusEta->Write();

  zeeDauMinusPhi->Draw(); 
  zeeDauMinusPhi->Write();
  
  // Jets

  LeadingJetPt->Draw(); 
  LeadingJetPt->Write(); 
  
  SubLeadingJetPt->Draw(); 
  SubLeadingJetPt->Write(); 
  
  JetMultiplicity->Draw(); 
  JetMultiplicity->Write(); 
  
  dijetMass->Draw(); 
  dijetMass->Write(); 
  
  diElectronMass->Draw();
  diElectronMass->Write();
  
  // DeltaPhiLeptDau1Jet->Draw();
  // DeltaPhiLeptDau1Jet->Write();
  
  // DeltaPhiLeptDau2Jet->Draw();
  // DeltaPhiLeptDau2Jet->Write();
  
  DeltaEtaLeptDau1Jet->Draw();
  DeltaEtaLeptDau1Jet->Write();
  
  DeltaEtaLeptDau2Jet->Draw();
  DeltaEtaLeptDau2Jet->Write();
  
  DeltaptLeptDau1Jet->Draw();
  DeltaptLeptDau1Jet->Write();
  
  DeltaptLeptDau2Jet->Draw();
  DeltaptLeptDau2Jet->Write();
     
  h_JetSSVHP->Draw();
  h_JetSSVHP->Write();

  h_JetSSVHE->Draw();
  h_JetSSVHE->Write(); 
  
  h_MET->Draw();
  h_MET->Write();

  output_file->Close();

}
