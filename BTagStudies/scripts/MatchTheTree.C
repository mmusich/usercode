#include <cmath>
#include <iostream>
#include <vector>
#include "TBranch.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLeafI.h"
#include "TTree.h"
#include "TMath.h"
#include "TString.h"
#include "TH1F.h"
#include "TObjString.h"
#include "JetInfo.h" 

#include "TSystem.h"
#include "TROOT.h"

using namespace std;

Int_t GetJetID(Double_t nhf, Double_t nEF, Double_t nconstituents, Double_t chf, Double_t nch, Double_t cef){  

  /*
    Double_t nhf = ( jet.neutralHadronEnergy() + jet.HFHadronEnergy() ) / jet.energy();  [NB: in the ntuple HFHadronEnergy is absent!]
    Double_t nEF = jet.neutralEmEnergyFraction();
    Double_t nconstituents = jet.numberOfDaughters();
    Double_t chf = jet.chargedHadronEnergyFraction();
    Double_t nch = jet.chargedMultiplicity();
    Double_t cef = jet.chargedEmEnergyFraction();
  */
  
  Int_t goodjetID;
  if(nhf<0.90 && nEF<0.90 && nconstituents>1 && chf>0 && nch>0 && cef<0.99) {
    goodjetID=3;  //tight
  } else if(nhf<0.95 && nEF<0.95 && nconstituents>1 && chf>0 && nch>0 && cef<0.99) {
    goodjetID=2;   //medium
  } else if(nhf<0.99 && nEF<0.99 && nconstituents>1 && chf>0 && nch>0 && cef<0.99) {
    goodjetID=1;   //loose
  } else {
    goodjetID=0.;
  }
  return goodjetID; 
}

void MatchTheTree(bool doTree, const TString& matrix_filename = "MatrixOfMatches.root",Int_t maxEvents=-1){
  // input 
  TFile *file_of_matches = new TFile(matrix_filename,"READ");  
  cout << matrix_filename << " open" << endl;

  TGraphErrors* matchList = (TGraphErrors*)file_of_matches->Get("matchList");

  TObjString *tos_file_in1=(TObjString*)file_of_matches->Get("rootfile1");
  TObjString *tos_label1=(TObjString*)file_of_matches->Get("label1");

  TObjString *tos_file_in2=(TObjString*)file_of_matches->Get("rootfile2");
  TObjString *tos_label2=(TObjString*)file_of_matches->Get("label2");

  // keep it for debug
  cout <<"first input file:"<< tos_label1->GetString() << endl;
  cout <<"second input file:"<< tos_label2->GetString() << endl;

  TFile *file_in1 = TFile::Open(tos_file_in1->GetString());
  TTree *tree1 = (TTree*)file_in1->Get("t");

  TFile *file_in2 = TFile::Open(tos_file_in2->GetString());
  TTree *tree2 = (TTree*)file_in2->Get("t");
 
  const int nMaxjets_ = 50;
  
  Int_t njet1=0; 
  Int_t njet2=0;

  // Event infos  
  UInt_t m_runNumber1, m_lumiNumber1, m_eventNumber1;
  UInt_t m_runNumber2, m_lumiNumber2, m_eventNumber2;

  // Event infos [MC]
  Float_t         m_pthat1, m_mcweight1;
  Bool_t          m_isBGluonSplitting1, m_isCGluonSplitting1;

  Float_t         m_pthat2, m_mcweight2;
  Bool_t          m_isBGluonSplitting2, m_isCGluonSplitting2;

  // Event infos [vertexing]
  Float_t         m_PVx1, m_PVx2, m_PVy1, m_PVy2, m_PVz1, m_PVz2;
  Float_t         m_PVChiSq1, m_PVChiSq2, m_PVndof1, m_PVndof2, m_PVNormChiSq1, m_PVNormChiSq2; 

  // Jet infos: general
  Int_t   MCTrueFlavor1[nMaxjets_],   MCTrueFlavor2[nMaxjets_];
  Int_t   nTracks1[nMaxjets_],        nTracks2[nMaxjets_];
  Float_t pt1[nMaxjets_],             pt2[nMaxjets_];
  Float_t eta1[nMaxjets_],            eta2[nMaxjets_];
  Float_t phi1[nMaxjets_],            phi2[nMaxjets_];

  // Infos on the most significant track
  Float_t IP3dFirst1[nMaxjets_],                   IP3dFirst2[nMaxjets_];
  Float_t IP3dErrorFirst1[nMaxjets_],              IP3dErrorFirst2[nMaxjets_];
  Float_t IP3dDecayLengthFirst1[nMaxjets_],        IP3dDecayLengthFirst2[nMaxjets_];
  Float_t IP3dTransverseMomentumFirst1[nMaxjets_], IP3dTransverseMomentumFirst2[nMaxjets_];
  Float_t IP3dEtaFirst1[nMaxjets_],                IP3dEtaFirst2[nMaxjets_];
  Float_t IP3dPhiFirst1[nMaxjets_],                IP3dPhiFirst2[nMaxjets_];
  
  // Infos on the most significant track
  Float_t IP3dSecond1[nMaxjets_],                  IP3dSecond2[nMaxjets_];
  Float_t IP3dErrorSecond1[nMaxjets_],             IP3dErrorSecond2[nMaxjets_];
  Float_t IP3dDecayLengthSecond1[nMaxjets_],       IP3dDecayLengthSecond2[nMaxjets_];
  Float_t IP3dTransverseMomentumSecond1[nMaxjets_],IP3dTransverseMomentumSecond2[nMaxjets_];
  Float_t IP3dEtaSecond1[nMaxjets_],               IP3dEtaSecond2[nMaxjets_];
  Float_t IP3dPhiSecond1[nMaxjets_],               IP3dPhiSecond2[nMaxjets_];

  // Infos on the most significant track
  Float_t IP3dThird1[nMaxjets_],                   IP3dThird2[nMaxjets_];
  Float_t IP3dErrorThird1[nMaxjets_],              IP3dErrorThird2[nMaxjets_];
  Float_t IP3dDecayLengthThird1[nMaxjets_],        IP3dDecayLengthThird2[nMaxjets_];
  Float_t IP3dTransverseMomentumThird1[nMaxjets_], IP3dTransverseMomentumThird2[nMaxjets_];
  Float_t IP3dEtaThird1[nMaxjets_],                IP3dEtaThird2[nMaxjets_];
  Float_t IP3dPhiThird1[nMaxjets_],                IP3dPhiThird2[nMaxjets_];

  // Infos on the most significant track
  Float_t IP3dFourth1[nMaxjets_],                  IP3dFourth2[nMaxjets_];
  Float_t IP3dErrorFourth1[nMaxjets_],             IP3dErrorFourth2[nMaxjets_];
  Float_t IP3dDecayLengthFourth1[nMaxjets_],       IP3dDecayLengthFourth2[nMaxjets_];
  Float_t IP3dTransverseMomentumFourth1[nMaxjets_],IP3dTransverseMomentumFourth2[nMaxjets_];
  Float_t IP3dEtaFourth1[nMaxjets_],               IP3dEtaFourth2[nMaxjets_];
  Float_t IP3dPhiFourth1[nMaxjets_],               IP3dPhiFourth2[nMaxjets_];

  // Jet-id
  Float_t jetNeutralHadronEnergyFraction1[nMaxjets_],jetNeutralHadronEnergyFraction2[nMaxjets_];    // nhf
  Float_t jetNeutralEmEnergyFraction1[nMaxjets_],    jetNeutralEmEnergyFraction2[nMaxjets_];        // nEF
  Int_t   jetnConstituents1[nMaxjets_],              jetnConstituents2[nMaxjets_];                  // nconstituents
  Float_t jetChargedMultiplicity1[nMaxjets_],        jetChargedMultiplicity2[nMaxjets_];            // nch
  Float_t jetChargedHadronEnergyFraction1[nMaxjets_],jetChargedHadronEnergyFraction2[nMaxjets_];    // chf
  Float_t jetChargedEmEnergyFraction1[nMaxjets_],    jetChargedEmEnergyFraction2[nMaxjets_];        // cef            

  // Jet infos: SV
  // TODO: add some text to describe the variables
  Float_t SV3dDistance1[nMaxjets_],      SV3dDistance2[nMaxjets_];   
  Float_t SV3dDistanceError1[nMaxjets_], SV3dDistanceError2[nMaxjets_];
  Float_t SV2dDistance1[nMaxjets_],      SV2dDistance2[nMaxjets_];   
  Float_t SV2dDistanceError1[nMaxjets_], SV2dDistanceError2[nMaxjets_];
  Float_t SVChi21[nMaxjets_],            SVChi22[nMaxjets_];
  Float_t SVDegreesOfFreedom1[nMaxjets_],SVDegreesOfFreedom2[nMaxjets_];   
  Float_t SVNormChi21[nMaxjets_],        SVNormChi22[nMaxjets_];  
  Float_t SVMass1[nMaxjets_],            SVMass2[nMaxjets_];  
  Int_t   SVtotCharge1[nMaxjets_],       SVtotCharge2[nMaxjets_]; 
  Int_t   SVnVertices1[nMaxjets_],       SVnVertices2[nMaxjets_];
  Int_t   SVnVertexTracks1[nMaxjets_],   SVnVertexTracks2[nMaxjets_];
  Int_t   SVnVertexTracksAll1[nMaxjets_],SVnVertexTracksAll2[nMaxjets_];

  // Jet infos: b-tag discriminators
  Float_t discrcsvglobal1[nMaxjets_],  discrcsvglobal2[nMaxjets_];
  Float_t discrjpglobal1[nMaxjets_],   discrjpglobal2[nMaxjets_];
  Float_t discrjbpglobal1[nMaxjets_],  discrjbpglobal2[nMaxjets_];
  Float_t discrssvheglobal1[nMaxjets_],discrssvheglobal2[nMaxjets_];
  Float_t discrssvhpglobal1[nMaxjets_],discrssvhpglobal2[nMaxjets_];
  Float_t discrtcheglobal1[nMaxjets_], discrtcheglobal2[nMaxjets_];
  Float_t discrtchpglobal1[nMaxjets_], discrtchpglobal2[nMaxjets_];
 
  // TTree #1
  tree1->SetBranchAddress("runNumber",&m_runNumber1);
  tree1->SetBranchAddress("eventNumber",&m_eventNumber1);
  tree1->SetBranchAddress("lumiBlockNumber",&m_lumiNumber1);

  tree1->SetBranchAddress("pthat",&m_pthat1);
  tree1->SetBranchAddress("mcweight",&m_mcweight1);
  tree1->SetBranchAddress("isBGluonSplitting",&m_isBGluonSplitting1);
  tree1->SetBranchAddress("isCGluonSplitting",&m_isCGluonSplitting1);  
  
  tree1->SetBranchAddress("PVx",&m_PVx1);
  tree1->SetBranchAddress("PVy",&m_PVy1);
  tree1->SetBranchAddress("PVz",&m_PVz1);
  tree1->SetBranchAddress("PVChi2",&m_PVChiSq1);
  tree1->SetBranchAddress("PVndof",&m_PVndof1);
  tree1->SetBranchAddress("PVNormalizedChi2",&m_PVNormChiSq1);

  tree1->SetBranchAddress("jetPt",&pt1);
  tree1->SetBranchAddress("jetPhi",&phi1);
  tree1->SetBranchAddress("jetEta",&eta1);
  tree1->SetBranchAddress("nJets",&njet1);
  tree1->SetBranchAddress("jetnTracks",&nTracks1);
  tree1->SetBranchAddress("MCTrueFlavor",&MCTrueFlavor1);

  tree1->SetBranchAddress("jetNeutralHadronEnergyFraction",&jetNeutralHadronEnergyFraction1);    
  tree1->SetBranchAddress("jetNeutralEmEnergyFraction",    &jetNeutralEmEnergyFraction1);        
  tree1->SetBranchAddress("jetnConstituents",              &jetnConstituents1);                  
  tree1->SetBranchAddress("jetChargedMultiplicity",        &jetChargedMultiplicity1);            
  tree1->SetBranchAddress("jetChargedHadronEnergyFraction",&jetChargedHadronEnergyFraction1);    
  tree1->SetBranchAddress("jetChargedEmEnergyFraction",    &jetChargedEmEnergyFraction1);        

  tree1->SetBranchAddress("SV3dDistance",&SV3dDistance1);
  tree1->SetBranchAddress("SV3dDistanceError",&SV3dDistanceError1);
  tree1->SetBranchAddress("SV2dDistance",&SV2dDistance1);
  tree1->SetBranchAddress("SV2dDistanceError",&SV2dDistanceError1);
  tree1->SetBranchAddress("SVChi2",&SVChi21);
  tree1->SetBranchAddress("SVDegreesOfFreedom",&SVDegreesOfFreedom1);
  tree1->SetBranchAddress("SVNormChi2",&SVNormChi21);
  tree1->SetBranchAddress("SVMass",&SVMass1);
  tree1->SetBranchAddress("SVtotCharge",&SVtotCharge1);
  tree1->SetBranchAddress("SVnVertices",&SVnVertices1);
  tree1->SetBranchAddress("SVnVertexTracks",&SVnVertexTracks1);
  tree1->SetBranchAddress("SVnVertexTracksAll",&SVnVertexTracksAll1);

  tree1->SetBranchAddress("IP3d1",&IP3dFirst1);		     
  tree1->SetBranchAddress("IP3dError1",&IP3dErrorFirst1);		     
  tree1->SetBranchAddress("IP3dDecayLength1",&IP3dDecayLengthFirst1);	     
  tree1->SetBranchAddress("IP3dTransverseMomentum1",&IP3dTransverseMomentumFirst1); 
  tree1->SetBranchAddress("IP3dEta1",&IP3dEtaFirst1);		     
  tree1->SetBranchAddress("IP3dPhi1",&IP3dPhiFirst1);		     
			   			  		    							     
  tree1->SetBranchAddress("IP3d2",&IP3dSecond1);		     
  tree1->SetBranchAddress("IP3dError2",&IP3dErrorSecond1);	     
  tree1->SetBranchAddress("IP3dDecayLength2",&IP3dDecayLengthSecond1);	     
  tree1->SetBranchAddress("IP3dTransverseMomentum2",&IP3dTransverseMomentumSecond1);
  tree1->SetBranchAddress("IP3dEta2",&IP3dEtaSecond1);               
  tree1->SetBranchAddress("IP3dPhi2",&IP3dPhiSecond1);               
			   			  	   	       				     
  tree1->SetBranchAddress("IP3d3",&IP3dThird1);		     		                     
  tree1->SetBranchAddress("IP3dError3",&IP3dErrorThird1);              
  tree1->SetBranchAddress("IP3dDecayLength3",&IP3dDecayLengthThird1);        
  tree1->SetBranchAddress("IP3dTransverseMomentum3",&IP3dTransverseMomentumThird1); 
  tree1->SetBranchAddress("IP3dEta3",&IP3dEtaThird1);                
  tree1->SetBranchAddress("IP3dPhi3",&IP3dPhiThird1);                
			   	   			   			  		     				     
  tree1->SetBranchAddress("IP3d4",&IP3dFourth1);                  
  tree1->SetBranchAddress("IP3dError4",&IP3dErrorFourth1);             
  tree1->SetBranchAddress("IP3dDecayLength4",&IP3dDecayLengthFourth1);       
  tree1->SetBranchAddress("IP3dTransverseMomentum4",&IP3dTransverseMomentumFourth1);
  tree1->SetBranchAddress("IP3dEta4",&IP3dEtaFourth1);               
  tree1->SetBranchAddress("IP3dPhi4",&IP3dPhiFourth1);               

  tree1->SetBranchAddress("standardCombinedSecondaryVertexPFBJetTags",&discrcsvglobal1);
  tree1->SetBranchAddress("standardJetProbabilityPFBJetTags",&discrjpglobal1);
  tree1->SetBranchAddress("standardJetBProbabilityPFBJetTags",&discrjbpglobal1);
  tree1->SetBranchAddress("standardSimpleSecondaryVertexHighEffPFBJetTags",&discrssvheglobal1);
  tree1->SetBranchAddress("standardSimpleSecondaryVertexHighPurPFBJetTags",&discrssvhpglobal1);
  tree1->SetBranchAddress("standardTrackCountingHighEffPFBJetTags",&discrtcheglobal1);
  tree1->SetBranchAddress("standardTrackCountingHighPurPFBJetTags",&discrtchpglobal1);

  // TTree #2
  tree2->SetBranchAddress("runNumber",&m_runNumber2);
  tree2->SetBranchAddress("eventNumber",&m_eventNumber2);
  tree2->SetBranchAddress("lumiBlockNumber",&m_lumiNumber2);

  tree2->SetBranchAddress("pthat",&m_pthat2);
  tree2->SetBranchAddress("mcweight",&m_mcweight2);
  tree2->SetBranchAddress("isBGluonSplitting",&m_isBGluonSplitting2);
  tree2->SetBranchAddress("isCGluonSplitting",&m_isCGluonSplitting2);
  
  tree2->SetBranchAddress("PVx",&m_PVx2);
  tree2->SetBranchAddress("PVy",&m_PVy2);
  tree2->SetBranchAddress("PVz",&m_PVz2);
  tree2->SetBranchAddress("PVChi2",&m_PVChiSq2);
  tree2->SetBranchAddress("PVndof",&m_PVndof2);
  tree2->SetBranchAddress("PVNormalizedChi2",&m_PVNormChiSq2);

  tree2->SetBranchAddress("jetPt",&pt2);
  tree2->SetBranchAddress("jetPhi",&phi2);
  tree2->SetBranchAddress("jetEta",&eta2);
  tree2->SetBranchAddress("nJets",&njet2);
  tree1->SetBranchAddress("jetnTracks",&nTracks2);
  tree2->SetBranchAddress("MCTrueFlavor",&MCTrueFlavor2);

  tree2->SetBranchAddress("jetNeutralHadronEnergyFraction",&jetNeutralHadronEnergyFraction2);    
  tree2->SetBranchAddress("jetNeutralEmEnergyFraction",    &jetNeutralEmEnergyFraction2);        
  tree2->SetBranchAddress("jetnConstituents",              &jetnConstituents2);                  
  tree2->SetBranchAddress("jetChargedMultiplicity",        &jetChargedMultiplicity2);            
  tree2->SetBranchAddress("jetChargedHadronEnergyFraction",&jetChargedHadronEnergyFraction2);    
  tree2->SetBranchAddress("jetChargedEmEnergyFraction",    &jetChargedEmEnergyFraction2);        

  tree2->SetBranchAddress("SV3dDistance",&SV3dDistance2);
  tree2->SetBranchAddress("SV3dDistanceError",&SV3dDistanceError2);
  tree2->SetBranchAddress("SV2dDistance",&SV2dDistance2);
  tree2->SetBranchAddress("SV2dDistanceError",&SV2dDistanceError2);
  tree2->SetBranchAddress("SVChi2",&SVChi22);
  tree2->SetBranchAddress("SVDegreesOfFreedom",&SVDegreesOfFreedom2);
  tree2->SetBranchAddress("SVNormChi2",&SVNormChi22);
  tree2->SetBranchAddress("SVMass",&SVMass2);
  tree2->SetBranchAddress("SVtotCharge",&SVtotCharge2);
  tree2->SetBranchAddress("SVnVertices",&SVnVertices2);
  tree2->SetBranchAddress("SVnVertexTracks",&SVnVertexTracks2);
  tree2->SetBranchAddress("SVnVertexTracksAll",&SVnVertexTracksAll2);

  tree2->SetBranchAddress("IP3d1",&IP3dFirst2);		     
  tree2->SetBranchAddress("IP3dError1",&IP3dErrorFirst2);		     
  tree2->SetBranchAddress("IP3dDecayLength1",&IP3dDecayLengthFirst2);	     
  tree2->SetBranchAddress("IP3dTransverseMomentum1",&IP3dTransverseMomentumFirst2); 
  tree2->SetBranchAddress("IP3dEta1",&IP3dEtaFirst2);		     
  tree2->SetBranchAddress("IP3dPhi1",&IP3dPhiFirst2);		     
			   			  		    							     
  tree2->SetBranchAddress("IP3d2",&IP3dSecond2);		     
  tree2->SetBranchAddress("IP3dError2",&IP3dErrorSecond2);	     
  tree2->SetBranchAddress("IP3dDecayLength2",&IP3dDecayLengthSecond2);	     
  tree2->SetBranchAddress("IP3dTransverseMomentum2",&IP3dTransverseMomentumSecond2);
  tree2->SetBranchAddress("IP3dEta2",&IP3dEtaSecond2);               
  tree2->SetBranchAddress("IP3dPhi2",&IP3dPhiSecond2);               
			   			  	   	       				     
  tree2->SetBranchAddress("IP3d3",&IP3dThird2);		     		                     
  tree2->SetBranchAddress("IP3dError3",&IP3dErrorThird2);              
  tree2->SetBranchAddress("IP3dDecayLength3",&IP3dDecayLengthThird2);        
  tree2->SetBranchAddress("IP3dTransverseMomentum3",&IP3dTransverseMomentumThird2); 
  tree2->SetBranchAddress("IP3dEta3",&IP3dEtaThird2);                
  tree2->SetBranchAddress("IP3dPhi3",&IP3dPhiThird2);                
			   	   			   			  		     				     
  tree2->SetBranchAddress("IP3d4",&IP3dFourth2);                  
  tree2->SetBranchAddress("IP3dError4",&IP3dErrorFourth2);             
  tree2->SetBranchAddress("IP3dDecayLength4",&IP3dDecayLengthFourth2);       
  tree2->SetBranchAddress("IP3dTransverseMomentum4",&IP3dTransverseMomentumFourth2);
  tree2->SetBranchAddress("IP3dEta4",&IP3dEtaFourth2);               
  tree2->SetBranchAddress("IP3dPhi4",&IP3dPhiFourth2);

  tree2->SetBranchAddress("standardSimpleSecondaryVertexHighPurPFBJetTags",&discrssvhpglobal2);
  tree2->SetBranchAddress("standardJetProbabilityPFBJetTags",&discrjpglobal2);
  tree2->SetBranchAddress("standardCombinedSecondaryVertexPFBJetTags",&discrcsvglobal2);
  tree2->SetBranchAddress("standardJetBProbabilityPFBJetTags",&discrjbpglobal2);
  tree2->SetBranchAddress("standardSimpleSecondaryVertexHighEffPFBJetTags",&discrssvheglobal2);
  tree2->SetBranchAddress("standardSimpleSecondaryVertexHighPurPFBJetTags",&discrssvhpglobal2);
  tree2->SetBranchAddress("standardTrackCountingHighEffPFBJetTags",&discrtcheglobal2);
  tree2->SetBranchAddress("standardTrackCountingHighPurPFBJetTags",&discrtchpglobal2);


  // output 
  TString tupleType_;
  if (doTree) tupleType_= "Tree";
  else tupleType_= "Ntuple";
		
  TString append = "JetByJetComparison"+tupleType_+"_"+tos_label1->GetString()+"Vs"+tos_label2->GetString()+".root";

  TFile *file_out=new TFile(append,"RECREATE");  
  file_out->cd();

  TH1F *hDeltapt         = new TH1F("hDeltapt","#Delta p_{T};#Delta p_{T} [GeV];jets",1000,-50,50);
  TH1F *hDeltaR          = new TH1F("hDeltaR","#Delta R;#sqrt{(#Delta #eta)^{2}  + (#Delta #phi)^{2}};jets",1000,0.,0.01);
  TH1F* hDeltaDiscrTCHE  = new TH1F("hDeltaDiscrTCHE","#Delta DiscrTCHE ;#Delta DiscrTCHE;jets",1000,-1,1);
  TH1F* hDeltaDiscrTCHP  = new TH1F("hDeltaDiscrTCHP","#Delta DiscrTCHP ;#Delta DiscrTCHP ;jets",1000,-1,1);
  TH1F* hDeltaDiscrSSVHE = new TH1F("hDeltaDiscrSSVHE","#Delta DiscrSSVHE ;#Delta DiscrSSVHE ;jets",1000,-1,1);
  TH1F* hDeltaDiscrSSVHP = new TH1F("hDeltaDiscrSSVHP","#Delta DiscrSSVHP ;#Delta DiscrSSVHP ;jets",1000,-1,1);
  TH1F* hDeltaDiscrCSV   = new TH1F("hDeltaDiscrCSV","#Delta DiscrCSV ;#Delta DiscrCSV ;jets",1000,-1,1);
  TH1F* hDeltaDiscrJP    = new TH1F("hDeltaDiscrJP","#Delta DiscrJP ;#Delta DiscrJP ;jets",1000,-1,1);
  TH1F* hDeltaDiscrJBP   = new TH1F("hDeltaDiscrJBP","#Delta DiscrJBP ;#Delta DiscrJBP ;jets",1000,-1,1);

  // output TTree: event infos
  Int_t   Trun, Tlumi, Tevent;
  Float_t Tpthat, Tmcweight;
  Bool_t  TisBGluonSplitting, TisCGluonSplitting;

  // TTree
  TTree *tree_out = 0;
  JetInfo *TJetInfoA=0;
  JetInfo *TJetInfoB=0;

  // ntuple
  TTree *ntuple_out = 0;
  Float_t TPVx[2],TPVy[2],TPVz[2];
  Float_t TPVChiSq[2], TPVndof[2], TPVNormChiSq[2]; 
  Int_t   TnTracks[2];
  
  Float_t Tpt[2], Teta[2], Tphi[2];
  Int_t   TMCTrueFlavor[2];    
  Int_t   TjetId[2];
  Float_t Ttche[2], Ttchp[2], Tssvhe[2], Tssvhp[2], Tcsv[2], Tjp[2], Tjbp[2];
  
  Float_t TSV3dDistance[2], TSV3dDistanceError[2], TSV2dDistance[2], TSV2dDistanceError[2], TSVChi2[2], TSVDegreesOfFreedom[2], TSVNormChi2[2], TSVMass[2];    
  Int_t   TSVtotCharge[2], TSVnVertices[2], TSVnVertexTracks[2], TSVnVertexTracksAll[2];
  
  Float_t TIP3dFirst[2],TIP3dErrorFirst[2],TIP3dDecayLengthFirst[2],TIP3dTransverseMomentumFirst[2],TIP3dEtaFirst[2],TIP3dPhiFirst[2];                	   
  Float_t TIP3dSecond[2],TIP3dErrorSecond[2],TIP3dDecayLengthSecond[2],TIP3dTransverseMomentumSecond[2],TIP3dEtaSecond[2],TIP3dPhiSecond[2];               
  Float_t TIP3dThird[2],TIP3dErrorThird[2],TIP3dDecayLengthThird[2],TIP3dTransverseMomentumThird[2],TIP3dEtaThird[2],TIP3dPhiThird[2];                	
  Float_t TIP3dFourth[2],TIP3dErrorFourth[2],TIP3dDecayLengthFourth[2],TIP3dTransverseMomentumFourth[2],TIP3dEtaFourth[2],TIP3dPhiFourth[2];                   


  if ( doTree ) {
    tree_out = new TTree("JetByJetComparisonTree","tree file1 vs file2");
    tree_out->Branch("run",   &Trun,   "Trun/I");
    tree_out->Branch("lumi",  &Tlumi,  "Tlumi/I");
    tree_out->Branch("evt",   &Tevent, "Tevent/I"); 
    tree_out->Branch("pthat",            &Tpthat,      "Tpthat/F");		  
    tree_out->Branch("mcweight",         &Tmcweight,   "Tmcweight/F");	  
    tree_out->Branch("isBGluonSplitting",&TisBGluonSplitting,"TisBGluonSplitting/B");
    tree_out->Branch("isCGluonSplitting",&TisCGluonSplitting,"TisCGluonSplitting/B");
    
    // output TTree: jet infos 
    TJetInfoA = new JetInfo();
    TJetInfoB = new JetInfo();
    tree_out->Branch("JetInfoA","JetInfo",&TJetInfoA);
    tree_out->Branch("JetInfoB","JetInfo",&TJetInfoB);
  } else {

    // output ntuple
    ntuple_out = new TTree("JetByJetComparisonNtuple","ntuple file1 vs file2");
    ntuple_out->Branch("run",   &Trun,   "Trun/I");
    ntuple_out->Branch("lumi",  &Tlumi,  "Tlumi/I");
    ntuple_out->Branch("evt",   &Tevent, "Tevent/I"); 
    ntuple_out->Branch("pthat",            &Tpthat,      "Tpthat/F");		  
    ntuple_out->Branch("mcweight",         &Tmcweight,   "Tmcweight/F");	  
    ntuple_out->Branch("isBGluonSplitting",&TisBGluonSplitting,"TisBGluonSplitting/B");
    ntuple_out->Branch("isCGluonSplitting",&TisCGluonSplitting,"TisCGluonSplitting/B");
    ntuple_out->Branch("isCGluonSplitting",&TisCGluonSplitting,"TisCGluonSplitting/B");

    ntuple_out->Branch("pt",    &Tpt,    "Tpt[2]/F");
    ntuple_out->Branch("eta",   &Teta,   "Teta[2]/F");
    ntuple_out->Branch("phi",   &Tphi,   "Tphi[2]/F");
    ntuple_out->Branch("jetId", &TjetId, "TjetId[2]/I");
    ntuple_out->Branch("jetnTracks",&TnTracks, "TnTracks[2]/I");
    ntuple_out->Branch("MCTrueFlavor",&TMCTrueFlavor,"TMCTrueFlavor[2]/I");

    ntuple_out->Branch("SV3dDistance",      &TSV3dDistance,      "TSV3dDistance[2]/F");
    ntuple_out->Branch("SV3dDistanceError", &TSV3dDistanceError, "TSV3dDistanceError[2]/F");
    ntuple_out->Branch("SV2dDistance",      &TSV2dDistance,      "TSV2dDistance[2]/F");
    ntuple_out->Branch("SV2dDistanceError", &TSV2dDistanceError, "TSV2dDistanceError[2]/F");
    ntuple_out->Branch("SVChi2",            &TSVChi2,            "TSVChi2[2]/F");
    ntuple_out->Branch("SVDegreesOfFreedom",&TSVDegreesOfFreedom,"TSVDegreesOfFreedom[2]/F");
    ntuple_out->Branch("SVNormChi2",        &TSVNormChi2,        "TSVNormChi2[2]/F");
    ntuple_out->Branch("SVMass",            &TSVMass,            "TSVMass[2]/F");
    ntuple_out->Branch("SVtotCharge",       &TSVtotCharge,       "TSVtotCharge[2]/I");
    ntuple_out->Branch("SVnVertices",       &TSVnVertices,       "TSVnVertices[2]/I");
    ntuple_out->Branch("SVnVertexTracks",   &TSVnVertexTracks,   "TSVnVertexTracks[2]/I");
    ntuple_out->Branch("SVnVertexTracksAll",&TSVnVertexTracksAll,"TSVnVertexTracksAll[2]/I");

    ntuple_out->Branch("IP3d1",                  &TIP3dFirst		       ,"TIP3dFirst[2]/F");		     	      
    ntuple_out->Branch("IP3dError1",             &TIP3dErrorFirst	       ,"TIP3dErrorFirst[2]/F");		   	     
    ntuple_out->Branch("IP3dDecayLength1",       &TIP3dDecayLengthFirst	       ,"TIP3dDecayLengthFirst[2]/F");	     	   
    ntuple_out->Branch("IP3dTransverseMomentum1",&TIP3dTransverseMomentumFirst ,"TIP3dTransverseMomentumFirst[2]/F"); 	   
    ntuple_out->Branch("IP3dEta1",               &TIP3dEtaFirst		       ,"TIP3dEtaFirst[2]/F");		     	   
    ntuple_out->Branch("IP3dPhi1",               &TIP3dPhiFirst		       ,"TIP3dPhiFirst[2]/F");		     	   
    
    ntuple_out->Branch("IP3d2",                  &TIP3dSecond		       ,"TIP3dSecond[2]/F");		     	   
    ntuple_out->Branch("IP3dError2",             &TIP3dErrorSecond	       ,"TIP3dErrorSecond[2]/F");	     	   
    ntuple_out->Branch("IP3dDecayLength2",       &TIP3dDecayLengthSecond       ,"TIP3dDecayLengthSecond[2]/F");	   	     
    ntuple_out->Branch("IP3dTransverseMomentum2",&TIP3dTransverseMomentumSecond,"TIP3dTransverseMomentumSecond[2]/F");	   
    ntuple_out->Branch("IP3dEta2",               &TIP3dEtaSecond               ,"TIP3dEtaSecond[2]/F");               	   
    ntuple_out->Branch("IP3dPhi2",               &TIP3dPhiSecond               ,"TIP3dPhiSecond[2]/F");               	   

    ntuple_out->Branch("IP3d3",                  &TIP3dThird		       ,"TIP3dThird[2]/F");		     	   		                     
    ntuple_out->Branch("IP3dError3",             &TIP3dErrorThird              ,"TIP3dErrorThird[2]/F");              	   
    ntuple_out->Branch("IP3dDecayLength3",       &TIP3dDecayLengthThird        ,"TIP3dDecayLengthThird[2]/F");            
    ntuple_out->Branch("IP3dTransverseMomentum3",&TIP3dTransverseMomentumThird ,"TIP3dTransverseMomentumThird[2]/F"); 	   
    ntuple_out->Branch("IP3dEta3",               &TIP3dEtaThird                ,"TIP3dEtaThird[2]/F");                    
    ntuple_out->Branch("IP3dPhi3",               &TIP3dPhiThird                ,"TIP3dPhiThird[2]/F");                	   
    
    ntuple_out->Branch("IP3d4",                  &TIP3dFourth                  ,"TIP3dFourth[2]/F");                  	   
    ntuple_out->Branch("IP3dError4",             &TIP3dErrorFourth             ,"TIP3dErrorFourth[2]/F");                 
    ntuple_out->Branch("IP3dDecayLength4",       &TIP3dDecayLengthFourth       ,"TIP3dDecayLengthFourth[2]/F");       	   
    ntuple_out->Branch("IP3dTransverseMomentum4",&TIP3dTransverseMomentumFourth,"TIP3dTransverseMomentumFourth[2]/F");
    ntuple_out->Branch("IP3dEta4",               &TIP3dEtaFourth               ,"TIP3dEtaFourth[2]/F");               	   
    ntuple_out->Branch("IP3dPhi4",               &TIP3dPhiFourth               ,"TIP3dPhiFourth[2]/F");                        

    ntuple_out->Branch("tche", &Ttche, "Ttche[2]/F");
    ntuple_out->Branch("tchp", &Ttchp, "Ttchp[2]/F");
    ntuple_out->Branch("ssvhe",&Tssvhe,"Tssvhe[2]/F");
    ntuple_out->Branch("ssvhp",&Tssvhp,"Tssvhp[2]/F");
    ntuple_out->Branch("csv",  &Tcsv,  "Tcsv[2]/F"); 
    ntuple_out->Branch("jp",   &Tjp,   "Tjp[2]/F"); 
    ntuple_out->Branch("jbp",  &Tjbp,  "Tjbp[2]/F");
    
  }

  const Float_t dRmatch(0.1);
  
  Int_t nMaxEvents_(matchList->GetN());
  if(maxEvents!=-1) nMaxEvents_ = maxEvents; 
  
  for (Int_t i=0; i<nMaxEvents_; i++) {
      
    Int_t entry1 = (Int_t)matchList->GetX()[i];
    Int_t entry2 = (Int_t)matchList->GetY()[i];
    
    tree1->GetEntry(entry1);
    tree2->GetEntry(entry2);

    // DEBUG
    // cout<<"---------------------------------------"<<endl;
    // cout<<"iRun1:"<<m_runNumber1<<" iRun2:"<<m_runNumber2<<" | iLumi1:"<<m_lumiNumber1<<" iLumi2:"<<m_lumiNumber2 <<" | iEvent1:"<<m_eventNumber1 <<" iEvent2:"<<m_eventNumber2<<endl;
    // cout <<"entry1: "<<entry1<<" entry2: "<<entry2<< endl;
    
    if ( njet1 > nMaxjets_ || njet2 > nMaxjets_ ) {
      cout << "more than " << nMaxjets_ << " jets " << max(njet1,njet2) << endl;
      continue;
    }

    for(Int_t j1 = 0; j1<njet1; j1++){ // loop on the jets in tree #1
      // match the jets in dR as the order of the jets in the two files can be different 
      Float_t dRmin=FLT_MAX;
      Int_t j2min(-1);
      for(Int_t j2 = 0; j2<njet2; j2++){ // loop on the jets in tree #2
	Float_t deta=eta1[j1]-eta2[j2];
	Float_t dphi=phi1[j1]-phi2[j2];
	Float_t dRtmp=TMath::Sqrt(deta*deta+dphi*dphi);
	if (dRtmp<dRmin) {
	  dRmin=dRtmp;
	  j2min=j2;
	}
      } // end of loop on the jets in tree #2
      hDeltaR->Fill(dRmin); 

      if (j2min > -1) hDeltapt->Fill(pt1[j1]-pt2[j2min]);
      if (dRmin<dRmatch){	// jet is matched in dR
 	hDeltaDiscrTCHE->Fill(discrtcheglobal1[j1]-discrtcheglobal2[j2min]);
 	hDeltaDiscrTCHP->Fill(discrtchpglobal1[j1]-discrtchpglobal2[j2min]);
 	hDeltaDiscrSSVHE->Fill(discrssvheglobal1[j1]-discrssvheglobal2[j2min]);
 	hDeltaDiscrSSVHP->Fill(discrssvhpglobal1[j1]-discrssvhpglobal2[j2min]); 
 	hDeltaDiscrCSV->Fill(discrcsvglobal1[j1]-discrcsvglobal2[j2min]);
 	hDeltaDiscrJP->Fill(discrjpglobal1[j1]-discrjpglobal2[j2min]);
 	hDeltaDiscrJBP->Fill(discrjbpglobal1[j1]-discrjbpglobal2[j2min]);

	Trun   =m_runNumber1;
	Tlumi  =m_lumiNumber1;
	Tevent =m_eventNumber1;

	Tpthat            =m_pthat1;		  
	Tmcweight	  =m_mcweight1;	  
	TisBGluonSplitting=m_isBGluonSplitting1;
	TisCGluonSplitting=m_isCGluonSplitting1;

	// sanity check about MC truth flavor
	if ( MCTrueFlavor1[j1]!=MCTrueFlavor2[j2min] ) 
	  cout << "MC True flavor mismatch: " << MCTrueFlavor1[j1] << ", " << MCTrueFlavor2[j2min] << endl;

	// tree 
	if ( doTree ) {
	  TJetInfoA->pv.PVx = m_PVx1;		  
	  TJetInfoA->pv.PVy = m_PVy1;		  
	  TJetInfoA->pv.PVz = m_PVz1;		  
	  TJetInfoA->pv.PVChiSq = m_PVChiSq1;	  
	  TJetInfoA->pv.PVndof = m_PVndof1;		  
	  TJetInfoA->pv.PVNormChiSq = m_PVNormChiSq1; 
	  
	  TJetInfoB->pv.PVx = m_PVx2;		  
	  TJetInfoB->pv.PVy = m_PVy2;		  
	  TJetInfoB->pv.PVz = m_PVz2;		  
	  TJetInfoB->pv.PVChiSq = m_PVChiSq2;	  
	  TJetInfoB->pv.PVndof = m_PVndof2;		  
	  TJetInfoB->pv.PVNormChiSq = m_PVNormChiSq2; 
	  

	  // matched jet infos from tree #1
	  TJetInfoA->pt    =pt1[j1];
	  TJetInfoA->eta   =eta1[j1];
	  TJetInfoA->phi   =phi1[j1];
	  TJetInfoA->nTracks = nTracks1[j1];
	
	  TJetInfoA->jetId = GetJetID(
				     jetNeutralHadronEnergyFraction1[j1], // nhf	       
				     jetNeutralEmEnergyFraction1[j1],     // nEF	       
				     jetnConstituents1[j1],               // nconstituents  
				     jetChargedMultiplicity1[j1],         // nch	       
				     jetChargedHadronEnergyFraction1[j1], // chf	       
				     jetChargedEmEnergyFraction1[j1]	  // cef            
				     );   

	  TJetInfoA->sv.SV3dDistance      =SV3dDistance1[j1];
	  TJetInfoA->sv.SV3dDistanceError =SV3dDistanceError1[j1];
	  TJetInfoA->sv.SV2dDistance      =SV2dDistance1[j1];      
	  TJetInfoA->sv.SV2dDistanceError =SV2dDistanceError1[j1];
	  TJetInfoA->sv.SVChi2            =SVChi21[j1];            
	  TJetInfoA->sv.SVDegreesOfFreedom=SVDegreesOfFreedom1[j1];
	  TJetInfoA->sv.SVNormChi2        =SVNormChi21[j1];
	  TJetInfoA->sv.SVMass            =SVMass1[j1];            
	  TJetInfoA->sv.SVtotCharge       =SVtotCharge1[j1];
	  TJetInfoA->sv.SVnVertices       =SVnVertices1[j1];
	  TJetInfoA->sv.SVnVertexTracks   =SVnVertexTracks1[j1];
	  TJetInfoA->sv.SVnVertexTracksAll=SVnVertexTracksAll1[j1];
	  
 	  TJetInfoA->trk[0].IP3d		      =IP3dFirst1[j1];                  
 	  TJetInfoA->trk[0].IP3dError	      =IP3dErrorFirst1[j1];             
 	  TJetInfoA->trk[0].IP3dDecayLength      =IP3dDecayLengthFirst1[j1];       	
 	  TJetInfoA->trk[0].IP3dransverseMomentum=IP3dTransverseMomentumFirst1[j1];	
 	  TJetInfoA->trk[0].IP3dEta	      =IP3dEtaFirst1[j1];               	
 	  TJetInfoA->trk[0].IP3dPhi              =IP3dPhiFirst1[j1];               	
	  
 	  TJetInfoA->trk[1].IP3d		      =IP3dSecond1[j1];                 
 	  TJetInfoA->trk[1].IP3dError	      =IP3dErrorSecond1[j1];            
 	  TJetInfoA->trk[1].IP3dDecayLength      =IP3dDecayLengthSecond1[j1];      	
 	  TJetInfoA->trk[1].IP3dransverseMomentum=IP3dTransverseMomentumSecond1[j1];	
 	  TJetInfoA->trk[1].IP3dEta	      =IP3dEtaSecond1[j1];              	
 	  TJetInfoA->trk[1].IP3dPhi              =IP3dPhiSecond1[j1];              	
	
 	  TJetInfoA->trk[2].IP3d		      =IP3dThird1[j1];                  
 	  TJetInfoA->trk[2].IP3dError	      =IP3dErrorThird1[j1];             	
 	  TJetInfoA->trk[2].IP3dDecayLength      =IP3dDecayLengthThird1[j1];       
 	  TJetInfoA->trk[2].IP3dransverseMomentum=IP3dTransverseMomentumThird1[j1];	
 	  TJetInfoA->trk[2].IP3dEta	      =IP3dEtaThird1[j1];               
 	  TJetInfoA->trk[2].IP3dPhi              =IP3dPhiThird1[j1];               	
	  
 	  TJetInfoA->trk[3].IP3d		      =IP3dFourth1[j1];                 	
 	  TJetInfoA->trk[3].IP3dError	      =IP3dErrorFourth1[j1];            
 	  TJetInfoA->trk[3].IP3dDecayLength      =IP3dDecayLengthFourth1[j1];      	
 	  TJetInfoA->trk[3].IP3dransverseMomentum=IP3dTransverseMomentumFourth1[j1];	
 	  TJetInfoA->trk[3].IP3dEta	      =IP3dEtaFourth1[j1];              	
 	  TJetInfoA->trk[3].IP3dPhi              =IP3dPhiFourth1[j1];              

	  TJetInfoA->tche =discrtcheglobal1[j1];	    
	  TJetInfoA->tchp =discrtchpglobal1[j1];	    
	  TJetInfoA->ssvhe=discrssvheglobal1[j1];    
	  TJetInfoA->ssvhp=discrssvhpglobal1[j1];    
	  TJetInfoA->csv  =discrcsvglobal1[j1];	    
	  TJetInfoA->jp   =discrjpglobal1[j1];        
	  TJetInfoA->jbp  =discrjbpglobal1[j1];               
	  
	  TJetInfoA->MCTrueFlavor = MCTrueFlavor1[j1];

	// -------------------------------------------------------------------------------------------------------------------------------
	// matched jet infos from tree #2
	  TJetInfoB->pt    =pt2[j2min];
	  TJetInfoB->eta   =eta2[j2min];
	  TJetInfoB->phi   =phi2[j2min];
	  TJetInfoB->nTracks = nTracks2[j2min];
	  
	  TJetInfoB->jetId = GetJetID(
				     jetNeutralHadronEnergyFraction2[j2min], // nhf	       
				     jetNeutralEmEnergyFraction2[j2min],     // nEF	       
				     jetnConstituents2[j2min],               // nconstituents  
				     jetChargedMultiplicity2[j2min],         // nch	       
				     jetChargedHadronEnergyFraction2[j2min], // chf	       
				     jetChargedEmEnergyFraction2[j2min]	     // cef            
				     );   
	
	  TJetInfoB->sv.SV3dDistance      =SV3dDistance2[j2min];
	  TJetInfoB->sv.SV3dDistanceError =SV3dDistanceError2[j2min];
	  TJetInfoB->sv.SV2dDistance      =SV2dDistance2[j2min];      
	  TJetInfoB->sv.SV2dDistanceError =SV2dDistanceError2[j2min];
	  TJetInfoB->sv.SVChi2            =SVChi22[j2min];            
	  TJetInfoB->sv.SVDegreesOfFreedom=SVDegreesOfFreedom2[j2min];
	  TJetInfoB->sv.SVNormChi2        =SVNormChi22[j2min];
	  TJetInfoB->sv.SVMass            =SVMass2[j2min];            
	  TJetInfoB->sv.SVtotCharge       =SVtotCharge2[j2min];
	  TJetInfoB->sv.SVnVertices       =SVnVertices2[j2min];
	  TJetInfoB->sv.SVnVertexTracks   =SVnVertexTracks2[j2min];
	  TJetInfoB->sv.SVnVertexTracksAll=SVnVertexTracksAll2[j2min];
	
 	  TJetInfoB->trk[0].IP3d		      =IP3dFirst2[j2min];                  
 	  TJetInfoB->trk[0].IP3dError	      =IP3dErrorFirst2[j2min];             
 	  TJetInfoB->trk[0].IP3dDecayLength      =IP3dDecayLengthFirst2[j2min];       	
 	  TJetInfoB->trk[0].IP3dransverseMomentum=IP3dTransverseMomentumFirst2[j2min];	
 	  TJetInfoB->trk[0].IP3dEta	      =IP3dEtaFirst2[j2min];               	
 	  TJetInfoB->trk[0].IP3dPhi              =IP3dPhiFirst2[j2min];               	
	  
 	  TJetInfoB->trk[1].IP3d		      =IP3dSecond2[j2min];                 
 	  TJetInfoB->trk[1].IP3dError	      =IP3dErrorSecond2[j2min];            
 	  TJetInfoB->trk[1].IP3dDecayLength      =IP3dDecayLengthSecond2[j2min];      	
 	  TJetInfoB->trk[1].IP3dransverseMomentum=IP3dTransverseMomentumSecond2[j2min];	
 	  TJetInfoB->trk[1].IP3dEta	      =IP3dEtaSecond2[j2min];              	
 	  TJetInfoB->trk[1].IP3dPhi              =IP3dPhiSecond2[j2min];              	
	  
 	  TJetInfoB->trk[2].IP3d		      =IP3dThird2[j2min];                  
 	  TJetInfoB->trk[2].IP3dError	      =IP3dErrorThird2[j2min];             	
 	  TJetInfoB->trk[2].IP3dDecayLength      =IP3dDecayLengthThird2[j2min];       
 	  TJetInfoB->trk[2].IP3dransverseMomentum=IP3dTransverseMomentumThird2[j2min];	
 	  TJetInfoB->trk[2].IP3dEta	      =IP3dEtaThird2[j2min];               
 	  TJetInfoB->trk[2].IP3dPhi              =IP3dPhiThird2[j2min];               	
	
 	  TJetInfoB->trk[3].IP3d		      =IP3dFourth2[j2min];                 	
 	  TJetInfoB->trk[3].IP3dError	      =IP3dErrorFourth2[j2min];            
 	  TJetInfoB->trk[3].IP3dDecayLength      =IP3dDecayLengthFourth2[j2min];      	
 	  TJetInfoB->trk[3].IP3dransverseMomentum=IP3dTransverseMomentumFourth2[j2min];	
 	  TJetInfoB->trk[3].IP3dEta	      =IP3dEtaFourth2[j2min];              	
 	  TJetInfoB->trk[3].IP3dPhi              =IP3dPhiFourth2[j2min];              

	  TJetInfoB->tche =discrtcheglobal2[j2min];	    
	  TJetInfoB->tchp =discrtchpglobal2[j2min];	    
	  TJetInfoB->ssvhe=discrssvheglobal2[j2min];    
	  TJetInfoB->ssvhp=discrssvhpglobal2[j2min];    
	  TJetInfoB->csv  =discrcsvglobal2[j2min];	    
	  TJetInfoB->jp   =discrjpglobal2[j2min];        
	  TJetInfoB->jbp  =discrjbpglobal2[j2min];               
	  
	  TJetInfoB->MCTrueFlavor = MCTrueFlavor2[j2min];

	  tree_out->Fill();

	} else {
	  // Flat ntuple
        TPVx[0] = m_PVx1;		  
	TPVy[0] = m_PVy1;		  
	TPVz[0] = m_PVz1;		  
	TPVChiSq[0] = m_PVChiSq1;	  
	TPVndof[0] = m_PVndof1;		  
	TPVNormChiSq[0] = m_PVNormChiSq1; 
       
	TPVx[1] = m_PVx2;		  
	TPVy[1] = m_PVy2;		  
	TPVz[1] = m_PVz2;		  
	TPVChiSq[1] = m_PVChiSq2;	  
	TPVndof[1] = m_PVndof2;		  
	TPVNormChiSq[1] = m_PVNormChiSq2; 

	// matched jet infos from tree #1
	Tpt[0]    =pt1[j1];
	Teta[0]   =eta1[j1];
	Tphi[0]   =phi1[j1];
	TnTracks[0] = nTracks1[j1];
       
	TjetId[0] = GetJetID(
			     jetNeutralHadronEnergyFraction1[j1], // nhf	       
			     jetNeutralEmEnergyFraction1[j1],     // nEF	       
			     jetnConstituents1[j1],               // nconstituents  
			     jetChargedMultiplicity1[j1],         // nch	       
			     jetChargedHadronEnergyFraction1[j1], // chf	       
			     jetChargedEmEnergyFraction1[j1]	    // cef            
			     );   

 	TSV3dDistance[0]      =SV3dDistance1[j1];
 	TSV3dDistanceError[0] =SV3dDistanceError1[j1];
 	TSV2dDistance[0]      =SV2dDistance1[j1];      
 	TSV2dDistanceError[0] =SV2dDistanceError1[j1];
 	TSVChi2[0]            =SVChi21[j1];            
 	TSVDegreesOfFreedom[0]=SVDegreesOfFreedom1[j1];
 	TSVNormChi2[0]        =SVNormChi21[j1];
 	TSVMass[0]            =SVMass1[j1];            
 	TSVtotCharge[0]       =SVtotCharge1[j1];
 	TSVnVertices[0]       =SVnVertices1[j1];
 	TSVnVertexTracks[0]   =SVnVertexTracks1[j1];
 	TSVnVertexTracksAll[0]=SVnVertexTracksAll1[j1];

	TIP3dFirst[0]		     	     =IP3dFirst1[j1];                  
	TIP3dErrorFirst[0]		     =IP3dErrorFirst1[j1];             
	TIP3dDecayLengthFirst[0]	     =IP3dDecayLengthFirst1[j1];       	
	TIP3dTransverseMomentumFirst[0]      =IP3dTransverseMomentumFirst1[j1];	
	TIP3dEtaFirst[0]		     =IP3dEtaFirst1[j1];               	
	TIP3dPhiFirst[0]		     =IP3dPhiFirst1[j1];               		 		    	   	     						

	TIP3dSecond[0]		     	     =IP3dSecond1[j1];                 
	TIP3dErrorSecond[0]	     	     =IP3dErrorSecond1[j1];            
	TIP3dDecayLengthSecond[0]	     =IP3dDecayLengthSecond1[j1];      	
	TIP3dTransverseMomentumSecond[0]     =IP3dTransverseMomentumSecond1[j1];	
	TIP3dEtaSecond[0]                    =IP3dEtaSecond1[j1];              	
	TIP3dPhiSecond[0]                    =IP3dPhiSecond1[j1];              	

	TIP3dThird[0]		     	     =IP3dThird1[j1];                  
	TIP3dErrorThird[0]                   =IP3dErrorThird1[j1];             	
	TIP3dDecayLengthThird[0]             =IP3dDecayLengthThird1[j1];       
	TIP3dTransverseMomentumThird[0]      =IP3dTransverseMomentumThird1[j1];	
	TIP3dEtaThird[0]                     =IP3dEtaThird1[j1];               
	TIP3dPhiThird[0]                     =IP3dPhiThird1[j1];               			   		   	     						

	TIP3dFourth[0]                       =IP3dFourth1[j1];                 	
	TIP3dErrorFourth[0]                  =IP3dErrorFourth1[j1];            
	TIP3dDecayLengthFourth[0]            =IP3dDecayLengthFourth1[j1];      	
	TIP3dTransverseMomentumFourth[0]     =IP3dTransverseMomentumFourth1[j1];	
	TIP3dEtaFourth[0]                    =IP3dEtaFourth1[j1];              	
	TIP3dPhiFourth[0]                    =IP3dPhiFourth1[j1];              

 	Ttche[0] =discrtcheglobal1[j1];	    
 	Ttchp[0] =discrtchpglobal1[j1];	    
 	Tssvhe[0]=discrssvheglobal1[j1];    
 	Tssvhp[0]=discrssvhpglobal1[j1];    
 	Tcsv[0]  =discrcsvglobal1[j1];	    
 	Tjp[0]   =discrjpglobal1[j1];        
 	Tjbp[0]  =discrjbpglobal1[j1];    

	TMCTrueFlavor[0] = MCTrueFlavor1[j1];

	// -------------------------------------------------------------------------------------------------------------------------------
	// matched jet infos from tree #2
	Tpt[1]    =pt2[j2min];
	Teta[1]   =eta2[j2min];
	Tphi[1]   =phi2[j2min];
	TnTracks[1] = nTracks2[j2min];

	TjetId[1] = GetJetID(
			     jetNeutralHadronEnergyFraction2[j2min], // nhf	       
			     jetNeutralEmEnergyFraction2[j2min],     // nEF	       
			     jetnConstituents2[j2min],               // nconstituents  
			     jetChargedMultiplicity2[j2min],         // nch	       
			     jetChargedHadronEnergyFraction2[j2min], // chf	       
			     jetChargedEmEnergyFraction2[j2min]	     // cef            
			     );   
	
 	TSV3dDistance[1]      =SV3dDistance2[j2min];
 	TSV3dDistanceError[1] =SV3dDistanceError2[j2min];
 	TSV2dDistance[1]      =SV2dDistance2[j2min];      
 	TSV2dDistanceError[1] =SV2dDistanceError2[j2min];
 	TSVChi2[1]            =SVChi22[j2min];            
 	TSVDegreesOfFreedom[1]=SVDegreesOfFreedom2[j2min];
 	TSVNormChi2[1]        =SVNormChi22[j2min];
 	TSVMass[1]            =SVMass2[j2min];            
 	TSVtotCharge[1]       =SVtotCharge2[j2min];
 	TSVnVertices[1]       =SVnVertices2[j2min];
 	TSVnVertexTracks[1]   =SVnVertexTracks2[j2min];
 	TSVnVertexTracksAll[1]=SVnVertexTracksAll2[j2min];	

	TIP3dFirst[1]		     	     =IP3dFirst2[j2min];                  
	TIP3dErrorFirst[1]		     =IP3dErrorFirst2[j2min];             
	TIP3dDecayLengthFirst[1]	     =IP3dDecayLengthFirst2[j2min];       	
	TIP3dTransverseMomentumFirst[1]      =IP3dTransverseMomentumFirst2[j2min];	
	TIP3dEtaFirst[1]		     =IP3dEtaFirst2[j2min];               	
	TIP3dPhiFirst[1]		     =IP3dPhiFirst2[j2min];               	

        TIP3dSecond[1]		     	     =IP3dSecond2[j2min];                 
	TIP3dErrorSecond[1]	     	     =IP3dErrorSecond2[j2min];            
	TIP3dDecayLengthSecond[1]	     =IP3dDecayLengthSecond2[j2min];      	
	TIP3dTransverseMomentumSecond[1]     =IP3dTransverseMomentumSecond2[j2min];	
	TIP3dEtaSecond[1]                    =IP3dEtaSecond2[j2min];              	
	TIP3dPhiSecond[1]                    =IP3dPhiSecond2[j2min];              	

	TIP3dThird[1]		     	     =IP3dThird2[j2min];                  
	TIP3dErrorThird[1]                   =IP3dErrorThird2[j2min];             	
	TIP3dDecayLengthThird[1]             =IP3dDecayLengthThird2[j2min];       
	TIP3dTransverseMomentumThird[1]      =IP3dTransverseMomentumThird2[j2min];	
	TIP3dEtaThird[1]                     =IP3dEtaThird2[j2min];               
	TIP3dPhiThird[1]                     =IP3dPhiThird2[j2min];               	

	TIP3dFourth[1]                       =IP3dFourth2[j2min];                 	
	TIP3dErrorFourth[1]                  =IP3dErrorFourth2[j2min];            
	TIP3dDecayLengthFourth[1]            =IP3dDecayLengthFourth2[j2min];      	
	TIP3dTransverseMomentumFourth[1]     =IP3dTransverseMomentumFourth2[j2min];	
	TIP3dEtaFourth[1]                    =IP3dEtaFourth2[j2min];              	
	TIP3dPhiFourth[1]                    =IP3dPhiFourth2[j2min];              

 	Ttche[1] =discrtcheglobal2[j2min];	    
 	Ttchp[1] =discrtchpglobal2[j2min];	    
 	Tssvhe[1]=discrssvheglobal2[j2min];    
 	Tssvhp[1]=discrssvhpglobal2[j2min];    
 	Tcsv[1]  =discrcsvglobal2[j2min];	    
 	Tjp[1]   =discrjpglobal2[j2min];        
 	Tjbp[1]  =discrjbpglobal2[j2min];

	TMCTrueFlavor[1] = MCTrueFlavor2[j2min];

	ntuple_out->Fill();

	}


      } // end of if jet matched clause
    } // end of loop on the jets in tree #1

  } // end of the loop on the TGraphErrors

  file_of_matches->Close();
  file_in1->Close();
  file_in2->Close();

  file_out->cd();
  file_out->ls();
  file_out->Write();
  file_out->Close();

}





	  





