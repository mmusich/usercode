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

void MatchTheTree(const TString& matrix_filename = "MatrixOfMatches.root"){
  
  // input 
  TFile *file_of_matches = new TFile(matrix_filename,"READ");  
  cout << matrix_filename << " open" << endl;

  TGraphErrors* matchList = (TGraphErrors*)file_of_matches->Get("matchList");

  TObjString *tos_file_in1=(TObjString*)file_of_matches->Get("rootfile1");
  TObjString *tos_label1=(TObjString*)file_of_matches->Get("label1");

  TObjString *tos_file_in2=(TObjString*)file_of_matches->Get("rootfile2");
  TObjString *tos_label2=(TObjString*)file_of_matches->Get("label2");

  // keep it for debug
  cout << tos_label1->GetString() << endl;
  cout << tos_label2->GetString() << endl;


  TFile *file_in1 = new TFile(tos_file_in1->GetString(), "READ");
  TTree *tree1 = (TTree*)file_in1->Get("t");

  TFile *file_in2 = new TFile(tos_file_in2->GetString(), "READ");
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

  // Jet infos: general
  Int_t   MCTrueFlavor1[nMaxjets_], MCTrueFlavor2[nMaxjets_];
  Float_t pt1[nMaxjets_],           pt2[nMaxjets_];
  Float_t eta1[nMaxjets_],          eta2[nMaxjets_];
  Float_t phi1[nMaxjets_],          phi2[nMaxjets_];


  // Jet-id
  Float_t  jetNeutralHadronEnergyFraction1[nMaxjets_],jetNeutralHadronEnergyFraction2[nMaxjets_];    // nhf
  Float_t  jetNeutralEmEnergyFraction1[nMaxjets_],    jetNeutralEmEnergyFraction2[nMaxjets_];        // nEF
  Int_t    jetnConstituents1[nMaxjets_],              jetnConstituents2[nMaxjets_];                  // nconstituents
  Float_t  jetChargedMultiplicity1[nMaxjets_],        jetChargedMultiplicity2[nMaxjets_];            // nch
  Float_t  jetChargedHadronEnergyFraction1[nMaxjets_],jetChargedHadronEnergyFraction2[nMaxjets_];    // chf
  Float_t  jetChargedEmEnergyFraction1[nMaxjets_],    jetChargedEmEnergyFraction2[nMaxjets_];        // cef            


  // Jet infos: SV
  // TODO: add some text to describe the variables
  Float_t         SV3dDistance1[nMaxjets_],      SV3dDistance2[nMaxjets_];   
  Float_t         SV3dDistanceError1[nMaxjets_], SV3dDistanceError2[nMaxjets_];
  Float_t         SV2dDistance1[nMaxjets_],      SV2dDistance2[nMaxjets_];   
  Float_t         SV2dDistanceError1[nMaxjets_], SV2dDistanceError2[nMaxjets_];
  Float_t         SVChi21[nMaxjets_],            SVChi22[nMaxjets_];
  Float_t         SVDegreesOfFreedom1[nMaxjets_],SVDegreesOfFreedom2[nMaxjets_];   
  Float_t         SVNormChi21[nMaxjets_],        SVNormChi22[nMaxjets_];  
  Float_t         SVMass1[nMaxjets_],            SVMass2[nMaxjets_];  
  Int_t           SVtotCharge1[nMaxjets_],       SVtotCharge2[nMaxjets_]; 
  Int_t           SVnVertices1[nMaxjets_],       SVnVertices2[nMaxjets_];
  Int_t           SVnVertexTracks1[nMaxjets_],   SVnVertexTracks2[nMaxjets_];
  Int_t           SVnVertexTracksAll1[nMaxjets_],SVnVertexTracksAll2[nMaxjets_];

  // Jet infos: b-tag discriminators
  Float_t discrcsvglobal1[nMaxjets_], discrcsvglobal2[nMaxjets_];
  Float_t discrjpglobal1[nMaxjets_], discrjpglobal2[nMaxjets_];
  Float_t discrjbpglobal1[nMaxjets_], discrjbpglobal2[nMaxjets_];
  Float_t discrssvheglobal1[nMaxjets_], discrssvheglobal2[nMaxjets_];
  Float_t discrssvhpglobal1[nMaxjets_], discrssvhpglobal2[nMaxjets_];
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

  tree1->SetBranchAddress("jetPt",&pt1);
  tree1->SetBranchAddress("jetPhi",&phi1);
  tree1->SetBranchAddress("jetEta",&eta1);
  tree1->SetBranchAddress("nJets",&njet1);
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
  
  tree2->SetBranchAddress("jetPt",&pt2);
  tree2->SetBranchAddress("jetPhi",&phi2);
  tree2->SetBranchAddress("jetEta",&eta2);
  tree2->SetBranchAddress("nJets",&njet2);
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

  tree2->SetBranchAddress("standardSimpleSecondaryVertexHighPurPFBJetTags",&discrssvhpglobal2);
  tree2->SetBranchAddress("standardJetProbabilityPFBJetTags",&discrjpglobal2);
  tree2->SetBranchAddress("standardCombinedSecondaryVertexPFBJetTags",&discrcsvglobal2);
  tree2->SetBranchAddress("standardJetBProbabilityPFBJetTags",&discrjbpglobal2);
  tree2->SetBranchAddress("standardSimpleSecondaryVertexHighEffPFBJetTags",&discrssvheglobal2);
  tree2->SetBranchAddress("standardSimpleSecondaryVertexHighPurPFBJetTags",&discrssvhpglobal2);
  tree2->SetBranchAddress("standardTrackCountingHighEffPFBJetTags",&discrtcheglobal2);
  tree2->SetBranchAddress("standardTrackCountingHighPurPFBJetTags",&discrtchpglobal2);


  // output 
  TFile *file_out=new TFile("JetByJetComparison.root","RECREATE");  
  file_out->cd();

  TH1F *hDeltapt = new TH1F("hDeltapt","#Delta p_{T};#Delta p_{T} [GeV];jets",1000,-50,50);
  TH1F *hDeltaR =  new TH1F("hDeltaR","#Delta R;#sqrt{(#Delta #eta)^{2}  + (#Delta #phi)^{2}};jets",1000,0.,0.01);
  TH1F* hDeltaDiscrTCHE = new TH1F("hDeltaDiscrTCHE","#Delta DiscrTCHE ;#Delta DiscrTCHE;jets",1000,-1,1);
  TH1F* hDeltaDiscrTCHP = new TH1F("hDeltaDiscrTCHP","#Delta DiscrTCHP ;#Delta DiscrTCHP ;jets",1000,-1,1);
  TH1F* hDeltaDiscrSSVHE = new TH1F("hDeltaDiscrSSVHE","#Delta DiscrSSVHE ;#Delta DiscrSSVHE ;jets",1000,-1,1);
  TH1F* hDeltaDiscrSSVHP = new TH1F("hDeltaDiscrSSVHP","#Delta DiscrSSVHP ;#Delta DiscrSSVHP ;jets",1000,-1,1);
  TH1F* hDeltaDiscrCSV = new TH1F("hDeltaDiscrCSV","#Delta DiscrCSV ;#Delta DiscrCSV ;jets",1000,-1,1);
  TH1F* hDeltaDiscrJP = new TH1F("hDeltaDiscrJP","#Delta DiscrJP ;#Delta DiscrJP ;jets",1000,-1,1);
  TH1F* hDeltaDiscrJBP = new TH1F("hDeltaDiscrJBP","#Delta DiscrJBP ;#Delta DiscrJBP ;jets",1000,-1,1);

  // output TTree: event infos
  Int_t   Trun, Tlumi, Tevent;
  Float_t Tpthat, Tmcweight;
  Bool_t  TisBGluonSplitting, TisCGluonSplitting;

  // output TTree: jet infos
  Float_t Tpt[2], Teta[2], Tphi[2];
  Int_t TMCTrueFlavor[2];
  Int_t TjetId[2];
  Float_t Ttche[2], Ttchp[2], Tssvhe[2], Tssvhp[2], Tcsv[2], Tjp[2], Tjbp[2];

  Float_t  TSV3dDistance[2], TSV3dDistanceError[2], TSV2dDistance[2], TSV2dDistanceError[2], TSVChi2[2], TSVDegreesOfFreedom[2], TSVNormChi2[2], TSVMass[2];
  Int_t  TSVtotCharge[2], TSVnVertices[2], TSVnVertexTracks[2], TSVnVertexTracksAll[2];


  TTree *tout = new TTree("JetByJetComparison","file1 vs file2");

  tout->Branch("run",   &Trun,   "Trun/I");
  tout->Branch("lumi",  &Tlumi,  "Tlumi/I");
  tout->Branch("evt",   &Tevent, "Tevent/I");

  tout->Branch("pthat",            &Tpthat,            "Tpthat/F");		  
  tout->Branch("mcweight",         &Tmcweight,	       "Tmcweight/F");	  
  tout->Branch("isBGluonSplitting",&TisBGluonSplitting,"TisBGluonSplitting/B");
  tout->Branch("isCGluonSplitting",&TisCGluonSplitting,"TisCGluonSplitting/B");

  tout->Branch("pt",    &Tpt,    "Tpt[2]/F");
  tout->Branch("eta",   &Teta,   "Teta[2]/F");
  tout->Branch("phi",   &Tphi,   "Tphi[2]/F");
  tout->Branch("jetId", &TjetId, "TjetId[2]/I");
  tout->Branch("MCTrueFlavor",&TMCTrueFlavor,"TMCTrueFlavor[2]/I");

  tout->Branch("SV3dDistance",      &TSV3dDistance,      "TSV3dDistance[2]/F");
  tout->Branch("SV3dDistanceError", &TSV3dDistanceError, "TSV3dDistanceError[2]/F");
  tout->Branch("SV2dDistance",      &TSV2dDistance,      "TSV2dDistance[2]/F");
  tout->Branch("SV2dDistanceError", &TSV2dDistanceError, "TSV2dDistanceError[2]/F");
  tout->Branch("SVChi2",            &TSVChi2,            "TSVChi2[2]/F");
  tout->Branch("SVDegreesOfFreedom",&TSVDegreesOfFreedom,"TSVDegreesOfFreedom[2]/F");
  tout->Branch("SVNormChi2",        &TSVNormChi2,        "TSVNormChi2[2]/F");
  tout->Branch("SVMass",            &TSVMass,            "TSVMass[2]/F");
  tout->Branch("SVtotCharge",       &TSVtotCharge,       "TSVtotCharge[2]/I");
  tout->Branch("SVnVertices",       &TSVnVertices,       "TSVnVertices[2]/I");
  tout->Branch("SVnVertexTracks",   &TSVnVertexTracks,   "TSVnVertexTracks[2]/I");
  tout->Branch("SVnVertexTracksAll",&TSVnVertexTracksAll,"TSVnVertexTracksAll[2]/I");
  
  tout->Branch("tche", &Ttche, "Ttche[2]/F");
  tout->Branch("tchp", &Ttchp, "Ttchp[2]/F");
  tout->Branch("ssvhe",&Tssvhe,"Tssvhe[2]/F");
  tout->Branch("ssvhp",&Tssvhp,"Tssvhp[2]/F");
  tout->Branch("csv",  &Tcsv,  "Tcsv[2]/F"); 
  tout->Branch("jp",   &Tjp,   "Tjp[2]/F"); 
  tout->Branch("jbp",  &Tjbp,  "Tjbp[2]/F");

  const Float_t dRmatch(0.1);

  for (Int_t i=0; i<matchList->GetN(); i++) {
      
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

	// matched jet infos from tree #1
	Tpt[0]    =pt1[j1];
	Teta[0]   =eta1[j1];
	Tphi[0]   =phi1[j1];

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

 	Ttche[0] =discrtcheglobal1[j1];	    
 	Ttchp[0] =discrtchpglobal1[j1];	    
 	Tssvhe[0]=discrssvheglobal1[j1];    
 	Tssvhp[0]=discrssvhpglobal1[j1];    
 	Tcsv[0]  =discrcsvglobal1[j1];	    
 	Tjp[0]   =discrjpglobal1[j1];        
 	Tjbp[0]  =discrjbpglobal1[j1];               

	TMCTrueFlavor[0] = MCTrueFlavor1[j1];

	// matched jet infos from tree #2
	Tpt[1]    =pt2[j2min];
	Teta[1]   =eta2[j2min];
	Tphi[1]   =phi2[j2min];

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

 	Ttche[1] =discrtcheglobal2[j2min];	    
 	Ttchp[1] =discrtchpglobal2[j2min];	    
 	Tssvhe[1]=discrssvheglobal2[j2min];    
 	Tssvhp[1]=discrssvhpglobal2[j2min];    
 	Tcsv[1]  =discrcsvglobal2[j2min];	    
 	Tjp[1]   =discrjpglobal2[j2min];        
 	Tjbp[1]  =discrjbpglobal2[j2min];               

	TMCTrueFlavor[1] = MCTrueFlavor2[j2min];

	// sanity check about MC truth flavor
	if ( MCTrueFlavor1[j1]!=MCTrueFlavor2[j2min] ) 
	  cout << "MC True flavor mismatch: " << MCTrueFlavor1[j1] << ", " << MCTrueFlavor2[j2min] << endl;

	tout->Fill();

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

