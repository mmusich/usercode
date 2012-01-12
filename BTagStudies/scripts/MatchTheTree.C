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

void MatchTheTree(const TString& matrix_filename = "MatrixOfMatches.root"){

  // input 
  TFile *file_of_matches = new TFile(matrix_filename,"open");  
  cout << matrix_filename << " open" << endl;

  TObjString *tos_file_in1=(TObjString*)file_of_matches->Get("rootfile1");
  TObjString *tos_label1=(TObjString*)file_of_matches->Get("label1");

  TObjString *tos_file_in2=(TObjString*)file_of_matches->Get("rootfile2");
  TObjString *tos_label2=(TObjString*)file_of_matches->Get("label2");

  // keep it for debug
  cout << tos_label1->GetString() << endl;
  cout << tos_label2->GetString() << endl;

  TFile *file_in1 = new TFile(tos_file_in1->GetString(), "open");
  TTree *tree1 = (TTree*)file_in1->Get("t");

  TFile *file_in2 = new TFile(tos_file_in2->GetString(), "open");
  TTree *tree2 = (TTree*)file_in2->Get("t");
 
  const int nMaxjets_ = 100000;
  
  Int_t njet1=0; 
  Int_t njet2=0;
  
  UInt_t m_runNumber1, m_lumiNumber1, m_eventNumber1;
  UInt_t m_runNumber2, m_lumiNumber2, m_eventNumber2;

  Float_t discrcsvglobal1[nMaxjets_], discrcsvglobal2[nMaxjets_];
  Float_t discrjpglobal1[nMaxjets_], discrjpglobal2[nMaxjets_];
  Float_t discrjbpglobal1[nMaxjets_], discrjbpglobal2[nMaxjets_];
  Float_t discrssvheglobal1[nMaxjets_], discrssvheglobal2[nMaxjets_];
  Float_t discrssvhpglobal1[nMaxjets_], discrssvhpglobal2[nMaxjets_];
  Float_t discrtcheglobal1[nMaxjets_], discrtcheglobal2[nMaxjets_];
  Float_t discrtchpglobal1[nMaxjets_], discrtchpglobal2[nMaxjets_];
  Float_t pt1[nMaxjets_], pt2[nMaxjets_];
  Float_t eta1[nMaxjets_], eta2[nMaxjets_];
  Float_t phi1[nMaxjets_], phi2[nMaxjets_];
  
  //  Double_t chisquared1[nMaxtracks_], chisquared2[nMaxtracks_];
  //  Double_t ndof1[nMaxtracks_], ndof2[nMaxtracks_]; 
  //  Double_t normchisq1[nMaxtracks_], normchisq2[nMaxtracks_];
  //  tree1->SetBranchAddress("nTracks",&ntrack1);
  //  tree2->SetBranchAddress("nTracks",&ntrack2);
  
  tree1->SetBranchAddress("jetPt",&pt1);
  tree2->SetBranchAddress("jetPt",&pt2);
  tree1->SetBranchAddress("jetPhi",&phi1);
  tree2->SetBranchAddress("jetPhi",&phi2);
  tree1->SetBranchAddress("jetEta",&eta1);
  tree2->SetBranchAddress("jetEta",&eta2);
  tree1->SetBranchAddress("nJets",&njet1);
  tree2->SetBranchAddress("nJets",&njet2);
  
  tree1->SetBranchAddress("runNumber",&m_runNumber1);
  tree1->SetBranchAddress("eventNumber",&m_eventNumber1);
  tree1->SetBranchAddress("lumiBlockNumber",&m_lumiNumber1);
  
  tree2->SetBranchAddress("runNumber",&m_runNumber2);
  tree2->SetBranchAddress("eventNumber",&m_eventNumber2);
  tree2->SetBranchAddress("lumiBlockNumber",&m_lumiNumber2);

  tree1->SetBranchAddress("standardCombinedSecondaryVertexPFBJetTags",&discrcsvglobal1);
  tree2->SetBranchAddress("standardCombinedSecondaryVertexPFBJetTags",&discrcsvglobal2);
  tree1->SetBranchAddress("standardJetProbabilityPFBJetTags",&discrjpglobal1);
  tree2->SetBranchAddress("standardJetProbabilityPFBJetTags",&discrjpglobal2);
  tree1->SetBranchAddress("standardJetBProbabilityPFBJetTags",&discrjbpglobal1);
  tree2->SetBranchAddress("standardJetBProbabilityPFBJetTags",&discrjbpglobal2);
  tree1->SetBranchAddress("standardSimpleSecondaryVertexHighEffPFBJetTags",&discrssvheglobal1);
  tree2->SetBranchAddress("standardSimpleSecondaryVertexHighEffPFBJetTags",&discrssvheglobal2);
  tree1->SetBranchAddress("standardSimpleSecondaryVertexHighPurPFBJetTags",&discrssvhpglobal1);
  tree2->SetBranchAddress("standardSimpleSecondaryVertexHighPurPFBJetTags",&discrssvhpglobal2);
  tree1->SetBranchAddress("standardTrackCountingHighEffPFBJetTags",&discrtcheglobal1);
  tree2->SetBranchAddress("standardTrackCountingHighEffPFBJetTags",&discrtcheglobal2);
  tree1->SetBranchAddress("standardTrackCountingHighPurPFBJetTags",&discrtchpglobal1);
  tree2->SetBranchAddress("standardTrackCountingHighPurPFBJetTags",&discrtchpglobal2);

  TGraphErrors* matchList = (TGraphErrors*)file_of_matches->Get("matchList");

  // output 
  TFile *file_out=new TFile("JetByJetComparison.root","recreate");  
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

  Int_t Trun, Tlumi, Tevent;
  Float_t Tpt, Teta, Tphi;
  Float_t Ttche1, Ttchp1, Tssvhe1, Tssvhp1, Tcsv1, Tjp1, Tjbp1;
  Float_t Ttche2, Ttchp2, Tssvhe2, Tssvhp2, Tcsv2, Tjp2, Tjbp2;

  TTree *tout = new TTree("JetByJetComparison","file1 vs file2");

  tout->Branch("run",   &Trun,   "Trun/I");
  tout->Branch("lumi",  &Tlumi,  "Tlumi/I");
  tout->Branch("evt",   &Tevent, "Tevent/I");
  tout->Branch("pt",    &Tpt,    "Tpt/F");
  tout->Branch("eta",   &Teta,   "Teta/F");
  tout->Branch("phi",   &Tphi,   "Tphi/F");

  tout->Branch("tche1", &Ttche1, "Ttche1/F");
  tout->Branch("tchp1", &Ttchp1, "Ttchp1/F");
  tout->Branch("ssvhe1",&Tssvhe1,"Tssvhe1/F");
  tout->Branch("ssvhp1",&Tssvhp1,"Tssvhp1/F");
  tout->Branch("csv1",  &Tcsv1,  "Tcsv1/F"); 
  tout->Branch("jp1",   &Tjp1,   "Tjp1/F"); 
  tout->Branch("jbp1",  &Tjbp1,  "Tjbp1/F");

  tout->Branch("tche2", &Ttche2, "Ttche2/F");
  tout->Branch("tchp2", &Ttchp2, "Ttchp2/F");
  tout->Branch("ssvhe2",&Tssvhe2,"Tssvhe2/F");
  tout->Branch("ssvhp2",&Tssvhp2,"Tssvhp2/F");
  tout->Branch("csv2",  &Tcsv2,  "Tcsv2/F"); 
  tout->Branch("jp2",   &Tjp2,   "Tjp2/F"); 
  tout->Branch("jbp2",  &Tjbp2,  "Tjbp2/F");

  const Float_t dRmatch(0.1);
  for (Int_t i=0; i<matchList->GetN(); i++) {
      
    Int_t entry1 = (Int_t)matchList->GetX()[i];
    Int_t entry2 = (Int_t)matchList->GetY()[i];
    
    tree1->GetEntry(entry1);
    tree2->GetEntry(entry2);
    
    // cout <<"ntracks1:"<<ntrack1<<" ntracks2 "<<ntrack2<< endl;
    // cout<<"---------------------------------------"<<endl;
    // cout<<"iRun1:"<<m_runNumber1<<" iRun2:"<<m_runNumber2<<" | iLumi1:"<<m_lumiNumber1<<" iLumi2:"<<m_lumiNumber2 <<" | iEvent1:"<<m_eventNumber1 <<" iEvent2:"<<m_eventNumber2<<endl;
    // cout <<"entry1: "<<entry1<<" entry2: "<<entry2<< endl;
    // cout <<"pt1:    "<<pt1[0]<<" pt2:    "<<pt2[0]<< endl;   
    
   
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

      if (dRmin<dRmatch){	// jet is matched in dR
	hDeltapt->Fill(pt1[j1]-pt2[j2min]);
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
	Tpt    =pt1[j1];
	Teta   =eta1[j1];
	Tphi   =phi1[j1];
	Ttche1 =discrtcheglobal1[j1];	    
	Ttchp1 =discrtchpglobal1[j1];	    
	Tssvhe1=discrssvheglobal1[j1];    
	Tssvhp1=discrssvhpglobal1[j1];    
	Tcsv1  =discrcsvglobal1[j1];	    
	Tjp1   =discrjpglobal1[j1];        
	Tjbp1  =discrjbpglobal1[j1];               
	Ttche2 =discrtcheglobal2[j2min];	    
	Ttchp2 =discrtchpglobal2[j2min];	    
	Tssvhe2=discrssvheglobal2[j2min];    
	Tssvhp2=discrssvhpglobal2[j2min];    
	Tcsv2  =discrcsvglobal2[j2min];	    
	Tjp2   =discrjpglobal2[j2min];        
	Tjbp2  =discrjbpglobal2[j2min];               
	tout->Fill();
      }

    }// end of loop on the jets in tree #1    

  }

  file_in1->Close();
  file_in2->Close();
  file_of_matches->Close();

  file_out->cd();
  file_out->ls();
  file_out->Write();
  file_out->Close();

}

