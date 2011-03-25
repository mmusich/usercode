/*  PreselAnalyzer
 *  Analyzer of the preselection in PAT Muons
 */

#include "PreselAnalyzer.h"

// Collaborating Class Header
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

// for "luminosity"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

using namespace std;
using namespace edm;

/// Constructor
PreselAnalyzer::PreselAnalyzer(const ParameterSet& pset) : 
  theMuonLabel_(pset.getUntrackedParameter<string>("MuonCollection")),
  theZllLabel_(pset.getUntrackedParameter<string>("ZllCollection"))
{
}

/// Destructor
PreselAnalyzer::~PreselAnalyzer(){
}

void PreselAnalyzer::beginJob(){

  // Book histograms
  edm::Service<TFileService> fs;

  //Creating sub directories
  TFileDirectory spectra = fs->mkdir("Spectra");
  TFileDirectory subspectra1 = spectra.mkdir("Comparing Presel Mu with other mu");
  TFileDirectory subspectra2 = spectra.mkdir("Studying second mu");
  TFileDirectory subspectra3 = spectra.mkdir("Track Quality");
  TFileDirectory scatter = fs->mkdir("Scatter plots");
  TFileDirectory counters = fs->mkdir("Counters");

  //Spectra histograms
  h1dr_ = subspectra1.make<TH1F>("dr","DeltaR between patMuons",50,0,10);
  h1dphi_ = subspectra1.make<TH1F>("dphi","DeltaPhi between patMuons",50,0,3.15);
  h1deta_ = subspectra1.make<TH1F>("deta","DeltaEta between patMuons",50,0,10);

  //Scatter plots
  h2drpt_ = scatter.make<TH2F>("drpt","DeltaR vs Pt",50,0,10,50,0,100);
  h2dphideta_ = scatter.make<TH2F>("dphideta","Delta Eta vs Delta Phi",50,0,3.15,50,0,5);
  h2dphipt_ = scatter.make<TH2F>("dphipt","DeltaPhi vs Pt",250,0,120,50,0,3.15);
  h2detapt_ = scatter.make<TH2F>("detapt","DeltaEta vs Pt",250,0,120,100,0,10);
  h2drmu_ = scatter.make<TH2F>("drmu","DeltaR vs Muons type esclusive",3,0,3,100,0,10);

  // not selected
  h1Pt_NotSel_ = subspectra2.make<TH1F>("Pt_NotSel","Pt Not selected",50,0,100);
  h1Chi2norm_NotSel_ = subspectra2.make<TH1F>("Chi2norm_NotSel","Chi2norm Not selected",50,0,30);
  h1NValTkHits_NotSel_ = subspectra2.make<TH1F>("NValTkHits_NotSel","NValTkHits Not selected",20,-0.5,19.5);
  h1NValPxHits_NotSel_ = subspectra2.make<TH1F>("NValPxHits_NotSel","NValPxHits Not selected",10,-0.5, 9.5);
  h1NValMuHits_NotSel_ = subspectra2.make<TH1F>("NValMuHits_NotSel","NValMuHits Not selected",20,-0.5,19.5);
  h1DB_NotSel_ = subspectra2.make<TH1F>("DB_NotSel","dB Not selected",20,0,2);
  h1TotIsoOverpT_NotSel_ = subspectra2.make<TH1F>("TotIsoOverlap_NotSel","TotIsoOverlap Not selected",20,0,1.0);
  h1NofMatch_NotSel_ = subspectra2.make<TH1F>("NofMatch_NotSel","NofMatch Not selected",10,-0.5,9.5);
  h1Eta_NotSel_ = subspectra2.make<TH1F>("Eta_NotSel","Eta Not selected",20,0,4);

  // Track quality cut
  h1Pt_TkQuality_ = subspectra3.make<TH1F>("Pt_TkQuality","Pt (TkQuality cut)",50,0,100);
  h1DB_TkQuality_ = subspectra3.make<TH1F>("DB_TkQuality","dB (TkQuality cut)",20,0,2);
  h1Eta_TkQuality_ = subspectra3.make<TH1F>("Eta_TkQuality","Eta (TkQuality cut)",20,0,4);
  h1Phi_TkQuality_ = subspectra3.make<TH1F>("Phi_TkQuality","Phi (TkQuality cut)",50,0,3.15);
  h1Dxy_TkQuality_ = subspectra3.make<TH1F>("Dxy_TkQuality","Dxy (TkQuality cut)",400,-200,200);
  h1Dz_TkQuality_ = subspectra3.make<TH1F>("Dz_TkQuality","Dz (TkQuality cut)",400,-200,200);

  //Presel failure histogram
  h1PreselFail_ = counters.make<TH1F>("PreselFail","Reason for preselection cut failure",11,-0.5,10.5);
  h1PreselFail_->GetXaxis()->SetBinLabel(1,"Pt");
  h1PreselFail_->GetXaxis()->SetBinLabel(2,"non Global");
  h1PreselFail_->GetXaxis()->SetBinLabel(3,"non Tracker");
  h1PreselFail_->GetXaxis()->SetBinLabel(4,"Chi2");
  h1PreselFail_->GetXaxis()->SetBinLabel(5,"Tracker hits");
  h1PreselFail_->GetXaxis()->SetBinLabel(6,"Pixel hits");
  h1PreselFail_->GetXaxis()->SetBinLabel(7,"Muon hits");
  h1PreselFail_->GetXaxis()->SetBinLabel(8,"dB");
  h1PreselFail_->GetXaxis()->SetBinLabel(9,"Isolation");
  h1PreselFail_->GetXaxis()->SetBinLabel(10,"Matches");
  h1PreselFail_->GetXaxis()->SetBinLabel(11,"Eta");
  h1Cut_ = counters.make<TH1F>("Cut","Various counters",3,-0.5,2.5);
  h1Cut_->GetXaxis()->SetBinLabel(1,"All Muons");
  h1Cut_->GetXaxis()->SetBinLabel(2,"Presel Muons");
  h1Cut_->GetXaxis()->SetBinLabel(3,"TkQuality Muons");
}

void PreselAnalyzer::endJob(){

}

void PreselAnalyzer::analyze(const Event & iEvent, const EventSetup& eventSetup)
{
  // Get the Muon collection
  Handle<pat::MuonCollection> muons;
  iEvent.getByLabel(theMuonLabel_, muons);
  pat::MuonCollection::const_iterator muon,muon2;

  // Get the Z collection
  Handle<reco::CompositeCandidateCollection> zll;
  iEvent.getByLabel(theZllLabel_, zll);  

  // Let's look inside the muon collection.
  for (muon = muons->begin();  muon != muons->end(); ++muon){

    h1Cut_->Fill(0);
    //Fill histograms to study Presel and NonPresel Muons
    if (zll->size()==1) { //if1
    
      //Let muon pass the preselection (if2)  
      if ( muon->pt() > 10. && 
	   muon->isGlobalMuon()==true && muon->isTrackerMuon()==true && 
	   muon->normChi2() < 15 && muon->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 && muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0  &&
	   muon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && 
	   muon->dB() < 0.2 && (muon->trackIso()+muon->caloIso()) < 0.15*muon->pt() && muon->numberOfMatches() > 1 && abs(muon->eta()) < 2.1){
	
	int charge_ref = muon->charge();
	h1Cut_->Fill(1);
		
	//Compare muon with any muon (muon2) with same charge
	for (muon2 = muons->begin();  muon2 != muons->end(); ++muon2){

	  if ( muon != muon2 && (charge_ref*(muon2->charge()))>0) { //if 3

	    float deta= muon->eta()-muon2->eta();
	    float dphi= muon->phi()-muon2->phi();
	    float dr= TMath::Sqrt(deta*deta+dphi*dphi);

	    //Fill histograms to compare muon and muon2
	    h1dr_->Fill(dr);
	    h1deta_->Fill(fabs(deta));
	    h1dphi_->Fill(fabs(dphi));
	    h2dphideta_->Fill(fabs(dphi),fabs(deta));
	    h2drpt_->Fill(muon2->pt(),dr);
	    h2dphipt_->Fill(muon2->pt(),fabs(dphi));
	    h2detapt_->Fill(muon2->pt(),fabs(deta));
	    
	    //not selected
	    h1Pt_NotSel_->Fill(muon2->pt());
	    h1Eta_NotSel_->Fill(abs(muon2->eta()));
	    
	    if ( muon2->isTrackerMuon() ) {	     
	      h1NValTkHits_NotSel_->Fill(muon2->innerTrack()->hitPattern().numberOfValidTrackerHits());
	      h1NValPxHits_NotSel_->Fill(muon2->innerTrack()->hitPattern().numberOfValidPixelHits());
	    }
	    if (muon2->isGlobalMuon()){
	      h1Chi2norm_NotSel_->Fill(muon2->normChi2());
	      h1NValMuHits_NotSel_->Fill(muon2->globalTrack()->hitPattern().numberOfValidMuonHits());
	    }
	    
	    h1DB_NotSel_->Fill(muon2->dB());
	    h1TotIsoOverpT_NotSel_->Fill( (muon2->trackIso()+muon2->caloIso())/muon2->pt());
	    h1NofMatch_NotSel_->Fill(muon2->numberOfMatches());
	    
	  }//end if3
	  
	}// end for loop
	
      }//end if2

      //Let muon pass a TkQuality cut 

      if ( muon->isGlobalMuon() && muon->isTrackerMuon() ){ //for GlobalMuon
	if ( muon->normChi2() < 15 && 
	     muon->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 && 
	     muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
	     muon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
	     (muon->trackIso()+muon->caloIso()) < 0.15*muon->pt() && 
	     muon->numberOfMatches() > 1
	     ){
	  h1Cut_->Fill(2);
	  h1Pt_TkQuality_->Fill(muon->pt());
	  h1DB_TkQuality_->Fill(muon->dB());
	  h1Phi_TkQuality_->Fill(abs(muon->phi()));
	  h1Eta_TkQuality_->Fill(abs(muon->eta()));
	  h1Dxy_TkQuality_->Fill(muon->globalTrack()->d0());
	  h1Dz_TkQuality_->Fill(muon->globalTrack()->dz());
	}
      }// end if GlobalMuon
      
      if ( !(muon->isStandAloneMuon()) && muon->isTrackerMuon() ){ //for TrackerMuon
	if ( muon->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 && 
	     muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
	     (muon->trackIso()+muon->caloIso()) < 0.15*muon->pt() && 
	     muon->numberOfMatches() > 1
	     ){
	  h1Cut_->Fill(2);
	  h1Pt_TkQuality_->Fill(muon->pt());
	  h1DB_TkQuality_->Fill(muon->dB());
	  h1Phi_TkQuality_->Fill(abs(muon->phi()));
	  h1Eta_TkQuality_->Fill(abs(muon->eta()));
	  h1Dxy_TkQuality_->Fill(muon->innerTrack()->d0());
	  h1Dz_TkQuality_->Fill(muon->innerTrack()->dz());
	}
      }// end if TrackerMuon
      
      
      if ( !(muon->isTrackerMuon()) ){ //for StandAloneMuon
	if ( muon->outerTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
	     (muon->trackIso()+muon->caloIso()) < 0.15*muon->pt() && 
	     muon->numberOfMatches() > 1
	     ){
	  h1Cut_->Fill(2);
	  h1Pt_TkQuality_->Fill(muon->pt());
	  h1DB_TkQuality_->Fill(muon->dB());
	  h1Phi_TkQuality_->Fill(abs(muon->phi()));
	  h1Eta_TkQuality_->Fill(abs(muon->eta()));
	  h1Dxy_TkQuality_->Fill(muon->outerTrack()->d0());
	  h1Dz_TkQuality_->Fill(muon->outerTrack()->dz());
	}
      }// end if StandAloneMuon

      //Study of preselction failure

      if ( muon->pt() <= 10. ) h1PreselFail_->Fill(0);
      if ( !(muon->isGlobalMuon()) ) h1PreselFail_->Fill(1);
      if ( !(muon->isTrackerMuon()) ) h1PreselFail_->Fill(2);
      if ( muon->isGlobalMuon() && muon->normChi2() >= 15 ) h1PreselFail_->Fill(3);
      if ( muon->isTrackerMuon() && muon->innerTrack()->hitPattern().numberOfValidTrackerHits() <= 10 ) h1PreselFail_->Fill(4);
      if ( muon->isTrackerMuon() && muon->innerTrack()->hitPattern().numberOfValidPixelHits() == 0 ) h1PreselFail_->Fill(5);
      if ( muon->isGlobalMuon() && muon->globalTrack()->hitPattern().numberOfValidMuonHits() == 0 ) h1PreselFail_->Fill(6);
      if ( muon->dB() >= 0.2 ) h1PreselFail_->Fill(7);
      if ( (muon->trackIso()+muon->caloIso()) >= 0.15*muon->pt() ) h1PreselFail_->Fill(8);
      if ( muon->numberOfMatches() <= 1 ) h1PreselFail_->Fill(9);
      if ( abs(muon->eta()) >= 2.1 ) h1PreselFail_->Fill(10);

    }//end if1

  }//end loop on MuonCollection

}

DEFINE_FWK_MODULE(PreselAnalyzer);







