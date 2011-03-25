//class JetLeptonCleanAnalyzer

#include "JetLeptonCleanAnalyzer.h"

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


#include "DataFormats/PatCandidates/interface/Jet.h"

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
JetLeptonCleanAnalyzer::JetLeptonCleanAnalyzer(const ParameterSet& pset) : 
  theMuonLabel_(pset.getUntrackedParameter<string>("MuonCollection")),
  theZllLabel_(pset.getUntrackedParameter<string>("ZllCollection")),
  theJetsLabel_(pset.getUntrackedParameter<string>("JetCollection")),
  theJetsPULabel_(pset.getUntrackedParameter<string>("JetCollectionPU")),
  theJetsJPTLabel_(pset.getUntrackedParameter<string>("JetCollectionJPT")),
  theJetsPFLabel_(pset.getUntrackedParameter<string>("JetCollectionPF"))
{
}

/// Destructor
JetLeptonCleanAnalyzer::~JetLeptonCleanAnalyzer(){
}

void JetLeptonCleanAnalyzer::beginJob(){

  // Book histograms
  edm::Service<TFileService> fs;

  h1dr_ = fs->make<TH1F>("dr","DeltaR Lepton/Jet",50,0,5);
  h2drpT_ = fs->make<TH2F>("drpT","pT vs DeltaR Lepton/Jet",50,0,5,75,0,75);

  

}

void JetLeptonCleanAnalyzer::endJob(){

}

void JetLeptonCleanAnalyzer::analyze(const Event & iEvent, const EventSetup& eventSetup)
{
  
  // Get the Muon collection
  Handle<pat::MuonCollection> muons;
  iEvent.getByLabel(theMuonLabel_, muons);

  // Get the Z collection
  Handle<reco::CompositeCandidateCollection> zll;
  iEvent.getByLabel(theZllLabel_, zll);

  //Get Jet collections
  Handle<pat::JetCollection> jets;
  iEvent.getByLabel(theJetsLabel_, jets);

  Handle<pat::JetCollection> jetsPU;
  iEvent.getByLabel(theJetsPULabel_, jetsPU);

  Handle<pat::JetCollection> jetsJPT;
  iEvent.getByLabel(theJetsJPTLabel_, jetsJPT);

  Handle<pat::JetCollection> jetsPF;
  iEvent.getByLabel(theJetsPFLabel_, jetsPF);




  pat::MuonCollection::const_iterator muon;
  pat::JetCollection::const_iterator jet;
  
  if( muons->size() > 2 ){
    
    //Loop on muon collection
    for (muon = muons->begin();  muon != muons->end(); ++muon){
      
      if ( muon->pt() > 10. && 
	   muon->isGlobalMuon()==true && muon->isTrackerMuon()==true && 
	   muon->normChi2() < 15 && muon->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 && muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0  &&
	   muon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && 
	   muon->dB() < 0.2 && (muon->trackIso()+muon->caloIso()) < 0.15*muon->pt() && muon->numberOfMatches() > 1 && abs(muon->eta()) < 2.1){}
      else {
	float etaMu,phiMu,etaJet,phiJet,deta,dphi,dr;
	etaMu = muon->eta();
	phiMu = muon->phi();

	float drMin=FLT_MAX;
	for(jet = jetsPF->begin(); jet != jetsPF->end(); ++jet){

	  etaJet = jet->eta();
	  phiJet = jet->phi();
	  deta = etaMu - etaJet;
	  dphi = phiMu - phiJet;
	  dr = TMath::Sqrt(deta*deta + dphi*dphi);
	  if ( dr<drMin ) drMin=dr; 
	}
	h1dr_->Fill(drMin);
	h2drpT_->Fill(drMin,muon->pt());
	
      }
      
         
    }//end loop on MuonCollection

  }//end if muons size

  
}



DEFINE_FWK_MODULE(JetLeptonCleanAnalyzer);







