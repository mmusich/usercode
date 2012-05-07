// -*- C++ -*-
//
// Package:    LeptFromDecayFilter
// Class:      LeptFromDecayFilter
// 
/**\class LeptFromDecayFilter LeptFromDecayFilter.cc UserCode/LeptFromDecayFilter/src/LeptFromDecayFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Stefano Casasso,,,
//         Created:  Thu Mar 31 09:15:46 CEST 2011
// $Id: LeptFromDecayFilter.cc,v 1.3 2011/05/04 14:17:41 emiglior Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Parse.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

// includes from Alberto's code
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Common/interface/GetProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimGeneral/HepPDTRecord/interface/PDTRecord.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
//
// class declaration
//
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

using namespace edm;
using namespace std;
using namespace reco;



class LeptFromDecayFilter : public edm::EDFilter {
public:
  explicit LeptFromDecayFilter(const edm::ParameterSet&);
  ~LeptFromDecayFilter();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  //InputTag
  edm::InputTag MuonCollection_;
  edm::InputTag ElectronCollection_;
  edm::InputTag JetCollection_;
  edm::InputTag ZmmCollection_;

  //Histograms

  //GEN
  TH1F *h1PtFromZ_;
  TH1F *h1PtFromB_;
  TH1F *h1PtFromC_;
  TH1F *h1PtFromBC_;
 
  //RECO
  TH1F *h1DrMJMin_;
  TH1F *h1DrMJMinVBTFCut_;
  TH1F *h1PtSoftMuon_;
  TH1F *h1PtLeadingMuon_;
  TH1I *h1Nmuon_;
  TH1I *h1Nelectron_;
  TH1I *h1Njet_;
  TH1I *h1Nzmm_;
  edm::InputTag genParticles_;
  string decay_chain_string;
  bool hasMu_,hasEle_,hasB_,hasC_,hasBC_;
  int nPassedEvts_,nNotPassedEvts_,nTotEvts_;
 
      
      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
LeptFromDecayFilter::LeptFromDecayFilter(const edm::ParameterSet& iConfig)

{
  genParticles_=iConfig.getParameter<InputTag>("src");
  MuonCollection_=iConfig.getParameter<InputTag>("MuonCollection");
  ElectronCollection_=iConfig.getParameter<InputTag>("ElectronCollection");
  JetCollection_=iConfig.getParameter<InputTag>("JetCollection");
  ZmmCollection_=iConfig.getParameter<InputTag>("ZmmCollection");
  decay_chain_string = iConfig.getUntrackedParameter<string>("DecayChainSelection");
  //now do what ever initialization is needed

}


LeptFromDecayFilter::~LeptFromDecayFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
LeptFromDecayFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //Define and Set counters to 0
  int nbMu_(0),ncMu_(0),nbcMu_(0),nbEle_(0),ncEle_(0),nbcEle_(0);

  nTotEvts_++;
  
  // get gen particle candidates
  edm::Handle<GenParticleCollection> genParticlesCollection;
  iEvent.getByLabel(genParticles_, genParticlesCollection);

  for( GenParticleCollection::const_iterator genp = genParticlesCollection->begin();genp != genParticlesCollection->end(); ++ genp ){
    
    //do a copy of the genParticle, not to skip anything in the loop
    const reco::Candidate* p = &(*genp);

    //************************
    //Histograms at GEN level
    //************************
    
    //if final-state muons
    if( abs( p->pdgId() )==13 && ( p->status())==1 ){
      
      //Check muon's mother
      while (p->mother()!=0 && abs(p->mother()->pdgId()) == 13) {
	p = p->mother();
      }


      if ( JetMCTagUtils::decayFromCHadron(*p) && JetMCTagUtils::decayFromBHadron(*p) ) {
	nbcMu_++;
	h1PtFromBC_->Fill(p->pt());
      } else if ( JetMCTagUtils::decayFromBHadron(*p) ) {
	nbMu_++;
	h1PtFromB_->Fill(p->pt());
      } else if ( JetMCTagUtils::decayFromCHadron(*p) ) {
	ncMu_++;
	h1PtFromC_->Fill(p->pt());       
      } else if ( abs(p->mother()->pdgId()) == 23 ) {
	h1PtFromZ_->Fill(p->pt());
      }
    } //end if final state muons

    //if final-state electrons
    if( abs( p->pdgId() )==11 && ( p->status())==1 ){
      
      //Check electron's mother
      while (p->mother()!=0 && abs(p->mother()->pdgId()) == 11) {
	p = p->mother();
      }

      if ( JetMCTagUtils::decayFromCHadron(*p) && JetMCTagUtils::decayFromBHadron(*p) ) {
	nbcEle_++;
	h1PtFromBC_->Fill(p->pt());
      } else if ( JetMCTagUtils::decayFromBHadron(*p) ) {
	nbEle_++;
	h1PtFromB_->Fill(p->pt());
      } else if ( JetMCTagUtils::decayFromCHadron(*p) ) {
	ncEle_++;
	h1PtFromC_->Fill(p->pt());       
      } else if ( abs(p->mother()->pdgId()) == 23 ) {
	h1PtFromZ_->Fill(p->pt());
      }
    }//end if final state electrons
    
  }//end loop on genParticles

  
  //***************
  //Histograms RECO
  //***************

  // get patMuonsWithTrigger
  edm::Handle<pat::MuonCollection> muonsHandle;
  iEvent.getByLabel(MuonCollection_, muonsHandle);
  h1Nmuon_->Fill(muonsHandle->size());
  
  // get patElectrons
  edm::Handle<pat::ElectronCollection> electronsHandle;
  iEvent.getByLabel(ElectronCollection_, electronsHandle);
  h1Nelectron_->Fill(electronsHandle->size());
  
  // get patJets
  edm::Handle<pat::JetCollection> jetsHandle;
  iEvent.getByLabel(JetCollection_, jetsHandle);
  h1Njet_->Fill(jetsHandle->size());
  
  // get Zmm candidates
  Handle<reco::CompositeCandidateCollection> zmmHandle;
  iEvent.getByLabel(ZmmCollection_, zmmHandle);
  h1Nzmm_->Fill(zmmHandle->size());

  
  pat::MuonCollection::const_iterator SoftMuon,LeadingMuon,muon;
  pat::JetCollection::const_iterator jet;

  //Define SoftMuon and LeadingMuon
  if( muonsHandle->size() != 0 ){
    for( muon = muonsHandle->begin();muon != muonsHandle->end(); ++ muon ){
      if (muon == muonsHandle->begin() || muon->pt() < SoftMuon->pt())SoftMuon = muon;
      if (muon == muonsHandle->begin() || muon->pt() > LeadingMuon->pt())LeadingMuon = muon;
    }//end loop on MuonCollection
    h1PtSoftMuon_->Fill(SoftMuon->pt());
    h1PtLeadingMuon_->Fill(LeadingMuon->pt());
  }

  //DeltaR soft muon/jets
  if( muonsHandle->size() != 0 ){
    float DrMJMin_(0);
    if( jetsHandle->size() >= 1 ){
    for( jet = jetsHandle->begin();jet != jetsHandle->end();++jet ){
      float DphiMJ_,DetaMJ_,DrMJ_;
      DphiMJ_ = abs(jet->phi()-SoftMuon->phi());
      DetaMJ_ = abs(jet->eta()-SoftMuon->eta());
      DrMJ_ = TMath::Sqrt(DphiMJ_*DphiMJ_ + DetaMJ_*DetaMJ_);
      if( DrMJMin_ == 0 || DrMJ_ < DrMJMin_ )DrMJMin_ = DrMJ_;
    }//end loop on jet collection
    
    h1DrMJMin_->Fill(DrMJMin_);
    if ( SoftMuon->pt() > 20. && //VBTF cut again
	 SoftMuon->isGlobalMuon()==true && SoftMuon->isTrackerMuon()==true && 
	 SoftMuon->normChi2() < 10 && SoftMuon->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 && SoftMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0  &&
	 SoftMuon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && 
	 SoftMuon->dB() < 0.2 && (SoftMuon->trackIso()+SoftMuon->caloIso()) < 0.15*SoftMuon->pt() && SoftMuon->numberOfMatches() > 1 && abs(SoftMuon->eta()) < 2.1 )h1DrMJMinVBTFCut_->Fill(DrMJMin_);
    else{}
    
    }//end if jets size >= 1
  }//end if muon size != 0

  //************
  //FILTER
  //************
  
  
  //It's time decide whether the event pass the filter or not

  bool decision_(false);

  if ( hasMu_ && !(hasEle_)){//if interested in muons
    if ( hasB_  && nbMu_ != 0 ) decision_ = true;
    if ( hasC_  && ncMu_ != 0 ) decision_ = true;
    if ( hasBC_ && nbcMu_ != 0) decision_ = true;
  }

  if ( hasEle_ && !(hasMu_)){//if interested in electrons
    if ( hasB_  && nbEle_ != 0 ) decision_ = true;
    if ( hasC_  && ncEle_ != 0 ) decision_ = true;
    if ( hasBC_ && nbcEle_ != 0) decision_ = true;
  }
  
  if ( decision_ ){
    nPassedEvts_++;
    LogInfo("Event selected") << "Run: " << iEvent.id().run() << ", Evt: " << iEvent.id().event();
    return true;
  }
  else {
    nNotPassedEvts_++;
    LogInfo("Event skipped") << "Run: " << iEvent.id().run() << ", Evt: " << iEvent.id().event();
    return false;
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
LeptFromDecayFilter::beginJob()
{
  //Directories
  edm::Service<TFileService> fs;
  TFileDirectory HistoPreFilter = fs->mkdir("Histograms before filter");
  TFileDirectory Gen = HistoPreFilter.mkdir("GEN");
  TFileDirectory Reco = HistoPreFilter.mkdir("RECO");

  //GEN
  h1PtFromZ_ = Gen.make<TH1F>("PtFromZ","Pt leptons from Z",100,0,100);
  h1PtFromB_ = Gen.make<TH1F>("PtFromB","Pt leptons from b-hadrons",100,0,100);
  h1PtFromC_ = Gen.make<TH1F>("PtFromC","Pt leptons from c-hadrons",100,0,100);
  h1PtFromBC_ = Gen.make<TH1F>("PtFromBC","Pt leptons from b-c-hadrons",100,0,100);
 
  //RECO
  h1DrMJMin_ = Reco.make<TH1F>("DrMJMin_","DeltaR min soft muon/jet",500,0,10);
  h1DrMJMinVBTFCut_ = Reco.make<TH1F>("DrMJMinVBTFCut","DeltaR min soft muon pt > 20 Gev/jet",500,0,10);
  h1PtSoftMuon_ = Reco.make<TH1F>("PtSoftMuon","Pt soft muon",100,0,100);
  h1PtLeadingMuon_ = Reco.make<TH1F>("PtLeadingMuon","Pt leading muon",100,0,100);
  h1Nmuon_ = Reco.make<TH1I>("Nmuon","Number of muons in event",20,-0.5,19.5);
  h1Nelectron_ = Reco.make<TH1I>("Nelectron","Number of electrons in event",20,-0.5,19.5);
  h1Nzmm_ = Reco.make<TH1I>("Nzmm","Number of zmms in event",20,-0.5,19.5);
  h1Njet_ = Reco.make<TH1I>("Njet","Number of jets in event",20,-0.5,19.5);
  
  
  //Initialize counters
  nTotEvts_ = 0;
  nPassedEvts_ = 0;
  nNotPassedEvts_ = 0;


  //Set booleans to false
  hasMu_  = false;
  hasEle_ = false;
  hasB_   = false;
  hasC_   = false;
  hasBC_  = false;

  std::vector<std::string> tokens = edm::tokenize(decay_chain_string, ">");

  if( tokens.size() != 2 ){
    cout<<"Invalid decay chain string"<<endl;
  }
  if( tokens[1] == "m" ) {
    hasMu_ = true;
  }
  if( tokens[1] == "e" ) {
    hasEle_ = true;
  }
  if ( tokens[0] == "b" ) {
    hasB_ = true;
  }
  else if ( tokens[0] == "c" ) {
    hasC_ = true;
  }
  else if ( tokens[0] == "bc" ) {
    hasBC_ = true;
  }


}

// ------------ method called once each job just after ending the event loop  ------------
void 
LeptFromDecayFilter::endJob() {
  //Print counters results
  cout<<endl<<"Total number of events: "<<nTotEvts_<<endl;
  //cout<<endl<<"Number of events that passed the filter: "<<nPassedEvts_<<endl;
  //cout<<endl<<"Number of events that didn't pass the filter: "<<nNotPassedEvts_<<endl<<endl;

}

//define this as a plug-in
DEFINE_FWK_MODULE(LeptFromDecayFilter);
