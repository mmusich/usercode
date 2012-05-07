#include <map>
#include <string>

#include "TH1.h"
#include <Math/VectorUtil.h>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

// for "evt filters"
#include "DataFormats/Common/interface/MergeableCounter.h"

class PATAnalyzer : public edm::EDAnalyzer {

public:
  /// default constructor
  explicit PATAnalyzer(const edm::ParameterSet&);
  /// default destructor
  ~PATAnalyzer();
  
private:
  /// everything that needs to be done before the event loop
  virtual void beginJob() ;
  /// everything that needs to be done during the event loop
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  /// check if histogram was booked
  bool booked(const std::string histName) const { return histContainer_.find(histName.c_str())!=histContainer_.end(); };
  /// fill histogram if it had been booked before
  void fill(const std::string histName, double value) const { if(booked(histName.c_str())) histContainer_.find(histName.c_str())->second->Fill(value); };
  /// sort the candidates by difference
  std::vector<reco::CompositeCandidate> sortCandidatesByDifference(std::vector<reco::CompositeCandidate> unsortedCands);

  /// check the jet ID
  Int_t jetId(pat::Jet jet);

  /// for Event yields
  void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup&);

  /// everything that needs to be done after the event loop
  virtual void endJob() ;
 
  // simple map to contain all histograms; 
  // histograms are booked in the beginJob() 
  // method
  std::map<std::string,TH1F*> histContainer_; 
  std::vector<TH1F*> h_ptJetsCorrected;

  // Filter studies
  std::vector<std::string> numEventsNames_;

  // input tags  
  edm::InputTag metSrc_;
  edm::InputTag elecSrc_;
  edm::InputTag muonSrc_;
  edm::InputTag jetSrc_;
  edm::InputTag tauSrc_;
  edm::InputTag photonSrc_;
  edm::InputTag theZmmSrc_;
  edm::InputTag theZeeSrc_;

  /// correction levels for pat jet
  std::vector<std::string> str_corrLevels;
  std::vector<std::string> corrections_;
  Int_t corrsize;

  // for trend plots 	 
  unsigned int runNumberFirst_, runNumberLast_;

};

PATAnalyzer::PATAnalyzer(const edm::ParameterSet& iConfig):
  histContainer_(),
  //photonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc")),
  //tauSrc_(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc" )),
  numEventsNames_(iConfig.getUntrackedParameter< std::vector<std::string> >("numEventsNames")),
  metSrc_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc")),                 // TO BE FIXED
  elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc")),
  muonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc")),
  jetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc" )),
  theZmmSrc_(iConfig.getUntrackedParameter<edm::InputTag>("ZmmSrc")),
  theZeeSrc_(iConfig.getUntrackedParameter<edm::InputTag>("ZeeSrc")),
  runNumberFirst_(iConfig.getUntrackedParameter<unsigned int>("firstRunNumber",0)), 	 
  runNumberLast_(iConfig.getUntrackedParameter<unsigned int>("lastRunNumber",0))
{
  // read and parse corrLevels
  std::vector<std::string> str_corrLevels = iConfig.getParameter<std::vector<std::string> >("corrLevels");
  corrections_.reserve(str_corrLevels.size());
  for (std::vector<std::string>::const_iterator it = str_corrLevels.begin(), ed = str_corrLevels.end(); it != ed; ++it) {
    corrections_.push_back(*it);
  }
  corrsize=corrections_.size();
}

PATAnalyzer::~PATAnalyzer()
{
}

void
PATAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  // get met collection  
  edm::Handle<edm::View<pat::MET> > mets;
  iEvent.getByLabel(metSrc_,mets);

  // get electron collection
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByLabel(elecSrc_,electrons);

  // get muon collection
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(muonSrc_,muons);

  // get jet collection
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jetSrc_,jets);
  
  // get tau collection  
  // edm::Handle<edm::View<pat::Tau> > taus;
  //iEvent.getByLabel(tauSrc_,taus);
  
  // get photon collection  
  //edm::Handle<edm::View<pat::Photon> > photons;
  //iEvent.getByLabel(photonSrc_,photons);


  edm::Handle<edm::View<reco::Vertex> > vertexCollection;
  iEvent.getByLabel("offlinePrimaryVertices", vertexCollection);

  unsigned int vertexCollectionSize = vertexCollection->size();  
  int nvvertex = 0;
  for (unsigned int i=0; i<vertexCollectionSize; i++) {
    const reco::Vertex& vertex = vertexCollection->at(i);
    if (vertex.isValid()) nvvertex++;
  }

  fill("nVvertex",nvvertex);
  
  // get rho from event
//   double _rho=-999.;
//   edm::Handle<double> rhoH;
//   iEvent.getByLabel("kt6PFJetsForIso","rho",rhoH);
//   _rho = *rhoH;
  
//      __  __                      _     _      
//     |  \/  |_  _ ___ _ _    _ __| |___| |_ ___
//     | |\/| | || / _ \ ' \  | '_ \ / _ \  _(_-<
//     |_|  |_|\_,_\___/_||_| | .__/_\___/\__/__/
//                            |_|                

  // loop over muons
  for (edm::View<pat::Muon>::const_iterator muon = muons->begin();  muon != muons->end(); ++muon){
    fill("muon_pt",muon->pt());  
    fill("muon_eta",muon->eta());
    fill("muon_phi",muon->phi());
    fill("muon_dB",muon->dB());
    double muIP      = fabs(muon->dB(pat::Muon::PV3D));
    double muIPError = muon->edB(pat::Muon::PV3D); 
    double SIP =  muIP/muIPError;
    fill("muon_IP",muIP);
    fill("muon_SIP",SIP);
    if(muon->isGlobalMuon()==true || muon->isTrackerMuon()==true){
      fill("muon_iso",(muon->trackIso() + muon->caloIso())/muon->pt());
      fill("isGlobalMuon",muon->isGlobalMuon()); 
      fill("isTrackerMuon",muon->isTrackerMuon()); 
      fill("muon_trackerhits",muon->innerTrack()->hitPattern().numberOfValidTrackerHits());
      fill("muon_pixelhits",muon->innerTrack()->hitPattern().numberOfValidPixelHits());
      if(muon->isGlobalMuon()==true){
	fill("muon_chi2",muon->globalTrack()->normalizedChi2());
	fill("muon_muonhits",muon->globalTrack()->hitPattern().numberOfValidMuonHits());
	fill("muon_numberOfMatches",muon->numberOfMatches());
      }
    }
  }
  
  //std::cout<<"after muon plots"<<std::endl;

  edm::View<pat::Muon>::const_iterator muon1,muon2;
  reco::Candidate::LorentzVector p1,p2,sum;

  //Fill histogram AllDiMuMass
  if(muons->size()>1){
    
    for (muon1 = muons->begin(); muon1 != muons->end(); muon1++){
      for (muon2 = muon1+1; muon2 != muons->end(); muon2++) {
	p1 = muon1->p4();
	p2 = muon2->p4();  
	sum = p1 + p2;
	if((muon1->charge()*muon2->charge())<0){ 
	  fill("AllDiMuMass",sum.mass());
	  fill("AllDiMuPt",sum.pt());
	  fill("AllDiMuEta",sum.eta());
	  fill("AllDiMuPhi",sum.phi());

	  double deltaR = ROOT::Math::VectorUtil::DeltaR(muon1->momentum(),muon2->momentum());
	  double deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(muon1->momentum(),muon2->momentum());
	  fill("AllDiMuDeltaPhi",deltaPhi);
	  fill("AllDiMuDeltaR",deltaR);
	}
      }						 
    }
  }
  
  //Fill histogram PreselMuMass
  if(muons->size()>1){
    
    for (muon1 = muons->begin(); muon1 != muons->end(); muon1++){
      
      if ( muon1->pt() > 5. && 
	   muon1->isGlobalMuon()==true && muon1->isTrackerMuon()==true && 
	   muon1->normChi2() < 10 && muon1->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 && muon1->innerTrack()->hitPattern().numberOfValidPixelHits() > 0  &&
	   muon1->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && 
	   fabs(muon1->dB()) < 0.02 && (muon1->trackIso()+muon1->caloIso()) < 0.15*muon1->pt() && muon1->numberOfMatches() > 1 && abs(muon1->eta()) < 2.1){
	
	for (muon2 = muon1+1; muon2 != muons->end(); muon2++) {
	
	  if ( muon2->pt() > 5. && 
	       muon2->isGlobalMuon()==true && muon2->isTrackerMuon()==true && 
	       muon2->normChi2() < 10 && muon2->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 && muon2->innerTrack()->hitPattern().numberOfValidPixelHits() > 0  &&
	       muon2->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && 
	       fabs(muon2->dB()) < 0.02 && (muon2->trackIso()+muon2->caloIso()) < 0.15*muon2->pt() && muon2->numberOfMatches() > 1 && abs(muon2->eta()) < 2.1){
	    p1 = muon1->p4();
	    p2 = muon2->p4();  
	    sum = p1 + p2;
	    if((muon1->charge()*muon2->charge())<0){ 
	      fill("PreselDiMuMass",sum.mass());
	      fill("PreselDiMuPt",sum.pt());
	      fill("PreselDiMuEta",sum.eta());
	      fill("PreselDiMuPhi",sum.phi());
	      
	      // check the Delta R between the two muons
	      // (Delta_R^2 = Delta_Eta^2 + Delta_Phi^2)
	      double deltaR = ROOT::Math::VectorUtil::DeltaR(muon1->momentum(),muon2->momentum());
	      double deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(muon1->momentum(),muon2->momentum());
	      fill("PreselDiMuDeltaPhi",deltaPhi);
	      fill("PreselDiMuDeltaR",deltaR);
	      
	      if(sum.mass()>50){
		fill("MassCutDiMuMass",sum.mass());
		fill("MassCutDiMuPt",sum.pt());
		fill("MassCutDiMuEta",sum.eta());
		fill("MassCutDiMuPhi",sum.phi());
		fill("MassCutDiMuDeltaPhi",deltaPhi);
		fill("MassCutDiMuDeltaR",deltaR); 
		
	      } //if mass cut 
	    } // if charge ok
	  } //is second muon presel
	} // loop on second muon 
      } //if first muon presel
    } //loop on first muon
  } //if muon collection size > 1
 
  //std::cout<<"after muon mass plots"<<std::endl;

//      _     _          _     _      
//   _ | |___| |_   _ __| |___| |_ ___
//  | || / -_)  _| | '_ \ / _ \  _(_-<
//   \__/\___|\__| | .__/_\___/\__/__/
//                 |_|                

  // loop over jets
  size_t nJets=0;
  for(edm::View<pat::Jet>::const_iterator jet = jets->begin();  jet!=jets->end(); ++jet){
    // fill basic kinematics for all jets
    fill( "jet_pt" , jet->pt());
    for(Int_t j = 0; j <  corrsize; j++) {
      h_ptJetsCorrected[j]->Fill(jet->correctedJet(corrections_[j]).pt());
    }
    fill("jet_eta",jet->eta());
    fill("jet_phi",jet->phi());
    fill("jet_nef",jet->neutralEmEnergyFraction());
    fill("jet_nhf",(jet->neutralHadronEnergy() + jet->HFHadronEnergy() ) / jet->energy());
    fill("jet_chf",jet->chargedHadronEnergyFraction());
    fill("jet_nch",jet->chargedMultiplicity());
    fill("jet_cef",jet->chargedEmEnergyFraction()); 
    fill("jet_NConstituents",jet->numberOfDaughters());
    fill("jet_MuonMultiplicity",jet->muonMultiplicity()); 
    fill("jet_id",this->jetId(*jet));

    // special for CaloJets 
    // fill("jet_emf",jet->emEnergyFraction() );
    // fill("jet_Area",jet->jetArea ());
    // fill("jet_NCarrying60",jet->n60());
    // fill("jet_NCarrying90",jet->n90());

    double discrTC = jet->bDiscriminator("trackCountingHighEffBJetTags");
    double discrSSVHP = jet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
    double discrSSVHE = jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    double discrCSV = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
    fill("jet_discrTC", discrTC);
    fill("jet_discrSSVHP", discrSSVHP);
    fill("jet_discrSSVHE", discrSSVHE);
    fill("jet_discrCSV", discrCSV);
  
    // multiplicity of pt>15 jets
    if(jet->pt()>15){ 
      ++nJets;
    } 
  }

  //std::cout<<"after jet loop"<<std::endl;

  if(jets->size()>0){
    // fill basic kinematics for leading jets
    fill( "leadingjet_pt" , (*jets)[0].pt());
    fill( "leadingjet_eta", (*jets)[0].eta());
    fill( "leadingjet_phi", (*jets)[0].phi());
    
    //std::cout<<"after leading jet plots"<<std::endl;
    
    // fill basic kinematics for sub-leading jets (if any)
    if(jets->size()>1){
      fill( "di_jet_mass", ((*jets)[0].p4()+(*jets)[1].p4()).mass());  
      fill( "subleadingjet_pt" , (*jets)[1].pt());
      fill( "subleadingjet_eta", (*jets)[1].eta());
      fill( "subleadingjet_phi", (*jets)[1].phi());
    }
  }
  
  //std::cout<<"after jet plots"<<std::endl;

  // do something similar for the other candidates
  //histContainer_["photons"]->Fill(photons->size() );
  //histContainer_["taus" ]->Fill(taus->size()  );
  histContainer_["met"]->Fill(mets->empty() ? 0 : (*mets)[0].et());             // TO BE FIXED!!!
  if(mets->size()>0){
    histContainer_["met_phi"]->Fill((*mets)[0].phi()); 
  }
  histContainer_["jets"]->Fill(nJets);
  histContainer_["elecs" ]->Fill(electrons->size());
  histContainer_["muons"]->Fill(muons->size() );

  //std::cout<<"after collection plots"<<std::endl;

//   ___ _        _                       _     _      
//  | __| |___ __| |_ _ _ ___ _ _    _ __| |___| |_ ___
//  | _|| / -_) _|  _| '_/ _ \ ' \  | '_ \ / _ \  _(_-<
//  |___|_\___\__|\__|_| \___/_||_| | .__/_\___/\__/__/
//                                  |_|                

  // loop over electrons
  for(edm::View<pat::Electron>::const_iterator elec=electrons->begin(); elec !=electrons->end(); ++elec){
    // fill simple histograms
    histContainer_["electron_pt"]->Fill(elec->pt());
    histContainer_["electron_eta"]->Fill(elec->eta());
    histContainer_["electron_phi"]->Fill(elec->phi());
    histContainer_["electron_iso"]->Fill((elec->trackIso()+elec->caloIso())/elec->pt() );
    histContainer_["electron_eop"]->Fill(elec->eSeedClusterOverP());
    histContainer_["electron_clus"]->Fill(elec->e1x5()/elec->e5x5());
    fill("electron_dB",elec->dB());
    double eleIP      = fabs(elec->dB(pat::Electron::PV3D));
    double eleIPError = elec->edB(pat::Electron::PV3D); 
    double SIP =  eleIP/eleIPError;
    fill("electron_IP",eleIP);
    fill("electron_SIP",SIP);
    
    //fill electron id histograms
    if( elec->electronID("eidVBTFCom80") > 0.5 )
      histContainer_["electron_eIDs" ]->Fill(0);
    if( elec->electronID("eidVBTFCom95") > 0.5 )
      histContainer_["electron_eIDs" ]->Fill(1);
    if( elec->electronID("eidVBTFRel80") > 0.5 )
      histContainer_["electron_eIDs" ]->Fill(2);
    if( elec->electronID("eidVBTFRel95") > 0.5 )
      histContainer_["electron_eIDs" ]->Fill(3);
  }
  
  //std::cout<<"after electron plots"<<std::endl;

  edm::View<pat::Electron>::const_iterator electron1,electron2;
  reco::Candidate::LorentzVector ep1,ep2,esum;

  //Fill histogram AllDiEleMass
  if(electrons->size()>1){
    
    for (electron1 = electrons->begin(); electron1 != electrons->end(); electron1++){
      for (electron2 = electron1+1; electron2 != electrons->end(); electron2++) {
	ep1 = electron1->p4();
	ep2 = electron2->p4();  
	esum = ep1 + ep2;
	if((electron1->charge()*electron2->charge())<0){
	  fill("AllDiEleMass",esum.mass());
	  fill("AllDiElePt",esum.pt());
	  fill("AllDiEleEta",esum.eta());
	  fill("AllDiElePhi",esum.phi());

	  double deltaR = ROOT::Math::VectorUtil::DeltaR(electron1->momentum(),electron2->momentum());
	  double deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(electron1->momentum(),electron2->momentum());
	  fill("AllDiEleDeltaPhi",deltaPhi);
	  fill("AllDiEleDeltaR",deltaR);
	}
      }						 
    }
  }
  
  //Fill histogram PreselEleMass
  if(electrons->size()>1){
    
    for (electron1 = electrons->begin(); electron1 != electrons->end(); electron1++){
      
      if ( electron1->pt() > 5.0 && abs(electron1->eta()) < 2.5 && (electron1->isEE() || electron1->isEB()) && !electron1->isEBEEGap() && electron1->electronID("eidVBTFRel85") == 7  && (fabs(electron1->superCluster()->eta())<1.444 || fabs(electron1->superCluster()->eta())>1.566 )){
	
	for (electron2 = electron1+1; electron2 != electrons->end(); electron2++) {
	
	  if (  electron2->pt() > 5.0 && abs(electron2->eta()) < 2.5 && (electron2->isEE() || electron2->isEB()) && !electron2->isEBEEGap() && electron2->electronID("eidVBTFRel85") == 7 && (fabs(electron2->superCluster()->eta())<1.444 || fabs(electron2->superCluster()->eta())>1.566 )){
	    ep1 = electron1->p4();
	    ep2 = electron2->p4();  
	    esum = ep1 + ep2;
	    if((electron1->charge()*electron2->charge())<0){ 
	      fill("PreselDiEleMass",esum.mass());
	      fill("PreselDiElePt",esum.pt());
	      fill("PreselDiEleEta",esum.eta());
	      fill("PreselDiElePhi",esum.phi());

	      double deltaR = ROOT::Math::VectorUtil::DeltaR(electron1->momentum(),electron2->momentum());
	      double deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(electron1->momentum(),electron2->momentum());
	      fill("PreselDiEleDeltaPhi",deltaPhi);
	      fill("PreselDiEleDeltaR",deltaR);
	      
	      if(esum.mass()>50){
		fill("MassCutDiEleMass",esum.mass());
		fill("MassCutDiElePt",esum.pt());
		fill("MassCutDiEleEta",esum.eta());
		fill("MassCutDiElePhi",esum.phi());
		fill("MassCutDiEleDeltaPhi",deltaPhi);
		fill("MassCutDiEleDeltaR",deltaR); 
	      }
	    }
	  }
	}
      }
    }
  }

  //std::cout<<"after electron mass plots"<<std::endl;

//   ___              _ _    _      _                      _      _    _        
//  / __|__ _ _ _  __| (_)__| |__ _| |_ ___  __ ____ _ _ _(_)__ _| |__| |___ ___
// | (__/ _` | ' \/ _` | / _` / _` |  _/ -_) \ V / _` | '_| / _` | '_ \ / -_|_-<
//  \___\__,_|_||_\__,_|_\__,_\__,_|\__\___|  \_/\__,_|_| |_\__,_|_.__/_\___/__/

  unsigned int runNumber=iEvent.id().run(); 	 
  histContainer_["trendplot"]->Fill(runNumber);

  // Get the Z->mm collection
  edm::Handle<reco::CompositeCandidateCollection> zmmHandle;
  iEvent.getByLabel(theZmmSrc_, zmmHandle);
  const reco::CompositeCandidateCollection & zmm = *(zmmHandle.product());
  fill("NZmm",zmm.size()); 
  if ( zmm.size()>0 ) histContainer_["trendplotMu"]->Fill(runNumber);

  std::vector<reco::CompositeCandidate> UnsortedZmm;
  
  for(reco::CompositeCandidateCollection::const_iterator ZmmCandidate=zmm.begin(); ZmmCandidate != zmm.end(); ++ZmmCandidate){
    UnsortedZmm.push_back(*ZmmCandidate);
  }

  if(UnsortedZmm.size()>0){

    std::vector<reco::CompositeCandidate> ZmmCandidatesSorted = this->sortCandidatesByDifference(UnsortedZmm);
    
    fill("ZmmMass",ZmmCandidatesSorted[0].mass());  
    fill("ZmmPt",ZmmCandidatesSorted[0].pt());            
    fill("ZmmEta",ZmmCandidatesSorted[0].eta());                
    fill("ZmmPhi",ZmmCandidatesSorted[0].phi());            
    fill("ZmmLeadingLepPt",ZmmCandidatesSorted[0].daughter(0)->pt());      
    fill("ZmmLeadingLepEta",ZmmCandidatesSorted[0].daughter(0)->eta());      
    fill("ZmmLeadingLepPhi",ZmmCandidatesSorted[0].daughter(0)->phi());      
    fill("ZmmSubLeadingLepPt",ZmmCandidatesSorted[0].daughter(1)->pt());     
    fill("ZmmSubLeadingLepEta",ZmmCandidatesSorted[0].daughter(1)->eta());     
    fill("ZmmSubLeadingLepPhi",ZmmCandidatesSorted[0].daughter(1)->phi());    
  }

  //   // Get the Z->ee collection
  edm::Handle<reco::CompositeCandidateCollection> zeeHandle;
  iEvent.getByLabel(theZeeSrc_, zeeHandle);
  const reco::CompositeCandidateCollection & zee = *(zeeHandle.product());
  fill("NZee",zee.size());
  if ( zee.size()>0 ) histContainer_["trendplotEle"]->Fill(runNumber);

  std::vector<reco::CompositeCandidate> UnsortedZee;

  for(reco::CompositeCandidateCollection::const_iterator ZeeCandidate=zee.begin(); ZeeCandidate != zee.end(); ++ZeeCandidate){
    UnsortedZee.push_back(*ZeeCandidate);
  }
  
  if(UnsortedZee.size()>0){

    std::vector<reco::CompositeCandidate> ZeeCandidatesSorted = this->sortCandidatesByDifference(UnsortedZee);
    
    fill("ZeeMass",ZeeCandidatesSorted[0].mass());  
    fill("ZeePt",ZeeCandidatesSorted[0].pt());            
    fill("ZeeEta",ZeeCandidatesSorted[0].eta());                
    fill("ZeePhi",ZeeCandidatesSorted[0].phi());            
    fill("ZeeLeadingLepPt",ZeeCandidatesSorted[0].daughter(0)->pt());      
    fill("ZeeLeadingLepEta",ZeeCandidatesSorted[0].daughter(0)->eta());      
    fill("ZeeLeadingLepPhi",ZeeCandidatesSorted[0].daughter(0)->phi());      
    fill("ZeeSubLeadingLepPt",ZeeCandidatesSorted[0].daughter(1)->pt());     
    fill("ZeeSubLeadingLepEta",ZeeCandidatesSorted[0].daughter(1)->eta());     
    fill("ZeeSubLeadingLepPhi",ZeeCandidatesSorted[0].daughter(1)->phi());     
  }
  
}

std::vector<reco::CompositeCandidate> PATAnalyzer::sortCandidatesByDifference(std::vector<reco::CompositeCandidate> unsortedCands){

  std::vector<reco::CompositeCandidate> sortedCands = unsortedCands;
  Double_t mZ=91.1876;
  std::vector<Double_t> diffZmass; 
  unsigned int ZCandSize=unsortedCands.size();
  

  if(unsortedCands.size()!=1){
    std::cout<<"Z candidates: "<<ZCandSize<<std::endl;
  }
  
  for(reco::CompositeCandidateCollection::const_iterator ZllCandidate1=unsortedCands.begin(); ZllCandidate1 != unsortedCands.end(); ++ZllCandidate1){
    diffZmass.push_back(fabs(mZ - ZllCandidate1->mass()));
  }

  for (unsigned int i = 0; i < ZCandSize; i++) {   
    for (unsigned int j = i+1; j < ZCandSize; j++) {

      if(diffZmass[i] > diffZmass[j]) {
	
	reco::CompositeCandidate auxCand = sortedCands[i];
	sortedCands[i] = sortedCands[j];
	sortedCands[j] = auxCand;  
  
      }
    }    
  }
 
  return sortedCands;
}

void 
PATAnalyzer::beginJob()
{
  // register to the TFileService
  edm::Service<TFileService> fs;

  // book histograms:

  // General Object Multiplicity
  // histContainer_["photons"]=fs->make<TH1F>("photons","photon multiplicity",10, 0,10);
  // histContainer_["taus"]=fs->make<TH1F>("taus","tau multiplicity",10,0, 10);

  TFileDirectory SelEvents = fs->mkdir("Selevents");
  
  unsigned int nbin = numEventsNames_.size();
  if ( nbin > 0 ) {
    histContainer_["SelectedEvtsAfterFilter_"]=SelEvents.make<TH1F>("SelectedEvts","Selected events",nbin,-0.5,-0.5+nbin);  
    for (unsigned int i = 0;i < nbin; i++) histContainer_["SelectedEvtsAfterFilter_"]->GetXaxis()->SetBinLabel(i+1,(numEventsNames_[i].c_str()));
  }
  
  TFileDirectory EventVars = fs->mkdir("EventVars");

  int nBin=runNumberLast_-runNumberFirst_+1; // assume that TH1F for trend plots are rebinned in the ROOT script used to display them...
  histContainer_["trendplot"]   =EventVars.make<TH1F>("evtvsrun","n evt vs. run; run number",nBin,(float)runNumberFirst_-0.5,(float)runNumberLast_+0.5);
  histContainer_["trendplotMu"] =EventVars.make<TH1F>("evtvsrun_mu","n #mu #mu evt vs. run; run number",nBin,(float)runNumberFirst_-0.5,(float)runNumberLast_+0.5);
  histContainer_["trendplotEle"]=EventVars.make<TH1F>("evtvsrun_ele","n #ele #ele evt vs. run; run number",nBin,(float)runNumberFirst_-0.5,(float)runNumberLast_+0.5);

  histContainer_["met"]     =EventVars.make<TH1F>("met","missing E_{T}; #slash{E}_{T} (GeV)",30,0,150);                             // TO BE FIXED!!!!
  histContainer_["met_phi"] =EventVars.make<TH1F>("METphi","MET #phi; #phi_{#slash{E}_{T}} (rad)",70,-TMath::Pi(),TMath::Pi());
  histContainer_["elecs"]   =EventVars.make<TH1F>("n_elecs","electron multiplicity; n_{electrons};events",15,-0.5,14.5);
  histContainer_["muons"]   =EventVars.make<TH1F>("n_muons","muon multiplicity; n_{muons};events",15,-0.5,14.5);
  histContainer_["jets"]    =EventVars.make<TH1F>("n_jets","jet multiplicity; n_{jets};events",15,-0.5,14.5);
  histContainer_["nVvertex"]=EventVars.make<TH1F>("NVvertex","Nb. of valid vertices",26,-0.5,25.5);

  TFileDirectory ElectronVars = fs->mkdir("ElectronVars");

  // electron variables (for all electrons)
  histContainer_["electron_pt"  ]=ElectronVars.make<TH1F>("electron_pt","electron pt;p_{T} (GeV)",150,0.,150.);
  histContainer_["electron_eta" ]=ElectronVars.make<TH1F>("electron_eta","electron eta; electron #eta",30,-3.,3.);
  histContainer_["electron_phi" ]=ElectronVars.make<TH1F>("electron_phi","electron phi; electron #phi (rad);",60,-TMath::Pi(),TMath::Pi());
  histContainer_["electron_iso" ]=ElectronVars.make<TH1F>("electron_iso","electron iso; (#sum_{trk} p_{T} + #sum_{cal} E_{T})/p_{T} (combRelIso);",60,0.,20.);
  histContainer_["electron_eop" ]=ElectronVars.make<TH1F>("electron_eop","electron eop; electron E/p;",40,0.,1.);
  histContainer_["electron_clus"]=ElectronVars.make<TH1F>("electron_clus","electron clus; electron E_{1x5}/E_{5x5}",40,0.,1.);
  histContainer_["electron_eIDs"]=ElectronVars.make<TH1F>("electron_eIDs","electron eIDS; electron eleID;",4,0.,4.);
  histContainer_["electron_IP"]  =ElectronVars.make<TH1F>("electron_IP","electron IP (cm);electron IP (cm)",100,0.,5.);
  histContainer_["electron_dB"]  =ElectronVars.make<TH1F>("electron_dB","electron dB; dB [cm]",100,0.,5.);
  histContainer_["electron_SIP"] =ElectronVars.make<TH1F>("electron_SIP","electron IP/#sigma_{IP};electron IP/#sigma_{IP}",100,0.,10.);
  histContainer_["AllDiEleMass"] = ElectronVars.make<TH1F>("AllDiEleMass","Invariant mass spectrum of all di-electrons; M(e^{+}e^{-}) (GeV)",100,0,120);
  histContainer_["AllDiElePt"] = ElectronVars.make<TH1F>("AllDiElePt","p_{T} of all di-electrons; p_{T}(e^{+}e^{-}) (GeV)",30,0,150);
  histContainer_["AllDiEleEta"] = ElectronVars.make<TH1F>("AllDiEleEta","#eta of all di-electrons; #eta(e^{+}e^{-})",50,-5,5);
  histContainer_["AllDiElePhi"] = ElectronVars.make<TH1F>("AllDiElePhi","#phi of all di-electrons; #phi(e^{+}e^{-}) (rad)",30,-TMath::Pi(),TMath::Pi());  
  histContainer_["AllDiEleDeltaPhi"] = ElectronVars.make<TH1F>("AllDiEleDeltaPhi","#Delta #phi of all di-electrons; #Delta#phi(e^{+}e^{-}) (rad)",60,0.,TMath::Pi());  
  histContainer_["AllDiEleDeltaR"] = ElectronVars.make<TH1F>("AllDiEleDeltaR","#Delta R of all di-electrons; #DeltaR(e^{+}e^{-}) (rad)",100,0,10); 

  histContainer_["PreselDiEleMass"] = ElectronVars.make<TH1F>("PreselDiEleMass","Invariant mass spectrum of Presel di-electrons; M(e^{+}e^{-}) (GeV)",100,0,120);
  histContainer_["PreselDiElePt"] = ElectronVars.make<TH1F>("PreselDiElePt","p_{T} of Presel di-electrons; p_{T}(e^{+}e^{-}) (GeV)",30,0,150);
  histContainer_["PreselDiEleEta"] = ElectronVars.make<TH1F>("PreselDiEleEta","#eta of Presel di-electrons; #eta(e^{+}e^{-})",50,-5,5);
  histContainer_["PreselDiElePhi"] = ElectronVars.make<TH1F>("PreselDiElePhi","#phi of Presel di-electrons; #phi(e^{+}e^{-}) (rad)",30,-TMath::Pi(),TMath::Pi());  
  histContainer_["PreselDiEleDeltaPhi"] = ElectronVars.make<TH1F>("PreselDiEleDeltaPhi","#Delta #phi of Presel di-electrons; #Delta#phi(e^{+}e^{-}) (rad)",60,0.,TMath::Pi());  
  histContainer_["PreselDiEleDeltaR"] = ElectronVars.make<TH1F>("PreselDiEleDeltaR","#Delta R of Presel di-electrons; #DeltaR(e^{+}e^{-}) (rad)",100,0,10);  

  histContainer_["MassCutDiEleMass"] = ElectronVars.make<TH1F>("MassCutDiEleMass","Invariant mass spectrum of M(ll)> 50 GeV di-electrons; M(e^{+}e^{-}) (GeV)",100,50.,150.);
  histContainer_["MassCutDiElePt"] = ElectronVars.make<TH1F>("MassCutDiElePt","p_{T} of M(ll)> 50 GeV di-electrons; p_{T}(e^{+}e^{-}) (GeV)",30,0,150);
  histContainer_["MassCutDiEleEta"] = ElectronVars.make<TH1F>("MassCutDiEleEta","#eta of M(ll)> 50 GeV di-electrons; #eta(e^{+}e^{-})",50,-5,5);
  histContainer_["MassCutDiElePhi"] = ElectronVars.make<TH1F>("MassCutDiElePhi","#phi of M(ll)> 50 GeV di-electrons; #phi(e^{+}e^{-}) (rad)",30,-TMath::Pi(),TMath::Pi());  
  histContainer_["MassCutDiEleDeltaPhi"] = ElectronVars.make<TH1F>("MassCutDiEleDeltaPhi","#Delta #phi of M(ll)> 50 GeV di-electrons; #Delta#phi(e^{+}e^{-}) (rad)",60,0.,TMath::Pi());  
  histContainer_["MassCutDiEleDeltaR"] = ElectronVars.make<TH1F>("MassCutDiEleDeltaR","#Delta R of M(ll)> 50 GeV di-electrons; #DeltaR(e^{+}e^{-}) (rad)",100,0,10);  

  TFileDirectory JetVars = fs->mkdir("JetVars");

  // jet variables (for all jets)
  histContainer_["jet_pt"]=JetVars.make<TH1F>("alljet_pt","p_{T}(Jet) all jets;p_{T}(Jet) (GeV)",60,0.,300.);
  histContainer_["jet_eta"]=JetVars.make<TH1F>("alljet_eta","#eta (Jet) all jets;#eta (Jet)",60,-3.,3.);
  histContainer_["jet_phi"]=JetVars.make<TH1F>("alljet_phi","#phi (Jet) all jets;#phi (Jet)",60,-TMath::Pi(),TMath::Pi());
  histContainer_["jet_nhf"]=JetVars.make<TH1F>("nhf","nhf;jet neutral Hadron fraction",40,0.,1.);
  histContainer_["jet_nef"]=JetVars.make<TH1F>("nef","nhf;jet neutral EM fraction",40,0.,1.);
  histContainer_["jet_nch"]=JetVars.make<TH1F>("nch","nhf;charged multiplicity",40,-0.5,39.5);
  histContainer_["jet_chf"]=JetVars.make<TH1F>("chf","nhf;charged Hadron fraction",40,0.,1.);
  histContainer_["jet_cef"]=JetVars.make<TH1F>("cef","nhf;charged EM fraction",40,0.,1.);
  histContainer_["jet_NConstituents"]= JetVars.make<TH1F>("JetNConstituents","Jet n. of constituents;number of jet constituents",100,-0.5,99.5);
  histContainer_["jet_MuonMultiplicity"]= JetVars.make<TH1F>("JetMuonMultiplicity","Muon multiplicity;jet muon multiplicity",30,-0.5,29.5);
  histContainer_["jet_id"]=JetVars.make<TH1F>("jetid","Jet Id level (none, loose, medium, tight);",4,-1.5,2.5);
  
  TString jetIdBinLabels[4] ={"none","loose","medium","tight"};

  for(UInt_t bin=1; bin<=4; bin++){
    histContainer_["jet_id"]->GetXaxis()->SetBinLabel(bin,jetIdBinLabels[bin-1]);    
  }

  // only CALO jets
  // histContainer_["jet_emf"]=JetVars.make<TH1F>("emf","emf;jet EM fraction",40,0.,1.);
  // histContainer_["jet_Area"]= JetVars.make<TH1F>("JetArea","Area;jet Area",40,0.,1.);
  // histContainer_["jet_NCarrying60"]= JetVars.make<TH1F>("JetNCarrying60","n. of constituents with 60% of Jet E; n_{const}(0.6)",40,-0.5,39.5);
  // histContainer_["jet_NCarrying90"]= JetVars.make<TH1F>("JetNCarrying90","n. of constituents with 90% of Jet E; n_{const}(0.9)",40,-0.5,39.5);

  // jet pt corrected

  for(Int_t i = 0; i < corrsize ; i++) {
    h_ptJetsCorrected.push_back(JetVars.make<TH1F>(("jet_pt"+corrections_[i]).c_str(),("p_{T}(Jet);p_{T} (Jet "+corrections_[i]+ ")(GeV)").c_str(),60,0.,300.));
  }
 
  // leading jet variables 
  histContainer_["leadingjet_pt"]=JetVars.make<TH1F>("leadingjet_pt","p_{T}(Jet);p_{T}(Leading Jet) (GeV)",60,0.,300.);
  histContainer_["leadingjet_eta"]=JetVars.make<TH1F>("leadingjet_eta","#eta (Jet);#eta (Leading Jet)",60,-3.,3.);
  histContainer_["leadingjet_phi"]=JetVars.make<TH1F>("leadingjet_phi","#phi (Jet);#phi (Leading Jet)",60,-TMath::Pi(),TMath::Pi());

  // leading jet variables (if any)
  histContainer_["subleadingjet_pt"]=JetVars.make<TH1F>("subleadingjet_pt","p_{T}(Jet);p_{T}(Subleading Jet) (GeV)",60,0.,300.);
  histContainer_["subleadingjet_eta"]=JetVars.make<TH1F>("subleadingjet_eta","#eta (Jet);#eta (Subleading Jet)",60,-3.,3.);
  histContainer_["subleadingjet_phi"]=JetVars.make<TH1F>("subleadingjet_phi","#phi (Jet);#phi (Subleading Jet)",60,-TMath::Pi(),TMath::Pi());
  histContainer_["jet_discrTC"]=JetVars.make<TH1F>("jet_discrTC", "TC discriminant;TC discriminant",100,0.,20.); 
  histContainer_["jet_discrSSVHP"]=JetVars.make<TH1F>("jet_discrSSVHP","SSVHP discriminant; SSVHP discriminant",100,-2.,10.);
  histContainer_["jet_discrSSVHE"]=JetVars.make<TH1F>("jet_discrSSVHE","SSVHE discriminant;SSVHE discriminant",100,-2.,10.);
  histContainer_["jet_discrCSV"]=JetVars.make<TH1F>("jet_discrCSV","CSV discriminant;CSV discriminant",100,0.,1.);
						
  // dijet mass (if available)
  histContainer_["di_jet_mass" ]=JetVars.make<TH1F>("di_jet_mass","M_{jj}; M_{jj} (GeV)",50,0.,500.);

  TFileDirectory MuonVars = fs->mkdir("MuonVars");

  // muon variables (for all muons)
  histContainer_["muon_pt"]=MuonVars.make<TH1F>("muon_pt","muon pt;muon p_{T} (GeV)",150,0.,150.);
  histContainer_["muon_eta"]=MuonVars.make<TH1F>("muon_eta","muon eta; muon #eta",30,-3.,3.);
  histContainer_["muon_phi"]=MuonVars.make<TH1F>("muon_phi","muon phi; muon #phi (rad);",60,-TMath::Pi(),TMath::Pi());
  histContainer_["muon_iso"]=MuonVars.make<TH1F>("muon_iso","muon iso;  (#sum_{trk} p_{T} + #sum_{cal} E_{T}) / p_{T} (combRelIso)",60,0.,20.);
  histContainer_["muon_IP"]=MuonVars.make<TH1F>("muon_IP","muon IP (cm);muon IP (cm)",100,0.,5.);
  histContainer_["muon_SIP"]=MuonVars.make<TH1F>("muon_SIP","muon IP/#sigma_{IP};muon IP/#sigma_{IP}",100,0.,10.);
  histContainer_["isGlobalMuon"]=MuonVars.make<TH1F>("isGlobalMuon","isGlobalMuon",2,-0.5,1.5); 
  histContainer_["isTrackerMuon"]=MuonVars.make<TH1F>("isTrackerMuon","isTrackerMuon",2,-0.5,1.5);  
  histContainer_["muon_chi2"]= MuonVars.make<TH1F>("muon_chi2","muon #chi^{2}/ndof; muon #chi^{2}/ndof", 100,0.,10.);
  histContainer_["muon_trackerhits"]=MuonVars.make<TH1F>("muon_trackerhits","muon Trk hits;tracker hits",40,-0.5,39.5);
  histContainer_["muon_pixelhits"]=MuonVars.make<TH1F>("muon_pixelhits","muon Pxl hits;pixel hits",10,-0.5,9.5);
  histContainer_["muon_muonhits"]=   MuonVars.make<TH1F>("muon_muonhits","muon Muon hits;muon hits",100,-0.,99.5);
  histContainer_["muon_dB"]= MuonVars.make<TH1F>("muon_dB","muon dB; dB [cm]",100,0.,5.);
  histContainer_["muon_numberOfMatches"]= MuonVars.make<TH1F>("muon_numberOfMatches","muon number of matches;n.of matches",10,-0.5,9.5);
  
  histContainer_["AllDiMuMass"] = MuonVars.make<TH1F>("AllDiMuMass","Invariant mass spectrum of all di-muons; M(#mu^{+}#mu^{-}) (GeV)",100,0,120);
  histContainer_["AllDiMuPt"] = MuonVars.make<TH1F>("AllDiMuPt","p_{T} of all di-muons; p_{T}(#mu^{+}#mu^{-}) (GeV)",30,0,150);
  histContainer_["AllDiMuEta"] = MuonVars.make<TH1F>("AllDiMuEta","#eta of all di-muons; #eta(#mu^{+}#mu^{-})",50,-5,5);
  histContainer_["AllDiMuPhi"] = MuonVars.make<TH1F>("AllDiMuPhi","#phi of all di-muons; #phi(#mu^{+}#mu^{-}) (rad)",30,-TMath::Pi(),TMath::Pi());
  histContainer_["AllDiMuDeltaPhi"] = MuonVars.make<TH1F>("AllDiMuDeltaPhi","#Delta #phi of all di-muons; #Delta#phi(#mu^{+}#mu^{-}) (rad)",60,0.,TMath::Pi());  
  histContainer_["AllDiMuDeltaR"] = MuonVars.make<TH1F>("AllDiMuDeltaR","#Delta R of all di-muons; #DeltaR(#mu^{+}#mu^{-})",100,0,10); 

  histContainer_["PreselDiMuMass"] = MuonVars.make<TH1F>("PreselDiMuMass","Invariant mass spectrum of Presel di-muons; M(#mu^{+}#mu^{-}) (GeV)",100,0,120);
  histContainer_["PreselDiMuPt"] = MuonVars.make<TH1F>("PreselDiMuPt","p_{T} of Presel di-muons; p_{T}(#mu^{+}#mu^{-}) (GeV)",30,0,150);
  histContainer_["PreselDiMuEta"] = MuonVars.make<TH1F>("PreselDiMuEta","#eta of Presel di-muons; #eta(#mu^{+}#mu^{-})",50,-5,5);
  histContainer_["PreselDiMuPhi"] = MuonVars.make<TH1F>("PreselDiMuPhi","#phi of Presel di-muons; #phi(#mu^{+}#mu^{-}) (rad)",30,-TMath::Pi(),TMath::Pi()); 
  histContainer_["PreselDiMuDeltaPhi"] = MuonVars.make<TH1F>("PreselDiMuDeltaPhi","#Delta #phi of Presel di-muons; #Delta#phi(#mu^{+}#mu^{-}) (rad)",60,0.,TMath::Pi());  
  histContainer_["PreselDiMuDeltaR"] = MuonVars.make<TH1F>("PreselDiMuDeltaR","#Delta R of Presel di-muons; #DeltaR(#mu^{+}#mu^{-})",100,0,10);  

  histContainer_["MassCutDiMuMass"] = MuonVars.make<TH1F>("MassCutDiMuMass","Invariant mass spectrum of M(ll)> 50 GeV di-muons; M(#mu^{+}#mu^{-}) (GeV)",100,50.,150.);
  histContainer_["MassCutDiMuPt"] = MuonVars.make<TH1F>("MassCutDiMuPt","p_{T} of M(ll)> 50 GeV di-muons; p_{T}(#mu^{+}#mu^{-}) (GeV)",30,0,150);
  histContainer_["MassCutDiMuEta"] = MuonVars.make<TH1F>("MassCutDiMuEta","#eta of M(ll)> 50 GeV di-muons; #eta(#mu^{+}#mu^{-})",50,-5,5);
  histContainer_["MassCutDiMuPhi"] = MuonVars.make<TH1F>("MassCutDiMuPhi","#phi of M(ll)> 50 GeV di-muons; #phi(#mu^{+}#mu^{-}) (rad)",30,-TMath::Pi(),TMath::Pi());  
  histContainer_["MassCutDiMuDeltaPhi"] = MuonVars.make<TH1F>("MassCutDiMuDeltaPhi","#Delta #phi of M(ll)> 50 GeV di-muons; #Delta#phi(#mu^{+}#mu^{-}) (rad)",60,0.,TMath::Pi());  
  histContainer_["MassCutDiMuDeltaR"] = MuonVars.make<TH1F>("MassCutDiMuDeltaR","#Delta R of M(ll)> 50 GeV di-muons; #DeltaR(#mu^{+}#mu^{-})",100,0,10); 
  
  TFileDirectory CandidateVars = fs->mkdir("CandidateVars");

  // Candidate collections histograms
  histContainer_["NZmm"] =CandidateVars.make<TH1F>("NZmm","Number of Z #rightarrow #mu#mu per event; N^{#mu^{+}#mu^{-}}_{cand}",20,-0.5,19.5);

  histContainer_["ZmmMass"]             =CandidateVars.make<TH1F>("ZmmMass","Invariant mass of Z #rightarrow #mu#mu candidates ; M(#mu^{+}#mu^{-}) (GeV)",100,0.,120.);
  histContainer_["ZmmPt"]               =CandidateVars.make<TH1F>("ZmmPt","p_{T} of Z #rightarrow #mu#mu per event; p_{T}(#mu^{+}#mu^{-}) (GeV)",50,0.,150.);
  histContainer_["ZmmEta"]              =CandidateVars.make<TH1F>("ZmmEta","#eta of Z #rightarrow #mu#mu per event; #eta(#mu^{+}#mu^{-})",50,-5,5);
  histContainer_["ZmmPhi"]              =CandidateVars.make<TH1F>("ZmmPhi","#phi of Z #rightarrow #mu#mu per event; #phi(#mu^{+}#mu^{-}) (rad)",30,-TMath::Pi(),TMath::Pi());
  histContainer_["ZmmLeadingLepPt"]     =CandidateVars.make<TH1F>("ZmmLeadingLepPt","p_{T} of leading lepton; p^{lead #mu}_{T} (GeV)",100,0.,100.);
  histContainer_["ZmmLeadingLepEta"]    =CandidateVars.make<TH1F>("ZmmLeadingLepEta","#eta of leading lepton; #eta^{lead #mu}",60,-3,3);
  histContainer_["ZmmLeadingLepPhi"]    =CandidateVars.make<TH1F>("ZmmLeadingLepPhi","#phi of leading lepton; #phi^{lead #mu} (rad)",60,-TMath::Pi(),TMath::Pi());
  histContainer_["ZmmSubLeadingLepPt"]  =CandidateVars.make<TH1F>("ZmmSubLeadingLepPt","p_{T} of sub-leading lepton; p^{2nd #mu}_{T} (GeV)",100,0.,100.);
  histContainer_["ZmmSubLeadingLepEta"] =CandidateVars.make<TH1F>("ZmmSubLeadingLepEta","#eta of sub-leading lepton; #eta^{2nd #mu}",60,-3,3);
  histContainer_["ZmmSubLeadingLepPhi"] =CandidateVars.make<TH1F>("ZmmSubLeadingLepPhi","#phi of sub-leading lepton; #phi^{2nd #mu} (rad)",60,-TMath::Pi(),TMath::Pi());

  histContainer_["NZee"] =CandidateVars.make<TH1F>("NZee","Number of Z #rightarrow ee per event; N^{e^{+}e^{-}}_{cand}",20,-0.5,19.5);

  histContainer_["ZeeMass"]             =CandidateVars.make<TH1F>("ZeeMass","Invariant mass of Z #rightarrow #mu#mu candidates ; M(e^{+}e^{-}) (GeV)",100,0.,120.);
  histContainer_["ZeePt"]               =CandidateVars.make<TH1F>("ZeePt","p_{T} of Z #rightarrow ee per event; p_{T}(e^{+}e^{-}) (GeV)",50,0.,150.);
  histContainer_["ZeeEta"]              =CandidateVars.make<TH1F>("ZeeEta","#eta of Z #rightarrow ee per event; #eta(e^{+}e^{-})",50,-5,5);
  histContainer_["ZeePhi"]              =CandidateVars.make<TH1F>("ZeePhi","#phi of Z #rightarrow ee per event; #phi(e^{+}e^{-}) (rad)",30,-TMath::Pi(),TMath::Pi());
  histContainer_["ZeeLeadingLepPt"]     =CandidateVars.make<TH1F>("ZeeLeadingLepPt","p_{T} of leading lepton; p^{lead e}_{T} (GeV)",100,0.,100.);
  histContainer_["ZeeLeadingLepEta"]    =CandidateVars.make<TH1F>("ZeeLeadingLepEta","#eta of leading lepton; #eta^{lead e}",60,-3,3);
  histContainer_["ZeeLeadingLepPhi"]    =CandidateVars.make<TH1F>("ZeeLeadingLepPhi","#phi of leading lepton; #phi^{lead e} (rad)",60,-TMath::Pi(),TMath::Pi());
  histContainer_["ZeeSubLeadingLepPt"]  =CandidateVars.make<TH1F>("ZeeSubLeadingLepPt","p_{T} of sub-leading lepton; p^{2nd e}_{T} (GeV)",100,0.,100.);
  histContainer_["ZeeSubLeadingLepEta"] =CandidateVars.make<TH1F>("ZeeSubLeadingLepEta","#eta of sub-leading lepton; #eta^{2nd e}",60,-3,3);
  histContainer_["ZeeSubLeadingLepPhi"] =CandidateVars.make<TH1F>("ZeeSubLeadingLepPhi","#phi of sub-leading lepton; #phi^{2nd e} (rad)",60,-TMath::Pi(),TMath::Pi());
   
}

// based on https://hypernews.cern.ch/HyperNews/CMS/get/muon/433/1/1.html
void PATAnalyzer::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup) 
{

  for (size_t i = 0; i < numEventsNames_.size(); i++) {
    
    std::string name = numEventsNames_[i];
    edm::Handle<edm::MergeableCounter> numEventsCounter;
    iLumi.getByLabel(name, numEventsCounter);
    
    if (numEventsCounter.isValid()) 
      histContainer_["SelectedEvtsAfterFilter_"]->AddBinContent(i + 1, numEventsCounter->value); 
  }
}


void 
PATAnalyzer::endJob() 
{
}

Int_t 
PATAnalyzer::jetId(pat::Jet jet){

  Double_t nhf = ( jet.neutralHadronEnergy() + jet.HFHadronEnergy() ) / jet.energy();
  Double_t nEF = jet.neutralEmEnergyFraction();
  Double_t nconstituents = jet.numberOfDaughters();
  Double_t chf = jet.chargedHadronEnergyFraction();
  Double_t nch = jet.chargedMultiplicity();
  Double_t cef = jet.chargedEmEnergyFraction();

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

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATAnalyzer);
