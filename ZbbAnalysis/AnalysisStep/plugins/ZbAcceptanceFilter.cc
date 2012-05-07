/////////////////////////////////////////////////////////
// Filter to define the acceptance for Zb analysis     //
//                                                     //
// - lepton acceptance defined based on GEN level info //
//   lepton flavour, eta, pT                           //
//                                                     //
// - b acceptance defined based on b-parton matched    //
//   gen jets matched to reco jets                     //
///////////////////////////////////////////////////////// 

#include <iostream>
#include <string>
#include <algorithm>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbUtils.h"
#include <Math/VectorUtil.h>

class ZbAcceptanceFilter : public edm::EDFilter {

public:
  /// default constructor
  explicit ZbAcceptanceFilter(const edm::ParameterSet&);
  /// default destructor
  ~ZbAcceptanceFilter();
  
private:

  Int_t nGoodBGEN;
  Int_t nGoodBReco;

  // lepton acceptance
  std::string zll_decay_;
  int lepton_pdgid1_;
  int lepton_pdgid2_;
  double muon_etamax_;
  double muon_ptmin_;
  double ele_etamax_;
  double ele_ptmin_;

  // b/c acceptance
  std::string genjet_pdgid_ ;
  int hfparton_pdgid_;
  double genjet_etamax_;
  double genjet_ptmin_;
  bool isExclusiveMeasurement_;
  int ngenjet_good_;  
  edm::InputTag genPSrc_;
  edm::InputTag recoJetSrc_;
  edm::InputTag genJetSrc_;
  edm::InputTag elecSrc_;
  edm::InputTag muonSrc_;

  // methods
  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  std::vector<math::XYZTLorentzVectorD> sort4MomentaByPt(std::vector<math::XYZTLorentzVectorD> leptonMomenta_);
  virtual void endJob();

  // ---------- member data ---------------------------
  double wasRun;
  double wasAccept;
  std::string theDecayChannel_;
};


ZbAcceptanceFilter::ZbAcceptanceFilter(const edm::ParameterSet& iConfig)
{

  // gen particles
  genPSrc_ = iConfig.getUntrackedParameter<edm::InputTag>("genPSrc");

  //===================================================================
  // lepton selection
  zll_decay_ = iConfig.getUntrackedParameter<std::string>("flavLept");
  if ( zll_decay_=="Zee" ){
    lepton_pdgid1_ = 11;
    lepton_pdgid2_ = 11;
    theDecayChannel_ = "Z->ee";
  } else if ( zll_decay_=="Zmm"){
    lepton_pdgid1_ = 13;
    lepton_pdgid2_ = 13;
    theDecayChannel_ = "Z->mumu";
  } else if ( zll_decay_=="all"){
    lepton_pdgid1_ = 11;
    lepton_pdgid2_ = 13;
    theDecayChannel_ = "Z->ll";
  }
  else
    std::cout << "UNKNOWN ZLL FINAL STATE! " << std::endl;

  elecSrc_ = iConfig.getUntrackedParameter<edm::InputTag>("electronSrc");
  muonSrc_ =iConfig.getUntrackedParameter<edm::InputTag>("muonSrc");

  muon_etamax_  = iConfig.getUntrackedParameter<double>("etaMuMax");
  muon_ptmin_   = iConfig.getUntrackedParameter<double>("ptMuMin");
  
  ele_etamax_  = iConfig.getUntrackedParameter<double>("etaEleMax");
  ele_ptmin_   = iConfig.getUntrackedParameter<double>("ptEleMin");
  
  //===================================================================
  // b/c parton/jet selection
  genjet_pdgid_   = iConfig.getUntrackedParameter<std::string>("flavJet");
  if (genjet_pdgid_=="b"){
    hfparton_pdgid_ = 5;
  } else if (genjet_pdgid_=="c"){
    hfparton_pdgid_ = 4;
  } else if (genjet_pdgid_=="all"){
    hfparton_pdgid_ = 0;
  }
  else {
    std::cout << "WARNING: Not a valid FINAL STATE (b/c/all) is specified! " << std::endl;
    hfparton_pdgid_ = -1;  
  }

  genjet_etamax_  = iConfig.getUntrackedParameter<double>("etaGenJetMax");
  genjet_ptmin_   = iConfig.getUntrackedParameter<double>("ptGenJetMin");
  ngenjet_good_   = iConfig.getUntrackedParameter<int>("ngoodGenJet");
  isExclusiveMeasurement_ = iConfig.getUntrackedParameter<bool>("isExclusive");
  recoJetSrc_ = iConfig.getUntrackedParameter<edm::InputTag>("recoJetSrc");
  genJetSrc_ = iConfig.getUntrackedParameter<edm::InputTag>("genJetSrc");

  // --------- initialize counters --------------------
  wasRun=0;
  wasAccept=0;
} 

ZbAcceptanceFilter::~ZbAcceptanceFilter()
{
} 

// ------------------------------- filter method ---------------------------------------------------

bool 
ZbAcceptanceFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 
  wasRun++;
  bool isFilter_      = false;

  //----- lepton section
  bool isSelectedZll_ = false; 
  bool isLLInAcceptance_ = false;
  
  //----- heavy quarks acceptance
  bool isQQInAcceptance_ = false;

  //-----handles
  edm::Handle<reco::GenParticleCollection> genParticlesCollection;
  iEvent.getByLabel(genPSrc_,genParticlesCollection);

  edm::Handle<reco::GenJetCollection> genJetsHandle;
  iEvent.getByLabel(genJetSrc_,genJetsHandle);
  const reco::GenJetCollection & genJets = *(genJetsHandle.product());
  
  edm::Handle<edm::View<pat::Jet> > recoJets;
  iEvent.getByLabel(recoJetSrc_,recoJets);

  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByLabel(elecSrc_,electrons);

  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(muonSrc_,muons);
  
  // stored vector of lepton momenta
  std::vector<math::XYZTLorentzVectorD> leptonMomenta_; 
  
  // stored vector of genJet momenta
  std::vector<math::XYZTLorentzVectorD> selGenJetMomenta_; 

  bool isStatus1_,isStatus2_,isStatus3_,isHF_,hasHdaughter_;
  bool isHFPartonGenJetMatched_;
  bool isGenJetRecoMatched_;
  bool isDiMuonEvent_;
  
  for(reco::GenParticleCollection::const_iterator genp = genParticlesCollection->begin();genp != genParticlesCollection->end(); ++ genp ) {

    //Set booleans to false
    isStatus1_=false;
    isStatus2_=false;
    isStatus3_=false;
    isHF_=false;
    hasHdaughter_=false;
    isHFPartonGenJetMatched_=false;
    isGenJetRecoMatched_=false;
    
    // reset booleans
    if (genp->status()==1)isStatus1_=true;
    if (genp->status()==2)isStatus2_=true;
    if (genp->status()==3)isStatus3_=true;    
    if (fabs(genp->pdgId())==hfparton_pdgid_) isHF_=true;
    
    for (size_t i=0;i<genp->numberOfDaughters();i++){
      if (fabs(genp->daughter(i)->pdgId()) == 91 || fabs(genp->daughter(i)->pdgId()) == 92) hasHdaughter_=true;
    }

    // lepton selection
    if( (fabs(genp->pdgId()) == lepton_pdgid1_ || fabs(genp->pdgId()) == lepton_pdgid2_) && (isStatus3_==true) && (ZbbUtils::getParentCode(&(*genp))==23) ){  // lepton of selected flavour and daughter of a Z
      isSelectedZll_ = true;
      // dump the 4-momentum into a LorentzVector
      math::XYZTLorentzVectorD p4LeptonGEN(genp->px(),genp->py(),genp->pz(),genp->energy());
    
      // gen lepton to reco matching
      if( fabs(genp->pdgId())==13){
	isDiMuonEvent_=true;
	for(edm::View<pat::Muon>::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon){
	  if(ROOT::Math::VectorUtil::DeltaR(muon->momentum(),p4LeptonGEN.Vect())<0.1){
	    leptonMomenta_.push_back(p4LeptonGEN);
	  } // if gen-reco are matched
	} // loop on reco muons
      } else if( fabs(genp->pdgId())==11) {
	for(edm::View<pat::Electron>::const_iterator electron=electrons->begin(); electron!=electrons->end(); ++electron){
	  if(ROOT::Math::VectorUtil::DeltaR(electron->momentum(),p4LeptonGEN.Vect())<0.1){
	    leptonMomenta_.push_back(p4LeptonGEN);
	  } // if gen-reco are matched
	} // loop on reco electrons
      } // flavour switch
    
    } else if( (fabs(genp->pdgId()) == 11 || fabs(genp->pdgId()) == 13) && isStatus3_==true && (ZbbUtils::getParentCode(&(*genp))!=23) ) {
      std::cout<<"Ptl: "<<fabs(genp->pdgId())<<" parent code:"<<ZbbUtils::getParentCode(&(*genp))<<std::endl;
    }
    
    // heavy flavour selection
    if (isHF_){
      if ( (isStatus2_ || isStatus3_) && hasHdaughter_){
	nGoodBGEN++;
	math::XYZTLorentzVectorD p4HFGEN(genp->px(),genp->py(),genp->pz(),genp->energy());
	std::pair<int,double> GenGenAssociation = std::make_pair(-1,9999.);
	double minDeltaRGenGen(9999.);
	int i(0);
	for(std::vector<reco::GenJet>::const_iterator genjet_it=genJets.begin(); genjet_it!=genJets.end(); ++genjet_it){ 
	  if(ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4HFGEN.Vect())<0.5 ) { 
	    isHFPartonGenJetMatched_=true;
	    if(ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4HFGEN.Vect())< minDeltaRGenGen){
	      minDeltaRGenGen = ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4HFGEN.Vect());
	      GenGenAssociation.first  = i;
	      GenGenAssociation.second = minDeltaRGenGen;
	    } // end if minimum distance gen jet-parton 
	  } // ends if parton matched to gen jet
	  //std::cout<<"i: "<<i<<" minDeltaR: "<<minDeltaRGenGen<<std::endl;
	  i++;
	}// ends loop on gen jets

	if(!isHFPartonGenJetMatched_ || GenGenAssociation.first == -1 ) continue;
	
	//std::cout<<"*init********************************************************"<<std::endl;
	//std::cout<<"sel i: "<<GenGenAssociation.first<<" sel minDeltaR: "<<GenGenAssociation.second<<std::endl;
	//std::cout<<"*end********************************************************"<<std::endl;

	for(edm::View<pat::Jet>::const_iterator recojet_it=recoJets->begin(); recojet_it!=recoJets->end(); ++recojet_it){  
	  if(ROOT::Math::VectorUtil::DeltaR(genJets.at(GenGenAssociation.first).momentum(),recojet_it->momentum())<0.5) {
	    isGenJetRecoMatched_=true;
	    nGoodBReco++;
	  } // end if distance gen jet- reco jet < 0.5
	} //ends loop on reco jets
	
	if(isHFPartonGenJetMatched_ && isGenJetRecoMatched_){
	  selGenJetMomenta_.push_back(genJets.at(GenGenAssociation.first).p4());  
	}	
      } // if gen partron "last before hadronization"
    } // is heavy flavour parton
  } // ends loop on genParticles

  // lepton acceptance
  if( leptonMomenta_.size()==2 ) {   
    if(isDiMuonEvent_){
      isLLInAcceptance_ = ( std::max(fabs(leptonMomenta_[0].eta()),fabs(leptonMomenta_[1].eta())) < muon_etamax_ ) && 
	( std::min(leptonMomenta_[0].pt(),leptonMomenta_[1].pt()) > muon_ptmin_ );
    } else {
      isLLInAcceptance_ = ( std::max(fabs(leptonMomenta_[0].eta()),fabs(leptonMomenta_[1].eta())) < ele_etamax_ ) && 
	( std::min(leptonMomenta_[0].pt(),leptonMomenta_[1].pt()) > ele_ptmin_ );
    }
  } 
  
  // heavy flavour acceptance
  int nSelGenJet = 0; 
  for(std::vector<math::XYZTLorentzVectorD>::const_iterator it = selGenJetMomenta_.begin(); it!=selGenJetMomenta_.end(); it++){
    if ( (it->pt()>genjet_ptmin_) && (fabs(it->eta())<genjet_etamax_) ) nSelGenJet++;
  }

  if ( isExclusiveMeasurement_ ) {
    isQQInAcceptance_ =  nSelGenJet == ngenjet_good_ ;
  } else {
    isQQInAcceptance_ =  nSelGenJet >= ngenjet_good_ ;
  }

  // do not throw away events in case no selection on flavour ('all') is done
  if ( hfparton_pdgid_ == 0 ) isQQInAcceptance_ = true;

  isFilter_= isSelectedZll_ && isLLInAcceptance_ && isQQInAcceptance_;
  
  if(isFilter_)
    wasAccept++;
  
  return isFilter_;
}

// ------------------------------- beginJob ---------------------------------------------------

void ZbAcceptanceFilter::beginJob()
{  
 nGoodBGEN=0;
 nGoodBReco=0;
}

// ------------------------------- sort leptons by pt -----------------------------------------

std::vector<math::XYZTLorentzVectorD> 
ZbAcceptanceFilter::sort4MomentaByPt(std::vector<math::XYZTLorentzVectorD> leptonMomenta_)
{
  std::vector<math::XYZTLorentzVectorD> sorted4Momenta_ = leptonMomenta_;
  std::vector<double> ptParticles_;

  for(std::vector<math::XYZTLorentzVectorD>::const_iterator it = leptonMomenta_.begin(); it!=leptonMomenta_.end(); it++){
    ptParticles_.push_back(it->pt());
  }
  
  for (unsigned int i = 0; i < ptParticles_.size(); i++) {   
    for (unsigned int j = i+1; j < ptParticles_.size(); j++) {
      if(ptParticles_[i] > ptParticles_[j]) {
	math::XYZTLorentzVectorD auxMomentum = sorted4Momenta_[i];
	sorted4Momenta_[i] = sorted4Momenta_[j];
	sorted4Momenta_[j] = auxMomentum;   
      }// if inverted
    }// loop on second
  }//loop on first 
  
  return sorted4Momenta_;
}

// ------------------------------- endJob -----------------------------------------------------

void ZbAcceptanceFilter::endJob()
{

  //std::cout<<"nGoodBGEN: "<<nGoodBGEN<<std::endl;
  //std::cout<<"nGoodBReco:"<<nGoodBReco<<std::endl;

  std::cout<<"***********************************"<<std::endl;
  std::cout<<"* ZbAcceptanceFilter::endJob() "<<std::endl;
  std::cout<<"***********************************"<<std::endl;
  std::cout<<"* Was run on       : "<<wasRun<<" events"<<std::endl;
  if(zll_decay_=="Zmm")
    std::cout<<"* Muon acceptance defined as |eta(l)|<"<<muon_etamax_<<" && min(pt1(l),pt2(l)>"<<muon_ptmin_<<std::endl;
  else if(zll_decay_=="Zee")
    std::cout<<"* Electron acceptance defined as |eta(l)|<"<<ele_etamax_<<" && min(pt1(l),pt2(l)>"<<ele_ptmin_<<std::endl;
  else {
    std::cout<<"* Muon acceptance defined as |eta(l)|<"<<muon_etamax_<<" && min(pt1(l),pt2(l)>"<<muon_ptmin_<<std::endl;
    std::cout<<"* Electron acceptance defined as |eta(l)|<"<<ele_etamax_<<" && min(pt1(l),pt2(l)>"<<ele_ptmin_<<std::endl;
  }
  std::cout<<"* HF acceptance defined as |eta(HF)|<"<<genjet_etamax_<<" && pt(HF)>"<<genjet_ptmin_<<std::endl;
  std::cout<<"* Was accepted on  : "<<wasAccept<<" "<<theDecayChannel_ <<" + "<<ngenjet_good_<<" "<<genjet_pdgid_<<" events";
  if(isExclusiveMeasurement_){
    std::cout<<" (exclusive final state)"<<std::endl;
  }else{
    std::cout<<" (inclusive final state)"<<std::endl;
  }
  std::cout<<"* Filter efficiency: "<<(double(wasAccept/wasRun))*100<<"%"<<std::endl;
  std::cout<<"***********************************"<<std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ZbAcceptanceFilter);
