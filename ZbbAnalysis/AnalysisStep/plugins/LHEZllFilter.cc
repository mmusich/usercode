#include <iostream>
#include <string>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"


class LHEZllFilter : public edm::EDFilter {

public:
  /// default constructor
  explicit LHEZllFilter(const edm::ParameterSet&);
  /// default destructor
  ~LHEZllFilter();
  
private:
  std::string zll_decay_;
  int lepton_pdgid_;
  double lepton_etamax_;
  std::vector<double> lepton_ptmin_;

  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  std::vector<math::XYZTLorentzVectorD> sort4MomentaByPt(std::vector<math::XYZTLorentzVectorD> leptonMomenta_);
  virtual void endJob();

  // ---------- member data ---------------------------
  double wasRun;
  double wasAccept;
  std::string theDecayChannel_;
};


LHEZllFilter::LHEZllFilter(const edm::ParameterSet& iConfig)
{
  zll_decay_ = iConfig.getUntrackedParameter<std::string>("ZllDecaySelection");
  if ( zll_decay_=="Zee" ){
    lepton_pdgid_ = 11;
    theDecayChannel_ = "Z->ee";
  } else if ( zll_decay_=="Zmm"){
    lepton_pdgid_ = 13;
    theDecayChannel_ = "Z->mumu";
  }
  else
    std::cout << "UNKNOWN ZLL FINAL STATE! " << std::endl;

  lepton_etamax_  = iConfig.getUntrackedParameter<double>("etaMax");
  lepton_ptmin_   = iConfig.getUntrackedParameter<std::vector<double> >("ptMins");

  // --------- initialize counters --------------------
  wasRun=0;
  wasAccept=0;
} 

LHEZllFilter::~LHEZllFilter()
{
} 

// ------------------------------- filter method ---------------------------------------------------

bool 
LHEZllFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
// based on http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/Validation/EventGenerator/plugins/DuplicationChecker.cc
// http://arxiv.org/pdf/hep-ph/0609017v1
{
  wasRun++;
  bool isFilter_      = false;
  bool isSelectedZll_ = false; 
  bool isInAcceptance_ = false;

  edm::Handle<LHEEventProduct> evt;
  iEvent.getByType( evt );

  const lhef::HEPEUP hepeup_ = evt->hepeup();
  const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP;
  const std::vector<int> idup_ = hepeup_.IDUP;
  const std::vector<int> istup_ = hepeup_.ISTUP;
  
  // loop on the particles in the event
  int nup = hepeup_.NUP;

  // stored vector of lepton momenta
  std::vector<math::XYZTLorentzVectorD> leptonMomenta_; 
  
  for ( int iup=0; iup<nup; iup++) {
    if ( idup_[iup]==lepton_pdgid_ ||  idup_[iup]==-lepton_pdgid_ ) { // is a lepton of the selected flavour
      // std::cout <<"idup_[hepeup_.MOTHUP[iup].second]: "<<idup_[hepeup_.MOTHUP[iup].second] << " idup_[hepeup.MOTHUP[iup].first]: " <<idup_[hepeup_.MOTHUP[iup].first]<<std::endl;
      if(istup_[iup]==1){
	isSelectedZll_ = true;
	// dump the 4-momentum into a LorentzVector
	math::XYZTLorentzVectorD p4((pup_[iup])[0],(pup_[iup])[1],(pup_[iup])[2],(pup_[iup])[3]);  
	leptonMomenta_.push_back(p4);
      } // if status = outgoing
    } // closes if lepton
  } // closes for on LHE ptls
  
  //std::cout<<"lepton vector size: "<<leptonMomenta_.size()<<std::endl;
  
  // sort leptons momenta by pt
  std::vector<math::XYZTLorentzVectorD> sortedLeptonMomenta_ = sort4MomentaByPt(leptonMomenta_); 
  if(sortedLeptonMomenta_.size()==2){
    if( (fabs(sortedLeptonMomenta_[0].eta())<lepton_etamax_ && fabs(sortedLeptonMomenta_[1].eta())<lepton_etamax_) && 
	(sortedLeptonMomenta_[0].pt()>lepton_ptmin_[0] && sortedLeptonMomenta_[1].pt()>lepton_ptmin_[1]) ) isInAcceptance_=true;
  }

  isFilter_= isSelectedZll_ && isInAcceptance_;
  
  if(isFilter_){
    wasAccept++;
  }
  
  return isFilter_;
}

// ------------------------------- beginJob ---------------------------------------------------

void LHEZllFilter::beginJob()
{  
}

// ------------------------------- sort leptons by pt -----------------------------------------

std::vector<math::XYZTLorentzVectorD> 
LHEZllFilter::sort4MomentaByPt(std::vector<math::XYZTLorentzVectorD> leptonMomenta_)
{
  std::vector<math::XYZTLorentzVectorD> sorted4Momenta_ = leptonMomenta_;
  std::vector<double> ptLeptons_;

  for(std::vector<math::XYZTLorentzVectorD>::const_iterator it = leptonMomenta_.begin(); it!=leptonMomenta_.end(); it++){
    ptLeptons_.push_back(it->pt());
  }
  
  for (unsigned int i = 0; i < ptLeptons_.size(); i++) {   
    for (unsigned int j = i+1; j < ptLeptons_.size(); j++) {
      if(ptLeptons_[i] > ptLeptons_[j]) {
	math::XYZTLorentzVectorD auxMomentum = sorted4Momenta_[i];
	sorted4Momenta_[i] = sorted4Momenta_[j];
	sorted4Momenta_[j] = auxMomentum;   
      }// if inverted
    }// loop on second
  }//loop on first 
  
  return sorted4Momenta_;
}

// ------------------------------- endJob -----------------------------------------------------

void LHEZllFilter::endJob()
{
  std::cout<<"***********************************"<<std::endl;
  std::cout<<"* LHEZllFilter::endJob() "<<std::endl;
  std::cout<<"***********************************"<<std::endl;
  std::cout<<"* Was run on       : "<<wasRun<<" events"<<std::endl;
  std::cout<<"* Acceptance defined as |eta(l)|<"<<lepton_etamax_<<" && pt1(l)>"<<lepton_ptmin_[0]<<" && pt2(l)>"<<lepton_ptmin_[1]<<std::endl;
  std::cout<<"* Was accepted on  : "<<wasAccept<<" "<<theDecayChannel_ <<" events"<<std::endl; 
  std::cout<<"* Filter efficiency: "<<(double(wasAccept/wasRun))*100<<"%"<<std::endl;
  std::cout<<"***********************************"<<std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LHEZllFilter);
