// -*- C++ -*-
//
// Package:    LeadingJetsProducer
// Class:      LeadingJetsProducer
// 
/**\class LeadingJetsProducer LeadingJetsProducer.cc ZbbAnalysis/Tools/plugins/LeadingJetsProducer.cc

*/
//
// Original Authors:  Alberto Traverso, Marco Musich
//          Created:  Wed May 16 21:12:25 CEST 2012
// $Id: LeadingJetsProducer.cc,v 1.1 2012/05/24 22:33:52 musich Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbUtils.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Utilities/interface/EDMException.h"


using namespace edm;
using namespace std;
using namespace reco;
//
// class declaration
//

class LeadingJetsProducer : public edm::EDProducer {
public:
  explicit LeadingJetsProducer(const edm::ParameterSet&);
  ~LeadingJetsProducer() {}
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() {}
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() {}

  float beta(pat::Jet const&, edm::Handle<reco::VertexCollection> const&) const;
  float betaStar(pat::Jet const&, edm::Handle<reco::VertexCollection> const&) const;
  template<class T> std::vector<T> sortObjectsBypT(std::vector<T> const& unsortedObj);

  // ----------member data ---------------------------
  edm::InputTag  jetCollection_,genJetCollection_,primaryVertices_;
  bool   verbose_;
  double etaCut_;
  double betaStarCut_;
  double betaCut_;
  int collectionLength_;
};

LeadingJetsProducer::LeadingJetsProducer(const edm::ParameterSet& iConfig):
  jetCollection_(iConfig.getParameter<edm::InputTag>("jets")),
  genJetCollection_(iConfig.getParameter<edm::InputTag>("genJets")),
  primaryVertices_(iConfig.getParameter<edm::InputTag>("primaryVertices")),
  verbose_(iConfig.getParameter<bool>("verbose")),
  etaCut_(iConfig.getParameter<double>("etaCut")),
  betaStarCut_(iConfig.getParameter<double>("betastarCut")),
  betaCut_(iConfig.getParameter<double>("betaCut")),
  collectionLength_(iConfig.getParameter<int>("collectionLength"))
{
  produces<std::vector<pat::Jet> >("patLeadingJets");
  produces<std::vector<reco::GenJet> >("genLeadingJets");
}

float
LeadingJetsProducer::beta(pat::Jet const& jet, edm::Handle<VertexCollection> const& vertices) const
{
  // definition of beta: ratio of charged pT from first vertex over the sum of all the charged pT in jet. 

  LogDebug("LeadingJetsProducer") << "Computing beta for jet with Pt " << jet.pt()
                              << " in an event with " << vertices->size() << " vertices."; 

  // by definition, 0 if there is no primary vertex.
  if(vertices->size()<1) return 0.;
  // loop over the tracks making the PV, and store the keys
  list<int> trackrefs;
  const reco::Vertex pv = vertices.product()->operator[](0);
  for( reco::Vertex::trackRef_iterator tk = pv.tracks_begin(); tk < pv.tracks_end(); ++tk) {
    trackrefs.push_back(tk->key());
    LogDebug("LeadingJetsProducer") << "Key from PV: " << trackrefs.back();
  }
  // now loop over the jet charged constituents, and count pt
  float ptsum = 0.;
  float ptsumall = 0.;
  const std::vector< reco::PFCandidatePtr > constituents = jet.getPFConstituents();
  for(std::vector< reco::PFCandidatePtr >::const_iterator jetconst = constituents.begin(); jetconst < constituents.end(); ++jetconst) {
    if((*jetconst)->trackRef().isNull()) continue;
    list<int>::iterator found = find(trackrefs.begin(), trackrefs.end(), (*jetconst)->trackRef().key());
    LogDebug("LeadingJetsProducer") << "found charged component in jet with key " << (*jetconst)->trackRef().key() << ". Found = " << (found!=trackrefs.end());
    if(found!=trackrefs.end()) {
      ptsum += (*jetconst)->pt();
    }
    ptsumall += (*jetconst)->pt();
    LogDebug("LeadingJetsProducer") << "ptsum= " << ptsum << " ptsumall= " << ptsumall << " ratio= " << ptsum/ptsumall;
  }
  if(ptsumall<0.001) // non-null: 0.001 is much lower than any pt cut at reco level.
    return -1.;
  return ptsum/ptsumall;
  
}

// ------------ method called to produce the beta*  -----------
float
LeadingJetsProducer::betaStar(pat::Jet const& jet, edm::Handle<VertexCollection> const& vertices) const
{
  // defined as the ratio of charged pT coming from good PU vertices over the sum of all charged pT in jet. 

  LogDebug("LeadingJetsProducer") << "Computing beta* for jet with Pt " << jet.pt()
                              << " in an event with " << vertices->size() << " vertices."; 

  // by definition, 0 if there is no PU vertex.
  if(vertices->size()<2) return 0.;
  // loop over the tracks making the PU vertices, and store the keys
  list<int> trackrefs;
  for(VertexCollection::const_iterator PUvertex = vertices->begin()+1; PUvertex<vertices->end(); ++PUvertex) {
    for( reco::Vertex::trackRef_iterator tk = PUvertex->tracks_begin(); tk < PUvertex->tracks_end(); ++tk) {
      trackrefs.push_back(tk->key());
      LogDebug("LeadingJetsProducer") << "Key from PU: " << trackrefs.back();
    }
  }
  // now loop over the jet charged constituents, and count pt
  float ptsum = 0.;
  float ptsumall = 0.;
  const std::vector< reco::PFCandidatePtr > constituents = jet.getPFConstituents();
  for(std::vector< reco::PFCandidatePtr >::const_iterator jetconst = constituents.begin(); jetconst < constituents.end(); ++jetconst) {
    if((*jetconst)->trackRef().isNull()) continue;
    list<int>::iterator found = find(trackrefs.begin(), trackrefs.end(), (*jetconst)->trackRef().key());
    LogDebug("LeadingJetsProducer") << "found charged component in jet with key " << (*jetconst)->trackRef().key() << ". Found = " << (found!=trackrefs.end());
    if(found!=trackrefs.end()) ptsum += (*jetconst)->pt();
    ptsumall += (*jetconst)->pt();
    LogDebug("LeadingJetsProducer") << "ptsum= " << ptsum << " ptsumall= " << ptsumall << " ratio= " << ptsum/ptsumall;
  }
  if(ptsumall<0.001) // non-null: 0.001 is much lower than any pt cut at reco level.
    return -1.;
  return ptsum/ptsumall;

}

//______________________________________________________________________________
template<class T>
std::vector<T> LeadingJetsProducer::sortObjectsBypT(std::vector<T> const& unsortedObj){

  std::vector<T> sortedObj= unsortedObj;
  std::vector<Float_t> pTObj;
  UInt_t ObjSize=unsortedObj.size();
  for (UInt_t n=0;n <ObjSize; ++n){
    pTObj.push_back(unsortedObj[n].pt());
  }

  for (UInt_t i = 0; i < ObjSize; ++i){
    for (UInt_t j = i+1; j < ObjSize; ++j){
      if (pTObj[i] < pTObj[j]){
	T auxObj = sortedObj[i];
	sortedObj[i] = sortedObj[j];
	sortedObj[j] = auxObj;
      }// if
    }// j loop
  }// i loop

  return sortedObj;

}
		   
// ------------ method called to produce the data  ------------

void
LeadingJetsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 
  
  edm::Handle<std::vector<pat::Jet> > jets; // Jet collection
  iEvent.getByLabel(jetCollection_,jets);   //load it
 
  edm::Handle<reco::GenJetCollection> genJetsHandle;
  iEvent.getByLabel(genJetCollection_,genJetsHandle);
  const reco::GenJetCollection & genJets = *(genJetsHandle.product());

  edm::Handle<VertexCollection> vertices;       // vertex collection
  iEvent.getByLabel(primaryVertices_,vertices); // load it
  
  //                      _     _               _     
  //  _ _ ___ __ ___   _ | |___| |_   __ ___ __| |___ 
  // | '_/ -_) _/ _ \ | || / -_)  _| / _/ _ | _` / -_)
  // |_| \___\__\___/  \__/\___|\__| \__\___|__,_\___|

  std::vector<pat::Jet> selectedJets;
  for(std::vector<pat::Jet>::const_iterator recojet_it=jets->begin(); recojet_it!=jets->end(); ++recojet_it){
    
    if ( !recojet_it->isPFJet()  ){ 
      throw cms::Exception("InvalidInput") 
	<< "Input pat::Jet is not of PF-type !!\n";
    }
    
    // calculate the associations
    float _betastar = betaStar(*recojet_it,vertices);
    float _beta     = beta(*recojet_it,vertices);
    
    if(verbose_)  std::cout<< "jet Pt = " << recojet_it->pt() << " beta= " << _beta << " beta*= " << _betastar<<std::endl;
    
    // starting beta or beta* selection (set beta < 0 default not to apply , b* big not to apply)
    if(fabs(recojet_it->eta())<etaCut_ && _beta > betaCut_ && _betastar < betaStarCut_){
      selectedJets.push_back(*recojet_it);
    }  // if beta cut 
  } // closing for
  
  std::vector<pat::Jet> sortedJetsByPt = sortObjectsBypT(selectedJets);
  std::vector<pat::Jet> leadingJets;
  if(selectedJets.size()){
    int njets(0);
    for(std::vector<pat::Jet>::const_iterator index_it=selectedJets.begin(); index_it <selectedJets.end(); index_it++){
      njets++;
      if(njets<=collectionLength_) leadingJets.push_back(*index_it);
    }
  }
  
  // the output
  std::auto_ptr<std::vector<pat::Jet> > jetColl( new std::vector<pat::Jet> (leadingJets) );

  for (unsigned int i = 0; i< jetColl->size();++i){
    pat::Jet & j = (*jetColl)[i];
    if ( !j.isPFJet() )
      throw cms::Exception("InvalidInput") 
	<< "Input pat::Jet is not of PF-type !!\n";
    j.addUserFloat("beta",beta(j,vertices));
    j.addUserFloat("betaStar",betaStar(j,vertices));
    LogDebug("LeadingJetsProducerSummary") << "jet Pt = " << j.pt() << " beta= " << j.userFloat("beta") << " beta*= " << j.userFloat("betaStar");
  } //closing for 
  
  //                  _     _               _     
  //  __ _ ___ _ _ _ | |___| |_   __ ___ __| |___ 
  // / _` / -_) ' \ || / -_)  _| / _/ _ | _` / -_)
  // \__, \___|_||_\__/\___|\__| \__\___|__,_\___|
  // |___/                                        
  
  std::vector<reco::GenJet> selectedGenJets;
  for(std::vector<reco::GenJet>::const_iterator genjet_it=genJets.begin(); genjet_it!=genJets.end(); ++genjet_it){
    if(fabs(genjet_it->eta())<etaCut_+0.5){
      selectedGenJets.push_back(*genjet_it);
    }  // if eta cut 
  } // ends for on the genJets
  
  std::vector<reco::GenJet> sortedGenJetsByPt = sortObjectsBypT(selectedGenJets);
  std::vector<reco::GenJet> leadingGenJets;
  if(selectedJets.size()){
    int ngenjets(0);
    for(std::vector<reco::GenJet>::const_iterator index_it=selectedGenJets.begin(); index_it <selectedGenJets.end(); index_it++){
      ngenjets++;
      if(ngenjets<=collectionLength_) leadingJets.push_back(*index_it);
    }
  }

  // the output
  std::auto_ptr<std::vector<reco::GenJet> > genjetColl( new std::vector<reco::GenJet> (leadingGenJets) );

  // add to the event
  iEvent.put(jetColl,"patLeadingJets");
  iEvent.put(genjetColl,"genLeadingJets");

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LeadingJetsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LeadingJetsProducer);










