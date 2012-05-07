// -*- C++ -*-
//
// Package:    ZVertexProducer
// Class:      ZVertexProducer
// 
/**\class ZVertexProducer ZVertexProducer.cc ZbbAnalysis/ZVertexProducer/src/ZVertexProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Musich,40 2-A16,+41227671519,
//         Created:  Wed May 25 14:42:59 CEST 2011
// $Id: ZVertexProducer.cc,v 1.4 2012/02/07 14:17:24 musich Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

//
// class declaration
//

class ZVertexProducer : public edm::EDProducer {
   public:
      explicit ZVertexProducer(const edm::ParameterSet&);
      ~ZVertexProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      std::vector<reco::CompositeCandidate> sortCandidatesByDifference(std::vector<reco::CompositeCandidate>);

      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  // ----------member data ---------------------------
  
  edm::InputTag theVertexSrc_;
  edm::InputTag theZmmSrc_;
  edm::InputTag theZeeSrc_;
  std::string alias, iName_; 

};

//
// constructors and destructor
//
ZVertexProducer::ZVertexProducer(const edm::ParameterSet& iConfig):
  theVertexSrc_(iConfig.getUntrackedParameter<edm::InputTag>("VertexSrc")),
  theZmmSrc_(iConfig.getUntrackedParameter<edm::InputTag>("ZmmSrc")),
  theZeeSrc_(iConfig.getUntrackedParameter<edm::InputTag>("ZeeSrc"))
{
  iName_="ZVertexProducer"; 
  produces< std::vector<reco::Vertex> >(alias = iName_ + "ZVertex").setBranchAlias(alias);
  produces<pat::MuonCollection>("patMuonsFromZ");
  produces<pat::ElectronCollection>("patElectronsFromZ");
}


ZVertexProducer::~ZVertexProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
ZVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  
  std::auto_ptr< std::vector<reco::Vertex> > vtxVect ( new std::vector<reco::Vertex> );
  std::auto_ptr< pat::MuonCollection> MuonsFromZ(new pat::MuonCollection());
  std::auto_ptr< pat::ElectronCollection> ElectronsFromZ(new pat::ElectronCollection());


  //  std::cout << "ZVertexProducer::produce() INIT" << std::endl; 

  // get vertex collection
  edm::Handle<edm::View<reco::Vertex> > vertices;
  iEvent.getByLabel(theVertexSrc_,vertices);
  
  // Get the Z->mm collection
  edm::Handle<reco::CompositeCandidateCollection> zmmHandle;
  iEvent.getByLabel(theZmmSrc_, zmmHandle);

  // Get the Z->ee collection
  edm::Handle<reco::CompositeCandidateCollection> zeeHandle;
  iEvent.getByLabel(theZeeSrc_, zeeHandle);

  //  std::cout << "ZVertexProducer::produce() AFTER getByLabel" << std::endl; 

  bool hasZCand_(false);

  reco::CompositeCandidate ZCand;
  if ( zmmHandle.isValid() && !zeeHandle.isValid() ) {
    const reco::CompositeCandidateCollection & zmm = *(zmmHandle.product());   
    if ( zmm.size()>0 ) {
      ZCand=sortCandidatesByDifference(zmm)[0];
      hasZCand_ = true;
    }

  } else if ( !zmmHandle.isValid() && zeeHandle.isValid() ) {
    const reco::CompositeCandidateCollection & zee = *(zeeHandle.product());   
    if ( zee.size()>0 ) {
      ZCand=sortCandidatesByDifference(zee)[0];  
      hasZCand_ = true;
    }

  } else if ( zmmHandle.isValid() && zeeHandle.isValid()) {  
    //    cout << "ZVertexProducer::produce Mixed Flavour" << endl;
    const reco::CompositeCandidateCollection & zmm = *(zmmHandle.product());   
    const reco::CompositeCandidateCollection & zee = *(zeeHandle.product());   

    std::vector<reco::CompositeCandidate> unsortedMixedFlavourCands;

    reco::CompositeCandidate ZmmCand; 
    if ( zmm.size() > 0 ) { 
      hasZCand_ = true;
      ZmmCand =sortCandidatesByDifference(zmm)[0];
      unsortedMixedFlavourCands.push_back(ZmmCand);
    }
    reco::CompositeCandidate ZeeCand; 
    if ( zee.size() > 0 ) { 
      hasZCand_ = true;
      ZeeCand =sortCandidatesByDifference(zee)[0];
      unsortedMixedFlavourCands.push_back(ZeeCand);
    }

    if ( zee.size() > 0 || zmm.size() > 0 ) {  
	   ZCand=sortCandidatesByDifference(unsortedMixedFlavourCands)[0];       
    }
  }

  //  std::cout << "ZVertexProducer::produce() AFTER sortCandidates" << std::endl;   

  if ( hasZCand_ ) {  
    const reco::Candidate* lepton1 = ZCand.daughter(0);
    const reco::Candidate* lepton2 = ZCand.daughter(1);
  
    if(lepton1->isMuon()){
    
      const pat::Muon* muon0 = dynamic_cast<const pat::Muon*>(&(*lepton1->masterClone()));
      MuonsFromZ->push_back(*muon0);
      const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(&(*lepton2->masterClone()));   
      MuonsFromZ->push_back(*muon1);
      
    } else if (lepton1->isElectron()) {
      
      const pat::Electron* ele0 = dynamic_cast<const pat::Electron*>(&(*lepton1->masterClone()));
      ElectronsFromZ->push_back(*ele0);
      const pat::Electron* ele1 = dynamic_cast<const pat::Electron*>(&(*lepton2->masterClone()));
      ElectronsFromZ->push_back(*ele1);
    
    } else {
      std::cout<<"%MSG-e ZVertexProducer::produce: !We have a problem here!"<<std::endl;
      return;
    }

    //    std::cout << "ZVertexProducer::produce() push VTX valid handle" << std::endl;   

    Double_t vz = (lepton1->vz()+lepton2->vz())/2.;
    Double_t mindz = 1000;
    reco::Vertex closestVertex = vertices->at(0);
    unsigned int vertexCollectionSize = vertices->size();  
    for (unsigned int i=0; i<vertexCollectionSize; i++) {
      reco::Vertex vertex = vertices->at(i);
      if( fabs(vertex.z()-vz)<mindz ){
	mindz = fabs(vertex.z()-vz);
	// cout<<"ZVertexProducer::produce vertex found"<<endl;
	closestVertex = vertex;
      }
    }
    vtxVect->push_back(closestVertex);
  } else {

    //    std::cout << "ZVertexProducer::produce() push VTX unvalid handles" << std::endl;   
    unsigned int vertexCollectionSize = vertices->size();  
    for (unsigned int i=0; i<vertexCollectionSize; i++) {
      reco::Vertex vertex = vertices->at(i);
      vtxVect->push_back(vertex);
    }
  }


    
  //  std::cout << "ZVertexProducer::produce() NOW PUT" << std::endl;   
  iEvent.put(vtxVect, iName_+"ZVertex");
  //  std::cout << "ZVertexProducer::produce() PUT VTX OK" << std::endl;   
  iEvent.put(MuonsFromZ,"patMuonsFromZ");
  //  std::cout << "ZVertexProducer::produce() PUT MU OK" << std::endl;   
  iEvent.put(ElectronsFromZ,"patElectronsFromZ");
  //  std::cout << "ZVertexProducer::produce() PUT ELE OK" << std::endl;   

  return;
}

// ------------ method to sort the Z  ------------

std::vector<reco::CompositeCandidate>
ZVertexProducer::sortCandidatesByDifference(std::vector<reco::CompositeCandidate> unsortedCands){

  // if(debug_) std::cout<<"ZbbEventContentAnalyzer::sortCandidatesByDifference"<<std::endl;

  std::vector<reco::CompositeCandidate> sortedCands = unsortedCands;
  Double_t mZ=91.1876;
  std::vector<Double_t> diffZmass; 
  unsigned int ZCandSize=unsortedCands.size();
  
  if(unsortedCands.size()>1){
    // if(debug_) std::cout<<"ZbbEventContentAnalyzer::sortCandidatesByDifference: Z candidates= "<<ZCandSize<<std::endl;

    for(reco::CompositeCandidateCollection::const_iterator ZllCandidate=unsortedCands.begin(); ZllCandidate != unsortedCands.end(); ++ZllCandidate){
      diffZmass.push_back(fabs(mZ - ZllCandidate->mass()));
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
      
    for(reco::CompositeCandidateCollection::const_iterator ZllCandSorted=sortedCands.begin(); ZllCandSorted != sortedCands.end(); ++ZllCandSorted){
      // if(debug_) std::cout<<"Candidate mass:"<<ZllCandSorted->mass()<<std::endl;
    }
  }
  
  return sortedCands;
}

// ------------ method called once each job just before starting event loop  ------------
void 
ZVertexProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZVertexProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
ZVertexProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ZVertexProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ZVertexProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ZVertexProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZVertexProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZVertexProducer);
