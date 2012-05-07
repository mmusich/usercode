// -*- C++ -*-
//
// Package:    ZTauTauFilter
// Class:      ZTauTauFilter
// 
/**\class ZTauTauFilter ZTauTauFilter.cc ZbbAnalysis/AnalysisStep/plugins/ZTauTauFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Musich , Reviser: Stefano Casasso 40 2-A16,+41227671519,
//         Created:  Thu Apr 14 11:37:27 CEST 2011
// $Id: ZTauTauFilter.cc,v 1.1 2012/03/08 14:49:58 musich Exp $
//

// system include files
#include <memory>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
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
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"

//
// class declaration
//
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

//
// class declaration
//

class ZTauTauFilter : public edm::EDFilter {
   public:
      explicit ZTauTauFilter(const edm::ParameterSet&);
      ~ZTauTauFilter();

   private:
       virtual void beginJob() ;
       virtual bool filter(edm::Event&, const edm::EventSetup&);
       virtual void endJob() ;

  // print event info
  void PrintEvent(const edm::Event& iEvent,Bool_t verbose=false);

  
      // ----------member data ---------------------------
  edm::InputTag genParticles_;
  bool status2_,status3_;
  double ptcut_;
  double wasRun;
  double wasAccept;

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
ZTauTauFilter::ZTauTauFilter(const edm::ParameterSet& iConfig)
{
  genParticles_=iConfig.getParameter<edm::InputTag>("src");
  //now do what ever initialization is needed
}


ZTauTauFilter::~ZTauTauFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ZTauTauFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   wasRun++;
   using namespace edm;
   using namespace std;
   using namespace reco;

   Bool_t theFilter=false;
   Int_t ntausFromZ(0);
   Int_t ntaus(0);

   std::vector<std::pair<int,math::XYZTLorentzVectorD> > idAndPTaus;
   math::XYZTLorentzVectorD pTauPos, pTauNeg, pZTauTau;
   // initialize the gen 4vectors
   pTauPos.SetXYZT(0,0,0,0);
   pTauNeg.SetXYZT(0,0,0,0);
   pZTauTau.SetXYZT(0,0,0,0);

   edm::Handle<GenParticleCollection> genParticlesCollection;
   iEvent.getByLabel(genParticles_, genParticlesCollection);

   for(GenParticleCollection::const_iterator genp = genParticlesCollection->begin();genp != genParticlesCollection->end(); ++ genp ) {  // loop over GEN particles
  
     if(abs(genp->pdgId()==15)){
       ntaus++;
     }

     //------------------------------------------- Z boson -------------------------------------------
     if(genp->pdgId()==23){    
       for(UInt_t i=0; i<genp->numberOfDaughters(); i++){
	 if(abs(genp->daughter(i)->pdgId())==15 && genp->daughter(i)->status()==3){
	   //  std::cout<<"found a tau from Z"<<std::endl;
	   ntausFromZ++;
	   std::pair<int,math::XYZTLorentzVectorD> theTauInfo;
	   theTauInfo = make_pair(0,pTauPos);
	   theTauInfo.first = genp->pdgId();
	   theTauInfo.second = genp->p4();
	   idAndPTaus.push_back(theTauInfo);
	   if(theTauInfo.first > 0) pTauPos = theTauInfo.second;
	   else pTauNeg = theTauInfo.second; 
	 }
       }
     }
   }//end loop over genParticles
   
   if(ntausFromZ==2){
     theFilter=true;
     pZTauTau=pTauNeg+pTauPos;
     // std::cout<<"M(Z->tautau)="<<pZTauTau.M()<<std::endl;
   }

   if(theFilter) {
     // PrintEvent(iEvent);
     wasAccept++;
     return true;
   }
   
   else return false;
   
}

// ------------ method called once each job just before starting event loop  ------------
void 
ZTauTauFilter::beginJob()
{

  wasRun=0;
  wasAccept=0;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZTauTauFilter::endJob() {
  
  std::cout<<"***********************************"<<std::endl;
  std::cout<<"* ZTauTauFilter::endJob() "<<std::endl;
  std::cout<<"***********************************"<<std::endl;
  std::cout<<"* Was run on       : "<<wasRun<<" events"<<std::endl;
  std::cout<<"* Was accepted on  : "<<wasAccept<<" events"<<std::endl;
  std::cout<<"* Filter efficiency: "<<(double(wasAccept/wasRun))*100<<"%"<<std::endl;
  std::cout<<"***********************************"<<std::endl;

}

void 
ZTauTauFilter::PrintEvent(const edm::Event& iEvent, Bool_t verbose){
  if(verbose){
    std::cout<< "================================================================="<<std::endl;
    std::cout<< "Run: "<<iEvent.eventAuxiliary().run()<<"\t  Lumi:"<<iEvent.eventAuxiliary().luminosityBlock()<<"\t Event: "<< iEvent.eventAuxiliary().id().event()<<std::endl;
    time_t rawtime = iEvent.eventAuxiliary().time().unixTime();
    struct tm * timeinfo;
    timeinfo = localtime(&rawtime); 
    std::cout<< "Recorded on:"<< asctime (timeinfo) <<std::endl;
    std::cout<< "-----------------------------------------------------------------"<<std::endl;
  } else {
    std::cout<<iEvent.eventAuxiliary().run()<<":"<<iEvent.eventAuxiliary().luminosityBlock()<<":"<< iEvent.eventAuxiliary().id().event()<<std::endl;
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZTauTauFilter);

