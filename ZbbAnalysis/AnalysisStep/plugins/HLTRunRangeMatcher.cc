// -*- C++ -*-
//
// Package:    HLTRunRangeMatcher
// Class:      HLTRunRangeMatcher
// 
/**\class HLTRunRangeMatcher HLTRunRangeMatcher.cc ZbbAnalysis/HLTRunRangeMatcher/src/HLTRunRangeMatcher.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Musich,40 2-A16,+41227671519,
//         Created:  Thu Apr 14 11:37:27 CEST 2011
// $Id: HLTRunRangeMatcher.cc,v 1.3 2011/08/12 20:11:10 musich Exp $
//
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

// Trigger Event includes
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include <DataFormats/Common/interface/RefVector.h>

//
// class declaration
//

class HLTRunRangeMatcher : public edm::EDFilter {
   public:
      explicit HLTRunRangeMatcher(const edm::ParameterSet&);
      ~HLTRunRangeMatcher();

   private:
       virtual void beginJob() ;
       virtual bool filter(edm::Event&, const edm::EventSetup&);
       virtual void endJob() ;
  
      // ----------member data ---------------------------
       std::map<std::vector<std::string>,std::vector<uint> > theAssociativeMap_;
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
HLTRunRangeMatcher::HLTRunRangeMatcher(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

  edm::VParameterSet hltrunrangeselection;
 
  hltrunrangeselection = iConfig.getParameter<edm::VParameterSet>("HLTRunRangeList");
  for(edm::VParameterSet::const_iterator pset = hltrunrangeselection.begin(); pset!=hltrunrangeselection.end(); pset++){
    std::vector<std::string> theString_ = pset->getUntrackedParameter<std::vector<std::string> >("hltpath");
    std::vector<uint> theRunRange_ = pset->getUntrackedParameter<std::vector<uint> >("runrange");
    //std::cout <<"path: "<< theString_ <<" in range:"<< theRunRange_[0] << "-"<< theRunRange_[1] << std::endl;
    theAssociativeMap_[theString_]=theRunRange_; 
  }
  wasRun=0;
  wasAccept=0;
}


HLTRunRangeMatcher::~HLTRunRangeMatcher()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HLTRunRangeMatcher::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   wasRun++;
   using namespace edm;
   using namespace std;

   Bool_t theFilter_=false;

   // Handle to the TriggerEvent
   Handle<pat::TriggerEvent> theTriggerEvent;
   iEvent.getByLabel("patTriggerEvent",theTriggerEvent);
 
   // cout<<"Trigger object match size:  "<<theTriggerEvent->triggerMatchers().size()<<endl; 
  
   // get the fired trigger paths 
   const pat::TriggerPathRefVector acceptedTrigPaths(theTriggerEvent->acceptedPaths());

   // loop over selected trigger objects
   for ( pat::TriggerPathRefVector::const_iterator iTrig = acceptedTrigPaths.begin(); iTrig != acceptedTrigPaths.end(); ++iTrig ) {
     //cout<<"fired path: "<<(*iTrig)->name()<<endl;
     //loop over the map
     map<vector<string>,vector<uint> >::iterator it;  
     for ( it=theAssociativeMap_.begin() ; it !=theAssociativeMap_.end(); it++ ){
       //cout << (*it).first << " => " << (*it).second[0] << "-"<<(*it).second[1]<<endl;
       // if event in the right run-range
       if( iEvent.run () >=(*it).second[0] && iEvent.run () <=(*it).second[1]){
	 //cout<<"the trigger to be matched is: "<<(*it).first[0]<<endl;	 
	 if(((*iTrig)->name())==((*it).first[0])){
	   theFilter_=true;   // if there is a match theFilter_ becomes true and breaks nested for
	   //cout<<"found a matched trigger: "<<(*iTrig)->name()<<endl;
	   break;
	 }
       }// ends if on run-range 
     }// ends loop over map
   }//  ends loop over fired triggerpaths
 
   if(theFilter_) wasAccept++;

   return theFilter_;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HLTRunRangeMatcher::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HLTRunRangeMatcher::endJob() {
  
  std::cout<<"***********************************"<<std::endl;
  std::cout<<"* HLTRunRangeMatcher::endJob() "<<std::endl;
  std::cout<<"***********************************"<<std::endl;
  std::cout<<"* Was run on       : "<<wasRun<<" events"<<std::endl;
  std::cout<<"* Was accepted on  : "<<wasAccept<<" events"<<std::endl;
  std::cout<<"* Filter efficiency: "<<(double(wasAccept/wasRun))*100<<"%"<<std::endl;
  std::cout<<"***********************************"<<std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTRunRangeMatcher);
