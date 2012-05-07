// -*- C++ -*-
//
// Package:    EventListFilter
// Class:      EventListFilter
// 
/**\class EventListFilter EventListFilter.cc ZbbAnalysis/EventListFilter/src/EventListFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Musich,40 2-A16,+41227671519,
//         Created:  Sun May 15 20:08:03 CEST 2011
// $Id: EventListFilter.cc,v 1.1 2011/05/16 13:21:10 musich Exp $
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
#include "FWCore/Utilities/interface/Parse.h"


#include <iostream>
#include <fstream>

//
// class declaration
//

class EventListFilter : public edm::EDFilter {
   public:
      explicit EventListFilter(const edm::ParameterSet&);
      ~EventListFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
      std::string infilename_;
      std::ifstream infile;
      std::vector<unsigned int> the_runnumber, the_luminumber, the_eventnumber;
      double wasRun, wasAccepted;
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
EventListFilter::EventListFilter(const edm::ParameterSet& iConfig)
{
  infilename_= iConfig.getUntrackedParameter<std::string>("Inputfile");
  //now do what ever initialization is needed

}


EventListFilter::~EventListFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
EventListFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  wasRun++;
  using namespace edm;
  using namespace std;
  
  bool thefilter_= false;
  
  for(unsigned int i=0; i<the_runnumber.size(); i++){
    if( iEvent.eventAuxiliary().run() ==  the_runnumber[i] ) { 
      if( iEvent.eventAuxiliary().luminosityBlock() == the_luminumber[i] ){
	if( iEvent.eventAuxiliary().id().event() ==  the_eventnumber[i] )  {
	  cout << "Event matched! Run " <<iEvent.eventAuxiliary().run() << " Event "<<  iEvent.eventAuxiliary().id().event() << " LumiSection "<<  iEvent.eventAuxiliary().luminosityBlock() <<endl;
	  thefilter_=true;
	  wasAccepted++;
	  break;
	} // if event id ok
      } // if lumi ok
    } // if run ok
  } //closes loop on vector
  return thefilter_;
}

// ------------ method called once each job just before starting event loop  ------------
void 
EventListFilter::beginJob()
{

  wasRun=0;
  wasAccepted=0;
  
  using namespace edm;

  infile.open (infilename_.c_str(), std::ifstream::in);

  try {
    if (!infile) 
      throw std::runtime_error("Couldn't open file");
    std::string theline_; 
    while (true) { 
      infile >> theline_;
      std::vector<std::string> tokens = edm::tokenize(theline_, ":");
      the_runnumber.push_back(atoi(tokens[0].c_str())); 
      the_luminumber.push_back(atoi(tokens[1].c_str())); 
      the_eventnumber.push_back(atoi(tokens[2].c_str()));
      if (infile.eof()) break; 
      if (infile.fail())
	throw std::runtime_error("Error while reading from file");
    }
  } catch (std::exception& e) {
    infile.close();
    throw;
  }
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EventListFilter::endJob() {
  infile.close();
  std::cout<<"***********************************"<<std::endl;
  std::cout<<"* EventListFilter::endJob() "<<std::endl;
  std::cout<<"***********************************"<<std::endl;
  std::cout<<"* Was run on       : "<<wasRun<<" events"<<std::endl;
  std::cout<<"* Was accepted on  : "<<wasAccepted<<" events"<<std::endl;
  std::cout<<"* Filter efficiency: "<<(double(wasAccepted/wasRun))*100<<"%"<<std::endl;
  std::cout<<"***********************************"<<std::endl;
}

// ------------ method called when starting to processes a run  ------------
bool 
EventListFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
EventListFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
EventListFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
EventListFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EventListFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(EventListFilter);
