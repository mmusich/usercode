/** \class LumiBlockAnalyzer
 *  Analyzer of the lumi block infos
 *
 *  $Date: 2012/03/05 09:55:22 $
 *  $Revision: 1.5 $
 *  \author E.Migliore - Torino <migliore@to.infn.it>
 */


#include "LumiBlockAnalyzer.h"

// Collaborating Class Header
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// for "luminosity"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "TH1F.h"

using namespace std;
using namespace edm;

/// Constructor
LumiBlockAnalyzer::LumiBlockAnalyzer(const ParameterSet& pset) : 
  numEventsNames_(pset.getUntrackedParameter< vector<string> >("numEventsNames")),
  theLumiSummaryTag_(pset.getUntrackedParameter<string>("LumiSummaryTag"))
{
  edm::VParameterSet runrangeselection = pset.getParameter<edm::VParameterSet>("runRangeList");
  for(edm::VParameterSet::const_iterator itpset = runrangeselection.begin(); itpset!=runrangeselection.end(); itpset++){
    std::vector<uint> aRunRange = (*itpset).getUntrackedParameter<std::vector<uint> >("runrange");
    std::pair<uint, uint> aRunRangePair (aRunRange[0],aRunRange[1]);
    runRangeList_.push_back(aRunRangePair);
  }

}

/// Destructor
LumiBlockAnalyzer::~LumiBlockAnalyzer(){
}

void LumiBlockAnalyzer::beginJob(){

  totLS_ = 0;  
  totLumiRecorded_=0; 
  totLumiDelivered_=0;

  invalidLS_totLS_ = 0;
  invalidLS_totLumiRecorded_=0;
  invalidLS_totLumiDelivered_=0;

  reclumibyrun_.clear();
  // Book histograms
  edm::Service<TFileService> fileService;

  unsigned int nbin = numEventsNames_.size();
  if ( nbin > 0 ) {
    h1SelectedEvts_ = fileService->make<TH1F>("SelectedEvts","Selected events",nbin,-0.5,-0.5+nbin);  
    for (unsigned int i = 0;i < nbin; i++) h1SelectedEvts_->GetXaxis()->SetBinLabel(i+1,(numEventsNames_[i].c_str()));
  }

}

void LumiBlockAnalyzer::endJob(){
  cout << "====================================================" << endl;
  cout << "==== LUMI VALUES BELOW SHOULD BE DIVIDED BY 10! ====" << endl;
  cout << "====================================================" << endl;

  double ubarnsTopbarns = 1000000;
  double corrFactor     = 10;
  
  double CF = 1/(ubarnsTopbarns*corrFactor);
  
  map<int,pair<int,double> >::iterator it1;

  for ( it1=reclumibyrun_.begin() ; it1 !=reclumibyrun_.end(); it1++ )
    cout << (*it1).first << " => " << (*it1).second.first << " / "  
	 << setiosflags(ios::fixed) << setprecision(2) 
	 << (*it1).second.second * CF <<"/pb" << endl;
  
  // print out info for valid lumi summaries
  cout << "LumiBlockAnalyzer::endJob() Total number of LS:   "<< totLS_ << endl; 
  cout << "LumiBlockAnalyzer::endJob() Total delivered lumi: "<< setiosflags(ios::fixed) << setprecision(2) 
       << totLumiDelivered_*CF << "/pb "<<endl; 
  cout << "LumiBlockAnalyzer::endJob() Total recorded lumi:  "<< setiosflags(ios::fixed) << setprecision(2) 
       << totLumiRecorded_*CF <<  "/pb"<< endl; 
  

  // print out info for invalid lumi summaries
  cout << "LumiBlockAnalyzer::endJob() Total number of INVALID LS:   "<< invalidLS_totLS_ << endl;
  cout << "LumiBlockAnalyzer::endJob() Total INVALID delivered lumi: "<< setiosflags(ios::fixed) << setprecision(2)
       << invalidLS_totLumiDelivered_*CF << "/pb" << endl;
  cout << "LumiBlockAnalyzer::endJob() Total INVALID recorded lumi:  "<< setiosflags(ios::fixed) << setprecision(2)
       << invalidLS_totLumiRecorded_*CF << "/pb" << endl;


}

 

void LumiBlockAnalyzer::analyze(const Event & event, const EventSetup& eventSetup)
{
}


// based on https://hypernews.cern.ch/HyperNews/CMS/get/muon/433/1/1.html
void LumiBlockAnalyzer::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup) 
{

//  std::string ss="lumiProducer::RECO";
  edm::InputTag theLumiSummaryInputTag_(theLumiSummaryTag_);
  Handle<LumiSummary> lumiSummary;
  iLumi.getByLabel(theLumiSummaryInputTag_, lumiSummary);

  // collect info only for the selected runs
  bool runInTheRange(false);
  std::vector<std::pair<uint, uint> >::iterator it;
  for (it=runRangeList_.begin(); it!=runRangeList_.end(); it++){
    if ( (*it).first<=iLumi.run() && iLumi.run()<=(*it).second ) {
      runInTheRange=true;
      break;
    }
  }
  if (!runInTheRange) return;


  if(lumiSummary->isValid()){    

    float dellumi = lumiSummary->intgDelLumi();
    float reclumi = dellumi*lumiSummary->liveFrac();

    totLS_ ++;
    totLumiDelivered_+=(double)dellumi;
    totLumiRecorded_+=(double)reclumi;     
    reclumibyrun_[iLumi.run()].first++;
    reclumibyrun_[iLumi.run()].second+=(double)reclumi;
  }else{
    cout << "no valid lumi summary data run " << iLumi.run() << endl;

    float dellumi = lumiSummary->intgDelLumi();
    float reclumi = dellumi*lumiSummary->liveFrac();
    invalidLS_totLS_ ++;
    invalidLS_totLumiDelivered_+=(double)dellumi;
    invalidLS_totLumiRecorded_+=(double)reclumi;     
  }  

  for (size_t i = 0; i < numEventsNames_.size(); i++) {
    
    string name = numEventsNames_[i];
    Handle<MergeableCounter> numEventsCounter;
    iLumi.getByLabel(name, numEventsCounter);
    
    if (numEventsCounter.isValid()) 
      if ( h1SelectedEvts_  ) h1SelectedEvts_->AddBinContent(i + 1, numEventsCounter->value); 
  }
  return;
}



DEFINE_FWK_MODULE(LumiBlockAnalyzer);







