#ifndef LumiBlockAnalyzer_H
#define LumiBlockAnalyzer_H

/** \class LumiBlockAnalyzer
 *  Analyzer of the lumi block information
 *
 *  $Date: 2011/04/21 08:18:27 $
 *  $Revision: 1.3 $
 *  \author E.Migliore - Torino <migliore@to.infn.it>
 */

// Base Class Headers
#include "FWCore/Framework/interface/EDAnalyzer.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
  class LuminosityBlock;
}

class TH1F;


class LumiBlockAnalyzer: public edm::EDAnalyzer {
public:
  /// Constructor
  LumiBlockAnalyzer(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~LumiBlockAnalyzer();

  // Operations

  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup&);


  virtual void beginJob() ;
  virtual void endJob() ;
protected:

private:

  std::vector<std::pair<uint, uint> > runRangeList_ ;

  // Filter studies
  std::vector<std::string> numEventsNames_;
  TH1F *h1SelectedEvts_;
  int totLS_;
  double totLumiRecorded_, totLumiDelivered_;

  int invalidLS_totLS_;
  double invalidLS_totLumiRecorded_, invalidLS_totLumiDelivered_;

  // Luminosity studies
  std::string theLumiSummaryTag_;
  std::map<int,std::pair<int,double> > reclumibyrun_;

};
#endif

