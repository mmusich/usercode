#ifndef JetLeptonCleanAnalyzer_H
#define JetLeptonCleanAnalyzer_H

//class JetLeptonCleanAnalyzer

// Base Class Headers
#include "FWCore/Framework/interface/EDAnalyzer.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
  class LuminosityBlock;
}

class TH1I;
class TH1F;
class TH2F;

class JetLeptonCleanAnalyzer: public edm::EDAnalyzer {
public:
  /// Constructor
  JetLeptonCleanAnalyzer(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~JetLeptonCleanAnalyzer();

  // Operations

  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  
  virtual void beginJob() ;
  virtual void endJob() ;
protected:

private:

  // Muon studies
  std::string theMuonLabel_;
  std::string theZllLabel_;
  std::string theJetsLabel_;
  std::string theJetsPULabel_;
  std::string theJetsJPTLabel_;
  std::string theJetsPFLabel_;

  //Histograms
  TH1F *h1dr_;
  TH2F *h2drpT_;
  
  
  
};
#endif

