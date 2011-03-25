#ifndef PreselAnalyzer_H
#define PreselAnalyzer_H

/** PreselAnalyzer
 */

// Base Class Headers
#include "FWCore/Framework/interface/EDAnalyzer.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

class TH1I;
class TH1F;
class TH2F;

class PreselAnalyzer: public edm::EDAnalyzer {
public:
  /// Constructor
  PreselAnalyzer(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~PreselAnalyzer();

  // Operations

  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);


  virtual void beginJob() ;
  virtual void endJob() ;
protected:

private:

  std::string theMuonLabel_;
  std::string theZllLabel_;
  
  //Spectra
  TH1F *h1dr_;
  TH1F *h1dphi_;
  TH1F *h1deta_;

  //Scatter plots
  TH2F *h2dphideta_;
  TH2F *h2dphipt_;
  TH2F *h2detapt_;
  TH2F *h2drpt_;
  TH2F *h2drmu_;

  // for not selected
  TH1F* h1Pt_NotSel_;
  TH1F* h1Chi2norm_NotSel_;
  TH1F* h1NValTkHits_NotSel_;
  TH1F* h1NValPxHits_NotSel_;
  TH1F* h1NValMuHits_NotSel_;
  TH1F* h1DB_NotSel_;
  TH1F* h1TotIsoOverpT_NotSel_;
  TH1F* h1NofMatch_NotSel_;
  TH1F* h1Eta_NotSel_;

  // Track quality cut
  TH1F* h1Pt_TkQuality_;
  TH1F* h1DB_TkQuality_;
  TH1F* h1Eta_TkQuality_;
  TH1F* h1Phi_TkQuality_;
  TH1F *h1Dxy_TkQuality_;
  TH1F *h1Dz_TkQuality_;

  //Counters
  TH1F *h1PreselFail_;
  TH1F *h1Cut_;

};
#endif

