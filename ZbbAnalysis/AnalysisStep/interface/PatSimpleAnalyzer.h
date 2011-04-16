#ifndef PATSIMPLEANALYZER_H
#define PATSIMPLEANALYZER_H

#include <map>
#include <string>

#include "TH1.h"
#include "TFile.h"
#include <Math/VectorUtil.h>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ServiceRegistry/interface/ServiceRegistry.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

namespace patba {
  class PatSimpleAnalyzer : public edm::EDAnalyzer {
    
  public:
    /// default constructor
    explicit PatSimpleAnalyzer();
    //  explicit PatSimpleAnalyzer(const edm::ParameterSet&);
    /// default destructor
    ~PatSimpleAnalyzer();
  
    //private:
    /// everything that needs to be done before the event loop
    virtual void beginJob() ;
    /// everything that needs to be done during the event loop
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void analyzeEventOnly(const edm::EventBase&);
    /// check if histogram was booked
    bool booked(const std::string histName) const { return histContainer_.find(histName.c_str())!=histContainer_.end(); };
    /// fill histogram if it had been booked before
    void fill(const std::string histName, double value) const { if(booked(histName.c_str())) histContainer_.find(histName.c_str())->second->Fill(value); };
    /// everything that needs to be done after the event loop
    virtual void endJob() ;
    
    // simple map to contain all histograms; 
    // histograms are booked in the beginJob() 
    // method

  private:
    std::map<std::string,TH1F*> histContainer_; 
    std::vector<TH1F*> h_ptJetsCorrected;
    
    // input tags  
    edm::InputTag metSrc_;
    edm::InputTag elecSrc_;
    edm::InputTag muonSrc_;
    edm::InputTag jetSrc_;
    edm::InputTag tauSrc_;
    edm::InputTag photonSrc_;
    
    /// correction levels for pat jet
    std::vector<std::string> str_corrLevels;
    std::vector<std::string> corrections_;
    Int_t corrsize;
    
    // register to the TFileService
    // edm::Service<TFileService> fs;  

  };
}
#endif
