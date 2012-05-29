#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "ZbbAnalysis/AnalysisStep/interface/EventCategory.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbTypeDefs.h"
#include "ZbbAnalysis/AnalysisStep/interface/AcceptanceCuts.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <vector>
#include <TH1D.h>
#include <TString.h>

#ifndef ZbbAnalysis_AnalysisStep_ZbbUtils_h
#define ZbbAnalysis_AnalysisStep_ZbbUtils_h

namespace ZbbUtils {

  // jet flavour constants
  enum FLAVOUR {
    UDSG_JETS = 0,
    C_JETS,
    B_JETS,
    NONID_JETS,
    N_JET_TYPES
  };

  // check if jet id ok
  Bool_t isJetIdOk(const pat::Jet& jet,TString jetcategory);
  // check if jet is ok
  Bool_t isGoodJet(const pat::Jet& jet, const reco::CompositeCandidate& ZCand,const AcceptanceCuts& lCuts);
  // check if jet is btagged
  Bool_t isBJet(const pat::Jet& jet,TString theAlgoWP);
  // check if b-tagged jet is MC matched
  Bool_t isBJetMCMatched(const pat::Jet& jet,edm::Handle<reco::GenParticleCollection> genParticlesCollection,const reco::GenJetCollection & genJets);
  // check if the Z->ll is MC matched
  Bool_t isZLLCandidateMCMatched(const reco::CompositeCandidate& ZCand,edm::Handle<reco::GenParticleCollection> genParticlesCollection);
  // check met in the event
  Bool_t isGoodMet(const pat::MET& met, Double_t theMetCut);
  // helper functions
  Bool_t hasOverlap(Double_t etaJ, Double_t phiJ, Double_t etaLetpt,Double_t phiLept,Double_t theDeltaR);
  // sort the Z candidates by difference with M_Z(PDG)
  std::vector<reco::CompositeCandidate> sortCandidatesByDifference(std::vector<reco::CompositeCandidate> unsortedCands);
  // sort jets by pT
  std::vector<pat::Jet> sortJetsBypT(std::vector<pat::Jet> unsortedJets);
  // sort generic objects by pT
  template<class T> std::vector<T> sortObjectsBypT(std::vector<T> const& unsortedObj);

  //Presel muon cut
  Bool_t isPreselMu(const pat::Muon& muon, const reco::Vertex& vertex,const AcceptanceCuts& lCuts);
  //Presel electron cut
  Bool_t isPreselEle(const pat::Electron& electron, const reco::Vertex& vertex,const AcceptanceCuts& lCuts);
  // check if tight Z
  Bool_t isTightZCandidate(reco::CompositeCandidate ZCand, const reco::BeamSpot& beamSpot, Bool_t isMuChannel,const AcceptanceCuts& lCuts);
  // produce the cut
  Bool_t ProduceTheCut(const EventCategory& ec,TString categoryname,TString extractThisCut);
  // corrects the weight for exclusive final state
  Double_t getTheWeight(const EventCategory& ec,TString categoryname);

  // produce ok deltaphi
  Double_t myDeltaPhi(Double_t phi1, Double_t phi2);
   // Calculation the Collins-Soper angle (adapted from code by R. Arnaldi)
  Double_t CostCS(reco::CompositeCandidate ZCand);
  // Calculation of the cosTheta
  Double_t CosThetaStar(reco::CompositeCandidate ZCand);

  // compute SF for lepton-efficiency
  Double_t getElectronTrgScaleFactor_L( double pt, double eta,const std::string& runPeriod,bool uncert);
  Double_t getElectronTrgScaleFactor_H( double pt, double eta,const std::string& runPeriod,bool uncert);
  Double_t getElectronOfflineScaleFactor( double pt, double eta, bool uncert);  
  
  Double_t getMuonOfflineScaleFactor(double pt, double  eta, bool uncert);
  Double_t getMuonTrgScaleFactor_L( double pt, double eta,const std::string& runPeriod,bool uncert);
  Double_t getMuonTrgScaleFactor_H( double pt, double eta,const std::string& runPeriod,bool uncert);
 
  // compute SF for b-efficiency
  Double_t getbEffScaleFactor(FLAVOUR jetFlavour, const pat::Jet& jet, const edm::EventSetup& iSetup, std::string beffcalibmethod, std::string bmistagcalibmethod);
  Double_t getbEffScaleFactorAR(std::map<EffMCpair,std::vector<double> > theAssociativeMapEffbMC_,std::map<EffMCpair,std::vector<double> > theAssociativeMapEffcMC_,std::vector<pat::Jet> theBtaggedJets, const edm::EventSetup& iSetup,std::string beffcalibmethod, std::string bmistagcalibmethod,bool isExclusive,int minNBtag,int theSFkFactor_,int theLEffkFactor_);

  // compute weight for PU
  std::vector<Double_t> generate_pu_weights(TH1D* data_npu_estimated); 

  // calculates efficiencies 
  std::pair< Double_t, Double_t > effCalc(const double& num, const double& den);

  // get parton flavour
  FLAVOUR getPartonFlavour(const pat::Jet& jet);

  // get the parent code of a GenParticle
  int getParentCode(const reco::GenParticle* genLep);
  
  // function to get the distribution of true expected vertices (for PU reweighting)
  std::vector<float> getPU(Bool_t getData);
}
#endif
