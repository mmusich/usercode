#ifndef ZbbAnalysis_AnalysisStep_ZbbMultiLhistos_h
#define ZbbAnalysis_AnalysisStep_ZbbMultiLhistos_h

#include "ZbbAnalysis/AnalysisStep/interface/ZbbUtils.h"
#include "ZbbAnalysis/AnalysisStep/interface/EventCategory.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

/** \class ZbbMultiLhistos
 *
 * A set of histograms for llbb candidates.
 *
 *  $Date: 2012/02/14 10:50:14 $
 *  $Revision: 1.4 $
 *  \author S. Casasso - Torino 
 */

#include <iostream>
#include <TString.h>
#include <string>
class TH1F;
class TH2F;
class TH2I;

using namespace ZbbUtils;

class ZbbAllLeptons{

  AcceptanceCuts lCuts_;
  TString N_;

  //Muons
  TH1F* h_nmu_;
  TH1F* h_allmu_pt_;
  TH1F* h_allmu_eta_;
  TH1F* h_allmu_trackIso_;
  TH1F* h_allmu_caloIso_;
  TH1F *h_allmu_IPT_;
  TH1F *h_allmu_IPL_;
  TH1F *h_allmu_SIP_;
  TH1F *h_allmu_dr_lj_;
  TH1F *h_allmu_dr_lbj_;
  TH1F *h_allmu_ptrel_lj_;
  TH1F *h_allmu_ptrel_lbj_;

  //Electrons
  TH1F* h_nele_;
  TH1F* h_allele_pt_;
  TH1F* h_allele_eta_;
  TH1F* h_allele_trackIso_;
  TH1F* h_allele_caloIso_;
  TH1F *h_allele_IPT_;
  TH1F *h_allele_IPL_;
  TH1F *h_allele_SIP_;
  TH1F *h_allele_dr_lj_;
  TH1F *h_allele_dr_lbj_;
  TH1F *h_allele_ptrel_lj_;
  TH1F *h_allele_ptrel_lbj_;


public:
  ZbbAllLeptons(const TString& name, AcceptanceCuts theCuts) : 
    lCuts_(theCuts),
    N_(name){}

  void book();
  void fillMu(const EventCategory& ec_mu, const EventCategory& ec_ele, edm::Handle<edm::View<pat::Muon> > muons, edm::Handle<edm::View<pat::Jet> > jets,std::string bTagAlgoWP);
  void fillEle(const EventCategory& ec_mu, const EventCategory& ec_ele,edm::Handle<edm::View<pat::Electron> > electrons, edm::Handle<edm::View<pat::Jet> > jets,std::string bTagAlgoWP);
  ~ZbbAllLeptons(){} 

};

class ZbbExtraLeptons{

  AcceptanceCuts lCuts_;
  TString N_;
 
  //Muons
  TH1F* h_muNotFromZ_pt_;
  TH1F* h_muNotFromZ_eta_;
  TH1F *h_muNotFromZ_phi_;
  TH1F* h_muNotFromZ_IPT_;
  TH1F* h_muNotFromZ_IPL_;
  TH1F* h_muNotFromZ_SIP_;
  TH1F* h_muNotFromZ_dr_lj_;
  TH1F* h_muNotFromZ_dr_lbj_;
  TH1F* h_muNotFromZ_ptrel_lj_;
  TH1F* h_muNotFromZ_ptrel_lbj_;
  TH1F* h_muNotFromZ_trackIso_;
  TH1F* h_muNotFromZ_caloIso_;
  TH1F* h_muNotFromZ_combIso_;
  TH1F* h_addDiMu_mass_;

  //Electrons
  TH1F* h_eleNotFromZ_pt_;
  TH1F* h_eleNotFromZ_eta_;
  TH1F *h_eleNotFromZ_phi_;
  TH1F* h_eleNotFromZ_IPT_;
  TH1F* h_eleNotFromZ_IPL_;
  TH1F* h_eleNotFromZ_SIP_;
  TH1F* h_eleNotFromZ_dr_lj_;
  TH1F* h_eleNotFromZ_dr_lbj_;
  TH1F* h_eleNotFromZ_ptrel_lj_;
  TH1F* h_eleNotFromZ_ptrel_lbj_;
  TH1F* h_eleNotFromZ_trackIso_;
  TH1F* h_eleNotFromZ_caloIso_;
  TH1F* h_eleNotFromZ_combIso_;
  TH1F* h_addDiEle_mass_;

  //Dileptons
  TH1F *h_addDiLept_mass_;

  //Counters
  TH1F* h_addmu_;
  TH1F* h_addele_;
  TH1F* h_addlept_;
  TH1F* h_addleptflav_3levts_;
  TH1F* h_addleptflav_4levts_;

 public:
  ZbbExtraLeptons(const TString& name, AcceptanceCuts theCuts) :
    lCuts_(theCuts),
    N_(name){}
    
  void book();
  void fillLeptons(const EventCategory& ec_mu, const EventCategory& ec_ele, edm::Handle<edm::View<pat::Muon> > muons, edm::Handle<edm::View<pat::Electron> > electrons, edm::Handle<edm::View<pat::Jet> > jets,std::string bTagAlgoWP);
  ~ZbbExtraLeptons(){} 
};

class MCTruthLFromHF{

  AcceptanceCuts lCuts_;
  TString N_;
 
  //Generic lepton histograms
  TH1F *h_pt_;
  TH1F *h_eta_;
  TH1F *h_phi_;
  TH1F *h_dr_lj_;
  TH1F *h_dr_lbj_;
  TH1F *h_ptrel_lj_;
  TH1F *h_ptrel_lbj_;
  TH1F *h_IPT_;
  TH1F *h_IPL_;
  TH1F *h_CombIso_;

  //Specific muon histograms
  TH1F *h_Mu_RecoAlgo_;
  TH1F *h_Mu_MuonSelector_;
  TH1F *h_Mu_NormChi2_;
  TH1F *h_Mu_TrkHits_;
  TH1F *h_Mu_PixHits_;
  TH1F *h_Mu_MuHits_;
  TH1F *h_Mu_NOfMatches_;
  TH1F *h_Mu_parents_;

  //Specific electron histograms
  TH1F *h_Ele_eleID_;
  TH1F *h_Ele_parents_;

  //GEN level histograms
  TH1F *h_pt_GENmatchedRECO_;
  TH1F *h_pt_rebin_GENmatchedRECO_;
  TH1F *h_eta_GENmatchedRECO_;
  TH1F *h_phi_GENmatchedRECO_;
  TH1F *h_pt_GEN_;
  TH1F *h_pt_rebin_GEN_;
  TH1F *h_eta_GEN_;
  TH1F *h_phi_GEN_;
  
 public:

  MCTruthLFromHF(const TString& name, AcceptanceCuts theCuts) : 
    lCuts_(theCuts),
    N_(name){}
  void book();
  void fillMu(const EventCategory& ec_mu,const EventCategory& ec_ele,const pat::Muon &muon,edm::Handle<edm::View<pat::Jet> > jets,std::string bTagAlgoWP);
  void fillEle(const EventCategory& ec_mu,const EventCategory& ec_ele,const pat::Electron &electron,edm::Handle<edm::View<pat::Jet> > jets,std::string bTagAlgoWP);
  void fillGenParticle(const reco::Candidate &genp);
  ~MCTruthLFromHF(){}    
};


#endif


