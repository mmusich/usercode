#ifndef ZbbAnalysis_AnalysisStep_ZbbHistos_h
#define ZbbAnalysis_AnalysisStep_ZbbHistos_h

#include "ZbbAnalysis/AnalysisStep/interface/ZbbUtils.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbTypeDefs.h"
#include "ZbbAnalysis/AnalysisStep/interface/Zbbstruct4JEC.h"
#include "ZbbAnalysis/AnalysisStep/interface/EventCategory.h"
#include "ZbbAnalysis/AnalysisStep/interface/AcceptanceCuts.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

/** \class ZbbHistos
 *
 * A set of histograms for llbb candidates.
 *
 *  $Date: 2012/02/16 10:21:11 $
 *  $Revision: 1.33 $
 *  \author M. Musich - Torino
 */

#include <iostream>
#include <TString.h>
#include <string>
#include "TTree.h"
#include "TNtuple.h"
class TH1F;
class TH2F;
class TH2I;
class TProfile; 
class XYZTLorentzVectorD;

using namespace ZbbUtils;

//______________________________________________________________________________
class ZbbCandidate{

  AcceptanceCuts lCuts_;
  TString N_;

  // event weight
  TH1F* h_eventWeight_;

  // z plots
  TH1F* h_zllmass_; 
  TH1F* h_zllpt_; 
  TH1F* h_zlleta_; 
  TH1F* h_lldeltaPhi_;
  TH1F* h_lldeltaR_;  
  TH1F* h_zllCosThetaCS_;
  TH1F* h_zllCosThetaStar_;
  TH1F* h_zllAcop_;
  TH1F* h_zllMt_;   
    
  // muons
  TH1F* h_zmumumass_; 
  TH1F* h_zmumupt_; 
  TH1F* h_zmumueta_; 
  TH1F* h_mu1pt_;     
  TH1F* h_mu2pt_;     
  TH1F* h_mu1eta_;    
  TH1F* h_mu2eta_;    
  TH1F* h_mumudeltaPhi_;
  TH1F* h_mumudeltaR_;  
  TH1F* h_zmumuCosThetaCS_;
  TH1F* h_zmumuCosThetaStar_;
  TH1F* h_zmumuAcop_;
  TH1F* h_zmumuMt_;

  // electrons
  TH1F* h_zeemass_; 
  TH1F* h_zeept_;   
  TH1F* h_zeeeta_; 
  TH1F* h_ele1pt_;     
  TH1F* h_ele2pt_;     
  TH1F* h_ele1eta_;    
  TH1F* h_ele2eta_;    
  TH1F* h_eedeltaPhi_; 
  TH1F* h_eedeltaR_;  
  TH1F* h_zeeCosThetaCS_;
  TH1F* h_zeeCosThetaStar_;
  TH1F* h_zeeAcop_;
  TH1F* h_zeeMt_;

  //jet specific
  TH1F* h_nj_;            
  TH1F* h_SSVHEdisc_;     
  TH1F* h_SSVHPdisc_;    
  TH1F* h_TCHEdisc_;   
  TH1F* h_TCHPdisc_;
  TH1F* h_CSVdisc_;
  TH1F* h_JPdisc_;
  TH1F* h_met_;           
  TH1F* h_phimet_;        
  TH1F* h_jetpt_;         
  TH1F* h_jeteta_;        
  TH1F* h_jetphi_;        
  TH1F* h_jetoverlapLept_;  
  TH1F* h_jet1pt_;        
  TH1F* h_jet1eta_;       
  TH1F* h_jet2pt_;        
  TH1F* h_jet2eta_;       
  TH1F* h_nhf_;           
  TH1F* h_nef_;           
  TH1F* h_nconstituents_; 
  TH1F* h_chf_;           
  TH1F* h_nch_;           
  TH1F* h_cef_;           
  TH1F* h_jetid_;         
  TH1F* h_alljetDeltaR_;
  TH1F* h_alljetDeltaPhi_;
  
  //b jet stuff   
  TH1F* h_nb_;          
  TH2I* h_njb_;     
  TH1F* h_sumbHtOversumHt_; 
  TH1F* h_bjetpt_;     
  TH1F* h_bjeteta_; 
  TH1F* h_bjetmass_;

  //sv stuff
  TH1F* h_SVnVertices_; 
  TH1F* h_SVJetdeltaR_; 
  TH1F* h_SVmass_;
  TH1F* h_SVdist_;      
  TH1F* h_SVdistErr_;   
  TH1F* h_SVdistSig_;   
  TH1F* h_SVnTracks_;   
  TH1F* h_SVchi2_;      

  TH1F* h_bjet1pt_;     
  TH1F* h_bjet1eta_; 
  TH1F* h_bjet1mass_;
  
  TH1F* h_SSVHEdisc1_;
  TH1F* h_SSVHPdisc1_;
  TH1F* h_TCHEdisc1_; 
  TH1F* h_TCHPdisc1_;          		 
  TH1F* h_CSVdisc1_;
  TH1F* h_JPdisc1_;

  //TProfile* p_bTagAlgoWP_SFb_pt_;
  //TProfile* p_bTagAlgoWP_EffbMC_pt_;
  TH1F*  h_bTagAlgoWP_MCflav_;
    
  TH1F* h_bjet2pt_;     
  TH1F* h_bjet2eta_;  
  TH1F* h_bjet2mass_;

  TH1F* h_SSVHEdisc2_;
  TH1F* h_SSVHPdisc2_;
  TH1F* h_TCHEdisc2_; 
  TH1F* h_TCHPdisc2_; 
  TH1F* h_CSVdisc2_;
  TH1F* h_JPdisc2_;

  TH1F* h_bjetOtherJetsDeltaR_;
  TH1F* h_bjetOtherJetsDeltaPhi_;

  TH1F* h_scaldptZbj1_; 
  TH1F* h_drZbj1_;      
  TH1F* h_vecdptZbj1_;  
  TH1F* h_dphiZbj1_;    
  TH1F* h_bbM_;      
  TH1F* h_bbPt_;
  TH1F* h_bbDeltaR_;
  TH1F* h_ZbM_;         
  TH1F* h_ZbPt_;        
  TH1F* h_ZbbM_;        
  TH1F* h_ZbbPt_;       
  TH2F* h_ZbbM2D_;      

public:

  ZbbCandidate(const TString& name, AcceptanceCuts theCuts): 
    lCuts_(theCuts),
    N_(name){}
  void book();
  void fillMu(const EventCategory& ec);
  void fillEle(const EventCategory& ec);
  void fillJets(const EventCategory& ec_mu, const EventCategory& ec_ele, edm::Handle<edm::View<pat::Jet> > jets, edm::Handle<edm::View<pat::MET> > mets, std::string bTagAlgoWP, 
  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffbMC,  
  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffcMC);
  ~ZbbCandidate(){}    
};


//______________________________________________________________________________
class ZbbCandidateNtuple{

  AcceptanceCuts lCuts_;
  TString N_;

  //Structures and TTree
  Event_info Event; Z_info Z; bjet_info bjet; jet2_info jet2; MET_info MET; 
  btagjets_info btagjets; std::vector<bjet_info> btagjets2;
  TTree* ZbbNtuple;

 public:

 ZbbCandidateNtuple(const TString& name, AcceptanceCuts theCuts):
  lCuts_(theCuts),
  N_(name){}
  void book();
  void fill(const EventCategory& ec_mu, const EventCategory& ec_ele, edm::Handle<edm::View<pat::Jet> > jets, std::string bTagAlgoWP, edm::Handle<edm::View<pat::MET> > mets, const edm::Event& iEvent, Int_t nvvertex); //,edm::Handle<reco::GenParticleCollection> genParticlesCollection, const reco::GenJetCollection& genJets);
  ~ZbbCandidateNtuple(){}    
  
};


//______________________________________________________________________________
class ZbbMatchedCandidate{

  AcceptanceCuts lCuts_;
  TString N_;

  //--------------- all leptons
  // z plots
  TH1F* h_zllmass_; 
  TH1F* h_zllpt_; 
  TH1F* h_zllrapidity_;  
  TH1F* h_zllCosThetaStar_;  
  TH1F* h_lept1pt_;     
  TH1F* h_lept2pt_;     
  TH1F* h_lept1eta_;    
  TH1F* h_lept2eta_;    
  
  //b jet stuff   
  TH1F* h_nb_; 
  TH1F* h_sumbHtOversumHt_;
  TH1F* h_bjetmass_;
  TH1F* h_bjet1pt_;     
  TH1F* h_bjet1eta_; 
  TH1F* h_bjet2pt_;     
  TH1F* h_bjet2eta_;  
  TH1F* h_bbDeltaR_;
  TH1F* h_bbM_;

  //--------------- muons
  // z plots
  TH1F* h_zmumumass_; 
  TH1F* h_zmumupt_; 
  TH1F* h_zmumurapidity_;  
  TH1F* h_zmumuCosThetaStar_;  
  TH1F* h_mu1pt_;     
  TH1F* h_mu2pt_;     
  TH1F* h_mu1eta_;    
  TH1F* h_mu2eta_;    
  
  //b jet stuff   
  TH1F* h_mmnb_; 
  TH1F* h_mmsumbHtOversumHt_;
  TH1F* h_mmbjetmass_;
  TH1F* h_mmbjet1pt_;     
  TH1F* h_mmbjet1eta_; 
  TH1F* h_mmbjet2pt_;     
  TH1F* h_mmbjet2eta_;  
  TH1F* h_mmbbDeltaR_;
  TH1F* h_mmbbM_;

  //--------------- electrons
  // z plots
  TH1F* h_zeleelemass_; 
  TH1F* h_zeleelept_; 
  TH1F* h_zeleelerapidity_;  
  TH1F* h_zeleeleCosThetaStar_;  
  TH1F* h_ele1pt_;     
  TH1F* h_ele2pt_;     
  TH1F* h_ele1eta_;    
  TH1F* h_ele2eta_;    
  
  //b jet stuff   
  TH1F* h_eenb_; 
  TH1F* h_eesumbHtOversumHt_;
  TH1F* h_eebjetmass_;
  TH1F* h_eebjet1pt_;     
  TH1F* h_eebjet1eta_; 
  TH1F* h_eebjet2pt_;     
  TH1F* h_eebjet2eta_;  
  TH1F* h_eebbDeltaR_;
  TH1F* h_eebbM_;

public:

  ZbbMatchedCandidate(const TString& name, AcceptanceCuts theCuts): 
    lCuts_(theCuts),
    N_(name){}
  void book();
  void fillMC(const EventCategory& ec_mu, const EventCategory& ec_ele, edm::Handle<edm::View<pat::Jet> > jets, std::string bTagAlgoWP,
	      edm::Handle<reco::GenParticleCollection> genParticlesCollection,const reco::GenJetCollection & genJets);
  void fillDATA(const EventCategory& ec_mu, const EventCategory& ec_ele, edm::Handle<edm::View<pat::Jet> > jets, std::string bTagAlgoWP);
  
  ~ZbbMatchedCandidate(){}    
};

//______________________________________________________________________________
class ZbbCollections{

  AcceptanceCuts lCuts_;
  TString N_;

  TH1F* h_eventWeight_;
  
  // z plots
  TH1F* h_zllmass_; 
  TH1F* h_zllpt_; 
  TH1F* h_zlleta_; 
  TH1F* h_lldeltaPhi_;
  TH1F* h_lldeltaR_;  
  TH1F* h_zllCosThetaCS_;
  TH1F* h_zllCosThetaStar_;

   // Z -> muons
  TH1F* h_zmumumass_; 
  TH1F* h_zmumupt_;   
  TH1F* h_zmumueta_; 
  TH1F* h_mu1pt_;     
  TH1F* h_mu2pt_;     
  TH1F* h_mu1eta_;    
  TH1F* h_mu2eta_;   
  TH1F* h_mumudeltaPhi_;
  TH1F* h_mumudeltaR_;  
  TH1F* h_zmumuCosThetaCS_;
  TH1F* h_zmumuCosThetaStar_;
  
  // Z -> electrons
  TH1F* h_zeemass_; 
  TH1F* h_zeept_;   
  TH1F* h_zeeeta_; 
  TH1F* h_ele1pt_;     
  TH1F* h_ele2pt_;     
  TH1F* h_ele1eta_;    
  TH1F* h_ele2eta_;    
  TH1F* h_eedeltaPhi_; 
  TH1F* h_eedeltaR_; 
  TH1F* h_zeeCosThetaCS_;
  TH1F* h_zeeCosThetaStar_;

  //jet specific
  TH1F* h_nj_;            
  TH1F* h_SSVHEdisc_;     
  TH1F* h_SSVHPdisc_;   
  TH1F* h_TCHEdisc_;   
  TH1F* h_TCHPdisc_;
  TH1F* h_CSVdisc_;
  TH1F* h_JPdisc_;
  TH1F* h_met_;           
  TH1F* h_phimet_;        
  TH1F* h_jetpt_;         
  TH1F* h_jeteta_;        
  TH1F* h_jetphi_;        
  TH1F* h_jetoverlapLept_;  
  TH1F* h_jet1pt_;        
  TH1F* h_jet1eta_;       
  TH1F* h_jet2pt_;        
  TH1F* h_jet2eta_;       
  TH1F* h_nhf_;           
  TH1F* h_nef_;           
  TH1F* h_nconstituents_; 
  TH1F* h_chf_;           
  TH1F* h_nch_;           
  TH1F* h_cef_;           
  TH1F* h_jetid_;  
  TH1F* h_alljetDeltaR_;
  TH1F* h_alljetDeltaPhi_;
  
  TH1F* h_rawjetspt_;  
  TH1F* h_rawjetseta_; 
  TH1F* h_jetMCflav_;

public:
  ZbbCollections(const TString& name, AcceptanceCuts theCuts) : 
    lCuts_(theCuts),
    N_(name){}
  
  void book();
  void fillMu(const EventCategory& ec);
  void fillEle(const EventCategory& ec);
  void fillJets(const EventCategory& ec_mu, const EventCategory& ec_ele, edm::Handle<edm::View<pat::Jet> > jets, edm::Handle<edm::View<pat::MET> > mets);
  void fillLeptons(const EventCategory& ec_mu,const EventCategory& ec_ele, std::string bTagAlgoWP, edm::Handle<edm::View<pat::Muon> > muons, edm::Handle<edm::View<pat::Electron> > electrons, edm::Handle<edm::View<pat::Jet> > jets);

  ~ZbbCollections(){} 
};

//______________________________________________________________________________
class ZbbMCinfo{

  AcceptanceCuts lCuts_;
  TString N_;


  // from genParticles  
  TH1F* h_nbGENP_;
  TH1F* h_ncGENP_;
  
  // Z histos

  TH1F* h_GENP_pt_Z_;  
  TH1F* h_GENP_eta_Z_;
  TH1F* h_GENP_mass_Z_;
  		
  TH1F* h_GENP_pt_ZfromLept_;   
  TH1F* h_GENP_eta_ZfromLept_;  
  TH1F* h_GENP_mass_ZfromLept_; 
  
  TH1F* h_GENP_pt_ZfromElectrons_;   
  TH1F* h_GENP_eta_ZfromElectrons_;  
  TH1F* h_GENP_mass_ZfromElectrons_; 

  TH1F* h_GENP_pt_ZfromMuons_;   
  TH1F* h_GENP_eta_ZfromMuons_;  
  TH1F* h_GENP_mass_ZfromMuons_; 

  // leading b histos
  TH1F* h_GENP_pt_b1_;
  TH1F* h_GENP_eta_b1_;
  TH1F* h_GENP_mass_b1_;

  TH1F* h_GENP_recoMatched_pt_b1_;
  TH1F* h_GENP_recoMatched_eta_b1_;
  
  TH1F* h_GENP_recoBtagged_pt_b1_;
  TH1F* h_GENP_recoBtagged_eta_b1_;

  // from LHEvent
  TH1F* h_nbLHE_;

  //for matching studies
  TH1F* h_GENPDeltaR_b1LBH_RECOJet_;

  // special histos for b/c efficiency studies
  TH1F* h_GENP_pt_b1_etaBrl_;
  TH1F* h_GENP_pt_b1_etaFwd_;
  TH1F* h_GENP_pt_c1_etaBrl_;
  TH1F* h_GENP_pt_c1_etaFwd_;

  TH1F* h_RECO_pt_b1_etaBrl_;
  TH1F* h_RECO_pt_b1_etaFwd_;
  TH1F* h_RECO_pt_c1_etaBrl_;
  TH1F* h_RECO_pt_c1_etaFwd_;

  // Z histos
  TH1F* h_LHE_pt_Z_;   
  TH1F* h_LHE_eta_Z_;  
  TH1F* h_LHE_mass_Z_; 

  // leading b histos
  TH1F* h_LHE_pt_b1_;
  TH1F* h_LHE_eta_b1_;
  TH1F* h_LHE_mass_b1_;

  TH1F* h_LHE_recoMatched_pt_b1_;  
  TH1F* h_LHE_recoMatched_eta_b1_; 

  TH1F* h_LHE_recoBtagged_pt_b1_;  
  TH1F* h_LHE_recoBtagged_eta_b1_; 
  
  // correlations
  
  TH2F* h2_LHEpt_vs_GENpt;
  TH2F* h2_LHEpt_vs_RECOJetpt;
  TH2F* h2_GENpt_vs_RECOJetpt;

public:
  ZbbMCinfo(const TString& name, AcceptanceCuts theCuts) : 
    lCuts_(theCuts),
    N_(name){}

  void book();
  void fill(const EventCategory& ec_mu,const EventCategory& ec_ele,std::string bTagAlgoWP, edm::Handle<LHEEventProduct> lh_evt, edm::Handle<reco::GenParticleCollection> genParticlesCollection, const reco::GenJetCollection& genJets, edm::Handle<edm::View<pat::Jet> > recoJets);
  ~ZbbMCinfo(){}
};

//______________________________________________________________________________
class ZbbBasicComponents{

  AcceptanceCuts lCuts_;
  TString N_;

  // electron specific
  TH1F* h_electron_SIP_;  
  TH1F* h_electron_iso_;     
  TH1F* h_electron_eop_;     
  TH1F* h_electron_clus_;    
  TH1F* h_electron_eIDs_;  
  TH1F* h_electron_eIDcut_; 
  TH1F* h_electronpair_type_;
  TH1F* h_electron_type_;

  // muon specific
  TH1F* h_muon_iso_;         
  TH1F* h_muon_SIP_;         
  TH1F* h_muon_chi2_;        
  TH1F* h_muon_trackerhits_; 
  TH1F* h_muon_pixelhits_;   
  TH1F* h_muon_numberOfMatches_;
  TH1F* h_muon_muonhits_;  

  TH1F* h_muon_muonhitsCSCOnly_;
  TH1F* h_muon_muonhitsEndcaps_;
  TH1F* h_muon_muonhitsOverlap_; 
  TH1F* h_muon_muonhitsBarrel_;  

  TH1F* h_muonpair_type_;
  TH1F* h_muon_type_;

public:
  ZbbBasicComponents(const TString& name, AcceptanceCuts theCuts) : 
    lCuts_(theCuts),
    N_(name){}

  void book();
  void fillMu(const EventCategory& ec);
  void fillEle(const EventCategory& ec);
  ~ZbbBasicComponents(){}
};

//______________________________________________________________________________
class ZbbABCDMatrix{
  TString N_;
 
  TH2F* h_MllVsMET_; 
  TH1F* h_MllVsNOMASSCUT_;
  
  TH2F* h_MmumuVsMET_; 
  TH1F* h_MmumuVsNOMASSCUT_;

  TH2F* h_MeeVsMET_; 
  TH1F* h_MeeVsNOMASSCUT_;
  
public:
  ZbbABCDMatrix(const TString& name) : 
    N_(name){}
    
  void book();
  void fill(const EventCategory& ec_mu, const EventCategory& ec_ele,edm::Handle<edm::View<pat::MET> > mets);
  void fillMu(const EventCategory& ec_mu,edm::Handle<edm::View<pat::MET> > mets);
  void fillEle(const EventCategory& ec_ele,edm::Handle<edm::View<pat::MET> > mets);
  ~ZbbABCDMatrix(){}  
};

#endif



