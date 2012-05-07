// -*- C++ -*-
//
// Package:    AnalysisStep
// Class:      ZbbAcceptanceCalculator
// 
/**\class ZbbAcceptanceCalculator ZbbAcceptanceCalculator.cc ZbbAnalysis/AnalysisStep/plugins/ZbbAcceptanceCalculator.cc
   
 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Musich,40 2-A16,+41227671519,
//         Created:  Tue Nov 22 10:30:38 CET 2011
// $Id: ZbbAcceptanceCalculator.cc,v 1.4 2011/12/07 09:49:44 emiglior Exp $
//
//

// system include files
#include <memory>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/Parse.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Common/interface/GetProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimGeneral/HepPDTRecord/interface/PDTRecord.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "ZbbAnalysis/AnalysisStep/interface/LumiReWeighting.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbUtils.h"
#include "ZbbAnalysis/AnalysisStep/interface/AcceptanceCuts.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbTypeDefs.h"
//#include "ZbbAnalysis/AnalysisStep/interface/PU.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"
#include <boost/foreach.hpp>
#include <boost/bind.hpp>

//
// Root utils
//
#include "TEfficiency.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMath.h"

struct EventProperties
{
  
  Bool_t is_right_flavour_,is_gen_Z_yes_,is_gen_j_yes_,is_gen_b_yes_,is_gen_Zkin_yes_,is_rec_Zkin_yes_,is_rec_b_yes_,is_rec_j_yes_;
  Bool_t is_rec_lep_idiso_yes_,is_rec_b_HE_yes_,is_rec_b_HP_yes_;            
  Int_t nBFlavourJets_ , nOkBFlavouredJets_,nRecoJets_,nOkRecoJets_;
  Float_t mLL_,evtWeight_,evtWeightSq_,lepIdIsoWeight_,lepIdIsoWeightSq_,bTagWeightHE_,bTagWeightHESq_,bTagWeightHP_,bTagWeightHPSq_;
								
  EventProperties() {
    // std::cout<<"Event properties exist!"<<std::endl;
    this->reset();
  }

  void reset(){
    
    //initialize bools and counters
    is_right_flavour_=false;
    is_gen_Z_yes_=false;
    is_gen_b_yes_=false;
    is_gen_j_yes_=false;
    is_gen_Zkin_yes_=false;
    is_rec_Zkin_yes_=false;
    is_rec_b_yes_=false;  
    is_rec_j_yes_=false;  
    is_rec_lep_idiso_yes_=false;
    is_rec_b_HE_yes_=false;
    is_rec_b_HP_yes_=false;

    nBFlavourJets_=0; 
    nOkBFlavouredJets_=0;
    nOkRecoJets_=0;
    nRecoJets_=0;
    mLL_ = -1.;
    evtWeightSq_=1.;
    evtWeight_=1.;
    lepIdIsoWeight_=1;
    lepIdIsoWeightSq_=1;
    bTagWeightHE_=1;
    bTagWeightHESq_=1;
    bTagWeightHP_=1;
    bTagWeightHPSq_=1;

    //  if(!is_gen_Z_yes_)    std::cout<<"is_gen_Z_yes_    = false"<<std::endl;
    //  if(!is_gen_b_yes_)    std::cout<<"is_gen_b_yes_    = false"<<std::endl;
    //  if(!is_gen_Zkin_yes_) std::cout<<"is_gen_Zkin_yes_ = false"<<std::endl;
    //  if(!is_rec_Zkin_yes_) std::cout<<"is_rec_Zkin_yes_ = false"<<std::endl;
    //  if(!is_rec_b_yes_)    std::cout<<"is_rec_b_yes_    = false"<<std::endl;
  }

  // getter methods
  Bool_t isRightFlavour() const {return is_right_flavour_;}
  Bool_t isGenZyes() const {return is_gen_Z_yes_;}
  Bool_t isGenbyes() const {return is_gen_b_yes_;}
  Bool_t isGenjetyes()    const {return is_gen_j_yes_;}
  Bool_t isGen_Zkin_yes() const {return is_gen_Zkin_yes_;}
  Bool_t isRec_Zkin_yes() const {return is_rec_Zkin_yes_;}
  Bool_t isRecjet_yes()   const {return is_rec_j_yes_;}
  Bool_t isRecb_yes()    const {return is_rec_b_yes_;}
  Bool_t isRecLepIdIso_yes() const {return is_rec_lep_idiso_yes_;}
  Bool_t isRecbHE_yes() const {return is_rec_b_HE_yes_;}
  Bool_t isRecbHP_yes() const {return is_rec_b_HP_yes_;}

};

//
// class declaration
//

class ZbbAcceptanceCalculator : public edm::EDAnalyzer {
   public:
      explicit ZbbAcceptanceCalculator(const edm::ParameterSet&);
      ~ZbbAcceptanceCalculator();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
 
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual std::vector<math::XYZTLorentzVectorD> sort4MomentaByPt(std::vector<math::XYZTLorentzVectorD> leptonMomenta_);
      virtual std::vector<const reco::Candidate*> sortLeptonsByPt(std::vector<const reco::Candidate*> leptons_);
      virtual std::pair<bool,bool> WhichFlavour(std::vector<const reco::Candidate*> theGenMuons,std::vector<const reco::Candidate*> theGenElectrons);
      virtual std::pair<const reco::Candidate*,const reco::Candidate*> getTheZLeptons(std::vector<const reco::Candidate*> theGenMuons,std::vector<const reco::Candidate*> theGenElectrons);
      virtual std::vector< std::pair<const pat::Jet*,const reco::GenJet*> > matchJetsByDr( std::vector<const pat::Jet*> theRecoJets, std::vector<const reco::GenJet*> theGenJets, double const& maxDR, bool const& uniqueFirst,  bool const& uniqueSecond);
      template<class T, class U>  std::vector< std::pair<T,U> > MatchByDR(std::vector<T> const& c1,std::vector<U> const& c2,double const& maxDR,bool const& uniqueFirst,bool const& uniqueSecond);
      template<class T, class U>  bool TemplDRCompare(std::pair<T,U> p1, std::pair<T,U> p2);
      virtual bool DRCompare(std::pair<const pat::Jet*,const reco::GenJet*>  p1, std::pair<const pat::Jet*,const reco::GenJet*>  p2);					    
      virtual std::vector<reco::CompositeCandidate> sortCandidatesByDifference(std::vector<reco::CompositeCandidate> unsortedCands);
      virtual bool hasBottom(const reco::Candidate &c);
      virtual bool hasBottomIC(const reco::Candidate &c);
      virtual bool hasCharm(const reco::Candidate &c); 
      virtual Bool_t isTightZCandidate(reco::CompositeCandidate ZCand, const reco::BeamSpot& beamSpot, Bool_t isMuChannel, const AcceptanceCuts& lCuts);
      virtual Double_t getMuonScaleFactor(reco::CompositeCandidate &theRecoZcand);
      virtual Double_t getElectronScaleFactor(reco::CompositeCandidate &theRecoZcand);

  // ----------member data ---------------------------

  ofstream outfile_;

  edm::InputTag genJetSrc_;
  edm::InputTag genParticles_;
  edm::InputTag MuonCollection_;
  edm::InputTag ElectronCollection_;
  edm::InputTag JetCollection_;
  edm::InputTag ZmmCollection_;
  edm::InputTag ZeeCollection_;
  
  bool unLockDefaults_;
  Double_t minMassCut_;     
  Double_t maxMassCut_;     
  Double_t jetEtaCut_;      
  Double_t muonEtaCut_;     
  Double_t eleEtaCut_;      
  Double_t jetPtCut_;       
  Double_t muonPtCut_;      
  Double_t elePtCut_;       
  Double_t dRLeptons_;    
  Double_t dRJets_; 
  Double_t theBparticlePtCut_;

  //defining acceptance cuts
  AcceptanceCuts lCuts_;

  // some bools
  bool verbose_;
  bool saveNTuple_;
  bool isMCatNLO_; 
  bool applypucorrection_;
  bool useClopperPearsonErrors_;
  bool partonLevel_;
  
  // flavour of decay
  std::string decay_chain_string_;

  // weights for PU reweighting
  edm::LumiReWeighting* theLumiWeights_;
  
  // vpsets for effb (MC) maps	 
  edm::VParameterSet  effbcMC_SSVHE;
  edm::VParameterSet  effbcMC_SSVHP;
  // associative map for b-tag efficiency
  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffbMCSSVHE_;
  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffbMCSSVHP_;
  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffcMCSSVHE_;
  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffcMCSSVHP_;
  
  // unweighted counters
  EventProperties myBools_;

  // unweighted counters
  Double_t N_tot_events, N_right_flavour; 
  Double_t N_gen_Z_yes,N_gen_j_yes,N_gen_b_yes, N_gen_Zkin_yes;	
  Double_t Al_Num,Al_Den,CHad_Num,Al_Num_alljets,CHad_Num_alljets;
  Double_t N_rec_Zkin_yes,N_rec_b_yes,N_rec_j_yes;      
  Double_t N_rec_y_gen_n, N_rec_y_gen_y, N_rec_n_gen_n, N_rec_n_gen_y;
  Double_t Epsilon_l_Num, Epsilon_l_Den, Epsilon_l_Num_alljets, Epsilon_l_Den_alljets;
  Double_t Epsilon_bHE_Num, Epsilon_bHE_Den; 
  Double_t Epsilon_bHP_Num, Epsilon_bHP_Den;  

  // weighted counters
  Double_t sumW_tot_events,sumW_right_flavour;
  Double_t sumW_gen_Z_yes,sumW_gen_b_yes,sumW_gen_j_yes,sumW_gen_Zkin_yes;
  Double_t sumW_Al_Num,sumW_Al_Den,sumW_CHad_Num,sumW_Al_Num_alljets,sumW_CHad_Num_alljets; 		    
  Double_t sumW_rec_Zkin_yes,sumW_rec_b_yes,sumW_rec_j_yes;
  Double_t sumW_rec_y_gen_n,sumW_rec_y_gen_y,sumW_rec_n_gen_n,sumW_rec_n_gen_y;
  Double_t sumW_Epsilon_l_Num, sumW_Epsilon_l_Den, sumW_Epsilon_l_Num_alljets, sumW_Epsilon_l_Den_alljets;
  Double_t sumW_Epsilon_bHE_Num, sumW_Epsilon_bHE_Den; 
  Double_t sumW_Epsilon_bHP_Num, sumW_Epsilon_bHP_Den;  

  // weighted square counters
  Double_t sumW2_tot_events,sumW2_right_flavour;
  Double_t sumW2_gen_Z_yes,sumW2_gen_b_yes,sumW2_gen_j_yes,sumW2_gen_Zkin_yes;
  Double_t sumW2_Al_Num,sumW2_Al_Den,sumW2_CHad_Num,sumW2_Al_Num_alljets,sumW2_CHad_Num_alljets;
  Double_t sumW2_rec_Zkin_yes,sumW2_rec_b_yes,sumW2_rec_j_yes;
  Double_t sumW2_rec_y_gen_n,sumW2_rec_y_gen_y,sumW2_rec_n_gen_n,sumW2_rec_n_gen_y;
  Double_t sumW2_Epsilon_l_Num, sumW2_Epsilon_l_Den, sumW2_Epsilon_l_Num_alljets, sumW2_Epsilon_l_Den_alljets;
  Double_t sumW2_Epsilon_bHE_Num, sumW2_Epsilon_bHE_Den;
  Double_t sumW2_Epsilon_bHP_Num, sumW2_Epsilon_bHP_Den;  

  // flavour selection
  bool isZEE, isZMM, isZTT;
  
  // histograms
  TH1F *h1UnwCategories_;
  TH1F *h1WgtCategories_;
  
  TH1F *h_GEN_mass_ZfromMuons_, *h_GEN_mass_ZfromElectrons_, *h_GEN_nBflavouredJets_, *h_GEN_nOkBflavouredJets_;  
  TH1F *h_GEN_GenbGenJetMinDeltaR_;

  TH1F *h_selbjet_GenRecoResponse_, *h_bjet_GenRecoResponse_;     
  TH2F *h_bjet_GenPtvsRecoPt_,*h_GEN_GenbGenJetMinDeltaRVsEta_,*h_GEN_GenbGenJetMinDeltaRVsPt_;
  TH2F *h_selbjet_GenPtvsRecoPt_ ,*h_selbjet_GenReco_DeltaR_vsEta_,*h_selbjet_GenReco_DeltaR_vsPt_;

  TH1F *h_bparticlePt_, *h_bparticleEta_, *h_bparticleEnergyFraction_;           

  TH1F *h_bjet_GenPt_,*h_bjet_GenEta_,*h_bjet_GenMass_; 
  TH1F *h_bjet_GenMatchedPt_, *h_bjet_GenMatchedEta_, *h_bjet_GenMatchedMass_;
  TH1F *h_bjet_RecoMatchedPt_,*h_bjet_RecoMatchedEta_,*h_bjet_RecoMatchedMass_;     
  TH1F *h_selbjet_RecoMatchedPt_,*h_selbjet_RecoMatchedEta_,*h_selbjet_RecoMatchedMass_;  

  TH1F *h_bjet_GenReco_DeltaR_,*h_selbjet_GenReco_DeltaR_,*h_mu_GenReco_DeltaR_,*h_ele_GenReco_DeltaR_;  

  TH2F *h_bjet_GenReco_DeltaR_vsPt_, *h_bjet_GenReco_DeltaR_vsEta_;
  TH2F *h_mu_GenReco_DeltaR_vsPt_,*h_mu_GenReco_DeltaR_vsEta_;  
  TH2F *h_ele_GenReco_DeltaR_vsPt_, *h_ele_GenReco_DeltaR_vsEta_; 

  TH1F *h_ele_GenPt_,*h_ele_GenEta_,*h_ele_RecoPt_,*h_ele_RecoEta_;              
  TH1F *h_mu_GenPt_,*h_mu_GenEta_,*h_mu_RecoPt_,*h_mu_RecoEta_;              

  TH1F *h_eleMatched_GenPt_,*h_eleMatched_GenEta_,*h_eleMatched_RecoPt_,*h_eleMatched_RecoEta_;              
  TH1F *h_muMatched_GenPt_,*h_muMatched_GenEta_,*h_muMatched_RecoPt_,*h_muMatched_RecoEta_;    

  TH1F *h_RECO_mass_ZfromMuons_,*h_RECO_mass_ZfromElectrons_;  
 
  TH2F * h_GENMatrix_;                
  TH2F * h_GENKinMatrix_;             
  TH2F * h_RECOMatrix_;            

  TTree *myNTuple_;

  // define the 4vector of the Z and of the leptons
  math::XYZTLorentzVectorD pLpos, pLneg, pZLL;

};

//
// constructors and destructor
//

ZbbAcceptanceCalculator::ZbbAcceptanceCalculator(const edm::ParameterSet& iConfig){
  
  //get the collections
  genParticles_      =iConfig.getParameter<edm::InputTag>("GenSrc");
  MuonCollection_    =iConfig.getParameter<edm::InputTag>("MuonCollection");
  ElectronCollection_=iConfig.getParameter<edm::InputTag>("ElectronCollection");
  JetCollection_     =iConfig.getParameter<edm::InputTag>("JetCollection");
  ZmmCollection_     =iConfig.getParameter<edm::InputTag>("ZmmCollection");
  ZeeCollection_     =iConfig.getParameter<edm::InputTag>("ZeeCollection");
  genJetSrc_         =iConfig.getParameter<edm::InputTag>("genJetSrc");

  // define the cuts
  unLockDefaults_    =iConfig.getParameter<bool>("unLockDefaultCuts");
  minMassCut_        =iConfig.getParameter<double>("minMassCut");
  maxMassCut_        =iConfig.getParameter<double>("maxMassCut");
  jetEtaCut_         =iConfig.getParameter<double>("jetEtaCut");
  muonEtaCut_        =iConfig.getParameter<double>("muonEtaCut");
  eleEtaCut_         =iConfig.getParameter<double>("eleEtaCut");
  jetPtCut_          =iConfig.getParameter<double>("jetPtCut");
  muonPtCut_         =iConfig.getParameter<double>("muonPtCut");
  elePtCut_          =iConfig.getParameter<double>("elePtCut");
  dRLeptons_         =iConfig.getParameter<double>("dRLeptonMatch");
  dRJets_            =iConfig.getParameter<double>("dRJetMatch");
  theBparticlePtCut_ =iConfig.getParameter<double>("BparticlePtCut");

  if(!unLockDefaults_){
    minMassCut_= 60.;    
    maxMassCut_= 120.;    
    jetEtaCut_ = 2.1;     
    muonEtaCut_= 2.1;    
    eleEtaCut_ = 2.5;     
    jetPtCut_  = 25.;      
    muonPtCut_ = 20.;     
    elePtCut_  = 25.;      
    dRLeptons_ = 0.3; 
    dRJets_    = 0.5;   
  }

  decay_chain_string_      = iConfig.getUntrackedParameter<std::string>("DecayChainSelection");
  isMCatNLO_               = iConfig.getUntrackedParameter<bool>("isMCatNLO", false);
  applypucorrection_       = iConfig.getUntrackedParameter<bool>("applyPUCorr", true); 
  verbose_                 = iConfig.getUntrackedParameter<bool>("verbose", false);
  saveNTuple_              = iConfig.getUntrackedParameter<bool>("saveNTuple", false);
  useClopperPearsonErrors_ = iConfig.getUntrackedParameter<bool>("useClopperPearsonErrors", false);
  partonLevel_             = iConfig.getUntrackedParameter<bool>("PartonLevel", false);
  effbcMC_SSVHE            = iConfig.getParameter<edm::VParameterSet>("EffbcMCSSVHE");
  effbcMC_SSVHP            = iConfig.getParameter<edm::VParameterSet>("EffbcMCSSVHP");

  // unweighted counters
  N_tot_events=0.; N_right_flavour=0.;
  N_gen_Z_yes=0.; N_gen_b_yes=0.; N_gen_Zkin_yes=0.;N_gen_j_yes=0.;	
  Al_Num=0.; Al_Den=0.; CHad_Num=0.; Al_Num_alljets=0; CHad_Num_alljets=0;
  N_rec_Zkin_yes=0.; N_rec_b_yes=0.; N_rec_j_yes=0.;    
  N_rec_y_gen_n=0.; N_rec_y_gen_y=0.; N_rec_n_gen_n=0.; N_rec_n_gen_y=0.;
  Epsilon_l_Num=0; Epsilon_l_Den=0; Epsilon_l_Num_alljets=0; Epsilon_l_Den_alljets=0;
  Epsilon_bHE_Num=0; Epsilon_bHE_Den=0;
  Epsilon_bHP_Num=0; Epsilon_bHP_Den=0;

  // weighted counters
  sumW_tot_events=0.; sumW_right_flavour=0.;
  sumW_gen_Z_yes=0.; sumW_gen_b_yes=0.;	 sumW_gen_Zkin_yes=0.; sumW_gen_j_yes=0.;	
  sumW_Al_Num=0.; sumW_Al_Den=0.; sumW_CHad_Num=0.; sumW_Al_Num_alljets=0.; sumW_CHad_Num_alljets=0.; 
  sumW_rec_Zkin_yes=0.; sumW_rec_b_yes=0.; sumW_rec_j_yes=0.; 
  sumW_rec_y_gen_n=0.; sumW_rec_y_gen_y=0.; sumW_rec_n_gen_n=0.; sumW_rec_n_gen_y=0.;
  sumW_Epsilon_l_Num=0; sumW_Epsilon_l_Den=0; sumW_Epsilon_l_Num_alljets=0; sumW_Epsilon_l_Den_alljets=0;
  sumW_Epsilon_bHE_Num=0; sumW_Epsilon_bHE_Den=0;
  sumW_Epsilon_bHP_Num=0; sumW_Epsilon_bHP_Den=0;

  // squared weighted counters
  sumW2_tot_events=0.; sumW2_right_flavour=0.;
  sumW2_gen_Z_yes=0.; sumW2_gen_b_yes=0.; sumW2_gen_Zkin_yes=0.; sumW2_gen_j_yes=0.;	
  sumW2_Al_Num=0.; sumW2_Al_Den=0.; sumW2_CHad_Num=0.; sumW2_Al_Num_alljets=0.; sumW2_CHad_Num_alljets=0.; 
  sumW2_rec_Zkin_yes=0.; sumW2_rec_b_yes=0.; sumW2_rec_j_yes=0.;
  sumW2_rec_y_gen_n=0.; sumW2_rec_y_gen_y=0.; sumW2_rec_n_gen_n=0.; sumW2_rec_n_gen_y=0.;
  sumW2_Epsilon_l_Num=0; sumW2_Epsilon_l_Den=0; sumW2_Epsilon_l_Num_alljets=0; sumW2_Epsilon_l_Den_alljets=0;
  sumW2_Epsilon_bHE_Num=0; sumW2_Epsilon_bHE_Den=0;
  sumW2_Epsilon_bHP_Num=0; sumW2_Epsilon_bHP_Den=0;

  //initialize the lumiweight
  theLumiWeights_=0;

  //select the flavour of leptons from Z
  isZEE=false; 
  isZMM=false; 
  isZTT=false;

  std::vector<std::string> tokens = edm::tokenize(decay_chain_string_, "_");
  if( tokens.size() != 2 ){
    std::cout<<"Invalid decay chain string"<<std::endl;
  } else {
    if( tokens[1] == "m" ) {
      std::cout<<"Selected Z->mm channel"<<std::endl;
      isZMM=true;
    } else if( tokens[1] == "e" ) {
      std::cout<<"Selected Z->ee channel"<<std::endl;
      isZEE=true; 
    } else {
      std::cout<<"Invalid decay chain string"<<std::endl;
      isZTT=true;
    }
  }

  // b-tag MC efficiency SSVHE tagger
  EffMCpair theEffbSSVHE_;
  EffMCpair theEffcSSVHE_;
  for(edm::VParameterSet::const_iterator pset =  effbcMC_SSVHE.begin(); pset != effbcMC_SSVHE.end(); pset++){
    std::vector<double> thePtRange_ = pset->getUntrackedParameter<std::vector<double> >("ptrange");
    double theEffb_barrel_  = pset->getUntrackedParameter<double>("effb_barrel");
    double theEffb_forward_ = pset->getUntrackedParameter<double>("effb_forward");
    double theEffc_barrel_  = pset->getUntrackedParameter<double>("effc_barrel");
    double theEffc_forward_ = pset->getUntrackedParameter<double>("effc_forward");
    theEffbSSVHE_ = std::make_pair(theEffb_barrel_,theEffb_forward_);
    theEffcSSVHE_ = std::make_pair(theEffc_barrel_,theEffc_forward_);
    theAssociativeMapEffbMCSSVHE_[theEffbSSVHE_]=thePtRange_; 
    theAssociativeMapEffcMCSSVHE_[theEffcSSVHE_]=thePtRange_; 
  }

  // b-tag MC efficiency SSVHP tagger
  EffMCpair theEffbSSVHP_;
  EffMCpair theEffcSSVHP_;
  for(edm::VParameterSet::const_iterator pset =  effbcMC_SSVHP.begin(); pset != effbcMC_SSVHP.end(); pset++){
    std::vector<double> thePtRange_ = pset->getUntrackedParameter<std::vector<double> >("ptrange");
    double theEffb_barrel_  = pset->getUntrackedParameter<double>("effb_barrel");
    double theEffb_forward_ = pset->getUntrackedParameter<double>("effb_forward");
    double theEffc_barrel_  = pset->getUntrackedParameter<double>("effc_barrel");
    double theEffc_forward_ = pset->getUntrackedParameter<double>("effc_forward");
    theEffbSSVHP_ = std::make_pair(theEffb_barrel_,theEffb_forward_);
    theEffcSSVHP_ = std::make_pair(theEffc_barrel_,theEffc_forward_);
    theAssociativeMapEffbMCSSVHP_[theEffbSSVHP_]=thePtRange_; 
    theAssociativeMapEffcMCSSVHP_[theEffcSSVHP_]=thePtRange_; 
  }
  
  //setting the offline cuts
  lCuts_.setDefault();
  
  if(unLockDefaults_){
    lCuts_.clear();
    lCuts_.set(jetEtaCut_,jetPtCut_,muonEtaCut_,muonPtCut_,eleEtaCut_,elePtCut_);
  }
}

ZbbAcceptanceCalculator::~ZbbAcceptanceCalculator()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void
ZbbAcceptanceCalculator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
   
  //------------------------ Handles -------------------------------
  
  // get gen particle candidates
  edm::Handle<reco::GenParticleCollection> genParticlesCollection;
  iEvent.getByLabel(genParticles_, genParticlesCollection);
  
  // get patMuons
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(MuonCollection_,muons);
  
  // get patElectrons
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByLabel(ElectronCollection_,electrons);
  
  // get reco patJets
  edm::Handle<edm::View<pat::Jet> > recoJetsHandle;
  iEvent.getByLabel(JetCollection_,recoJetsHandle);
  const edm::View<pat::Jet> & recoJets = *(recoJetsHandle.product());
  
  // get genJets
  edm::Handle<reco::GenJetCollection> genJetsHandle;
  iEvent.getByLabel(genJetSrc_,genJetsHandle);
  const reco::GenJetCollection & genJets = *(genJetsHandle.product());
  
  // get reco Zmm candidates
  Handle<reco::CompositeCandidateCollection> zmmHandle;
  iEvent.getByLabel(ZmmCollection_, zmmHandle);
  
  // get reco Zee candidates
  Handle<reco::CompositeCandidateCollection> zeeHandle;
  iEvent.getByLabel(ZeeCollection_, zeeHandle);
  
  // get event weight 
  edm::Handle<GenEventInfoProduct> genEventInfoHandle;
  iEvent.getByLabel("generator",genEventInfoHandle);
  double wMC = genEventInfoHandle->weight();

  // get beamSpot
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  if ( beamSpotHandle.isValid() ) {
    beamSpot = *beamSpotHandle;
    if( verbose_ )  std::cout<<"Beamspot x: "<< beamSpot.x0() <<" y: "<< beamSpot.y0() <<" z: "<< beamSpot.z0() <<std::endl;
  } else {
    if( verbose_ )  std::cout << "No BeamSpot found!" << std::endl;
  }

  //------------------------ reweighting ----------------------------------

  // reweighting rule for NLO
  // if(isMCatNLO_) {
  //     wMC > 0 ?  wMC=1. : wMC=-1.;
  //   }

  // applying PU corrections
  if(applypucorrection_){
    double thePUWeight_(1.);
    const edm::EventBase* iEventB = dynamic_cast<const edm::EventBase*>(&iEvent);
    if (theLumiWeights_) thePUWeight_ = theLumiWeights_->weight(*iEventB);
    wMC*=thePUWeight_;
  }

  double wLepIdIso(1.);
  double wBtagHE(1.);
  double wBtagHP(1.);

  myBools_.reset();

  N_tot_events++;
  sumW_tot_events+=wMC;
  sumW2_tot_events+=(wMC*wMC);

  myBools_.evtWeight_=wMC;
  myBools_.evtWeightSq_=(wMC*wMC);
  
  //------------------------ define leptons -------------------------------

  int pdgIdL(0), statusL(0); 
  Double_t etaCut(0.), ptCut(0.), jetEtaCut(0.), jetPtCut(0.), minMassCut(0.), maxMassCut(0.);

  minMassCut = minMassCut_;     
  maxMassCut = maxMassCut_;     
  jetEtaCut  = jetEtaCut_; 
  jetPtCut   = jetPtCut_;
 
  if ( isZEE ) {
    pdgIdL = 11;
    statusL = 3;
    etaCut = eleEtaCut_;
    ptCut  = elePtCut_ ;
  } else if ( isZMM ) {
    pdgIdL = 13;
    statusL = 3;
    etaCut = muonEtaCut_;
    ptCut =  muonPtCut_;
  } else if ( isZTT ) {
    pdgIdL = 15;
    statusL = 2;
  }

  std::vector<const reco::Candidate*> theGenMuons;
  std::vector<const reco::Candidate*> theGenMuons_Status3;
  std::vector<math::XYZTLorentzVectorD> theGenMuonsMomenta;
  
  std::vector<const reco::Candidate*> theGenElectrons;
  std::vector<math::XYZTLorentzVectorD> theGenElectronsMomenta;

  std::vector<const reco::GenJet*> the_b_flavour_genJets;
  std::vector<const reco::GenJet*> the_any_flavour_genJets;
  std::vector<const reco::GenJet*> the_filtered_any_flavour_genJets;
  std::vector<const reco::GenJet*> the_filtered_b_flavour_genJets;

  //------------------------ master loop on genParticles -------------------------------

  // initialize the gen 4vectors
  pLpos.SetXYZT(0,0,0,0);
  pLneg.SetXYZT(0,0,0,0);
  pZLL.SetXYZT(0,0,0,0);

  // lepton matching
  for( reco::GenParticleCollection::const_iterator genp = genParticlesCollection->begin();genp != genParticlesCollection->end(); ++ genp ) { 
    
    // make a copy of genp for better handling
    const reco::Candidate* p = &(*genp);
    
    if(isZMM){
      if( abs( p->pdgId() )==pdgIdL && ( p->status())==statusL ){
	theGenMuons.push_back(p);
	math::XYZTLorentzVectorD p4Lepttemp(p->px(),p->py(),p->pz(),p->energy());  
	theGenMuonsMomenta.push_back(p4Lepttemp);
      }
      if( abs( p->pdgId() )==pdgIdL && ( p->status())==3 ){
	theGenMuons_Status3.push_back(p);
      }

    } else if(isZEE) {
      if( abs( p->pdgId() )==pdgIdL && ( p->status())==statusL ){
	theGenElectrons.push_back(p);
	math::XYZTLorentzVectorD p4Lepttemp(p->px(),p->py(),p->pz(),p->energy());  
	theGenElectronsMomenta.push_back(p4Lepttemp);
      }
    } else {  
      std::cout<<"This case is not contemplated at the moment!"<<std::endl;
    } // ends swithces on flavour
  } // ends lepton master loop on genParticles

  // b-jet matching

  // parton level matching
  if ( partonLevel_ ){

    // loop on genJets
    for(std::vector<reco::GenJet>::const_iterator genjet_it=genJets.begin(); genjet_it!=genJets.end(); ++genjet_it){ 
      // loop on genJets constituents to check for a b
      std::vector <const reco::GenParticle*> mcparts = genjet_it->getGenConstituents();
      // for each MC particle in turn  

      Bool_t foundAb(false);
      double dR_BpartBjet_(9999.);

      for (unsigned i = 0; i < mcparts.size (); i++) {
	const reco::GenParticle* mcpart = mcparts[i];
	if ( mcpart->pdgId() == -5 || mcpart->pdgId() == 5 ) {
	  dR_BpartBjet_ = ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(), mcpart->momentum()); 
	  foundAb=true;
	  if(verbose_) std::cout<< "b parton found" <<std::endl;	  
	  h_bparticlePt_->Fill(mcpart->pt(),wMC);
	  h_bparticleEta_->Fill(mcpart->eta(),wMC);
	  h_bparticleEnergyFraction_->Fill(mcpart->pt()/genjet_it->pt(),wMC); 
	}
      }

      if(foundAb){
	h_GEN_GenbGenJetMinDeltaR_->Fill(dR_BpartBjet_,wMC);
	h_GEN_GenbGenJetMinDeltaRVsEta_->Fill(genjet_it->eta(),dR_BpartBjet_,wMC); 
	h_GEN_GenbGenJetMinDeltaRVsPt_->Fill(genjet_it->pt(),dR_BpartBjet_,wMC);
	h_bjet_GenPt_->Fill(genjet_it->pt(),wMC);
	h_bjet_GenEta_->Fill(genjet_it->eta(),wMC);
	h_bjet_GenMass_->Fill(genjet_it->mass(),wMC);
	the_b_flavour_genJets.push_back(&(*genjet_it));
      }
    }
  } else {
    // hadron level matching
    for( reco::GenParticleCollection::const_iterator genp = genParticlesCollection->begin();genp != genParticlesCollection->end(); ++ genp ) { 
      
      // make a copy of genp for better handling
      const reco::Candidate* p = &(*genp);
      
      // if b-hadrons
      if( hasBottomIC(*p) ){ 
	
	// checks if all daughter are bottomless
	Bool_t hasBottomedDaughter = false;
	for(UInt_t i=0; i<p->numberOfDaughters(); i++){
	  if(hasBottomIC(*(p->daughter(i)))) hasBottomedDaughter=true;
	}
	
	// if last-bhadron
	if(!hasBottomedDaughter &&  p->pt() > theBparticlePtCut_){
	  
	  h_bparticlePt_->Fill(p->pt(),wMC);
	  h_bparticleEta_->Fill(p->eta(),wMC);

	  math::XYZTLorentzVectorD p4HFGEN(p->px(),p->py(),p->pz(),p->energy());
	  
	  // indices for b-hadron association
	  std::pair<int,double>  GenGenAssociation = std::make_pair(-1,9999.);
	  std::pair<const reco::GenJet*,double> GenJetAssociation;
	  if(genJets.size()>0) GenJetAssociation = std::make_pair(&genJets.at(0),9999.);
	  
	  double minDeltaRGenGen(9999.);
	  int i(0);
	  
	  for(std::vector<reco::GenJet>::const_iterator genjet_it=genJets.begin(); genjet_it!=genJets.end(); ++genjet_it){ 
	    double dR_tmp = ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4HFGEN.Vect()); 
	    //if( dR_tmp<dRJets_ ) { 
	    if( dR_tmp<minDeltaRGenGen ){
	      minDeltaRGenGen = dR_tmp;
	      GenGenAssociation.first  = i;
	      GenGenAssociation.second = minDeltaRGenGen;
	      
	      GenJetAssociation.first  = &(*genjet_it);
	      GenJetAssociation.second = minDeltaRGenGen;
	    } // end if minimum distance gen jet-parton
	    if(verbose_) std::cout<<"i: "<<i<<" minDeltaR: "<<minDeltaRGenGen<<std::endl;
	    //} // ends if b-hadron loosely matched to gen jet
	    i++;
	  }// ends loop on gen jets
	  
	  if( GenGenAssociation.first != -1 ){
	    if(GenJetAssociation.second<dRJets_ ) { 
	      the_b_flavour_genJets.push_back(GenJetAssociation.first);
	      h_bparticleEnergyFraction_->Fill(p4HFGEN.pt()/(GenJetAssociation.first)->pt(),wMC); 
	    }
	    h_GEN_GenbGenJetMinDeltaR_->Fill(GenJetAssociation.second,wMC);
	    h_GEN_GenbGenJetMinDeltaRVsEta_->Fill((GenJetAssociation.first)->eta(),GenJetAssociation.second,wMC); 
	    h_GEN_GenbGenJetMinDeltaRVsPt_->Fill((GenJetAssociation.first)->pt(),GenJetAssociation.second,wMC);
	    h_bjet_GenPt_->Fill((GenJetAssociation.first)->pt(),wMC);
	    h_bjet_GenEta_->Fill((GenJetAssociation.first)->eta(),wMC);
	    h_bjet_GenMass_->Fill((GenJetAssociation.first)->mass(),wMC);
	  } // pushes back vector of b-flavoured genJets	
	} // if last b-hadron 
      } // if bottomed particle
    } // end second master loop on genParticles
  } // end case parton/hadron level corrections
  

  if(verbose_) std::cout<<"before muon or electron case"<<std::endl;

  // ========================= Muons Case ===============================
  //Muon fix: we want the status 1 muons that we actually detect.  At this
  //point it is given that two status 3 muons were found, we replace these
  //with the status 1 pair having a mass closest to the Z.  If for some reason
  //two status 1 muons are not found, discard this event.
  
  if (isZMM && theGenMuons_Status3.size()==2){

    myBools_.is_right_flavour_=true;
    N_right_flavour++;
    sumW_right_flavour+=wMC;
    sumW2_right_flavour+=(wMC*wMC);
    
    if(theGenMuons.size()>1){
      if(verbose_) std::cout<<"theGenMuon.size()="<<theGenMuons.size()<<std::endl;
      
      std::vector<math::XYZTLorentzVectorD>  sortedMuonsMomenta = sort4MomentaByPt(theGenMuonsMomenta);
      std::vector<const reco::Candidate*>    sortedGenMuons     = sortLeptonsByPt(theGenMuons);
      
      std::vector<reco::CompositeCandidate> unsortedCands;
      
      // prepare vector of unsorted composite candidates
      for (std::vector<const reco::Candidate*>::const_iterator muonsBegin = sortedGenMuons.begin(),muonsEnd = sortedGenMuons.end(),imuon = muonsBegin;imuon != muonsEnd; ++imuon ) {
	for ( std::vector<const reco::Candidate*>::const_iterator jmuon = imuon + 1; jmuon != muonsEnd; ++jmuon ) {
	  if ( (*imuon)->charge() * (*jmuon)->charge() < 0 ) {
	    
	    reco::CompositeCandidate ZCand;
	    ZCand.addDaughter( **imuon, "mu1");
	    ZCand.addDaughter( **jmuon, "mu2");
	    
	    AddFourMomenta addp4;
	    addp4.set( ZCand );      
	    
	    unsortedCands.push_back(ZCand);
	  } // if opposite charge
	} // loop on second muon
      } // loop on first muon
      
      if(unsortedCands.size() >0){
	
	reco::CompositeCandidate theBestZCand = sortCandidatesByDifference(unsortedCands)[0];
	
	if ( theBestZCand.daughter(0)->pdgId() < 0 ) {
	  pLneg=theBestZCand.daughter(0)->p4();
	  pLpos=theBestZCand.daughter(1)->p4();
	} else if ( theBestZCand.daughter(0)->pdgId() > 0  ) {
	  pLneg=theBestZCand.daughter(1)->p4();
	  pLpos=theBestZCand.daughter(0)->p4();
	}
	
	pZLL=pLpos+pLneg;
	Double_t mLL=pZLL.mass();
	myBools_.mLL_=mLL;
	h_GEN_mass_ZfromMuons_->Fill(mLL,wMC);
	
	h_mu_GenPt_->Fill(pLpos.pt(),wMC);                
	h_mu_GenEta_->Fill(pLpos.eta(),wMC);   
	h_mu_GenPt_->Fill(pLneg.pt(),wMC);                
	h_mu_GenEta_->Fill(pLneg.eta(),wMC);  
	
	if(mLL>minMassCut && mLL<maxMassCut){
	  myBools_.is_gen_Z_yes_=true;
	} 
      }
    } // everything need to be incapsulated in the if there at least 2 genMuons (status 1)
    
    // ========================= Electrons Case ===============================      
  } else if (isZEE && theGenElectrons.size()==2 ) {

    myBools_.is_right_flavour_=true;
    N_right_flavour++;

    sumW_right_flavour+=wMC;
    sumW2_right_flavour+=(wMC*wMC);

    std::vector<math::XYZTLorentzVectorD>  sortedElectronsMomenta = sort4MomentaByPt(theGenElectronsMomenta);
    std::vector<const reco::Candidate*>    sortedGenElectrons     = sortLeptonsByPt(theGenElectrons);
    
    if( (theGenElectrons[0]->pdgId()*theGenElectrons[1]->pdgId() )< 0) {
      if (  theGenElectrons[0]->pdgId() < 0 ) {
	pLneg=theGenElectrons[0]->p4();
 	pLpos=theGenElectrons[1]->p4();
      } else if ( theGenElectrons[0]->pdgId() > 0  ) {
	pLneg=theGenElectrons[1]->p4();
 	pLpos=theGenElectrons[0]->p4();
      }
    } else {
      std::cout<<"This should never happen!"<<std::endl;
    }

    pZLL=pLpos+pLneg;
    Double_t mLL=pZLL.mass();
    myBools_.mLL_=mLL;
    h_GEN_mass_ZfromElectrons_->Fill(mLL,wMC);
    
    h_ele_GenPt_->Fill(pLpos.pt(),wMC);                
    h_ele_GenEta_->Fill(pLpos.eta(),wMC);   
    h_ele_GenPt_->Fill(pLneg.pt(),wMC);                
    h_ele_GenEta_->Fill(pLneg.eta(),wMC);     

    if(mLL>minMassCut && mLL<maxMassCut){
      myBools_.is_gen_Z_yes_=true;
    }
  }

  // ========================= study any flavour jets ===========================
  if(verbose_) std::cout<<"before any flavour matching"<<std::endl;
  Int_t nOkJets(0);
  Int_t nJets(0);
  for(std::vector<reco::GenJet>::const_iterator genjet_it=genJets.begin(); genjet_it!=genJets.end(); ++genjet_it){
    nJets++;
    the_any_flavour_genJets.push_back(&(*genjet_it));
    if(genjet_it->pt()>jetPtCut && TMath::Abs(genjet_it->eta())<jetEtaCut &&
       ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),pLpos.Vect())>dRJets_ &&
       ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),pLneg.Vect())>dRJets_
       ){
      nOkJets++;
      the_filtered_any_flavour_genJets.push_back(&(*genjet_it));
    }
  }

  if(nOkJets>=1){  
    myBools_.is_gen_j_yes_=true;
  }

  // ========================= study b-flavoured jets ===========================
  if(verbose_) std::cout<<"before b-flavour matching"<<std::endl;
  Int_t nOkBFlavouredJets(0);
  Int_t nBFlavourJets(0);
  for(std::vector<const reco::GenJet*>::const_iterator genjet_it=the_b_flavour_genJets.begin(); genjet_it!=the_b_flavour_genJets.end(); ++genjet_it){
    nBFlavourJets++;
    if(pLpos.mass()!=0. && pLneg.mass()!=0.){
      if((*genjet_it)->pt()>jetPtCut && TMath::Abs((*genjet_it)->eta())<jetEtaCut &&
	 ROOT::Math::VectorUtil::DeltaR((*genjet_it)->momentum(),pLpos.Vect())>dRJets_ &&
	 ROOT::Math::VectorUtil::DeltaR((*genjet_it)->momentum(),pLneg.Vect())>dRJets_
	 ){
	nOkBFlavouredJets++;
	the_filtered_b_flavour_genJets.push_back(*genjet_it);
      }
    }
  }

  h_GEN_nOkBflavouredJets_->Fill(nOkBFlavouredJets,wMC);
  h_GEN_nBflavouredJets_->Fill(nBFlavourJets,wMC);
  myBools_.nBFlavourJets_= nBFlavourJets;
  myBools_.nOkBFlavouredJets_ = nOkBFlavouredJets;

  if(nOkBFlavouredJets>=1 && myBools_.isRightFlavour()) {
    myBools_.is_gen_b_yes_=true;
  }
  
  // ========================= study kinematics of Z leptons ====================
  Int_t nOkLeptons(0);
  
  if(verbose_) std::cout<<"before Z kin ok"<<std::endl;
  if(pLpos.pt() > ptCut  && TMath::Abs(pLpos.eta()) < etaCut ) nOkLeptons++;
  if(pLneg.pt() > ptCut  && TMath::Abs(pLneg.eta()) < etaCut ) nOkLeptons++;

  if(nOkLeptons ==2) myBools_.is_gen_Zkin_yes_=true;

  // ========================= reco leptons matching ===========================
  // stored vector of lepton momenta
  std::vector<math::XYZTLorentzVectorD> recoLeptonMomenta_; 
  reco::CompositeCandidate recoZCand;
  Int_t nRecoMatchedLeptons(0);
  Int_t compositeCharge(1);

  if(verbose_) std::cout<<"before reco leptons matching"<<std::endl;
  if(isZMM && myBools_.isRightFlavour()){
    
    // indices for b-hadron association
    std::pair<int,double> posGenMuonRecoMuonAssociation = std::make_pair(-1,9999.);
    std::pair<int,double> negGenMuonRecoMuonAssociation = std::make_pair(-1,9999.);
    std::pair<pat::Muon,double> posGenMuonAssociation;
    std::pair<pat::Muon,double> negGenMuonAssociation;
    if((*(muons.product())).size()>0){ 
      posGenMuonAssociation = std::make_pair((*(muons.product())).at(0),9999.);
      negGenMuonAssociation = std::make_pair((*(muons.product())).at(0),9999.);
    }
    
    double posMinDeltaRGenMuon(9999.), negMinDeltaRGenMuon(9999.);
    int p(0),n(0);
    
    //========================== matching for positive muon =======================
    for(edm::View<pat::Muon>::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon){
      h_mu_RecoPt_->Fill(muon->pt(),wMC);   
      h_mu_RecoEta_->Fill(muon->eta(),wMC);  
      double dR_tmp = ROOT::Math::VectorUtil::DeltaR(muon->momentum(),pLpos.Vect());
      //  if( dR_tmp<dRLeptons_ ){
      if( dR_tmp<posMinDeltaRGenMuon ){
	posMinDeltaRGenMuon = dR_tmp;
	posGenMuonRecoMuonAssociation.first  = p;
	posGenMuonRecoMuonAssociation.second = posMinDeltaRGenMuon;
	
	posGenMuonAssociation.first  = (*muon);
	posGenMuonAssociation.second = posMinDeltaRGenMuon;
      } // if gen-reco are matched
	// } // loop on reco muons
      p++;
    }

    if( posGenMuonRecoMuonAssociation.first != -1 ){
      
      if(posGenMuonAssociation.second<dRLeptons_ ){
	recoLeptonMomenta_.push_back( posGenMuonAssociation.first.p4() );
	nRecoMatchedLeptons++;
	compositeCharge*=posGenMuonAssociation.first.charge();   
	h_muMatched_GenPt_->Fill(pLpos.pt(),wMC);                
	h_muMatched_GenEta_->Fill(pLpos.eta(),wMC);               
	h_muMatched_RecoPt_->Fill(posGenMuonAssociation.first.pt(),wMC);               
	h_muMatched_RecoEta_->Fill(posGenMuonAssociation.first.eta(),wMC); 
      }
      
      h_mu_GenReco_DeltaR_->Fill(posGenMuonAssociation.second,wMC);        
      h_mu_GenReco_DeltaR_vsPt_->Fill(pLpos.pt(),posGenMuonAssociation.second,wMC);   
      h_mu_GenReco_DeltaR_vsEta_->Fill(pLpos.eta(),posGenMuonAssociation.second,wMC);  

    }

    //========================== matching for negative muon =============================
    for(edm::View<pat::Muon>::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon){
      double dR_tmp = ROOT::Math::VectorUtil::DeltaR(muon->momentum(),pLneg.Vect());
      //  if( dR_tmp<dRLeptons_ ){
	if( dR_tmp<negMinDeltaRGenMuon ){
	  negMinDeltaRGenMuon = ROOT::Math::VectorUtil::DeltaR(muon->momentum(),pLneg.Vect());
	  negGenMuonRecoMuonAssociation.first  = n;
	  negGenMuonRecoMuonAssociation.second = negMinDeltaRGenMuon;
	  
	  negGenMuonAssociation.first  = (*muon);
	  negGenMuonAssociation.second = negMinDeltaRGenMuon;
	} // if gen-reco are matched
	//} // loop on reco muons
      n++;
    }

    if( negGenMuonRecoMuonAssociation.first != -1 ){
      
      if(negGenMuonAssociation.second<dRLeptons_){
	recoLeptonMomenta_.push_back( negGenMuonAssociation.first.p4() );
	nRecoMatchedLeptons++;
	compositeCharge*=negGenMuonAssociation.first.charge();
	h_muMatched_GenPt_->Fill(pLneg.pt(),wMC);                	      
	h_muMatched_GenEta_->Fill(pLneg.eta(),wMC);               	      
	h_muMatched_RecoPt_->Fill(negGenMuonAssociation.first.pt(),wMC);   
	h_muMatched_RecoEta_->Fill(negGenMuonAssociation.first.eta(),wMC); 
      }

      h_mu_GenReco_DeltaR_->Fill(negGenMuonAssociation.second,wMC);        
      h_mu_GenReco_DeltaR_vsPt_->Fill(pLneg.pt(),negGenMuonAssociation.second,wMC);   
      h_mu_GenReco_DeltaR_vsEta_->Fill(pLneg.eta(),negGenMuonAssociation.second,wMC);  

    }

    //==========================
    if(nRecoMatchedLeptons==2 && compositeCharge<0){
      if(recoLeptonMomenta_[0].pt() > ptCut  && TMath::Abs(recoLeptonMomenta_[0].eta()) < etaCut && 
	 recoLeptonMomenta_[1].pt() > ptCut  && TMath::Abs(recoLeptonMomenta_[1].eta()) < etaCut ){
	
	math::XYZTLorentzVectorD recoZMomentum = recoLeptonMomenta_[0] + recoLeptonMomenta_[1];
	h_RECO_mass_ZfromMuons_->Fill(recoZMomentum.mass(),wMC);
	
	recoZCand.addDaughter(negGenMuonAssociation.first,"elereco1");
	recoZCand.addDaughter(posGenMuonAssociation.first,"elereco2");	
	AddFourMomenta addp4;
	addp4.set( recoZCand );     

	if( recoZMomentum.mass() > minMassCut &&  recoZMomentum.mass() < maxMassCut){
	  myBools_.is_rec_Zkin_yes_=true;
	} // if reco Z in mass window 
      } // if ok reco muons 
    } // if ok reco Z
   
  } else if (isZEE && myBools_.isRightFlavour()) {

    // indices for b-hadron association
    std::pair<int,double> posGenElectronRecoElectronAssociation = std::make_pair(-1,9999.);
    std::pair<int,double> negGenElectronRecoElectronAssociation = std::make_pair(-1,9999.);
    std::pair<pat::Electron,double> posGenElectronAssociation;
    std::pair<pat::Electron,double> negGenElectronAssociation;
    if((*(electrons.product())).size()>0){ 
      posGenElectronAssociation = std::make_pair((*(electrons.product())).at(0),9999.);
      negGenElectronAssociation = std::make_pair((*(electrons.product())).at(0),9999.);
    }
    
    double posMinDeltaRGenElectron(9999.), negMinDeltaRGenElectron(9999.);
    int p(0),n(0);
    
    //========================== matching for positive electron =============================
    for(edm::View<pat::Electron>::const_iterator electron=electrons->begin(); electron!=electrons->end(); ++electron){
      h_ele_RecoPt_->Fill(electron->pt(),wMC);  
      h_ele_RecoEta_->Fill(electron->eta(),wMC); 
      double dR_tmp=ROOT::Math::VectorUtil::DeltaR(electron->momentum(),pLpos.Vect());
      //if( dR_tmp<dRLeptons_){
      if( dR_tmp<posMinDeltaRGenElectron){
	posMinDeltaRGenElectron = dR_tmp;
	posGenElectronRecoElectronAssociation.first  = p;
	posGenElectronRecoElectronAssociation.second = posMinDeltaRGenElectron;
	
	posGenElectronAssociation.first  = (*electron);
	posGenElectronAssociation.second = posMinDeltaRGenElectron;
      } // if gen-reco are matched
      // } // loop on reco electrons
      p++;
    }

    if( posGenElectronRecoElectronAssociation.first != -1 ){

      if(posGenElectronAssociation.second<dRLeptons_){
	recoLeptonMomenta_.push_back( posGenElectronAssociation.first.p4() );
	nRecoMatchedLeptons++;
	compositeCharge*=posGenElectronAssociation.first.charge();
	h_eleMatched_GenPt_->Fill(pLpos.pt(),wMC);                
	h_eleMatched_GenEta_->Fill(pLpos.eta(),wMC);               
	h_eleMatched_RecoPt_->Fill(posGenElectronAssociation.first.pt(),wMC);               
	h_eleMatched_RecoEta_->Fill(posGenElectronAssociation.first.eta(),wMC);   
      }      

      h_ele_GenReco_DeltaR_->Fill(posGenElectronAssociation.second,wMC);        
      h_ele_GenReco_DeltaR_vsPt_->Fill(pLpos.pt(),posGenElectronAssociation.second,wMC);   
      h_ele_GenReco_DeltaR_vsEta_->Fill(pLpos.eta(),posGenElectronAssociation.second,wMC);  
                 
    }

    //========================== matching for negative electron =============================
    for(edm::View<pat::Electron>::const_iterator electron=electrons->begin(); electron!=electrons->end(); ++electron){
      double dR_tmp=ROOT::Math::VectorUtil::DeltaR(electron->momentum(),pLneg.Vect());
      //   if( dR_tmp<dRLeptons_ ){
      if( dR_tmp<negMinDeltaRGenElectron ){
	negMinDeltaRGenElectron = dR_tmp;
	negGenElectronRecoElectronAssociation.first  = n;
	negGenElectronRecoElectronAssociation.second = negMinDeltaRGenElectron;
	
	negGenElectronAssociation.first  = (*electron);
	negGenElectronAssociation.second = negMinDeltaRGenElectron;
      } // if gen-reco are matched
	// } // loop on reco electrons
      n++;
    }

    if( negGenElectronRecoElectronAssociation.first != -1 ){

      if(negGenElectronAssociation.second<dRLeptons_ ){
	recoLeptonMomenta_.push_back( negGenElectronAssociation.first.p4() );
	nRecoMatchedLeptons++;
	compositeCharge*=negGenElectronAssociation.first.charge();
	h_eleMatched_GenPt_->Fill(pLneg.pt(),wMC);                
	h_eleMatched_GenEta_->Fill(pLneg.eta(),wMC);               
	h_eleMatched_RecoPt_->Fill(negGenElectronAssociation.first.pt(),wMC);               
	h_eleMatched_RecoEta_->Fill(negGenElectronAssociation.first.eta(),wMC); 
      }
      
      h_ele_GenReco_DeltaR_->Fill(negGenElectronAssociation.second,wMC);        
      h_ele_GenReco_DeltaR_vsPt_->Fill(pLneg.pt(),negGenElectronAssociation.second,wMC);   
      h_ele_GenReco_DeltaR_vsEta_->Fill(pLneg.eta(),negGenElectronAssociation.second,wMC);  
   
    }
   
    //==========================
    if(nRecoMatchedLeptons==2 && compositeCharge<0){
      if(recoLeptonMomenta_[0].pt() > ptCut  && TMath::Abs(recoLeptonMomenta_[0].eta()) < etaCut && 
	 recoLeptonMomenta_[1].pt() > ptCut  && TMath::Abs(recoLeptonMomenta_[1].eta()) < etaCut &&	 
	 ( abs(negGenElectronAssociation.first.superCluster()->eta())<1.444 || fabs(posGenElectronAssociation.first.superCluster()->eta())>1.566 ) // veto on ECAL crack
	 ){
	
	math::XYZTLorentzVectorD recoZMomentum = recoLeptonMomenta_[0] + recoLeptonMomenta_[1];
	h_RECO_mass_ZfromElectrons_->Fill(recoZMomentum.mass(),wMC);

	recoZCand.addDaughter(negGenElectronAssociation.first,"elereco1");
	recoZCand.addDaughter(posGenElectronAssociation.first,"elereco2");	
	AddFourMomenta addp4;
	addp4.set( recoZCand );     
	
	if( recoZMomentum.mass() > minMassCut &&  recoZMomentum.mass() < maxMassCut){
	  myBools_.is_rec_Zkin_yes_=true;
	} // if reco Z in mass window
      } // if ok reco electrons
    } // if ok reco Z
  } // ends electron case

  // ========================= filter the reco jets =================================

  if(verbose_)std::cout<<"before reco jet matching: recoLeptonMomenta_.size()="<<recoLeptonMomenta_.size()<<" nRecoMatchedLeptons:"<<nRecoMatchedLeptons<<" compositeCharge: "<<compositeCharge<<std::endl;  
  std::vector<const pat::Jet*> the_filteredRecoJets; 
  std::vector<const pat::Jet*> all_the_RecoJets; 

  int nRecoJets(0);    
  int nOkRecoJets(0);  

  if(verbose_) std::cout<<"is_rec_Zkin_yes!"<<std::endl;
  for(edm::View<pat::Jet>::const_iterator recojet_it=recoJets.begin(); recojet_it!=recoJets.end(); ++recojet_it){ 
    all_the_RecoJets.push_back(&(*recojet_it));
    nRecoJets++;
    if(recojet_it->pt()>jetPtCut && TMath::Abs(recojet_it->eta())<jetEtaCut) {

      if(recoLeptonMomenta_.size()>1){
	if( ROOT::Math::VectorUtil::DeltaR(recojet_it->momentum(), recoLeptonMomenta_[0].Vect() )>dRJets_ && 
	    ROOT::Math::VectorUtil::DeltaR(recojet_it->momentum(), recoLeptonMomenta_[1].Vect() )>dRJets_ ){
	  
	  the_filteredRecoJets.push_back(&(*recojet_it));
	  nOkRecoJets++;
	}  
	//     } else { // if at least 2 matched leptons 
	// 	the_filteredRecoJets.push_back(&(*recojet_it));
	// 	nOkRecoJets++;
      }
    } // ends if reco jet is ok
  } // for loop on reco jets 

  myBools_.nRecoJets_=nRecoJets;
  myBools_.nOkRecoJets_=nOkRecoJets;

  //================================================================
  // matching of any-flavour genJets
  std::vector<std::pair<const pat::Jet*,const reco::GenJet*> > recAnyFlavJGenAnyFlavJMatch = matchJetsByDr(the_filteredRecoJets,the_any_flavour_genJets,dRJets_,true,true);
  if (recAnyFlavJGenAnyFlavJMatch.size() > 0) {
    myBools_.is_rec_j_yes_=true;
  }

  if(the_b_flavour_genJets.size()!=0 && verbose_){
    std::cout<<"b genJets: "<<the_b_flavour_genJets.size()<<" filtered b genJets: "<< nOkBFlavouredJets<<" recoJets: "<<recoJets.size()<<" filtered recoJets:"<<the_filteredRecoJets.size()<<std::endl;
  }

  // ========================= reco jets matching ====================================

  std::vector<std::pair<const pat::Jet*,const reco::GenJet*> > filtRecJGenJMatch = matchJetsByDr(the_filteredRecoJets,the_b_flavour_genJets,dRJets_,true,true);
  if (filtRecJGenJMatch.size() > 0 ) {
    myBools_.is_rec_b_yes_=true;
    if(verbose_) std::cout<<"recJGenJMatch.size(): "<<filtRecJGenJMatch.size()<<" deltaR(j,j):"<<ROOT::Math::VectorUtil::DeltaR((*filtRecJGenJMatch[0].first).momentum(),(*filtRecJGenJMatch[0].second).momentum())<<std::endl;
  }

  std::vector<const pat::Jet*> filtVecRecJets;
  filtVecRecJets.resize(filtRecJGenJMatch.size());
  std::vector<const reco::GenJet*> filtVecHadJets;
  filtVecHadJets.resize(filtRecJGenJMatch.size());
  
  for (unsigned i = 0; i < filtVecHadJets.size(); ++i){ 
    
    h_bjet_GenMatchedPt_->Fill((filtRecJGenJMatch[i].second)->pt(),wMC);   
    h_bjet_GenMatchedEta_->Fill((filtRecJGenJMatch[i].second)->eta(),wMC);  
    h_bjet_GenMatchedMass_->Fill((filtRecJGenJMatch[i].second)->mass(),wMC); 

    h_selbjet_RecoMatchedPt_->Fill((filtRecJGenJMatch[i].first)->pt(),wMC);  
    h_selbjet_RecoMatchedEta_->Fill((filtRecJGenJMatch[i].first)->eta(),wMC); 
    h_selbjet_RecoMatchedMass_->Fill((filtRecJGenJMatch[i].first)->mass(),wMC);
    
    double dR_tmp= ROOT::Math::VectorUtil::DeltaR((filtRecJGenJMatch[i].second)->momentum(),(filtRecJGenJMatch[i].first)->momentum());
    
    h_selbjet_GenRecoResponse_->Fill((filtRecJGenJMatch[i].first)->pt()/(filtRecJGenJMatch[i].second)->pt(),wMC);
    h_selbjet_GenPtvsRecoPt_->Fill((filtRecJGenJMatch[i].first)->pt(),(filtRecJGenJMatch[i].second)->pt(),wMC);
    h_selbjet_GenReco_DeltaR_->Fill(dR_tmp,wMC);        
    h_selbjet_GenReco_DeltaR_vsPt_->Fill((filtRecJGenJMatch[i].second)->pt(),dR_tmp,wMC);   
    h_selbjet_GenReco_DeltaR_vsEta_->Fill((filtRecJGenJMatch[i].second)->eta(),dR_tmp,wMC);  
    
    filtVecRecJets[i] = (filtRecJGenJMatch[i].first);
    filtVecHadJets[i] = (filtRecJGenJMatch[i].second);
  }

  // ========================= match all reco jets ====================================

  std::vector<std::pair<const pat::Jet*,const reco::GenJet*> > recJGenJMatch = matchJetsByDr(all_the_RecoJets,the_b_flavour_genJets,10.,true,true);
  std::vector<const pat::Jet*> vecRecJets;
  vecRecJets.resize(recJGenJMatch.size());
  std::vector<const reco::GenJet*> vecHadJets;
  vecHadJets.resize(recJGenJMatch.size());
  
  for (unsigned i = 0; i < vecHadJets.size(); ++i){ 
    
    double dR_tmp= ROOT::Math::VectorUtil::DeltaR((recJGenJMatch[i].second)->momentum(),(recJGenJMatch[i].first)->momentum());
    
    h_bjet_RecoMatchedPt_->Fill((recJGenJMatch[i].first)->pt(),wMC);    
    h_bjet_RecoMatchedEta_->Fill((recJGenJMatch[i].first)->eta(),wMC);  
    h_bjet_RecoMatchedMass_->Fill((recJGenJMatch[i].first)->mass(),wMC);

    h_bjet_GenRecoResponse_->Fill(((recJGenJMatch[i].first)->pt())/((recJGenJMatch[i].second)->pt()),wMC);
    h_bjet_GenPtvsRecoPt_->Fill((recJGenJMatch[i].first)->pt(),(recJGenJMatch[i].second)->pt(),wMC);
    h_bjet_GenReco_DeltaR_->Fill(dR_tmp,wMC);        
    h_bjet_GenReco_DeltaR_vsPt_->Fill((recJGenJMatch[i].second)->pt(),dR_tmp,wMC);   
    h_bjet_GenReco_DeltaR_vsEta_->Fill((recJGenJMatch[i].second)->eta(),dR_tmp,wMC);  
    
    vecRecJets[i] = (recJGenJMatch[i].first);
    vecHadJets[i] = (recJGenJMatch[i].second);
  }

  // ========================= match all reco jets ====================================

  //  for(edm::View<pat::Jet>::const_iterator recojet_it=recoJets.begin(); recojet_it!=recoJets.end(); ++recojet_it){  
  //     for(std::vector<const reco::GenJet*>::const_iterator genjet_it=the_b_flavour_genJets.begin(); genjet_it!=the_b_flavour_genJets.end(); ++genjet_it){
  //       double dR_tmp= ROOT::Math::VectorUtil::DeltaR((*recojet_it)->momentum(),(*genjet_it)->momentum());
  //       if( dR_tmp < dRJets_  && myBools_.isRightFlavour() ) {
  //       } // if reco jet matched to gen jet
  //     } // loop on gen b-flavoured jets 
  //   } // if ok reco jet
  
  if(verbose_) std::cout<<"before access to struct "<<std::endl;
  
  // ========================= Counters for acceptance and C_had ======================== 

  // is Gen Z ok
  if(myBools_.isGenZyes()){
    N_gen_Z_yes++;
    sumW_gen_Z_yes+=wMC;
    sumW2_gen_Z_yes+=(wMC*wMC);
  } 
 
  // is Gen B ok + gen Z ok
  if(myBools_.isGenbyes()){
    N_gen_b_yes++;
    sumW_gen_b_yes+=wMC;
    sumW2_gen_b_yes+=(wMC*wMC);

    if(myBools_.isGenZyes()){
      Al_Den++;
      sumW_Al_Den+=wMC;
      sumW2_Al_Den+=(wMC*wMC);
    }
  }

  // is Gen Z kin ok + gen b ok
  if(myBools_.isGen_Zkin_yes()){
    N_gen_Zkin_yes++;
    sumW_gen_Zkin_yes+=wMC;
    sumW2_gen_Zkin_yes+=(wMC*wMC);
    
    if(myBools_.isGenbyes()){
      Al_Num++;
      sumW_Al_Num+=wMC;
      sumW2_Al_Num+=(wMC*wMC);
    }
  }
  
  // is Rec Z ok 
  if(myBools_.isRec_Zkin_yes()){
    N_rec_Zkin_yes++;
    sumW_rec_Zkin_yes+=wMC;
    sumW2_rec_Zkin_yes+=(wMC*wMC);
  }

  // is Rec b ok + Z kin ok
  if(myBools_.isRecb_yes()){
    N_rec_b_yes++;
    sumW_rec_b_yes+=wMC;
    sumW2_rec_b_yes+=(wMC*wMC);

    if(myBools_.isRec_Zkin_yes()) {
      CHad_Num++;
      sumW_CHad_Num+=wMC;
      sumW2_CHad_Num+=(wMC*wMC);
    }
  }
  
  // ========================= reco level corrections  ======================== 
   
  // ------------------------- epsilon_l
  Bool_t rec_yes_alljets = (myBools_.isRecjet_yes() && myBools_.isRec_Zkin_yes());
  Bool_t rec_yes = (myBools_.isRecb_yes() && myBools_.isRec_Zkin_yes());
  if(rec_yes){
    
    Epsilon_l_Den++;
    sumW_Epsilon_l_Den+=wMC;
    sumW2_Epsilon_l_Den+=(wMC*wMC);

    if (isZMM){
      if(isTightZCandidate(recoZCand,beamSpot,true,lCuts_)){
	myBools_.is_rec_lep_idiso_yes_=true;
	Epsilon_l_Num++;
	wLepIdIso*=getMuonScaleFactor(recoZCand);
	myBools_.lepIdIsoWeight_=wLepIdIso;
	myBools_.lepIdIsoWeightSq_=wLepIdIso*wLepIdIso;
	sumW_Epsilon_l_Num+=(wMC*wLepIdIso);
	sumW2_Epsilon_l_Num+=(wMC*wMC*wLepIdIso*wLepIdIso);
      }
    } else if (isZEE){
      if(isTightZCandidate(recoZCand,beamSpot,false,lCuts_)){
	myBools_.is_rec_lep_idiso_yes_=true;
	Epsilon_l_Num++;
	wLepIdIso*=getElectronScaleFactor(recoZCand);
	myBools_.lepIdIsoWeight_=wLepIdIso;
	myBools_.lepIdIsoWeightSq_=wLepIdIso*wLepIdIso;
	sumW_Epsilon_l_Num+=(wMC*wLepIdIso);
	sumW2_Epsilon_l_Num+=(wMC*wMC*wLepIdIso*wLepIdIso);
      }
    }
  }

  // all jets sample
  if(rec_yes_alljets){
    
    Epsilon_l_Den_alljets++;
    sumW_Epsilon_l_Den_alljets+=wMC;
    sumW2_Epsilon_l_Den_alljets+=(wMC*wMC);

    if (isZMM){
      if(isTightZCandidate(recoZCand,beamSpot,true,lCuts_)){
	myBools_.is_rec_lep_idiso_yes_=true;
	Epsilon_l_Num_alljets++;
	wLepIdIso*=getMuonScaleFactor(recoZCand);
	myBools_.lepIdIsoWeight_=wLepIdIso;
	myBools_.lepIdIsoWeightSq_=wLepIdIso*wLepIdIso;
	sumW_Epsilon_l_Num_alljets+=(wMC*wLepIdIso);
	sumW2_Epsilon_l_Num_alljets+=(wMC*wMC*wLepIdIso*wLepIdIso);
      }
    } else if (isZEE){
      if(isTightZCandidate(recoZCand,beamSpot,false,lCuts_)){
	myBools_.is_rec_lep_idiso_yes_=true;
	Epsilon_l_Num_alljets++;
	wLepIdIso*=getElectronScaleFactor(recoZCand);
	myBools_.lepIdIsoWeight_=wLepIdIso;
	myBools_.lepIdIsoWeightSq_=wLepIdIso*wLepIdIso;
	sumW_Epsilon_l_Num_alljets+=(wMC*wLepIdIso);
	sumW2_Epsilon_l_Num_alljets+=(wMC*wMC*wLepIdIso*wLepIdIso);
      }
    }
  }

  // ------------------------- epsilon_b
  
  std::vector<pat::Jet> theBtaggedJetsHE;
  std::vector<pat::Jet> theBtaggedJetsHP;  

  if(myBools_.isRecLepIdIso_yes()){
    
    Epsilon_bHE_Den++;
    sumW_Epsilon_bHE_Den+=(wMC*wLepIdIso);
    sumW2_Epsilon_bHE_Den+=(wMC*wMC*wLepIdIso*wLepIdIso);

    Epsilon_bHP_Den++;
    sumW_Epsilon_bHP_Den+=(wMC*wLepIdIso);
    sumW2_Epsilon_bHP_Den+=(wMC*wMC*wLepIdIso*wLepIdIso);

    for (unsigned i = 0; i < filtVecHadJets.size(); ++i){ 
      if(ZbbUtils::isBJet(*filtVecRecJets[i],"SSVHEM")){
	theBtaggedJetsHE.push_back(*filtVecRecJets[i]);
	myBools_.is_rec_b_HE_yes_=true;
      }
      if(ZbbUtils::isBJet(*filtVecRecJets[i],"SSVHPT")){
	theBtaggedJetsHP.push_back(*filtVecRecJets[i]);
	myBools_.is_rec_b_HP_yes_=true;
      }
    }

    // if SSVHE tag in the event
    if(myBools_. isRecbHE_yes()){
      Epsilon_bHE_Num++;
      wBtagHE=ZbbUtils::getbEffScaleFactorAR(theAssociativeMapEffbMCSSVHE_,theAssociativeMapEffcMCSSVHE_,theBtaggedJetsHE,iSetup,"BTAGSSVHEM","MISTAGSSVHEM",false,1,0,0); 
      myBools_.bTagWeightHE_=wBtagHE;
      myBools_.bTagWeightHESq_=wBtagHE*wBtagHE;
      sumW_Epsilon_bHE_Num+=(wMC*wBtagHE*wLepIdIso);
      sumW2_Epsilon_bHE_Num+=(wMC*wMC*wLepIdIso*wLepIdIso*wBtagHE*wBtagHE);
    } 

    // if SSVHP tag in the event
    if(myBools_. isRecbHP_yes()){
      Epsilon_bHP_Num++;
      wBtagHP=ZbbUtils::getbEffScaleFactorAR(theAssociativeMapEffbMCSSVHP_,theAssociativeMapEffcMCSSVHP_,theBtaggedJetsHP,iSetup,"BTAGSSVHPT","MISTAGSSVHPT",false,1,0,0); 
      myBools_.bTagWeightHP_=wBtagHP;
      myBools_.bTagWeightHPSq_=wBtagHP*wBtagHP;
      sumW_Epsilon_bHP_Num+=(wMC*wBtagHP*wLepIdIso);
      sumW2_Epsilon_bHP_Num+=(wMC*wMC*wLepIdIso*wLepIdIso*wBtagHP*wBtagHP);
    }
  }

  
  Bool_t gen_yes = (myBools_.isGenZyes()  && myBools_.isGenbyes());

  //We now have everything we need to calculate c_hadron!
  if((rec_yes && !gen_yes)) {
    N_rec_y_gen_n++;
    sumW_rec_y_gen_n+=wMC;
    sumW2_rec_y_gen_n+=(wMC*wMC);
  }
  if((rec_yes && gen_yes)) {
    N_rec_y_gen_y++;
    sumW_rec_y_gen_y+=wMC;
    sumW2_rec_y_gen_y+=(wMC*wMC);
  }
  if((!rec_yes && !gen_yes)) {
    N_rec_n_gen_n++;
    sumW_rec_n_gen_n+=wMC;
    sumW2_rec_n_gen_n+=(wMC*wMC);
  }
  if((!rec_yes && gen_yes)){
    N_rec_n_gen_y++;
    sumW_rec_n_gen_y+=wMC;
    sumW2_rec_n_gen_y+=(wMC*wMC);
  }

  //========================= Filling the matrices ====================
  
  if(myBools_.isRightFlavour()){

    // filling the gen matrix
    if(myBools_.isGenZyes()){
      if(myBools_.isGenbyes()){
	h_GENMatrix_->Fill(1.,1.,wMC);
      } else {
	h_GENMatrix_->Fill(1.,0.,wMC);
      }
    } else {
      if(myBools_.isGenbyes()){
	h_GENMatrix_->Fill(0.,1.,wMC);
      } else {
	h_GENMatrix_->Fill(0.,0.,wMC);
      }
    } 

    // filling the gen kin matrix
    if(myBools_.isGen_Zkin_yes()){
      if(myBools_.isGenbyes()){
	h_GENKinMatrix_->Fill(1.,1.,wMC);
      } else {
	h_GENKinMatrix_->Fill(1.,0.,wMC);
      }
    } else {
      if(myBools_.isGenbyes()){
	h_GENKinMatrix_->Fill(0.,1.,wMC);
      } else {
	h_GENKinMatrix_->Fill(0.,0.,wMC);
      }
    } 

    // filling the reco matrix
    if(myBools_.isRec_Zkin_yes()){
      if(myBools_.isRecb_yes()){
	h_RECOMatrix_->Fill(1.,1.,wMC);
      } else {
	h_RECOMatrix_->Fill(1.,0.,wMC);
      }
    } else {
      if(myBools_.isRecb_yes()){
	h_RECOMatrix_->Fill(0.,1.,wMC);
      } else {
	h_RECOMatrix_->Fill(0.,0.,wMC);
      }
    } 
  }

  if(saveNTuple_){   
    myNTuple_->Fill();
  }
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
ZbbAcceptanceCalculator::beginJob()
{
  // from https://twiki.cern.ch/twiki/bin/view/CMS/VectorBosonPlusHeavyFlavor#Pile_up_reweighting
  float data_2011_to_run173692_v3[36] ={
    1.39818e+07,
    5.70465e+07,
    1.35101e+08,
    2.2739e+08,
    3.0301e+08,
    3.39416e+08,
    3.32104e+08,
    2.91552e+08,
    2.34186e+08,
    1.74644e+08,
    1.2225e+08,
    8.0982e+07,
    5.10706e+07,
    3.07955e+07,
    1.7812e+07,
    9.90513e+06,
    5.3052e+06,
    2.74071e+06,
    1.36736e+06,
    659563,
    307928,
    139290,
    61111.2,
    26031.4,
    10776.8,
    4340.47,
    1702.41,
    650.858,
    242.78,
    88.4374,
    31.4874,
    10.9669,
    3.73959,
    1.24937,
    0.600761,
    0
  };
  
  // from https://twiki.cern.ch/twiki/bin/view/CMS/VectorBosonPlusHeavyFlavor#Pile_up_reweighting
  float MC_summer11_PUS4[36] ={
    0.145346,   
    0.0642802,		
    0.0695255,	
    0.0696747,
    0.0692955,
    0.0684997,
    0.0669528,
    0.0645515,
    0.0609865,
    0.0563323,
    0.0507322,
    0.0444681,
    0.0379205,
    0.0315131,
    0.025422,
    0.0200184,
    0.0153776,
    0.0115387,
    0.00847608,
    0.00608715,
    0.00428255,
    0.00297185,
    0.00201918,
    0.0013449,
    0.000881587,
    0.000569954,
    0.000361493,
    0.000228692,
    0.000140791,
    8.44606e-05,
    5.10204e-05,
    3.07802e-05,
    1.81401e-05,
    1.00201e-05,
    5.80004e-06,
    0
  };
  
  std::vector<float> dataPU(36);
  std::vector<float> MCPU(36);
  // for EPS analysis from ZZ analysis
  // copy(data_EPS11, data_EPS11+25, dataPU.begin());
  // copy(PoissonOneXDist, PoissonOneXDist+25, MCPU.begin());
  // for 2/fb analysis from https://twiki.cern.ch/twiki/bin/view/CMS/VectorBosonPlusHeavyFlavor#Pile_up_reweighting
  copy(data_2011_to_run173692_v3,data_2011_to_run173692_v3+36,dataPU.begin());
  copy(MC_summer11_PUS4,MC_summer11_PUS4+36,MCPU.begin());
  theLumiWeights_ = new edm::LumiReWeighting(MCPU, dataPU);
  
  // ==================== Summary Text File =========================

  const char* outfile_name="";
  if(isZEE){
    outfile_name="CorrectionFactors_Zee.txt";
  } else if (isZMM){
    outfile_name="CorrectionFactors_Zmm.txt";
  }

  outfile_.open (outfile_name);

  // ==================== Summary Histograms =========================

  edm::Service<TFileService> fs;
  TH1F::SetDefaultSumw2(kTRUE);

  TString particle_type_;

  if ( partonLevel_ ){
    particle_type_="parton";
  } else {
    particle_type_="hadron"; 
  }

  h1UnwCategories_ =fs->make<TH1F>("h1UnwCategories","Events for each category; Category",10,-0.5,9.5);
  h1WgtCategories_ =fs->make<TH1F>("h1WgtCategories","Events for each category; Category",10,-0.5,9.5);

  h_GEN_mass_ZfromMuons_        = fs->make<TH1F>("GEN_mass_ZfromMuons","mass of #mu#mu (genParticles); GEN M_{#mu^{+}#mu^{-}} (GeV)",100.,0.,150.);
  h_GEN_mass_ZfromElectrons_    = fs->make<TH1F>("GEN_mass_ZfromElectrons","mass of #it{ee} (genParticles); GEN M_{e^{+}e^{-}} (GeV)",100.,0.,150.);

  h_GEN_nBflavouredJets_        = fs->make<TH1F>("GEN_nBflavouredJets","number of b flav genJets; n_{b-jets}",10,-0.5,9.5);
  h_GEN_nOkBflavouredJets_      = fs->make<TH1F>("GEN_nOkBflavouredJets","number of b flav genJets; n^{sel}_{b-jets}",10,-0.5,9.5);

  h_GEN_GenbGenJetMinDeltaR_      = fs->make<TH1F>("GEN_GenbGenJetMinDeltaR","#DeltaR(B,j); #DeltaR(B,j)",100,0.,1.);
  h_GEN_GenbGenJetMinDeltaRVsEta_ = fs->make<TH2F>("GEN_GenbGenJetMinDeltaR_vsJetEta","#DeltaR(B,j) vs jet #eta;"+particle_type_+" Jet #eta; #DeltaR(B,j)",60,-3.,3.,100,0.,1.);
  h_GEN_GenbGenJetMinDeltaRVsPt_  = fs->make<TH2F>("GEN_GenbGenJetMinDeltaR_vsJetPt","#DeltaR(B,j) vs jet p_{T};"+particle_type_+" Jet p_{T} (GeV); #DeltaR(B,j)",100,0.,100.,100,0.,1.);
  h_bparticlePt_                = fs->make<TH1F>("GEN_bparticle_Pt","b-particle p_{T};b-"+particle_type_+" p_{T} (GeV)",100,0.,100.);
  h_bparticleEta_               = fs->make<TH1F>("GEN_bparticle_Eta","b-particle #eta;b-"+particle_type_+" #eta",100,-5.,5.);
  h_bparticleEnergyFraction_    = fs->make<TH1F>("GEN_bparticle_EF","fraction of "+particle_type_+" jet p_{T} carried by b-particle #eta; p_{T}(b-particle)/p_{T} "+particle_type_+"jet",100,0.,1.);

  h_bjet_GenPt_                 = fs->make<TH1F>("GEN_bjet_Pt",particle_type_+" Jet p_{T};"+particle_type_+" Jet p_{T} (GeV)",100,0.,100.);
  h_bjet_GenEta_                = fs->make<TH1F>("GEN_bjet_Eta",particle_type_+" Jet #eta;"+particle_type_+" Jet #eta",100,-5.,5.);
  h_bjet_GenMass_               = fs->make<TH1F>("GEN_bjet_Mass",particle_type_+" Jet mass;"+particle_type_+" Jet mass (GeV)",100,0.,40.);

  h_bjet_GenMatchedPt_          = fs->make<TH1F>("GENMatched_bjet_Pt","reco matched "+particle_type_+" Jet p_{T};reco matched "+particle_type_+" Jet p_{T} (GeV)",100,0.,100.);
  h_bjet_GenMatchedEta_         = fs->make<TH1F>("GENMatched_bjet_Eta","reco matched "+particle_type_+" Jet #eta;reco matched "+particle_type_+" Jet #eta",60,-3.,3.);
  h_bjet_GenMatchedMass_        = fs->make<TH1F>("GENMatched_bjet_Mass","reco matched "+particle_type_+" Jet mass;reco matched "+particle_type_+" Jet mass (GeV)",100,0.,40.);
    														  
  h_bjet_RecoMatchedPt_         = fs->make<TH1F>("RECOMatched_bjet_Pt","gen matched reco Jet p_{T};gen-mathed reco Jet p_{T} (GeV)",100,0.,100.);
  h_bjet_RecoMatchedEta_        = fs->make<TH1F>("RECOMatched_bjet_Eta","gen matched reco Jet #eta;gen-matched reco Jet #eta",60,-3.,3.);	  
  h_bjet_RecoMatchedMass_       = fs->make<TH1F>("RECOMatched_bjet_Mass","gen matched reco Jet mass;gen-matched reco Jet mass (GeV)",100,0.,40.);
  
  h_selbjet_RecoMatchedPt_      = fs->make<TH1F>("RECOMatched_selbjet_Pt","gen matched reco Jet p_{T};selected reco Jet p_{T} (GeV)",100,0.,100.);
  h_selbjet_RecoMatchedEta_     = fs->make<TH1F>("RECOMatched_selbjet_Eta","gen matched reco Jet #eta;selected reco Jet #eta",60,-3.,3.);	     
  h_selbjet_RecoMatchedMass_    = fs->make<TH1F>("RECOMatched_selbjet_Mass","gen matched reco Jet mass;selected reco Jet mass (GeV)",100,0.,40.);
  
  h_bjet_GenRecoResponse_       = fs->make<TH1F>("bjetGenRecoResponse","Ratio of p_{T} of genJet/ p_{T} of recoJet; reco Jet p_{T} / "+particle_type_+" Jet p_{T};jets",40,0.,2.);
  h_bjet_GenPtvsRecoPt_         = fs->make<TH2F>("GEN_bjet_GenPtvsRecoPt",particle_type_+" Jet p_{T} (GeV); reco Jet p_{T} (GeV);"+particle_type_+" Jet p_{T} (GeV)",100,0.,100.,100,0.,100);
  h_bjet_GenReco_DeltaR_        = fs->make<TH1F>("GEN_bjet_GenReco_DeltaR","#DeltaR(bjet_{gen},bjet_{reco}) ; #DeltaR(bjet_{gen},bjet_{reco})",100,0.,1.);
  h_bjet_GenReco_DeltaR_vsPt_   = fs->make<TH2F>("GEN_bjet_GenReco_DeltaR_vsPt","#DeltaR(bjet_{gen},bjet_{reco}) vs bjet p_{T};gen bjet p_{T} (GeV);#DeltaR(bjet_{gen},bjet_{reco})",200,0.,200.,100,0.,1.);
  h_bjet_GenReco_DeltaR_vsEta_  = fs->make<TH2F>("GEN_bjet_GenReco_DeltaR_vsEta","#DeltaR(bjet_{gen},bjet_{reco}) vs bjet #eta;gen  bjet #eta (GeV);#DeltaR(bjet_{gen},bjet_{reco})",100,-3.,3.,100,0.,1.);

  h_selbjet_GenRecoResponse_       = fs->make<TH1F>("selbjetGenRecoResponse","Ratio of p_{T} of genJet/ p_{T} of recoJet; reco Jet p_{T} / "+particle_type_+" Jet p_{T};jets",40,0.,2.);
  h_selbjet_GenPtvsRecoPt_         = fs->make<TH2F>("GEN_selbjet_GenPtvsRecoPt",particle_type_+" Jet p_{T} (GeV); reco Jet p_{T} (GeV);"+particle_type_+" Jet p_{T} (GeV)",100,0.,100.,100,0.,100);
  h_selbjet_GenReco_DeltaR_        = fs->make<TH1F>("GEN_selbjet_GenReco_DeltaR","#DeltaR(bjet_{gen},bjet_{reco}) ; #DeltaR(bjet_{gen},bjet_{reco})",100,0.,1.);
  h_selbjet_GenReco_DeltaR_vsPt_   = fs->make<TH2F>("GEN_selbjet_GenReco_DeltaR_vsPt","#DeltaR(bjet_{gen},bjet_{reco}) vs bjet p_{T};gen bjet p_{T} (GeV);#DeltaR(bjet_{gen},bjet_{reco})",200,0.,200.,100,0.,1.);
  h_selbjet_GenReco_DeltaR_vsEta_  = fs->make<TH2F>("GEN_selbjet_GenReco_DeltaR_vsEta","#DeltaR(bjet_{gen},bjet_{reco}) vs bjet #eta;gen  bjet #eta (GeV);#DeltaR(bjet_{gen},bjet_{reco})",100,-3.,3.,100,0.,1.);

  h_mu_GenReco_DeltaR_          = fs->make<TH1F>("GEN_muon_GenReco_DeltaR","#DeltaR(#mu_{gen},#mu_{reco}) ; #DeltaR(#mu_{gen},#mu_{reco})",100,0.,0.5);
  h_mu_GenReco_DeltaR_vsPt_     = fs->make<TH2F>("GEN_muon_GenReco_DeltaR_vsPt","#DeltaR(#mu_{gen},#mu_{reco}) vs muon p_{T}; gen muon p_{T} (GeV);#DeltaR(#mu_{gen},#mu_{reco})",200,0.,200.,100,0.,0.5);
  h_mu_GenReco_DeltaR_vsEta_    = fs->make<TH2F>("GEN_muon_GenReco_DeltaR_vsEta","#DeltaR(#mu_{gen},#mu_{reco}) vs muon #eta; gen muon #eta (GeV);#DeltaR(#mu_{gen},#mu_{reco})",100,-3.,3.,100,0.,0.5);
  
  h_ele_GenReco_DeltaR_         = fs->make<TH1F>("GEN_electron_GenReco_DeltaR","#DeltaR(e_{gen},e_{reco}); #DeltaR(e_{gen},e_{reco})",100,0.,0.5);
  h_ele_GenReco_DeltaR_vsPt_    = fs->make<TH2F>("GEN_electron_GenReco_DeltaR_vsPt","#DeltaR(e_{gen},e_{reco}) vs electron p_{T}; gen electron p_{T} (GeV);#DeltaR(e_{gen},e_{reco})",200,0.,200.,100,0.,0.5);
  h_ele_GenReco_DeltaR_vsEta_   = fs->make<TH2F>("GEN_electron_GenReco_DeltaR_vsEta","#DeltaR(e_{gen},e_{reco}) vs electron #eta; gen electron #eta (GeV);#DeltaR(e_{gen},e_{reco})",100,-3.,3.,100,0.,0.5);

  h_ele_GenPt_                  = fs->make<TH1F>("GEN_electron_Pt","gen Electron p_{T};gen Electron p_{T} (GeV)",100,0.,100.);
  h_ele_GenEta_                 = fs->make<TH1F>("GEN_electron_Eta","gen Electron #eta;gen Electron #eta (GeV)",80,-4.,4.);
  h_ele_RecoPt_                 = fs->make<TH1F>("RECO_electron_Pt","reco Electron p_{T};reco Electron p_{T} (GeV)",100,0.,100.);
  h_ele_RecoEta_                = fs->make<TH1F>("RECO_electron_Eta","reco Electron #eta;reco Electron #eta (GeV)",70,-3.5,3.5);
  
  h_mu_GenPt_                   = fs->make<TH1F>("GEN_muon_Pt","gen Muon p_{T};gen Muon p_{T} (GeV)",100,0.,100.);				     	                             
  h_mu_GenEta_                  = fs->make<TH1F>("GEN_muon_Eta","gen Muon #eta;gen Muon #eta (GeV)",70,-3.5,3.5);	  			   	
  h_mu_RecoPt_                  = fs->make<TH1F>("RECO_muon_Pt","reco Muon p_{T};reco Muon p_{T} (GeV)",100,0.,100.);		       
  h_mu_RecoEta_                 = fs->make<TH1F>("RECO_muon_Eta","reco Muon #eta;reco Muon #eta (GeV)",70,-3.5,3.5);

  h_eleMatched_GenPt_           = fs->make<TH1F>("GENMatched_electron_Pt","gen matched Electron p_{T};gen matched Electron p_{T} (GeV)",100,0.,100.);
  h_eleMatched_GenEta_          = fs->make<TH1F>("GENMatched_electron_Eta","gen matched Electron #eta;gen matched Electron #eta (GeV)",70,-3.5,3.5);
  h_eleMatched_RecoPt_          = fs->make<TH1F>("RECOMatched_electron_Pt","reco matched Electron p_{T};reco matched Electron p_{T} (GeV)",100,0.,100.);
  h_eleMatched_RecoEta_         = fs->make<TH1F>("RECOMatched_electron_Eta","reco matched Electron #eta;reco matched Electron #eta (GeV)",70,-3.5,3.5);
  
  h_muMatched_GenPt_            = fs->make<TH1F>("GENMatched_muon_Pt","gen matched Muon p_{T};gen matched Muon p_{T} (GeV)",100,0.,100.);				     	                             
  h_muMatched_GenEta_           = fs->make<TH1F>("GENMatched_muon_Eta","gen matched Muon #eta;gen matched Muon #eta (GeV)",70,-3.5,3.5);	  			   	
  h_muMatched_RecoPt_           = fs->make<TH1F>("RECOMatched_muon_Pt","reco matched Muon p_{T};reco matched Muon p_{T} (GeV)",100,0.,100.);		       
  h_muMatched_RecoEta_          = fs->make<TH1F>("RECOMatched_muon_Eta","reco matched Muon #eta;reco matched Muon #eta (GeV)",70,-3.5,3.5);

  h_RECO_mass_ZfromMuons_       = fs->make<TH1F>("RECO_mass_ZfromMuons","mass of #mu#mu; RECO M_{#mu^{+}#mu^{-}} (GeV)",100.,0.,150.);
  h_RECO_mass_ZfromElectrons_   = fs->make<TH1F>("RECO_mass_ZfromElectrons","mass of #it{ee}; RECO M_{e^{+}e^{-}} (GeV)",100.,0.,150.);

  TString genMatrixTitle    = "gen Acceptance Matrix";
  TString genKinMatrixTitle = "kin Acceptance Matrix";
  TString recoMatrixTitle   = "reco Acceptance Matrix";
  
  if(isZEE){
    genMatrixTitle    =  genMatrixTitle+" (Z #rightarrow ee)";  
    genKinMatrixTitle =	 genKinMatrixTitle+" (Z #rightarrow ee)"; 
    recoMatrixTitle   =  recoMatrixTitle+" (Z #rightarrow ee)";   
  } else {
    genMatrixTitle    =  genMatrixTitle+" (Z #rightarrow #mu#mu)"; 
    genKinMatrixTitle =	 genKinMatrixTitle+" (Z #rightarrow #mu#mu)"; 
    recoMatrixTitle   =	 recoMatrixTitle+" (Z #rightarrow #mu#mu)";   
  }

  // matrices
  h_GENMatrix_                      = fs->make<TH2F>("GENMatrix",genMatrixTitle,2,-0.5,1.5,2,-0.5,1.5);
  h_GENKinMatrix_                   = fs->make<TH2F>("GENKinMatrix",genKinMatrixTitle,2,-0.5,1.5,2,-0.5,1.5);
  h_RECOMatrix_                     = fs->make<TH2F>("RECOMatrix",recoMatrixTitle,2,-0.5,1.5,2,-0.5,1.5);   

  TString genMatrixXBinLabels[2]    = {"gen Z N","gen Z Y"};
  TString genKinMatrixXBinLabels[2] = {"gen Z kin N","gen Z kin Y"};
  TString genMatrixYBinLabels[2]    = {"gen b N","gen b Y"}; 
  TString recoMatrixXBinLabels[2]   = {"reco Z kin N","reco Z kin Y"};
  TString recoMatrixYBinLabels[2]   = {"reco b N","reco b Y"};

  for(UInt_t bin=1;bin<=2; bin++){
    h_GENMatrix_->GetXaxis()->SetBinLabel(bin,genMatrixXBinLabels[bin-1]); 
    h_GENMatrix_->GetYaxis()->SetBinLabel(bin,genMatrixYBinLabels[bin-1]); 
    h_GENKinMatrix_->GetXaxis()->SetBinLabel(bin,genKinMatrixXBinLabels[bin-1]);  
    h_GENKinMatrix_->GetYaxis()->SetBinLabel(bin,genMatrixYBinLabels[bin-1]);
    h_RECOMatrix_->GetXaxis()->SetBinLabel(bin,recoMatrixXBinLabels[bin-1]);  
    h_RECOMatrix_->GetYaxis()->SetBinLabel(bin,recoMatrixYBinLabels[bin-1]);  
  }

  TString BinLabels[10]  ={"Tot evts","Sel Flav","genZ","genb","genZkin","genZkin+genb","genZ+genb","recZkin","recb","recZkin+recb"};
  
  for(UInt_t bin=1; bin<=10; bin++){
    h1UnwCategories_->GetXaxis()->SetBinLabel(bin,BinLabels[bin-1]); 
    h1WgtCategories_->GetXaxis()->SetBinLabel(bin,BinLabels[bin-1]); 
  }

  // ==================== Summary Ntuple ============================

  if(saveNTuple_){
    myNTuple_ = fs->make<TTree>("T","Data for Acceptance");
    //  myNTuple_->Branch("Acceptance",&myBools_.is_gen_Z_yes,"isgenZyes/O:isgenbyes:isgenZkinyes:isrecZkinyes:isrecbyes");
    myNTuple_->Branch("isRightFlavour",&(myBools_.is_right_flavour_),"isRightFlavour/O");
    myNTuple_->Branch("isgenZyes",&(myBools_.is_gen_Z_yes_),"isgenZyes/O");
    myNTuple_->Branch("isgenbyes",&(myBools_.is_gen_b_yes_),"isgenbyes/O");
    myNTuple_->Branch("isgenjetyes",&(myBools_.is_gen_j_yes_),"isgenjetyes/O");
    myNTuple_->Branch("isgenZkinyes",&(myBools_.is_gen_Zkin_yes_),"isgenZkinyes/O");
    myNTuple_->Branch("isrecZkinyes",&(myBools_.is_rec_Zkin_yes_),"isrecZkinyes/O");
    myNTuple_->Branch("isrecbyes",&(myBools_.is_rec_b_yes_),"isrecbyes/O");
    myNTuple_->Branch("isrecjetyes",&(myBools_.is_rec_j_yes_),"isrecjetyes/O");
    myNTuple_->Branch("isreclepidiso_yes",&(myBools_.is_rec_lep_idiso_yes_),"isreclepidiso_yes/O");
    myNTuple_->Branch("nbgenJets",&(myBools_.nBFlavourJets_),"nbgenJets/I");
    myNTuple_->Branch("nFilteredBgenJets", &(myBools_.nOkBFlavouredJets_),"nFilteredBgenJets/I");
    myNTuple_->Branch("nrecoJets",&(myBools_.nRecoJets_),"nrecoJets/I");
    myNTuple_->Branch("nFilteredRecoJets",&(myBools_.nOkRecoJets_),"nFilteredRecoJets/I"); 
    myNTuple_->Branch("llMass",&(myBools_.mLL_),"llMass/F");
    myNTuple_->Branch("evtweight",&(myBools_.evtWeight_),"evtweight/F");
    myNTuple_->Branch("evtweightSq",&(myBools_.evtWeightSq_),"evtweightSq/F");
    myNTuple_->Branch("lepIdIsoWeight",&(myBools_.lepIdIsoWeight_),"lepIdIsoWeight/F");
    myNTuple_->Branch("lepIdIsoWeightSq",&(myBools_.lepIdIsoWeightSq_),"lepIdIsoWeightSq/F");
    myNTuple_->Branch("bTagWeightHE",&(myBools_.bTagWeightHE_)," bTagWeightHE/F");  
    myNTuple_->Branch("bTagWeightHESq",&(myBools_.bTagWeightHESq_),"bTagWeightHESq/F");
    myNTuple_->Branch("bTagWeightHP",&(myBools_.bTagWeightHP_ )," bTagWeightHP/F");  
    myNTuple_->Branch("bTagWeightHPSq",&(myBools_.bTagWeightHPSq_),"bTagWeightHPSq/F");
    myNTuple_->Branch("isreclepidiso_yes",&(myBools_.is_rec_lep_idiso_yes_),"isreclepidiso_yes/O");
    myNTuple_->Branch("isrecbHE_yes",&(myBools_.is_rec_b_HE_yes_),"isrec_b_HE_yes/O");
    myNTuple_->Branch("isrecbHP_yes",&(myBools_.is_rec_b_HP_yes_),"isrec_b_HP_yes/O");
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZbbAcceptanceCalculator::endJob() 
{

  h1UnwCategories_->SetBinContent(1,N_tot_events);
  h1WgtCategories_->SetBinContent(1,sumW_tot_events);

  h1UnwCategories_->SetBinContent(2,N_right_flavour);
  h1WgtCategories_->SetBinContent(2,sumW_right_flavour);

  h1UnwCategories_->SetBinContent(3,N_gen_Z_yes);
  h1WgtCategories_->SetBinContent(3,sumW_gen_Z_yes);

  h1UnwCategories_->SetBinContent(4,N_gen_b_yes);
  h1WgtCategories_->SetBinContent(4,sumW_gen_b_yes);

  h1UnwCategories_->SetBinContent(5,N_gen_Zkin_yes);
  h1WgtCategories_->SetBinContent(5,sumW_gen_Zkin_yes);

  h1UnwCategories_->SetBinContent(6,Al_Num);
  h1WgtCategories_->SetBinContent(6,sumW_Al_Num);

  h1UnwCategories_->SetBinContent(7,Al_Den);
  h1WgtCategories_->SetBinContent(7,sumW_Al_Den);

  h1UnwCategories_->SetBinContent(8,N_rec_Zkin_yes);
  h1WgtCategories_->SetBinContent(8,sumW_rec_Zkin_yes);

  h1UnwCategories_->SetBinContent(9,N_rec_b_yes);
  h1WgtCategories_->SetBinContent(9,sumW_rec_b_yes);

  h1UnwCategories_->SetBinContent(10,CHad_Num);
  h1WgtCategories_->SetBinContent(10,sumW_CHad_Num);

  // Definition of correction factors

  Double_t UnwAl    =  (Al_Num/Al_Den);
  Double_t UnwCHad  =  (CHad_Num/Al_Num);
  Double_t UnwEpsilon_l = (Epsilon_l_Num/Epsilon_l_Den);
  Double_t UnwEpsilon_bHE = (Epsilon_bHE_Num/Epsilon_bHE_Den); 
  Double_t UnwEpsilon_bHP = (Epsilon_bHP_Num/Epsilon_bHP_Den); 

  Double_t Al   = (sumW_Al_Num/sumW_Al_Den); 
  Double_t CHad = (sumW_CHad_Num/sumW_Al_Num);
  Double_t Epsilon_l = (sumW_Epsilon_l_Num/sumW_Epsilon_l_Den);
  Double_t Epsilon_bHE = (sumW_Epsilon_bHE_Num/sumW_Epsilon_bHE_Den); 
  Double_t Epsilon_bHP = (sumW_Epsilon_bHP_Num/sumW_Epsilon_bHP_Den); 


  Double_t ErrUnwAl(-1.), ErrAl(-1.);  
  Double_t ErrUnwCHad(-1.), ErrCHad(-1.);
  Double_t ErrUnwEpsilonl(-1.), ErrEpsilonl(-1.);
  Double_t ErrUnwEpsilonbHE(-1.), ErrEpsilonbHE(-1.);
  Double_t ErrUnwEpsilonbHP(-1.), ErrEpsilonbHP(-1.);

  // trivial error (not taking into account weights)

  //Double_t Al =  sumW_Al_Num/sumW_Al_Den;
  //Double_t CHad =  sumW_CHad_Num/sumW_Al_Num;

  //Double_t ErrAl   = TMath::Sqrt((Al*(1-Al))/sumW_Al_Den);
  //Double_t ErrCHad = TMath::Sqrt((CHad*(1-CHad))/sumW_Al_Num);

  if(useClopperPearsonErrors_){
    
    TEfficiency eff;
    
    ErrUnwAl =   (eff.ClopperPearson(Al_Den+0.5,Al_Num+0.5,0.683 ,1)-UnwAl);
    ErrUnwCHad = (eff.ClopperPearson(Al_Num+0.5,CHad_Num+0.5,0.683 ,1)-UnwCHad);
    ErrUnwEpsilonl = (eff.ClopperPearson(Epsilon_l_Den+0.5,Epsilon_l_Num+0.5,0.683,1)-UnwEpsilon_l);
    ErrUnwEpsilonbHE = (eff.ClopperPearson(Epsilon_bHE_Den+0.5,Epsilon_bHE_Num+0.5,0.683,1)-UnwEpsilon_bHE);
    ErrUnwEpsilonbHP = (eff.ClopperPearson(Epsilon_bHP_Den+0.5,Epsilon_bHP_Num+0.5,0.683,1)-UnwEpsilon_bHP);

    ErrAl =   (eff.ClopperPearson(sumW_Al_Den+0.5,sumW_Al_Num+0.5,0.683 ,1)-Al);
    ErrCHad = (eff.ClopperPearson(sumW_Al_Num+0.5,sumW_CHad_Num+0.5,0.683 ,1)-CHad);
    ErrEpsilonl = (eff.ClopperPearson(sumW_Epsilon_l_Den+0.5,sumW_Epsilon_l_Num+0.5,0.683,1)-Epsilon_l); 
    ErrEpsilonbHE = (eff.ClopperPearson(sumW_Epsilon_bHE_Den+0.5,sumW_Epsilon_bHE_Num+0.5,0.683,1)-Epsilon_bHE);
    ErrEpsilonbHP = (eff.ClopperPearson(sumW_Epsilon_bHP_Den+0.5,sumW_Epsilon_bHP_Num+0.5,0.683,1)-Epsilon_bHP);

  } else{
    
    // simple binomial errors
    
    ErrUnwAl   = TMath::Sqrt((UnwAl*(1-UnwAl))/Al_Den);
    ErrUnwCHad = TMath::Sqrt((UnwCHad*(1-UnwCHad))/Al_Num);
    ErrUnwEpsilonl = TMath::Sqrt((UnwEpsilon_l*(1-UnwEpsilon_l))/Epsilon_l_Den);
    ErrUnwEpsilonbHE = TMath::Sqrt((UnwEpsilon_bHE*(1-UnwEpsilon_bHE))/Epsilon_bHE_Den);
    ErrUnwEpsilonbHP = TMath::Sqrt((UnwEpsilon_bHP*(1-UnwEpsilon_bHP))/Epsilon_bHP_Den);

    // actual formula taking into account the weights in the errors detailed in here:  www.desy.de/~blist/notes/effic.ps.gz
    // errEff = sqrt(sumw2_pass * (sumW_fail)^2 + (sumW2_fail) * (sumW_pass)^2) / (sumW_all)^2 

    ErrAl   =  TMath::Sqrt(sumW2_Al_Num*(sumW_Al_Den-sumW_Al_Num)*(sumW_Al_Den-sumW_Al_Num) + (sumW2_Al_Den - sumW2_Al_Num)*sumW_Al_Num*sumW_Al_Num)/(sumW_Al_Den*sumW_Al_Den);
    ErrCHad =  TMath::Sqrt(sumW2_CHad_Num*(sumW_Al_Num-sumW_CHad_Num)*(sumW_Al_Num-sumW_CHad_Num) + (sumW2_Al_Num - sumW2_CHad_Num)*sumW_CHad_Num*sumW_CHad_Num)/(sumW_Al_Num*sumW_Al_Num);
    ErrEpsilonl = TMath::Sqrt(sumW2_Epsilon_l_Num*(sumW_Epsilon_l_Den-sumW_Epsilon_l_Num)*(sumW_Epsilon_l_Den-sumW_Epsilon_l_Num) + (sumW2_Epsilon_l_Den - sumW2_Epsilon_l_Num)*sumW_Epsilon_l_Num*sumW_Epsilon_l_Num)/(sumW_Epsilon_l_Den*sumW_Epsilon_l_Den);
    ErrEpsilonbHE = TMath::Sqrt(sumW2_Epsilon_bHE_Num*(sumW_Epsilon_bHE_Den-sumW_Epsilon_bHE_Num)*(sumW_Epsilon_bHE_Den-sumW_Epsilon_bHE_Num) + (sumW2_Epsilon_bHE_Den - sumW2_Epsilon_bHE_Num)*sumW_Epsilon_bHE_Num*sumW_Epsilon_bHE_Num)/(sumW_Epsilon_bHE_Den*sumW_Epsilon_bHE_Den);
    ErrEpsilonbHP = TMath::Sqrt(sumW2_Epsilon_bHP_Num*(sumW_Epsilon_bHP_Den-sumW_Epsilon_bHP_Num)*(sumW_Epsilon_bHP_Den-sumW_Epsilon_bHP_Num) + (sumW2_Epsilon_bHP_Den - sumW2_Epsilon_bHP_Num)*sumW_Epsilon_bHP_Num*sumW_Epsilon_bHP_Num)/(sumW_Epsilon_bHP_Den*sumW_Epsilon_bHP_Den);
  }

  // ================================= write on  screen ======================

  std::cout<<"================================================================="<<std::endl;
  std::cout<<"Summary of selected events:  ";
  if( isZMM ) {
    std::cout<<"| Selected Z->mm channel"<<std::endl;
  } else if( isZEE ) {
    std::cout<<"| Selected Z->ee channel"<<std::endl;
  }
  std::cout<<"---------------------------Selection cuts-----------------------"<<std::endl;
  if( isZEE ) {
    std::cout<<"ele pT<" << elePtCut_  <<" ele |eta|< " << eleEtaCut_;
  } else if ( isZMM ){
    std::cout<<"mu pT< " << muonPtCut_ <<" mu  |eta|< " << muonEtaCut_;
  }
  std::cout<<" jet pT<"  << jetPtCut_  <<" jet |eta|< " << jetEtaCut_  <<std::endl;
  std::cout<< minMassCut_ <<" < m(ll)< " << maxMassCut_<<" dR(l,l) <" << dRLeptons_ <<" dR(j,j) < "<< dRJets_ <<std::endl;
  std::cout<<"================================================================="<<std::endl;
  std::cout<<"N(Tot events):                "<<N_tot_events<<std::endl;	  
  std::cout<<"SumW(Tot events):             "<<sumW_tot_events<<std::endl;	  
  std::cout<<"SumW2(Tot events):            "<<sumW2_tot_events<<std::endl;
  std::cout<<"================================================================="<<std::endl;
  std::cout<<"N(Sel Flav):                  "<<N_right_flavour<<"\t | "<< (N_right_flavour/N_tot_events)*100.<<"% of total"<<std::endl;	  
  std::cout<<"SumW(Sel Flav):               "<<sumW_right_flavour<<"\t | "<< (sumW_right_flavour/sumW_tot_events)*100.<<"% of total"<<std::endl;	   
  std::cout<<"SumW2(Sel Flav):              "<<sumW2_right_flavour<<std::endl;
  std::cout<<"================================================================="<<std::endl;  
  std::cout<<"N(Gen Z):                     "<<N_gen_Z_yes<<" \t | "<< (N_gen_Z_yes/N_tot_events)*100.<<"% of total"<<" \t | "<< (N_gen_Z_yes/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;  
  std::cout<<"SumW(Gen Z):                  "<<sumW_gen_Z_yes<<"\t | "<< (sumW_gen_Z_yes/sumW_tot_events)*100.<<"% of total"<<std::endl;
  std::cout<<"SumW2(Gen Z):                 "<<sumW2_gen_Z_yes<<std::endl;	  
  std::cout<<"================================================================="<<std::endl;
  std::cout<<"N(Gen b):                     "<<N_gen_b_yes<<" \t | "<< (N_gen_b_yes/N_tot_events)*100.<<"% of total"<<" \t | "<< (N_gen_b_yes/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;
  std::cout<<"SumW(Gen b):                  "<<sumW_gen_b_yes<<"\t | "<< (sumW_gen_b_yes/sumW_tot_events)*100.<<"% of total"<<std::endl;
  std::cout<<"SumW2(Gen b):                 "<<sumW2_gen_b_yes<<std::endl;	  
  std::cout<<"================================================================="<<std::endl;
  std::cout<<"N(Gen Z kin):                 "<<N_gen_Zkin_yes<<" \t | "<< (N_gen_Zkin_yes/N_tot_events)*100.<<"% of total"<<" \t | "<< (N_gen_Zkin_yes/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;
  std::cout<<"SumW(Gen Z kin):              "<<sumW_gen_Zkin_yes<<"\t | "<< (sumW_gen_Zkin_yes/sumW_tot_events)*100.<<"% of total"<<std::endl;	  
  std::cout<<"SumW2(Gen Z kin):             "<<sumW2_gen_Zkin_yes<<std::endl;	  
  std::cout<<"================================================================="<<std::endl;
  std::cout<<"N(Gen Z kin && Gen b ):       "<<Al_Num<<" \t | "<< (Al_Num/N_tot_events)*100.<<"% of total"<<" \t | "<< (Al_Num/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;	  
  std::cout<<"SumW(Gen Z kin && Gen b):     "<<sumW_Al_Num<<"\t | "<< (sumW_Al_Num/sumW_tot_events)*100.<<"% of total"<<std::endl;
  std::cout<<"SumW2(Gen Z kin && Gen b):    "<<sumW2_Al_Num<<std::endl;	  
  std::cout<<"================================================================="<<std::endl;
  std::cout<<"N(Gen Z && Gen b):            "<<Al_Den<<" \t | "<< (Al_Den/N_tot_events)*100.<<"% of total"<<" \t | "<< (Al_Den/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;
  std::cout<<"Sumw(Gen Z && Gen b):         "<<sumW_Al_Den<<"\t | "<< (sumW_Al_Den/sumW_tot_events)*100.<<"% of total"<<std::endl;
  std::cout<<"Sumw2(Gen Z && Gen b):        "<<sumW2_Al_Den<<std::endl;
  std::cout<<"================================================================="<<std::endl;
  std::cout<<"N(Reco Z kin):                "<<N_rec_Zkin_yes<<" \t | "<< (N_rec_Zkin_yes/N_tot_events)*100.<<"% of total"<<" \t | "<< (N_rec_Zkin_yes/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;  	  
  std::cout<<"Sumw(Reco Z kin):             "<<sumW_rec_Zkin_yes<<"\t | "<< (sumW_rec_Zkin_yes/sumW_tot_events)*100.<<"% of total"<<std::endl; 
  std::cout<<"Sumw2(Reco Z kin):            "<<sumW2_rec_Zkin_yes<<std::endl;
  std::cout<<"================================================================="<<std::endl;
  std::cout<<"N(Reco b):                    "<<N_rec_b_yes<<" \t | "<< (N_rec_b_yes/N_tot_events)*100.<<"% of total"<<" \t | "<< (N_rec_b_yes/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;
  std::cout<<"SumW(Reco b):                 "<<sumW_rec_b_yes<<"\t | "<< (sumW_rec_b_yes/sumW_tot_events)*100.<<"% of total"<<std::endl;
  std::cout<<"SumW2(Reco b):                "<<sumW2_rec_b_yes<<std::endl;
  std::cout<<"================================================================="<<std::endl;
  std::cout<<"N(Reco Z kin && Reco b):      "<<CHad_Num<<" \t | "<< (CHad_Num/N_tot_events)*100.<<"% of total"<<" \t | "<< (CHad_Num/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;
  std::cout<<"SumW(Reco Z kin && Reco b) :  "<<sumW_CHad_Num<<"\t | "<< (sumW_CHad_Num/sumW_tot_events)*100.<<"% of total"<<std::endl;
  std::cout<<"SumW2(Reco Z kin && Reco b):  "<<sumW2_CHad_Num<<std::endl;
  if ( !partonLevel_ ) {
	 std::cout<<"================================================================="<<std::endl;
         std::cout<<"Final acceptance (A_l) & hadron correction factor (C_had) "<<std::endl;
 	 std::cout<<"================================================================="<<std::endl;
         std::cout<<"Unweighted A_l   :          "<< std::setprecision (3) << UnwAl*100.          <<" +/- "<< std::setprecision (2) << ErrUnwAl*100.  <<"%"<<std::endl;
         std::cout<<"Weighted A_l     :          "<< std::setprecision (3) << Al*100.             <<" +/- "<< std::setprecision (2) << ErrAl*100.     <<"%"<<std::endl;
         std::cout<<"Unweighted C_had :          "<< std::setprecision (3) << UnwCHad*100.        <<" +/- "<< std::setprecision (2) << ErrUnwCHad*100.<<"%"<<std::endl;
         std::cout<<"Weighted C_had   :          "<< std::setprecision (3) << CHad*100.           <<" +/- "<< std::setprecision (2) << ErrCHad*100.   <<"%"<<std::endl;
	 std::cout<<"Unweighted Epsilon_l :      "<< std::setprecision (3) << UnwEpsilon_l*100.   <<" +/- "<< std::setprecision (2) << ErrUnwEpsilonl*100.<<"%"<<std::endl;
         std::cout<<"Weighted Epsilon_l   :      "<< std::setprecision (3) << Epsilon_l*100.      <<" +/- "<< std::setprecision (2) << ErrEpsilonl*100.   <<"%"<<std::endl; 
	 std::cout<<"Unweighted Epsilon_b (HE) : "<< std::setprecision (3) << UnwEpsilon_bHE*100. <<" +/- "<< std::setprecision (2) << ErrUnwEpsilonbHE*100.<<"%"<<std::endl;
         std::cout<<"Weighted Epsilon_b (HE)   : "<< std::setprecision (3) << Epsilon_bHE*100.    <<" +/- "<< std::setprecision (2) << ErrEpsilonbHE*100.   <<"%"<<std::endl;
	 std::cout<<"Unweighted Epsilon_b (HP) : "<< std::setprecision (3) << UnwEpsilon_bHP*100. <<" +/- "<< std::setprecision (2) << ErrUnwEpsilonbHP*100.<<"%"<<std::endl;
         std::cout<<"Weighted Epsilon_b (HP)   : "<< std::setprecision (3) << Epsilon_bHP*100.    <<" +/- "<< std::setprecision (2) << ErrEpsilonbHP*100.   <<"%"<<std::endl; 

  } else {
         std::cout<<"================================================================="<<std::endl;
         std::cout<<"Final acceptance (A_l) & parton correction factor (C_part) "<<std::endl;
         std::cout<<"================================================================="<<std::endl;
         std::cout<<"Unweighted A_l    :          "<< std::setprecision (3) << UnwAl*100.          <<" +/- "<< std::setprecision (2) << ErrUnwAl*100.  <<"%"<<std::endl;
         std::cout<<"Weighted A_l      :          "<< std::setprecision (3) << Al*100.             <<" +/- "<< std::setprecision (2) << ErrAl*100.     <<"%"<<std::endl;
         std::cout<<"Unweighted C_part :          "<< std::setprecision (3) << UnwCHad*100.        <<" +/- "<< std::setprecision (2) << ErrUnwCHad*100.<<"%"<<std::endl;
         std::cout<<"Weighted C_part   :          "<< std::setprecision (3) << CHad*100.           <<" +/- "<< std::setprecision (2) << ErrCHad*100.   <<"%"<<std::endl;
	 std::cout<<"Unweighted Epsilon_l :       "<< std::setprecision (3) << UnwEpsilon_l*100.   <<" +/- "<< std::setprecision (2) << ErrUnwEpsilonl*100.<<"%"<<std::endl;
         std::cout<<"Weighted Epsilon_l   :       "<< std::setprecision (3) << Epsilon_l*100.      <<" +/- "<< std::setprecision (2) << ErrEpsilonl*100.   <<"%"<<std::endl;
	 std::cout<<"Unweighted Epsilon_b (HE) :  "<< std::setprecision (3) << UnwEpsilon_bHE*100. <<" +/- "<< std::setprecision (2) << ErrUnwEpsilonbHE*100.<<"%"<<std::endl;
         std::cout<<"Weighted Epsilon_b (HE)   :  "<< std::setprecision (3) << Epsilon_bHE*100.    <<" +/- "<< std::setprecision (2) << ErrEpsilonbHE*100.   <<"%"<<std::endl;
	 std::cout<<"Unweighted Epsilon_b (HP) :  "<< std::setprecision (3) << UnwEpsilon_bHP*100. <<" +/- "<< std::setprecision (2) << ErrUnwEpsilonbHP*100.<<"%"<<std::endl;
         std::cout<<"Weighted Epsilon_b (HP)   :  "<< std::setprecision (3) << Epsilon_bHP*100.    <<" +/- "<< std::setprecision (2) << ErrEpsilonbHP*100.   <<"%"<<std::endl; 
  } 
  std::cout<<"================================================================="<<std::endl;
 

  // ========================== write on the file ================================

  outfile_<<"================================================================="<<std::endl;
  outfile_<<"Summary of selected events:  ";
  if( isZMM ) {
    outfile_<<"| Selected Z->mm channel"<<std::endl;
  } else if( isZEE ) {
    outfile_<<"| Selected Z->ee channel"<<std::endl;
  }
  outfile_<<"---------------------------Selection cuts-----------------------"<<std::endl;
  if( isZEE ) {
    outfile_<<"ele pT<" << elePtCut_  <<" ele |eta|< " << eleEtaCut_;
  } else if ( isZMM ){
    outfile_<<"mu pT< " << muonPtCut_ <<" mu  |eta|< " << muonEtaCut_;
  }
  outfile_<<" jet pT<"  << jetPtCut_  <<" jet |eta|< " << jetEtaCut_  <<std::endl;
  outfile_<< minMassCut_ <<" < m(ll)< " << maxMassCut_<<" dR(l,l) <" << dRLeptons_ <<" dR(j,j) < "<< dRJets_ <<std::endl;
  outfile_<<"================================================================="<<std::endl;
  outfile_<<"N(Tot events):                "<<N_tot_events<<std::endl;	  
  outfile_<<"SumW(Tot events):             "<<sumW_tot_events<<std::endl;	  
  outfile_<<"SumW2(Tot events):            "<<sumW2_tot_events<<std::endl;
  outfile_<<"================================================================="<<std::endl;
  outfile_<<"N(Sel Flav):                  "<<N_right_flavour<<"\t | "<< (N_right_flavour/N_tot_events)*100.<<"% of total"<<std::endl;	  
  outfile_<<"SumW(Sel Flav):               "<<sumW_right_flavour<<"\t | "<< (sumW_right_flavour/sumW_tot_events)*100.<<"% of total"<<std::endl;	   
  outfile_<<"SumW2(Sel Flav):              "<<sumW2_right_flavour<<std::endl;
  outfile_<<"================================================================="<<std::endl;  
  outfile_<<"N(Gen Z):                     "<<N_gen_Z_yes<<" \t | "<< (N_gen_Z_yes/N_tot_events)*100.<<"% of total"<<" \t | "<< (N_gen_Z_yes/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;  
  outfile_<<"SumW(Gen Z):                  "<<sumW_gen_Z_yes<<"\t | "<< (sumW_gen_Z_yes/sumW_tot_events)*100.<<"% of total"<<std::endl;
  outfile_<<"SumW2(Gen Z):                 "<<sumW2_gen_Z_yes<<std::endl;	  
  outfile_<<"================================================================="<<std::endl;
  outfile_<<"N(Gen b):                     "<<N_gen_b_yes<<" \t | "<< (N_gen_b_yes/N_tot_events)*100.<<"% of total"<<" \t | "<< (N_gen_b_yes/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;
  outfile_<<"SumW(Gen b):                  "<<sumW_gen_b_yes<<"\t | "<< (sumW_gen_b_yes/sumW_tot_events)*100.<<"% of total"<<std::endl;
  outfile_<<"SumW2(Gen b):                 "<<sumW2_gen_b_yes<<std::endl;	  
  outfile_<<"================================================================="<<std::endl;
  outfile_<<"N(Gen Z kin):                 "<<N_gen_Zkin_yes<<" \t | "<< (N_gen_Zkin_yes/N_tot_events)*100.<<"% of total"<<" \t | "<< (N_gen_Zkin_yes/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;
  outfile_<<"SumW(Gen Z kin):              "<<sumW_gen_Zkin_yes<<"\t | "<< (sumW_gen_Zkin_yes/sumW_tot_events)*100.<<"% of total"<<std::endl;	  
  outfile_<<"SumW2(Gen Z kin):             "<<sumW2_gen_Zkin_yes<<std::endl;	  
  outfile_<<"================================================================="<<std::endl;
  outfile_<<"N(Gen Z kin && Gen b ):       "<<Al_Num<<" \t | "<< (Al_Num/N_tot_events)*100.<<"% of total"<<" \t | "<< (Al_Num/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;	  
  outfile_<<"SumW(Gen Z kin && Gen b):     "<<sumW_Al_Num<<"\t | "<< (sumW_Al_Num/sumW_tot_events)*100.<<"% of total"<<std::endl;
  outfile_<<"SumW2(Gen Z kin && Gen b):    "<<sumW2_Al_Num<<std::endl;	  
  outfile_<<"================================================================="<<std::endl;
  outfile_<<"N(Gen Z && Gen b):            "<<Al_Den<<" \t | "<< (Al_Den/N_tot_events)*100.<<"% of total"<<" \t | "<< (Al_Den/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;
  outfile_<<"Sumw(Gen Z && Gen b):         "<<sumW_Al_Den<<"\t | "<< (sumW_Al_Den/sumW_tot_events)*100.<<"% of total"<<std::endl;
  outfile_<<"Sumw2(Gen Z && Gen b):        "<<sumW2_Al_Den<<std::endl;
  outfile_<<"================================================================="<<std::endl;
  outfile_<<"N(Reco Z kin):                "<<N_rec_Zkin_yes<<" \t | "<< (N_rec_Zkin_yes/N_tot_events)*100.<<"% of total"<<" \t | "<< (N_rec_Zkin_yes/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;  	  
  outfile_<<"Sumw(Reco Z kin):             "<<sumW_rec_Zkin_yes<<"\t | "<< (sumW_rec_Zkin_yes/sumW_tot_events)*100.<<"% of total"<<std::endl; 
  outfile_<<"Sumw2(Reco Z kin):            "<<sumW2_rec_Zkin_yes<<std::endl;
  outfile_<<"================================================================="<<std::endl;
  outfile_<<"N(Reco b):                    "<<N_rec_b_yes<<" \t | "<< (N_rec_b_yes/N_tot_events)*100.<<"% of total"<<" \t | "<< (N_rec_b_yes/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;
  outfile_<<"SumW(Reco b):                 "<<sumW_rec_b_yes<<"\t | "<< (sumW_rec_b_yes/sumW_tot_events)*100.<<"% of total"<<std::endl;
  outfile_<<"SumW2(Reco b):                "<<sumW2_rec_b_yes<<std::endl;
  outfile_<<"================================================================="<<std::endl;
  outfile_<<"N(Reco Z kin && Reco b):      "<<CHad_Num<<" \t | "<< (CHad_Num/N_tot_events)*100.<<"% of total"<<" \t | "<< (CHad_Num/N_right_flavour)*100.<<"% of Sel Flav"<<std::endl;
  outfile_<<"SumW(Reco Z kin && Reco b) :  "<<sumW_CHad_Num<<"\t | "<< (sumW_CHad_Num/sumW_tot_events)*100.<<"% of total"<<std::endl;
  outfile_<<"SumW2(Reco Z kin && Reco b):  "<<sumW2_CHad_Num<<std::endl;
  if ( !partonLevel_ ) {
    outfile_<<"================================================================="<<std::endl;
    outfile_<<"Final acceptance (A_l) & hadron correction factor (C_had) "<<std::endl;
    outfile_<<"================================================================="<<std::endl;
    outfile_<<"Unweighted A_l   :          "<< std::setprecision (3) << UnwAl*100.          <<" +/- "<< std::setprecision (2) << ErrUnwAl*100.  <<"%"<<std::endl;
    outfile_<<"Weighted A_l     :          "<< std::setprecision (3) << Al*100.             <<" +/- "<< std::setprecision (2) << ErrAl*100.     <<"%"<<std::endl;
    outfile_<<"Unweighted C_had :          "<< std::setprecision (3) << UnwCHad*100.        <<" +/- "<< std::setprecision (2) << ErrUnwCHad*100.<<"%"<<std::endl;
    outfile_<<"Weighted C_had   :          "<< std::setprecision (3) << CHad*100.           <<" +/- "<< std::setprecision (2) << ErrCHad*100.   <<"%"<<std::endl;
    outfile_<<"Unweighted Epsilon_l :      "<< std::setprecision (3) << UnwEpsilon_l*100.   <<" +/- "<< std::setprecision (2) << ErrUnwEpsilonl*100.<<"%"<<std::endl;
    outfile_<<"Weighted Epsilon_l   :      "<< std::setprecision (3) << Epsilon_l*100.      <<" +/- "<< std::setprecision (2) << ErrEpsilonl*100.   <<"%"<<std::endl; 
    outfile_<<"Unweighted Epsilon_b (HE) : "<< std::setprecision (3) << UnwEpsilon_bHE*100. <<" +/- "<< std::setprecision (2) << ErrUnwEpsilonbHE*100.<<"%"<<std::endl;
    outfile_<<"Weighted Epsilon_b (HE)   : "<< std::setprecision (3) << Epsilon_bHE*100.    <<" +/- "<< std::setprecision (2) << ErrEpsilonbHE*100.   <<"%"<<std::endl;
    outfile_<<"Unweighted Epsilon_b (HP) : "<< std::setprecision (3) << UnwEpsilon_bHP*100. <<" +/- "<< std::setprecision (2) << ErrUnwEpsilonbHP*100.<<"%"<<std::endl;
    outfile_<<"Weighted Epsilon_b (HP)   : "<< std::setprecision (3) << Epsilon_bHP*100.    <<" +/- "<< std::setprecision (2) << ErrEpsilonbHP*100.   <<"%"<<std::endl; 
  } else {
    outfile_<<"================================================================="<<std::endl;
    outfile_<<"Final acceptance (A_l) & parton correction factor (C_part) "<<std::endl;
    outfile_<<"================================================================="<<std::endl;
    outfile_<<"Unweighted A_l    :          "<< std::setprecision (3) << UnwAl*100.          <<" +/- "<< std::setprecision (2) << ErrUnwAl*100.  <<"%"<<std::endl;
    outfile_<<"Weighted A_l      :          "<< std::setprecision (3) << Al*100.             <<" +/- "<< std::setprecision (2) << ErrAl*100.     <<"%"<<std::endl;
    outfile_<<"Unweighted C_part :          "<< std::setprecision (3) << UnwCHad*100.        <<" +/- "<< std::setprecision (2) << ErrUnwCHad*100.<<"%"<<std::endl;
    outfile_<<"Weighted C_part   :          "<< std::setprecision (3) << CHad*100.           <<" +/- "<< std::setprecision (2) << ErrCHad*100.   <<"%"<<std::endl;
    outfile_<<"Unweighted Epsilon_l :       "<< std::setprecision (3) << UnwEpsilon_l*100.   <<" +/- "<< std::setprecision (2) << ErrUnwEpsilonl*100.<<"%"<<std::endl;
    outfile_<<"Weighted Epsilon_l   :       "<< std::setprecision (3) << Epsilon_l*100.      <<" +/- "<< std::setprecision (2) << ErrEpsilonl*100.   <<"%"<<std::endl;
    outfile_<<"Unweighted Epsilon_b (HE) :  "<< std::setprecision (3) << UnwEpsilon_bHE*100. <<" +/- "<< std::setprecision (2) << ErrUnwEpsilonbHE*100.<<"%"<<std::endl;
    outfile_<<"Weighted Epsilon_b (HE)   :  "<< std::setprecision (3) << Epsilon_bHE*100.    <<" +/- "<< std::setprecision (2) << ErrEpsilonbHE*100.   <<"%"<<std::endl;
    outfile_<<"Unweighted Epsilon_b (HP) :  "<< std::setprecision (3) << UnwEpsilon_bHP*100. <<" +/- "<< std::setprecision (2) << ErrUnwEpsilonbHP*100.<<"%"<<std::endl;
    outfile_<<"Weighted Epsilon_b (HP)   :  "<< std::setprecision (3) << Epsilon_bHP*100.    <<" +/- "<< std::setprecision (2) << ErrEpsilonbHP*100.   <<"%"<<std::endl; 
  } 

  outfile_.close();

  h_GENMatrix_->Scale(100./sumW_right_flavour);
  h_GENKinMatrix_->Scale(100./sumW_right_flavour);
  h_RECOMatrix_->Scale(100./sumW_right_flavour);

}

// ------------ method called when starting to processes a run  ------------
void 
ZbbAcceptanceCalculator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ZbbAcceptanceCalculator::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ZbbAcceptanceCalculator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ZbbAcceptanceCalculator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZbbAcceptanceCalculator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// ----------- method to decide flavour of Z in the event
std::pair<bool,bool> 
ZbbAcceptanceCalculator::WhichFlavour(std::vector<const reco::Candidate*> theGenMuons,std::vector<const reco::Candidate*> theGenElectrons){
  
  std::pair<bool,bool> theRes = std::make_pair(false,false);

  bool isDiMu_(false), isDiEle_(false);

  // checks that there is at least a di-lepton
  if(theGenElectrons.size()>1 || theGenMuons.size()>1){

    if(theGenElectrons.size()>1 && theGenMuons.size()<=1){
      isDiEle_=true;
    } else if(theGenMuons.size()>1 && theGenElectrons.size()<=1){
      isDiMu_=true;
    } else {
      
      Double_t mZ_=91.1876;
      
      std::vector<const reco::Candidate*> sortedGenElectrons = sortLeptonsByPt(theGenElectrons);
      std::vector<const reco::Candidate*> sortedGenMuons     = sortLeptonsByPt(theGenMuons);
      
      Double_t diffmassMu_(0.),diffmassEle_(0.);
      
      diffmassMu_=TMath::Abs((sortedGenMuons[0]->p4()+sortedGenMuons[1]->p4()).M() - mZ_);
      diffmassEle_=TMath::Abs((sortedGenElectrons[0]->p4()+sortedGenElectrons[1]->p4()).M() - mZ_);
      
      if(diffmassMu_<diffmassEle_){
	isDiMu_=true;
      } else {
	isDiEle_=true;
      }
    }
  }
  
  theRes.first  = isDiMu_;
  theRes.second = isDiEle_;

  return theRes;
}

// ------------ method to return the Z leptons from vector of candidates ------------------
std::pair<const reco::Candidate*,const reco::Candidate*> 
ZbbAcceptanceCalculator::getTheZLeptons(std::vector<const reco::Candidate*> theGenMuons,std::vector<const reco::Candidate*> theGenElectrons){
  
  std::pair<const reco::Candidate*,const reco::Candidate*> theZLeptons;
  
  std::vector<const reco::Candidate*>   sortedGenMuons     = sortLeptonsByPt(theGenMuons);
  std::vector<const reco::Candidate*>   sortedGenElectrons = sortLeptonsByPt(theGenElectrons);
  
  if(WhichFlavour(sortedGenMuons,sortedGenElectrons).first){
    theZLeptons.first  = sortedGenMuons[0];
    theZLeptons.second = sortedGenMuons[1];
  } else if(WhichFlavour(sortedGenMuons,sortedGenElectrons).second) {
    theZLeptons.first  = sortedGenElectrons[0];
    theZLeptons.second = sortedGenElectrons[1];
  }

  return theZLeptons;
  
}

// ------------ method to tag b-hadrons ---------------------------------------------------
bool 
ZbbAcceptanceCalculator::hasBottom(const reco::Candidate &c) 
// PDG-like implementation
{
  int code1;
  int code2;
  bool tmpHasBottom = false;
  code1 = (int)( ( abs(c.pdgId() ) / 100)%10 );
  code2 = (int)( ( abs(c.pdgId() ) /1000)%10 );
  if ( code1 == 5 || code2 == 5) tmpHasBottom = true;
  return tmpHasBottom;
}

bool 
ZbbAcceptanceCalculator::hasBottomIC(const reco::Candidate &c) 
// Imperial College implementation
{
  
  bool tmpHasBottom = false;
  bool code1 =  ( abs(c.pdgId())> 500 && abs(c.pdgId())< 600);
  bool code2 =  ( abs(c.pdgId())> 5000 && abs(c.pdgId())< 6000);
  if ( code1 || code2 ) tmpHasBottom = true;
  return tmpHasBottom;
}


// ------------ method to tag c-hadrons ---------------------------------------------------
bool 
ZbbAcceptanceCalculator::hasCharm(const reco::Candidate &c) 
{
  int code1;
  int code2;
  bool tmpHasCharm = false;
  code1 = (int)( ( abs(c.pdgId() ) / 100)%10 );
  code2 = (int)( ( abs(c.pdgId() ) /1000)%10 );
  if ( code1 == 4 || code2 == 4) tmpHasCharm = true;
  return tmpHasCharm;
}

// --------------  sort leptons 4mom by pt ------------------------------------------------

std::vector<math::XYZTLorentzVectorD>
ZbbAcceptanceCalculator::sort4MomentaByPt(std::vector<math::XYZTLorentzVectorD> leptonMomenta_)
{
  // std::cout<<"Sorting"<<std::endl;
  std::vector<math::XYZTLorentzVectorD> sorted4Momenta_ = leptonMomenta_;
  std::vector<double> ptParticles_;
  
  for(std::vector<math::XYZTLorentzVectorD>::const_iterator it = leptonMomenta_.begin(); it!=leptonMomenta_.end(); it++){
    ptParticles_.push_back(it->pt());
    //std::cout<<it->pt()<<std::endl;
  }
  
  for (unsigned int i = 0; i < ptParticles_.size(); i++) {   
    for (unsigned int j = i+1; j < ptParticles_.size(); j++) {
      if(ptParticles_[j] > ptParticles_[i]) {
	math::XYZTLorentzVectorD auxMomentum = sorted4Momenta_[j];
	sorted4Momenta_[j] = sorted4Momenta_[i];
	sorted4Momenta_[i] = auxMomentum;
      }// if inverted
    }// loop on second
  }//loop on first 
  
  return sorted4Momenta_;
}

// ------------------------------- sort leptons by pt -----------------------------------------
std::vector<const reco::Candidate*>
ZbbAcceptanceCalculator::sortLeptonsByPt(std::vector<const reco::Candidate*> leptons_)
{
  //std::cout<<"Sorting"<<std::endl;
  std::vector<const reco::Candidate*> sortedLeptons_ = leptons_;
  std::vector<double> ptParticles_;

  for(std::vector<const reco::Candidate*>::const_iterator it = leptons_.begin(); it!=leptons_.end(); it++){
    ptParticles_.push_back((*it)->pt());
    //std::cout<<it->pt()<<std::endl;
  }
  
  for (unsigned int i = 0; i < ptParticles_.size(); i++) {   
    for (unsigned int j = i+1; j < ptParticles_.size(); j++) {
      if(ptParticles_[j] > ptParticles_[i]) {
	const reco::Candidate* auxLepton = sortedLeptons_[j];
	sortedLeptons_[j] = sortedLeptons_[i];
	sortedLeptons_[i] = auxLepton;
      }// if inverted
    }// loop on second
  }//loop on first 
  
  return sortedLeptons_;
}

//______________________________________________________________________________
std::vector<reco::CompositeCandidate> 
ZbbAcceptanceCalculator::sortCandidatesByDifference(std::vector<reco::CompositeCandidate> unsortedCands){

  // if(debug_) std::cout<<"ZbbEventContentAnalyzer::sortCandidatesByDifference"<<std::endl;
  std::vector<reco::CompositeCandidate> sortedCands = unsortedCands;
  Double_t mZ=91.1876;
  std::vector<Double_t> diffZmass; 
  unsigned int ZCandSize=unsortedCands.size();
  
  if(unsortedCands.size()>1){
    // if(debug_) std::cout<<"ZbbEventContentAnalyzer::sortCandidatesByDifference: Z candidates= "<<ZCandSize<<std::endl;

    for(reco::CompositeCandidateCollection::const_iterator ZllCandidate=unsortedCands.begin(); ZllCandidate != unsortedCands.end(); ++ZllCandidate){
      diffZmass.push_back(fabs(mZ - ZllCandidate->mass()));
    }
    
    for (unsigned int i = 0; i < ZCandSize; i++) {   
      for (unsigned int j = i+1; j < ZCandSize; j++) {
	if(diffZmass[i] > diffZmass[j]) {
	  reco::CompositeCandidate auxCand = sortedCands[i];
	  sortedCands[i] = sortedCands[j];
	  sortedCands[j] = auxCand;   
	}
      }    
    } 
      
    for(reco::CompositeCandidateCollection::const_iterator ZllCandSorted=sortedCands.begin(); ZllCandSorted != sortedCands.end(); ++ZllCandSorted){
      // if(debug_) std::cout<<"Candidate mass:"<<ZllCandSorted->mass()<<std::endl;
    }
  }
  
  return sortedCands;
}

//______________________________________________________________________________
// powerful function to match all jets
//______________________________________________________________________________
std::vector< std::pair<const pat::Jet*,const reco::GenJet*> > 
ZbbAcceptanceCalculator::matchJetsByDr( std::vector<const pat::Jet*> theRecoJets, 
					std::vector<const reco::GenJet*> theGenJets,
					double const& maxDR,
					bool const& uniqueFirst,
					bool const& uniqueSecond ){
  
  unsigned n = theRecoJets.size();
  unsigned m = theGenJets.size();
  std::vector< std::pair<const pat::Jet*,const reco::GenJet*> > pairVec;
  std::vector< std::pair<const pat::Jet*,const reco::GenJet*> > pairVec_raw;
  pairVec_raw.resize(n*m);
  unsigned vecIndex = 0;
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < m; ++j, ++vecIndex) {
      pairVec_raw[vecIndex] = (std::pair<const pat::Jet*,const reco::GenJet*>(theRecoJets[i],theGenJets[j]));
    }
  }

  //std::cout<<"----------------inside matching method; pairVec.size()="<<pairVec.size()<<std::endl;
  // delete all matches outside the cut
  for(unsigned int k = 0; k <pairVec_raw.size(); k++ ){
    //std::cout<<"pair["<<k<<"] DeltaR ="<<ROOT::Math::VectorUtil::DeltaR( (*pairVec[k].first).momentum(),(*pairVec[k].second).momentum())<<std::endl;
    if(ROOT::Math::VectorUtil::DeltaR( (*pairVec_raw[k].first).momentum(),(*pairVec_raw[k].second).momentum()) < maxDR ) pairVec.push_back(pairVec_raw[k]); //pairVec.erase(pairVec.begin()+k); 
  }
  
  //  std::cout<<"after cleaning"<<std::endl;
  //   for(unsigned int i = 0; i <pairVec.size(); i++ ){
  //     std::cout<<"pair["<<i<<"] DeltaR ="<<ROOT::Math::VectorUtil::DeltaR( (*pairVec[i].first).momentum(),(*pairVec[i].second).momentum())<<std::endl;
  //   }
  //   std::cout<<"----------------outside matching method"<<std::endl;  

  // sort the pairs in ascending dR
  for(unsigned int i = 0; i <pairVec.size(); i++ ){
    for(unsigned int j = i+1; j<pairVec.size(); j++  ){
      if(DRCompare(pairVec[i],pairVec[j])){
	std::pair<const pat::Jet*,const reco::GenJet*> auxPair = pairVec[i];
	pairVec[i] = pairVec[j];
	pairVec[j] = auxPair;
      }
    }
  }
  
  // starts the selection of the case
  if (uniqueFirst && uniqueSecond) {
    std::vector< std::pair<const pat::Jet*,const reco::GenJet*> > matchedPairs;
    std::vector<const pat::Jet*> fVec;
    std::vector<const reco::GenJet*> sVec;
    for(std::vector<std::pair<const pat::Jet*,const reco::GenJet*> >::const_iterator aPair=pairVec.begin(); aPair!=pairVec.end(); aPair++){
      bool inFVec = std::count(fVec.begin(),fVec.end(),(*aPair).first);
      bool inSVec = std::count(sVec.begin(),sVec.end(),(*aPair).second);
      if (!inFVec && !inSVec) {
	matchedPairs.push_back(*aPair);
	fVec.push_back((*aPair).first);
	sVec.push_back((*aPair).second);
      }
    }
    return matchedPairs;
  } else if (uniqueFirst) {
    std::vector< std::pair<const pat::Jet*,const reco::GenJet*> > matchedPairs;
    std::vector<const pat::Jet*> fVec;
    for(std::vector<std::pair<const pat::Jet*,const reco::GenJet*> >::const_iterator aPair=pairVec.begin(); aPair!=pairVec.end(); aPair++){
      bool inFVec = std::count(fVec.begin(),fVec.end(),(*aPair).first);
      if (!inFVec) {
	matchedPairs.push_back(*aPair);
	fVec.push_back((*aPair).first);
      }
    }
    return matchedPairs;
  } else if (uniqueSecond) {
    std::vector< std::pair<const pat::Jet*,const reco::GenJet*> > matchedPairs;
    std::vector<const reco::GenJet*> sVec;
    for(std::vector<std::pair<const pat::Jet*,const reco::GenJet*> >::const_iterator aPair=pairVec.begin(); aPair!=pairVec.end(); aPair++){
      bool inSVec = std::count(sVec.begin(),sVec.end(),(*aPair).second);
      if (!inSVec) {
	matchedPairs.push_back(*aPair);
	sVec.push_back((*aPair).second);
      }
    }
    return matchedPairs;
  }
  return pairVec;
}

//______________________________________________________________________________
bool 
ZbbAcceptanceCalculator::DRCompare(std::pair<const pat::Jet*,const reco::GenJet*> p1, std::pair<const pat::Jet*,const reco::GenJet*> p2) {
  double dR1 = ROOT::Math::VectorUtil::DeltaR((*p1.first).momentum(),(*p1.second).momentum());
  double dR2 = ROOT::Math::VectorUtil::DeltaR((*p2.first).momentum(),(*p2.second).momentum());
  return (dR1 > dR2);
}

//______________________________________________________________________________
// even more powerful function to match all type of objects
//______________________________________________________________________________
template<class T, class U> 
std::vector< std::pair<T,U> > ZbbAcceptanceCalculator::MatchByDR(std::vector<T> const& c1,
								 std::vector<U> const& c2,
								 double const& maxDR,
								 bool const& uniqueFirst,
								 bool const& uniqueSecond) {
  
  unsigned n = c1.size();
  unsigned m = c2.size();
  std::vector< std::pair<T,U> > pairVec;
  std::vector< std::pair<T,U> > pairVec_raw;  
  pairVec_raw.resize(n*m);
  unsigned vecIndex = 0;
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < m; ++j, ++vecIndex) {
      pairVec_raw[vecIndex] = (std::pair<T,U>(c1[i],c2[j]));
    }
  }

  // fills the the vector only if deltaR condition is satisfied
  for(unsigned int k = 0; k <pairVec_raw.size(); k++ ){
    if(ROOT::Math::VectorUtil::DeltaR( (*pairVec_raw[k].first).momentum(),(*pairVec_raw[k].second).momentum()) < maxDR ) pairVec.push_back(pairVec_raw[k]); 
  }

  // sort the pairs in ascending dR
  for(unsigned int i = 0; i <pairVec.size(); i++ ){
    for(unsigned int j = i+1; j<pairVec.size(); j++  ){
      if(TemplDRCompare(pairVec[i],pairVec[j])){
	std::pair<T,U> auxPair = pairVec[i];
	pairVec[i] = pairVec[j];
	pairVec[j] = auxPair;
      }
    }
  }

  if (uniqueFirst && uniqueSecond) {
    std::vector< std::pair<T,U> > uPairVec;
    std::vector<T> fVec;
    std::vector<U> sVec;
    std::pair<T,U> aPair;
    BOOST_FOREACH(aPair, pairVec) {
      bool inFVec = std::count(fVec.begin(),fVec.end(),aPair.first);
      bool inSVec = std::count(sVec.begin(),sVec.end(),aPair.second);
      if (!inFVec && !inSVec) {
	uPairVec.push_back(aPair);
	fVec.push_back(aPair.first);
	sVec.push_back(aPair.second);
      }
    }
    return uPairVec;
  } else if (uniqueFirst) {
    std::vector< std::pair<T,U> > uPairVec;
    std::vector<T> fVec;
    std::pair<T,U> aPair;
    BOOST_FOREACH(aPair, pairVec) {
      bool inFVec = std::count(fVec.begin(),fVec.end(),aPair.first);
      if (!inFVec) {
	uPairVec.push_back(aPair);
	fVec.push_back(aPair.first);
      }
    }
    return uPairVec;
  } else if (uniqueSecond) {
    std::vector< std::pair<T,U> > uPairVec;
    std::vector<U> sVec;
    std::pair<T,U> aPair;
    BOOST_FOREACH(aPair, pairVec) {
      bool inSVec = std::count(sVec.begin(),sVec.end(),aPair.second);
      if (!inSVec) {
	uPairVec.push_back(aPair);
	sVec.push_back(aPair.second);
      }
    }
    return uPairVec;
  }
  return pairVec;
}

//______________________________________________________________________________
template<class T, class U> 
bool ZbbAcceptanceCalculator::TemplDRCompare(std::pair<T,U> p1, std::pair<T,U> p2) {
  double dR1 = ROOT::Math::VectorUtil::DeltaR((*p1.first).momentum(),(*p1.second).momentum());
  double dR2 = ROOT::Math::VectorUtil::DeltaR((*p2.first).momentum(),(*p2.second).momentum());
  return (dR1 > dR2);
}

//______________________________________________________________________________
Bool_t  ZbbAcceptanceCalculator::isTightZCandidate(reco::CompositeCandidate ZCand, const reco::BeamSpot& beamSpot, Bool_t isMuChannel, const AcceptanceCuts& lCuts){

  Bool_t istightZcandidate = false;
  const reco::Candidate* lep0 = ZCand.daughter(0);
  const reco::Candidate* lep1 = ZCand.daughter(1);

  if(isMuChannel){
    const pat::Muon* muon0 = dynamic_cast<const pat::Muon*>(&(*lep0));
    const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(&(*lep1));    

    if( 
       ( muon0->pt()>lCuts.muPtMin_ &&
	 muon0->isGlobalMuon() && muon0->isTrackerMuon() && 
	 muon0->globalTrack()->normalizedChi2() < 10 && 
	 muon0->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 &&                      
	 muon0->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&                         
	 muon0->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&                         
	 fabs(muon0->dB()) < 0.02 &&                                                                     
	 ((muon0->trackIso() + muon0->caloIso()) <0.15*muon0->pt()) &&                                             
	 muon0->numberOfMatches() > 1 && 
	 abs(muon0->eta()) < lCuts.muEtaMax_ ) &&
       ( muon1->pt()>lCuts.muPtMin_ &&
	 muon1->isGlobalMuon() && muon1->isTrackerMuon() && 
	 muon1->globalTrack()->normalizedChi2() < 10 && 
	 muon1->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 &&                      
	 muon1->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&                         
	 muon1->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&                         
	 fabs(muon1->dB()) < 0.02 &&                                                                     
	 ((muon1->trackIso() + muon1->caloIso()) <0.15*muon1->pt()) &&                                             
	 muon1->numberOfMatches() > 1 && 
	 abs(muon1->eta()) < lCuts.muEtaMax_ )
       ){ 
      istightZcandidate=true;
    }
  } else {
    const pat::Electron* ele0 = dynamic_cast<const pat::Electron*>(&(*lep0));
    const pat::Electron* ele1 = dynamic_cast<const pat::Electron*>(&(*lep1));
    if(
       (fabs(ele0->gsfTrack()->dxy(beamSpot))<0.02 && ele0->pt() > lCuts.elePtMin_ && abs(ele0->eta()) < lCuts.eleEtaMax_ &&                         
	(ele0->isEE() || ele0->isEB()) && !ele0->isEBEEGap() && 
	(fabs(ele0->superCluster()->eta())<1.444 || fabs(ele0->superCluster()->eta())>1.566 ) &&
	(ele0->electronID("eidVBTFRel85") == 7) ) &&    
       (fabs(ele1->gsfTrack()->dxy(beamSpot))<0.02 && ele1->pt() > lCuts.elePtMin_ && abs(ele1->eta()) < lCuts.eleEtaMax_ &&                           
	(ele1->isEE() || ele1->isEB()) && !ele1->isEBEEGap() &&    
	(fabs(ele1->superCluster()->eta())<1.444 || fabs(ele1->superCluster()->eta())>1.566 ) &&
	(ele1->electronID("eidVBTFRel85") == 7) )
       ){
      istightZcandidate=true;
    }
  }
  
  return istightZcandidate;
}

//______________________________________________________________________________
Double_t ZbbAcceptanceCalculator::getMuonScaleFactor(reco::CompositeCandidate &theRecoZcand){

  Double_t muonScaleFactor_(1.);

  Double_t ptmu1= (theRecoZcand.daughter(0))->pt();
  Double_t etamu1=(theRecoZcand.daughter(0))->eta();
  Double_t ptmu2= (theRecoZcand.daughter(1))->pt();
  Double_t etamu2=(theRecoZcand.daughter(1))->eta();
  
  // muon offline scale factors
  Double_t wmuOFF1 =ZbbUtils::getMuonOfflineScaleFactor(ptmu1,etamu1,false);
  Double_t wmuOFF2 =ZbbUtils::getMuonOfflineScaleFactor(ptmu2,etamu2,false);
  
  // HLT Double_TMu7 SF (triggered 2011A1 period)	
  Double_t wmuTRGMu7_1 =ZbbUtils::getMuonTrgScaleFactor_H(ptmu1,etamu1,"2011A1",false);
  Double_t wmuTRGMu7_2 =ZbbUtils::getMuonTrgScaleFactor_L(ptmu2,etamu2,"2011A1",false);
  
  // HLT Mu13Mu8 SF   (triggered 2011A2 period)
  Double_t wmu1TRGMu138_high =ZbbUtils::getMuonTrgScaleFactor_H(ptmu1,etamu1,"2011A2",false);
  Double_t wmu1TRGMu138_low  =ZbbUtils::getMuonTrgScaleFactor_L(ptmu1,etamu1,"2011A2",false);	
  Double_t wmu2TRGMu138_high =ZbbUtils::getMuonTrgScaleFactor_H(ptmu2,etamu2,"2011A2",false);
  Double_t wmu2TRGMu138_low  =ZbbUtils::getMuonTrgScaleFactor_L(ptmu2,etamu2,"2011A2",false); 
  
  // Weighted average: valid for ~ 2 fb^-1 statistics
  Double_t triggerW=0.11*(wmuTRGMu7_1*wmuTRGMu7_2)+0.89*(wmu1TRGMu138_high*wmu2TRGMu138_low + wmu1TRGMu138_low*wmu2TRGMu138_high - wmu1TRGMu138_high*wmu2TRGMu138_high); 
  
  muonScaleFactor_*=(wmuOFF1*wmuOFF2*triggerW);  

  return muonScaleFactor_;

}

//______________________________________________________________________________
Double_t ZbbAcceptanceCalculator::getElectronScaleFactor(reco::CompositeCandidate &theRecoZcand){
  
  Double_t electronScaleFactor_(1.);

  Double_t ptele1= (theRecoZcand.daughter(0))->pt();
  Double_t etaele1=(theRecoZcand.daughter(0))->eta();
  Double_t ptele2= (theRecoZcand.daughter(1))->pt();
  Double_t etaele2=(theRecoZcand.daughter(1))->eta();
  
  // electron offline scale factors
  Double_t weleOFF1 =ZbbUtils::getElectronOfflineScaleFactor(ptele1,etaele1,false);
  Double_t weleOFF2 =ZbbUtils::getElectronOfflineScaleFactor(ptele2,etaele2,false);
  
  // HLT HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5  SF (triggered 2011A1 period)
  Double_t wele1TRG178_v5_high =ZbbUtils::getElectronTrgScaleFactor_H(ptele1,etaele1,"2011A1",false);
  Double_t wele1TRG178_v5_low  =ZbbUtils::getElectronTrgScaleFactor_L(ptele1,etaele1,"2011A1",false);	
  Double_t wele2TRG178_v5_high =ZbbUtils::getElectronTrgScaleFactor_H(ptele2,etaele2,"2011A1",false);
  Double_t wele2TRG178_v5_low  =ZbbUtils::getElectronTrgScaleFactor_L(ptele2,etaele2,"2011A1",false);  
  
  // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6  SF (triggered 2011A2 period)
  Double_t wele1TRG178_v6_high =ZbbUtils::getElectronTrgScaleFactor_H(ptele1,etaele1,"2011A2",false);
  Double_t wele1TRG178_v6_low  =ZbbUtils::getElectronTrgScaleFactor_L(ptele1,etaele1,"2011A2",false);	
  Double_t wele2TRG178_v6_high =ZbbUtils::getElectronTrgScaleFactor_H(ptele2,etaele2,"2011A2",false);
  Double_t wele2TRG178_v6_low  =ZbbUtils::getElectronTrgScaleFactor_L(ptele2,etaele2,"2011A2",false);  
  
  Double_t triggerW = 0.50*(wele1TRG178_v5_high*wele2TRG178_v5_low + wele1TRG178_v5_low*wele2TRG178_v5_high - wele1TRG178_v5_high*wele2TRG178_v5_high)+0.50*(wele1TRG178_v6_high*wele2TRG178_v6_low + wele1TRG178_v6_low*wele2TRG178_v6_high - wele1TRG178_v6_high*wele2TRG178_v6_high);
  
  electronScaleFactor_*=(weleOFF1*weleOFF2*triggerW);
  
  return electronScaleFactor_;

}

//define this as a plug-in
DEFINE_FWK_MODULE(ZbbAcceptanceCalculator);
