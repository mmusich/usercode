// -*- C++ -*-
//
// Package:    SimpleCFAnalyzer
// Class:      SimpleCFAnalyzer
// 
/**\class SimpleCFAnalyzer SimpleCFAnalyzer.cc ZbbAnalysis/AnalysisStep/plugins/SimpleCFAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Musich,40 2-A16,+41227671519,
//         Created:  Wed May  9 18:06:07 CEST 2012
// $Id: SimpleCFAnalyzer.cc,v 1.2 2012/05/14 16:04:20 musich Exp $
//
//


// system include files
#include <memory>

// user include files
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
#include "ZbbAnalysis/AnalysisStep/interface/ZbbTypeDefs.h"
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

struct EventProperties {
  
  bool is_right_flavour_;
  bool is_gen_Z_yes_;
  bool is_gen_j_yes_;
  bool is_gen_b_yes_;
  bool is_gen_Zkin_yes_;
  bool is_rec_Zkin_yes_;
  bool is_rec_b_yes_;
  bool is_rec_j_yes_;
  bool is_rec_lep_idiso_yes_;
  bool is_rec_b_HE_yes_;
  bool is_rec_b_HP_yes_;
  bool is_rec_b_CSVM_yes_;
  bool is_rec_b_CSVT_yes_;            
  Float_t evtWeight_,evtWeightSq_;
  Float_t lepIdIsoWeight_,lepIdIsoWeightSq_;
  Float_t bTagWeightSSVHE_,bTagWeightSSVHESq_;
  Float_t bTagWeightSSVHP_,bTagWeightSSVHPSq_; 
  Float_t bTagWeightCSVM_,bTagWeightCSVMSq_;
  Float_t bTagWeightCSVT_,bTagWeightCSVTSq_;
								
  EventProperties() {
    this->reset();
  }

  void reset(){
    
    //initialize bools and counters
    is_right_flavour_=false;
    is_gen_Z_yes_=false;
    is_gen_j_yes_=false;
    is_gen_b_yes_=false;
    is_gen_Zkin_yes_=false;
    is_rec_Zkin_yes_=false;
    is_rec_b_yes_=false;  
    is_rec_j_yes_=false;  
    is_rec_lep_idiso_yes_=false;
    is_rec_b_HE_yes_=false;
    is_rec_b_HP_yes_=false;
    is_rec_b_CSVM_yes_=false;
    is_rec_b_CSVT_yes_=false;

    evtWeightSq_=1.;
    evtWeight_=1.;
    lepIdIsoWeight_=1.;
    lepIdIsoWeightSq_=1.;
    bTagWeightSSVHE_=1.;
    bTagWeightSSVHESq_=1.;
    bTagWeightSSVHP_=1.;
    bTagWeightSSVHPSq_=1.;
    bTagWeightCSVM_=1.;
    bTagWeightCSVMSq_=1.;
    bTagWeightCSVT_=1.;
    bTagWeightCSVTSq_=1.;
  }

  // getter methods
  bool isRightFlavour()    const {return is_right_flavour_;}
  bool isGenZyes()         const {return is_gen_Z_yes_;}
  bool isGenbyes()         const {return is_gen_b_yes_;}
  bool isGenjetyes()       const {return is_gen_j_yes_;}
  bool isGen_Zkin_yes()    const {return is_gen_Zkin_yes_;}
  bool isRec_Zkin_yes()    const {return is_rec_Zkin_yes_;}
  bool isRecjet_yes()      const {return is_rec_j_yes_;}
  bool isRecb_yes()        const {return is_rec_b_yes_;}
  bool isRecLepIdIso_yes() const {return is_rec_lep_idiso_yes_;}
  bool isRecbHE_yes()      const {return is_rec_b_HE_yes_;}
  bool isRecbHP_yes()      const {return is_rec_b_HP_yes_;}
  bool isRecbCSVM_yes()    const {return is_rec_b_CSVM_yes_;}
  bool isRecbCSVT_yes()    const {return is_rec_b_CSVT_yes_;}

};

//
// class declaration
//

class SimpleCFAnalyzer : public edm::EDAnalyzer {
public:
  explicit SimpleCFAnalyzer(const edm::ParameterSet&);
  ~SimpleCFAnalyzer();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  // auxiliary methods to reconstruct best Z
  reco::CompositeCandidate theBestZCand(std::vector<const reco::Candidate*> theLeptons);
  std::vector<const reco::Candidate*>   sortLeptonsByPt(std::vector<const reco::Candidate*> leptons_);
  std::vector<reco::CompositeCandidate> sortCandidatesByDifference(std::vector<reco::CompositeCandidate> unsortedCands);

  // method to tag b in the particle
  bool hasBottom(const reco::Candidate &c);
 
  // auxiliary methods to match any two objects
  template<class T, class U> bool dRCompare(std::pair<T,U> p1, std::pair<T,U> p2); 
  template<class T, class U> std::vector< std::pair<T,U> > matchByDR(std::vector<T> const& c1, std::vector<U> const& c2, double const& maxDR, bool const& uFirst, bool const& uSecond);

  // function for espilon_l
  bool   isTightZCandidate(reco::CompositeCandidate ZCand, const reco::BeamSpot& beamSpot, bool isMuChannel);
  Double_t getMuonScaleFactor(reco::CompositeCandidate &theRecoZcand);
  Double_t getElectronScaleFactor(reco::CompositeCandidate &theRecoZcand);

  virtual void endJob() ;
  
  // ----------member data ---------------------------

  // inputtags
  edm::InputTag genParticlesSrc_;
  edm::InputTag genJetSrc_;
  edm::InputTag recoJetSrc_;
  edm::InputTag recoMuonsSrc_;
  edm::InputTag recoElectronsSrc_;
  
  // applied cuts
  Double_t minMassCut_;     
  Double_t maxMassCut_;    
  
  Double_t jetEtaCut_;      
  Double_t jetPtCut_;  

  Double_t muonEtaCut_;     
  Double_t muonPtCut_;      

  Double_t eleEtaCut_;      
  Double_t elePtCut_;       
  
  Double_t dRLeptons_;    
  Double_t dRJets_; 
  
  Double_t theBparticlePtCut_;

  // some bools
  bool verbose_;
  bool saveNTuple_;
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
  edm::VParameterSet  effbcMC_CSVM;
  edm::VParameterSet  effbcMC_CSVT;

  // associative map for b-tag efficiency
  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffbMCSSVHE_;
  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffbMCSSVHP_;

  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffcMCSSVHE_;
  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffcMCSSVHP_;

  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffbMCCSVM_;
  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffbMCCSVT_;

  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffcMCCSVM_;
  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffcMCCSVT_;

  // event properties stored here
  EventProperties myBools_;

  // flavour selection
  bool isZEE, isZMM;

  // histograms
  TH1F *_h_gen_Z_mass;
  TH1F *_h_gen_Zdau1_eta;
  TH1F *_h_gen_Zdau2_eta;
  TH1F *_h_gen_Zdau1_pt;
  TH1F *_h_gen_Zdau2_pt;

  TH1F *_h_gen_bjet_eta;
  TH1F *_h_gen_bjet_pt;
  TH1F *_h_gen_bhad_eta;
  TH1F *_h_gen_bhad_pt;
  
  TH1F *_h_reco_Z_mass;
  TH1F *_h_reco_Zdau1_eta;
  TH1F *_h_reco_Zdau2_eta;
  TH1F *_h_reco_Zdau1_pt;
  TH1F *_h_reco_Zdau2_pt;

  TH1F *_h_reco_bjet_eta;
  TH1F *_h_reco_bjet_pt;
  TH1F *_h_bjet_energyresponse;

  // output ntuple
  TTree *myNTuple_;

};

SimpleCFAnalyzer::SimpleCFAnalyzer(const edm::ParameterSet& iConfig)

{

  //get the collections
  genParticlesSrc_  =iConfig.getParameter<edm::InputTag>("GenparticleCollection");
  recoJetSrc_       =iConfig.getParameter<edm::InputTag>("JetCollection");
  recoMuonsSrc_     =iConfig.getParameter<edm::InputTag>("MuonCollection");
  recoElectronsSrc_ =iConfig.getParameter<edm::InputTag>("ElectronCollection");
  genJetSrc_        =iConfig.getParameter<edm::InputTag>("GenJetCollection");

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
  
  decay_chain_string_      = iConfig.getUntrackedParameter<std::string>("DecayChainSelection");
  applypucorrection_       = iConfig.getUntrackedParameter<bool>("applyPUCorr", true); 
  verbose_                 = iConfig.getUntrackedParameter<bool>("verbose", false);
  saveNTuple_              = iConfig.getUntrackedParameter<bool>("saveNTuple", false);
  useClopperPearsonErrors_ = iConfig.getUntrackedParameter<bool>("useClopperPearsonErrors", false);
  partonLevel_             = iConfig.getUntrackedParameter<bool>("PartonLevel", false);
  
  effbcMC_SSVHE            = iConfig.getParameter<edm::VParameterSet>("EffbcMCSSVHE");
  effbcMC_SSVHP            = iConfig.getParameter<edm::VParameterSet>("EffbcMCSSVHP");
  effbcMC_CSVM             = iConfig.getParameter<edm::VParameterSet>("EffbcMCCSVM");
  effbcMC_CSVT             = iConfig.getParameter<edm::VParameterSet>("EffbcMCCSVT");

  //initialize the lumiweight
  theLumiWeights_=0;

  //select the flavour of leptons from Z
  isZEE=false; isZMM=false;

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
  
  // b-tag MC efficiency CSVM tagger
  EffMCpair theEffbCSVM_;
  EffMCpair theEffcCSVM_;
  for(edm::VParameterSet::const_iterator pset =  effbcMC_CSVM.begin(); pset != effbcMC_CSVM.end(); pset++){
    std::vector<double> thePtRange_ = pset->getUntrackedParameter<std::vector<double> >("ptrange");
    double theEffb_barrel_  = pset->getUntrackedParameter<double>("effb_barrel");
    double theEffb_forward_ = pset->getUntrackedParameter<double>("effb_forward");
    double theEffc_barrel_  = pset->getUntrackedParameter<double>("effc_barrel");
    double theEffc_forward_ = pset->getUntrackedParameter<double>("effc_forward");
    theEffbCSVM_ = std::make_pair(theEffb_barrel_,theEffb_forward_);
    theEffcCSVM_ = std::make_pair(theEffc_barrel_,theEffc_forward_);
    theAssociativeMapEffbMCCSVM_[theEffbCSVM_]=thePtRange_; 
    theAssociativeMapEffcMCCSVM_[theEffcCSVM_]=thePtRange_; 
  }

  // b-tag MC efficiency CSVT tagger
  EffMCpair theEffbCSVT_;
  EffMCpair theEffcCSVT_;
  for(edm::VParameterSet::const_iterator pset =  effbcMC_CSVT.begin(); pset != effbcMC_CSVT.end(); pset++){
    std::vector<double> thePtRange_ = pset->getUntrackedParameter<std::vector<double> >("ptrange");
    double theEffb_barrel_  = pset->getUntrackedParameter<double>("effb_barrel");
    double theEffb_forward_ = pset->getUntrackedParameter<double>("effb_forward");
    double theEffc_barrel_  = pset->getUntrackedParameter<double>("effc_barrel");
    double theEffc_forward_ = pset->getUntrackedParameter<double>("effc_forward");
    theEffbCSVT_ = std::make_pair(theEffb_barrel_,theEffb_forward_);
    theEffcCSVT_ = std::make_pair(theEffc_barrel_,theEffc_forward_);
    theAssociativeMapEffbMCCSVT_[theEffbCSVT_]=thePtRange_; 
    theAssociativeMapEffcMCCSVT_[theEffcCSVT_]=thePtRange_; 
 }
}


SimpleCFAnalyzer::~SimpleCFAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SimpleCFAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   // define the 4vector of the Z and of the gen leptons
   math::XYZTLorentzVectorD pLpos, pLneg, pZLL;
  
   // initialize the gen 4vectors
   pLpos.SetXYZT(0,0,0,0);
   pLneg.SetXYZT(0,0,0,0);
   pZLL.SetXYZT(0,0,0,0);

   //------------------------ Handles -------------------------------
  
   // get gen particle candidates
   edm::Handle<reco::GenParticleCollection> genParticlesCollection;
   iEvent.getByLabel(genParticlesSrc_, genParticlesCollection);

   // get patMuons
   edm::Handle<edm::View<pat::Muon> > recoMuonsHandle;
   iEvent.getByLabel(recoMuonsSrc_,recoMuonsHandle);
   const edm::View<pat::Muon> & recoMuons = *(recoMuonsHandle.product());

   // get patElectrons
   edm::Handle<edm::View<pat::Electron> > recoElectronsHandle;
   iEvent.getByLabel(recoElectronsSrc_,recoElectronsHandle);
   const edm::View<pat::Electron> & recoElectrons = *(recoElectronsHandle.product());
 
   // get reco patJets
   edm::Handle<edm::View<pat::Jet> > recoJetsHandle;
   iEvent.getByLabel(recoJetSrc_,recoJetsHandle);
   const edm::View<pat::Jet> & recoJets = *(recoJetsHandle.product());
   
   // get genJets
   edm::Handle<reco::GenJetCollection> genJetsHandle;
   iEvent.getByLabel(genJetSrc_,genJetsHandle);
   const reco::GenJetCollection & genJets = *(genJetsHandle.product());
   
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
   } else {
     if( verbose_ )  std::cout << "SimpleCFAnalyzer::analyze() No BeamSpot found!" << std::endl;
   }
   
   // applying PU corrections
   if(applypucorrection_){
     double thePUWeight_(1.);
     const edm::EventBase* iEventB = dynamic_cast<const edm::EventBase*>(&iEvent);
     if (theLumiWeights_) thePUWeight_ = theLumiWeights_->weight(*iEventB);
     wMC*=thePUWeight_;
   }

   myBools_.reset();

   myBools_.evtWeight_=wMC;
   myBools_.evtWeightSq_=(wMC*wMC);
   
   double wLepIdIso(1.);
   double wBtagHE(1.);
   double wBtagHP(1.);
   double wBtagCSVM(1.);
   double wBtagCSVT(1.);
   
   //----------------------------------------------------------------------------------------------------
   //2) Discard event if the gen. Z does not decay to the desired flavour of lepton.

   std::vector<const reco::Candidate*> theGenMuons;
   std::vector<math::XYZTLorentzVectorD> theGenMuonsMomenta;

   std::vector<const reco::Candidate*> theGenElectrons;
   std::vector<math::XYZTLorentzVectorD> theGenElectronsMomenta;

   for( reco::GenParticleCollection::const_iterator genp = genParticlesCollection->begin();genp != genParticlesCollection->end(); ++ genp ) { 
     
     // make a copy of genp for better handling
     const reco::Candidate* p = &(*genp);
     
     if(isZMM){
       if( abs( p->pdgId() )==13 && ( p->status())==3 ){
	 theGenMuons.push_back(p);
	 math::XYZTLorentzVectorD p4Lepttemp(p->px(),p->py(),p->pz(),p->energy());  
	 theGenMuonsMomenta.push_back(p4Lepttemp);
       }
     } else if(isZEE) {
       if( abs( p->pdgId() )==11 && ( p->status())==3 ){
	 theGenElectrons.push_back(p);
	 math::XYZTLorentzVectorD p4Lepttemp(p->px(),p->py(),p->pz(),p->energy());  
	 theGenElectronsMomenta.push_back(p4Lepttemp);
       }
     } else {  
       std::cout<<"SimpleCFAnalyzer::analyze() this case is not contemplated at the moment!"<<std::endl;
     } // ends swithces on flavour
   } // ends lepton master loop on genParticles
   
   
   if (isZMM && theGenMuons.size()==2){
     myBools_.is_right_flavour_=true;
   } else if (isZEE && theGenElectrons.size()==2){
     myBools_.is_right_flavour_=true;
   }

   if(!myBools_.isRightFlavour()) return;  
   
   //----------------------------------------------------------------------------------------------------
   // 3) For Z->ee,  take the two status 3 daughter leptons as the Z candidate pair. For Z->mumu, take the status 1 pair most compatible with the Z mass.
   // 4) Label the event as gen_Z_yes if the gen. lepton pair has mass 60 < m_ll < 120 and the leptons have opposite charge.

   reco::CompositeCandidate myBestZCand;

   if(isZMM){
     myBestZCand = theBestZCand(theGenMuons);
   } else if(isZEE){
     myBestZCand = theBestZCand(theGenElectrons); 
   }

   if (myBestZCand.daughter(0)->pdgId() < 0 ) {
     pLneg=myBestZCand.daughter(0)->p4();
     pLpos=myBestZCand.daughter(1)->p4();
   } else if ( myBestZCand.daughter(0)->pdgId() > 0  ) {
     pLneg=myBestZCand.daughter(1)->p4();
     pLpos=myBestZCand.daughter(0)->p4();
   }
   
   pZLL=pLpos+pLneg;
   Double_t mLL=pZLL.mass();
   
   if(mLL>minMassCut_ && mLL<maxMassCut_){
     myBools_.is_gen_Z_yes_=true;
     _h_gen_Z_mass->Fill(mLL);
   } 
   
   //----------------------------------------------------------------------------------------------------
   //5) Out of the lepton pair, count how many would pass the cuts pT > 25 (20), |Eta| < 2.5 (2.1) for the electrons (muons).  If both leptons pass, label the event gen_Zkin_yes.

   Int_t nOkLeptons(0);
   if(isZMM){
     if(pLpos.pt() > muonPtCut_ && TMath::Abs(pLpos.eta()) < muonEtaCut_) nOkLeptons++;
     if(pLneg.pt() > muonPtCut_ && TMath::Abs(pLneg.eta()) < muonEtaCut_) nOkLeptons++;
   } else if (isZEE){
     if(pLpos.pt() > elePtCut_  && TMath::Abs(pLpos.eta()) < eleEtaCut_ ) nOkLeptons++;
     if(pLneg.pt() > elePtCut_  && TMath::Abs(pLneg.eta()) < eleEtaCut_ ) nOkLeptons++;
   }

   if(nOkLeptons ==2) {
     myBools_.is_gen_Zkin_yes_=true;
     _h_gen_Zdau1_eta->Fill(pLpos.eta()); 
     _h_gen_Zdau2_eta->Fill(pLpos.eta()); 
     _h_gen_Zdau1_pt->Fill(pLneg.pt());  
     _h_gen_Zdau2_pt->Fill(pLpos.pt());  
   } 

   //----------------------------------------------------------------------------------------------------
   //5) Match gen-jets to the last B hadrons before decay. To identify these, we take any B hadron that does not have a daughter that is also
   //a B hadron. In the matching algorithm, the maximum Delta_R is 0.5 and each B hadron can be matched to at most one gen-jet (the best match is taken).  Call the matched gen-jets

   std::vector<const reco::Candidate*> the_last_bhadrons;
   std::vector<const reco::GenJet*>    the_b_flavour_genJets;
   std::vector<const reco::GenJet*>    the_any_flavour_genJets;
   std::vector<const reco::GenJet*>    the_filtered_any_flavour_genJets;
   std::vector<const reco::GenJet*>    the_filtered_b_flavour_genJets;
   
   for( reco::GenParticleCollection::const_iterator genp = genParticlesCollection->begin();genp != genParticlesCollection->end(); ++ genp ) { 
      
     const reco::Candidate* p = &(*genp);
    
     if(hasBottom(*p)){ 
       
       bool hasBottomedDaughter = false;
       for(UInt_t i=0; i<p->numberOfDaughters(); i++){
	 if(hasBottom(*(p->daughter(i)))) hasBottomedDaughter=true;
       }
       
       if(!hasBottomedDaughter && p->pt() > theBparticlePtCut_){  
	 the_last_bhadrons.push_back(p);
       }
     }
   }
      
   for(std::vector<reco::GenJet>::const_iterator genjet_it=genJets.begin(); genjet_it!=genJets.end(); ++genjet_it){ 
     the_any_flavour_genJets.push_back(&(*genjet_it));
   }   
   
   // do the actual matching
   std::vector<std::pair<const reco::GenJet*,const reco::Candidate*> > bHad_genJetMatch = matchByDR(the_any_flavour_genJets,the_last_bhadrons,dRJets_,true,true);

   if(verbose_) std::cout<<"SimpleCFAnalyzer::analyze() size of bHad_genJetMatch= "<<bHad_genJetMatch.size()<<std::endl;

   for (unsigned i = 0; i < bHad_genJetMatch.size(); ++i){ 
     the_b_flavour_genJets.push_back(bHad_genJetMatch[i].first);
     _h_gen_bjet_eta->Fill(bHad_genJetMatch[i].first->eta()); 
     _h_gen_bjet_pt->Fill(bHad_genJetMatch[i].first->pt());  
     _h_gen_bhad_eta->Fill(bHad_genJetMatch[i].second->eta()); 
     _h_gen_bhad_pt->Fill(bHad_genJetMatch[i].second->pt()); 
   }

   //----------------------------------------------------------------------------------------------------
   //6) Count the number of these b-flavour gen-jets with pT > 25, |Eta| < 2.1 and Delta_R > 0.5 to each of the two leptons. If >= 1, label events as gen_b_yes.
   Int_t nOkBFlavouredJets(0);
   
   for(std::vector<const reco::GenJet*>::const_iterator bgenjet_it=the_b_flavour_genJets.begin(); bgenjet_it!=the_b_flavour_genJets.end(); ++bgenjet_it){
     if((*bgenjet_it)->pt()>jetPtCut_ && TMath::Abs((*bgenjet_it)->eta())<jetEtaCut_  &&
	ROOT::Math::VectorUtil::DeltaR((*bgenjet_it)->momentum(),pLpos.Vect())>dRJets_ &&
	ROOT::Math::VectorUtil::DeltaR((*bgenjet_it)->momentum(),pLneg.Vect())>dRJets_
	){
       nOkBFlavouredJets++;
       the_filtered_b_flavour_genJets.push_back(*bgenjet_it);
     }
   }

   if(nOkBFlavouredJets>=1) {
     myBools_.is_gen_b_yes_=true;
   }
   
   // ========================= study any flavour jets ===========================
   Int_t nOkJets(0);
   
   for(std::vector<reco::GenJet>::const_iterator genjet_it=genJets.begin(); genjet_it!=genJets.end(); ++genjet_it){
     if(genjet_it->pt()>jetPtCut_ && TMath::Abs(genjet_it->eta())<jetEtaCut_ &&
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

   //----------------------------------------------------------------------------------------------------
   //9) Match the reco leptons to the two gen leptons (regardless of whether the gen leptons passed/failed any of the cuts above). We used a Delta_R limit of 0.3 for this matching, but this may be excessive. If two reco matches are found, apply the pT and |Eta| cuts, and if both leptons survive, determine if the reco pair has 60 < m_ll < 120 and opposite charge.  If it does, label the event rec_Zkin_yes, and discard all other reco leptons
     
   std::vector<const pat::Muon*>      the_recoMuons; 
   std::vector<const pat::Electron*>  the_recoElectrons; 

   reco::CompositeCandidate myRecoZCand;
   std::vector<math::XYZTLorentzVectorD> recoLeptonMomenta; 

   for(edm::View<pat::Muon>::const_iterator muon=recoMuons.begin(); muon!=recoMuons.end(); ++muon){
     the_recoMuons.push_back(&(*muon));
   }
   
   for(edm::View<pat::Electron>::const_iterator electron= recoElectrons.begin(); electron!= recoElectrons.end(); ++electron){
     the_recoElectrons.push_back(&(*electron));
   }
   
   if(isZMM){

     std::vector<const reco::Candidate*> the_matchedMuons; 

     // do the actual matching
     std::vector<std::pair<const pat::Muon*,const reco::Candidate*> > recoMu_genMuMatch = matchByDR(the_recoMuons,theGenMuons,dRLeptons_,true,true);
     for (unsigned i = 0; i < recoMu_genMuMatch.size(); ++i){ 
       the_matchedMuons.push_back(recoMu_genMuMatch[i].first);
     }
     
     if(verbose_) std::cout<<"SimpleCFAnalyzer::analyze() size of recoMu_genMuMatch= "<<recoMu_genMuMatch.size()<<std::endl;

     // check 2 matches
     if(the_matchedMuons.size()==2){
       Int_t nOkRecoMuons(0.);
       if(the_matchedMuons[0]->pt() > muonPtCut_ && TMath::Abs(the_matchedMuons[0]->eta()) < muonEtaCut_) nOkRecoMuons++;
       if(the_matchedMuons[1]->pt() > muonPtCut_ && TMath::Abs(the_matchedMuons[1]->eta()) < muonEtaCut_) nOkRecoMuons++;

       if(verbose_) std::cout<<"SimpleCFAnalyzer::analyze() nOkLept= "<<nOkRecoMuons<<std::endl;
       if(verbose_) std::cout<<"SimpleCFAnalyzer::analyze() q(LL)=   "<<the_matchedMuons[0]->charge()*the_matchedMuons[1]->charge()<<std::endl;

       // check 2 passes + opposite sign
       if( nOkRecoMuons==2 && (the_matchedMuons[0]->charge()*the_matchedMuons[1]->charge())<0.){
	 
	 math::XYZTLorentzVectorD p1reco(the_matchedMuons[0]->px(),the_matchedMuons[0]->py(),the_matchedMuons[0]->pz(),the_matchedMuons[0]->energy()); 
	 math::XYZTLorentzVectorD p2reco(the_matchedMuons[1]->px(),the_matchedMuons[1]->py(),the_matchedMuons[1]->pz(),the_matchedMuons[1]->energy()); 
	 
	 recoLeptonMomenta.push_back(p1reco);
	 recoLeptonMomenta.push_back(p2reco);

	 myRecoZCand = theBestZCand(the_matchedMuons);
 
	 math::XYZTLorentzVectorD pZreco = p1reco+p2reco;
	 Double_t mLLreco=pZreco.mass();

	 if(verbose_) std::cout<<"SimpleCFAnalyzer::analyze() mLL= "<<mLLreco<<std::endl;
	 
	 if(mLLreco>minMassCut_ && mLLreco<maxMassCut_){
	   myBools_.is_rec_Zkin_yes_=true;
	   _h_reco_Z_mass->Fill(mLLreco);    
	   _h_reco_Zdau1_eta->Fill(p1reco.eta()); 
	   _h_reco_Zdau2_eta->Fill(p2reco.eta()); 
	   _h_reco_Zdau1_pt->Fill(p1reco.pt());  
	   _h_reco_Zdau2_pt->Fill(p2reco.pt());  
	 } 
       }
     } 
   } else if (isZEE){

     std::vector<const reco::Candidate*>  the_matchedElectrons;

     // do the actual matching
     std::vector<std::pair<const pat::Electron*,const reco::Candidate*> > recoEle_genEleMatch = matchByDR(the_recoElectrons,theGenElectrons,dRLeptons_,true,true);
     for (unsigned i = 0; i < recoEle_genEleMatch.size(); ++i){ 
       the_matchedElectrons.push_back(recoEle_genEleMatch[i].first);
     }

     if(verbose_) std::cout<<"SimpleCFAnalyzer::analyze() size of recoEle_genEleMatch= "<<recoEle_genEleMatch.size()<<std::endl;

      // check 2 matches
     if(the_matchedElectrons.size()==2){
       Int_t nOkRecoElectrons(0.);

       const pat::Electron* ele0 = dynamic_cast<const pat::Electron*>(&(*the_matchedElectrons[0]));
       const pat::Electron* ele1 = dynamic_cast<const pat::Electron*>(&(*the_matchedElectrons[0]));

       if(the_matchedElectrons[0]->pt() > elePtCut_ && TMath::Abs(the_matchedElectrons[0]->eta()) < eleEtaCut_ && 
	  (fabs(ele0->superCluster()->eta())<1.444 || fabs(ele0->superCluster()->eta())>1.566 )) nOkRecoElectrons++;
       if(the_matchedElectrons[1]->pt() > elePtCut_ && TMath::Abs(the_matchedElectrons[1]->eta()) < eleEtaCut_ &&
	  (fabs(ele1->superCluster()->eta())<1.444 || fabs(ele1->superCluster()->eta())>1.566 )) nOkRecoElectrons++;
       
       // check 2 passes + opposite sign
       if( nOkRecoElectrons==2 && (the_matchedElectrons[0]->charge()*the_matchedElectrons[1]->charge())<0.){
	 
	 math::XYZTLorentzVectorD p1reco(the_matchedElectrons[0]->px(),the_matchedElectrons[0]->py(),the_matchedElectrons[0]->pz(),the_matchedElectrons[0]->energy()); 
	 math::XYZTLorentzVectorD p2reco(the_matchedElectrons[1]->px(),the_matchedElectrons[1]->py(),the_matchedElectrons[1]->pz(),the_matchedElectrons[1]->energy()); 

	 recoLeptonMomenta.push_back(p1reco);
	 recoLeptonMomenta.push_back(p2reco);

	 myRecoZCand = theBestZCand(the_matchedElectrons);
	 
	 math::XYZTLorentzVectorD pZreco = p1reco+p2reco;
	 Double_t mLLreco=pZreco.mass();
   
	 if(mLLreco>minMassCut_ && mLLreco<maxMassCut_){
	   myBools_.is_rec_Zkin_yes_=true;
	   _h_reco_Z_mass->Fill(mLLreco);    
	   _h_reco_Zdau1_eta->Fill(p1reco.eta()); 
	   _h_reco_Zdau2_eta->Fill(p2reco.eta()); 
	   _h_reco_Zdau1_pt->Fill(p1reco.pt());  
	   _h_reco_Zdau2_pt->Fill(p2reco.pt()); 
	 } 
       }
     } 
   }

   //----------------------------------------------------------------------------------------------------
   // 10)  Filter the reco jets on pT (25.0), |Eta| (2.1) and Delta_R to the two matched reco leptons (0.5) (if they survived the cuts in step 9).
   std::vector<const pat::Jet*> the_filteredRecoJets; 

   for(edm::View<pat::Jet>::const_iterator recojet_it=recoJets.begin(); recojet_it!=recoJets.end(); ++recojet_it){ 
     if(recojet_it->pt()>jetPtCut_ && TMath::Abs(recojet_it->eta())<jetEtaCut_) {
       if(myBools_.is_rec_Zkin_yes_){
	 if(ROOT::Math::VectorUtil::DeltaR(recojet_it->momentum(), recoLeptonMomenta[0].Vect() )>dRJets_ && 
	    ROOT::Math::VectorUtil::DeltaR(recojet_it->momentum(), recoLeptonMomenta[1].Vect() )>dRJets_ ){
	   
	   the_filteredRecoJets.push_back(&(*recojet_it));
	 }  
       } else { 
	 the_filteredRecoJets.push_back(&(*recojet_it));
       }
     } // ends if reco jet is ok
   } // for loop on reco jets 
   
   //----------------------------------------------------------------------------------------------------
   //11) Match the remaining reco jets to the collection of b-flavour gen-jets with Delta_R (0.5).  If at least one reco jet is matched, label the event rec_b_yes.
   std::vector<std::pair<const pat::Jet*,const reco::GenJet*> > filtRecoJet_bGenJetMatch = matchByDR(the_filteredRecoJets,the_b_flavour_genJets,dRJets_,true,true);
   if (filtRecoJet_bGenJetMatch.size() > 0 ) {
     myBools_.is_rec_b_yes_=true; 
   }
   
   std::vector<const pat::Jet*> the_filtRecoMatchedJets; 
   for (unsigned i = 0; i < filtRecoJet_bGenJetMatch.size(); ++i){ 
     the_filtRecoMatchedJets.push_back(filtRecoJet_bGenJetMatch[i].first);
     _h_reco_bjet_eta->Fill(filtRecoJet_bGenJetMatch[i].first->eta());     
     _h_reco_bjet_pt->Fill(filtRecoJet_bGenJetMatch[i].first->pt());      
     _h_bjet_energyresponse->Fill((filtRecoJet_bGenJetMatch[i].first)->pt()/(filtRecoJet_bGenJetMatch[i].second)->pt());
   }

   // matching of any-flavour genJets
   std::vector<std::pair<const pat::Jet*,const reco::GenJet*> > filtRecoJet_GenJetMatch  = matchByDR(the_filteredRecoJets,the_any_flavour_genJets,dRJets_,true,true);
   if (filtRecoJet_GenJetMatch.size() > 0) {
     myBools_.is_rec_j_yes_=true;
   }

   if(verbose_) std::cout<<"SimpleCFAnalyzer::analyze() size of filtRecoJet_GenJetMatch= "<<filtRecoJet_GenJetMatch.size()<<std::endl;


   //----------------------------------------------------------------------------------------------------
   //12) Calculate lepton efficiencies

   bool rec_yes_alljets = (myBools_.isRecjet_yes() && myBools_.isRec_Zkin_yes());
   bool rec_yes         = (myBools_.isRecb_yes() && myBools_.isRec_Zkin_yes());
   
   if(rec_yes || rec_yes_alljets){
     
     if (isZMM){
       if(isTightZCandidate(myRecoZCand,beamSpot,true)){
	 myBools_.is_rec_lep_idiso_yes_=true;	
	 wLepIdIso*=getMuonScaleFactor(myRecoZCand);
	 myBools_.lepIdIsoWeight_=wLepIdIso;
	 myBools_.lepIdIsoWeightSq_=wLepIdIso*wLepIdIso;
       }
     } else if (isZEE){
       if(isTightZCandidate(myRecoZCand,beamSpot,false)){
	 myBools_.is_rec_lep_idiso_yes_=true;
	 wLepIdIso*=getElectronScaleFactor(myRecoZCand);
	 myBools_.lepIdIsoWeight_=wLepIdIso;
	 myBools_.lepIdIsoWeightSq_=wLepIdIso*wLepIdIso;
       }
     }
   }
    
   //----------------------------------------------------------------------------------------------------
   //12) Calculate b-tagging efficiencies
  
  std::vector<pat::Jet> theBtaggedJetsHE;
  std::vector<pat::Jet> theBtaggedJetsHP;  
  std::vector<pat::Jet> theBtaggedJetsCSVM;
  std::vector<pat::Jet> theBtaggedJetsCSVT;

  if(myBools_.isRecLepIdIso_yes()){
    
    for (unsigned i = 0; i < the_filtRecoMatchedJets.size(); ++i){ 
      if(ZbbUtils::isBJet(*the_filtRecoMatchedJets[i],"SSVHEM")){
	theBtaggedJetsHE.push_back(*the_filtRecoMatchedJets[i]);
	myBools_.is_rec_b_HE_yes_=true;
      }
      if(ZbbUtils::isBJet(*the_filtRecoMatchedJets[i],"SSVHPT")){
	theBtaggedJetsHP.push_back(*the_filtRecoMatchedJets[i]);
	myBools_.is_rec_b_HP_yes_=true;
      }
      if(ZbbUtils::isBJet(*the_filtRecoMatchedJets[i],"CSVM")){
	theBtaggedJetsCSVM.push_back(*the_filtRecoMatchedJets[i]);
	myBools_.is_rec_b_CSVM_yes_=true;
      }
      if(ZbbUtils::isBJet(*the_filtRecoMatchedJets[i],"CSVT")){
	theBtaggedJetsCSVT.push_back(*the_filtRecoMatchedJets[i]);
	myBools_.is_rec_b_CSVT_yes_=true;
      }
    }

    // if SSVHE tag in the event
    if(myBools_.isRecbHE_yes()){
      wBtagHE=ZbbUtils::getbEffScaleFactorAR(theAssociativeMapEffbMCSSVHE_,theAssociativeMapEffcMCSSVHE_,theBtaggedJetsHE,iSetup,"BTAGSSVHEM","MISTAGSSVHEM",false,1,0,0); 
      myBools_.bTagWeightSSVHE_=wBtagHE;
      myBools_.bTagWeightSSVHESq_=wBtagHE*wBtagHE;
    } 

    // if SSVHP tag in the event
    if(myBools_.isRecbHP_yes()){
      wBtagHP=ZbbUtils::getbEffScaleFactorAR(theAssociativeMapEffbMCSSVHP_,theAssociativeMapEffcMCSSVHP_,theBtaggedJetsHP,iSetup,"BTAGSSVHPT","MISTAGSSVHPT",false,1,0,0); 
      myBools_.bTagWeightSSVHP_=wBtagHP;
      myBools_.bTagWeightSSVHPSq_=wBtagHP*wBtagHP;
    }

    // if CSVM tag in the event
    if(myBools_.isRecbCSVM_yes()){
      wBtagCSVM=ZbbUtils::getbEffScaleFactorAR(theAssociativeMapEffbMCCSVM_,theAssociativeMapEffcMCCSVM_,theBtaggedJetsCSVM,iSetup,"BTAGCSVM","MISTAGCSVM",false,1,0,0); 
      myBools_.bTagWeightCSVM_=wBtagCSVM;
      myBools_.bTagWeightCSVMSq_=wBtagCSVM*wBtagCSVM;
    }
    
    // if CSVT tag in the event
    if(myBools_.isRecbCSVT_yes()){  
      wBtagCSVT=ZbbUtils::getbEffScaleFactorAR(theAssociativeMapEffbMCCSVT_,theAssociativeMapEffcMCCSVT_,theBtaggedJetsCSVT,iSetup,"BTAGCSVT","MISTAGCSVT",false,1,0,0); 
      myBools_.bTagWeightCSVT_=wBtagCSVT;
      myBools_.bTagWeightCSVTSq_=wBtagCSVT*wBtagCSVT;
    }
  }

  if(saveNTuple_){   
    myNTuple_->Fill();
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
SimpleCFAnalyzer::beginJob()
{
  
  using namespace std;
 
  std::vector<float> dataPU = ZbbUtils::getPU(true);
  std::vector<float> MCPU   = ZbbUtils::getPU(false);
  
  theLumiWeights_ = new edm::LumiReWeighting(MCPU,dataPU);

  edm::Service<TFileService> fs;
  TH1F::SetDefaultSumw2(kTRUE);

  _h_gen_Z_mass    = fs->make<TH1F>("gen_Z_mass","gen Z mass; gen m_{ll} (GeV); events",100,50.,150.);	       
  _h_gen_Zdau1_eta = fs->make<TH1F>("gen_Z_dau1_eta","gen Z daughter #eta; gen #eta_{l}; events",80,-2.8,2.8);      
  _h_gen_Zdau2_eta = fs->make<TH1F>("gen_Z_dau2_eta","gen Z daughter #eta; gen #eta_{l}; events",80,-2.8,2.8);       
  _h_gen_Zdau1_pt  = fs->make<TH1F>("gen_Z_dau1_pt","gen Z daughter p_{T}; gen p_{T}^{l} (GeV); events",150.,0.,150.);        
  _h_gen_Zdau2_pt  = fs->make<TH1F>("gen_Z_dau2_pt","gen Z daughter p_{T}; gen p_{T}^{l} (GeV); events",150.,0.,150.);    
  			       
  _h_gen_bjet_eta  = fs->make<TH1F>("gen_bjet_eta","gen b-jet #eta; gen #eta_{b jet}; events",80,-2.8,2.8);        
  _h_gen_bjet_pt   = fs->make<TH1F>("gen_bjet_pt","gen b-jet p_{T}; gen p_{T}^{b jet} (GeV); events",150.,0.,150); 	       
  _h_gen_bhad_eta  = fs->make<TH1F>("gen_bhad_eta","gen b-hadron #eta; gen #eta_{b had}; events",80,-2.8,2.8);      
  _h_gen_bhad_pt   = fs->make<TH1F>("gen_bhad_pt","gen b-hadron p_{T}; gen p_{T}^{b had} (GeV); events",150.,0.,150); 	       
  			       
  _h_reco_Z_mass     = fs->make<TH1F>("reco_Z_mass","reco Z mass; reco m_{ll} (GeV); events",100,50.,150.);		       
  _h_reco_Zdau1_eta  = fs->make<TH1F>("reco_Z_dau1_eta","reco Z daughter #eta; reco #eta_{l}; events",80,-2.8,2.8);      
  _h_reco_Zdau2_eta  = fs->make<TH1F>("reco_Z_dau2_eta","reco Z daughter #eta; reco #eta_{l}; events",80,-2.8,2.8);       
  _h_reco_Zdau1_pt   = fs->make<TH1F>("reco_Z_dau1_pt","reco Z daughter p_{T}; reco p_{T}^{l} (GeV); events",150.,0.,150);      
  _h_reco_Zdau2_pt   = fs->make<TH1F>("reco_Z_dau2_pt","reco Z daughter p_{T}; reco p_{T}^{l} (GeV); events",150.,0.,150);      
  			       
  _h_reco_bjet_eta        = fs->make<TH1F>("reco_bjet_eta","reco b-jet #eta; reco #eta_{b jet}; events",80,-2.8,2.8);     
  _h_reco_bjet_pt         = fs->make<TH1F>("reco_bjet_pt","reco b-jet p_{T}; reco p_{T}^{b jet} (GeV); events",150.,0.,150);       
  _h_bjet_energyresponse  = fs->make<TH1F>("reco_bjet_energyresponse","reco b-jet p_{T}/gen b-jet p_{T}; p_{T}^{reco}/p_{T}^{gen} for b-jets; events",100,0.,2.);

  // ==================== Summary Ntuple ============================
  
  if(saveNTuple_){
    myNTuple_ = fs->make<TTree>("T","Data for Acceptance");
    myNTuple_->Branch("isRightFlavour",&(myBools_.is_right_flavour_),"isRightFlavour/O");
    myNTuple_->Branch("isgenZyes",&(myBools_.is_gen_Z_yes_),"isgenZyes/O");
    myNTuple_->Branch("isgenjetyes",&(myBools_.is_gen_j_yes_),"isgenjetyes/O");
    myNTuple_->Branch("isgenbyes",&(myBools_.is_gen_b_yes_),"isgenbyes/O");
    myNTuple_->Branch("isgenZkinyes",&(myBools_.is_gen_Zkin_yes_),"isgenZkinyes/O");
    myNTuple_->Branch("isrecZkinyes",&(myBools_.is_rec_Zkin_yes_),"isrecZkinyes/O");
    myNTuple_->Branch("isrecbyes",&(myBools_.is_rec_b_yes_),"isrecbyes/O");
    myNTuple_->Branch("isrecjetyes",&(myBools_.is_rec_j_yes_),"isrecjetyes/O");
    myNTuple_->Branch("isreclepidiso_yes",&(myBools_.is_rec_lep_idiso_yes_),"isreclepidiso_yes/O");
    myNTuple_->Branch("isrecbHE_yes",&(myBools_.is_rec_b_HE_yes_),"isrec_b_HE_yes/O");
    myNTuple_->Branch("isrecbHP_yes",&(myBools_.is_rec_b_HP_yes_),"isrec_b_HP_yes/O");
    myNTuple_->Branch("isrecbCSVM_yes",&(myBools_.is_rec_b_CSVM_yes_),"isrec_b_CSVM_yes/O");
    myNTuple_->Branch("isrecbCSVT_yes",&(myBools_.is_rec_b_CSVT_yes_),"isrec_b_CSVT_yes/O");
    myNTuple_->Branch("evtweight",&(myBools_.evtWeight_),"evtweight/F");
    myNTuple_->Branch("evtweightSq",&(myBools_.evtWeightSq_),"evtweightSq/F");
    myNTuple_->Branch("lepIdIsoWeight",&(myBools_.lepIdIsoWeight_),"lepIdIsoWeight/F");
    myNTuple_->Branch("lepIdIsoWeightSq",&(myBools_.lepIdIsoWeightSq_),"lepIdIsoWeightSq/F");
    myNTuple_->Branch("bTagWeightHE",&(myBools_.bTagWeightSSVHE_)," bTagWeightHE/F");  
    myNTuple_->Branch("bTagWeightHESq",&(myBools_.bTagWeightSSVHESq_),"bTagWeightHESq/F");
    myNTuple_->Branch("bTagWeightHP",&(myBools_.bTagWeightSSVHP_ )," bTagWeightHP/F");  
    myNTuple_->Branch("bTagWeightHPSq",&(myBools_.bTagWeightSSVHPSq_),"bTagWeightHPSq/F");
    myNTuple_->Branch("bTagWeightCSVM",&(myBools_.bTagWeightCSVM_)," bTagWeightCSVM/F");  
    myNTuple_->Branch("bTagWeightCSVMSq",&(myBools_.bTagWeightCSVMSq_),"bTagWeightCSVMSq/F");
    myNTuple_->Branch("bTagWeightCSVT",&(myBools_.bTagWeightCSVT_ )," bTagWeightCSVT/F");  
    myNTuple_->Branch("bTagWeightCSVTSq",&(myBools_.bTagWeightCSVTSq_),"bTagWeightCSVTSq/F");
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
SimpleCFAnalyzer::endJob() 
{
}

// ------------ method called to produce the best Z candidate -----------------------
reco::CompositeCandidate 
SimpleCFAnalyzer::theBestZCand(std::vector<const reco::Candidate*> theLeptons){

  reco::CompositeCandidate theBestZCand;

  std::vector<const reco::Candidate*>   sortedLeptons = sortLeptonsByPt(theLeptons);
  std::vector<reco::CompositeCandidate> unsortedCands;
  
  // prepare vector of unsorted composite candidates
  for (std::vector<const reco::Candidate*>::const_iterator lBegin = sortedLeptons.begin(),lEnd = sortedLeptons.end(),ilept = lBegin;ilept != lEnd; ++ilept ) {
    for ( std::vector<const reco::Candidate*>::const_iterator jlept = ilept + 1; jlept != lEnd; ++jlept ) {
      if ( (*ilept)->charge() * (*jlept)->charge() < 0 ) {
	
	reco::CompositeCandidate ZCand;
	ZCand.addDaughter( **ilept,"l1");
	ZCand.addDaughter( **jlept,"l2");
	
	AddFourMomenta addp4;
	addp4.set(ZCand);      
	
	unsortedCands.push_back(ZCand);
      } // if opposite charge
    } // loop on second muon
  } // loop on first muon
  
  if(unsortedCands.size() >0){
    theBestZCand  = sortCandidatesByDifference(unsortedCands)[0];
  }

  return theBestZCand;

}

// ------------ method called to sort the leptons in pT -----------------------------
std::vector<const reco::Candidate*>
SimpleCFAnalyzer::sortLeptonsByPt(std::vector<const reco::Candidate*> leptons_)
{
  
  std::vector<const reco::Candidate*> sortedLeptons_ = leptons_;
  std::vector<double> ptParticles_;

  for(std::vector<const reco::Candidate*>::const_iterator it = leptons_.begin(); it!=leptons_.end(); it++){
    ptParticles_.push_back((*it)->pt());
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

// ------------ method called to Z candidates by difference abs(Mll-MZ) --------------
std::vector<reco::CompositeCandidate> 
SimpleCFAnalyzer::sortCandidatesByDifference(std::vector<reco::CompositeCandidate> unsortedCands)
{


  std::vector<reco::CompositeCandidate> sortedCands = unsortedCands;
  Double_t mZ=91.1876;
  std::vector<Double_t> diffZmass; 
  unsigned int ZCandSize=unsortedCands.size();
  
  if(unsortedCands.size()>1){
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
  } 
  return sortedCands;
}

// ------------ method called to tag b-hadrons -------------------------------------
bool
SimpleCFAnalyzer::hasBottom(const reco::Candidate &c) 
{
  bool tmpHasBottom = false;
  bool code1 =  ( abs(c.pdgId())> 500 && abs(c.pdgId())< 600);
  bool code2 =  ( abs(c.pdgId())> 5000 && abs(c.pdgId())< 6000);
  if ( code1 || code2 ) tmpHasBottom = true;
  return tmpHasBottom;
}

// ------------ method called to match any type of two objects ---------------------
template<class T, class U> 
std::vector< std::pair<T,U> > SimpleCFAnalyzer::matchByDR(std::vector<T> const& c1,
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
      if(dRCompare(pairVec[i],pairVec[j])){
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

// ------------ method called to match any type of two objects ---------------------
template<class T, class U> 
bool SimpleCFAnalyzer::dRCompare(std::pair<T,U> p1, std::pair<T,U> p2) {
  double dR1 = ROOT::Math::VectorUtil::DeltaR((*p1.first).momentum(),(*p1.second).momentum());
  double dR2 = ROOT::Math::VectorUtil::DeltaR((*p2.first).momentum(),(*p2.second).momentum());
  return (dR1 > dR2);
}

// ------------ method called to make lepton offline cuts --------------------------
bool SimpleCFAnalyzer::isTightZCandidate(reco::CompositeCandidate ZCand, const reco::BeamSpot& beamSpot, bool isMuChannel){

  bool istightZcandidate = false;
  const reco::Candidate* lep0 = ZCand.daughter(0);
  const reco::Candidate* lep1 = ZCand.daughter(1);

  if(isMuChannel){
    const pat::Muon* muon0 = dynamic_cast<const pat::Muon*>(&(*lep0));
    const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(&(*lep1));    

    if( ( muon0->isGlobalMuon() && muon0->isTrackerMuon() && 
	  muon0->globalTrack()->normalizedChi2() < 10 && 
	  muon0->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 &&                      
	  muon0->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&                         
	  muon0->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&                         
	  fabs(muon0->dB()) < 0.02 &&                                                                     
	  ((muon0->trackIso() + muon0->caloIso()) <0.15*muon0->pt()) &&                                             
	  muon0->numberOfMatches() > 1
	  ) &&
	( muon1->isGlobalMuon() && muon1->isTrackerMuon() && 
	  muon1->globalTrack()->normalizedChi2() < 10 && 
	  muon1->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 &&                      
	  muon1->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&                         
	  muon1->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&                         
	  fabs(muon1->dB()) < 0.02 &&                                                                     
	  ((muon1->trackIso() + muon1->caloIso()) <0.15*muon1->pt()) &&                                             
	  muon1->numberOfMatches() > 1 
	  )
	){ 
      istightZcandidate=true;
    }
  } else {
    
    const pat::Electron* ele0 = dynamic_cast<const pat::Electron*>(&(*lep0));
    const pat::Electron* ele1 = dynamic_cast<const pat::Electron*>(&(*lep1));
    
    if( (fabs(ele0->gsfTrack()->dxy(beamSpot))<0.02 && (ele0->electronID("eidVBTFRel85") == 7) ) &&    
	(fabs(ele1->gsfTrack()->dxy(beamSpot))<0.02 && (ele1->electronID("eidVBTFRel85") == 7) ) ){
      istightZcandidate=true;
    }
  }

  
  return istightZcandidate;
}

// ------------ method called to calculate muon scale factors --------------------
Double_t 
SimpleCFAnalyzer::getMuonScaleFactor(reco::CompositeCandidate &theRecoZcand){

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

// ------------ method called to calculate electron scale factors --------------------
Double_t 
SimpleCFAnalyzer::getElectronScaleFactor(reco::CompositeCandidate &theRecoZcand){
  
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
DEFINE_FWK_MODULE(SimpleCFAnalyzer);
