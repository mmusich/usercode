#include <utility>
#include <map>
#include <list>
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <time.h>

#include "TFile.h"
#include "TH1.h"
#include <Math/VectorUtil.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"
#include "RecoBTag/PerformanceDB/interface/BtagPerformance.h" 
#include "RecoBTag/Records/interface/BTagPerformanceRecord.h"

#include "ZbbAnalysis/AnalysisStep/interface/ZbbUtils.h"
#include "ZbbAnalysis/AnalysisStep/interface/AcceptanceCuts.h"
#include "ZbbAnalysis/AnalysisStep/interface/EventCategory.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbHistos.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbMultiLhistos.h"
#include "ZbbAnalysis/AnalysisStep/interface/VtxHistos.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbTypeDefs.h"
#include "ZbbAnalysis/AnalysisStep/interface/VtxAssociatorsUtils.h"

#include "ZbbAnalysis/AnalysisStep/interface/LumiReWeighting.h"
#include "ZbbAnalysis/AnalysisStep/interface/PU.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

class ZbbEventContentAnalyzer : public edm::EDAnalyzer {

public:
  // default constructor
  explicit ZbbEventContentAnalyzer(const edm::ParameterSet&);
  // default destructor
  ~ZbbEventContentAnalyzer();
  
private:
  // everything that needs to be done before the event loop
  virtual void beginJob();

  // everything that needs to be done during the event loop
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  
  // check if histogram was booked
  bool booked(const std::string histName) const { return histContainer_.find(histName.c_str())!=histContainer_.end(); };
  
  // fill histogram if it had been booked before
  void fill(const std::string histName, double value, double w) const { if(booked(histName.c_str())) histContainer_.find(histName.c_str())->second->Fill(value,w); };

  // fill trend histograms of number of event as a function of run number
  void fillTrendPlot(EventCategory& iEvtCategory, unsigned int theRunNumber_,Bool_t isMuChannel);

  // decide the category of the event
  void setEventCategory(const edm::Event& iEvent, const edm::EventSetup& iSetup, EventCategory& iEvtCategory, Bool_t isMuChannel, Bool_t AndOrSwitch) ;

  // check if trigger is ok
  Bool_t isTriggerOK(pat::TriggerEvent theTriggerEvent,UInt_t theRunNumber, std::map<std::vector<std::string>,std::vector<uint> > theAssociativeMap_);
   
  // check if trigger match is ok  
  Bool_t isTriggerMatchOK(const reco::CompositeCandidate& ZCand,UInt_t theRunNumber, std::map<std::vector<std::string>,std::vector<uint> > theAssociativeMap_, Bool_t isMuChannel, Bool_t AndOrSwitch);

  // print event info
  void PrintEvent(const edm::Event& iEvent,Int_t category, Bool_t verbose);
  
  // for Lumi Calculation
  void endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup); 

  /// everything that needs to be done after the event loop
  virtual void endJob() ;
  
  // simple map to contain all histograms; histograms are booked in the beginJob() method
  std::map<std::string,TH1F*> histContainer_; 
  std::vector<TH1F*> h_ptJetsCorrected;
  
  // utils
  ofstream outfile_;
  TString  outfilename_;

  //bits of selection
  double wasRun,wasTrigOk,wasLL,wasLLTight,wasLLTightAndTrigMatch,wasZLL,wasGoodJet,wasJetVertexAssoc,wasJetBTag,wasJetBTagExcl,wasJetBTagPlusMET,wasJetBTagPlusl,wasJetBTagPlus2l,wasJetBTagMCMatched;

  // switch to run all the plots
  bool doAllThePlotting_;

  //mc switch
  bool ismc_;

  //b-tag algo and WP
  std::string bTagAlgoWP_;

  //defining acceptance cuts
  AcceptanceCuts lCuts_;
  bool unLockDefaults_; // to unlock default cuts
  double jetEtaCut_,muonEtaCut_,eleEtaCut_,jetPtCut_,muonPtCut_,elePtCut_,betaCut_,betaStarCut_;

  //debug switch
  bool debug_;

  //switches trigger selection (true = OR | double triggers , false = AND | single object triggers) 
  bool andOr_;
  
  // input tags 
  edm::InputTag genPSrc_;
  edm::InputTag genJetSrc_;
  edm::InputTag metSrc_;
  edm::InputTag elecSrc_;
  edm::InputTag muonSrc_;
  edm::InputTag jetSrc_;
  edm::InputTag tauSrc_;
  edm::InputTag photonSrc_;
  edm::InputTag theVertexSrc_;
  edm::InputTag theZVertexSrc_;
  edm::InputTag theZmmSrc_;
  edm::InputTag theZeeSrc_;
  edm::InputTag theTriggerEventSrc_;

  // vpsets for mu and ele trigger maps
  edm::VParameterSet muhltrunrangeselection_;
  edm::VParameterSet elehltrunrangeselection_;  

  // correction levels for pat jet
  std::vector<std::string> str_corrLevels;
  std::vector<std::string> corrections_;
  Int_t corrsize;
  
  // MC weight from GenEventInfoProduct
  bool applymcweight_;
  bool isMCatNLO_;
  double theMCweight_;

  // pileup corrections
  bool applypucorrection_;
  edm::LumiReWeighting* theLumiWeights_;

  // vpsets for effb (MC) maps	 
  edm::VParameterSet effbcMC_ptrange_selection_;	
  // associative map for b-tag efficiency
  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffbMC_;
  std::map<EffMCpair,std::vector<double> > theAssociativeMapEffcMC_;

  // lepton efficiency 
  bool applyLeptonEff_;
  
  // calibration for b efficiency
  bool applybeffcalib_;  
  std::string beffcalibmethod_; 
  std::string bmistagcalibmethod_; 
  int theSFkFactor_,theLEffkFactor_;

  // definition of final state
  int  minBtags_;
  
  // associative maps for trigger-run ranges
  std::map<std::vector<std::string>,std::vector<uint> > theAssociativeMapMu_;
  std::map<std::vector<std::string>,std::vector<uint> > theAssociativeMapEle_;

  // evtcategory structs
  EventCategory theMMEventCategory_;
  EventCategory theEEEventCategory_;

  // histogram containers (for Z+b candidates)
  ZbbCandidate* ZbbCandidate_evtcat8;
  ZbbCandidate* ZbbCandidate_evtcat9;
  ZbbCandidate* ZbbCandidate_evtcat10;
  ZbbCandidate* ZbbCandidate_evtcat11;

  //n-tuple for Z+b candidates
  ZbbCandidateNtuple* ZbbCandidateNtuple_evtcat8;
  //Z_info Z; jet1_info jet1; jet2_info jet2; MET_info MET;
  //TTree ZbbNtuple4JEC;

  // histogram containers (for Z+b candidates, matched)
  ZbbMatchedCandidate* ZbbMatchedCandidate_evtcat12;

  // histogram containers (for intermediate levels)
  ZbbCollections* ZbbCollections_evtcat4;
  ZbbCollections* ZbbCollections_evtcat5;
  ZbbCollections* ZbbCollections_evtcat6;
  ZbbCollections* ZbbCollections_evtcat7;

  // histograms containers (for ZbbAllLeptons)
  ZbbAllLeptons* ZbbAllLeptons_evtcat4;
  ZbbAllLeptons* ZbbAllLeptons_evtcat5;
  ZbbAllLeptons* ZbbAllLeptons_evtcat6;
  ZbbAllLeptons* ZbbAllLeptons_evtcat8;
  ZbbAllLeptons* ZbbAllLeptons_evtcat9;
  ZbbAllLeptons* ZbbAllLeptons_evtcat10;
  ZbbAllLeptons* ZbbAllLeptons_evtcat11;

  // histograms containers (for ZbbExtraLeptons)
  ZbbExtraLeptons* ZbbExtraLeptons_evtcat4;
  ZbbExtraLeptons* ZbbExtraLeptons_evtcat5;
  ZbbExtraLeptons* ZbbExtraLeptons_evtcat6;
  ZbbExtraLeptons* ZbbExtraLeptons_evtcat8;
  ZbbExtraLeptons* ZbbExtraLeptons_evtcat9;
  ZbbExtraLeptons* ZbbExtraLeptons_evtcat10;
  ZbbExtraLeptons* ZbbExtraLeptons_evtcat11;
  
  // histogram containers (for object components)
  ZbbBasicComponents* ZbbBasicComponents_evtcat4;
  ZbbBasicComponents* ZbbBasicComponents_evtcat5;
  ZbbBasicComponents* ZbbBasicComponents_evtcat6;
  ZbbBasicComponents* ZbbBasicComponents_evtcat8;
  
  // histogram containers (for MC truth) 
  // If not booked for it crashes Execute.C 
  ZbbMCinfo* ZbbMCinfo_evtcat0;
  ZbbMCinfo* ZbbMCinfo_evtcat5;
  ZbbMCinfo* ZbbMCinfo_evtcat6;
  ZbbMCinfo* ZbbMCinfo_evtcat8;
  
  //Lepton categories
  MCTruthLFromHF *MCTruthL_Mu_all;
  MCTruthLFromHF *MCTruthL_Mu_fromHF;
  MCTruthLFromHF *MCTruthL_Mu_fromB;
  MCTruthLFromHF *MCTruthL_Mu_fromC;
  MCTruthLFromHF *MCTruthL_Mu_fromBC;

  MCTruthLFromHF *MCTruthL_Ele_all;
  MCTruthLFromHF *MCTruthL_Ele_fromHF;
  MCTruthLFromHF *MCTruthL_Ele_fromB;
  MCTruthLFromHF *MCTruthL_Ele_fromC;
  MCTruthLFromHF *MCTruthL_Ele_fromBC;

  //Lepton from HF categories
  MCTruthLFromHF *MCTruthL_evtcat4;
  MCTruthLFromHF *MCTruthL_evtcat5;
  MCTruthLFromHF *MCTruthL_evtcat6;
  MCTruthLFromHF *MCTruthL_evtcat8;
  MCTruthLFromHF *MCTruthL_evtcat10;
  MCTruthLFromHF *MCTruthL_evtcat11;

  //histograms for ABCD 
  ZbbABCDMatrix*  ZbbABCDMatrix_evtcat0;
  ZbbABCDMatrix*  ZbbABCDMatrix_evtcat5;
  ZbbABCDMatrix*  ZbbABCDMatrix_evtcat6;
  ZbbABCDMatrix*  ZbbABCDMatrix_evtcat7;
  ZbbABCDMatrix*  ZbbABCDMatrix_evtcat8;
  
  // vertex quantities histos
  VtxHistos* VtxHistos_evcat5;
  VtxHistos* VtxHistos_evcat6;
  VtxHistos* VtxHistos_evcat7;
  VtxHistos* VtxHistos_evcat8;
  
  // string for trend plots
  TString evtVsRunMuName[12],evtVsRunEleName[12],evtVsRunAllName[12];
  TString xsecVsRunMuName[12],xsecVsRunEleName[12],xsecVsRunAllName[12];

  //++++++++++++++++++++++++++++++++++++++++++++++
  // For luminosity plots
  //++++++++++++++++++++++++++++++++++++++++++++++

  // Luminosity studies
  std::string theLumiSummaryTag_; 
  //std::vector<std::pair<uint, uint> > runRangeList_ ;
  std::map<int,std::pair<int,double> > reclumibyrun_;

  int totLS_;
  double totLumiRecorded_, totLumiDelivered_;
  int invalidLS_totLS_;
  double invalidLS_totLumiRecorded_, invalidLS_totLumiDelivered_;

};

ZbbEventContentAnalyzer::ZbbEventContentAnalyzer(const edm::ParameterSet& iConfig):
  histContainer_(),
  outfilename_(iConfig.getParameter<std::string>("OutfileName")), 
  doAllThePlotting_(iConfig.getParameter<bool>("doAllThePlotting")),
  ismc_(iConfig.getParameter<bool>("isMC")),
  bTagAlgoWP_(iConfig.getParameter<std::string>("bTagAlgoWP")),
  unLockDefaults_(iConfig.getParameter<bool>("unLockDefaultCuts")),
  jetEtaCut_(iConfig.getParameter<double>("jetEtaCut")),
  muonEtaCut_(iConfig.getParameter<double>("muonEtaCut")),
  eleEtaCut_(iConfig.getParameter<double>("eleEtaCut")),
  jetPtCut_(iConfig.getParameter<double>("jetPtCut")),
  muonPtCut_(iConfig.getParameter<double>("muonPtCut")),
  elePtCut_(iConfig.getParameter<double>("elePtCut")),
  betaCut_(iConfig.getParameter<double>("betaCut")),
  betaStarCut_(iConfig.getParameter<double>("betaStarCut")), 
  debug_(iConfig.getParameter<bool>("Debug")),  
  andOr_(iConfig.getParameter<bool>("andOr")),
  genPSrc_(iConfig.getUntrackedParameter<edm::InputTag>("genPSrc")),
  genJetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("genJetSrc")),
  metSrc_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc")),        
  elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc")),
  muonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc")),
  jetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc" )),
  theVertexSrc_(iConfig.getUntrackedParameter<edm::InputTag>("VertexSrc")),
  theZVertexSrc_(iConfig.getUntrackedParameter<edm::InputTag>("ZVertexSrc")),
  theZmmSrc_(iConfig.getUntrackedParameter<edm::InputTag>("ZmmSrc")),
  theZeeSrc_(iConfig.getUntrackedParameter<edm::InputTag>("ZeeSrc")),
  theTriggerEventSrc_(iConfig.getUntrackedParameter<edm::InputTag>("TriggerEventSrc")),
  muhltrunrangeselection_(iConfig.getParameter<edm::VParameterSet>("muHLTRunRangeList")),
  elehltrunrangeselection_(iConfig.getParameter<edm::VParameterSet>("eleHLTRunRangeList")),
  applymcweight_(iConfig.getParameter<bool>("applyMCweight")),
  isMCatNLO_(iConfig.getParameter<bool>("isMCatNLO")),
  applypucorrection_(iConfig.getParameter<bool>("applyPUcorrection")),
  effbcMC_ptrange_selection_(iConfig.getParameter<edm::VParameterSet>("EffbcMCPtRangeList")),
  applyLeptonEff_(iConfig.getParameter<bool>("applyLeptonEfficiency")),  
  applybeffcalib_(iConfig.getParameter<bool>("applybEffCalibration")),  
  beffcalibmethod_(iConfig.getParameter<std::string>("bEffCalibrationMethod")),
  bmistagcalibmethod_(iConfig.getParameter<std::string>("bMistagCalibrationMethod")), 
  theSFkFactor_(iConfig.getParameter<int>("bSFkFactor")),
  theLEffkFactor_(iConfig.getParameter<int>("lEffkFactor")),
  minBtags_(iConfig.getParameter<int>("minBtags"))
{
  //intialize the associative maps

  // muon triggers
  for(edm::VParameterSet::const_iterator pset1 = muhltrunrangeselection_.begin(); pset1!=muhltrunrangeselection_.end(); pset1++){
    std::vector<std::string> theString_ = pset1->getUntrackedParameter<std::vector<std::string> >("hltpath");
    std::vector<uint> theRunRange_ = pset1->getUntrackedParameter<std::vector<uint> >("runrange");
    theAssociativeMapMu_[theString_]=theRunRange_; 
  }
  
  // electron triggers
  for(edm::VParameterSet::const_iterator pset2 = elehltrunrangeselection_.begin(); pset2!=elehltrunrangeselection_.end(); pset2++){
    std::vector<std::string> theString_ = pset2->getUntrackedParameter<std::vector<std::string> >("hltpath");
    std::vector<uint> theRunRange_ = pset2->getUntrackedParameter<std::vector<uint> >("runrange");
    theAssociativeMapEle_[theString_]=theRunRange_; 
  }

  //initialize the MC weight
  theMCweight_=1;

  // b-tag MC efficiency 
  EffMCpair theEffb_;
  EffMCpair theEffc_;
  for(edm::VParameterSet::const_iterator pset3 = effbcMC_ptrange_selection_.begin(); pset3 != effbcMC_ptrange_selection_.end(); pset3++){
    std::vector<double> thePtRange_ = pset3->getUntrackedParameter<std::vector<double> >("ptrange");
    double theEffb_barrel_  = pset3->getUntrackedParameter<double>("effb_barrel");
    double theEffb_forward_ = pset3->getUntrackedParameter<double>("effb_forward");
    double theEffc_barrel_  = pset3->getUntrackedParameter<double>("effc_barrel");
    double theEffc_forward_ = pset3->getUntrackedParameter<double>("effc_forward");
    theEffb_ = std::make_pair(theEffb_barrel_,theEffb_forward_);
    theEffc_ = std::make_pair(theEffc_barrel_,theEffc_forward_);
    theAssociativeMapEffbMC_[theEffb_]=thePtRange_; 
    theAssociativeMapEffcMC_[theEffc_]=thePtRange_; 
  }
  //initialize the lumiweight
  theLumiWeights_=0;

  //setting the offline cuts
  lCuts_.setDefault();
  
  if(unLockDefaults_){
    lCuts_.clear();
    lCuts_.set(jetEtaCut_,jetPtCut_,muonEtaCut_,muonPtCut_,eleEtaCut_,elePtCut_,betaCut_,betaStarCut_);
  }

  // Preparing the names of trend histograms
  char a[128],b[128],c[128],d[128],e[128],f[128];
  for(Int_t nCut=0;nCut<12;nCut++){

    sprintf(a,"evtsVsRunMu_after%i_cut",nCut);
    sprintf(b,"evtsVsRunEle_after%i_cut",nCut);
    sprintf(c,"evtsVsRunAll_after%i_cut",nCut);
    evtVsRunMuName[nCut] =a;
    evtVsRunEleName[nCut]=b;
    evtVsRunAllName[nCut]=c;
   
    sprintf(d,"xsecVsRunMu_after%i_cut",nCut);
    sprintf(e,"xsecVsRunEle_after%i_cut",nCut);
    sprintf(f,"xsecVsRunAll_after%i_cut",nCut);
    xsecVsRunMuName[nCut]=d;
    xsecVsRunEleName[nCut]=e;
    xsecVsRunAllName[nCut]=f;
    
  }
  
  //setting the lumiproducer
  theLumiSummaryTag_ = "lumiProducer::RECO";

}

ZbbEventContentAnalyzer::~ZbbEventContentAnalyzer()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//ZbbEventContentAnalyzer::Analyze
///////////////////////////////////////////////////////////////////////////////////////////////////

void
ZbbEventContentAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 
  wasRun++;

  // get the run number
  unsigned int theRunNumber_=iEvent.id().run(); 

  // get met collection  
  edm::Handle<edm::View<pat::MET> > mets;
  iEvent.getByLabel(metSrc_,mets);

  // get electron collection
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByLabel(elecSrc_,electrons);

  // get muon collection
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(muonSrc_,muons);

  // get jet collection
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jetSrc_,jets);
  
  // get vertex collection
  edm::Handle<edm::View<reco::Vertex> > vertexCollection;
  iEvent.getByLabel(theVertexSrc_, vertexCollection);
  edm::View<reco::Vertex> vertices = *(vertexCollection.product());
  unsigned int vertexCollectionSize = vertexCollection->size();  
  int nvvertex = 0;
  for (unsigned int i=0; i<vertexCollectionSize; i++) {
    const reco::Vertex& vertex = vertexCollection->at(i);
    if (vertex.isValid()) nvvertex++;
  }

  // get Z vertex collection
  edm::Handle<edm::View<reco::Vertex> > ZvertexCollection;
  iEvent.getByLabel(theZVertexSrc_, ZvertexCollection);
  const reco::Vertex& theZVertexFromProd = ZvertexCollection->at(0);

  // get the Z->mm collection
  edm::Handle<reco::CompositeCandidateCollection> zmmHandle;
  iEvent.getByLabel(theZmmSrc_, zmmHandle);
  
  // get the Z->ee collection
  edm::Handle<reco::CompositeCandidateCollection> zeeHandle;
  iEvent.getByLabel(theZeeSrc_, zeeHandle);
 
  //#################################
  // fills eventCategory
  //#################################

  // clear the struct
  theMMEventCategory_.reset();
  theEEEventCategory_.reset();
  
  // fill the struct
  setEventCategory(iEvent,iSetup,theMMEventCategory_,true,andOr_);
  setEventCategory(iEvent,iSetup,theEEEventCategory_,false,andOr_);
  
  //#################################
  // resolve ambiguity if there are two dileptons
  //#################################

  if(theMMEventCategory_.isLL() && theEEEventCategory_.isLL()){
    
    const reco::CompositeCandidateCollection & zmm = *(zmmHandle.product());  
    const reco::CompositeCandidateCollection & zee = *(zeeHandle.product());
    
    std::vector<reco::CompositeCandidate> unsortedMixedFlavourCands;
    
    reco::CompositeCandidate ZmmCand =ZbbUtils::sortCandidatesByDifference(zmm)[0];
    reco::CompositeCandidate ZeeCand =ZbbUtils::sortCandidatesByDifference(zee)[0];
    
    unsortedMixedFlavourCands.push_back(ZmmCand);
    unsortedMixedFlavourCands.push_back(ZeeCand);
    
    reco::CompositeCandidate ZCand = ZbbUtils::sortCandidatesByDifference(unsortedMixedFlavourCands)[0];   
    
    const reco::Candidate* lepton1 = ZCand.daughter(0);
    const reco::Candidate* lepton2 = ZCand.daughter(1);

    if(lepton1->isMuon() && lepton2->isMuon()){
      theEEEventCategory_.LL_=false;
    } else if(lepton1->isElectron() && lepton2->isElectron()) {
      theMMEventCategory_.LL_=false;
    } else {
      std::cout<<"ZbbEventContentAnalyzer::analyze: Shouldn't be here!"<<std::endl;
      return;
    }
  }

  //#################################
  // get the event weight   
  //#################################

  Double_t weight=1;
  Double_t bTagEffWeightCF=1;
  
  if (theMMEventCategory_.isLLTight()) {
    weight=theMMEventCategory_.getWeight();
    bTagEffWeightCF=theMMEventCategory_.getBTagEffWeightCF();
  } else if(theEEEventCategory_.isLLTight()) {
    weight=theEEEventCategory_.getWeight();
    bTagEffWeightCF=theEEEventCategory_.getBTagEffWeightCF();
  }
  
  //#################################
  // Get vertex infos for PU reweight
  //#################################

  // fill the number of vertices (w & w/o reweighted) histo
  fill("nVvertex",nvvertex,1);
  fill("nVvertexReweight",nvvertex,weight);
  
  // Get the pile-up summary info (for MC)
  if(ismc_ && applypucorrection_){
    edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel("addPileupInfo", PupInfo);
    
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    
    // list of available methods: 
    // int num_PU_vertices_;                // the number of pileup interactions that have been added to the event
    // int bunchCrossing_;                  // to which bunch crossing does this interaction belong?  New in 3_11_3 and 4_1_3 or later
    // std::vector<float> zpositions_;      // the true primary vertex position along the z axis for each added interaction
    // std::vector<float> sumpT_lowpT_;     // the sum of the transverse momentum of the tracks originating from each interaction, where track pT > low_cut
    // std::vector<float> sumpT_highpT_;    // the sum of the transverse momentum of the tracks originating from each interaction, where track pT > high_cut
    // std::vector<int> ntrks_lowpT_;       // the number of tracks originating from each interaction, where track pT > low_cut
    // std::vector<int> ntrks_highpT_;      // the number of tracks originating from each interaction, where track pT > high_cut
    
    int npv = -1;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(debug_) std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
      fill("nBunchCrossing",PVI->getBunchCrossing(),1);     
      int BX = PVI->getBunchCrossing();
      if(BX == 0) { 
	npv = PVI->getPU_NumInteractions();
	continue;
      } 
    }
    fill("nPileUp",npv,1); 
    fill("nPileUpReweight",npv,weight);
  }
  
  //#################################
  // Fill trend plots for Data only
  //#################################

  if(!ismc_){
    fillTrendPlot(theMMEventCategory_,theRunNumber_,true);
    fillTrendPlot(theEEEventCategory_,theRunNumber_,false);
  }

  //#################################
  // Fill cut flow table
  //#################################

  // fill the event weight histo for all events
  fill("eventWeight",weight,1);
  // filled at every event
  fill("eventCategoryMu",0,weight);
  fill("UnwEventCategoryMu",0,1);

  // fill event category histogram for di-muon channel
  if(theMMEventCategory_.isTrigger()){
    wasTrigOk++;
    fill("eventCategoryMu",1,weight);
    fill("UnwEventCategoryMu",1,1);
    if(theMMEventCategory_.isLL()){
      wasLL++;
      fill("eventCategoryMu",2,weight);
      fill("UnwEventCategoryMu",2,1);
      if(theMMEventCategory_.isLLTight()){
	wasLLTight++;
	fill("eventCategoryMu",3,weight);
	fill("UnwEventCategoryMu",3,1);
	if(theMMEventCategory_.isLLTightAndTriggerMatched()){
	  wasLLTightAndTrigMatch++;
	  fill("eventCategoryMu",4,weight);
	  fill("UnwEventCategoryMu",4,1);
	  if(theMMEventCategory_.isZLL()){ 
	    wasZLL++;
	    fill("eventCategoryMu",5,weight);
	    fill("UnwEventCategoryMu",5,1);
	    if(theMMEventCategory_.isGoodJet()){
	      wasGoodJet++;
	      fill("eventCategoryMu",6,weight);
	      fill("UnwEventCategoryMu",6,1);
	      if(theMMEventCategory_.isJetVertexAssoc()){
		wasJetVertexAssoc++;
		fill("eventCategoryMu",7,weight);
		fill("UnwEventCategoryMu",7,1);
	      }
	      if(theMMEventCategory_.isJetBTagged()){
		wasJetBTag++;
		this->PrintEvent(iEvent,8,false);
		fill("eventCategoryMu",8,weight);
		fill("UnwEventCategoryMu",8,1);
		if(theMMEventCategory_.isJetBTaggedMCMatched()){
		  wasJetBTagMCMatched++;
		} // if the b-tagged jet is MC matched
		if(theMMEventCategory_.isExclusive()){
		  wasJetBTagExcl++;
		  fill("eventCategoryMu",9,weight*bTagEffWeightCF);
		  fill("UnwEventCategoryMu",9,1);
		} // is exclusive final state 
		if (theMMEventCategory_.isThreeLeptons()){
		  wasJetBTagPlusl++;
		  fill("eventCategoryMu",10,weight);
		  fill("UnwEventCategoryMu",10,1);
		  if (theMMEventCategory_.isFourLeptons()){
		    wasJetBTagPlus2l++;
		    fill("eventCategoryMu",11,weight);
		    fill("UnwEventCategoryMu",11,1);
		  } //two extra lepton
		} //one extra lepton
	      } // has at least a good jet b tagged
	    } // has at least a good jet
	  } // the di-lepton is a Z 
	} // di-lepton is matched with HLT object
      } // tight di-lepton
    } // has at least a di-lepton
  } // is trigger ok

  // filled at every event
  fill("eventCategoryEle",0,weight);
  fill("UnwEventCategoryEle",0,1);
  
  // fill event category histogram for di-electron channel
  if(theEEEventCategory_.isTrigger()){
    wasTrigOk++;
    fill("eventCategoryEle",1,weight);
    fill("UnwEventCategoryEle",1,1);
    if(theEEEventCategory_.isLL()){
      wasLL++;
      fill("eventCategoryEle",2,weight);
      fill("UnwEventCategoryEle",2,1);
      if(theEEEventCategory_.isLLTight()){
	wasLLTight++;	
	fill("eventCategoryEle",3,weight);
	fill("UnwEventCategoryEle",3,1);
	if(theEEEventCategory_.isLLTightAndTriggerMatched()){
	  wasLLTightAndTrigMatch++;	
	  fill("eventCategoryEle",4,weight);
	  fill("UnwEventCategoryEle",4,1);
	  if(theEEEventCategory_.isZLL()){
	    wasZLL++;
	    fill("eventCategoryEle",5,weight);
	    fill("UnwEventCategoryEle",5,1);
	    if(theEEEventCategory_.isGoodJet()){
	      wasGoodJet++;
	      fill("eventCategoryEle",6,weight);
	      fill("UnwEventCategoryEle",6,1);
	      if(theEEEventCategory_.isJetVertexAssoc()){
		wasJetVertexAssoc++;
		fill("eventCategoryEle",7,weight);
		fill("UnwEventCategoryEle",7,1);
	      }
	      if(theEEEventCategory_.isJetBTagged()){
		wasJetBTag++;
		fill("eventCategoryEle",8,weight);
		fill("UnwEventCategoryEle",8,1);
		this->PrintEvent(iEvent,8,false);
		if(theEEEventCategory_.isJetBTaggedMCMatched()){
		  wasJetBTagMCMatched++;
		} // if the b-tagged jet is MC matched
		if(theEEEventCategory_.isExclusive()){
		  wasJetBTagExcl++;
		  fill("eventCategoryEle",9,weight*bTagEffWeightCF);
		  fill("UnwEventCategoryEle",9,1);
		} // is exclusive final state 
		if (theEEEventCategory_.isThreeLeptons()){
		  wasJetBTagPlusl++;
		  fill("eventCategoryEle",10,weight);
		  fill("UnwEventCategoryEle",10,1);
		  if (theEEEventCategory_.isFourLeptons()){
		    wasJetBTagPlus2l++;
		    fill("eventCategoryEle",11,weight);
		    fill("UnwEventCategoryEle",11,1);
		  } //two extra lepton
		} //one extra lepton
	      } // has at least a good jet b tagged
	    } // has at least a good jet
	  } // the di-lepton is a Z 
	} // at least one lepton is matched with HLT object
      } // was a tight di-lepton
    } // has at least a di-lepton
  } // is trigger ok
    
  //-------------------------- Basic Components Plots ---------------------------------
  
  // fill muons
  if( doAllThePlotting_){
    ZbbBasicComponents_evtcat4->fillMu(theMMEventCategory_);
    ZbbBasicComponents_evtcat5->fillMu(theMMEventCategory_);
    ZbbBasicComponents_evtcat6->fillMu(theMMEventCategory_);
  }
  ZbbBasicComponents_evtcat8->fillMu(theMMEventCategory_);
 
  // fill electrons
  if( doAllThePlotting_){
    ZbbBasicComponents_evtcat4->fillEle(theEEEventCategory_);
    ZbbBasicComponents_evtcat5->fillEle(theEEEventCategory_);
    ZbbBasicComponents_evtcat6->fillEle(theEEEventCategory_);
  }
  ZbbBasicComponents_evtcat8->fillEle(theEEEventCategory_);

  //-------------------------- Collection Plots ---------------------------------------
  
  if( doAllThePlotting_){
    // fill muons
    ZbbCollections_evtcat4->fillMu(theMMEventCategory_);
    ZbbCollections_evtcat5->fillMu(theMMEventCategory_);
    ZbbCollections_evtcat6->fillMu(theMMEventCategory_);
    ZbbCollections_evtcat7->fillMu(theMMEventCategory_);

    // fill electrons
    ZbbCollections_evtcat4->fillEle(theEEEventCategory_);
    ZbbCollections_evtcat5->fillEle(theEEEventCategory_);
    ZbbCollections_evtcat6->fillEle(theEEEventCategory_);
    ZbbCollections_evtcat7->fillEle(theEEEventCategory_);
    
    //fill jets
    ZbbCollections_evtcat4->fillJets(theMMEventCategory_,theEEEventCategory_,jets,mets);
    ZbbCollections_evtcat5->fillJets(theMMEventCategory_,theEEEventCategory_,jets,mets);
    ZbbCollections_evtcat6->fillJets(theMMEventCategory_,theEEEventCategory_,jets,mets);
    ZbbCollections_evtcat7->fillJets(theMMEventCategory_,theEEEventCategory_,jets,mets); 
  }

  //-------------------------- Candidate Plots ---------------------------------------

  // fill muons
  ZbbCandidate_evtcat8->fillMu(theMMEventCategory_);
  ZbbCandidate_evtcat9->fillMu(theMMEventCategory_);
  ZbbCandidate_evtcat10->fillMu(theMMEventCategory_);
  ZbbCandidate_evtcat11->fillMu(theMMEventCategory_);
  
  // fill electrons
  ZbbCandidate_evtcat8->fillEle(theEEEventCategory_);
  ZbbCandidate_evtcat9->fillEle(theEEEventCategory_);
  ZbbCandidate_evtcat10->fillEle(theEEEventCategory_);
  ZbbCandidate_evtcat11->fillMu(theMMEventCategory_);

  //fill jets
  ZbbCandidate_evtcat8->fillJets(theMMEventCategory_,theEEEventCategory_,jets,mets,bTagAlgoWP_,theAssociativeMapEffbMC_,theAssociativeMapEffcMC_);
  ZbbCandidate_evtcat9->fillJets(theMMEventCategory_,theEEEventCategory_,jets,mets,bTagAlgoWP_,theAssociativeMapEffbMC_,theAssociativeMapEffcMC_);
  ZbbCandidate_evtcat10->fillJets(theMMEventCategory_,theEEEventCategory_,jets,mets,bTagAlgoWP_,theAssociativeMapEffbMC_,theAssociativeMapEffcMC_);
  ZbbCandidate_evtcat11->fillJets(theMMEventCategory_,theEEEventCategory_,jets,mets,bTagAlgoWP_,theAssociativeMapEffbMC_,theAssociativeMapEffcMC_);

  //-------------------------- Candidate N-tuple --------------------------------------

  //fill n-tuple
  ZbbCandidateNtuple_evtcat8->fill(theMMEventCategory_,theEEEventCategory_,jets,bTagAlgoWP_,mets, iEvent, nvvertex); //,genParticlesCollection,genJets);

  //-------------------------- All leptons Plots --------------------------------------
  
  //fill muons
  if( doAllThePlotting_){
    ZbbAllLeptons_evtcat4->fillMu(theMMEventCategory_,theEEEventCategory_,muons,jets,bTagAlgoWP_);
    ZbbAllLeptons_evtcat5->fillMu(theMMEventCategory_,theEEEventCategory_,muons,jets,bTagAlgoWP_);
    ZbbAllLeptons_evtcat6->fillMu(theMMEventCategory_,theEEEventCategory_,muons,jets,bTagAlgoWP_);
  }
  ZbbAllLeptons_evtcat8->fillMu(theMMEventCategory_,theEEEventCategory_,muons,jets,bTagAlgoWP_);
  ZbbAllLeptons_evtcat9->fillMu(theMMEventCategory_,theEEEventCategory_,muons,jets,bTagAlgoWP_);
  ZbbAllLeptons_evtcat10->fillMu(theMMEventCategory_,theEEEventCategory_,muons,jets,bTagAlgoWP_);
  ZbbAllLeptons_evtcat11->fillMu(theMMEventCategory_,theEEEventCategory_,muons,jets,bTagAlgoWP_);

  //fill electrons
  if( doAllThePlotting_){
    ZbbAllLeptons_evtcat4->fillEle(theMMEventCategory_,theEEEventCategory_,electrons,jets,bTagAlgoWP_);
    ZbbAllLeptons_evtcat5->fillEle(theMMEventCategory_,theEEEventCategory_,electrons,jets,bTagAlgoWP_);
    ZbbAllLeptons_evtcat6->fillEle(theMMEventCategory_,theEEEventCategory_,electrons,jets,bTagAlgoWP_);
  }
  ZbbAllLeptons_evtcat8->fillEle(theMMEventCategory_,theEEEventCategory_,electrons,jets,bTagAlgoWP_);
  ZbbAllLeptons_evtcat9->fillEle(theMMEventCategory_,theEEEventCategory_,electrons,jets,bTagAlgoWP_);
  ZbbAllLeptons_evtcat10->fillEle(theMMEventCategory_,theEEEventCategory_,electrons,jets,bTagAlgoWP_);
  ZbbAllLeptons_evtcat11->fillEle(theMMEventCategory_,theEEEventCategory_,electrons,jets,bTagAlgoWP_);

  //-------------------------- Extra leptons Plots ------------------------------------
  
  //fill leptons
  if( doAllThePlotting_){
    ZbbExtraLeptons_evtcat4->fillLeptons(theMMEventCategory_,theEEEventCategory_,muons,electrons,jets,bTagAlgoWP_);
    ZbbExtraLeptons_evtcat5->fillLeptons(theMMEventCategory_,theEEEventCategory_,muons,electrons,jets,bTagAlgoWP_);
    ZbbExtraLeptons_evtcat6->fillLeptons(theMMEventCategory_,theEEEventCategory_,muons,electrons,jets,bTagAlgoWP_);
  } 
  ZbbExtraLeptons_evtcat8->fillLeptons(theMMEventCategory_,theEEEventCategory_,muons,electrons,jets,bTagAlgoWP_);
  ZbbExtraLeptons_evtcat9->fillLeptons(theMMEventCategory_,theEEEventCategory_,muons,electrons,jets,bTagAlgoWP_);
  ZbbExtraLeptons_evtcat10->fillLeptons(theMMEventCategory_,theEEEventCategory_,muons,electrons,jets,bTagAlgoWP_);
  ZbbExtraLeptons_evtcat11->fillLeptons(theMMEventCategory_,theEEEventCategory_,muons,electrons,jets,bTagAlgoWP_);

  //-------------------------- ABCD Plots ----------------------------------------------
  if( doAllThePlotting_){
    ZbbABCDMatrix_evtcat0->fill(theMMEventCategory_,theEEEventCategory_,mets);
    ZbbABCDMatrix_evtcat5->fill(theMMEventCategory_,theEEEventCategory_,mets);
    ZbbABCDMatrix_evtcat6->fill(theMMEventCategory_,theEEEventCategory_,mets);
    ZbbABCDMatrix_evtcat7->fill(theMMEventCategory_,theEEEventCategory_,mets);
  }
  ZbbABCDMatrix_evtcat8->fill(theMMEventCategory_,theEEEventCategory_,mets);

  if( doAllThePlotting_){
    ZbbABCDMatrix_evtcat0->fillMu(theMMEventCategory_,mets);
    ZbbABCDMatrix_evtcat5->fillMu(theMMEventCategory_,mets);
    ZbbABCDMatrix_evtcat6->fillMu(theMMEventCategory_,mets);
    ZbbABCDMatrix_evtcat7->fillMu(theMMEventCategory_,mets);
  }
  ZbbABCDMatrix_evtcat8->fillMu(theMMEventCategory_,mets);

  if( doAllThePlotting_){
    ZbbABCDMatrix_evtcat0->fillEle(theEEEventCategory_,mets);
    ZbbABCDMatrix_evtcat5->fillEle(theEEEventCategory_,mets);
    ZbbABCDMatrix_evtcat6->fillEle(theEEEventCategory_,mets);
    ZbbABCDMatrix_evtcat7->fillEle(theEEEventCategory_,mets);
  }
  ZbbABCDMatrix_evtcat8->fillEle(theEEEventCategory_,mets);

  //-------------------------- Vertex Histos -------------------------------------------
  
  if( doAllThePlotting_){
    VtxHistos_evcat5->fill(theMMEventCategory_,theEEEventCategory_,vertices,jets,bTagAlgoWP_);
    VtxHistos_evcat6->fill(theMMEventCategory_,theEEEventCategory_,vertices,jets,bTagAlgoWP_);
    VtxHistos_evcat7->fill(theMMEventCategory_,theEEEventCategory_,vertices,jets,bTagAlgoWP_);
  }
  VtxHistos_evcat8->fill(theMMEventCategory_,theEEEventCategory_,vertices,jets,bTagAlgoWP_);

  //------------------------- Matched Candidate Plots (fake matching for DATA) ---------
  if(!ismc_) ZbbMatchedCandidate_evtcat12->fillDATA(theMMEventCategory_,theEEEventCategory_,jets,bTagAlgoWP_);

  //--------------------------  MC truth Plots -----------------------------------------

  if ( ismc_ ) {
    
     // get LHE Event
    edm::Handle<LHEEventProduct> lhevtCollection;
    lhevtCollection.clear();
    if( lhevtCollection.isValid () ) iEvent.getByType( lhevtCollection );
  
    edm::Handle<reco::GenParticleCollection> genParticlesCollection;
    iEvent.getByLabel(genPSrc_, genParticlesCollection);

    edm::Handle<reco::GenJetCollection> genJetsHandle;
    iEvent.getByLabel(genJetSrc_,genJetsHandle);
    const reco::GenJetCollection & genJets = *(genJetsHandle.product());
   

    //------------------------- Matched Candidate Plots ---------------------------------
    ZbbMatchedCandidate_evtcat12->fillMC(theMMEventCategory_,theEEEventCategory_,jets,bTagAlgoWP_,genParticlesCollection,genJets);

    //-- next fill the histos
    if(doAllThePlotting_){
      ZbbMCinfo_evtcat0->fill(theMMEventCategory_,theEEEventCategory_,bTagAlgoWP_,lhevtCollection,genParticlesCollection,genJets,jets);
      ZbbMCinfo_evtcat5->fill(theMMEventCategory_,theEEEventCategory_,bTagAlgoWP_,lhevtCollection,genParticlesCollection,genJets,jets);
      ZbbMCinfo_evtcat6->fill(theMMEventCategory_,theEEEventCategory_,bTagAlgoWP_,lhevtCollection,genParticlesCollection,genJets,jets);
    }
    ZbbMCinfo_evtcat8->fill(theMMEventCategory_,theEEEventCategory_,bTagAlgoWP_,lhevtCollection,genParticlesCollection,genJets,jets);

    //Histogram containers (MC Truth leptons from HF)

    //Muons
    for(edm::View<pat::Muon>::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon){
      if( muon->genParticleById(13,1).isNonnull() ) {
	const reco::Candidate *GenMuon = &(*muon->genParticleById(13,1));
	while (GenMuon->mother()!=0 && abs(GenMuon->mother()->pdgId()) == 13) GenMuon = GenMuon->mother();
	if (ZbbUtils::isPreselMu(*muon,theZVertexFromProd,lCuts_)){
	  if( doAllThePlotting_){
	    MCTruthL_evtcat4->fillMu(theMMEventCategory_,theEEEventCategory_,*muon,jets,bTagAlgoWP_);
	    MCTruthL_evtcat5->fillMu(theMMEventCategory_,theEEEventCategory_,*muon,jets,bTagAlgoWP_);
	    MCTruthL_evtcat6->fillMu(theMMEventCategory_,theEEEventCategory_,*muon,jets,bTagAlgoWP_);
	  }
	  MCTruthL_evtcat8->fillMu(theMMEventCategory_,theEEEventCategory_,*muon,jets,bTagAlgoWP_);
	  MCTruthL_evtcat10->fillMu(theMMEventCategory_,theEEEventCategory_,*muon,jets,bTagAlgoWP_);
	  MCTruthL_evtcat11->fillMu(theMMEventCategory_,theEEEventCategory_,*muon,jets,bTagAlgoWP_);
	}// end if decay from B
      }// end if muon is matched with GenMuon
    }// end loop over muons
    
    //Electrons
    for(edm::View<pat::Electron>::const_iterator electron=electrons->begin(); electron!=electrons->end(); ++electron){
      if( electron->genParticleById(11,1).isNonnull() ) {
	const reco::Candidate *GenElectron = &(*electron->genParticleById(11,1));
	while (GenElectron->mother()!=0 && abs(GenElectron->mother()->pdgId()) == 11) GenElectron = GenElectron->mother();
	if (ZbbUtils::isPreselEle(*electron,theZVertexFromProd,lCuts_)){
	  if( doAllThePlotting_){
	    MCTruthL_evtcat4->fillEle(theMMEventCategory_,theEEEventCategory_,*electron,jets,bTagAlgoWP_);
	    MCTruthL_evtcat5->fillEle(theMMEventCategory_,theEEEventCategory_,*electron,jets,bTagAlgoWP_);
	    MCTruthL_evtcat6->fillEle(theMMEventCategory_,theEEEventCategory_,*electron,jets,bTagAlgoWP_);
	  }
	  MCTruthL_evtcat8->fillEle(theMMEventCategory_,theEEEventCategory_,*electron,jets,bTagAlgoWP_);
	  MCTruthL_evtcat10->fillEle(theMMEventCategory_,theEEEventCategory_,*electron,jets,bTagAlgoWP_);
	  MCTruthL_evtcat11->fillEle(theMMEventCategory_,theEEEventCategory_,*electron,jets,bTagAlgoWP_);
	}// end if decay from B
      }// end if electron is matched with GenElectron
    }// end loop over electrons

    //Fill genParticles
    for( reco::GenParticleCollection::const_iterator genp=genParticlesCollection->begin(); genp!=genParticlesCollection->end(); genp++){
      const reco::Candidate *GenParticle = &(*genp);
      
      //MUONS
      if ( abs(GenParticle->pdgId()==13) && GenParticle->status()==1 ){
	while (GenParticle->mother()!=0 && abs(GenParticle->mother()->pdgId()) == 13) GenParticle = GenParticle->mother();
	MCTruthL_Mu_all->fillGenParticle(*GenParticle);
	if (JetMCTagUtils::decayFromBHadron(*GenParticle) || JetMCTagUtils::decayFromCHadron(*GenParticle)) {
	  MCTruthL_Mu_fromHF->fillGenParticle(*GenParticle);
	  if (JetMCTagUtils::decayFromBHadron(*GenParticle) && !(JetMCTagUtils::decayFromCHadron(*GenParticle))){
	    MCTruthL_Mu_fromB->fillGenParticle(*GenParticle);
	  }
	  if (!(JetMCTagUtils::decayFromBHadron(*GenParticle)) && JetMCTagUtils::decayFromCHadron(*GenParticle)){
	    MCTruthL_Mu_fromC->fillGenParticle(*GenParticle);
	  }
	  if (JetMCTagUtils::decayFromBHadron(*GenParticle) && JetMCTagUtils::decayFromCHadron(*GenParticle)){
	    MCTruthL_Mu_fromBC->fillGenParticle(*GenParticle);
	  }
	}
      }//end if muon

      //ELECTRONS
      if ( abs(GenParticle->pdgId()==11) && GenParticle->status()==1 ){
	while (GenParticle->mother()!=0 && abs(GenParticle->mother()->pdgId()) == 11) GenParticle = GenParticle->mother();
	MCTruthL_Ele_all->fillGenParticle(*GenParticle);
	if (JetMCTagUtils::decayFromBHadron(*GenParticle) || JetMCTagUtils::decayFromCHadron(*GenParticle)) {
	  MCTruthL_Ele_fromHF->fillGenParticle(*GenParticle);
	  if (JetMCTagUtils::decayFromBHadron(*GenParticle) && !(JetMCTagUtils::decayFromCHadron(*GenParticle))){
	    MCTruthL_Ele_fromB->fillGenParticle(*GenParticle);
	  }
	  if (!(JetMCTagUtils::decayFromBHadron(*GenParticle)) && JetMCTagUtils::decayFromCHadron(*GenParticle)){
	    MCTruthL_Ele_fromC->fillGenParticle(*GenParticle);
	  }
	  if (JetMCTagUtils::decayFromBHadron(*GenParticle) && JetMCTagUtils::decayFromCHadron(*GenParticle)){
	    MCTruthL_Ele_fromBC->fillGenParticle(*GenParticle);
	  }
	}
      }//end if electron
    }//end loop over genParticles
  }// end if MC 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//Zbbtd::EventContentAnalyzer::beginJob
///////////////////////////////////////////////////////////////////////////////////////////////////

void 
ZbbEventContentAnalyzer::beginJob()
{

  // directly taken from UserCode/Torino/ZZAnalysis/plugins/ZZAnalyzer.cc
  if ( ismc_ ) {
    
    //std::vector<float> dataPU(36);
    //std::vector<float> MCPU(36);
    
    // for EPS analysis from ZZ analysis
    //copy(data_EPS11, data_EPS11+25, dataPU.begin());
    //copy(PoissonOneXDist, PoissonOneXDist+25, MCPU.begin());

    // for 2/fb analysis from https://twiki.cern.ch/twiki/bin/view/CMS/VectorBosonPlusHeavyFlavor#Pile_up_reweighting
    //copy(data_2011_to_run173692_v2,data_2011_to_run173692_v2+36,dataPU.begin());
    //copy(MC_summer11_PUS4,MC_summer11_PUS4+36,MCPU.begin());

    std::vector<float> dataPU = ZbbUtils::getPU(true);
    std::vector<float> MCPU   = ZbbUtils::getPU(false);

    theLumiWeights_ = new edm::LumiReWeighting(MCPU, dataPU);
  }

  // category counters                 name                      cat n.
  wasRun=0;                         // All evts                     0       
  wasTrigOk=0;                      // HLT                          1
  wasLL=0;                          // ll                           2 
  wasLLTight=0;                     // ll tight                     3
  wasLLTightAndTrigMatch=0;         // l-HLT match                  4
  wasZLL=0;                         // Z(ll)                        5
  wasGoodJet=0;                     // Z(ll) + jet                  6
  wasJetVertexAssoc=0;              // Z(ll) + jet vertex assoc     7
  wasJetBTag=0;                     // Z(ll) + b                    8   
  wasJetBTagExcl=0;                 // Z(ll) + nb (exclusive)       9
  wasJetBTagPlusl=0;                // Z(ll) + b + l               10 
  wasJetBTagPlus2l=0;               // Z(ll) + b + ll              11

  wasJetBTagMCMatched=0;            // extra category for MC studies

  outfile_.open (outfilename_);
  // register to the TFileService
  edm::Service<TFileService> fs;
  TFileDirectory EventYields = fs->mkdir("EventYields");

  // book histograms:
  unsigned int bin1(1),bin2(1); 
  histContainer_["muonTriggerBits"]  = EventYields.make<TH1F>("muonTriggerBits","muon trigger path firing the event,",theAssociativeMapMu_.size(),-0.5,theAssociativeMapMu_.size()-0.5);
  for (std::map<std::vector<std::string>,std::vector<uint> >::iterator it=theAssociativeMapMu_.begin() ; it !=theAssociativeMapMu_.end(); it++, ++bin1){
    histContainer_["muonTriggerBits"]->GetXaxis()->SetBinLabel(bin1,((*it).first[0]).c_str()); 
  }

  histContainer_["electronTriggerBits"]  = EventYields.make<TH1F>("electronTriggerBits","electron trigger path firing the event,",theAssociativeMapEle_.size(),-0.5,theAssociativeMapEle_.size()-0.5);
  for (std::map<std::vector<std::string>,std::vector<uint> >::iterator it=theAssociativeMapEle_.begin() ; it !=theAssociativeMapEle_.end(); it++, ++bin2){
    histContainer_["electronTriggerBits"]->GetXaxis()->SetBinLabel(bin2,((*it).first[0]).c_str());  
  }

  histContainer_["eventCategoryMu"]  = EventYields.make<TH1F>("eventCategoryMu","Cut-flow Summary (number of events AFTER each cut - muon channel)",12,-0.5,11.5);
  histContainer_["eventCategoryEle"] = EventYields.make<TH1F>("eventCategoryEle","Cut-flow Summary (number of events AFTER each cut - electron channel)",12,-0.5,11.5);
  histContainer_["eventCategoryAll"] = EventYields.make<TH1F>("eventCategoryAll","Cut-flow Summary (number of events AFTER each cut - ee + #mu#mu channels)",12,-0.5,11.5);

  histContainer_["UnwEventCategoryMu"]  = EventYields.make<TH1F>("UnwEventCategoryMu","Cut-flow Summary (muon channel)",12,-0.5,11.5);
  histContainer_["UnwEventCategoryEle"] = EventYields.make<TH1F>("UnwEventCategoryEle","Cut-flow Summary (electron channel)",12,-0.5,11.5);
  histContainer_["UnwEventCategoryAll"] = EventYields.make<TH1F>("UnwEventCategoryAll","Cut-flow Summary (ee + #mu#mu channels)",12,-0.5,11.5);

  histContainer_["eventWeight"]      = EventYields.make<TH1F>("eventWeight","event weight; w_{event}",150,0.0,1.5);
  histContainer_["nVvertexReweight"] = EventYields.make<TH1F>("NVvertexReweight","Nb. of valid vertices (reweighted); N_{VTX} (reweighted)",26,-0.5,25.5);
  histContainer_["nVvertex"]         = EventYields.make<TH1F>("NVvertex","Nb. of valid vertices; N_{VTX}",26,-0.5,25.5);
  histContainer_["nPileUp"]          = EventYields.make<TH1F>("NPileup","Nb. of pile-up interactions; n_{PU} (MC only)",26,-0.5,25.5);
  histContainer_["nPileUpReweight"]  = EventYields.make<TH1F>("NPileupReweight","Nb. of pile-up intereactions (reweighted); n_{PU} (reweighted)",26, -0.5, 25.5);
  histContainer_["nBunchCrossing"]   = EventYields.make<TH1F>("NBunchCrossing","Bunch crossing assignment to PU; bx_{PU} (MC only)",5,-2.5,2.5);
  
  TString MuEvtCategoryBinLabels[12]  = {"All evts","HLT","ll","ll tight","l-HLT match","Z(#mu#mu)","Z(#mu#mu)+j","Z(#mu#mu)+j_{vtx}","Z(#mu#mu)b","Z(#mu#mu)b Excl","Z(#mu#mu)b+l","Z(#mu#mu)b+2l"};
  TString EleEvtCategoryBinLabels[12] = {"All evts","HLT","ll","ll tight","l-HLT match","Z(ee)","Z(ee)+j","Z(ee)+j_{vtx}","Z(ee)b","Z(ee)b Excl","Z(ee)b+l","Z(ee)b+2l"};
  TString AllEvtCategoryBinLabels[12] = {"All evts","HLT","ll","ll tight","l-HLT match","Z(ll)","Z(ll)+j","Z(ll)+j_{vtx}","Z(ll)b","Z(ll)b Excl","Z(ll)b+l","Z(ll)b+2l"};
  
  for(UInt_t bin=1; bin<=12; bin++){
    histContainer_["eventCategoryMu"]->GetXaxis()->SetBinLabel(bin,MuEvtCategoryBinLabels[bin-1]); 
    histContainer_["eventCategoryEle"]->GetXaxis()->SetBinLabel(bin,EleEvtCategoryBinLabels[bin-1]);    
    histContainer_["eventCategoryAll"]->GetXaxis()->SetBinLabel(bin,AllEvtCategoryBinLabels[bin-1]);    
  }

  // trend plot events vs run
  
  TFileDirectory YieldsVsRun = fs->mkdir("YieldsVsRun");

  Float_t runNumberFirst_= 160410.-0.5;
  Float_t runNumberLast_ = 174000.+0.5;
  Int_t    nRunBins=runNumberLast_-runNumberFirst_;
  
  char evtVsRunMuTitle_[128],evtVsRunEleTitle_[128], evtVsRunAllTitle_[128];
  char evtVsRunMuName_[12],evtVsRunEleName_[12],evtVsRunAllName_[12];

  char xsecVsRunMuTitle_[128],xsecVsRunEleTitle_[128], xsecVsRunAllTitle_[128];
  char xsecVsRunMuName_[12],xsecVsRunEleName_[12],xsecVsRunAllName_[12];
  
  for(UInt_t nCut=0; nCut<12; nCut++){

    // number of events plots
    
    // muon channel
    sprintf(evtVsRunMuName_,"evtsVsRunMu_after%i_cut",nCut);
    sprintf(evtVsRunMuTitle_,"n.(#mu#mu) candidates after %s cut;run number; n.(#mu#mu)",MuEvtCategoryBinLabels[nCut].Data());
    histContainer_[evtVsRunMuName_]  = YieldsVsRun.make<TH1F>(evtVsRunMuName_,evtVsRunMuTitle_,nRunBins,runNumberFirst_,runNumberLast_);

    // electron channel
    sprintf(evtVsRunEleName_,"evtsVsRunEle_after%i_cut",nCut);
    sprintf(evtVsRunEleTitle_,"n.(#it{ee}) candidates after %s cut; run number; n.(#it{ee})",EleEvtCategoryBinLabels[nCut].Data());
    histContainer_[evtVsRunEleName_] = YieldsVsRun.make<TH1F>(evtVsRunEleName_,evtVsRunEleTitle_,nRunBins,runNumberFirst_,runNumberLast_);
    
    // combining channels 
    sprintf(evtVsRunAllName_,"evtsVsRunAll_after%i_cut",nCut);
    sprintf(evtVsRunAllTitle_,"n.(ll) candidates after %s cut; run number; n.(ll)",AllEvtCategoryBinLabels[nCut].Data());
    histContainer_[evtVsRunAllName_] = YieldsVsRun.make<TH1F>(evtVsRunAllName_,evtVsRunAllTitle_,nRunBins,runNumberFirst_,runNumberLast_);

    // cross-sections plots

    // muon channel
    sprintf(xsecVsRunMuName_,"xsecVsRunMu_after%i_cut",nCut);
    sprintf(xsecVsRunMuTitle_,"n.(#mu#mu)/L candidates after %s cut;run number; #sigma B(#mu#mu) [pb]",MuEvtCategoryBinLabels[nCut].Data());
    histContainer_[xsecVsRunMuName_]  = YieldsVsRun.make<TH1F>(xsecVsRunMuName_,xsecVsRunMuTitle_,nRunBins,runNumberFirst_,runNumberLast_);

    // electron channel
    sprintf(xsecVsRunEleName_,"xsecVsRunEle_after%i_cut",nCut);
    sprintf(xsecVsRunEleTitle_,"n.(#it{ee})/L candidates after %s cut; run number; #sigma B(#it{ee}) [pb]",EleEvtCategoryBinLabels[nCut].Data());
    histContainer_[xsecVsRunEleName_] = YieldsVsRun.make<TH1F>(xsecVsRunEleName_,xsecVsRunEleTitle_,nRunBins,runNumberFirst_,runNumberLast_);
    
    // combining channels 
    sprintf(xsecVsRunAllName_,"xsecVsRunAll_after%i_cut",nCut);
    sprintf(xsecVsRunAllTitle_,"n.(ll)/L candidates after %s cut; run number;  #sigma B(ll) [pb]",AllEvtCategoryBinLabels[nCut].Data());
    histContainer_[xsecVsRunAllName_] = YieldsVsRun.make<TH1F>(xsecVsRunAllName_,xsecVsRunAllTitle_,nRunBins,runNumberFirst_,runNumberLast_);

  }

  // store the SumW2 struct 
  histContainer_["eventCategoryMu"]->Sumw2(); 
  histContainer_["eventCategoryEle"]->Sumw2();   
  histContainer_["eventCategoryAll"]->Sumw2();   
  
  for(UInt_t nCut=0; nCut<12; nCut++){
    histContainer_[(evtVsRunMuName[nCut]).Data()]->Sumw2(); 
    histContainer_[(evtVsRunEleName[nCut]).Data()]->Sumw2();   
    histContainer_[(evtVsRunAllName[nCut]).Data()]->Sumw2();  

    histContainer_[(xsecVsRunMuName[nCut]).Data()]->Sumw2(); 
    histContainer_[(xsecVsRunEleName[nCut]).Data()]->Sumw2();   
    histContainer_[(xsecVsRunAllName[nCut]).Data()]->Sumw2(); 
  }

  histContainer_["eventWeight"]->Sumw2(); 
  histContainer_["nVvertexReweight"]->Sumw2(); 
  histContainer_["nVvertex"]->Sumw2(); 
  histContainer_["nPileUp"]->Sumw2(); 
  histContainer_["nPileUpReweight"]->Sumw2(); 
  histContainer_["nBunchCrossing"]->Sumw2(); 

  // histogram containers (basic objects)

  if( doAllThePlotting_){
    ZbbBasicComponents_evtcat4 = new ZbbBasicComponents("LLTightAndTrigMatch",lCuts_);
    ZbbBasicComponents_evtcat4->book();
    
    ZbbBasicComponents_evtcat5 = new ZbbBasicComponents("ZLL",lCuts_);
    ZbbBasicComponents_evtcat5->book();
    
    ZbbBasicComponents_evtcat6 = new ZbbBasicComponents("GoodJet",lCuts_);
    ZbbBasicComponents_evtcat6->book();
  }
  
  ZbbBasicComponents_evtcat8 = new ZbbBasicComponents("JetBTag",lCuts_);
  ZbbBasicComponents_evtcat8->book();

  // histogram containers (intermediate steps)

  if( doAllThePlotting_){
    ZbbCollections_evtcat4 = new ZbbCollections("LLTightAndTrigMatch",lCuts_);
    ZbbCollections_evtcat4->book();
    
    ZbbCollections_evtcat5 = new ZbbCollections("ZLL",lCuts_);
    ZbbCollections_evtcat5->book();

    ZbbCollections_evtcat6 = new ZbbCollections("GoodJet",lCuts_);
    ZbbCollections_evtcat6->book();

    ZbbCollections_evtcat7 = new ZbbCollections("GoodJetAndVertexAssoc",lCuts_);
    ZbbCollections_evtcat7->book();
  }  

  // histogram containers (for Z+b candidates)
  ZbbCandidate_evtcat8  = new ZbbCandidate("JetBTag",lCuts_);
  ZbbCandidate_evtcat8->book();

  ZbbCandidate_evtcat9  = new ZbbCandidate("JetBTagExcl",lCuts_);
  ZbbCandidate_evtcat9->book();

  ZbbCandidate_evtcat10 = new ZbbCandidate("JetBTagPlusl",lCuts_);
  ZbbCandidate_evtcat10->book();

  ZbbCandidate_evtcat11 = new ZbbCandidate("JetBTagPlus2l",lCuts_);
  ZbbCandidate_evtcat11->book();

  //n-tuple for Z+b candidates
  ZbbCandidateNtuple_evtcat8 = new ZbbCandidateNtuple("JetBTag",lCuts_); 
  ZbbCandidateNtuple_evtcat8->book();

  // histogram containers (for Z+b matched candidates)
  ZbbMatchedCandidate_evtcat12 = new ZbbMatchedCandidate("JetBTagMCMatched",lCuts_);
  ZbbMatchedCandidate_evtcat12->book();
  
  // histogram containers (all leptons)
  if( doAllThePlotting_){
    ZbbAllLeptons_evtcat4  = new ZbbAllLeptons("LLTightAndTrigMatch",lCuts_);
    ZbbAllLeptons_evtcat4->book();
    
    ZbbAllLeptons_evtcat5  = new ZbbAllLeptons("ZLL",lCuts_);
    ZbbAllLeptons_evtcat5->book();
    
    ZbbAllLeptons_evtcat6  = new ZbbAllLeptons("GoodJet",lCuts_);
    ZbbAllLeptons_evtcat6->book();
  }

  ZbbAllLeptons_evtcat8  = new ZbbAllLeptons("JetBTag",lCuts_);
  ZbbAllLeptons_evtcat8->book();
  
  ZbbAllLeptons_evtcat9  = new ZbbAllLeptons("JetBTagExcl",lCuts_);
  ZbbAllLeptons_evtcat9->book();

  ZbbAllLeptons_evtcat10 = new ZbbAllLeptons("JetBTagPlusl",lCuts_);
  ZbbAllLeptons_evtcat10->book();

  ZbbAllLeptons_evtcat11 = new ZbbAllLeptons("JetBTagPlus2l",lCuts_);
  ZbbAllLeptons_evtcat11->book();

  // histogram containers (extra leptons)
  if( doAllThePlotting_){
    ZbbExtraLeptons_evtcat4  = new ZbbExtraLeptons("LLTightAndTrigMatch",lCuts_);
    ZbbExtraLeptons_evtcat4->book();
    
    ZbbExtraLeptons_evtcat5  = new ZbbExtraLeptons("ZLL",lCuts_);
    ZbbExtraLeptons_evtcat5->book();
    
    ZbbExtraLeptons_evtcat6  = new ZbbExtraLeptons("GoodJet",lCuts_);
    ZbbExtraLeptons_evtcat6->book(); 
  }

  ZbbExtraLeptons_evtcat8  = new ZbbExtraLeptons("JetBTag",lCuts_);
  ZbbExtraLeptons_evtcat8->book();
  
  ZbbExtraLeptons_evtcat9  = new ZbbExtraLeptons("JetBTagExcl",lCuts_);
  ZbbExtraLeptons_evtcat9->book();
  
  ZbbExtraLeptons_evtcat10 = new ZbbExtraLeptons("JetBTagPlusl",lCuts_);
  ZbbExtraLeptons_evtcat10->book();
  
  ZbbExtraLeptons_evtcat11 = new ZbbExtraLeptons("JetBTagPlus2l",lCuts_);
  ZbbExtraLeptons_evtcat11->book();
  
  // histogram containers (for MC truth)
  if( doAllThePlotting_){
    ZbbMCinfo_evtcat0 = new ZbbMCinfo("AllEvts",lCuts_);
    ZbbMCinfo_evtcat0->book();
    
    ZbbMCinfo_evtcat5 = new ZbbMCinfo("ZLL",lCuts_);
    ZbbMCinfo_evtcat5->book();
    
    ZbbMCinfo_evtcat6 = new ZbbMCinfo("GoodJet",lCuts_);
    ZbbMCinfo_evtcat6->book(); 
  }
  
  ZbbMCinfo_evtcat8 = new ZbbMCinfo("JetBTag",lCuts_);
  ZbbMCinfo_evtcat8->book();

  //MCTruth Leptons from HF
  if( doAllThePlotting_){
    MCTruthL_evtcat4 = new MCTruthLFromHF("LLTightAndTrigMatch",lCuts_);
    MCTruthL_evtcat4->book();
    
    MCTruthL_evtcat5 = new MCTruthLFromHF("ZLL",lCuts_);
    MCTruthL_evtcat5->book();
    
    MCTruthL_evtcat6 = new MCTruthLFromHF("GoodJet",lCuts_);
    MCTruthL_evtcat6->book();
  }

  MCTruthL_evtcat8 = new MCTruthLFromHF("JetBTag",lCuts_);
  MCTruthL_evtcat8->book();
  
  MCTruthL_evtcat10 = new MCTruthLFromHF("JetBTagPlusl",lCuts_);
  MCTruthL_evtcat10->book();
  
  MCTruthL_evtcat11 = new MCTruthLFromHF("JetBTagPlus2l",lCuts_);
  MCTruthL_evtcat11->book();
  
  //MCTruth Leptons from HF

  MCTruthL_Mu_all    = new MCTruthLFromHF("MCTruthL_Mu_all",lCuts_);
  MCTruthL_Mu_all->book();
  
  MCTruthL_Mu_fromHF = new MCTruthLFromHF("MCTruthL_Mu_fromHF",lCuts_);
  MCTruthL_Mu_fromHF->book();
  
  MCTruthL_Mu_fromB  = new MCTruthLFromHF("MCTruthL_Mu_fromB",lCuts_);
  MCTruthL_Mu_fromB->book();
  
  MCTruthL_Mu_fromC  = new MCTruthLFromHF("MCTruthL_Mu_fromC",lCuts_);
  MCTruthL_Mu_fromC->book();
  
  MCTruthL_Mu_fromBC = new MCTruthLFromHF("MCTruthL_Mu_fromBC",lCuts_);
  MCTruthL_Mu_fromBC->book();
  
  MCTruthL_Ele_all    = new MCTruthLFromHF("MCTruthL_Ele_all",lCuts_);
  MCTruthL_Ele_all->book();
  
  MCTruthL_Ele_fromHF = new MCTruthLFromHF("MCTruthL_Ele_fromHF",lCuts_);
  MCTruthL_Ele_fromHF->book();

  MCTruthL_Ele_fromB  = new MCTruthLFromHF("MCTruthL_Ele_fromB",lCuts_);
  MCTruthL_Ele_fromB->book();

  MCTruthL_Ele_fromC  = new MCTruthLFromHF("MCTruthL_Ele_fromC",lCuts_);
  MCTruthL_Ele_fromC->book();

  MCTruthL_Ele_fromBC = new MCTruthLFromHF("MCTruthL_Ele_fromBC",lCuts_);
  MCTruthL_Ele_fromBC->book();

  // histograms for ABCD
  if( doAllThePlotting_){
    ZbbABCDMatrix_evtcat0 = new ZbbABCDMatrix("AllEvts");
    ZbbABCDMatrix_evtcat0->book();
    
    ZbbABCDMatrix_evtcat5 = new ZbbABCDMatrix("ZLL");
    ZbbABCDMatrix_evtcat5->book();
    
    ZbbABCDMatrix_evtcat6 = new ZbbABCDMatrix("GoodJet");
    ZbbABCDMatrix_evtcat6->book();
    
    ZbbABCDMatrix_evtcat7 = new ZbbABCDMatrix("GoodJetAndVertexAssoc");
    ZbbABCDMatrix_evtcat7->book();
  }

  ZbbABCDMatrix_evtcat8 = new ZbbABCDMatrix("JetBTag");
  ZbbABCDMatrix_evtcat8->book();

  // vertex histograms 
  if( doAllThePlotting_){
    VtxHistos_evcat5 = new VtxHistos("ZLL",lCuts_);
    VtxHistos_evcat5->book();
    
    VtxHistos_evcat6 = new VtxHistos("GoodJet",lCuts_);
    VtxHistos_evcat6->book();
    
    VtxHistos_evcat7 = new VtxHistos("GoodJetAndVertexAssoc",lCuts_);
    VtxHistos_evcat7->book();
  }
  
  VtxHistos_evcat8 = new VtxHistos("JetBTag",lCuts_);
  VtxHistos_evcat8->book();
  
  //++++++++++++++++++++++++++++++++++++
  // Lumi stuff
  //++++++++++++++++++++++++++++++++++++
  
  totLS_ = 0;  
  totLumiRecorded_=0; 
  totLumiDelivered_=0;

  invalidLS_totLS_ = 0;
  invalidLS_totLumiRecorded_=0;
  invalidLS_totLumiDelivered_=0;

  reclumibyrun_.clear();

}

///////////////////////////////////////////////////////////////////////////////////////////////////
//ZbbEventContentAnalyzer::endJob
///////////////////////////////////////////////////////////////////////////////////////////////////

void 
ZbbEventContentAnalyzer::endJob() 
{

  using namespace std;

  // adding up combined histograms

  histContainer_["eventCategoryAll"]->Add(histContainer_["eventCategoryEle"]);
  histContainer_["eventCategoryAll"]->Add(histContainer_["eventCategoryMu"]);

  histContainer_["UnwEventCategoryAll"]->Add(histContainer_["UnwEventCategoryEle"]);
  histContainer_["UnwEventCategoryAll"]->Add(histContainer_["UnwEventCategoryMu"]);
  
  //  for(Int_t nCut=0;nCut<12;nCut++){
  //     histContainer_[(evtVsRunAllName[nCut]).Data()]->Add(histContainer_[(evtVsRunMuName[nCut]).Data()]);
  //     histContainer_[(evtVsRunAllName[nCut]).Data()]->Add(histContainer_[(evtVsRunEleName[nCut]).Data()]); 
  //   }

  // special treatment since trigger is not emulated in MC
  if(ismc_ ) wasTrigOk= wasTrigOk/2.;

  if(debug_){

    cout<<"***********************************"<<endl;
    cout<<"* ZbbEventContentAnalyzer::endJob()"<<endl;
    cout<<"* acceptance cuts:                 "<<endl;  
    lCuts_.print();
    cout<<"***********************************"<<endl;
    cout<<"* "<<bTagAlgoWP_<<" Selection"<<endl;
    cout<<"***********************************"<<endl;
    cout<<"* Was run on        : "<<wasRun<<                " events"<<endl;
    cout<<"* Was trigger ok on : "<<wasTrigOk<<             " events"<<endl;
    cout<<"* Was di-lepton     : "<<wasLL<<                 " events"<<endl;
    cout<<"* Was tight di-lep  : "<<wasLLTight<<            " events"<<endl;
    cout<<"* Was di-lep match  : "<<wasLLTightAndTrigMatch<<" events"<<endl;
    cout<<"* Was a Z evt       : "<<wasZLL<<                " events"<<endl;
    cout<<"* Was good jet      : "<<wasGoodJet<<            " events"<<endl;
    cout<<"* Was jet-vtx assoc : "<<wasJetVertexAssoc<<     " events"<<endl;
    cout<<"* Was b-tag         : "<<wasJetBTag<<            " events"<<endl;
    cout<<"* Was b-tag Excl FS : "<<wasJetBTagExcl<<        " events"<<endl;
    cout<<"* Was b-tag + l     : "<<wasJetBTagPlusl<<       " events"<<endl;
    cout<<"* Was b-tag + 2l    : "<<wasJetBTagPlus2l<<      " events"<<endl;
    cout<<"***********************************"<<endl;
    
    if(ismc_){
      cout<<"\n";
      cout<<"* Was b-tag jetMC-Matched : "<<wasJetBTagMCMatched<<" events"<<endl;
      cout<<"\n";
    }
    
  }

  double all = wasRun;

  // absolute efficiencies

  pair<Double_t,Double_t> etrig      = ZbbUtils::effCalc(wasTrigOk,all);
  pair<Double_t,Double_t> epresel    = ZbbUtils::effCalc(wasLL,all);
  pair<Double_t,Double_t> eid        = ZbbUtils::effCalc(wasLLTight,all);
  pair<Double_t,Double_t> etrigmatch = ZbbUtils::effCalc(wasLLTightAndTrigMatch,all);
  pair<Double_t,Double_t> ecand      = ZbbUtils::effCalc(wasZLL,all);
  pair<Double_t,Double_t> ejet       = ZbbUtils::effCalc(wasGoodJet,all);
  pair<Double_t,Double_t> ebtag      = ZbbUtils::effCalc(wasJetBTag,all);
  pair<Double_t,Double_t> ebtagEx    = ZbbUtils::effCalc(wasJetBTagExcl,all);
  pair<Double_t,Double_t> ebtag1l    = ZbbUtils::effCalc(wasJetBTagPlusl,all);
  pair<Double_t,Double_t> ebtag2l    = ZbbUtils::effCalc(wasJetBTagPlus2l,all);

  // relative efficiencies

  pair<Double_t,Double_t> eTrigToPresel    = ZbbUtils::effCalc(wasLL,wasTrigOk);
  pair<Double_t,Double_t> ePreselToId      = ZbbUtils::effCalc(wasLLTight,wasLL);
  pair<Double_t,Double_t> eIdToTrigMatch   = ZbbUtils::effCalc(wasLLTightAndTrigMatch,wasLLTight);
  pair<Double_t,Double_t> eTrigMatchToCand = ZbbUtils::effCalc(wasZLL,wasLLTightAndTrigMatch);
  pair<Double_t,Double_t> eCandToJet       = ZbbUtils::effCalc(wasGoodJet,wasZLL);
  pair<Double_t,Double_t> eJetToBtag       = ZbbUtils::effCalc(wasJetBTag,wasGoodJet);
  pair<Double_t,Double_t> eJetToBtagEx     = ZbbUtils::effCalc(wasJetBTagExcl,wasGoodJet);
  pair<Double_t,Double_t> eBtagTo1l        = ZbbUtils::effCalc(wasJetBTagPlusl,wasJetBTag);
  pair<Double_t,Double_t> eBtagTo2l        = ZbbUtils::effCalc(wasJetBTagPlus2l,wasJetBTag);

  using namespace edm;
 
  // starts Report

  cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
  cout << "* Acceptance cuts selected with        "<<bTagAlgoWP_<<endl;  
  lCuts_.print();
  cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
  cout << "\n>>>>>> Z+ b-JETS  SELECTION SUMMARY BEGIN >>>>>>>>>>>>>>>"<<endl;  
  cout << "Number of events analyzed                          : " << wasRun           << " [events]"<<endl;
  cout << "Number of events triggered                         : " << wasTrigOk        << " [events]"<<endl;
  cout << "Number of events with di-leptons                   : " << wasLL            << " [events]"<<endl;
  cout << "Number of events after vector boson ID Cuts        : " << wasLLTight       << " [events]"<<endl;
  cout << "Number of events after HLT matching                : " << wasLLTightAndTrigMatch << " [events]"<<endl;
  cout << "Number of events after candidate inv. mass cut     : " << wasZLL           << " [events]"<<endl;
  cout << "Number of events after request of good jet         : " << wasGoodJet       << " [events]"<<endl;
  cout << "Number of events after request of jet vtx assoc    : " << wasJetVertexAssoc<< " [events]"<<endl;
  cout << "Number of events after b-tagging of jet            : " << wasJetBTag       << " [events]"<<endl;
  cout << "Number of events after request of exclusive F.S.   : " << wasJetBTagExcl   << " [events]"<<endl;
  cout << "Number of events after request of 1 extra lepton   : " << wasJetBTagPlusl  << " [events]"<<endl;
  cout << "Number of events after request of 2 extra leptons  : " << wasJetBTagPlus2l << " [events]"<<endl;

  if(ismc_){
    cout<<"\n";
    cout<<"Number of events after request of b-tagMC-Matched  : " <<wasJetBTagMCMatched<<" [events]"<<endl;
    cout<<"\n";
  }
  
  cout << ">>>>>>>>>>>>>> Efficiencies CUMULATIVE >>>>>>>>>>>>>>>>>>>>"<<endl;
  cout << "Trigger Efficiency:                      " << "(" << setw(4) << etrig.first*100.      <<" +/- "<< setw(2) <<etrig.second*100 << ")%"<<endl;
  cout << "Di-Lepton Candidate Building:            " << "(" << setw(4) << epresel.first*100.    <<" +/- "<< setw(2) <<epresel.second*100 << ")%"<<endl;
  cout << "Muon/electron ID Efficiency:             " << "(" << setw(4) << eid.first*100.        <<" +/- "<< setw(2) <<eid.second*100 << ")%"<<endl;
  cout << "HLT matching Efficiency:                 " << "(" << setw(4) << etrigmatch.first*100. <<" +/- "<< setw(2) <<etrigmatch.second*100 << ")%"<<endl;
  cout << "Invariant mass cut Efficiency:           " << "(" << setw(4) << ecand.first*100.      <<" +/- "<< setw(2) <<ecand.second*100 << ")%"<<endl;
  cout << "Jet selection Efficiency:                " << "(" << setw(4) << ejet.first*100.       <<" +/- "<< setw(2) <<ejet.second*100 << ")%"<<endl;
  cout << "b-Tagging Efficiency:                    " << "(" << setw(4) << ebtag.first*100.      <<" +/- "<< setw(2) <<ebtag.second*100 << ")%"<<endl;
  cout << "b-Tagging Efficiency (exclusive F.S.):   " << "(" << setw(4) << ebtagEx.first*100.    <<" +/- "<< setw(2) <<ebtagEx.second*100 << ")%"<<endl;
  cout << "Soft lepton Efficiency (after b-tag):    " << "(" << setw(4) << ebtag1l.first*100.    <<" +/- "<< setw(2) <<ebtag1l.second*100 << ")%"<<endl;
  cout << "2nd soft lepton Efficiency (after b-tag):" << "(" << setw(4) << ebtag2l.first*100.    <<" +/- "<< setw(2) <<ebtag2l.second*100 << ")%"<<endl;
  cout << ">>>>>>>>>>>>>>> Efficiencies RELATIVE >>>>>>>>>>>>>>>>>>>>>"<<endl;
  cout << "Di-Lepton from Trigger:                  " << "(" << setw(4) << eTrigToPresel.first*100.    <<" +/- "<<setw(2) <<eTrigToPresel.second*100. <<")%"<<endl;
  cout << "Lepton ID from Di-Lepton:                " << "(" << setw(4) << ePreselToId.first*100.      <<" +/- "<< setw(2) <<ePreselToId.second*100. << ")%"<<endl;
  cout << "HLT matching from Lepton ID              " << "(" << setw(4) << eIdToTrigMatch.first*100.   <<" +/- "<< setw(2) <<eIdToTrigMatch.second*100. << ")%"<<endl;
  cout << "Invariant mass from HLT matching:        " << "(" << setw(4) << eTrigMatchToCand.first*100. <<" +/- "<< setw(2) <<eTrigMatchToCand.second*100.<< ")%"<<endl;
  cout << "Jet selection from invariant mass:       " << "(" << setw(4) << eCandToJet.first*100.       <<" +/- "<< setw(2) <<eCandToJet.second*100. << ")%"<<endl;
  cout << "b-Tagging from jet selection:            " << "(" << setw(4) << eJetToBtag.first*100.       <<" +/- "<< setw(2) <<eJetToBtag.second*100. << ")%"<<endl;
  cout << "b-Tagging (excl F.S.) from jet selection:" << "(" << setw(4) << eJetToBtagEx.first*100.     <<" +/- "<< setw(2) <<eJetToBtagEx.second*100. << ")%"<<endl;
  cout << "Soft lepton from b-tag:                  " << "(" << setw(4) << eBtagTo1l.first*100.        <<" +/- "<< setw(2) <<eBtagTo1l.second*100. << ")%"<<endl;
  cout << "2nd soft lepton froma b-tag:             " << "(" << setw(4) << eBtagTo2l.first*100.        <<" +/- "<< setw(2) <<eBtagTo2l.second*100. << ")%"<<endl;
  cout << ">>>>>> Z+ b-JETS SELECTION SUMMARY END   >>>>>>>>>>>>>>>\n"<<endl;
  
  outfile_.close();

  //++++++++++++++++++++++++++++++++++++
  // Lumi stuff for data only
  //++++++++++++++++++++++++++++++++++++

  if(!ismc_){
    
    double ubarnsTopbarns = 1000000;
    double corrFactor     = 10;
    
    double CF = 1/(ubarnsTopbarns*corrFactor);
    
    map<int,pair<int,double> >::iterator it1;
    
    cout << ">>>>>> Begins Lumisummary  >>>>>>>>>>>>>>>"<<endl;
    for ( it1=reclumibyrun_.begin() ; it1 !=reclumibyrun_.end(); it1++ ){
      
      if(debug_){
	cout << (*it1).first << " => " << (*it1).second.first << " / "  
	     << setiosflags(ios::fixed) << setprecision(2) 
	     << (*it1).second.second * CF <<"/pb" << endl;
      }
      
      for(Int_t nCut=0;nCut<12;nCut++){
	Int_t    theBin        = histContainer_[(evtVsRunMuName[nCut]).Data()]->FindBin((*it1).first);
	Double_t theBinContentMu  = histContainer_[(evtVsRunMuName[nCut]).Data()]->GetBinContent(theBin);
	Double_t theBinContentEle = histContainer_[(evtVsRunEleName[nCut]).Data()]->GetBinContent(theBin);
	Double_t xsecMu(0.), xsecEle(0.);
	
	// there should be actaul lumi
	if((*it1).second.second>0.){ 
	  xsecMu  = theBinContentMu/((*it1).second.second * CF ); 
	  xsecEle = theBinContentEle/((*it1).second.second * CF ); 
	  histContainer_[(xsecVsRunMuName[nCut]).Data()]->SetBinContent(theBin,xsecMu);
	  histContainer_[(xsecVsRunEleName[nCut]).Data()]->SetBinContent(theBin,xsecEle);
	  histContainer_[(xsecVsRunAllName[nCut]).Data()]->SetBinContent(theBin,xsecMu+xsecEle);
	}
      }
    }
    
    // print out info for valid lumi summaries
    cout << "Total number of LS:   "<< totLS_ << endl; 
    cout << "Total delivered lumi: "<< setiosflags(ios::fixed) << setprecision(2) 
	 << totLumiDelivered_*CF << "/pb "<<endl; 
    cout << "Total recorded lumi:  "<< setiosflags(ios::fixed) << setprecision(2) 
	 << totLumiRecorded_*CF <<  "/pb"<< endl; 
    
    
    // print out info for invalid lumi summaries
    cout << "Total number of INVALID LS:   "<< invalidLS_totLS_ << endl;
    cout << "Total INVALID delivered lumi: "<< setiosflags(ios::fixed) << setprecision(2)
	 << invalidLS_totLumiDelivered_*CF << "/pb" << endl;
    cout << "Total INVALID recorded lumi:  "<< setiosflags(ios::fixed) << setprecision(2)
	 << invalidLS_totLumiRecorded_*CF << "/pb" << endl;
    cout << ">>>>>> ends Lumisummary  >>>>>>>>>>>>>>>\n"<<endl;
    
  }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
//ZbbEventContentAnalyzer::setEventCategory
///////////////////////////////////////////////////////////////////////////////////////////////////

void 
ZbbEventContentAnalyzer::setEventCategory(const edm::Event& iEvent, const edm::EventSetup& iSetup, EventCategory& iEvtCategory, Bool_t isMuChannel, Bool_t AndOrSwitch){

  if(debug_) std::cout<<"ZbbEventContentAnalyzer::setEventCategory"<<std::endl;

  // get beamSpot
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  if ( beamSpotHandle.isValid() ) {
    beamSpot = *beamSpotHandle;
    if( debug_ )  std::cout<<"Beamspot x: "<< beamSpot.x0() <<" y: "<< beamSpot.y0() <<" z: "<< beamSpot.z0() <<std::endl;
  } else {
    if( debug_ )  std::cout << "No BeamSpot found!" << std::endl;
  }

  // get met collection  
  edm::Handle<edm::View<pat::MET> > mets;
  iEvent.getByLabel(metSrc_,mets);

  // get electron collection
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByLabel(elecSrc_,electrons);

  // get muon collection
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(muonSrc_,muons);

  // get jet collection
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jetSrc_,jets);

  // get the Z->mm collection
  edm::Handle<reco::CompositeCandidateCollection> zmmHandle;
  iEvent.getByLabel(theZmmSrc_, zmmHandle);
  const reco::CompositeCandidateCollection & zmm = *(zmmHandle.product());

  // get the Z->ee collection
  edm::Handle<reco::CompositeCandidateCollection> zeeHandle;
  iEvent.getByLabel(theZeeSrc_, zeeHandle);
  const reco::CompositeCandidateCollection & zee = *(zeeHandle.product());
  
  // get offline primary vertex collection
  edm::Handle<edm::View<reco::Vertex> > vertexCollection;
  iEvent.getByLabel(theVertexSrc_, vertexCollection);
  edm::View<reco::Vertex> vertices = *(vertexCollection.product());

  // lists of track refs necessary for PU/PV beta/betastar
  std::list<int> trackrefs_PV = VtxAssociatorsUtils::buildTrackRefs(vertices,false);
  std::list<int> trackrefs_PU = VtxAssociatorsUtils::buildTrackRefs(vertices,true);

  // get Z vertexes collection
  edm::Handle<edm::View<reco::Vertex> > ZvertexCollection;
  iEvent.getByLabel(theZVertexSrc_, ZvertexCollection);
  const reco::Vertex& theZVertexFromProd = ZvertexCollection->at(0);
 
  // get the TriggerEvent
  edm::Handle<pat::TriggerEvent> theTriggerEventHandle;
  iEvent.getByLabel(theTriggerEventSrc_,theTriggerEventHandle);
  pat::TriggerEvent theTriggerEvent = *(theTriggerEventHandle.product()); 

  TString  jetcategory = "loose";
  Double_t themetcut = 40;
  
  iEvtCategory.MC_=ismc_;

  // get gen particle candidates
  edm::Handle<reco::GenParticleCollection> genParticlesCollection;
  iEvent.getByLabel(genPSrc_, genParticlesCollection);
  
  // get genJets
  edm::Handle<reco::GenJetCollection> genJetsHandle;
  iEvent.getByLabel(genJetSrc_,genJetsHandle);
  const reco::GenJetCollection & genJets = *(genJetsHandle.product());

  // applying PU corrections
  // again directly taken from UserCode/Torino/ZZAnalysis/plugins/ZZAnalyzer.cc
  if(ismc_ && applypucorrection_){
    double thePUWeight_(1.);
    const edm::EventBase* iEventB = dynamic_cast<const edm::EventBase*>(&iEvent);
    if (theLumiWeights_) thePUWeight_ = theLumiWeights_->weight(*iEventB);
    iEvtCategory.weight_ *=thePUWeight_;
  }

  // applies the MC weight corrections
  // for Sherpa (recipe from https://hypernews.cern.ch/HyperNews/CMS/get/ewk-vplusjets/317) and aMC@NLO
  if (ismc_ && applymcweight_) {  
    edm::Handle<GenEventInfoProduct> genEventInfoHandle;
    iEvent.getByLabel("generator",genEventInfoHandle);
    double wMCtmp = genEventInfoHandle->weight();

    // std::cout << "W MC " << wMCtmp<< std::endl; 
    
    // reweighting rule for NLO
    if(isMCatNLO_) {
      wMCtmp > 0 ?  wMCtmp=1. : wMCtmp=-1.;
    }
    iEvtCategory.weight_ *=wMCtmp;
  }

  // muon selection
  if ( isMuChannel ) {

    // check on the trigger
    if(ismc_){
      iEvtCategory.Trigger_=true;
    } else {
      if( this->isTriggerOK(theTriggerEvent,iEvent.run(),theAssociativeMapMu_)){
	iEvtCategory.Trigger_=true;
      }
    }
    
    if (zmm.size()>0) {
      iEvtCategory.LL_=true;
      std::vector<reco::CompositeCandidate> sortedZmmCands = ZbbUtils::sortCandidatesByDifference(zmm);
      std::vector<Bool_t> isTightZmmCand;
      std::pair<reco::CompositeCandidate, std::pair<Bool_t,Bool_t> > ZmmCandsAndQuals;

      for( std::vector<reco::CompositeCandidate>::const_iterator ZmmCand=sortedZmmCands.begin(); ZmmCand != sortedZmmCands.end(); ++ZmmCand){
	isTightZmmCand.push_back(ZbbUtils::isTightZCandidate(*ZmmCand,beamSpot,isMuChannel,lCuts_));
	ZmmCandsAndQuals.first  = *ZmmCand;
	ZmmCandsAndQuals.second.first = ZbbUtils::isTightZCandidate(*ZmmCand,beamSpot,isMuChannel,lCuts_);
	ZmmCandsAndQuals.second.second = this->isTriggerMatchOK(*ZmmCand,iEvent.run(),theAssociativeMapMu_,isMuChannel,AndOrSwitch);

	if( ZbbUtils::isTightZCandidate(*ZmmCand,beamSpot,isMuChannel,lCuts_)) {
	  iEvtCategory.bestZcandidate_= *ZmmCand;
	  iEvtCategory.LLTight_=true;

	  if(ismc_){
	    iEvtCategory.LLTightAndTriggerMatch_=true;
	  } else {
	    if(this->isTriggerMatchOK(*ZmmCand,iEvent.run(),theAssociativeMapMu_,isMuChannel,AndOrSwitch)){
	      iEvtCategory.LLTightAndTriggerMatch_=true;
	    }
	  }	  
	  break;
	}
      }
      
      // weight for lepton efficiency ------------------------------------------------------
      if(ismc_ && applyLeptonEff_ && iEvtCategory.LLTight_){
	
	double ptmu1= (iEvtCategory.bestZcandidate_.daughter(0))->pt();
	double etamu1=(iEvtCategory.bestZcandidate_.daughter(0))->eta();
	double ptmu2= (iEvtCategory.bestZcandidate_.daughter(1))->pt();
	double etamu2=(iEvtCategory.bestZcandidate_.daughter(1))->eta();
	
	// muon offline scale factors
	double wmuOFF1 =ZbbUtils::getMuonOfflineScaleFactor(ptmu1,etamu1,false);
	double wmuOFF2 =ZbbUtils::getMuonOfflineScaleFactor(ptmu2,etamu2,false);

	// HLT DoubleMu7 SF (triggered 2011A1 period)	
	double wmuTRGMu7_1 =ZbbUtils::getMuonTrgScaleFactor_H(ptmu1,etamu1,"2011A1",false);
	double wmuTRGMu7_2 =ZbbUtils::getMuonTrgScaleFactor_L(ptmu2,etamu2,"2011A1",false);

	// HLT Mu13Mu8 SF   (triggered 2011A2 period)
  	double wmu1TRGMu138_high =ZbbUtils::getMuonTrgScaleFactor_H(ptmu1,etamu1,"2011A2",false);
	double wmu1TRGMu138_low  =ZbbUtils::getMuonTrgScaleFactor_L(ptmu1,etamu1,"2011A2",false);	
	double wmu2TRGMu138_high =ZbbUtils::getMuonTrgScaleFactor_H(ptmu2,etamu2,"2011A2",false);
	double wmu2TRGMu138_low  =ZbbUtils::getMuonTrgScaleFactor_L(ptmu2,etamu2,"2011A2",false); 

	// Weighted average: valid for ~ 2 fb^-1 statistics
	double triggerW=0.11*(wmuTRGMu7_1*wmuTRGMu7_2)+0.89*(wmu1TRGMu138_high*wmu2TRGMu138_low + wmu1TRGMu138_low*wmu2TRGMu138_high - wmu1TRGMu138_high*wmu2TRGMu138_high); 
	
	iEvtCategory.weight_*=(wmuOFF1*wmuOFF2*triggerW);	
      }
      // more physical selection but non in the EWK box
      //if( fabs(iEvtCategory.bestZcandidate_.mass()-91.1876)<30 ) iEvtCategory.ZLL_=true;
      if( iEvtCategory.bestZcandidate_.mass() > 60 && iEvtCategory.bestZcandidate_.mass() < 120 ) iEvtCategory.ZLL_=true;
    }
  } else {
  // electron selection
    
    // check on the trigger
    if(ismc_){
      iEvtCategory.Trigger_=true;
    } else {
      if( this->isTriggerOK(theTriggerEvent,iEvent.run(),theAssociativeMapEle_)){
	iEvtCategory.Trigger_=true;
      }
    }
    
    if (zee.size()>0) {
      iEvtCategory.LL_=true;
      std::vector<reco::CompositeCandidate> sortedZeeCands = ZbbUtils::sortCandidatesByDifference(zee);
      std::vector<Bool_t> isTightZeeCand;
      std::pair<reco::CompositeCandidate, std::pair<Bool_t,Bool_t> > ZeeCandsAndQuals;

      for( std::vector<reco::CompositeCandidate>::const_iterator ZeeCand=sortedZeeCands.begin(); ZeeCand != sortedZeeCands.end(); ++ZeeCand){
	isTightZeeCand.push_back(ZbbUtils::isTightZCandidate(*ZeeCand,beamSpot,isMuChannel,lCuts_));
	ZeeCandsAndQuals.first  = *ZeeCand;
	ZeeCandsAndQuals.second.first = ZbbUtils::isTightZCandidate(*ZeeCand,beamSpot,isMuChannel,lCuts_);
	ZeeCandsAndQuals.second.second = this->isTriggerMatchOK(*ZeeCand,iEvent.run(),theAssociativeMapEle_,isMuChannel,AndOrSwitch);

	if( ZbbUtils::isTightZCandidate(*ZeeCand,beamSpot,isMuChannel,lCuts_)) {
	  iEvtCategory.bestZcandidate_= *ZeeCand;
	  iEvtCategory.LLTight_=true;

	  if(ismc_){
	    iEvtCategory.LLTightAndTriggerMatch_=true;
	  } else {
	    if(this->isTriggerMatchOK(*ZeeCand,iEvent.run(),theAssociativeMapEle_,isMuChannel,AndOrSwitch)){
	      iEvtCategory.LLTightAndTriggerMatch_=true;
	    }
	  }
	  break;
	}	
      }
      
      // weight for lepton efficiency ------------------------------------------------------
      if(ismc_ && applyLeptonEff_ && iEvtCategory.LLTight_){
	double ptele1= (iEvtCategory.bestZcandidate_.daughter(0))->pt();
	double etaele1=(iEvtCategory.bestZcandidate_.daughter(0))->eta();
	double ptele2= (iEvtCategory.bestZcandidate_.daughter(1))->pt();
	double etaele2=(iEvtCategory.bestZcandidate_.daughter(1))->eta();
	
	// electron offline scale factors
	double weleOFF1 =ZbbUtils::getElectronOfflineScaleFactor(ptele1,etaele1,false);
	double weleOFF2 =ZbbUtils::getElectronOfflineScaleFactor(ptele2,etaele2,false);
	
	// HLT HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5  SF (triggered 2011A1 period)
	double wele1TRG178_v5_high =ZbbUtils::getElectronTrgScaleFactor_H(ptele1,etaele1,"2011A1",false);
	double wele1TRG178_v5_low  =ZbbUtils::getElectronTrgScaleFactor_L(ptele1,etaele1,"2011A1",false);	
	double wele2TRG178_v5_high =ZbbUtils::getElectronTrgScaleFactor_H(ptele2,etaele2,"2011A1",false);
	double wele2TRG178_v5_low  =ZbbUtils::getElectronTrgScaleFactor_L(ptele2,etaele2,"2011A1",false);  

	// HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6  SF (triggered 2011A2 period)
	double wele1TRG178_v6_high =ZbbUtils::getElectronTrgScaleFactor_H(ptele1,etaele1,"2011A2",false);
	double wele1TRG178_v6_low  =ZbbUtils::getElectronTrgScaleFactor_L(ptele1,etaele1,"2011A2",false);	
	double wele2TRG178_v6_high =ZbbUtils::getElectronTrgScaleFactor_H(ptele2,etaele2,"2011A2",false);
	double wele2TRG178_v6_low  =ZbbUtils::getElectronTrgScaleFactor_L(ptele2,etaele2,"2011A2",false);  
	
	double triggerW = 0.50*(wele1TRG178_v5_high*wele2TRG178_v5_low + wele1TRG178_v5_low*wele2TRG178_v5_high - wele1TRG178_v5_high*wele2TRG178_v5_high)+0.50*(wele1TRG178_v6_high*wele2TRG178_v6_low + wele1TRG178_v6_low*wele2TRG178_v6_high - wele1TRG178_v6_high*wele2TRG178_v6_high);

	iEvtCategory.weight_*=(weleOFF1*weleOFF2*triggerW);
	
      }
      // more physical but not in the EWK box
      //if( fabs(iEvtCategory.bestZcandidate_.mass()-91.1876)<30 ) iEvtCategory.ZLL_=true;
      if( iEvtCategory.bestZcandidate_.mass() > 60 && iEvtCategory.bestZcandidate_.mass() < 120 ) iEvtCategory.ZLL_=true;
    }
  }
  
  if(iEvtCategory.LLTight_==true){
    iEvtCategory.theZvertex_=VtxAssociatorsUtils::findPrimaryVertex(iEvtCategory.bestZcandidate_,vertices);
  } else {
    iEvtCategory.theZvertex_=vertices.at(0);
  }

  Int_t nJets = 0;
  Int_t nAssocJets=0;
  Int_t nBjets = 0;
  Int_t nBMcMatchedJets=0;
 
  // vector of b-tagged jets 
  std::vector<pat::Jet> theBtaggedJets;

  if(iEvtCategory.LLTight_==true){
    for(edm::View<pat::Jet>::const_iterator jetcand=jets->begin(); jetcand!=jets->end(); ++jetcand){
      if(ZbbUtils::isJetIdOk((*jetcand),jetcategory) && ZbbUtils::isGoodJet((*jetcand),iEvtCategory.bestZcandidate_,lCuts_)){
	nJets++;
	// if(VtxAssociatorsUtils::jetVertex(iEvtCategory.theZvertex_,(*jetcand),4,4,0.15)){
	if(VtxAssociatorsUtils::beta(*jetcand,trackrefs_PV)<lCuts_.betaCut_){
	  nAssocJets++;  
	}
	if(ZbbUtils::isBJet((*jetcand),bTagAlgoWP_)){
	  nBjets++;
	  theBtaggedJets.push_back(*jetcand);
	  if(ismc_){
	    if(ZbbUtils::isBJetMCMatched((*jetcand),genParticlesCollection,genJets)){
	      nBMcMatchedJets++; 
	    }
	  }// if is mc makes the matching  
	} // close if b-tag
      } // close if good jet
    } // close the for on jets
  } // if there is at least a Z candidate

  if (nJets>0)      iEvtCategory.GoodJet_=true;
  if (nAssocJets>0) iEvtCategory.VertexAssoc_=true;

  // compute the weight related to b tagging efficiency (SF_b)
  // based on http://cmslxr.fnal.gov/lxr/source/RecoBTag/PerformanceDB/test/october/tests/PATJetAnalyzer.cc
  
  if(nBjets>0){
    
    if(nBjets>=minBtags_) {
      iEvtCategory.JetBTagged_=true; 

      if(ismc_){
	//std::cout<<"Number of matched b-jets/event: "<<nBMcMatchedJets<<std::endl;
	if(nBMcMatchedJets>=minBtags_){ 
	  iEvtCategory.JetBTaggedMCMatched_=true;
	  //std::cout<<"MC Matched b-jet event: "<<std::endl;
	}
      } else {
	iEvtCategory.JetBTaggedMCMatched_=true; 
      }

      if(nBjets==minBtags_) {
	iEvtCategory.ExclFinalState_=true; 
      }
    }
 
    iEvtCategory.nBTags_=nBjets;
    
    //compute SF only for true b jets
    //ZbbUtils::FLAVOUR flavour = ZbbUtils::getPartonFlavour(theBtaggedJets[0]);
    
    if ( ismc_ && applybeffcalib_ && iEvtCategory.isJetBTagged() ) {
      //iEvtCategory.weight_= ZbbUtils::getbEffScaleFactor(flavour,theBtaggedJets[0],iSetup,beffcalibmethod_,bmistagcalibmethod_);
      iEvtCategory.weight_*=ZbbUtils::getbEffScaleFactorAR(theAssociativeMapEffbMC_,theAssociativeMapEffcMC_,theBtaggedJets,iSetup,beffcalibmethod_,bmistagcalibmethod_,false,minBtags_,theSFkFactor_,theLEffkFactor_);
      
      if(iEvtCategory.isExclusive()){
	iEvtCategory.weightCorrectionFromInclToExcl_= ZbbUtils::getbEffScaleFactorAR(theAssociativeMapEffbMC_,theAssociativeMapEffcMC_,theBtaggedJets,iSetup,beffcalibmethod_,bmistagcalibmethod_,true,minBtags_,theSFkFactor_,theLEffkFactor_) / ZbbUtils::getbEffScaleFactorAR(theAssociativeMapEffbMC_,theAssociativeMapEffcMC_,theBtaggedJets,iSetup,beffcalibmethod_,bmistagcalibmethod_,false,minBtags_,theSFkFactor_,theLEffkFactor_);
	//std::cout<<"nBjets:"<<nBjets<<" changing factor:"<<iEvtCategory.getBTagEffWeightCF()<<std::endl; 
      } // ends if exclusive final state
    } // if b-eff calibration is switched on    
  } // close if at least 1 btag

  // vector of leptons not coming from Z
  std::vector<pat::Muon> addMu;
  std::vector<pat::Electron> addEle;
  
  if(iEvtCategory.ZLL_==true){
    //Additional muons
    if(muons->size()!=0){
      for(edm::View<pat::Muon>::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon){
	
	if (ZbbUtils::isPreselMu(*muon,iEvtCategory.theZvertex_,lCuts_)){
	  if (iEvtCategory.bestZcandidate_.daughter(0)->pt()!=0){
	    Double_t deltaR1   = ROOT::Math::VectorUtil::DeltaR(iEvtCategory.bestZcandidate_.daughter(0)->momentum(),muon->momentum());
	    Double_t deltaR2   = ROOT::Math::VectorUtil::DeltaR(iEvtCategory.bestZcandidate_.daughter(1)->momentum(),muon->momentum());
	    
	    //Veto on muons matched by deltaR with candidate's daughters
	    if(deltaR1 > 0.05 && deltaR2 > 0.05){
	      addMu.push_back(*muon);
	    }//end veto
	  }
	}//muon cut
      }//end loop over muons
    }//end if muons->size() != 0
    
    //Additional electrons
    if(electrons->size()!=0){
      for(edm::View<pat::Electron>::const_iterator electron=electrons->begin(); electron!=electrons->end(); ++electron){
	if (ZbbUtils::isPreselEle(*electron,iEvtCategory.theZvertex_,lCuts_)){
	  if (iEvtCategory.bestZcandidate_.daughter(0)->pt()!=0){
	    Double_t deltaR1   = ROOT::Math::VectorUtil::DeltaR(iEvtCategory.bestZcandidate_.daughter(0)->momentum(),electron->momentum());
	    Double_t deltaR2   = ROOT::Math::VectorUtil::DeltaR(iEvtCategory.bestZcandidate_.daughter(1)->momentum(),electron->momentum());
	    
	    //Veto on electrons matched by deltaR with candidate's daughters
	    if(deltaR1 > 0.05 && deltaR2 > 0.05){
	      addEle.push_back(*electron);
	    }//end veto
	  }
	}//electron cut
      }//end loop over electrons
    }//end if electrons->size()!=0
  } // if there is at least a Z candidate  
  
  if ( addMu.size()>0 || addEle.size()>0 ) {
    iEvtCategory.ThreeLeptons_=true;
    if ( (addMu.size()+addEle.size())>1 ) iEvtCategory.FourLeptons_=true;
  }
  
  if (ZbbUtils::isGoodMet((*mets)[0],themetcut)) iEvtCategory.MET_=true;

  return;
 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//ZbbEventContentAnalyzer::isTriggerOK
///////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t 
ZbbEventContentAnalyzer::isTriggerOK(pat::TriggerEvent theTriggerEvent,UInt_t theRunNumber, std::map<std::vector<std::string>,std::vector<uint> > theAssociativeMap_){
 
  if(debug_) std::cout<<"ZbbEventContentAnalyzer::isTriggerOk"<<std::endl;
 
  Bool_t triggerisok = false;

  // get the fired trigger paths 
  const pat::TriggerPathRefVector acceptedTrigPaths(theTriggerEvent.acceptedPaths());
  
  // loop over selected trigger objects
  for ( pat::TriggerPathRefVector::const_iterator iTrig = acceptedTrigPaths.begin(); iTrig != acceptedTrigPaths.end(); ++iTrig ) {
    //loop over the map
    std::map<std::vector<std::string>,std::vector<uint> >::iterator it;  
    int triggernumber(0);
    for (it=theAssociativeMap_.begin() ; it !=theAssociativeMap_.end(); it++){
      triggernumber++; 
      // if event in the right run-range
      if( theRunNumber >=(*it).second[0] && theRunNumber <=(*it).second[1]){	 
	if(((*iTrig)->name())==((*it).first[0])){
	  TString trigName = ((*it).first[0]);
	  //	  std::cout<<"trigName: "<<trigName<<" trigNumber: "<<triggernumber<<std::endl;
	  if(trigName.Contains("Mu")){
	    histContainer_["muonTriggerBits"]->Fill(triggernumber-1); 
	  } else {
	    histContainer_["electronTriggerBits"]->Fill(triggernumber-1);
	  }
	  triggerisok =true;   // if there is a match theFilter_ becomes true and breaks nested for
	  break;
	}
      }// ends if on run-range
    }// ends loop over map
  }//  ends loop over fired triggerpaths
  return triggerisok;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//ZbbEventContentAnalyzer::isTriggerMatckOK
///////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t 
ZbbEventContentAnalyzer::isTriggerMatchOK(const reco::CompositeCandidate& ZCand,UInt_t theRunNumber, std::map<std::vector<std::string>,std::vector<uint> > theAssociativeMap_,Bool_t isMuChannel, Bool_t AndOrSwitch){

  if(debug_) std::cout<<"ZbbEventContentAnalyzer::isTriggerMatchOk"<<std::endl;

  Bool_t triggermatchisok = false;
  const reco::Candidate* lep0 = ZCand.daughter(0);
  const reco::Candidate* lep1 = ZCand.daughter(1);

  if(isMuChannel){
  
    const pat::Muon* muon0 = dynamic_cast<const pat::Muon*>(&(*lep0->masterClone()));
    const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(&(*lep1->masterClone()));   
 
    if(debug_){
      std::cout<<"Trigger object match size Filter (0) "<< muon0->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").size()<<std::endl;  
      std::cout<<"Trigger object match size Filter (1) "<< muon1->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").size()<<std::endl;  
      std::cout<<"Trigger object match size (Mu6) "<< muon0->triggerObjectMatchesByPath("HLT_DoubleMu6_v1").size()<<std::endl;   
      std::cout<<"Trigger object match size (Mu7v2) "<< muon0->triggerObjectMatchesByPath("HLT_DoubleMu7_v*").size()<<std::endl;   
      std::cout<<"Trigger object match size (Mu7v3) "<< muon1->triggerObjectMatchesByPath("HLT_DoubleMu7_v*").size()<<std::endl;       
      std::cout<<"Trigger object match size (v1) "<< muon0->triggerObjectMatchesByPath("HLT_Mu13_Mu8_v1").size()<<std::endl;   
      std::cout<<"Trigger object match size (v2) "<< muon0->triggerObjectMatchesByPath("HLT_Mu13_Mu8_v2").size()<<std::endl;    
      std::cout<<"Trigger object match size (v3) "<< muon0->triggerObjectMatchesByPath("HLT_Mu13_Mu8_v3").size()<<std::endl;  
    }    

   //loop over the map
    std::map<std::vector<std::string>,std::vector<uint> >::iterator it;  

    /// changes in the trigger matching (R.C. 2011)

    for ( it=theAssociativeMap_.begin() ; it !=theAssociativeMap_.end(); it++){

      int nmuHLThigh=0;
      int nmuHLTlow=0;

      // if event in the right run-range
      // std::cout<<"Matched path is  "<<(*it).first[0] <<std::endl;   
      // std::cout<<"Run number  "<<theRunNumber <<std::endl;   
      // std::cout<<"Range is "<< (*it).second[0] <<std::endl;   
      if( theRunNumber >=(*it).second[0] && theRunNumber <=(*it).second[1]){
	if(AndOrSwitch){
	  if(!muon0->triggerObjectMatchesByFilter((*it).first[2]).empty())nmuHLThigh++;
	  if(!muon1->triggerObjectMatchesByFilter((*it).first[2]).empty())nmuHLThigh++;
	  if(!muon0->triggerObjectMatchesByFilter((*it).first[3]).empty())nmuHLTlow++;
	  if(!muon1->triggerObjectMatchesByFilter((*it).first[3]).empty())nmuHLTlow++;
	  // bool passDiMu= (nmuHLThigh>=2||(nmuHLThigh>=1 && nmuHLTlow>=2));	
	  // bool passDiMu= (nmuHLThigh>=2||(nmuHLThigh>0 && nmuHLTlow>0));                
	  bool passDiMu = (muon0->triggerObjectMatches().size() > 0 || muon1->triggerObjectMatches().size() > 0  );  // Louvain's Recipe

	  if(//!muon0->triggerObjectMatchesByPath((*it).first[1]).empty() && // trying to remove bug
	     //!muon1->triggerObjectMatchesByPath((*it).first[1]).empty() &&
	     passDiMu){
	    triggermatchisok=true;     
	    break;
	  }
	}else { 
	  //if(!muon0->triggerObjectMatchesByFilter((*it).first).empty() || !muon1->triggerObjectMatchesByFilter((*it).first).empty()){	 
	  if(!muon0->triggerObjectMatchesByPath((*it).first[1]).empty() || !muon1->triggerObjectMatchesByPath((*it).first[1]).empty()){
	    triggermatchisok=true;   // if there is a match theFilter_ becomes true and breaks nested for
	    break;
	  }
	}//ends AndOr switch
      }// ends if on run-range 
    }// ends loop over map
  } else {
  
    const pat::Electron* ele0 = dynamic_cast<const pat::Electron*>(&(*lep0->masterClone()));
    const pat::Electron* ele1 = dynamic_cast<const pat::Electron*>(&(*lep1->masterClone()));
    
    // std::cout<<"electron0 matches size: "<<ele0->triggerObjectMatches().size() <<"/n"; 
    // std::cout<<"electron1 matches size: "<<ele1->triggerObjectMatches().size() <<std::endl;  

    //loop over the map
    std::map<std::vector<std::string>,std::vector<uint> >::iterator it;  
    for ( it=theAssociativeMap_.begin(); it !=theAssociativeMap_.end(); it++ ){
      int neleHLThigh=0;
      int neleHLTlow=0;

      // if event in the right run-range
      if( theRunNumber >=(*it).second[0] && theRunNumber <=(*it).second[1]){	 
	if(AndOrSwitch){
	  if(!ele0->triggerObjectMatchesByFilter((*it).first[2]).empty())neleHLThigh++;
	  if(!ele1->triggerObjectMatchesByFilter((*it).first[2]).empty())neleHLThigh++;
	  if(!ele0->triggerObjectMatchesByFilter((*it).first[3]).empty())neleHLTlow++;
	  if(!ele1->triggerObjectMatchesByFilter((*it).first[3]).empty())neleHLTlow++;
	  //  bool passDiEle= (neleHLThigh>=2||(neleHLThigh>=1 && neleHLTlow>=2));
	  bool passDiEle= (ele0->triggerObjectMatches().size() > 0 && ele1->triggerObjectMatches().size() > 0);  // Louvain's Recipe
	  
	  // std::cout<<"neleHLThigh: "<<neleHLThigh<<"neleHLTlow: "<<neleHLTlow<<std::endl; 

	  if(//!ele0->triggerObjectMatchesByPath((*it).first[1]).empty() && 
	     //!ele1->triggerObjectMatchesByPath((*it).first[1]).empty() && 
	     passDiEle){
	    triggermatchisok=true;     
	    break;
	  }
	} else {
	  if( !ele0->triggerObjectMatchesByPath((*it).first[1]).empty() || !ele1->triggerObjectMatchesByPath((*it).first[1]).empty()){
	    triggermatchisok=true;   // if there is a match theFilter_ becomes true and breaks nested for
	    break;
	  }
	} // ends AndOr switch
      }// ends if on run-range 
    }// ends loop over map
  }

  return triggermatchisok;  

}
      
///////////////////////////////////////////////////////////////////////////////////////////////////
//ZbbEventContentAnalyzer::fillTrendPlot (fills yield after each cut as a function of run number)
///////////////////////////////////////////////////////////////////////////////////////////////////

void 
ZbbEventContentAnalyzer::fillTrendPlot(EventCategory& ec, unsigned int theRunNumber,Bool_t isMuChannel){

  std::map<UInt_t,Bool_t> theEvtCategoryCut;

  Bool_t  cuts[12]            = {true,                                // 0 
				 ec.isTrigger(),                      // 1
				 ec.isLL(),                           // 2
				 ec.isLLTight(),                      // 3
				 ec.isLLTightAndTriggerMatched(),     // 4
				 ec.isZLL(),                          // 5
				 ec.isGoodJet(),                      // 6
				 ec.isJetVertexAssoc(),               // 7
				 ec.isJetBTagged(),                   // 8
				 ec.isExclusive(),                    // 9
				 ec.isThreeLeptons(),                 // 10
				 ec.isFourLeptons()};                 // 11

  Bool_t theCut=true;
  
  for(UInt_t cat=0; cat<=7; cat++){
    
    theCut = theCut&&cuts[cat];
    theEvtCategoryCut[cat] = theCut;
    
  }
  
  //separate case for mc truth matching
  theEvtCategoryCut[8]  = (theEvtCategoryCut.find(6)->second && ec.isJetBTagged());
  theEvtCategoryCut[9]  = (theEvtCategoryCut.find(8)->second && ec.isExclusive());
  theEvtCategoryCut[10] = (theEvtCategoryCut.find(8)->second && ec.isThreeLeptons());
  theEvtCategoryCut[11] = (theEvtCategoryCut.find(8)->second && ec.isFourLeptons());

  // returns the cut
  
  for(UInt_t nCut=0;nCut<12;nCut++){
    if(theEvtCategoryCut.find(nCut)->second){     
      //std::cout<<"nCut: "<<nCut<<" | name:"<<histContainer_[(evtVsRunMuName[nCut]).Data()]->GetName()<<" | title: "<<histContainer_[(evtVsRunMuName[nCut]).Data()]->GetTitle()<<std::endl;
      if(isMuChannel) fill(evtVsRunMuName[nCut].Data(),(float)theRunNumber,1);
      else fill(evtVsRunEleName[nCut].Data(),(float)theRunNumber,1);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//ZbbEventContentAnalyzer::PrintEvent
///////////////////////////////////////////////////////////////////////////////////////////////////

void 
ZbbEventContentAnalyzer::PrintEvent(const edm::Event& iEvent,Int_t category, Bool_t verbose){

  if(verbose){
    outfile_<< "================================================================="<<std::endl;
    outfile_<< "Run: "<<iEvent.eventAuxiliary().run()<<"\t  Lumi:"<<iEvent.eventAuxiliary().luminosityBlock()<<"\t Event: "<< iEvent.eventAuxiliary().id().event()<<std::endl;
    outfile_<< "Passes selection level:" <<category<<std::endl;
    time_t rawtime = iEvent.eventAuxiliary().time().unixTime();
    struct tm * timeinfo;
    timeinfo = localtime(&rawtime); 
    outfile_<< "Recorded on:"<< asctime (timeinfo) <<std::endl;
    outfile_<< "-----------------------------------------------------------------"<<std::endl;
  } else {
    outfile_<<iEvent.eventAuxiliary().run()<<":"<<iEvent.eventAuxiliary().luminosityBlock()<<":"<< iEvent.eventAuxiliary().id().event()<<std::endl;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//ZbbEventContentAnalyzer::endLuminosityBlock based on https://hypernews.cern.ch/HyperNews/CMS/get/muon/433/1/1.html
///////////////////////////////////////////////////////////////////////////////////////////////////

void
ZbbEventContentAnalyzer::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup) 
{

  // std::string ss="lumiProducer::RECO";
  edm::InputTag theLumiSummaryInputTag_(theLumiSummaryTag_);
  edm::Handle<LumiSummary> lumiSummary;
  iLumi.getByLabel(theLumiSummaryInputTag_,lumiSummary);
  
  // collect info only for the selected runs
  //  bool runInTheRange(false);
  //   std::vector<std::pair<uint, uint> >::iterator it;
  //   for (it=runRangeList_.begin(); it!=runRangeList_.end(); it++){
  //     if ( (*it).first<=iLumi.run() && iLumi.run()<=(*it).second ) {
  //       runInTheRange=true;
  //       break;
  //     }
  //   }
  
  //   if (!runInTheRange) return;
  
  if(lumiSummary->isValid()){    
    
    float dellumi = lumiSummary->intgDelLumi();
    float reclumi = dellumi*lumiSummary->liveFrac();

    totLS_ ++;
    totLumiDelivered_+=(double)dellumi;
    totLumiRecorded_+=(double)reclumi;     
    reclumibyrun_[iLumi.run()].first++;
    reclumibyrun_[iLumi.run()].second+=(double)reclumi;
  } else {
    
    // cout << "no valid lumi summary data run " << iLumi.run() << endl;
    
    float dellumi = lumiSummary->intgDelLumi();
    float reclumi = dellumi*lumiSummary->liveFrac();
    invalidLS_totLS_ ++;
    invalidLS_totLumiDelivered_+=(double)dellumi;
    invalidLS_totLumiRecorded_+=(double)reclumi;     
  }  
  
  return;
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ZbbEventContentAnalyzer);
