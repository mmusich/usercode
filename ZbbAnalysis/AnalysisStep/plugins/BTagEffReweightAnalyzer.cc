#include <iostream>

#include "TH1.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"
#include "RecoBTag/PerformanceDB/interface/BtagPerformance.h" 
#include "RecoBTag/Records/interface/BTagPerformanceRecord.h"

#include "CondFormats/PhysicsToolsObjects/interface/PerformancePayloadFromTable.h"
#include "CondFormats/PhysicsToolsObjects/interface/PerformanceWorkingPoint.h"


class BTagEffReweightAnalyzer : public edm::EDAnalyzer {

public:
  /// default constructor
  explicit BTagEffReweightAnalyzer(const edm::ParameterSet&);
  /// default destructor
  virtual ~BTagEffReweightAnalyzer();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  bool isBTaggedJet(const pat::Jet& jet,std::string theWP, std::string theAlgo);

 // jet flavour constants
  enum FLAVOUR {
    UDSG_JETS = 0,
    C_JETS,
    B_JETS,
    NONID_JETS,
    N_JET_TYPES
  };

  double getbEffScaleFactorFromPAS(FLAVOUR jetFlavour, const pat::Jet& jet);
  double getbEffScaleFactorFromDB(FLAVOUR jetFlavour, const pat::Jet& jet, const edm::EventSetup& iSetup);

  void getbEffScaleFactorFromDB_demo(const edm::EventSetup& iSetup);

  //mc switch
  bool ismc_;
  
  // input tags 
  edm::InputTag jetSrc_;

  // associative maps
  std::map<double,std::vector<double> > theAssociativeMapSFb_;

  // b-tag calibrations
  std::string beff_; // efficiency
  std::string leff_; // mistag

  //  BtagPerformance* lPerf;  

  // few control histos
  TH1F* h_jetMCflav_;

  TH1F* h_DISCR_;
  TH1F* h_DISCR_SELCTD_jetMCflav_;

  TH1F* h_DISCR_SELCTD_SFb_PAS_;    
  TProfile* p_DISCR_SELCTD_SFb_PAS_pt_;    
  TProfile* p_DISCR_SELCTD_SFb_PAS_eta_;    

  TH1F* h_DISCR_SELCTD_SFb_DB_;    
  TProfile* p_DISCR_SELCTD_SFb_DB_pt_;    
  TProfile* p_DISCR_SELCTD_SFb_DB_eta_;    

  TH1F* h_DISCR_SELCTD_bjetpt_;     
  TH1F* h_DISCR_SELCTD_bjeteta_;  
  TH1F* h_DISCR_SELCTD_bjetpt_rw_PAS_;     
  TH1F* h_DISCR_SELCTD_bjeteta_rw_PAS_;  
  TH1F* h_DISCR_SELCTD_bjetpt_rw_DB_;     
  TH1F* h_DISCR_SELCTD_bjeteta_rw_DB_;  
};


BTagEffReweightAnalyzer::BTagEffReweightAnalyzer(const edm::ParameterSet& iConfig) :
  ismc_(iConfig.getParameter<bool>("isMC")),
  jetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc"))
{

  // vpsets for SFb maps
  edm::VParameterSet sfb_ptrange_selection;
 
  sfb_ptrange_selection = iConfig.getParameter<edm::VParameterSet>("SFbPtRangeList");
  //intialize the associative maps
  for(edm::VParameterSet::const_iterator pset = sfb_ptrange_selection.begin(); pset != sfb_ptrange_selection.end(); pset++){
    double theSFb_ = pset->getUntrackedParameter<double>("scale_factor");
    std::vector<double> thePtRange_ = pset->getUntrackedParameter<std::vector<double> >("ptrange");
    theAssociativeMapSFb_[theSFb_]=thePtRange_; 
  }

  // from Db
  beff_ =  iConfig.getParameter<std::string>("CalibrationForBEfficiency");
  leff_ =  iConfig.getParameter<std::string>("CalibrationForLEfficiency");
} 

BTagEffReweightAnalyzer::~BTagEffReweightAnalyzer()
{
  //  if ( lPerf ) delete lPerf;
} 


void BTagEffReweightAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;

  // all jets
  h_jetMCflav_               = fs->make<TH1F>("MCflav","Jet MC flavour",5,-0.5,4.5);
  h_DISCR_                   = fs->make<TH1F>("DISCR","DISCRIMINATOR; discriminant b-jet",100,-10,10);

  // b-tagged jets
  h_DISCR_SELCTD_jetMCflav_  = fs->make<TH1F>("DISCR_SELCTD_MCflav","DISCRIMINATOR SELECTED MC flavour; MC flavour",5,-0.5,4.5); 
  h_DISCR_SELCTD_SFb_PAS_       = fs->make<TH1F>("DISCR_SELCTD_SF_PAS_","SF_b DISCRIMINATOR from PAS",100,0.5,1.5);    
  h_DISCR_SELCTD_SFb_DB_        = fs->make<TH1F>("DISCR_SELCTD_SF_DB_","SF_b DISCRIMINATOR from DB",100,0.5,1.5);    

  float SFb_pTbins[12]={0, 15, 20, 30, 40, 50, 60, 70, 80, 100, 120, 150};
  p_DISCR_SELCTD_SFb_PAS_pt_     = fs->make<TProfile>("DISCR_SELCTD_SF_PAS_pt","SF_b vs b-jet pT for DISCRIMINATOR from PAS; pT b-jet; SF_b for DISCR",11,SFb_pTbins);
  p_DISCR_SELCTD_SFb_DB_pt_     = fs->make<TProfile>("DISCR_SELCTD_SF_DB_pt","SF_b vs b-jet pT for DISCRIMINATOR from DB; pT b-jet; SF_b for DISCR",11,SFb_pTbins);
  float SFb_etabins[3]={0,1.2,2.4};
  p_DISCR_SELCTD_SFb_PAS_eta_    = fs->make<TProfile>("DISCR_SELCTD_PAS_SF_eta","SF_b vs b-jet eta for DISCRIMINATOR from PAS; eta b-jet; SF_b for DISCR",2,SFb_etabins);
  p_DISCR_SELCTD_SFb_DB_eta_    = fs->make<TProfile>("DISCR_SELCTD_DB_SF_eta","SF_b vs b-jet eta for DISCRIMINATOR from DB; eta b-jet; SF_b for DISCR",2,SFb_etabins);

  h_DISCR_SELCTD_bjetpt_     = fs->make<TH1F>("DISCR_SELCTDjetpt" ,"jet p_{T};  jet p_{T} (GeV)",50,15,215);
  h_DISCR_SELCTD_bjeteta_    = fs->make<TH1F>("DISCR_SELCTDjeteta","jet #eta;   jet #eta",35,-3.5,3.5);
  h_DISCR_SELCTD_bjetpt_rw_PAS_  = fs->make<TH1F>("DISCR_SELCTDjetpt_rw_PAS" ,"jet p_{T} post reweight PAS; jet p_{T} (GeV)",50,15,215);
  h_DISCR_SELCTD_bjeteta_rw_PAS_ = fs->make<TH1F>("DISCR_SELCTDjeteta_rw_PAS","jet #eta post reweight DB; jet #eta",35,-3.5,3.5);
  h_DISCR_SELCTD_bjetpt_rw_DB_  = fs->make<TH1F>("DISCR_SELCTDjetpt_rw_DB" ,"jet p_{T} post reweight PAS; jet p_{T} (GeV)",50,15,215);
  h_DISCR_SELCTD_bjeteta_rw_DB_ = fs->make<TH1F>("DISCR_SELCTDjeteta_rw_DB","jet #eta post reweightt DB; jet #eta",35,-3.5,3.5);

}

void BTagEffReweightAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  getbEffScaleFactorFromDB_demo(iSetup);

  // get jet collection
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jetSrc_,jets);

  for(edm::View<pat::Jet>::const_iterator jetIt=jets->begin(); jetIt!=jets->end(); ++jetIt){
 
    FLAVOUR flavour(NONID_JETS);
    if ( ismc_  ) {	
      switch(std::abs(jetIt->partonFlavour())) {
      case 1:
      case 2:
      case 3:
      case 21:
	flavour = UDSG_JETS;
	break;
      case 4:
	flavour = C_JETS;
	break;
      case 5:
	flavour = B_JETS;
	break;
      default:
	flavour = NONID_JETS;
      }
      h_jetMCflav_->Fill(flavour); //all jets
    }

    h_DISCR_->Fill(jetIt->bDiscriminator("simpleSecondaryVertexHighEffBJetTags")); // all jets

    double weight[2]={0,0};
    if( isBTaggedJet((*jetIt),"HEM","SSV") ) {  // if b-tagged
      if ( ismc_ )
	{
	  h_DISCR_SELCTD_jetMCflav_->Fill(flavour);
	  // get the SFb
	  weight[0] = getbEffScaleFactorFromPAS(flavour,(*jetIt)); 
	  weight[1] = getbEffScaleFactorFromDB(flavour,(*jetIt),iSetup); 

	  p_DISCR_SELCTD_SFb_PAS_pt_->Fill(jetIt->pt(),weight[0]);
	  p_DISCR_SELCTD_SFb_PAS_eta_->Fill(fabs(jetIt->eta()),weight[0]);
	  h_DISCR_SELCTD_SFb_PAS_->Fill(weight[0]);

	  p_DISCR_SELCTD_SFb_DB_pt_->Fill(jetIt->pt(),weight[1]);
	  p_DISCR_SELCTD_SFb_DB_eta_->Fill(fabs(jetIt->eta()),weight[1]);
	  h_DISCR_SELCTD_SFb_DB_->Fill(weight[1]);
	}
      h_DISCR_SELCTD_bjetpt_->Fill(jetIt->pt());      
      h_DISCR_SELCTD_bjeteta_->Fill(jetIt->eta());	   
      h_DISCR_SELCTD_bjetpt_rw_PAS_->Fill(jetIt->pt(),weight[0]);      
      h_DISCR_SELCTD_bjeteta_rw_PAS_->Fill(jetIt->eta(),weight[0]);	   
      h_DISCR_SELCTD_bjetpt_rw_DB_->Fill(jetIt->pt(),weight[1]);      
      h_DISCR_SELCTD_bjeteta_rw_DB_->Fill(jetIt->eta(),weight[1]);	   

    } // close if b-tagged

  } // close loop on jets
}


void BTagEffReweightAnalyzer::endJob()
{
}



bool BTagEffReweightAnalyzer::isBTaggedJet(const pat::Jet& jet,std::string theWP, std::string theAlgo){  

  bool isbjet = false;

  if(theAlgo=="SSV"){
    if(theWP=="HEM") {
      if(jet.bDiscriminator("simpleSecondaryVertexHighEffBJetTags")>1.74) isbjet = true;
    } else if(theWP=="HPT") {
      if(jet.bDiscriminator("simpleSecondaryVertexHighPurBJetTags")>2.00) isbjet = true;
    } else {
      std::cout<<"isBTaggedJet Error unforeseen working point for b-tagging. Use HEM or HPT"<<std::endl;
    }
  } else if (theAlgo=="TC"){
    if(theWP=="HEL") {
      if(jet.bDiscriminator("trackCountingHighEffBJetTags")>1.70) isbjet = true;
    } else if(theWP=="HEM" ) {
      if(jet.bDiscriminator("trackCountingHighEffBJetTags")>3.30) isbjet = true;
    } else if(theWP=="HPL") {
      if(jet.bDiscriminator("trackCountingHighPurBJetTags")>1.19) isbjet = true;  
    } else if(theWP=="HPT" ) {
      if(jet.bDiscriminator("trackCountingHighPurBJetTags")>3.41) isbjet = true;
    } else {
      std::cout<<"isBTaggedJet Error: unforeseen working point for b-tagging. Use HEM or HPT"<<std::endl;
    }
  } else { 
    std::cout<<"isBTaggedJet Error: unforeseen algo for b-tagging. Use SSV or TC"<<std::endl;
  }
  return isbjet;
}

double BTagEffReweightAnalyzer::getbEffScaleFactorFromPAS(FLAVOUR jetFlavour, const pat::Jet& jet){

  double w(1);
  
  // it would be better determine the jet flavour inside this method....
  if ( jetFlavour != B_JETS ) return w; 
  double ETjet=jet.pt();
  double AbsEtajet=fabs(jet.eta());

  //loop over the map
  std::map<double,std::vector<double> >::iterator it;  
  for ( it=theAssociativeMapSFb_.begin() ; it !=theAssociativeMapSFb_.end(); it++ ){
    // if pT in the range
    if( ETjet >=(*it).second[0] && ETjet <=(*it).second[1]){
      w = (*it).first;
      break;
    }// ends if on pT range
  }// ends loop over map  


  return w;
}

double BTagEffReweightAnalyzer::getbEffScaleFactorFromDB(FLAVOUR jetFlavour, const pat::Jet& jet, const edm::EventSetup& iSetup){

  double w(1);

 
  double ETjet=jet.pt();
  double AbsEtajet=fabs(jet.eta());

  //  std::cout << "ET/eta " << ETjet << " / " << AbsEtajet << std::endl;
  edm::ESHandle<BtagPerformance> perfEFFH;
  iSetup.get<BTagPerformanceRecord>().get(beff_,perfEFFH);
  const BtagPerformance & pbeff = *(perfEFFH.product());

  edm::ESHandle<BtagPerformance> perfMISH;
  iSetup.get<BTagPerformanceRecord>().get(leff_,perfMISH);
  const BtagPerformance & pleff = *(perfMISH.product());

  //  std::cout <<" Discriminant is "<<pbeff.workingPoint().discriminantName()<<std::endl;
  // std::cout << "Working point: " << pbeff.workingPoint().cut() << std::endl;


  BinningPointByMap p;
  p.reset();
  p.insert(BinningVariables::JetAbsEta,AbsEtajet);
  p.insert(BinningVariables::JetEt,ETjet);


  // it would be better determine the jet flavour inside this method....
  if ( jetFlavour == B_JETS ) {
    if ( pbeff.isResultOk(PerformanceResult::BTAGBEFFCORR,p) )
      w=pbeff.getResult(PerformanceResult::BTAGBEFFCORR,p);
  } else {
    if ( pleff.isResultOk(PerformanceResult::BTAGLEFF,p) )
      w=pleff.getResult(PerformanceResult::BTAGLEFF,p);
  }

  return w;
}


//
void BTagEffReweightAnalyzer::getbEffScaleFactorFromDB_demo(const edm::EventSetup& iSetup){

  edm::ESHandle<BtagPerformance> perfEFFH;
  iSetup.get<BTagPerformanceRecord>().get(beff_,perfEFFH);
  const BtagPerformance & pbeff = *(perfEFFH.product());

  edm::ESHandle<BtagPerformance> perfMISH;
  iSetup.get<BTagPerformanceRecord>().get(leff_,perfMISH);
  const BtagPerformance & pleff = *(perfMISH.product());

  //  std::cout <<" Discriminant is "<<pbeff.workingPoint().discriminantName()<<std::endl;
  // std::cout << "Working point: " << pbeff.workingPoint().cut() << std::endl;


  BinningPointByMap p;


  double weffb(1.), weffl(1.);
  //eta  (0-0.8, 0.8-1.6, 1.6-2.4)
  double demoEta[3]={0.4, 1.2, 2.0};
  double demoEt[29]={15, 25, 35, 45, 55, 65, 75, 85, 95, 105,115,125,135,145,155,165,175,185,195, 205,215,225,235,245,255,265,275,285,295};
  for (int i=0; i<3; i++){
    for (int j=0; j<29; j++){
      p.reset();
      p.insert(BinningVariables::JetAbsEta,demoEta[i]);
      p.insert(BinningVariables::JetEt,demoEt[j]);      
      if ( pbeff.isResultOk(PerformanceResult::BTAGBEFF,p) ) weffb=pbeff.getResult(PerformanceResult::BTAGBEFF,p);
      if ( pleff.isResultOk(PerformanceResult::BTAGLEFF,p) )  weffl=pleff.getResult(PerformanceResult::BTAGLEFF,p);
      std::cout << "ET/eta " << demoEta[i] << " / " << demoEt[j] << " EFFb: " << weffb << " EFFl: " << weffl <<std::endl;	   	
    }
  }
  return;

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BTagEffReweightAnalyzer);



