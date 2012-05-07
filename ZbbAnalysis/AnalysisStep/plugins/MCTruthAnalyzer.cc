// -*- C++ -*-
//
// Package:    MCTruthAnalyzer
// Class:      MCTruthAnalyzer
// 
/**\class MCTruthAnalyzer MCTruthAnalyzer.cc ZbbAnalysis/AnalysisStep/src/MCTruthAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Authors:  Stefano Casasso, Marco Musich
//          Created:  Thu Mar 24 14:59:11 CET 2011
// $Id: MCTruthAnalyzer.cc,v 1.1 2011/12/05 13:34:46 musich Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Common/interface/GetProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
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
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

//
// class declaration
//
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMath.h"

using namespace edm;
using namespace std;
using namespace reco;

class MCTruthAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MCTruthAnalyzer(const edm::ParameterSet&);
      ~MCTruthAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual std::vector<math::XYZTLorentzVectorD> sort4MomentaByPt(std::vector<math::XYZTLorentzVectorD> leptonMomenta_);
      virtual std::vector<const reco::Candidate*> sortLeptonsByPt(std::vector<const reco::Candidate*> leptons_);
      virtual Double_t CosThetaStar(std::vector<const reco::Candidate*> sortedGenLeptons);
      virtual Double_t CostCS(std::vector<const reco::Candidate*> sortedGenLeptons);
      virtual std::pair<bool,bool> WhichFlavour(std::vector<const reco::Candidate*> theGenMuons,std::vector<const reco::Candidate*> theGenElectrons);
      virtual std::pair<const reco::Candidate*,const reco::Candidate*> getTheZLeptons(std::vector<const reco::Candidate*> theGenMuons,std::vector<const reco::Candidate*> theGenElectrons);
      virtual bool hasBottom(const reco::Candidate &c);
      virtual bool hasCharm(const reco::Candidate &c); 
      virtual void endJob() ;

  edm::InputTag genJetSrc_;
  edm::InputTag genParticles_;
  edm::InputTag MuonCollection_;
  edm::InputTag ElectronCollection_;
  edm::InputTag JetCollection_;
  edm::InputTag ZmmCollection_;
  bool isMCatNLO_; 

  //counters
  TH1F *h1NMu_;
  TH1F *h1NEle_;

  //MC weight
  TH1F *hMCWeight_;
  TH1F *hMCWeightZoom_;

  TH1F *h1NRecoJets_;
  TH2F *h1NRecoJetsVsWeight_;
  TH1F *h1NRecoJetsW_;

  TH1F *h1NGenJets_;
  TH2F *h1NGenJetsVsWeight_;
  TH1F *h1NGenJetsW_;

  //Mu mother Id
  TH1F *h1MotherIdMuon_;
  //Ele mother Id
  TH1F *h1MotherIdElectron_;
  
  // unweighted Spectra
  TH1F *h1PtSoftMuon_;
  TH1F *h1PtLeadingMuon_;
  TH1F *h1PtSoftElectron_;
  TH1F *h1PtLeadingElectron_;  
  TH1F *h1DrMinMJ_;
  TH1F *h1PtMuFromHF_;
  TH1F *h1PtAllJets_;    
  TH1F *h1PtSoftJet_;    
  TH1F *h1PtLeadingJet_; 

  //weighted spectra
  TH1F *h1PtSoftMuonW_;
  TH1F *h1PtLeadingMuonW_;
  TH1F *h1PtSoftElectronW_;
  TH1F *h1PtLeadingElectronW_;  
  TH1F *h1DrMinMJW_;
  TH1F *h1PtMuFromHFW_;
  TH1F *h1PtSoftJetW_;    
  TH1F *h1PtAllJetsW_;
  TH1F *h1PtLeadingJetW_; 
  
  //MC Truth for B-hadrons
  TH1F *h1FirstbHadronId_;
  TH1F *h1bHadronStatus_;
  TH1F *h1bHadronMotherId_;
  TH1F *h1FirstbHadronMotherId_;
  TH1F *h1nBHadrons_;

  TH1F *h1nLastBHadronsPerJet_;

  TH1F *h1BHadronsPt_;
  TH1F *h1BHadronsEta_;

  TH1F *h1nLastBHadrons_;

  TH1F *h1LastBHadronsPt_;
  TH1F *h1LastBHadronsEta_;
  TH1F *h1LastBHadronsMass_;

  TH1F *h1LastBHadronsPtTightAssoc_;
  TH1F *h1LastBHadronsEtaTightAssoc_;
  TH1F *h1LastBHadronsMassTightAssoc_;
  
  TH1F *h1PoorlyAssociatedLastBHadronsPt_;
  TH1F *h1PoorlyAssociatedLastBHadronsEta_;

  TH1F* hDeltaRLastBHadGenJet_;
  TH1F* hDeltaRLastBHadGoodGenJet_;

  TH1F* hDeltaRLastBHadGenJet_ptHat_020_;
  TH1F* hDeltaRLastBHadGenJet_ptHat_2025_;
  TH1F* hDeltaRLastBHadGenJet_ptHat_2530_;
  TH1F* hDeltaRLastBHadGenJet_ptHat_GT30_;

  //MC Truth for z
  TH1F* h_GENP_pt_Z_;  
  TH1F* h_GENP_eta_Z_; 
  TH1F* h_GENP_mass_Z_;
  			      
  TH1F* h_GENP_pt_ZfromLept_;  
  TH1F* h_GENP_eta_ZfromLept_; 
  TH1F* h_GENP_mass_ZfromLept_;
       		
  TH1F* h_GENP_pt_ZfromMuons_; 
  TH1F* h_GENP_eta_ZfromMuons_;
  TH1F* h_GENP_mass_ZfromMuons_;
  
  TH1F* h_GENP_pt_ZfromElectrons_;   
  TH1F* h_GENP_eta_ZfromElectrons_;  
  TH1F* h_GENP_mass_ZfromElectrons_; 

  TH1F* h_GENP_pt_Lept1_fromZ_;        
  TH1F* h_GENP_eta_Lept1_fromZ_;       
  				
  TH1F* h_GENP_pt_Muons1_fromZ_;       
  TH1F* h_GENP_eta_Muons1_fromZ_;      
  				
  TH1F* h_GENP_pt_Electrons1_fromZ_;   
  TH1F* h_GENP_eta_Electrons1_fromZ_;  
  				
  TH1F* h_GENP_pt_Lept2_fromZ_;        
  TH1F* h_GENP_eta_Lept2_fromZ_;       
  				
  TH1F* h_GENP_pt_Muons2_fromZ_;       
  TH1F* h_GENP_eta_Muons2_fromZ_;      
  				
  TH1F* h_GENP_pt_Electrons2_fromZ_;   
  TH1F* h_GENP_eta_Electrons2_fromZ_;  

  TH1F* h_GENP_zllCosThetaCS_;    
  TH1F* h_GENP_zllCosThetaStar_;
				      
  TH1F* h_GENP_mll50_pt_ZfromLept_;          
  TH1F* h_GENP_mll50_eta_ZfromLept_;         
  TH1F* h_GENP_mll50_mass_ZfromLept_;        
				      
  TH1F* h_GENP_mll50_pt_ZfromMuons_;         
  TH1F* h_GENP_mll50_eta_ZfromMuons_;        
  TH1F* h_GENP_mll50_mass_ZfromMuons_;       
  				      
  TH1F* h_GENP_mll50_pt_ZfromElectrons_;     
  TH1F* h_GENP_mll50_eta_ZfromElectrons_;    
  TH1F* h_GENP_mll50_mass_ZfromElectrons_;   
				      	      
  TH1F* h_GENP_mll50_pt_Lept1_fromZ_;        
  TH1F* h_GENP_mll50_eta_Lept1_fromZ_;       
				      
  TH1F* h_GENP_mll50_pt_Muons1_fromZ_;       
  TH1F* h_GENP_mll50_eta_Muons1_fromZ_;      
				      
  TH1F* h_GENP_mll50_pt_Electrons1_fromZ_;   
  TH1F* h_GENP_mll50_eta_Electrons1_fromZ_;  
				      
  TH1F* h_GENP_mll50_pt_Lept2_fromZ_;        
  TH1F* h_GENP_mll50_eta_Lept2_fromZ_;       
				      
  TH1F* h_GENP_mll50_pt_Muons2_fromZ_;       
  TH1F* h_GENP_mll50_eta_Muons2_fromZ_;      
				      
  TH1F* h_GENP_mll50_pt_Electrons2_fromZ_;   
  TH1F* h_GENP_mll50_eta_Electrons2_fromZ_;  
				      
  TH1F* h_GENP_mll50_zllCosThetaCS_;         
  TH1F* h_GENP_mll50_zllCosThetaStar_;       
 
  // genJet b plots	      
  TH1F* h_GENP_pt_b_;   
  TH1F* h_GENP_eta_b_;  
  TH1F* h_GENP_mass_b_; 
  TH1F* h_GENP_maxDistance_b_; 

  TH1F* h_GENP_pt_b_tightmatch_;       
  TH1F* h_GENP_eta_b_tightmatch_;      
  TH1F* h_GENP_mass_b_tightmatch_;    
  TH1F* h_GENP_maxDistance_b_tightmatch_;

  TH1F* h_GENP_mass_b_tightmatch_maxDistCut_;
  
  TH1F* h_GENP_pt_b_extratight_;     
  TH1F* h_GENP_eta_b_extratight_;    
  TH1F* h_GENP_mass_b_extratight_;   
  TH1F* h_GENP_maxDistance_b_extratight_; 

  TH1F* h1BgenJetNConst_;     
  TH1F* h1BgenJetConstPdgId_; 

  TH2F* h2_GENP_massb_vs_etab_;        
  TH2F* h2_GENP_massb_vs_ptb_;         
  TH2F* h2_GENP_massb_vs_maxDistanceb_;
			       
  TProfile* p_GENP_massb_vs_etab_;         
  TProfile* p_GENP_massb_vs_ptb_;          
  TProfile* p_GENP_massb_vs_maxDistanceb_;

  // all genJets
  TH1F* h1genJetNConst_;     
  TH1F* h1genJetConstPdgId_;
  TH1F* h_GENP_pt_all_;         
  TH1F* h_GENP_eta_all_;        
  TH1F* h_GENP_mass_all_;     
  TH1F* h_GENP_maxDistance_all_; 
  TH1F* h_GENP_MymaxDistance_all_;

  TH1F* h1NGoodGenJets_;        

  // total kinematics of the event
  TH1F* h_GENP_pxTot_;                 
  TH1F* h_GENP_pyTot_;                
  TH1F* h_GENP_pzTot_;                 
  TH1F* h_GENP_ETot_;                  
  TH1F* h_GENP_DeltaEcmTotalSqrtS_;    
  TH2F* h_GENP_DeltaEcmTotalSqrtS_Vs_NPart_;

  //Counters
  int theMuonCounter_,theElectronCounter_, theJetCounter_;
  int nTotEvts_;

  //anglular separation of b-hadrons
  TH1F *h1LastBhadronDeltaR_; 
  TH1F *h1LastBhadronGenJetMatchedDeltaR_;

  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MCTruthAnalyzer::MCTruthAnalyzer(const edm::ParameterSet& iConfig)
{
  genParticles_      =iConfig.getParameter<InputTag>("src");
  MuonCollection_    =iConfig.getParameter<InputTag>("MuonCollection");
  ElectronCollection_=iConfig.getParameter<InputTag>("ElectronCollection");
  JetCollection_     =iConfig.getParameter<InputTag>("JetCollection");
  ZmmCollection_     =iConfig.getParameter<InputTag>("ZmmCollection");
  genJetSrc_         =iConfig.getParameter<InputTag>("genJetSrc");
  isMCatNLO_         =iConfig.getUntrackedParameter<bool>("isMCatNLO");
  //now do what ever initialization is needed
}


MCTruthAnalyzer::~MCTruthAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
MCTruthAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  theMuonCounter_= 0;
  theElectronCounter_ = 0;
  theJetCounter_=0;
  nTotEvts_ ++;

   // get gen particle candidates
  edm::Handle<GenParticleCollection> genParticlesCollection;
  iEvent.getByLabel(genParticles_, genParticlesCollection);
  
  // get patMuonsWithTrigger
  edm::Handle<pat::MuonCollection> muonsHandle;
  iEvent.getByLabel(MuonCollection_, muonsHandle);
  
  // get patElectrons
  edm::Handle<pat::ElectronCollection> electronsHandle;
  iEvent.getByLabel(ElectronCollection_, electronsHandle);
  
  // get patJets
  edm::Handle<pat::JetCollection> jetsHandle;
  iEvent.getByLabel(JetCollection_, jetsHandle);
  
  // get Zmm candidates
  Handle<reco::CompositeCandidateCollection> zmmHandle;
  iEvent.getByLabel(ZmmCollection_, zmmHandle);
  
  // get event weight 
  edm::Handle<GenEventInfoProduct> genEventInfoHandle;
  iEvent.getByLabel("generator",genEventInfoHandle);
  double wMCtmp = genEventInfoHandle->weight();

  double wMC(1.);
  
  if(isMCatNLO_) {
    wMCtmp > 0 ?  wMC=1. : wMC=-1.;
  } else {
    wMC = wMCtmp;
  }

  //std::cout<<"genEventInfoHandle->weight(): "<<wMC<<" genEventInfoHandle->weights()[0]: "<<genEventInfoHandle->weights()[0]<<std::endl;

  // get LHE Event
  edm::Handle<LHEEventProduct> lhevtCollection;
  lhevtCollection.clear();

  if( lhevtCollection.isValid () ) iEvent.getByType( lhevtCollection );

  // get genJets
  edm::Handle<reco::GenJetCollection> genJetsHandle;
  iEvent.getByLabel(genJetSrc_,genJetsHandle);
  const reco::GenJetCollection & genJets = *(genJetsHandle.product());

  h1NGenJets_->Fill(genJets.size());  
  h1NGenJetsVsWeight_->Fill(genJets.size(),wMC);  
  h1NGenJetsW_->Fill(genJets.size(),wMC); 

  //**************
  //Event Weight
  //**************

  hMCWeight_->Fill(wMC);
  hMCWeightZoom_->Fill(wMC);

  //***************
  //Analyze muons
  //***************
  
  //float PtSoftMuon_ = 0;
  //float PtLeadingMuon_ = 0;
  pat::MuonCollection::const_iterator SoftMuon,LeadingMuon,muon;
  
  if( muonsHandle->size() > 1){
    //Loop on MuonCollection
    for( muon = muonsHandle->begin();muon != muonsHandle->end(); ++ muon ){
      if (muon == muonsHandle->begin() || muon->pt() < SoftMuon->pt())SoftMuon = muon;
      if (muon == muonsHandle->begin() || muon->pt() > LeadingMuon->pt())LeadingMuon = muon;
    }//end loop on MuonCollection
    
    h1PtSoftMuon_->Fill(std::min(SoftMuon->pt(),99.9));
    h1PtLeadingMuon_->Fill(std::min(LeadingMuon->pt(),149.9));
    
    h1PtSoftMuonW_->Fill(std::min(SoftMuon->pt(),99.9),wMC);
    h1PtLeadingMuonW_->Fill(std::min(LeadingMuon->pt(),149.9),wMC);

  }//end if muonsHandle size
  
  //***************
  //Analyze electrons
  //***************
  
  //float PtSoftElectron_ = 0;
  //float PtLeadingElectron_ = 0;
  pat::ElectronCollection::const_iterator SoftElectron,LeadingElectron,electron;
  
  if( electronsHandle->size() > 1){
    //Loop on ElectronCollection
    for( electron = electronsHandle->begin();electron != electronsHandle->end(); ++ electron ){
      if (electron == electronsHandle->begin() || electron->pt() < SoftElectron->pt())SoftElectron = electron;
      if (electron == electronsHandle->begin() || electron->pt() > LeadingElectron->pt())LeadingElectron = electron;
    }//end loop on ElectronCollection
    
    h1PtSoftElectron_->Fill(std::min(SoftElectron->pt(),99.9));
    h1PtLeadingElectron_->Fill(std::min(LeadingElectron->pt(),149.9));
    
    h1PtSoftElectronW_->Fill(std::min(SoftElectron->pt(),99.9),wMC);
    h1PtLeadingElectronW_->Fill(std::min(LeadingElectron->pt(),149.9),wMC);

  }//end if electronsHandle size


  //***************
  //Analyze jets
  //***************

  h1NRecoJets_->Fill(jetsHandle->size());
  h1NRecoJetsVsWeight_->Fill(jetsHandle->size(),wMC);  
  h1NRecoJetsW_->Fill(jetsHandle->size(),wMC);

  pat::JetCollection::const_iterator SoftJet,LeadingJet,jet;
  if( jetsHandle->size() > 0 ){
    for (jet = jetsHandle->begin();jet != jetsHandle->end(); ++jet){//loop on jet collection
      if (jet == jetsHandle->begin() || jet->pt() < SoftJet->pt())SoftJet = jet;
      if (jet == jetsHandle->begin() || jet->pt() > LeadingJet->pt())LeadingJet = jet;
      
      h1PtAllJets_->Fill(std::min(jet->pt(),199.9)); 
      h1PtAllJetsW_->Fill(std::min(jet->pt(),199.9),wMC);

    }//end loop on jet collection
    
    h1PtSoftJet_->Fill(std::min(SoftJet->pt(),99.9));    
    h1PtLeadingJet_->Fill(std::min(LeadingJet->pt(),199.9)); 
    
    h1PtSoftJetW_->Fill(std::min(SoftJet->pt(),99.9),wMC);    
    h1PtLeadingJetW_->Fill(std::min(LeadingJet->pt(),199.9),wMC); 
  }
  
  //***************
  //Analyze Z
  //***************

  pat::MuonCollection::const_iterator Muon1st,Muon2nd;
  if ( zmmHandle->size() >= 1 && muonsHandle->size() > 2 ){
    for ( muon = muonsHandle->begin();muon != muonsHandle->end(); ++ muon ){
      float PtMuon_ = muon->pt();
      if ( PtMuon_ < LeadingMuon->pt() ) {
	float DrMin_ = 0;
	for (jet = jetsHandle->begin();jet != jetsHandle->end(); ++jet){//loop on jet collection
	  float Dphi_ = abs(muon->phi()-jet->phi());
	  float Deta_ = abs(muon->eta()-jet->eta());
	  float Dr_ = TMath::Sqrt(Dphi_*Dphi_ + Deta_*Deta_);
	  if( DrMin_ == 0 || Dr_ < DrMin_)DrMin_ = Dr_;
	}//end loop on jets
	h1DrMinMJ_->Fill(DrMin_);
	h1DrMinMJW_->Fill(DrMin_,wMC);
      }//end if muon pt 
    }//end loop on muons
  }//end if


  //********************************
  // TotalKinematics analysis
  //********************************

  Bool_t verbose_ = false;

  float nEcms = 0.;
  unsigned int nInit = 0;

  std::vector<float> p4tot(4,0.);
  unsigned int nPart = 0;
  //Double_t tolerance_ = 0.5;

  //********************************
  // MC Truth Analysis
  //********************************

  // initialize pT and eta of the leading b/Z
  std::vector<const reco::Candidate*> theGenMuons;
  std::vector<const reco::Candidate*> theGenElectrons;
  
  math::XYZTLorentzVectorD p4ZGEN(0,0,0,0); 
  math::XYZTLorentzVectorD p4ZGENFromLept(0,0,0,0); 

  Int_t nBHad(0);
  Int_t nLastBHad(0);

  //Loop on genParticles
  for( GenParticleCollection::const_iterator genp = genParticlesCollection->begin();genp != genParticlesCollection->end(); ++ genp ) {  // loop over GEN particles
  
    //--------------------------------------- total kinematics ----------------------------------
    if ( nInit < 3 && (*genp).status() == 3 && (*genp).pdgId() == 2212 ) {
      nInit++;
      nEcms += (*genp).energy();
    }
    if ( (*genp).status() == 1) { 
      nPart++;
      if (verbose_) {
        std::cout << "Status 1 part # " << std::setw(4) << std::fixed << nPart 
                  << std::setw(14) << std::fixed << (*genp).pdgId() 
                  << std::setw(14) << std::fixed << (*genp).px() 
                  << std::setw(14) << std::fixed << (*genp).py() 
                  << std::setw(14) << std::fixed << (*genp).pz() << std::endl;
      }
      p4tot[0] += (*genp).px();
      p4tot[1] += (*genp).py();
      p4tot[2] += (*genp).pz();
      p4tot[3] += std::sqrt( (*genp).px()*(*genp).px() + 
			     (*genp).py()*(*genp).py() + 
			     (*genp).pz()*(*genp).pz() + 
			     (*genp).mass()*(*genp).mass()) ; 
    }
    
    //do a copy of the genParticle, not to skip anythin in the loop
    const reco::Candidate* p = &(*genp);
    
    //------------------------------------------- muons -------------------------------------------
    if( abs( p->pdgId() )==13 && ( p->status())==1 ){
      
      if(p->pt()>0.5){
	theGenMuons.push_back(p);
      }
      
      theMuonCounter_++;
      
      //Check muon's mother
      while (p->mother()!=0 && abs(p->mother()->pdgId()) == 13) {
	p = p->mother();
      }
      //if(p->mother()->pdgId() >= 1 && p->mother()->pdgId() <= 6)  cout<<"Hi, I'm a muon and my mother has a pdgId = "<<p->mother()->pdgId()<<endl;
      
      //Fill muon's mother id histogram
      h1MotherIdMuon_->Fill(p->mother()->pdgId());
           
      //Now get muons from heavy flavour
      if( (abs(p->mother()->pdgId()) > 411 && abs(p->mother()->pdgId()) < 435)   ||// c-mesons
	  (abs(p->mother()->pdgId()) > 441 && abs(p->mother()->pdgId()) < 445)   ||// ccbar
	  (abs(p->mother()->pdgId()) > 511 && abs(p->mother()->pdgId()) < 545)   ||// b-mesons
	  (abs(p->mother()->pdgId()) > 551 && abs(p->mother()->pdgId()) < 557)   ||// bbbar
	  (abs(p->mother()->pdgId()) > 4122 && abs(p->mother()->pdgId()) < 4444) ||// c-barions
	  (abs(p->mother()->pdgId()) > 5122 && abs(p->mother()->pdgId()) < 5554)   // b-barions
	  ) {
	h1PtMuFromHF_->Fill(p->pt());
	h1PtMuFromHFW_->Fill(p->pt(),wMC);
      }
    } //if muons status 1
    

    //------------------------------------------- electrons -------------------------------------------
    if( abs( p->pdgId() )==11 && ( p->status())==3 ){
      
      if(p->pt()>0.5){
	theGenElectrons.push_back(p);
      }
      
      theElectronCounter_++;
      
      //Check electron's mother
      while (p->mother()!=0 && abs(p->mother()->pdgId()) == 11) {
	p = p->mother();
      }
      //Fill electron's mother id histogram
      h1MotherIdElectron_->Fill(p->mother()->pdgId());
    }//end if electrons status 1   
  }
 
  //  if ( std::abs(p4tot[0]) > tolerance_ || std::abs(p4tot[1]) > tolerance_ || std::abs(p4tot[2]) > tolerance_ || std::abs(p4tot[3]-nEcms) > tolerance_ ) {
  //     // std::cout << "Initial sqrt(s) = " << nEcms << std::endl;
  //     //for (unsigned int i=0; i<4; i++) {
  //     // std::cout << "p4tot["<<i<<"] = " << p4tot[i] << std::endl;
  //     //}
  //   }
  
  h_GENP_pxTot_->Fill(p4tot[0]);            
  h_GENP_pyTot_->Fill(p4tot[1]);            
  h_GENP_pzTot_->Fill(p4tot[2]);            
  h_GENP_ETot_->Fill(p4tot[3]);             
  h_GENP_DeltaEcmTotalSqrtS_->Fill(p4tot[3]-nEcms);
  h_GENP_DeltaEcmTotalSqrtS_Vs_NPart_->Fill(nPart,p4tot[3]-nEcms);
 
  Bool_t isThereAZ(false);

   const reco::Candidate* LeadingZLepton;
   const reco::Candidate* SubLeadingZLepton;

  if(WhichFlavour(theGenMuons,theGenElectrons).first || WhichFlavour(theGenMuons,theGenElectrons).second ){
    LeadingZLepton    = getTheZLeptons(theGenMuons,theGenElectrons).first;
    SubLeadingZLepton = getTheZLeptons(theGenMuons,theGenElectrons).second;
    isThereAZ = true;
  }  

  std::vector<const reco::Candidate*> theLastBHadrons; 
  std::vector<bool> islastBHadronsAssoc;
  std::map<const reco::Candidate*,bool> lastBHadrons_and_Assoc;

  //  std::cout<<"before loop on genParticles"<<std::endl;

  //Loop on genParticles
  for( GenParticleCollection::const_iterator genp = genParticlesCollection->begin();genp != genParticlesCollection->end(); ++ genp ) {  // loop over GEN particles
    
    const reco::Candidate* p = &(*genp);

    //------------------------------------------- b hadrons -------------------------------------------
    if( hasBottom(*p) ){ //if b-hadrons
      nBHad++;

      // cout<<"catched b-hadron"<<endl;
      
      h1bHadronStatus_->Fill(p->status());
      h1bHadronMotherId_->Fill(p->mother()->pdgId());

      h1BHadronsPt_->Fill(std::min(p->pt(),49.9),wMC);  	         
      h1BHadronsEta_->Fill(std::min(TMath::Abs(p->eta()),4.9),wMC);

      // looks if all daughter are bottomless
      Bool_t hasBottomedDaughter = false;
      for(UInt_t i=0; i<p->numberOfDaughters(); i++){
	if(hasBottom(*(p->daughter(i)))){	
	  hasBottomedDaughter=true;
	}
      } 
      
      if(!hasBottomedDaughter && TMath::Abs(p->eta())<3.5 ){ //no bottomed daughter and b-hadron in |eta| < 3.5

	theLastBHadrons.push_back(p);

	h1LastBHadronsPt_->Fill(std::min(p->pt(),49.9),wMC);  
	h1LastBHadronsEta_->Fill(std::min(TMath::Abs(p->eta()),4.9),wMC); 
	h1LastBHadronsMass_->Fill(p->mass(),wMC); 

	nLastBHad++;
	math::XYZTLorentzVectorD p4HFGEN(p->px(),p->py(),p->pz(),p->energy());

	// indices for b-hadron association
	std::pair<int,double>  GenGenAssociation = std::make_pair(-1,9999.);
	std::pair<const reco::GenJet*,double> GenJetAssociation;
	if(genJets.size()>0) GenJetAssociation = std::make_pair(&genJets.at(0),9999.);

	double minDeltaRGenGen(9999.);
	int i(0);

	for(std::vector<reco::GenJet>::const_iterator genjet_it=genJets.begin(); genjet_it!=genJets.end(); ++genjet_it){ 
	  // first check that the gen-jet doesn't overlap with a Z
	  if(isThereAZ && ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),LeadingZLepton->momentum())>0.5 && ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),SubLeadingZLepton->momentum())>0.5){
	    if(ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4HFGEN.Vect())<3.5) { 
	      if(ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4HFGEN.Vect())< minDeltaRGenGen){
		minDeltaRGenGen = ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4HFGEN.Vect());
		GenGenAssociation.first  = i;
		GenGenAssociation.second = minDeltaRGenGen;
	     
		GenJetAssociation.first  = &(*genjet_it);
		GenJetAssociation.second = minDeltaRGenGen;
		
	      } // end if minimum distance gen jet-parton
	      //   std::cout<<"i: "<<i<<" minDeltaR: "<<minDeltaRGenGen<<std::endl;
	    } // ends if b-hadron loosely matched to gen jet
	  } // ends if for Z lepton removal
	  i++;
	}// ends loop on gen jets
	
	// if there is a loose association DeltaR<3.5
	if( GenGenAssociation.first != -1 ){

	  // if there is a poor association 0.5 < DeltaR < 3.5
	  if(minDeltaRGenGen>0.5){
	    h1PoorlyAssociatedLastBHadronsPt_->Fill(std::min(p->pt(),49.99),wMC);  	         
	    h1PoorlyAssociatedLastBHadronsEta_->Fill(std::min(TMath::Abs(p->eta()),4.99),wMC); 
	    lastBHadrons_and_Assoc[p] = false;
	    islastBHadronsAssoc.push_back(false);
	  } else if(minDeltaRGenGen <= 0.5) {
	    // there is a tight association DeltaR < 0.5
	    h_GENP_pt_b_tightmatch_  ->Fill(std::min(GenJetAssociation.first->pt(),149.9),wMC);   
	    h_GENP_eta_b_tightmatch_ ->Fill(GenJetAssociation.first->eta()>0 ? std::min(GenJetAssociation.first->eta(),4.99) : std::max(GenJetAssociation.first->eta(),-4.99) ,wMC);  
	    h_GENP_mass_b_tightmatch_->Fill(GenJetAssociation.first->mass(),wMC);  
	    
	    if(GenJetAssociation.first->maxDistance()<0.5) h_GENP_mass_b_tightmatch_maxDistCut_->Fill(GenJetAssociation.first->mass(),wMC); 
                                                                            	    
	    h2_GENP_massb_vs_etab_->Fill(GenJetAssociation.first->eta()>0 ? std::min(GenJetAssociation.first->eta(),4.99) : std::max(GenJetAssociation.first->eta(),-4.99),GenJetAssociation.first->mass(),wMC);   
	    h2_GENP_massb_vs_ptb_->Fill(std::min(GenJetAssociation.first->pt(),149.9),GenJetAssociation.first->mass(),wMC);    	    
	    h2_GENP_massb_vs_maxDistanceb_->Fill(GenJetAssociation.first->maxDistance(),GenJetAssociation.first->mass(),wMC);

	    p_GENP_massb_vs_etab_->Fill(GenJetAssociation.first->eta()>0 ? std::min(GenJetAssociation.first->eta(),4.99) : std::max(GenJetAssociation.first->eta(),-4.99),GenJetAssociation.first->mass(),wMC);   
	    p_GENP_massb_vs_ptb_->Fill(std::min(GenJetAssociation.first->pt(),149.9),GenJetAssociation.first->mass(),wMC);     
	    p_GENP_massb_vs_maxDistanceb_->Fill(GenJetAssociation.first->maxDistance(),GenJetAssociation.first->mass(),wMC);
	
	    h1LastBHadronsPtTightAssoc_->Fill(std::min(p->pt(),49.9),wMC);  
	    h1LastBHadronsEtaTightAssoc_->Fill(std::min(TMath::Abs(p->eta()),4.9),wMC); 
	    h1LastBHadronsMassTightAssoc_->Fill(p->mass(),wMC); 
	    
	    h_GENP_maxDistance_b_tightmatch_->Fill(GenJetAssociation.first->maxDistance(),wMC);

	    std::vector <const GenParticle*> mcparticles = GenJetAssociation.first->getGenConstituents();
	    h1BgenJetNConst_->Fill(mcparticles.size(),wMC);
	    for( std::vector <const GenParticle*>::const_iterator thepart =mcparticles.begin();thepart != mcparticles.end(); ++ thepart ) { 
	      //std::cout<<"particle id: "<<(**thepart).pdgId()<<std::endl;
	      switch (abs((**thepart).pdgId())) {
	      case 11:
		h1BgenJetConstPdgId_->Fill(0.,wMC); //e
		break;
	      case 12:
		h1BgenJetConstPdgId_->Fill(1.,wMC); //nu_e
		break;
	      case 13:
		h1BgenJetConstPdgId_->Fill(2.,wMC); //mu
		break;
	      case 14:
		h1BgenJetConstPdgId_->Fill(3.,wMC); //nu_mu
		break;
	      case 16:
		h1BgenJetConstPdgId_->Fill(4.,wMC); //nu_tau
		break;
	      case 22:
		h1BgenJetConstPdgId_->Fill(5.,wMC); //gamma
		break;
	      case 211:
		h1BgenJetConstPdgId_->Fill(6.,wMC); //pions
		break;
	      case 130:
		h1BgenJetConstPdgId_->Fill(7.,wMC); //K0L
		break;
	      case 310:
		h1BgenJetConstPdgId_->Fill(8.,wMC); //K0S
		break;  
	      case 321:
		h1BgenJetConstPdgId_->Fill(9.,wMC); //K+/-
		break;
	      case 2212:
		h1BgenJetConstPdgId_->Fill(10.,wMC); //proton
		break;	
	      case 2112:
		h1BgenJetConstPdgId_->Fill(11.,wMC); //neutron
		break;
	      default:
		h1BgenJetConstPdgId_->Fill(12.,wMC); //something else
		break;
	      }
	    } 
	 
	    if(minDeltaRGenGen <= 0.1){
	      h_GENP_pt_b_extratight_->Fill(std::min(GenJetAssociation.first->pt(),149.9),wMC);   													      
	      h_GENP_eta_b_extratight_->Fill(GenJetAssociation.first->eta()>0 ? std::min(GenJetAssociation.first->eta(),4.99) : std::max(GenJetAssociation.first->eta(),-4.99) ,wMC);       
	      h_GENP_mass_b_extratight_->Fill(GenJetAssociation.first->mass(),wMC);
	      h_GENP_maxDistance_b_extratight_->Fill(GenJetAssociation.first->maxDistance(),wMC);
	    }

	    islastBHadronsAssoc.push_back(true);
	    lastBHadrons_and_Assoc[p] = true;
	  } // ends if tight association
	  
	  // for all b-jets
	  h_GENP_pt_b_  ->Fill(std::min(GenJetAssociation.first->pt(),149.9),wMC);   
	  h_GENP_eta_b_ ->Fill(GenJetAssociation.first->eta()>0 ? std::min(GenJetAssociation.first->eta(),4.99) : std::max(GenJetAssociation.first->eta(),-4.99) ,wMC);  
	  h_GENP_mass_b_->Fill(GenJetAssociation.first->mass(),wMC);                                                                                                                              
	  h_GENP_maxDistance_b_->Fill(GenJetAssociation.first->maxDistance(),wMC);

	  // DeltaR analysis
	  hDeltaRLastBHadGenJet_->Fill(std::min(minDeltaRGenGen,3.49),wMC);

	  // emulate reco-cuts
	  if(GenJetAssociation.first->pt()>25 && abs(GenJetAssociation.first->eta())<2.1){
	    hDeltaRLastBHadGoodGenJet_->Fill(std::min(minDeltaRGenGen,3.49),wMC);
	  }
 
	  // as a function of pt
	  if(GenJetAssociation.first->pt()<20.){
	    hDeltaRLastBHadGenJet_ptHat_020_->Fill(std::min(minDeltaRGenGen,3.49),wMC);     
	  } else if(GenJetAssociation.first->pt()>=20. && GenJetAssociation.first->pt()<25.){
	    hDeltaRLastBHadGenJet_ptHat_2025_->Fill(std::min(minDeltaRGenGen,3.49),wMC);
	  } else if(GenJetAssociation.first->pt()>=25. && GenJetAssociation.first->pt()<30.){
	    hDeltaRLastBHadGenJet_ptHat_2530_->Fill(std::min(minDeltaRGenGen,3.49),wMC);
	  } else {
	    hDeltaRLastBHadGenJet_ptHat_GT30_->Fill(std::min(minDeltaRGenGen,3.49),wMC);
	  } // ends switch on pthats
	} // if b-hadron is associated to genJet
	else {
	  islastBHadronsAssoc.push_back(false);
	  lastBHadrons_and_Assoc[p] = true;
	}
      } // if all daughter are bottomless

      while (p->mother()!=0 && hasBottom(*(p->mother()))) {
	p = p->mother();
      }

      //cout<<"last b-hadron mother: "<<p->mother()->pdgId()<<endl;      
      h1FirstbHadronId_->Fill(p->pdgId());
      h1FirstbHadronMotherId_->Fill(p->mother()->pdgId());
    } // if particle has bottom
    
    //------------------------------------------- Z boson -------------------------------------------
    if(p->pdgId()==23){ //&& p->status()==2){
      math::XYZTLorentzVectorD p4ZGENtemp(p->px(),p->py(),p->pz(),p->energy());    
      p4ZGEN=p4ZGENtemp;
    }
  }//loop on genParticles
  
  h1nBHadrons_->Fill(nBHad,wMC);
  h1nLastBHadrons_->Fill(nLastBHad,wMC);

  // a little bit of b-hadrons
  for(std::vector<reco::GenJet>::const_iterator genjet_it=genJets.begin(); genjet_it!=genJets.end(); ++genjet_it){
    Int_t nBhadsInJet(0);
    for(std::vector<const reco::Candidate*>::const_iterator thebhadron_it = theLastBHadrons.begin();  thebhadron_it!= theLastBHadrons.end(); ++thebhadron_it){
      // first check that the gen-jet doesn't overlap with a Z
      if(isThereAZ && ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),LeadingZLepton->momentum())>0.5 && ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),SubLeadingZLepton->momentum())>0.5){
	if(ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),(*thebhadron_it)->momentum())<0.5 ) { 
	  nBhadsInJet++;
	} // ends if parton matched to gen jet
      }// if genJet doesn't overlap with Z leptons
    }// ends loop on gen jets
    h1nLastBHadronsPerJet_->Fill(nBhadsInJet);
  }// ends loop on "last" b-hadrons 

  // if(nLastBHad>2){
  //  std::cout<<"nBhads: "<<nLastBHad<<" Run: "<<iEvent.eventAuxiliary().run()<<"\t  Lumi:"<<iEvent.eventAuxiliary().luminosityBlock()<<"\t Event: "<< iEvent.eventAuxiliary().id().event()<<std::endl;
  //}
  
  //Double loop on b-hadrons
  //   for(std::vector<const reco::Candidate*>::const_iterator thebhadron_it1 = theLastBHadrons.begin();  thebhadron_it1!= theLastBHadrons.end(); ++thebhadron_it1){
  //     for(std::vector<const reco::Candidate*>::const_iterator thebhadron_it2 = thebhadron_it1+1;  thebhadron_it2!= theLastBHadrons.end(); ++thebhadron_it1){
  //       h1LastBhadronDeltaR_->Fill(ROOT::Math::VectorUtil::DeltaR((*thebhadron_it1)->momentum(),(*thebhadron_it2)->momentum()));
  //     }
  //   }

  // Double loop on b-hadrons
  std::map<const reco::Candidate*,bool>::const_iterator thebhadron_it1,thebhadron_it2;  
  for(thebhadron_it1 = lastBHadrons_and_Assoc.begin(); thebhadron_it1 != lastBHadrons_and_Assoc.end(); thebhadron_it1++){
    for(thebhadron_it2 = thebhadron_it1; thebhadron_it2 != lastBHadrons_and_Assoc.end(); thebhadron_it2++){
      if(thebhadron_it1!=thebhadron_it2){
  	h1LastBhadronDeltaR_->Fill(ROOT::Math::VectorUtil::DeltaR(((*thebhadron_it1).first)->momentum(),((*thebhadron_it2).first)->momentum()));
  	if((*thebhadron_it1).second && (*thebhadron_it2).second) {
  	  h1LastBhadronGenJetMatchedDeltaR_->Fill(ROOT::Math::VectorUtil::DeltaR(((*thebhadron_it1).first)->momentum(),((*thebhadron_it2).first)->momentum()));
  	}
      }
    }
  }
  
  // std::cout<<"theLastBHadrons.size()"<<theLastBHadrons.size()<<"islastBHadronsAssoc"<<islastBHadronsAssoc.size()<<std::endl;

  //Fill electron counter
  h1NEle_->Fill(theElectronCounter_,wMC);
  //Fill muon counter
  h1NMu_->Fill(theMuonCounter_,wMC);

  h_GENP_pt_Z_->Fill(p4ZGEN.pt(),wMC);     
  h_GENP_eta_Z_->Fill(p4ZGEN.eta(),wMC); 
  h_GENP_mass_Z_->Fill(p4ZGEN.mass(),wMC);
  
  std::vector<math::XYZTLorentzVectorD> theGenElectronsMomenta;
  std::vector<math::XYZTLorentzVectorD> theGenMuonsMomenta;

  //std::cout<<"genMuons.size(): "<<theGenMuons.size()<<" genElectrons.size(): "<<theGenElectrons.size()<<std::endl;
  
  if(WhichFlavour(theGenMuons,theGenElectrons).second){
  
    for (size_t i=0;i<theGenElectrons.size();i++){
      math::XYZTLorentzVectorD p4Lepttemp(theGenElectrons[i]->px(),theGenElectrons[i]->py(),theGenElectrons[i]->pz(),theGenElectrons[i]->energy());  
      theGenElectronsMomenta.push_back(p4Lepttemp);
    }
    
    std::vector<math::XYZTLorentzVectorD> sortedElectronsMomenta = sort4MomentaByPt(theGenElectronsMomenta);
    std::vector<const reco::Candidate*> sortedGenElectrons = sortLeptonsByPt(theGenElectrons);

    for (size_t i=0;i<2;i++){
      math::XYZTLorentzVectorD p4ZGENFromLepttemp = sortedElectronsMomenta[i];
      p4ZGENFromLept += p4ZGENFromLepttemp;
    }
    
    //for (size_t i=0;i<sortedElectronsMomenta.size();i++){
    //std::cout<<" sortedElectronsMomenta["<<i<<"].pt()= "<<sortedElectronsMomenta[i].pt()<<std::endl;;
    //}

    h_GENP_zllCosThetaStar_->Fill(CosThetaStar(sortedGenElectrons)); 
    h_GENP_zllCosThetaCS_->Fill(CostCS(sortedGenElectrons));        

    h_GENP_pt_Lept1_fromZ_->Fill(sortedElectronsMomenta[0].pt(),wMC);         
    h_GENP_eta_Lept1_fromZ_->Fill(sortedElectronsMomenta[0].eta(),wMC);        
        				
    h_GENP_pt_Electrons1_fromZ_->Fill(sortedElectronsMomenta[0].pt(),wMC);     
    h_GENP_eta_Electrons1_fromZ_->Fill(sortedElectronsMomenta[0].eta(),wMC);  
    
    h_GENP_pt_Lept2_fromZ_->Fill(sortedElectronsMomenta[1].pt(),wMC);          
    h_GENP_eta_Lept2_fromZ_->Fill(sortedElectronsMomenta[1].eta(),wMC);       
     
    h_GENP_pt_Electrons2_fromZ_->Fill(sortedElectronsMomenta[1].pt(),wMC);     
    h_GENP_eta_Electrons2_fromZ_->Fill(sortedElectronsMomenta[1].eta(),wMC);

    h_GENP_pt_ZfromLept_->Fill(p4ZGENFromLept.pt(),wMC);  
    h_GENP_eta_ZfromLept_->Fill(p4ZGENFromLept.eta(),wMC); 
    h_GENP_mass_ZfromLept_->Fill(p4ZGENFromLept.mass(),wMC);
    
    h_GENP_pt_ZfromElectrons_->Fill(p4ZGENFromLept.pt(),wMC);  
    h_GENP_eta_ZfromElectrons_->Fill(p4ZGENFromLept.eta(),wMC); 
    h_GENP_mass_ZfromElectrons_->Fill(p4ZGENFromLept.mass(),wMC);

    if(p4ZGENFromLept.mass()>50){
      
      h_GENP_mll50_zllCosThetaStar_->Fill(CosThetaStar(sortedGenElectrons)); 
      h_GENP_mll50_zllCosThetaCS_->Fill(CostCS(sortedGenElectrons));        
      
      h_GENP_mll50_pt_Lept1_fromZ_->Fill(sortedElectronsMomenta[0].pt(),wMC);         
      h_GENP_mll50_eta_Lept1_fromZ_->Fill(sortedElectronsMomenta[0].eta(),wMC);        
      
      h_GENP_mll50_pt_Electrons1_fromZ_->Fill(sortedElectronsMomenta[0].pt(),wMC);     
      h_GENP_mll50_eta_Electrons1_fromZ_->Fill(sortedElectronsMomenta[0].eta(),wMC);  
      
      h_GENP_mll50_pt_Lept2_fromZ_->Fill(sortedElectronsMomenta[1].pt(),wMC);          
      h_GENP_mll50_eta_Lept2_fromZ_->Fill(sortedElectronsMomenta[1].eta(),wMC);       
      
      h_GENP_mll50_pt_Electrons2_fromZ_->Fill(sortedElectronsMomenta[1].pt(),wMC);     
      h_GENP_mll50_eta_Electrons2_fromZ_->Fill(sortedElectronsMomenta[1].eta(),wMC);
      
      h_GENP_mll50_pt_ZfromLept_->Fill(p4ZGENFromLept.pt(),wMC);  
      h_GENP_mll50_eta_ZfromLept_->Fill(p4ZGENFromLept.eta(),wMC); 
      h_GENP_mll50_mass_ZfromLept_->Fill(p4ZGENFromLept.mass(),wMC);
      
      h_GENP_mll50_pt_ZfromElectrons_->Fill(p4ZGENFromLept.pt(),wMC);  
      h_GENP_mll50_eta_ZfromElectrons_->Fill(p4ZGENFromLept.eta(),wMC); 
      h_GENP_mll50_mass_ZfromElectrons_->Fill(p4ZGENFromLept.mass(),wMC);
    } 
  } else if(WhichFlavour(theGenMuons,theGenElectrons).first){
    
    for (size_t i=0;i<theGenMuons.size();i++){
      math::XYZTLorentzVectorD p4Lepttemp(theGenMuons[i]->px(),theGenMuons[i]->py(),theGenMuons[i]->pz(),theGenMuons[i]->energy());  
      theGenMuonsMomenta.push_back(p4Lepttemp);
    }

    std::vector<math::XYZTLorentzVectorD> sortedMuonsMomenta = sort4MomentaByPt(theGenMuonsMomenta);
    std::vector<const reco::Candidate*>   sortedGenMuons     = sortLeptonsByPt(theGenMuons);

    for (size_t i=0;i<2;i++){
      math::XYZTLorentzVectorD p4ZGENFromLepttemp = sortedMuonsMomenta[i];
      p4ZGENFromLept += p4ZGENFromLepttemp;
    }
    
    h_GENP_zllCosThetaStar_->Fill(CosThetaStar(sortedGenMuons)); 
    h_GENP_zllCosThetaCS_->Fill(CostCS(sortedGenMuons));  

    h_GENP_pt_Lept1_fromZ_->Fill(sortedMuonsMomenta[0].pt(),wMC);         
    h_GENP_eta_Lept1_fromZ_->Fill(sortedMuonsMomenta[0].eta(),wMC);        
        				
    h_GENP_pt_Muons1_fromZ_->Fill(sortedMuonsMomenta[0].pt(),wMC);     
    h_GENP_eta_Muons1_fromZ_->Fill(sortedMuonsMomenta[0].eta(),wMC);  
    
    h_GENP_pt_Lept2_fromZ_->Fill(sortedMuonsMomenta[1].pt(),wMC);          
    h_GENP_eta_Lept2_fromZ_->Fill(sortedMuonsMomenta[1].eta(),wMC);       
     
    h_GENP_pt_Muons2_fromZ_->Fill(sortedMuonsMomenta[1].pt(),wMC);     
    h_GENP_eta_Muons2_fromZ_->Fill(sortedMuonsMomenta[1].eta(),wMC);

    h_GENP_pt_ZfromLept_->Fill(p4ZGENFromLept.pt(),wMC);  
    h_GENP_eta_ZfromLept_->Fill(p4ZGENFromLept.eta(),wMC); 
    h_GENP_mass_ZfromLept_->Fill(p4ZGENFromLept.mass(),wMC);
    
    h_GENP_pt_ZfromMuons_->Fill(p4ZGENFromLept.pt(),wMC);  
    h_GENP_eta_ZfromMuons_->Fill(p4ZGENFromLept.eta(),wMC); 
    h_GENP_mass_ZfromMuons_->Fill(p4ZGENFromLept.mass(),wMC);

    if(p4ZGENFromLept.mass()>50.){
      
      h_GENP_mll50_zllCosThetaStar_->Fill(CosThetaStar(sortedGenMuons)); 
      h_GENP_mll50_zllCosThetaCS_->Fill(CostCS(sortedGenMuons));  
      
      h_GENP_mll50_pt_Lept1_fromZ_->Fill(sortedMuonsMomenta[0].pt(),wMC);         
      h_GENP_mll50_eta_Lept1_fromZ_->Fill(sortedMuonsMomenta[0].eta(),wMC);        
      
      h_GENP_mll50_pt_Muons1_fromZ_->Fill(sortedMuonsMomenta[0].pt(),wMC);     
      h_GENP_mll50_eta_Muons1_fromZ_->Fill(sortedMuonsMomenta[0].eta(),wMC);  
      
      h_GENP_mll50_pt_Lept2_fromZ_->Fill(sortedMuonsMomenta[1].pt(),wMC);          
      h_GENP_mll50_eta_Lept2_fromZ_->Fill(sortedMuonsMomenta[1].eta(),wMC);       
      
      h_GENP_mll50_pt_Muons2_fromZ_->Fill(sortedMuonsMomenta[1].pt(),wMC);     
      h_GENP_mll50_eta_Muons2_fromZ_->Fill(sortedMuonsMomenta[1].eta(),wMC);
      
      h_GENP_mll50_pt_ZfromLept_->Fill(p4ZGENFromLept.pt(),wMC);  
      h_GENP_mll50_eta_ZfromLept_->Fill(p4ZGENFromLept.eta(),wMC); 
      h_GENP_mll50_mass_ZfromLept_->Fill(p4ZGENFromLept.mass(),wMC);
      
      h_GENP_mll50_pt_ZfromMuons_->Fill(p4ZGENFromLept.pt(),wMC);  
      h_GENP_mll50_eta_ZfromMuons_->Fill(p4ZGENFromLept.eta(),wMC); 
      h_GENP_mll50_mass_ZfromMuons_->Fill(p4ZGENFromLept.mass(),wMC);
    }
  }

  //std::cout<<"before loop on genJets"<<std::endl;

  //loop on genJets
  Int_t ngoodGenJets(0);
  for(std::vector<reco::GenJet>::const_iterator genjet_it=genJets.begin(); genjet_it!=genJets.end(); ++genjet_it){ 
 
    if(isThereAZ && ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),LeadingZLepton->momentum())>0.5 && ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),SubLeadingZLepton->momentum())>0.5){
    
      h_GENP_pt_all_  ->Fill(std::min(genjet_it->pt(),149.9),wMC);   													       
      h_GENP_eta_all_ ->Fill(genjet_it->eta()>0 ? std::min(genjet_it->eta(),4.99) : std::max(genjet_it->eta(),-4.99) ,wMC);       
      h_GENP_mass_all_->Fill(genjet_it->mass(),wMC);    
      h_GENP_maxDistance_all_->Fill(genjet_it->maxDistance(),wMC);
      
      if(genjet_it->pt()>25. && fabs(genjet_it->eta())<2.1) ngoodGenJets++;
      
      std::vector <const GenParticle*> mcparts = genjet_it->getGenConstituents();
      h1genJetNConst_->Fill(mcparts.size(),wMC);

      Double_t maxDeltaRJetPart=0.;

      for( std::vector <const GenParticle*>::const_iterator thePart =mcparts.begin();thePart != mcparts.end(); ++ thePart ) { 
	//std::cout<<"particle id: "<<(**thePart).pdgId()<<std::endl;
	if( ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),(**thePart).momentum()) > maxDeltaRJetPart) maxDeltaRJetPart= ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),(**thePart).momentum());

	switch (abs((**thePart).pdgId())) {
	case 11:
	  h1genJetConstPdgId_->Fill(0.,wMC); //e
	  break;
	case 12:
	  h1genJetConstPdgId_->Fill(1.,wMC); //nu_e
	  break;
	case 13:
	  h1genJetConstPdgId_->Fill(2.,wMC); //mu
	  break;
	case 14:
	  h1genJetConstPdgId_->Fill(3.,wMC); //nu_mu
	  break;
	case 16:
	  h1genJetConstPdgId_->Fill(4.,wMC); //nu_tau
	  break;
	case 22:
	  h1genJetConstPdgId_->Fill(5.,wMC); //gamma
	  break;
	case 211:
	  h1genJetConstPdgId_->Fill(6.,wMC); //pions
	  break;
	case 130:
	  h1genJetConstPdgId_->Fill(7.,wMC); //K0L
	  break;
	case 310:
	  h1genJetConstPdgId_->Fill(8.,wMC); //K0S
	  break;  
	case 321:
	  h1genJetConstPdgId_->Fill(9.,wMC); //K+/-
	  break;
	case 2212:
	  h1genJetConstPdgId_->Fill(10.,wMC); //proton
	  break;	
	case 2112:
	  h1genJetConstPdgId_->Fill(11.,wMC); //neutron
	  break;
	default:
	  h1genJetConstPdgId_->Fill(12.,wMC); //something else
	  break;
	}
      }
      h_GENP_MymaxDistance_all_->Fill(maxDeltaRJetPart,wMC);
    }
    //std::cout<<"<<<<<<<<<<<<< end of jet <<<<<<<<<<<<<"<<std::endl;
  }

  //std::cout<<"after loop on genJets"<<std::endl;

  h1NGoodGenJets_->Fill(ngoodGenJets,wMC);
}

// ------------ method called once each job just before starting event loop  ------------
void 
MCTruthAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;
  
  TH1F::SetDefaultSumw2(kTRUE);

  // Counter histograms
  h1NMu_              = fs->make<TH1F>("NMu","Number of muons per event; number of muons N_{#mu}; events",20,-0.5,19.5);
  h1NEle_             = fs->make<TH1F>("NEle","Number of electrons per event; number of electrons N_{e}; events",20,-0.5,19.5);

  h1NRecoJets_        = fs->make<TH1F>("NRecoJets","Number of reco jets per event; N_{jets}; events",15,-0.5,14.5);
  h1NRecoJetsVsWeight_= fs->make<TH2F>("NRecoJetsVsWeight","Number of reco jets per event;Number of reco jets per event; MC Weight",15,-0.5,14.5,500.,0.,1.);
  h1NRecoJetsW_       = fs->make<TH1F>("NRecoJetsW","Number of reco jets per event (weighted); N_{jets}; weighted events",15,-0.5,14.5);

  h1NGenJets_         = fs->make<TH1F>("NGenJets","Number of gen jets per event; N_{jets}^{gen}; events",40,-0.5,39.5);
  h1NGenJetsVsWeight_ = fs->make<TH2F>("NGenJetsVsWeight","Number of gen jets per event; Number of gen jets per event; MC Weight",40,-0.5,39.5,500.,0.,1.);
  h1NGenJetsW_        = fs->make<TH1F>("NGenJetsW","Number of gen jets per event (weighted); N_{jets}^{gen}; events",40,-0.5,39.5);

  hMCWeight_          = fs->make<TH1F>("MCWeight","MC Weight; MC weight; events",500,0,50);
  hMCWeightZoom_      = fs->make<TH1F>("MCWeightZoom","MC Weight; MC weight; events",100,-2,2);

  // PDGId of lepton's mother
  h1MotherIdMuon_     = fs->make<TH1F>("MotherIdMuon","PDGId of muon's mother",500,-0.5,499.5);
  h1MotherIdElectron_ = fs->make<TH1F>("MotherIdElectron","PDGId of electron's mother",500,-0.5,499.5);

  // Spectra histograms
  h1PtSoftMuon_        = fs->make<TH1F>("PtSoftMuon","Pt of the less energetic muon; p_{T}(soft #mu) [GeV]; muons",100,0,100);
  h1PtLeadingMuon_     = fs->make<TH1F>("PtLeadingMuon","Pt of more energetic muon; p_{T}(hard #mu) [GeV]; muons",100,0,150);
  h1PtSoftElectron_    = fs->make<TH1F>("PtSoftElectron","Pt of the less energetic electron; p_{T}(soft e) [GeV]; electrons",100,0,100);
  h1PtLeadingElectron_ = fs->make<TH1F>("PtLeadingElectron","Pt of more energetic electron; p_{T}(hard e) [GeV]; electrons",100,0,150);
  h1DrMinMJ_           = fs->make<TH1F>("DrMinMJ","DeltaR min not-leading muon/jet; #DeltaR(#mu,jet); muons",10,0,10);
  h1PtMuFromHF_        = fs->make<TH1F>("PtMuFromHF","Pt muons from heavy flavour; p_{T}(HF-#mu) [GeV]; muons",50,0.,50.);
  h1PtSoftJet_         = fs->make<TH1F>("PtSoftJet","Pt of the less energetic jet; p_{T}(soft jet) [GeV]; jets",50,0,100);
  h1PtAllJets_         = fs->make<TH1F>("PtAllJet","Pt of all jets; p_{T}(all jets) [GeV]; jets",100,0,200);
  h1PtLeadingJet_      = fs->make<TH1F>("PtLeadingJet","Pt of more energetic jet; p_{T}(hard jet) [GeV]; jets",100,0,200);
  
  // weighted spectra histograms
  h1PtSoftMuonW_        = fs->make<TH1F>("PtSoftMuonW","Pt of the less energetic muon (weighted); p_{T}(soft #mu) [GeV];  weighted muons",100,0,100);
  h1PtLeadingMuonW_     = fs->make<TH1F>("PtLeadingMuonW","Pt of more energetic muon (weighted); p_{T}(hard #mu) [GeV]; weighted muons",100,0,150);
  h1PtSoftElectronW_    = fs->make<TH1F>("PtSoftElectronW","Pt of the less energetic electron (weighted); p_{T}(soft e) [GeV];  weighted electrons",100,0,100);
  h1PtLeadingElectronW_ = fs->make<TH1F>("PtLeadingElectronW","Pt of more energetic electron (weighted); p_{T}(hard e) [GeV]; weighted electrons",100,0,150);
  h1DrMinMJW_           = fs->make<TH1F>("DrMinMJW","DeltaR min not-leading muon/jet(weighted); #DeltaR(#mu,jet); weighted muons",10,0,10);
  h1PtMuFromHFW_        = fs->make<TH1F>("PtMuFromHFW","Pt muons from heavy flavour (weighted); p_{T}(HF-#mu) [GeV]; weighted muons",50,0.,50.);
  h1PtSoftJetW_         = fs->make<TH1F>("PtSoftJetW","Pt of the less energetic jet (weighted); p_{T}(soft jet) [GeV]; weighted jets",50,0,100);
  h1PtAllJetsW_         = fs->make<TH1F>("PtAllJetW","Pt of all jets; p_{T}(all jets) [GeV];weighted jets",100,0,200);
  h1PtLeadingJetW_      = fs->make<TH1F>("PtLeadingJetW","Pt of more energetic jet(weighted); p_{T}(hard jet) [GeV];weighted jets",100,0,200);
    
  //*************************************************************************************************
  // MC truth for b-hadrons
  //*************************************************************************************************

  h1bHadronStatus_         =  fs->make<TH1F>("bHadronStatus","particle status of b-hadron",5,-0.5,4.5);

  Double_t pdgIdBins[6] = {0.,10.,100.,1000.,10000.,100000.};
  h1bHadronMotherId_       = fs->make<TH1F>("bHadronMotherId","Pdg Id of b-hadron mother; pdg ID of b-hadron mother",5,pdgIdBins);
  h1FirstbHadronId_        = fs->make<TH1F>("FirstbHadronId","Pdg Id of first b-hadron in decay chain; Pdg Id of first b-hadron in decay chain",5,pdgIdBins);
  h1FirstbHadronMotherId_  = fs->make<TH1F>("FirstbHadronMotherId","Pdg Id of mother of first b-hadron in decay chain; Pdg Id of mother of first b-hadron in decay chain",10,-0.5,9.5);

  h1nBHadrons_             = fs->make<TH1F>("nBHadrons","Number of all b-hadrons; N_{B-Had}; b-hadrons",11,-0.5,10.5);
  h1BHadronsPt_            = fs->make<TH1F>("BHadronsPt","p_{T} of all b-hadron; p_{T}^{B-Had} [GeV]; b-hadrons",50,0.,50.);
  h1BHadronsEta_           = fs->make<TH1F>("BHadronsEta","#eta of all b-hadron; #eta^{B-Had}; b-hadrons",50,0.,5.);

  h1nLastBHadrons_         = fs->make<TH1F>("nLastBHadrons","Number of \"last \" b-hadrons; N_{B-Had}; b-hadrons",11,-0.5,10.5);
  h1nLastBHadronsPerJet_   = fs->make<TH1F>("nLastBHadronsPerJet","Number of \"last \" b-hadrons; N_{B-Had} per jet; gen Jets",11,-0.5,10.5);
    
  h1LastBHadronsPt_        = fs->make<TH1F>("LastBHadronsPt","p_{T} of \"last\" b-hadron; p_{T}^{B-Had} (GeV); b-hadrons",50,0.,50.);
  h1LastBHadronsEta_       = fs->make<TH1F>("LastBHadronsEta","Eta of \"last\" b-hadron; #eta^{B-Had}; b-hadrons",50,0.,5.); 
  h1LastBHadronsMass_      = fs->make<TH1F>("LastBHadronsMass","mass of \"last\" b-hadron; m_{B-Had} (GeV)",100,0.,10.);         

  h1LastBHadronsPtTightAssoc_   = fs->make<TH1F>("LastBHadronsPtTightAssoc","p_{T} of \"last\" b-hadron; p_{T}^{B-Had} (GeV); b-hadrons",50,0.,50.);
  h1LastBHadronsEtaTightAssoc_  = fs->make<TH1F>("LastBHadronsEtaTightAssoc","Eta of \"last\" b-hadron; #eta^{B-Had}; b-hadrons",50,0.,5.);   
  h1LastBHadronsMassTightAssoc_ = fs->make<TH1F>("LastBHadronsMassTightAssoc","mass of \"last\" b-hadron; m_{B-Had} (GeV)",100,0.,10.);        
  
  h1PoorlyAssociatedLastBHadronsPt_ = fs->make<TH1F>("PoorlyAssociatedLastBHadronsPt","p_{T} of poorly associated b-hadron; p_{T}^{B-Had} [GeV]; b-hadrons",50,0.,50.);
  h1PoorlyAssociatedLastBHadronsEta_= fs->make<TH1F>("PoorlyAssociatedLastBHadronsEta","#eta of poorly associated b-hadron; #eta^{B-Had}; b-hadrons",50,0.,5.);

  hDeltaRLastBHadGoodGenJet_        = fs->make<TH1F>("DeltaRLastBHadGoodGenJet","#Delta R (bHad,genJet) - p_{T}>25, |#eta|<2.1; #Delta R (bHad,genJet); b-hadrons",50,0.,3.5);
  hDeltaRLastBHadGenJet_            = fs->make<TH1F>("DeltaRLastBHadGenJet","#Delta R (bHad,genJet); #Delta R (bHad,genJet); b-hadrons",50,0.,3.5);
  hDeltaRLastBHadGenJet_ptHat_020_  = fs->make<TH1F>("DeltaRLastBHadGenJet_ptHat_020", "#Delta R (bHad,genJet) #hat{p}_{T}(genJet)<20 GeV; #Delta R (bHad,genJet); b-hadrons",50,0.,3.5);
  hDeltaRLastBHadGenJet_ptHat_2025_ = fs->make<TH1F>("DeltaRLastBHadGenJet_ptHat_2025","#Delta R (bHad,genJet) 20< #hat{p}_{T}(genJet)<25 GeV; #Delta R (bHad,genJet); b-hadrons",50,0.,3.5);
  hDeltaRLastBHadGenJet_ptHat_2530_ = fs->make<TH1F>("DeltaRLastBHadGenJet_ptHat_2530","#Delta R (bHad,genJet) 25< #hat{p}_{T}(genJet)<30 GeV; #Delta R (bHad,genJet); b-hadrons",50,0.,3.5);
  hDeltaRLastBHadGenJet_ptHat_GT30_ = fs->make<TH1F>("DeltaRLastBHadGenJet_ptHat_GT30","#Delta R (bHad,genJet) #hat{p}_{T}(genJet)>30 GeV; #Delta R (bHad,genJet); b-hadrons",50,0.,3.5);
  
  // histograms for b-hadron angular distance
  h1LastBhadronDeltaR_               = fs->make<TH1F>("LastBHadronDeltaR","Angular distance of b-hadrons in the event; #DeltaR(bHad,bHad)",100,0.,10.);
  h1LastBhadronGenJetMatchedDeltaR_  = fs->make<TH1F>("LastBHadronGenJetMatchedDeltaR","Angular distance of b-hadrons (GenJetMatched) in the event; #DeltaR(bHad,bHad)",100,0.,10.);

  //*************************************************************************************************
  // MC Truth for Z and leptons 
  //*************************************************************************************************

  Int_t nbin_pt=50;
  Double_t pt_min= 0;
  Double_t pt_max=150.;
  Double_t pt_Zmax=250.;
  Int_t nbin_eta=50;
  Double_t eta_min=-6;
  Double_t eta_max =6;

  // distributions without cut on M(ll)

  h_GENP_pt_Z_                  = fs->make<TH1F>("GENP_pt_Z","p_{T} of Z boson (genParticles); GEN p_{T} of the Z (GeV)",nbin_pt,pt_min,pt_Zmax);  
  h_GENP_eta_Z_                 = fs->make<TH1F>("GENP_eta_Z","#eta of Z boson (genParticles); GEN #eta of the Z",nbin_eta,eta_min,eta_max); 
  h_GENP_mass_Z_                = fs->make<TH1F>("GENP_mass_Z","mass of Z boson (genParticles); GEN M_{Z} (GeV)",100.,0.,150.); 

  h_GENP_pt_ZfromLept_          = fs->make<TH1F>("GENP_pt_ZfromLept","p_{T} of ll (genParticles); GEN p_{T}(l^{+}l^{-}) (GeV)",nbin_pt,pt_min,pt_Zmax);  
  h_GENP_eta_ZfromLept_         = fs->make<TH1F>("GENP_eta_ZfromLept","#eta of ll (genParticles); GEN #eta(l^{+}l^{-})",nbin_eta,eta_min,eta_max); 
  h_GENP_mass_ZfromLept_        = fs->make<TH1F>("GENP_mass_ZfromLept","mass of ll (genParticles); GEN M_{l^{+}l^{-}} (GeV)",100.,0.,150.); 

  h_GENP_pt_ZfromMuons_         = fs->make<TH1F>("GENP_pt_ZfromMuons","p_{T} of #mu#mu (genParticles); GEN p_{T}(#mu^{+}#mu^{-}) (GeV)",nbin_pt,pt_min,pt_Zmax);  
  h_GENP_eta_ZfromMuons_        = fs->make<TH1F>("GENP_eta_ZfromMuons","#eta of #mu#mu (genParticles); GEN #eta(#mu^{+}#mu^{-})",nbin_eta,eta_min,eta_max); 
  h_GENP_mass_ZfromMuons_       = fs->make<TH1F>("GENP_mass_ZfromMuons","mass of #mu#mu (genParticles); GEN M_{#mu^{+}#mu^{-}} (GeV)",100.,0.,150.);
  
  h_GENP_pt_ZfromElectrons_     = fs->make<TH1F>("GENP_pt_ZfromElectrons","p_{T} of #it{ee} (genParticles); GEN p_{T}(e ^{+}e ^{-}) (GeV)",nbin_pt,pt_min,pt_Zmax);  
  h_GENP_eta_ZfromElectrons_    = fs->make<TH1F>("GENP_eta_ZfromElectrons","#eta of #it{ee} (genParticles); GEN #eta(e^{+}e^{-})",nbin_eta,eta_min,eta_max); 
  h_GENP_mass_ZfromElectrons_   = fs->make<TH1F>("GENP_mass_ZfromElectrons","mass of #it{ee} (genParticles); GEN M_{e^{+}e^{-}} (GeV)",100.,0.,150.);

  // leptons distributions
  h_GENP_pt_Lept1_fromZ_        = fs->make<TH1F>("GENP_pt_Lept1_fromZ","p_{T} of l (genParticles); GEN p_{T}(l^{lead}) (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_eta_Lept1_fromZ_       = fs->make<TH1F>("GENP_eta_Lept1_fromZ","#eta of l (genParticles); GEN #eta(l^{lead})",nbin_eta,eta_min,eta_max); 

  h_GENP_pt_Muons1_fromZ_       = fs->make<TH1F>("GENP_pt_Muons1_fromZ","p_{T} of #mu (genParticles); GEN p_{T}(#mu^{lead}) (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_eta_Muons1_fromZ_      = fs->make<TH1F>("GENP_eta_Muons1_fromZ","#eta of #mu (genParticles); GEN #eta(#mu^{lead})",nbin_eta,eta_min,eta_max); 

  h_GENP_pt_Electrons1_fromZ_   = fs->make<TH1F>("GENP_pt_Electrons1_fromZ","p_{T} of #it{e} (genParticles); GEN p_{T}(e ^{lead}) (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_eta_Electrons1_fromZ_  = fs->make<TH1F>("GENP_eta_Electrons1_fromZ","#eta of #it{e} (genParticles); GEN #eta(e^{lead})",nbin_eta,eta_min,eta_max); 

  h_GENP_pt_Lept2_fromZ_        = fs->make<TH1F>("GENP_pt_Lept2_fromZ","p_{T} of l (genParticles); GEN p_{T}(l^{2nd}) (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_eta_Lept2_fromZ_       = fs->make<TH1F>("GENP_eta_Lept2_fromZ","#eta of l (genParticles); GEN #eta(l^{2nd})",nbin_eta,eta_min,eta_max); 

  h_GENP_pt_Muons2_fromZ_       = fs->make<TH1F>("GENP_pt_Muons2_fromZ","p_{T} of #mu (genParticles); GEN p_{T}(#mu^{2nd}) (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_eta_Muons2_fromZ_      = fs->make<TH1F>("GENP_eta_Muons2_fromZ","#eta of #mu (genParticles); GEN #eta(#mu^{2nd})",nbin_eta,eta_min,eta_max); 

  h_GENP_pt_Electrons2_fromZ_   = fs->make<TH1F>("GENP_pt_Electrons2_fromZ","p_{T} of #it{e} (genParticles); GEN p_{T}(e ^{2nd}) (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_eta_Electrons2_fromZ_  = fs->make<TH1F>("GENP_eta_Electrons2_fromZ","#eta of #it{e} (genParticles); GEN #eta(e^{2nd})",nbin_eta,eta_min,eta_max); 

  h_GENP_zllCosThetaCS_         = fs->make<TH1F>("GENP_llCosThetaCS","Cosine of Collins-Soper frame angle; cos(#theta_{CS,l^{+}})",20,-1.,1.); 
  h_GENP_zllCosThetaStar_       = fs->make<TH1F>("GENP_llCosThetaStar","Cosine of Z decay frame l^{+} angle; cos(#theta^{*}_{l^{+}})",20,-1.,1.); 

  // same distributions with a cut on M(ll)>50 (to make fair comparison aMC@NLO/Sherpa/MadGraph)

  h_GENP_mll50_pt_ZfromLept_          = fs->make<TH1F>("GENP_mll50_pt_ZfromLept","p_{T} of ll (genParticles); GEN p_{T}(l^{+}l^{-}) (GeV)",nbin_pt,pt_min,pt_Zmax);  
  h_GENP_mll50_eta_ZfromLept_         = fs->make<TH1F>("GENP_mll50_eta_ZfromLept","#eta of ll (genParticles); GEN #eta(l^{+}l^{-})",nbin_eta,eta_min,eta_max); 
  h_GENP_mll50_mass_ZfromLept_        = fs->make<TH1F>("GENP_mll50_mass_ZfromLept","mass of ll (genParticles); GEN M_{l^{+}l^{-}} (GeV)",60.,40.,160.); 

  h_GENP_mll50_pt_ZfromMuons_         = fs->make<TH1F>("GENP_mll50_pt_ZfromMuons","p_{T} of #mu#mu (genParticles); GEN p_{T}(#mu^{+}#mu^{-}) (GeV)",nbin_pt,pt_min,pt_Zmax);  
  h_GENP_mll50_eta_ZfromMuons_        = fs->make<TH1F>("GENP_mll50_eta_ZfromMuons","#eta of #mu#mu (genParticles); GEN #eta(#mu^{+}#mu^{-})",nbin_eta,eta_min,eta_max); 
  h_GENP_mll50_mass_ZfromMuons_       = fs->make<TH1F>("GENP_mll50_mass_ZfromMuons","mass of #mu#mu (genParticles); GEN M_{#mu^{+}#mu^{-}} (GeV)",60.,40.,160.);
  
  h_GENP_mll50_pt_ZfromElectrons_     = fs->make<TH1F>("GENP_mll50_pt_ZfromElectrons","p_{T} of #it{ee} (genParticles); GEN p_{T}(e ^{+}e ^{-}) (GeV)",nbin_pt,pt_min,pt_Zmax);  
  h_GENP_mll50_eta_ZfromElectrons_    = fs->make<TH1F>("GENP_mll50_eta_ZfromElectrons","#eta of #it{ee} (genParticles); GEN #eta(e^{+}e^{-})",nbin_eta,eta_min,eta_max); 
  h_GENP_mll50_mass_ZfromElectrons_   = fs->make<TH1F>("GENP_mll50_mass_ZfromElectrons","mass of #it{ee} (genParticles); GEN M_{e^{+}e^{-}} (GeV)",60.,40.,160.);

  // leptons distributions
  h_GENP_mll50_pt_Lept1_fromZ_        = fs->make<TH1F>("GENP_mll50_pt_Lept1_fromZ","p_{T} of l (genParticles); GEN p_{T}(l^{lead}) (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_mll50_eta_Lept1_fromZ_       = fs->make<TH1F>("GENP_mll50_eta_Lept1_fromZ","#eta of l (genParticles); GEN #eta(l^{lead})",nbin_eta,eta_min,eta_max); 

  h_GENP_mll50_pt_Muons1_fromZ_       = fs->make<TH1F>("GENP_mll50_pt_Muons1_fromZ","p_{T} of #mu (genParticles); GEN p_{T}(#mu^{lead}) (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_mll50_eta_Muons1_fromZ_      = fs->make<TH1F>("GENP_mll50_eta_Muons1_fromZ","#eta of #mu (genParticles); GEN #eta(#mu^{lead})",nbin_eta,eta_min,eta_max); 

  h_GENP_mll50_pt_Electrons1_fromZ_   = fs->make<TH1F>("GENP_mll50_pt_Electrons1_fromZ","p_{T} of #it{e} (genParticles); GEN p_{T}(e ^{lead}) (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_mll50_eta_Electrons1_fromZ_  = fs->make<TH1F>("GENP_mll50_eta_Electrons1_fromZ","#eta of #it{e} (genParticles); GEN #eta(e^{lead})",nbin_eta,eta_min,eta_max); 

  h_GENP_mll50_pt_Lept2_fromZ_        = fs->make<TH1F>("GENP_mll50_pt_Lept2_fromZ","p_{T} of l (genParticles); GEN p_{T}(l^{2nd}) (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_mll50_eta_Lept2_fromZ_       = fs->make<TH1F>("GENP_mll50_eta_Lept2_fromZ","#eta of l (genParticles); GEN #eta(l^{2nd})",nbin_eta,eta_min,eta_max); 

  h_GENP_mll50_pt_Muons2_fromZ_       = fs->make<TH1F>("GENP_mll50_pt_Muons2_fromZ","p_{T} of #mu (genParticles); GEN p_{T}(#mu^{2nd}) (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_mll50_eta_Muons2_fromZ_      = fs->make<TH1F>("GENP_mll50_eta_Muons2_fromZ","#eta of #mu (genParticles); GEN #eta(#mu^{2nd})",nbin_eta,eta_min,eta_max); 

  h_GENP_mll50_pt_Electrons2_fromZ_   = fs->make<TH1F>("GENP_mll50_pt_Electrons2_fromZ","p_{T} of #it{e} (genParticles); GEN p_{T}(e ^{2nd}) (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_mll50_eta_Electrons2_fromZ_  = fs->make<TH1F>("GENP_mll50_eta_Electrons2_fromZ","#eta of #it{e} (genParticles); GEN #eta(e^{2nd})",nbin_eta,eta_min,eta_max); 

  h_GENP_mll50_zllCosThetaCS_         = fs->make<TH1F>("GENP_mll50_llCosThetaCS","Cosine of Collins-Soper frame angle;GEN cos(#theta_{CS,l^{+}})",20,-1.,1.); 
  h_GENP_mll50_zllCosThetaStar_       = fs->make<TH1F>("GENP_mll50_llCosThetaStar","Cosine of Z decay frame l^{+} angle;GEN cos(#theta^{*}_{l^{+}})",20,-1.,1.); 

  //*************************************************************************************************
  // Total kinematics
  //*************************************************************************************************
  
  h_GENP_pxTot_                       = fs->make<TH1F>("GENP_px_Tot","Total p_{x} in the event; Total #sum p_{x} (GeV)",500,-2.5,2.5);    
  h_GENP_pyTot_                       = fs->make<TH1F>("GENP_py_Tot","Total p_{y} in the event; Total #sum p_{y} (GeV)",500,-2.5,2.5);   
  h_GENP_pzTot_                       = fs->make<TH1F>("GENP_pz_Tot","Total p_{z} in the event; Total #sum p_{z} (GeV)",500,-2.5,2.5);   
  h_GENP_ETot_                        = fs->make<TH1F>("GENP_E_Tot","Total Energy in the event; E_{total} (GeV)",500,6995.,7005.);  
  h_GENP_DeltaEcmTotalSqrtS_          = fs->make<TH1F>("GENP_DeltaEcmTotalSqrtS","#Delta(E_{tot},#sqrt{s}) in the event;#Delta(E_{tot},#sqrt{s}) (GeV)",500,-2.5,2.5);  
  h_GENP_DeltaEcmTotalSqrtS_Vs_NPart_ = fs->make<TH2F>("GENP_DeltaEcmTotalSqrtS_Vs_NPart","#Delta(E_{tot},#sqrt{s}) in the event vs N_{part};N_{part};#Delta(E_{tot},#sqrt{s}) (GeV)",500,-0.5,499.5,500,-2.5,2.5); 

  //*************************************************************************************************
  // GenJet b plots
  //*************************************************************************************************

  h_GENP_pt_b_                  = fs->make<TH1F>("GENP_pt_b","p_{T} b-jet in the event (genParticles); GEN p_{T} b-jets (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_eta_b_                 = fs->make<TH1F>("GENP_eta_b","#eta b-jet in the event (genParticles); GEN #eta b-jet",nbin_eta,-5,5); 
  h_GENP_mass_b_                = fs->make<TH1F>("GENP_mass_b","mass b-jet in the event (genParticles); GEN m_{b-jet} (GeV)",100,0.,50.); 
  h_GENP_maxDistance_b_         = fs->make<TH1F>("GENP_maxDistance_b","jet size = #Delta R_{max}(part,b-jet); #Delta R_{max}(part,b-jet)",100,0.,2.);  

  h2_GENP_massb_vs_etab_        = fs->make<TH2F>("h2_GENP_massb_vs_etab","mass b-genJet vs #eta b-genJet; GEN #eta b-jets; GEN m_{b-jet} (GeV) ",nbin_eta,-5,5,100,0.,50.);  
  h2_GENP_massb_vs_ptb_         = fs->make<TH2F>("h2_GENP_massb_vs_ptb","mass b-genJet vs p_{T} b-genJet; GEN p_{T} b-jets (GeV); GEN m_{b-jet} (GeV)",nbin_pt,pt_min,pt_max,100,0.,50.);  
  h2_GENP_massb_vs_maxDistanceb_= fs->make<TH2F>("h2_GENP_massb_vs_maxdist","mass b-genJet vs #Delta R_{max}(part,b-jet);#Delta R_{max}(part,b-jet); GEN m_{b-jet} (GeV)",100,0.,2.,100,0.,50.);  

  p_GENP_massb_vs_etab_         = fs->make<TProfile>("p_GENP_massb_vs_etab","mass b-genJet vs #eta b-genJet; GEN #eta b-jets; GEN m_{b-jet} (GeV)",nbin_eta,-5.,5.);  	  
  p_GENP_massb_vs_ptb_          = fs->make<TProfile>("p_GENP_massb_vs_ptb","mass b-genJet vs p_{T} b-genJet; GEN p_{T} b-jets (GeV); GEN m_{b-jet} (GeV)",nbin_pt,pt_min,pt_max);  
  p_GENP_massb_vs_maxDistanceb_ = fs->make<TProfile>("p_GENP_massb_vs_maxdist","mass b-genJet vs #Delta R_{max}(part,b-jet); GEN #Delta R_{max}(part,b-jet); GEN m_{b-jet} (GeV)",100,0.,2.);  

  h_GENP_pt_b_tightmatch_       = fs->make<TH1F>("GENP_pt_b_tightmatch","p_{T} b-genjet tightly matched to b-hadron; GEN p_{T} b-jets (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_eta_b_tightmatch_      = fs->make<TH1F>("GENP_eta_b_tightmatch","#eta b-jet tightly matched to b-hadron; GEN #eta b-jet",nbin_eta,-5,5); 
  h_GENP_mass_b_tightmatch_     = fs->make<TH1F>("GENP_mass_b_tightmatch","mass b-jet tightly matched to b-hadron; GEN m_{b-jet} (GeV)",100,0.,50.); 
  h_GENP_maxDistance_b_tightmatch_ =  fs->make<TH1F>("GENP_maxDistance_b_tightmatch","jet size = #Delta R_{max}(part,b-jet); #Delta R_{max}(part,b-jet)",100,0.,2.);  
  h_GENP_mass_b_tightmatch_maxDistCut_  = fs->make<TH1F>("GENP_mass_b_tightmatch_maxDistCut","mass b-jet tightly matched to b-hadron #Delta R_{max}(part,b-jet)<0.5; GEN m_{b-jet} (GeV)",100,0.,50.); 

  h_GENP_pt_b_extratight_       = fs->make<TH1F>("GENP_pt_b_extratight","p_{T} b-genjet tightly matched to b-hadron; GEN p_{T} b-jets (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_eta_b_extratight_      = fs->make<TH1F>("GENP_eta_b_extratight","#eta b-jet tightly matched to b-hadron; GEN #eta b-jet",nbin_eta,-5,5); 
  h_GENP_mass_b_extratight_     = fs->make<TH1F>("GENP_mass_b_extratight","mass b-jet tightly matched to b-hadron; GEN m_{b-jet} (GeV)",100,0.,50.); 
  h_GENP_maxDistance_b_extratight_ = fs->make<TH1F>("GENP_maxDistance_b_extratight","jet size = #Delta R_{max}(part,b-jet); #Delta R_{max}(part,b-jet)",100,0.,2.);  

  //*************************************************************************************************
  // All GenJets
  //*************************************************************************************************

  h_GENP_pt_all_                  = fs->make<TH1F>("GENP_pt_all","p_{T} all-jet in the event (genParticles); GEN p_{T} all-jets (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_eta_all_                 = fs->make<TH1F>("GENP_eta_all","#eta all-jet in the event (genParticles); GEN #eta all-jet",nbin_eta,-5,5); 
  h_GENP_mass_all_                = fs->make<TH1F>("GENP_mass_all","mass all-jet in the event (genParticles); GEN m_{all-jet} (GeV)",100,0.,50.); 
  h_GENP_maxDistance_all_         = fs->make<TH1F>("GENP_maxDistance_all","jet size = #Delta R_{max}(part,jet); #Delta R_{max}(part,jet)",100,0.,2.);  
  h_GENP_MymaxDistance_all_       = fs->make<TH1F>("GENP_MymaxDistance_all","jet size = #Delta R_{max}(part,jet); #Delta R_{max}(part,jet)",100,0.,2.);  

  h1NGoodGenJets_                 = fs->make<TH1F>("NGoodGenJets","Number of good genJets per event; N_{jets}^{gen}; events",40,-0.5,39.5);

  h1genJetNConst_                 =fs->make<TH1F>("genJetNConstituents","Number of constituents per jet; number of constituents; genJets",40,-0.5,39.5);
  h1genJetConstPdgId_             =fs->make<TH1F>("genJetConstPdgId","Chemical composition of genJets; particle type; genJets",13,-0.5,12.5);

  h1BgenJetNConst_                 =fs->make<TH1F>("BgenJetNConstituents","Number of constituents per b-genJet; number of constituents; b-genJets",40,-0.5,39.5);
  h1BgenJetConstPdgId_             =fs->make<TH1F>("BgenJetConstPdgId","Chemical composition of b-genJets; particle type; b-genJets",13,-0.5,12.5);

  TString PdgIdBinLabels[13]  ={"e","#nu_{e}","#mu","#nu_{#mu}","#nu_{#tau}","#gamma","#pi^{+/-}","K^{0}_{L}","K^{0}_{S}","K^{+/-}","p","n","#Lambda,#Sigma,#Xi"};
  
  for(UInt_t bin=1; bin<=13; bin++){
    h1genJetConstPdgId_->GetXaxis()->SetBinLabel(bin,PdgIdBinLabels[bin-1]); 
    h1BgenJetConstPdgId_->GetXaxis()->SetBinLabel(bin,PdgIdBinLabels[bin-1]); 
  }

  nTotEvts_ = 0;
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MCTruthAnalyzer::endJob() {
  cout<<endl<<"Number of total events processed by analyzer: "<<nTotEvts_<<endl;
}

// ----------- method to decide flavour of Z in the event
std::pair<bool,bool> 
MCTruthAnalyzer::WhichFlavour(std::vector<const reco::Candidate*> theGenMuons,std::vector<const reco::Candidate*> theGenElectrons){
  
  std::pair<bool,bool> theRes = make_pair(false,false);

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
MCTruthAnalyzer::getTheZLeptons(std::vector<const reco::Candidate*> theGenMuons,std::vector<const reco::Candidate*> theGenElectrons){
  
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
MCTruthAnalyzer::hasBottom(const reco::Candidate &c) 
{
  int code1;
  int code2;
  bool tmpHasBottom = false;
  code1 = (int)( ( abs(c.pdgId() ) / 100)%10 );
  code2 = (int)( ( abs(c.pdgId() ) /1000)%10 );
  if ( code1 == 5 || code2 == 5) tmpHasBottom = true;
  return tmpHasBottom;
}

// ------------ method to tag c-hadrons ---------------------------------------------------
bool 
MCTruthAnalyzer::hasCharm(const reco::Candidate &c) 
{
  int code1;
  int code2;
  bool tmpHasCharm = false;
  code1 = (int)( ( abs(c.pdgId() ) / 100)%10 );
  code2 = (int)( ( abs(c.pdgId() ) /1000)%10 );
  if ( code1 == 4 || code2 == 4) tmpHasCharm = true;
  return tmpHasCharm;
}

// ------------------------------- sort leptons 4mom by pt -----------------------------------------

std::vector<math::XYZTLorentzVectorD>
MCTruthAnalyzer::sort4MomentaByPt(std::vector<math::XYZTLorentzVectorD> leptonMomenta_)
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
MCTruthAnalyzer::sortLeptonsByPt(std::vector<const reco::Candidate*> leptons_)
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
Double_t MCTruthAnalyzer::CosThetaStar(std::vector<const reco::Candidate*> sortedGenLeptons){

  TLorentzVector p4_lp;
  TLorentzVector p4_lm;

  Double_t pla10 = sortedGenLeptons[0]->px();
  Double_t pla11 = sortedGenLeptons[0]->py();
  Double_t pla12 = sortedGenLeptons[0]->pz();
  Double_t e1 = sortedGenLeptons[0]->energy();
  Int_t lep1Charge = sortedGenLeptons[0]->charge();
  
  Double_t pla20 = sortedGenLeptons[1]->px();
  Double_t pla21 = sortedGenLeptons[1]->py();
  Double_t pla22 = sortedGenLeptons[1]->pz();
  Double_t e2 = sortedGenLeptons[1]->energy();

  if(lep1Charge>0.){
    p4_lp.SetPxPyPzE(pla10,pla11,pla12,e1);
    p4_lm.SetPxPyPzE(pla20,pla21,pla22,e2);
  } else {
    p4_lp.SetPxPyPzE(pla20,pla21,pla22,e2);
    p4_lm.SetPxPyPzE(pla10,pla11,pla12,e1);
  }

  TLorentzVector p4_Z;
  TVector3 boost;
  
  p4_Z=p4_lp+p4_lm;
  boost=p4_Z.BoostVector();
  TVector3 deboost =-1*boost; 
  p4_lp.Boost(deboost);
  
  TVector3 p3_lp =p4_lp.Vect();
  Double_t cost = TMath::Cos(p3_lp.Angle(boost));
  
  return cost;

}

Double_t MCTruthAnalyzer::CostCS(std::vector<const reco::Candidate*> sortedGenLeptons){
  // Calculation the Collins-Soper angle (adapted from code by R. Arnaldi)
  
  Double_t fMProton = 0.93827231;
  Double_t ebeam=3500.;  //temporary
  if(ebeam<=0){
    printf("Can not compute costCS with EBeam=%f\n",ebeam);
    return -999999999;
  }
  
  Double_t mp=fMProton;
  Double_t pbeam=TMath::Sqrt(ebeam*ebeam-mp*mp);
  
  Double_t pla10 = sortedGenLeptons[0]->px();
  Double_t pla11 = sortedGenLeptons[0]->py();
  Double_t pla12 = sortedGenLeptons[0]->pz();
  Double_t e1 = sortedGenLeptons[0]->energy();
  Int_t lep1Charge = sortedGenLeptons[0]->charge();
  
  Double_t pla20 = sortedGenLeptons[1]->px();
  Double_t pla21 = sortedGenLeptons[1]->py();
  Double_t pla22 = sortedGenLeptons[1]->pz();
  Double_t e2 = sortedGenLeptons[1]->energy();
  Int_t lep2Charge = sortedGenLeptons[1]->charge();

  // Fill the Lorentz vector for projectile and target
  // For the moment we do not consider the crossing angle
  TLorentzVector pProjCM(0.,0.,-pbeam,ebeam); // projectile
  TLorentzVector pTargCM(0.,0., pbeam,ebeam); // target
  //
  // --- Get the muons parameters in the CM frame
  //
  TLorentzVector pLep1CM(pla10,pla11,pla12,e1);
  TLorentzVector pLep2CM(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  TLorentzVector pDilepCM=pLep1CM+pLep2CM;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta=(-1./pDilepCM.E())*pDilepCM.Vect();
  TLorentzVector pLep1Dilep=pLep1CM;
  TLorentzVector pLep2Dilep=pLep2CM;
  TLorentzVector pProjDilep=pProjCM;
  TLorentzVector pTargDilep=pTargCM;
  pLep1Dilep.Boost(beta);
  pLep2Dilep.Boost(beta);
  pProjDilep.Boost(beta);
  pTargDilep.Boost(beta);
  //
  // --- Determine the z axis for the CS angle 
  //
  TVector3 zaxisCS=(((pProjDilep.Vect()).Unit())-((pTargDilep.Vect()).Unit())).Unit();
  //
  // --- Determine the CS angle (angle between lep+ and the z axis defined above)
  //
  Double_t cost;
  if(lep1Charge > 0) {
    cost = zaxisCS.Dot((pLep1Dilep.Vect()).Unit());
    // Theta CS is not properly defined for Like-Sign lepons
    if(lep2Charge > 0 && cost<0) cost=-cost;
  } else { 
    // Theta CS is not properly defined for Like-Sign lepons
    cost = zaxisCS.Dot((pLep2Dilep.Vect()).Unit());
    if(lep2Charge < 0 && cost<0) cost=-cost;
  }
  return cost;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MCTruthAnalyzer);
