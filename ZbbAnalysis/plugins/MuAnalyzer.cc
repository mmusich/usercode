//class MuAnalyzer

#include "MuAnalyzer.h"

// Collaborating Class Header
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

// for "luminosity"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

using namespace std;
using namespace edm;

/// Constructor
MuAnalyzer::MuAnalyzer(const ParameterSet& pset) : 
  numEventsNames_(pset.getUntrackedParameter< vector<string> >("numEventsNames")),
  theMuonLabel_(pset.getUntrackedParameter<string>("MuonCollection")),
  theJetPFLabel_(pset.getUntrackedParameter<string>("JetCollectionPF")),
  theZllLabel_(pset.getUntrackedParameter<string>("ZllCollection"))
{
}

/// Destructor
MuAnalyzer::~MuAnalyzer(){
}

void MuAnalyzer::beginJob(){

  // Book histograms
  edm::Service<TFileService> fs;

  //Creating subdirectories
  TFileDirectory counters = fs->mkdir("Counters");
  TFileDirectory spectra = fs->mkdir("Spectra");
  TFileDirectory DiMuInvMass = spectra.mkdir("DiMuInvMass");
  TFileDirectory All = spectra.mkdir("All muons");
  TFileDirectory Global = spectra.mkdir("Global muons");
  TFileDirectory GlobalMM = spectra.mkdir("Global muons many muons");
  TFileDirectory TK = spectra.mkdir("Only TK muons");
  TFileDirectory SA = spectra.mkdir("Only SA muons");
  TFileDirectory TKandSA = spectra.mkdir("!(Global) but TK && SA muons");

  //Histograms counters
  unsigned int nbin=numEventsNames_.size();
  h1SelectedEvts_ = counters.make<TH1F>("SelectedEvts","Selected events",nbin,-0.5,-0.5+nbin);  
  for (UInt_t i = 0;i < nbin; i++) h1SelectedEvts_->GetXaxis()->SetBinLabel(i+1,(numEventsNames_[i].c_str()));

  //N histograms
  h1NMuons_ = counters.make<TH1I>("NMuons","Number of muons per event",20,-0.5,19.5);
  h1NZll_ = counters.make<TH1I>("NZll","Number of Z->ll per event",20,-0.5,19.5);
  
  //Spectra histograms
  h1DiMuMassAll_ = DiMuInvMass.make<TH1F>("DiMuMassAll","Invariant mass spectrum of all dimuons",100,0,120);
  h1DiMuMassPresel_ = DiMuInvMass.make<TH1F>("DiMuMassPresel","Invariant mass spectrum of Presel dimuons",100,0,120);

  //All muons
  h1PtRecoAllMuon_ = All.make<TH1F>("pTAllMuons","p_{T}^{rec}",250,0,120);
  h1EtaRecoAllMuon_ = All.make<TH1F>("etaAllMuons","p_{T}^{rec}",20,0,4);
  h1PhiRecoAllMuon_ = All.make<TH1F>("phiAllMuons","p_{T}^{rec}",50,0,3.15);
  h1DxyRecoAllMuon_ = All.make<TH1F>("dxyAllMuons","d_{xy}^{rec}",400,-200,200);
  h1DzRecoAllMuon_ = All.make<TH1F>("dzAllMuons","d_{z}^{rec}",400,-200,200);
  h2EtaPtRecoAllMuon_ = All.make<TH2F>("EtaPtAllMuon","Eta vs Pt",20,0,4,250,0,120);
  //Global muons
  h1PtRecoGlobalMuon_ = Global.make<TH1F>("pTGlobalMuons","p_{T}^{rec}",250,0,120);
  h1EtaRecoGlobalMuon_ = Global.make<TH1F>("etaGlobalMuons","p_{T}^{rec}",20,0,4);
  h1PhiRecoGlobalMuon_ = Global.make<TH1F>("phiGlobalMuons","p_{T}^{rec}",50,0,3.15);
  h1DxyRecoGlobalMuon_ = Global.make<TH1F>("dxyGlobalMuons","d_{xy}^{rec}",400,-200,200);
  h1DzRecoGlobalMuon_ = Global.make<TH1F>("dzGlobalMuons","d_{z}^{rec}",400,-200,200);
  h2EtaPtRecoGlobalMuon_ = Global.make<TH2F>("EtaPtGlobalMuon","Eta vs Pt",20,0,4,250,0,120);

  //Global muons Many muons (MM)
  h1PtRecoGlobalMuonMM_ = GlobalMM.make<TH1F>("pTGlobalMuonsMM","p_{T}^{rec}",250,0,120);
  h1EtaRecoGlobalMuonMM_ = GlobalMM.make<TH1F>("etaGlobalMuonsMM","p_{T}^{rec}",20,0,4);
  h1PhiRecoGlobalMuonMM_ = GlobalMM.make<TH1F>("phiGlobalMuonsMM","p_{T}^{rec}",50,0,3.15);
  h1DxyRecoGlobalMuonMM_ = GlobalMM.make<TH1F>("dxyGlobalMuonsMM","d_{xy}^{rec}",400,-200,200);
  h1DzRecoGlobalMuonMM_ = GlobalMM.make<TH1F>("dzGlobalMuonsMM","d_{z}^{rec}",400,-200,200);
  h2EtaPtRecoGlobalMuonMM_ = GlobalMM.make<TH2F>("EtaPtGlobalMuonMM","Eta vs Pt",20,0,4,250,0,120);

  //Only Tracker muons
  h1PtRecoTrackerMuon_ = TK.make<TH1F>("pTTrackerMuons","p_{T}^{rec}",250,0,120);
  h1EtaRecoTrackerMuon_ = TK.make<TH1F>("etaTrackerMuons","p_{T}^{rec}",20,0,4);
  h1PhiRecoTrackerMuon_ = TK.make<TH1F>("phiTrackerMuons","p_{T}^{rec}",50,0,3.15);
  h1DxyRecoTrackerMuon_ = TK.make<TH1F>("dxyTrackerMuons","d_{xy}^{rec}",400,-200,200);
  h1DzRecoTrackerMuon_ = TK.make<TH1F>("dzTrackerMuons","d_{z}^{rec}",400,-200,200);
  h2EtaPtRecoTrackerMuon_ = TK.make<TH2F>("EtaPtTrackerMuon","Eta vs Pt",20,0,4,250,0,120);
  //Only SA muons
  h1PtRecoStandAloneMuon_ = SA.make<TH1F>("pTStandAloneMuons","p_{T}^{rec}",250,0,120);
  h1EtaRecoStandAloneMuon_ = SA.make<TH1F>("etaStandAloneMuons","p_{T}^{rec}",20,0,4);
  h1PhiRecoStandAloneMuon_ = SA.make<TH1F>("phiStandAloneMuons","p_{T}^{rec}",50,0,3.15);
  h1DxyRecoStandAloneMuon_ = SA.make<TH1F>("dxyStandAloneMuons","d_{xy}^{rec}",400,-200,200);
  h1DzRecoStandAloneMuon_ = SA.make<TH1F>("dzStandAloneMuons","d_{z}^{rec}",400,-200,200);
  h2EtaPtRecoStandAloneMuon_ = SA.make<TH2F>("EtaPtStandAloneMuon","Eta vs Pt",20,0,4,250,0,120);
  //!(Global) but TK && SA muons
  h1PtRecoTKandSAMuon_ = TKandSA.make<TH1F>("pTTKandSAMuons","p_{T}^{rec}",250,0,120);
  h1EtaRecoTKandSAMuon_ = TKandSA.make<TH1F>("etaTKandSAMuons","p_{T}^{rec}",20,0,4);
  h1PhiRecoTKandSAMuon_ = TKandSA.make<TH1F>("phiTKandSAMuons","p_{T}^{rec}",50,0,3.15);
  h1DxyRecoTKandSAMuon_ = TKandSA.make<TH1F>("dxyTKandSAMuons","d_{xy}^{rec}",400,-200,200);
  h1DzRecoTKandSAMuon_ = TKandSA.make<TH1F>("dzTKandSAMuons","d_{z}^{rec}",400,-200,200);
  h2EtaPtRecoTKandSAMuon_ = TKandSA.make<TH2F>("EtaPtTKandSAMuon","Eta vs Pt",20,0,4,250,0,120);

  //Check muon type histograms
  h1CheckMu_ = counters.make<TH1F>("CheckMu","Check if Mu are SA,Global or Tracker",4,-0.5,3.5);
  h1CheckMu_->GetXaxis()->SetBinLabel(1,"Global Muons");
  h1CheckMu_->GetXaxis()->SetBinLabel(2,"Only Tracker Muons");
  h1CheckMu_->GetXaxis()->SetBinLabel(3,"Only StandAlone Muons");
  h1CheckMu_->GetXaxis()->SetBinLabel(4,"!(Global) but SA && TK");
  h1MuType_ = counters.make<TH1F>("MuType","muon->type() response",25,-0.5,24.5);


}

void MuAnalyzer::endJob(){

}

void MuAnalyzer::analyze(const Event & iEvent, const EventSetup& eventSetup)
{
  

  // Get the Muon collection
  Handle<pat::MuonCollection> muonsHandle;
  iEvent.getByLabel(theMuonLabel_, muonsHandle);

  // Get the Z collection
  Handle<reco::CompositeCandidateCollection> zllHandle;
  iEvent.getByLabel(theZllLabel_, zllHandle);
  const reco::CompositeCandidateCollection & zll = *(zllHandle.product());

  //Get JetPF collection
  Handle<pat::JetCollection> jetsPFHandle;
  iEvent.getByLabel(theJetPFLabel_, jetsPFHandle);

  pat::MuonCollection::const_iterator muon1,muon2;
  reco::Candidate::LorentzVector p1,p2,sum;

  //Fill histogram DiMuMassAll
  if(muonsHandle->size()>=2){
    
    for (muon1 = muonsHandle->begin(); muon1 != muonsHandle->end(); muon1++){
      for (muon2 = muonsHandle->begin(); muon2 < muon1; muon2++) {
	p1 = muon1->p4();
	p2 = muon2->p4();  
	sum = p1 + p2;
	if((muon1->charge()*muon2->charge())<0) h1DiMuMassAll_->Fill(sum.mass());
      }
    }

  }
  
  //Fill histogram DiMuMassPresel
  if(muonsHandle->size()>=2){
    
    for (muon1 = muonsHandle->begin(); muon1 != muonsHandle->end(); muon1++){
      
      if ( muon1->pt() > 10. && 
	   muon1->isGlobalMuon()==true && muon1->isTrackerMuon()==true && 
	   muon1->normChi2() < 15 && muon1->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 && muon1->innerTrack()->hitPattern().numberOfValidPixelHits() > 0  &&
	   muon1->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && 
	   muon1->dB() < 0.2 && (muon1->trackIso()+muon1->caloIso()) < 0.15*muon1->pt() && muon1->numberOfMatches() > 1 && abs(muon1->eta()) < 2.1){
	
	for (muon2 = muonsHandle->begin(); muon2 < muon1; muon2++) {
	
	if ( muon2->pt() > 10. && 
	     muon2->isGlobalMuon()==true && muon2->isTrackerMuon()==true && 
	     muon2->normChi2() < 15 && muon2->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 && muon2->innerTrack()->hitPattern().numberOfValidPixelHits() > 0  &&
	     muon2->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && 
	     muon2->dB() < 0.2 && (muon2->trackIso()+muon2->caloIso()) < 0.15*muon2->pt() && muon2->numberOfMatches() > 1 && abs(muon2->eta()) < 2.1){
	  p1 = muon1->p4();
	  p2 = muon2->p4();  
	  sum = p1 + p2;
	  if((muon1->charge()*muon2->charge())<0) h1DiMuMassPresel_->Fill(sum.mass());
	}
      }
      }
    }
  }


  // Fill histogram NMuons
  h1NMuons_->Fill(muonsHandle->size());

  // Fill histogram NZll
  h1NZll_->Fill(zllHandle->size());



  // Fill some spectra histograms
  for (muon1 = muonsHandle->begin();  muon1 != muonsHandle->end(); ++muon1){

    //Fill histogram MuType
    h1MuType_->Fill(muon1->type());

    float etaMu,ptMu,phiMu;
    etaMu = muon1->eta();
    phiMu = muon1->phi();
    ptMu = muon1->pt();

    //Fill histogram spectra all muons
    h1PtRecoAllMuon_->Fill(ptMu);
    h1EtaRecoAllMuon_->Fill(etaMu);
    h1PhiRecoAllMuon_->Fill(phiMu);
    h2EtaPtRecoAllMuon_->Fill(etaMu,ptMu);
    

    //Spectra of Global muons
    if( muon1->isGlobalMuon() ){

      h1PtRecoGlobalMuon_->Fill(ptMu);
      h1EtaRecoGlobalMuon_->Fill(etaMu);
      h1PhiRecoGlobalMuon_->Fill(phiMu);
      h1DxyRecoGlobalMuon_->Fill(muon1->globalTrack()->d0());
      h1DzRecoGlobalMuon_->Fill(muon1->globalTrack()->dz());
      h1DxyRecoAllMuon_->Fill(muon1->globalTrack()->d0());
      h1DzRecoAllMuon_->Fill(muon1->globalTrack()->dz());
      h2EtaPtRecoGlobalMuon_->Fill(etaMu,ptMu);



      if ( jetsPFHandle->size() > 1 ) {
	h1PtRecoGlobalMuonMM_->Fill(ptMu);
	h1EtaRecoGlobalMuonMM_->Fill(etaMu);
	h1PhiRecoGlobalMuonMM_->Fill(phiMu);
	h1DxyRecoGlobalMuonMM_->Fill(muon1->globalTrack()->d0());
	h1DzRecoGlobalMuonMM_->Fill(muon1->globalTrack()->dz());
	h2EtaPtRecoGlobalMuonMM_->Fill(etaMu,ptMu);
      }

      
      h1CheckMu_->Fill(0);

    }// end if Global
    
    else {
      
      //Spectra of Only TK muons
      if ( muon1->isTrackerMuon() && !(muon1->isStandAloneMuon()) ){
	h1PtRecoTrackerMuon_->Fill(ptMu);
	h1EtaRecoTrackerMuon_->Fill(etaMu);
	h1PhiRecoTrackerMuon_->Fill(phiMu);
	h1DxyRecoTrackerMuon_->Fill(muon1->innerTrack()->d0());
	h1DzRecoTrackerMuon_->Fill(muon1->innerTrack()->dz());
	h1DxyRecoAllMuon_->Fill(muon1->innerTrack()->d0());
	h1DzRecoAllMuon_->Fill(muon1->innerTrack()->dz());
	h2EtaPtRecoTrackerMuon_->Fill(etaMu,ptMu);
	
	h1CheckMu_->Fill(1);
      }


      //Spectra of Only SA muons
      if ( muon1->isStandAloneMuon() && !(muon1->isTrackerMuon()) ){
	h1PtRecoStandAloneMuon_->Fill(ptMu);
	h1EtaRecoStandAloneMuon_->Fill(etaMu);
	h1PhiRecoStandAloneMuon_->Fill(phiMu);
	h1DxyRecoStandAloneMuon_->Fill(muon1->outerTrack()->d0());
	h1DzRecoStandAloneMuon_->Fill(muon1->outerTrack()->dz());
	h1DxyRecoAllMuon_->Fill(muon1->outerTrack()->d0());
	h1DzRecoAllMuon_->Fill(muon1->outerTrack()->dz());
	h2EtaPtRecoStandAloneMuon_->Fill(etaMu,ptMu);
	
	h1CheckMu_->Fill(2);
      }

      //Spectra of TK && SA but !(Global)  muons
      if ( muon1->isTrackerMuon() && muon1->isStandAloneMuon() ){
	h1PtRecoTKandSAMuon_->Fill(ptMu);
	h1EtaRecoTKandSAMuon_->Fill(etaMu);
	h1PhiRecoTKandSAMuon_->Fill(phiMu);
	h1DxyRecoTKandSAMuon_->Fill(muon1->innerTrack()->d0());
	h1DzRecoTKandSAMuon_->Fill(muon1->innerTrack()->dz());
	h1DxyRecoAllMuon_->Fill(muon1->innerTrack()->d0());
	h1DzRecoAllMuon_->Fill(muon1->innerTrack()->dz());
	h2EtaPtRecoTKandSAMuon_->Fill(etaMu,ptMu);
	
	h1CheckMu_->Fill(3);
      }



    }// end else !(Global)
    
    
  }//end loop on MuonCollection
  
}

// based on https://hypernews.cern.ch/HyperNews/CMS/get/muon/433/1/1.html
void MuAnalyzer::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup) 
{

  for (size_t i = 0; i < numEventsNames_.size(); i++) {
    
    string name = numEventsNames_[i];
    Handle<MergeableCounter> numEventsCounter;
    iLumi.getByLabel(name, numEventsCounter);
    
    if (numEventsCounter.isValid()) 
      h1SelectedEvts_->AddBinContent(i + 1, numEventsCounter->value); 
  }
}


DEFINE_FWK_MODULE(MuAnalyzer);







