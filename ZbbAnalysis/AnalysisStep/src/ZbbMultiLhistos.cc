#include "ZbbAnalysis/AnalysisStep/interface/EventCategory.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbMultiLhistos.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbUtils.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include "TCut.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <Math/VectorUtil.h>

#include <iostream>
#include <algorithm>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////
//ZbbAllLeptons methods
//////////////////////////////////////////////////////////////////////////////////

void ZbbAllLeptons::book(){

  edm::Service<TFileService> fs;
  std::string DirName = &(*N_);
  TFileDirectory HistoVars = fs->mkdir(DirName);

  TH1F::SetDefaultSumw2(kTRUE);

  //Muons
  h_nmu_                = HistoVars.make<TH1F>(N_+"_nmu","Number of muons in the event; N_{#mu}",15,-0.5,14.5);
  h_allmu_pt_           = HistoVars.make<TH1F>(N_+"_allmu_pt","p_{T} all muons;#mu p_{T} (GeV/c)",50,0,100); 
  h_allmu_eta_          = HistoVars.make<TH1F>(N_+"_allmu_eta","#eta all muons;#mu #eta",40,-4,4);
  h_allmu_IPT_          = HistoVars.make<TH1F>(N_+"_allmu_IPT","d_{xy} muons; d_{xy} (cm)",50,-2,2);
  h_allmu_IPL_          = HistoVars.make<TH1F>(N_+"_allmu_IPL","d_{z} muons; d_{z} (cm)",4000,-20,20);
  h_allmu_SIP_          = HistoVars.make<TH1F>(N_+"_allmu_SIP","muon IP/#sigma_{IP};muon IP/#sigma_{IP}",60,0.,20.); 
  h_allmu_dr_lj_        = HistoVars.make<TH1F>(N_+"_allmu_dr_lj","#Delta R_{min} muons / all jets; #Delta R (#mu,j)_{min}",100,0,10);
  h_allmu_dr_lbj_       = HistoVars.make<TH1F>(N_+"_allmu_dr_lbj","#Delta R_{min} muons / B-tagged jets; #Delta R (#mu,j_{b})_{min}",100,0,10);
  h_allmu_ptrel_lj_     = HistoVars.make<TH1F>(N_+"_allmu_ptrel_lj","p_{T}^{rel} min muons wrt jets; p_{T}^{rel} (GeV/c)",60,0,20);
  h_allmu_ptrel_lbj_    = HistoVars.make<TH1F>(N_+"_allmu_ptrel_lbj","p_{T}^{rel} min muons wrt b-jets; p_{T}^{rel} (GeV/c)",60,0,20);
  h_allmu_trackIso_     = HistoVars.make<TH1F>(N_+"_allmu_trackIso","trackIso/p_{T} muons; trackIso/p_{T}",50,0,5);
  h_allmu_caloIso_      = HistoVars.make<TH1F>(N_+"_allmu_caloIso","caloIso/p_{T} muons; caloIso/p_{T}",50,0,5);
  
  //Electrons
  h_nele_                = HistoVars.make<TH1F>(N_+"_nele","Number of electrons in the event; N_{e}",15,-0.5,14.5);
  h_allele_pt_           = HistoVars.make<TH1F>(N_+"_allele_pt","p_{T} electrons;e p_{T} (GeV/c)",50,0,100); 
  h_allele_eta_          = HistoVars.make<TH1F>(N_+"_allele_eta","#eta electrons;e #eta",40,-4,4);
  h_allele_IPT_          = HistoVars.make<TH1F>(N_+"_allele_IPT","d_{xy} electrons; d_{xy} (cm)",50,-2,2);
  h_allele_IPL_          = HistoVars.make<TH1F>(N_+"_allele_IPL","d_{z} electrons; d_{z} (cm)",4000,-20,20);
  h_allele_SIP_          = HistoVars.make<TH1F>(N_+"_allele_SIP","electron IP/#sigma_{IP};electron IP/#sigma_{IP}",60,0.,20.); 
  h_allele_dr_lj_        = HistoVars.make<TH1F>(N_+"_allele_drjmin","#Delta R_{min} electrons / all jets; #Delta R (e,j)_{min}",100,0,10);
  h_allele_dr_lbj_       = HistoVars.make<TH1F>(N_+"_allele_drBjmin","#Delta R_{min} electrons / B-tagged jets; #Delta R (e,j_{b})_{min}",100,0,10);
  h_allele_ptrel_lj_     = HistoVars.make<TH1F>(N_+"_allele_ptrel_lj","p_{T}^{rel} min electrons wrt jets; p_{T}^{rel} (GeV/c)",60,0,20);
  h_allele_ptrel_lbj_    = HistoVars.make<TH1F>(N_+"_allele_ptrel_lbj","p_{T}^{rel} min electrons wrt b-jets; p_{T}^{rel} (GeV/c)",60,0,20);
  h_allele_trackIso_     = HistoVars.make<TH1F>(N_+"_allele_trackIso","trackIso/p_{T} electrons; trackIso/p_{T}",50,0,5);
  h_allele_caloIso_      = HistoVars.make<TH1F>(N_+"_allele_caloIso","caloIso/p_{T} electrons; caloIso/p_{T}",50,0,5);
    
}

void ZbbAllLeptons::fillMu(const EventCategory& ec_mu, const EventCategory& ec_ele, edm::Handle<edm::View<pat::Muon> > muons, edm::Handle<edm::View<pat::Jet> > jets,std::string bTagAlgoWP){

  std::vector<pat::Jet> theGoodJets;
  std::vector<pat::Jet> theGoodBJets;
  
  Bool_t thecut_;   
  EventCategory ec; 

  if(ec_mu.isLLTightAndTriggerMatched()){
    thecut_= ZbbUtils::ProduceTheCut(ec_mu,N_,"");
    ec = ec_mu;
  } else {
    thecut_ = ZbbUtils::ProduceTheCut(ec_ele,N_,"");
    ec = ec_ele;
  }
  
  Double_t w = ZbbUtils::getTheWeight(ec,N_);

  // const reco::Vertex& theZVertex = Zvertexes->at(0);
  const reco::Vertex& theZVertex = ec.theZvertex_;

  for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
    if ( ZbbUtils::isJetIdOk(*jet,"loose") && jet->pt() > lCuts_.bjetPtMin_ && jet->eta() < lCuts_.bjetEtaMax_ ){
      theGoodJets.push_back(*jet);
      if (ZbbUtils::isBJet(*jet,bTagAlgoWP))theGoodBJets.push_back(*jet);
    }
  }//end loop over jets
  
  if (thecut_){
    
    if (muons->size()!=0){
      for(edm::View<pat::Muon>::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon){
	h_allmu_pt_->Fill(muon->pt(),w);
	h_allmu_eta_->Fill(muon->eta(),w);
	h_allmu_trackIso_->Fill(muon->trackIso()/muon->pt(),w);
	h_allmu_caloIso_->Fill(muon->caloIso()/muon->pt(),w);
	h_allmu_SIP_->Fill(fabs(muon->dB(pat::Muon::PV3D))/muon->edB(pat::Muon::PV3D),w); 

	if ( !(muon->track().isNull()) ){
	  reco::TrackRef theTrack = muon->track();
	  const math::XYZPoint myVertex(theZVertex.position().x(),theZVertex.position().y(),theZVertex.position().z());
	  double IPT = theTrack->dxy(myVertex);
	  double IPL = theTrack->dz(myVertex);
	  h_allmu_IPT_->Fill(IPT,w);
	  h_allmu_IPL_->Fill(IPL,w);
	}//if !(track().isNull())

	Float_t deltaRmin = 1000;
	pat::Jet closestJet;
	if(theGoodJets.size()!=0){
	  for(size_t j = 0;j < theGoodJets.size();j++){
	    Double_t deltaR   = ROOT::Math::VectorUtil::DeltaR(theGoodJets[j].momentum(),muon->momentum());
	    if (deltaR < deltaRmin) {
	      closestJet = theGoodJets[j];
	      deltaRmin = deltaR;
	    }
	  }//end loop over jets
	}//end GoodJets.size() veto
	h_allmu_dr_lj_->Fill(deltaRmin,w);
	
	Float_t deltaRBmin = 1000;
	pat::Jet closestBJet;
	if(theGoodBJets.size()!=0){
	  for(size_t j = 0;j < theGoodBJets.size();j++){
	    Double_t deltaRB   = ROOT::Math::VectorUtil::DeltaR(theGoodBJets[j].momentum(),muon->momentum());
	    if (deltaRB < deltaRBmin) {
	      closestBJet = theGoodBJets[j];
	      deltaRBmin = deltaRB;
	    }
	  }//end loop over B-jets
	}//end GoodBJets.size() veto
	h_allmu_dr_lbj_->Fill(deltaRBmin,w);	

	//PtRel l-jets
	math::XYZVector jetDir = closestJet.momentum().Unit();
	Double_t ptrel_lj = ROOT::Math::VectorUtil::Perp(muon->momentum(),jetDir);
	//TVector3 muvec(muon->momentum().X(),muon->momentum().Y(),muon->momentum().Z());
	//TVector3 jetvec(closestJet.momentum().X(),closestJet.momentum().Y(),closestJet.momentum().Z());
	//Float_t ptrel_lj = muvec.Perp(jetvec);
	if (deltaRmin < 0.5) h_allmu_ptrel_lj_->Fill(ptrel_lj,w);

	//PtRel l-bjets
	math::XYZVector bjetDir = closestBJet.momentum().Unit();
	Double_t ptrel_lbj = ROOT::Math::VectorUtil::Perp(muon->momentum(),bjetDir);
	//TVector3 bjetvec(closestBJet.momentum().X(),closestBJet.momentum().Y(),closestBJet.momentum().Z());
	//Float_t ptrel_lbj = muvec.Perp(bjetvec);
	if (deltaRBmin < 0.5) h_allmu_ptrel_lbj_->Fill(ptrel_lbj,w);
      }//end loop over muons
    }
    h_nmu_->Fill(muons->size());

  }//end if right event category

  return;
}

void ZbbAllLeptons::fillEle(const EventCategory& ec_mu, const EventCategory& ec_ele, edm::Handle<edm::View<pat::Electron> > electrons, edm::Handle<edm::View<pat::Jet> > jets,std::string bTagAlgoWP){
  
  std::vector<pat::Jet> theGoodJets;;
  std::vector<pat::Jet> theGoodBJets;

  Bool_t thecut_;   
  EventCategory ec; 
  
  if(ec_mu.isLLTightAndTriggerMatched()){
    thecut_= ZbbUtils::ProduceTheCut(ec_mu,N_,"");
    ec = ec_mu;
  } else {
    thecut_ = ZbbUtils::ProduceTheCut(ec_ele,N_,"");
    ec = ec_ele;
  }

  Double_t w = ZbbUtils::getTheWeight(ec,N_);

  // const reco::Vertex& theZVertex = Zvertexes->at(0);
  const reco::Vertex& theZVertex = ec.theZvertex_;
  
  for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
    if ( ZbbUtils::isJetIdOk(*jet,"loose") && jet->pt() > lCuts_.bjetPtMin_ && jet->eta() < lCuts_.bjetEtaMax_ ){
      theGoodJets.push_back(*jet);
      if (ZbbUtils::isBJet(*jet,bTagAlgoWP))theGoodBJets.push_back(*jet);
    }
  }//end loop over jets
  
  if (thecut_){
    
    if (electrons->size()!=0){
      for(edm::View<pat::Electron>::const_iterator electron=electrons->begin(); electron!=electrons->end(); ++electron){
	h_allele_pt_->Fill(electron->pt(),w);
	h_allele_eta_->Fill(electron->eta(),w);
	h_allele_trackIso_->Fill(electron->trackIso()/electron->pt(),w);
	h_allele_caloIso_->Fill(electron->caloIso()/electron->pt(),w);
	h_allele_SIP_->Fill(fabs(electron->dB(pat::Electron::PV3D))/electron->edB(pat::Electron::PV3D),w); 

	if ( !(electron->gsfTrack().isNull()) ){
	  reco::GsfTrackRef theTrack = electron->gsfTrack();
	  const math::XYZPoint myVertex(theZVertex.position().x(),theZVertex.position().y(),theZVertex.position().z());
	  double IPT = theTrack->dxy(myVertex);
	  double IPL = theTrack->dz(myVertex);
	  h_allele_IPT_->Fill(IPT,w);
	  h_allele_IPL_->Fill(IPL,w);
	}//if !(track().isNull())

	Float_t deltaRmin = 1000;
	pat::Jet closestJet;
	if(theGoodJets.size()!=0){
	  for(size_t j = 0;j < theGoodJets.size();j++){
	    Double_t deltaR   = ROOT::Math::VectorUtil::DeltaR(theGoodJets[j].momentum(),electron->momentum());
	    if (deltaR < deltaRmin) {
	      closestJet = theGoodJets[j];
	      deltaRmin = deltaR;
	    }
	  }//end loop over jets
	}//end GoodJets.size() veto
	h_allele_dr_lj_->Fill(deltaRmin,w);
	
	Float_t deltaRBmin = 1000;
	pat::Jet closestBJet;
	if(theGoodBJets.size()!=0){
	  for(size_t j = 0;j < theGoodBJets.size();j++){
	    Double_t deltaRB   = ROOT::Math::VectorUtil::DeltaR(theGoodBJets[j].momentum(),electron->momentum());
	    if (deltaRB < deltaRBmin) {
	      closestBJet = theGoodBJets[j];
	      deltaRBmin = deltaRB;
	    }
	  }//end loop over B-jets
	}//end GoodBJets.size() veto
	h_allele_dr_lbj_->Fill(deltaRBmin,w);	

	//PtRel l-jets
	math::XYZVector jetDir = closestJet.momentum().Unit();
	Double_t ptrel_lj = ROOT::Math::VectorUtil::Perp(electron->momentum(),jetDir);
	//TVector3 elevec(electron->momentum().X(),electron->momentum().Y(),electron->momentum().Z());
	//TVector3 jetvec(closestJet.momentum().X(),closestJet.momentum().Y(),closestJet.momentum().Z());
	//Float_t ptrel_lj = elevec.Perp(jetvec);
	if (deltaRmin < 0.5) h_allele_ptrel_lj_->Fill(ptrel_lj,w);

	//PtRel l-bjets
	math::XYZVector bjetDir = closestBJet.momentum().Unit();
	Double_t ptrel_lbj = ROOT::Math::VectorUtil::Perp(electron->momentum(),bjetDir);
	//TVector3 bjetvec(closestBJet.momentum().X(),closestBJet.momentum().Y(),closestBJet.momentum().Z());
	//Float_t ptrel_lbj = elevec.Perp(bjetvec);
	if (deltaRBmin < 0.5) h_allele_ptrel_lbj_->Fill(ptrel_lbj,w);
      }//end loop over electrons
    }
    h_nele_->Fill(electrons->size());

  }//end if right event category

  return;  
}

//////////////////////////////////////////////////////////////////////////////////
//ZbbExtraLeptons methods
//////////////////////////////////////////////////////////////////////////////////

void ZbbExtraLeptons::book(){
  
  edm::Service<TFileService> fs;
  std::string DirName = &(*N_);
  TFileDirectory HistoVars = fs->mkdir(DirName);

  TH1F::SetDefaultSumw2(kTRUE);

  //Muons
  h_muNotFromZ_pt_         = HistoVars.make<TH1F>(N_+"_muNotFromZ_pt","p_{T} muons not coming from Z; extra #mu p_{T} (GeV/c)",50,0,100); 
  h_muNotFromZ_eta_        = HistoVars.make<TH1F>(N_+"_muNotFromZ_eta","#eta muons not coming from Z; extra #mu #eta",40,-2.5,2.5);
  h_muNotFromZ_phi_        = HistoVars.make<TH1F>(N_+"_muNotFromZ_phi","#phi muons not coming from Z; #phi",40,-TMath::Pi(),TMath::Pi());
  h_muNotFromZ_IPT_        = HistoVars.make<TH1F>(N_+"_muNotFromZ_IPT","d_{xy} muons not coming from Z; extra #mu d_{xy} (cm)",1000,-2,2);
  h_muNotFromZ_IPL_        = HistoVars.make<TH1F>(N_+"_muNotFromZ_IPL","d_{z} muons not coming from Z; extra #mu d_{z} (cm)",4000,-20,20);
  h_muNotFromZ_SIP_        = HistoVars.make<TH1F>(N_+"_muNotFromZ_SIP","muon IP/#sigma_{IP};muon IP/#sigma_{IP}",120,0.,40.); 
  h_muNotFromZ_dr_lj_      = HistoVars.make<TH1F>(N_+"_muNotFromZ_dr_lj","#Delta R_{min} muons not coming from Z / all jets; #Delta R (#mu,j)_{min}",100,0,10);
  h_muNotFromZ_dr_lbj_     = HistoVars.make<TH1F>(N_+"_muNotFromZ_dr_lbj","#Delta R_{min} muons not coming from Z / B-tagged jets; #Delta R (#mu,j_{b})_{min}",100,0,10);
  h_muNotFromZ_ptrel_lj_   = HistoVars.make<TH1F>(N_+"_muNotFromZ_ptrel_lj","p_{T}^{rel} min muons not coming from Z wrt jets; p_{T}^{rel} (GeV/c)",60,0,5);
  h_muNotFromZ_ptrel_lbj_  = HistoVars.make<TH1F>(N_+"_muNotFromZ_ptrel_lbj","p_{T}^{rel} min muons not coming from Z wrt b-jets; p_{T}^{rel} (GeV/c)",60,0,5);
  h_muNotFromZ_trackIso_   = HistoVars.make<TH1F>(N_+"_muNotFromZ_trackIso","trackIso/p_{T} muons; trackIso/p_{T}",50,0,5);
  h_muNotFromZ_caloIso_    = HistoVars.make<TH1F>(N_+"_muNotFromZ_caloIso","caloIso/p_{T} muons; caloIso/p_{T}",50,0,5);
  h_muNotFromZ_combIso_    = HistoVars.make<TH1F>(N_+"_muNotFromZ_combIso","(trackIso+caloIso)/p_{T} muons; (trackIso+caloIso)/p_{T}",100,0,10);
  h_addDiMu_mass_          = HistoVars.make<TH1F>(N_+"_addDiMu_mass","Mass of additional di-#mu candidates; m(#mu^{+}#mu^{-}) (GeV/c^{2})",60,0,120);

  //Electrons
  h_eleNotFromZ_pt_        = HistoVars.make<TH1F>(N_+"_eleNotFromZ_pt","p_{T} electrons not coming from Z; extra e p_{T} (GeV/c)",50,0,100); 
  h_eleNotFromZ_eta_       = HistoVars.make<TH1F>(N_+"_eleNotFromZ_eta","#eta electrons not coming from Z; extra e #eta",40,-2.5,2.5);
  h_eleNotFromZ_phi_       = HistoVars.make<TH1F>(N_+"_eleNotFromZ_phi","#phi elctrons not coming from Z; #phi",40,-TMath::Pi(),TMath::Pi());
  h_eleNotFromZ_IPT_       = HistoVars.make<TH1F>(N_+"_eleNotFromZ_IPT","d_{xy} electrons not coming from Z; extra e d_{xy} (cm)",1000,-2,2);
  h_eleNotFromZ_IPL_       = HistoVars.make<TH1F>(N_+"_eleNotFromZ_IPL","d_{z} electrons not coming from Z; extra e d_{z} (cm)",4000,-20,20);
  h_eleNotFromZ_SIP_       = HistoVars.make<TH1F>(N_+"_eleNotFromZ_SIP","muon IP/#sigma_{IP};muon IP/#sigma_{IP}",120,0.,40.); 
  h_eleNotFromZ_dr_lj_     = HistoVars.make<TH1F>(N_+"_eleNotFromZ_dr_lj","#Delta R_{min} electrons not coming from Z / all jets; #Delta R (e,j)_{min}",100,0,10);
  h_eleNotFromZ_dr_lbj_    = HistoVars.make<TH1F>(N_+"_eleNotFromZ_dr_lbj","#Delta R_{min} electrons not coming from Z / B-tagged jets; #Delta R (e,j_{b})_{min}",100,0,10);
  h_eleNotFromZ_ptrel_lj_  = HistoVars.make<TH1F>(N_+"_eleNotFromZ_ptrel_lj","p_{T}^{rel} min electrons not coming from Z wrt jets; p_{T}^{rel} (GeV/c)",60,0,5);
  h_eleNotFromZ_ptrel_lbj_ = HistoVars.make<TH1F>(N_+"_eleNotFromZ_ptrel_lbj","p_{T}^{rel} min electrons not coming from Z wrt b-jets; p_{T}^{rel} (GeV/c)",60,0,5);
  h_eleNotFromZ_trackIso_  = HistoVars.make<TH1F>(N_+"_eleNotFromZ_trackIso","trackIso/p_{T} electrons; trackIso/p_{T}",50,0,5);
  h_eleNotFromZ_caloIso_   = HistoVars.make<TH1F>(N_+"_eleNotFromZ_caloIso","caloIso/p_{T} muons; caloIso/p_{T}",50,0,5);
  h_eleNotFromZ_combIso_   = HistoVars.make<TH1F>(N_+"_eleNotFromZ_combIso","(trackIso+caloIso)/p_{T} electrons; (trackIso+caloIso)/p_{T}",100,0,10);
  h_addDiEle_mass_         = HistoVars.make<TH1F>(N_+"_addDiEle_mass","Mass of additional di-e candidates; m(e^{+}e^{-}) (GeV/c^{2})",60,0,120);

  //Dileptons
  h_addDiLept_mass_    = HistoVars.make<TH1F>(N_+"_addDiLept_mass","Mass of additional di leptons candidates; m(l+ l-) (GeV/c^{2})",120,-120,120);

  //Counters
  h_addmu_                 = HistoVars.make<TH1F>(N_+"_addmu","Number of additional muons in the event; n_{#mu}^{extra}",10,-0.5,9.5);
  h_addele_                = HistoVars.make<TH1F>(N_+"_addele","Number of additional electrons in the event; n_{e}^{extra}",10,-0.5,9.5);
  h_addlept_               = HistoVars.make<TH1F>(N_+"_addlept","Number of additional leptons in the event; n_{l}^{extra}",10,-0.5,9.5);
  h_addleptflav_3levts_    = HistoVars.make<TH1F>(N_+"_h_addleptflav_3levt","3 leptons events / lepton flavor",2,-0.5,1.5);
  h_addleptflav_3levts_->GetXaxis()->SetBinLabel(1,"e");
  h_addleptflav_3levts_->GetXaxis()->SetBinLabel(2,"#mu");
  h_addleptflav_4levts_ = HistoVars.make<TH1F>(N_+"_h_addleptflav_4levt","4 leptons events / lepton flavor",3,-0.5,2.5);
  h_addleptflav_4levts_->GetXaxis()->SetBinLabel(1,"ee");
  h_addleptflav_4levts_->GetXaxis()->SetBinLabel(2,"e#mu");
  h_addleptflav_4levts_->GetXaxis()->SetBinLabel(3,"#mu#mu");

}

void ZbbExtraLeptons::fillLeptons(const EventCategory& ec_mu,const EventCategory& ec_ele, edm::Handle<edm::View<pat::Muon> > muons, edm::Handle<edm::View<pat::Electron> > electrons, edm::Handle<edm::View<pat::Jet> > jets,std::string bTagAlgoWP){
  
  // vector of leptons not coming from Z
  std::vector<pat::Muon> addMu;
  std::vector<pat::Electron> addEle;
  std::vector<pat::Jet> theGoodJets;;
  std::vector<pat::Jet> theGoodBJets;

  Bool_t thecut_;   
  EventCategory ec; 
  
  if(ec_mu.isLLTightAndTriggerMatched()){
    thecut_= ZbbUtils::ProduceTheCut(ec_mu,N_,"");
    ec = ec_mu;
  } else {
    thecut_ = ZbbUtils::ProduceTheCut(ec_ele,N_,"");
    ec = ec_ele;
  }

  Double_t w = ZbbUtils::getTheWeight(ec,N_);

  // const reco::Vertex& theZVertex = Zvertexes->at(0);
  const reco::Vertex& theZVertex = ec.theZvertex_;

  ///////////
  //JETS
  //////////
  for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
    if ( ZbbUtils::isJetIdOk(*jet,"loose") && jet->pt() > lCuts_.bjetPtMin_ && jet->eta() < lCuts_.bjetEtaMax_ ){
      theGoodJets.push_back(*jet);
      if (ZbbUtils::isBJet(*jet,bTagAlgoWP))theGoodBJets.push_back(*jet);
    }
  }//end loop over jets
  

  if (thecut_){
    
    ///////////
    //MUONS
    ///////////
    if(muons->size()!=0){
      for(edm::View<pat::Muon>::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon){
	
	if (ZbbUtils::isPreselMu(*muon,theZVertex,lCuts_)){
	  if (ec.bestZcandidate_.daughter(0)->pt()!=0){
	    Double_t deltaR1   = ROOT::Math::VectorUtil::DeltaR(ec.bestZcandidate_.daughter(0)->momentum(),muon->momentum());
	    Double_t deltaR2   = ROOT::Math::VectorUtil::DeltaR(ec.bestZcandidate_.daughter(1)->momentum(),muon->momentum());
	    
	    //Veto on muons matched by deltaR with candidate's daughters
	    if(deltaR1 > 0.05 && deltaR2 > 0.05){
	      addMu.push_back(*muon);
	    }//end veto
	  }
	}//muon cut
      }//end loop over muons
    }//end if muons->size() != 0
    
    h_addmu_->Fill(addMu.size());

    //Study muons not coming from Z
    if (addMu.size()!=0) {
      for (size_t i = 0;i < addMu.size();i++){
	h_muNotFromZ_pt_->Fill(addMu[i].pt(),w);
	h_muNotFromZ_eta_->Fill(addMu[i].eta(),w);
	h_muNotFromZ_phi_->Fill(addMu[i].phi(),w);
	h_muNotFromZ_SIP_->Fill(fabs(addMu[i].dB(pat::Muon::PV3D))/addMu[i].edB(pat::Muon::PV3D),w); 
	h_muNotFromZ_trackIso_->Fill(addMu[i].trackIso()/addMu[i].pt(),w);
	h_muNotFromZ_caloIso_->Fill(addMu[i].caloIso()/addMu[i].pt(),w);
	h_muNotFromZ_combIso_->Fill((addMu[i].caloIso()+addMu[i].trackIso())/addMu[i].pt(),w);

	if ( !(addMu[i].track().isNull()) ){
	  reco::TrackRef theTrack = addMu[i].track();
	  const math::XYZPoint myVertex(theZVertex.position().x(),theZVertex.position().y(),theZVertex.position().z());
	  double IPT = theTrack->dxy(myVertex);
	  double IPL  = theTrack->dz(myVertex);
	  h_muNotFromZ_IPT_->Fill(IPT,w);
	  h_muNotFromZ_IPL_->Fill(IPL,w);
	}

	Float_t deltaRmin = 1000;
	pat::Jet closestJet;
	if(theGoodJets.size()!=0){
	  for(size_t j = 0;j < theGoodJets.size();j++){
	    Double_t deltaR   = ROOT::Math::VectorUtil::DeltaR(theGoodJets[j].momentum(),addMu[i].momentum());
	    if (deltaR < deltaRmin) {
	      closestJet = theGoodJets[j];
	      deltaRmin = deltaR;
	    }
	  }//end loop over jets
	}//end GoodJets.size() veto
	h_muNotFromZ_dr_lj_->Fill(deltaRmin,w);

	Float_t deltaRBmin = 1000;
	pat::Jet closestBJet;
	if(theGoodBJets.size()!=0){
	  for(size_t j = 0;j < theGoodBJets.size();j++){
	    Double_t deltaRB   = ROOT::Math::VectorUtil::DeltaR(theGoodBJets[j].momentum(),addMu[i].momentum());
	    if (deltaRB < deltaRBmin) {
	      closestBJet = theGoodBJets[j];
	      deltaRBmin = deltaRB;
	    }
	  }//end loop over B-jets
	}//end GoodBJets.size() veto
	h_muNotFromZ_dr_lbj_->Fill(deltaRBmin,w);	

	//PtRel l-jets
	math::XYZVector jetDir = closestJet.momentum().Unit();
	Double_t ptrel_lj = ROOT::Math::VectorUtil::Perp(addMu[i].momentum(),jetDir);
	//TVector3 muvec(addMu[i].momentum().X(),addMu[i].momentum().Y(),addMu[i].momentum().Z());
	//TVector3 jetvec(closestJet.momentum().X(),closestJet.momentum().Y(),closestJet.momentum().Z());
	//Float_t ptrel_lj = muvec.Perp(jetvec);
	if (deltaRmin < 0.5) h_muNotFromZ_ptrel_lj_->Fill(ptrel_lj,w);

	//PtRel l-bjets
	math::XYZVector bjetDir = closestBJet.momentum().Unit();
	Double_t ptrel_lbj = ROOT::Math::VectorUtil::Perp(addMu[i].momentum(),bjetDir);
	//TVector3 bjetvec(closestBJet.momentum().X(),closestBJet.momentum().Y(),closestBJet.momentum().Z());
	//Float_t ptrel_lbj = muvec.Perp(bjetvec);
	if (deltaRBmin < 0.5) h_muNotFromZ_ptrel_lbj_->Fill(ptrel_lbj,w);

	//Invariant mass of additional dilepton pairs	
	if (addMu.size()>1){
	  for (size_t j = 0;j < i;j++){
	    TLorentzVector p4Mu1(addMu[i].px(),addMu[i].py(),addMu[i].pz(),addMu[i].energy());
	    TLorentzVector p4Mu2(addMu[j].px(),addMu[j].py(),addMu[j].pz(),addMu[j].energy());
	    TLorentzVector p4sum = p4Mu1 + p4Mu2;
	    if (addMu[i].charge()*addMu[j].charge()<0){
	      h_addDiMu_mass_->Fill(p4sum.M(),w); 
	      h_addDiLept_mass_->Fill(p4sum.M(),w);  
	    }
	    else h_addDiLept_mass_->Fill(-p4sum.M(),w); 
	  }
	}
	
      }//end loop over AddMu

    }//if addMu.size()!=0
    
    /////////////
    //ELECTRONS
    /////////////
    if(electrons->size()!=0){
      for(edm::View<pat::Electron>::const_iterator electron=electrons->begin(); electron!=electrons->end(); ++electron){
	if (ZbbUtils::isPreselEle(*electron,theZVertex, lCuts_)){
	  if (ec.bestZcandidate_.daughter(0)->pt()!=0){
	    Double_t deltaR1   = ROOT::Math::VectorUtil::DeltaR(ec.bestZcandidate_.daughter(0)->momentum(),electron->momentum());
	    Double_t deltaR2   = ROOT::Math::VectorUtil::DeltaR(ec.bestZcandidate_.daughter(1)->momentum(),electron->momentum());
	    
	    //Veto on electrons matched by deltaR with candidate's daughters
	    if(deltaR1 > 0.05 && deltaR2 > 0.05){
	      addEle.push_back(*electron);
	    }//end veto
	  }
	}//electron cut
      }//end loop over electrons
    }//end if electrons->size()!=0
    
    h_addele_->Fill(addEle.size(),w);

    //Study electrons not coming from Z
    if (addEle.size()!=0) {
      for (size_t i = 0;i < addEle.size();i++){
	h_eleNotFromZ_pt_->Fill(addEle[i].pt(),w);
	h_eleNotFromZ_eta_->Fill(addEle[i].eta(),w);
	h_eleNotFromZ_phi_->Fill(addEle[i].phi(),w);
	h_eleNotFromZ_SIP_->Fill(fabs(addEle[i].dB(pat::Electron::PV3D))/addEle[i].edB(pat::Electron::PV3D),w); 
	h_eleNotFromZ_trackIso_->Fill(addEle[i].trackIso()/addEle[i].pt(),w);
	h_eleNotFromZ_caloIso_->Fill(addEle[i].caloIso()/addEle[i].pt(),w);
	h_eleNotFromZ_combIso_->Fill((addEle[i].caloIso()+addEle[i].trackIso())/addEle[i].pt(),w);

	if ( !(addEle[i].gsfTrack().isNull()) ){
	  reco::GsfTrackRef theTrack = addEle[i].gsfTrack();
	  const math::XYZPoint myVertex(theZVertex.position().x(),theZVertex.position().y(),theZVertex.position().z());
	  double IPT = theTrack->dxy(myVertex);
	  double IPL = theTrack->dz(myVertex);
	  h_eleNotFromZ_IPT_->Fill(IPT,w);
	  h_eleNotFromZ_IPL_->Fill(IPL,w);
	}

	Float_t deltaRmin = 1000;
	pat::Jet closestJet;
	if(theGoodJets.size()!=0){
	  for(size_t j = 0;j < theGoodJets.size();j++){
	    Double_t deltaR   = ROOT::Math::VectorUtil::DeltaR(theGoodJets[j].momentum(),addEle[i].momentum());
	    if (deltaR < deltaRmin) {
	      closestJet = theGoodJets[j];
	      deltaRmin = deltaR;
	    }
	  }//end loop over jets
	}//end GoodJets.size() veto
	h_eleNotFromZ_dr_lj_->Fill(deltaRmin,w);
	
	Float_t deltaRBmin = 1000;
	pat::Jet closestBJet;
	if(theGoodBJets.size()!=0){
	  for(size_t j = 0;j < theGoodBJets.size();j++){
	    Double_t deltaRB   = ROOT::Math::VectorUtil::DeltaR(theGoodBJets[j].momentum(),addEle[i].momentum());
	    if (deltaRB < deltaRBmin) {
	      closestBJet = theGoodBJets[j];
	      deltaRBmin = deltaRB;
	    }
	  }//end loop over B-jets
	}//end GoodBJets.size() veto
	h_eleNotFromZ_dr_lbj_->Fill(deltaRBmin,w);
	
	//PtRel l-jets
	math::XYZVector jetDir = closestJet.momentum().Unit();
	Double_t ptrel_lj = ROOT::Math::VectorUtil::Perp(addEle[i].momentum(),jetDir);
	//TVector3 elevec(addEle[i].momentum().X(),addEle[i].momentum().Y(),addEle[i].momentum().Z());
	//TVector3 jetvec(closestJet.momentum().X(),closestJet.momentum().Y(),closestJet.momentum().Z());
	//Float_t ptrel_lj = elevec.Perp(jetvec);
	if (deltaRmin < 0.5) h_eleNotFromZ_ptrel_lj_->Fill(ptrel_lj,w);

	//PtRel l-bjets
	math::XYZVector bjetDir = closestBJet.momentum().Unit();
	Double_t ptrel_lbj = ROOT::Math::VectorUtil::Perp(addEle[i].momentum(),bjetDir);
	//TVector3 bjetvec(closestBJet.momentum().X(),closestBJet.momentum().Y(),closestBJet.momentum().Z());
	//Float_t ptrel_lbj = elevec.Perp(bjetvec);
	if (deltaRBmin < 0.5) h_eleNotFromZ_ptrel_lbj_->Fill(ptrel_lbj,w);
	
	//Invariant mass of additional dilepton pairs	
	if (addEle.size()>1){
	  for (size_t j = 0;j < i;j++){
	    TLorentzVector p4Ele1(addEle[i].px(),addEle[i].py(),addEle[i].pz(),addEle[i].energy());
	    TLorentzVector p4Ele2(addEle[j].px(),addEle[j].py(),addEle[j].pz(),addEle[j].energy());
	    TLorentzVector p4sum = p4Ele1 + p4Ele2;
	    if (addEle[i].charge()*addEle[j].charge()<0){
	      h_addDiEle_mass_->Fill(p4sum.M(),w);  
	      h_addDiLept_mass_->Fill(p4sum.M(),w); 
	    }
	    else h_addDiLept_mass_->Fill(-p4sum.M(),w); 
	  }
	}
      }//end loop over AddEle
      
    }//if AddEle.size()!=0

    for (size_t i = 0;i < addEle.size();i++){
      for (size_t j = 0;j < addMu.size();j++){
	TLorentzVector p4Ele(addEle[i].px(),addEle[i].py(),addEle[i].pz(),addEle[i].energy());
	TLorentzVector p4Mu(addMu[j].px(),addMu[j].py(),addMu[j].pz(),addMu[j].energy());
	TLorentzVector p4sum = p4Ele + p4Mu;
	h_addDiLept_mass_->Fill(-p4sum.M(),w); 
      }
    }

    h_addlept_->Fill(addMu.size()+addEle.size());
    
    if (addMu.size()==0 && addEle.size()==1) h_addleptflav_3levts_->Fill(0);
    if (addMu.size()==1 && addEle.size()==0) h_addleptflav_3levts_->Fill(1);
    if (addMu.size()==0 && addEle.size()==2) h_addleptflav_4levts_->Fill(0);
    if (addMu.size()==1 && addEle.size()==1) h_addleptflav_4levts_->Fill(1);
    if (addEle.size()==0 && addMu.size()==2) h_addleptflav_4levts_->Fill(2);
    
  }//if right category
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////
//MCTruthLFromHF methods
//////////////////////////////////////////////////////////////////////////////////

void MCTruthLFromHF::book(){

  edm::Service<TFileService> fs;

  //Book histograms this category
  Int_t nbin_eta = 80;
  Int_t nbin_pt = 75;
  Int_t nbin_phi = 80;
  Int_t nbin_dr = 50;
  Int_t nbin_ptrel = 60;
  Double_t pt_min = 0;
  Double_t pt_max = 150;
  Double_t eta_min = -4;
  Double_t eta_max = 4;
  Double_t phi_min = - TMath::Pi();
  Double_t phi_max = TMath::Pi();
  Double_t dr_min = 0;
  Double_t dr_max = 10;
  Double_t ptrel_min = 0;
  Double_t ptrel_max = 20;

  std::string cat = &(*N_); 
  TFileDirectory MainDir = fs->mkdir("MCTruthLFromHF");  
  TFileDirectory dir = MainDir.mkdir(cat); 

  Float_t xbinpt[15] = {0,1,2,3,4,5,6,7,8,9,10,15,20,40,80};

  TH1F::SetDefaultSumw2(kTRUE);

  //Generic lepton histograms
  h_pt_ = dir.make<TH1F>(N_+"_pT","p_{T} "+N_,nbin_pt,pt_min,pt_max);
  h_eta_ = dir.make<TH1F>(N_+"_eta","#eta "+N_,nbin_eta,eta_min,eta_max);
  h_phi_ = dir.make<TH1F>(N_+"_phi","#phi "+N_,nbin_phi,phi_min,phi_max);
  h_dr_lj_ = dir.make<TH1F>(N_+"_dr_lj","Minimum #DeltaR any lepton/any jet  "+N_,nbin_dr,dr_min,dr_max);
  h_dr_lbj_ = dir.make<TH1F>(N_+"_dr_lbj","Minimum #DeltaR any lepton/any B-tagged jet "+N_,nbin_dr,dr_min,dr_max);
  h_ptrel_lj_ = dir.make<TH1F>(N_+"_ptrel_lj","p_{T}^{rel} lepton - jet min "+N_,nbin_ptrel,ptrel_min,ptrel_max);
  h_ptrel_lbj_ = dir.make<TH1F>(N_+"_ptrel_lbj","p_{T}^{rel} lepton - b-jet min "+N_,nbin_ptrel,ptrel_min,ptrel_max);
  h_IPT_ = dir.make<TH1F>(N_+"_IPT","d_{xy} "+N_,50,-2,2);  
  h_IPL_ = dir.make<TH1F>(N_+"_IPL","d_{z} "+N_,120,-60,60);
  h_CombIso_    = dir.make<TH1F>(N_+"_CombIso","trackIso/p_{T}",50,0,10);
  

  //Specific muon histograms
  h_Mu_RecoAlgo_ = dir.make<TH1F>(N_+"_Mu_RecoAlgo","Muon reconstruction algorithm"+N_,4,-0.5,3.5);  
  h_Mu_RecoAlgo_->GetXaxis()->SetBinLabel(1,"Global");
  h_Mu_RecoAlgo_->GetXaxis()->SetBinLabel(2,"Tracker");
  h_Mu_RecoAlgo_->GetXaxis()->SetBinLabel(3,"Standalone");
  h_Mu_RecoAlgo_->GetXaxis()->SetBinLabel(4,"Track.+Stand.");
  h_Mu_MuonSelector_ = dir.make<TH1F>(N_+"_Mu_MuonSelector","Muon's selector type (see twiki about muon ID)",24,-0.5,23.5);
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(1,"All");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(2,"AllGlobalMuons");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(3,"AllStandAloneMuons");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(4,"AllTrackerMuons");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(5,"TrackerMuonArbitrated");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(6,"AllArbitrated");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(7,"GlobalMuonPromptTight");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(8,"TMLastStationLoose");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(9,"TMLastStationTight");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(10,"TM2DCompatibilityLoose");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(11,"TM2DCompatibilityTight");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(12,"TMOneStationLoose");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(13,"TMOneStationTight");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(14,"TMLastStationOptimizedLowPtLoose");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(15,"TMLastStationOptimizedLowPtTight");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(16,"GMTkChiCompatibility");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(17,"GMStaChiCompatibility");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(18,"GMTkKinkTight");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(19,"TMLastStationAngLoose");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(20,"TMLastStationAngTight");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(21,"TMOneStationAngLoose");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(22,"TMOneStationAngTight");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(23,"TMLastStationOptimizedBarrelLowPtLoose");
  h_Mu_MuonSelector_->GetXaxis()->SetBinLabel(24,"TMLastStationOptimizedBarrelLowPtTight");
  h_Mu_NormChi2_  = dir.make<TH1F>(N_+"_Mu_NormChi2","Normalized Chi2 track",200,0,50);
  h_Mu_TrkHits_  = dir.make<TH1F>(N_+"_Mu_TrkHits","Number of valid tracker hits",20,-0.5,19.5);
  h_Mu_PixHits_  = dir.make<TH1F>(N_+"_Mu_PixHits","Number of valid pixel hits",10,-0.5,9.5);
  h_Mu_MuHits_  = dir.make<TH1F>(N_+"_Mu_MuHits","Number of valid muon hits",10,-0.5,9.5);
  h_Mu_NOfMatches_  = dir.make<TH1F>(N_+"_Mu_NOfMatches","Number of matches",10,-0.5,9.5);
  h_Mu_parents_ = dir.make<TH1F>(N_+"_Mu_parents","Parents of additional selected muons",3,0.5,3.5);
  h_Mu_parents_->GetXaxis()->SetBinLabel(1,"from B");
  h_Mu_parents_->GetXaxis()->SetBinLabel(2,"from C");
  h_Mu_parents_->GetXaxis()->SetBinLabel(3,"other");

  //Specific electron histograms
  h_Ele_eleID_  = dir.make<TH1F>(N_+"_Ele_eleID","Electron ID cut",8,-0.5,7.5);
  h_Ele_eleID_->GetXaxis()->SetBinLabel(1,"Fails");
  h_Ele_eleID_->GetXaxis()->SetBinLabel(2,"ele ID only");
  h_Ele_eleID_->GetXaxis()->SetBinLabel(3,"ele iso only");
  h_Ele_eleID_->GetXaxis()->SetBinLabel(4,"ele ID+iso");
  h_Ele_eleID_->GetXaxis()->SetBinLabel(5,"conv rej");
  h_Ele_eleID_->GetXaxis()->SetBinLabel(6,"conv rej + ID");
  h_Ele_eleID_->GetXaxis()->SetBinLabel(7,"conv rej + iso");
  h_Ele_eleID_->GetXaxis()->SetBinLabel(8,"all");
  h_Ele_parents_ = dir.make<TH1F>(N_+"_Ele_parents","Parents of additional selected electrons",3,0.5,3.5);
  h_Ele_parents_->GetXaxis()->SetBinLabel(1,"from B");
  h_Ele_parents_->GetXaxis()->SetBinLabel(2,"from C");
  h_Ele_parents_->GetXaxis()->SetBinLabel(3,"other");

  //GEN level histograms
  h_pt_GEN_ = dir.make<TH1F>(N_+"_p_{T}_GEN","GEN p_{T} "+N_,nbin_pt,pt_min,pt_max);
  h_pt_rebin_GEN_ = dir.make<TH1F>(N_+"_p_{T}_rebin_GEN","p_{T} "+N_,14,xbinpt);
  h_eta_GEN_ = dir.make<TH1F>(N_+"_eta_GEN","GEN #eta "+N_,nbin_eta,eta_min,eta_max);
  h_phi_GEN_ = dir.make<TH1F>(N_+"_phi_GEN","GEN #phi "+N_,nbin_phi,phi_min,phi_max);
  h_pt_GENmatchedRECO_ = dir.make<TH1F>(N_+"_p_{T}_GENmatchedRECO","GENmatchedRECO p_{T} "+N_,nbin_pt,pt_min,pt_max);
  h_pt_rebin_GENmatchedRECO_ = dir.make<TH1F>(N_+"_p_{T}_rebin_GENmatchedRECO","p_{T} "+N_,14,xbinpt);
  h_eta_GENmatchedRECO_ = dir.make<TH1F>(N_+"_eta_GENmatchedRECO","GENmatchedRECO #eta "+N_,nbin_eta,eta_min,eta_max);
  h_phi_GENmatchedRECO_ = dir.make<TH1F>(N_+"_phi_GENmatchedRECO","GENmatchedRECO #phi "+N_,nbin_phi,phi_min,phi_max);

}//end book()

void MCTruthLFromHF::fillMu(const EventCategory& ec_mu,const EventCategory& ec_ele,const pat::Muon &muon, edm::Handle<edm::View<pat::Jet> > jets, std::string bTagAlgoWP){

  EventCategory ec; 
  Bool_t thecut_;

  if(ec_mu.isLLTightAndTriggerMatched()){
    thecut_= ZbbUtils::ProduceTheCut(ec_mu,N_,"");
    ec = ec_mu;
  } else {
    thecut_ = ZbbUtils::ProduceTheCut(ec_ele,N_,"");
    ec = ec_ele;
  }
  Double_t w = ZbbUtils::getTheWeight(ec,N_);

  // const reco::Vertex& theZVertex = Zvertexes->at(0);
  const reco::Vertex& theZVertex = ec.theZvertex_;

  if (thecut_){

    Float_t pt = muon.pt();
    Float_t eta = muon.eta();
    Float_t phi = muon.phi();
  
    h_pt_->Fill(pt);
    h_eta_->Fill(eta);
    h_phi_->Fill(phi);
    
    const reco::Candidate *GenMuon = &(*muon.genParticleById(13,1));
    while (GenMuon->mother()!=0 && abs(GenMuon->mother()->pdgId()) == 13) GenMuon = GenMuon->mother();
    
    h_pt_GENmatchedRECO_->Fill(GenMuon->pt());
    h_pt_rebin_GENmatchedRECO_->Fill(GenMuon->pt());
    h_eta_GENmatchedRECO_->Fill(GenMuon->eta());
    h_phi_GENmatchedRECO_->Fill(GenMuon->phi());
    
    std::vector<pat::Jet> theGoodJets;;
    std::vector<pat::Jet> theGoodBJets;
    
    for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
      if ( ZbbUtils::isJetIdOk(*jet,"loose") && jet->pt() > lCuts_.bjetPtMin_ && jet->eta() < lCuts_.bjetEtaMax_ ){
	theGoodJets.push_back(*jet);
	if (ZbbUtils::isBJet(*jet,bTagAlgoWP))theGoodBJets.push_back(*jet);
      }
    }//end loop over jets

    //Mu parents
    std::vector<pat::Muon> addMu;
    if (ec.bestZcandidate_.daughter(0)->pt()!=0){
      Double_t deltaR1   = ROOT::Math::VectorUtil::DeltaR(ec.bestZcandidate_.daughter(0)->momentum(),muon.momentum());
      Double_t deltaR2   = ROOT::Math::VectorUtil::DeltaR(ec.bestZcandidate_.daughter(1)->momentum(),muon.momentum());
      
      //Veto on muons matched by deltaR with candidate's daughters
      if(deltaR1 > 0.05 && deltaR2 > 0.05){
	addMu.push_back(muon);
      }//end veto
    }

    if (addMu.size()!=0) {
      for (size_t i = 0;i < addMu.size();i++){
	if( addMu[i].genParticleById(13,1).isNonnull() ) {
	  const reco::Candidate *genp = &(*addMu[i].genParticleById(13,1));
	  while (genp->mother()!=0 && abs(genp->mother()->pdgId()) == 13) genp = genp->mother();
	  if (JetMCTagUtils::decayFromBHadron(*genp))h_Mu_parents_->Fill(1,w);
	  else if (JetMCTagUtils::decayFromCHadron(*genp))h_Mu_parents_->Fill(2,w);
	  else h_Mu_parents_->Fill(3,w);
	}
      }
    }
   
    if ( !(muon.track().isNull()) ){
      reco::TrackRef theTrack = muon.track();
      const math::XYZPoint myVertex(theZVertex.position().x(),theZVertex.position().y(),theZVertex.position().z());
      double IPT = theTrack->dxy(myVertex);
      double IPL  = theTrack->dz(myVertex);
      h_IPT_->Fill(IPT);
      h_IPL_->Fill(IPL);
    }
    
    Float_t deltaRmin = 1000;
    pat::Jet closestJet;
    if(theGoodJets.size()!=0){
      for(size_t j = 0;j < theGoodJets.size();j++){
	Double_t deltaR   = ROOT::Math::VectorUtil::DeltaR(theGoodJets[j].momentum(),muon.momentum());
	if (deltaR < deltaRmin) {
	  closestJet = theGoodJets[j];
	  deltaRmin = deltaR;
	}
      }//end loop over jets
    }//end GoodJets.size() veto
    h_dr_lj_->Fill(deltaRmin);
    
    Float_t deltaRBmin = 1000;
    pat::Jet closestBJet;
    if(theGoodBJets.size()!=0){
      for(size_t j = 0;j < theGoodBJets.size();j++){
	Double_t deltaRB   = ROOT::Math::VectorUtil::DeltaR(theGoodBJets[j].momentum(),muon.momentum());
	if (deltaRB < deltaRBmin) {
	  closestBJet = theGoodBJets[j];
	  deltaRBmin = deltaRB;
	}
      }//end loop over B-jets
    }//end GoodBJets.size() veto
    h_dr_lbj_->Fill(deltaRBmin);    
    
    //PtRel l-jets
    math::XYZVector jetDir = closestJet.momentum().Unit();
    Double_t ptrel_lj = ROOT::Math::VectorUtil::Perp(muon.momentum(),jetDir);
    //TVector3 muvec(muon.momentum().X(),muon.momentum().Y(),muon.momentum().Z());
    //TVector3 jetvec(closestJet.momentum().X(),closestJet.momentum().Y(),closestJet.momentum().Z());
    //Float_t ptrel_lj = muvec.Perp(jetvec);
    h_ptrel_lj_->Fill(ptrel_lj);
    
    //PtRel l-bjets
    math::XYZVector bjetDir = closestBJet.momentum().Unit();
    Double_t ptrel_lbj = ROOT::Math::VectorUtil::Perp(muon.momentum(),bjetDir);
    //TVector3 bjetvec(closestBJet.momentum().X(),closestBJet.momentum().Y(),closestBJet.momentum().Z());
    //Float_t ptrel_lbj = muvec.Perp(bjetvec);
    h_ptrel_lbj_->Fill(ptrel_lbj);
    
    h_CombIso_->Fill( (muon.trackIso()+muon.caloIso())/muon.pt() );
  
    //Specific muon histograms
    
    if ( muon.isGlobalMuon() ){
      h_Mu_RecoAlgo_->Fill(0);
      h_Mu_NormChi2_->Fill( muon.normChi2() );
      h_Mu_TrkHits_->Fill( muon.innerTrack()->hitPattern().numberOfValidTrackerHits() );
      h_Mu_PixHits_->Fill( muon.innerTrack()->hitPattern().numberOfValidPixelHits() );
      h_Mu_MuHits_->Fill( muon.globalTrack()->hitPattern().numberOfValidMuonHits() );
      h_Mu_NOfMatches_->Fill( muon.numberOfMatches() );
    }//is global
    else {
      if ( muon.isTrackerMuon() && !(muon.isStandAloneMuon()) ){
	h_Mu_RecoAlgo_->Fill(1);
	h_Mu_TrkHits_->Fill( muon.innerTrack()->hitPattern().numberOfValidTrackerHits() );
	h_Mu_PixHits_->Fill( muon.innerTrack()->hitPattern().numberOfValidPixelHits() );
      }
      if ( !(muon.isTrackerMuon()) && muon.isStandAloneMuon() )h_Mu_RecoAlgo_->Fill(2);
      if ( muon.isTrackerMuon() && muon.isStandAloneMuon() ){
	h_Mu_RecoAlgo_->Fill(3);
	h_Mu_TrkHits_->Fill( muon.innerTrack()->hitPattern().numberOfValidTrackerHits() );
	h_Mu_PixHits_->Fill( muon.innerTrack()->hitPattern().numberOfValidPixelHits() );
      }
      h_Mu_NOfMatches_->Fill( muon.numberOfMatches() );
    }//is not global
    
    if ( muon.isGood("All") )h_Mu_MuonSelector_->Fill( 0 );
    if ( muon.isGood("AllGlobalMuons") )h_Mu_MuonSelector_->Fill( 1 );
    if ( muon.isGood("AllStandAloneMuons") )h_Mu_MuonSelector_->Fill( 2 );
    if ( muon.isGood("AllTrackerMuons") )h_Mu_MuonSelector_->Fill( 3 );
    if ( muon.isGood("TrackerMuonArbitrated") )h_Mu_MuonSelector_->Fill( 4 );
    if ( muon.isGood("AllArbitrated") )h_Mu_MuonSelector_->Fill( 5 );
    if ( muon.isGood("GlobalMuonPromptTight") )h_Mu_MuonSelector_->Fill( 6 );
    if ( muon.isGood("TMLastStationLoose") )h_Mu_MuonSelector_->Fill( 7 );
    if ( muon.isGood("TMLastStationTight") )h_Mu_MuonSelector_->Fill( 8 );
    if ( muon.isGood("TM2DCompatibilityLoose") )h_Mu_MuonSelector_->Fill( 9 );
    if ( muon.isGood("TM2DCompatibilityTight") )h_Mu_MuonSelector_->Fill( 10 );
    if ( muon.isGood("TMOneStationLoose") )h_Mu_MuonSelector_->Fill( 11 );
    if ( muon.isGood("TMOneStationTight") )h_Mu_MuonSelector_->Fill( 12 );
    if ( muon.isGood("TMLastStationOptimizedLowPtLoose") )h_Mu_MuonSelector_->Fill( 13 );
    if ( muon.isGood("TMLastStationOptimizedLowPtTight") )h_Mu_MuonSelector_->Fill( 14 );
    if ( muon.isGood("GMTkChiCompatibility") )h_Mu_MuonSelector_->Fill( 15 );
    if ( muon.isGood("GMStaChiCompatibility") )h_Mu_MuonSelector_->Fill( 16 );
    if ( muon.isGood("GMTkKinkTight") )h_Mu_MuonSelector_->Fill( 17 );
    if ( muon.isGood("TMLastStationAngLoose") )h_Mu_MuonSelector_->Fill( 18 );
    if ( muon.isGood("TMLastStationAngTight") )h_Mu_MuonSelector_->Fill( 19 );
    if ( muon.isGood("TMOneStationAngLoose") )h_Mu_MuonSelector_->Fill( 20 );
    if ( muon.isGood("TMOneStationAngTight") )h_Mu_MuonSelector_->Fill( 21 );
    if ( muon.isGood("TMLastStationOptimizedBarrelLowPtLoose") )h_Mu_MuonSelector_->Fill( 22 );
    if ( muon.isGood("TMLastStationOptimizedBarrelLowPtTight") )h_Mu_MuonSelector_->Fill( 23 );
    
  }//end if (thecut_)
  
  return;
}//end fillMu


void MCTruthLFromHF::fillEle(const EventCategory& ec_mu,const EventCategory& ec_ele,const pat::Electron &electron, edm::Handle<edm::View<pat::Jet> > jets,std::string bTagAlgoWP){

  //cout<<"debug0"<<endl;
  EventCategory ec; 
  Bool_t thecut_;
 
  if(ec_mu.isLLTightAndTriggerMatched()){
    thecut_= ZbbUtils::ProduceTheCut(ec_mu,N_,"");
    ec = ec_mu;
  } else {
    thecut_ = ZbbUtils::ProduceTheCut(ec_ele,N_,"");
    ec = ec_ele;
  }

  // const reco::Vertex& theZVertex = Zvertexes->at(0);
  const reco::Vertex& theZVertex = ec.theZvertex_;

  Double_t w = ZbbUtils::getTheWeight(ec,N_);
  
  if (thecut_){
    
    Float_t pt = electron.pt();
    Float_t eta = electron.eta();
    Float_t phi = electron.phi();
    
    h_pt_->Fill(pt);
    h_eta_->Fill(eta);
    h_phi_->Fill(phi);

    const reco::Candidate *GenElectron = &(*electron.genParticleById(11,1));
    while (GenElectron->mother()!=0 && abs(GenElectron->mother()->pdgId()) == 11) GenElectron = GenElectron->mother();
    
    h_pt_GENmatchedRECO_->Fill(GenElectron->pt());
    h_pt_rebin_GENmatchedRECO_->Fill(GenElectron->pt());
    h_eta_GENmatchedRECO_->Fill(GenElectron->eta());
    h_phi_GENmatchedRECO_->Fill(GenElectron->phi());
    
    std::vector<pat::Jet> theGoodJets;;
    std::vector<pat::Jet> theGoodBJets;
    
    for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
      if ( ZbbUtils::isJetIdOk(*jet,"loose") && jet->pt() > lCuts_.bjetPtMin_ && jet->eta() < lCuts_.bjetEtaMax_ ){
	theGoodJets.push_back(*jet);
	if (ZbbUtils::isBJet(*jet,bTagAlgoWP))theGoodBJets.push_back(*jet);
      }
    }//end loop over jets

    //Ele parents
    std::vector<pat::Electron> addEle;
    if (ec.bestZcandidate_.daughter(0)->pt()!=0){
      Double_t deltaR1   = ROOT::Math::VectorUtil::DeltaR(ec.bestZcandidate_.daughter(0)->momentum(),electron.momentum());
      Double_t deltaR2   = ROOT::Math::VectorUtil::DeltaR(ec.bestZcandidate_.daughter(1)->momentum(),electron.momentum());
      
      //Veto on electrons matched by deltaR with candidate's daughters
      if(deltaR1 > 0.05 && deltaR2 > 0.05){
	addEle.push_back(electron);
      }//end veto
    }

    if (addEle.size()!=0) {
      for (size_t i = 0;i < addEle.size();i++){
	if( addEle[i].genParticleById(11,1).isNonnull() ) {
	  const reco::Candidate *genp = &(*addEle[i].genParticleById(11,1));
	  while (genp->mother()!=0 && abs(genp->mother()->pdgId()) == 11) genp = genp->mother();
	  if (JetMCTagUtils::decayFromBHadron(*genp))h_Ele_parents_->Fill(1,w);
	  else if (JetMCTagUtils::decayFromCHadron(*genp))h_Ele_parents_->Fill(2,w);
	  else h_Ele_parents_->Fill(3,w);
	}
      }
    }

    if ( !(electron.track().isNull()) ){
      reco::TrackRef theTrack = electron.track();
      const math::XYZPoint myVertex(theZVertex.position().x(),theZVertex.position().y(),theZVertex.position().z());
      double IPT = theTrack->dxy(myVertex);
      double IPL  = theTrack->dz(myVertex);
      h_IPT_->Fill(IPT);
      h_IPL_->Fill(IPL);
    }
   
    Float_t deltaRmin = 1000;
    pat::Jet closestJet;
    if(theGoodJets.size()!=0){
      for(size_t j = 0;j < theGoodJets.size();j++){
	Double_t deltaR   = ROOT::Math::VectorUtil::DeltaR(theGoodJets[j].momentum(),electron.momentum());
	if (deltaR < deltaRmin) {
	  closestJet = theGoodJets[j];
	  deltaRmin = deltaR;
	}
      }//end loop over jets
    }//end GoodJets.size() veto
    h_dr_lj_->Fill(deltaRmin);
    
    Float_t deltaRBmin = 1000;
    pat::Jet closestBJet;
    if(theGoodBJets.size()!=0){
      for(size_t j = 0;j < theGoodBJets.size();j++){
	Double_t deltaRB   = ROOT::Math::VectorUtil::DeltaR(theGoodBJets[j].momentum(),electron.momentum());
	if (deltaRB < deltaRBmin) {
	  closestBJet = theGoodBJets[j];
	  deltaRBmin = deltaRB;
	}
      }//end loop over B-jets
    }//end GoodBJets.size() veto
    h_dr_lbj_->Fill(deltaRBmin);    
    
    //PtRel l-jets
    math::XYZVector jetDir = closestJet.momentum().Unit();
    Double_t ptrel_lj = ROOT::Math::VectorUtil::Perp(electron.momentum(),jetDir);
    //TVector3 elevec(electron.momentum().X(),electron.momentum().Y(),electron.momentum().Z());
    //TVector3 jetvec(closestJet.momentum().X(),closestJet.momentum().Y(),closestJet.momentum().Z());
    //Float_t ptrel_lj = elevec.Perp(jetvec);
    h_ptrel_lj_->Fill(ptrel_lj);
    
    //PtRel l-bjets
    math::XYZVector bjetDir = closestBJet.momentum().Unit();
    Double_t ptrel_lbj = ROOT::Math::VectorUtil::Perp(electron.momentum(),bjetDir);
    //TVector3 bjetvec(closestBJet.momentum().X(),closestBJet.momentum().Y(),closestBJet.momentum().Z());
    //Float_t ptrel_lbj = elevec.Perp(bjetvec);
    h_ptrel_lbj_->Fill(ptrel_lbj);
    
    h_CombIso_->Fill( (electron.trackIso()+electron.caloIso())/electron.pt() );
    
    //Specific electron histograms
    if ( electron.electronID("eidVBTFRel95")==0 )h_Ele_eleID_->Fill(0);
    if ( electron.electronID("eidVBTFRel95")==1 )h_Ele_eleID_->Fill(1);
    if ( electron.electronID("eidVBTFRel95")==2 )h_Ele_eleID_->Fill(2);
    if ( electron.electronID("eidVBTFRel95")==3 )h_Ele_eleID_->Fill(3);
    if ( electron.electronID("eidVBTFRel95")==4 )h_Ele_eleID_->Fill(4);
    if ( electron.electronID("eidVBTFRel95")==5 )h_Ele_eleID_->Fill(5);
    if ( electron.electronID("eidVBTFRel95")==6 )h_Ele_eleID_->Fill(6);
    if ( electron.electronID("eidVBTFRel95")==7 )h_Ele_eleID_->Fill(7);
    
  }//end if (thecut_)

  return;

}//end fillEle


void MCTruthLFromHF::fillGenParticle(const reco::Candidate &genp){

  h_pt_GEN_->Fill(genp.pt());
  h_pt_rebin_GEN_->Fill(genp.pt());
  h_eta_GEN_->Fill(genp.eta());
  h_phi_GEN_->Fill(genp.phi());

}//end fillGenParticle()
