#include <map>
#include <string>

#include "TH1.h"
#include <Math/VectorUtil.h>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

class PatBasicAnalyzer : public edm::EDAnalyzer {

public:
  /// default constructor
  explicit PatBasicAnalyzer(const edm::ParameterSet&);
  /// default destructor
  ~PatBasicAnalyzer();
  
private:
  /// everything that needs to be done before the event loop
  virtual void beginJob() ;
  /// everything that needs to be done during the event loop
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  /// check if histogram was booked
  bool booked(const std::string histName) const { return histContainer_.find(histName.c_str())!=histContainer_.end(); };
  /// fill histogram if it had been booked before
  void fill(const std::string histName, double value) const { if(booked(histName.c_str())) histContainer_.find(histName.c_str())->second->Fill(value); };
  /// everything that needs to be done after the event loop
  virtual void endJob() ;
  
  // simple map to contain all histograms; 
  // histograms are booked in the beginJob() 
  // method
  std::map<std::string,TH1F*> histContainer_; 
  std::vector<TH1F*> h_ptJetsCorrected;

  // input tags  
  edm::InputTag photonSrc_;
  edm::InputTag elecSrc_;
  edm::InputTag muonSrc_;
  edm::InputTag tauSrc_;
  edm::InputTag jetSrc_;
  edm::InputTag metSrc_;

  /// correction levels for pat jet
  std::vector<std::string> str_corrLevels;
  std::vector<std::string> corrections_;
  Int_t corrsize;

};

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

PatBasicAnalyzer::PatBasicAnalyzer(const edm::ParameterSet& iConfig):
  histContainer_(),
  //photonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc")),
  //tauSrc_(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc" )),
  //metSrc_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc" )),                 // TO BE FIXED
  elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc")),
  muonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc")),
  jetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc" ))
{
  // read and parse corrLevels
  std::vector<std::string> str_corrLevels = iConfig.getParameter<std::vector<std::string> >("corrLevels");
  corrections_.reserve(str_corrLevels.size());
  for (std::vector<std::string>::const_iterator it = str_corrLevels.begin(), ed = str_corrLevels.end(); it != ed; ++it) {
    corrections_.push_back(*it);
  }
  corrsize=corrections_.size();
}

PatBasicAnalyzer::~PatBasicAnalyzer()
{
}

void
PatBasicAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // get electron collection
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByLabel(elecSrc_,electrons);

  // get muon collection
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(muonSrc_,muons);

  // get jet collection
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jetSrc_,jets);

  // get tau collection  
  // edm::Handle<edm::View<pat::Tau> > taus;
  //iEvent.getByLabel(tauSrc_,taus);
  
  // get met collection  
  //edm::Handle<edm::View<pat::MET> > mets;
  //iEvent.getByLabel(metSrc_,mets);
  
  // get photon collection  
  //edm::Handle<edm::View<pat::Photon> > photons;
  //iEvent.getByLabel(photonSrc_,photons);
  
  //std::cout<<"before muon plots"<<std::endl;

//      __  __                      _     _      
//     |  \/  |_  _ ___ _ _    _ __| |___| |_ ___
//     | |\/| | || / _ \ ' \  | '_ \ / _ \  _(_-<
//     |_|  |_|\_,_\___/_||_| | .__/_\___/\__/__/
//                            |_|                

  // loop over muons
  for (edm::View<pat::Muon>::const_iterator muon = muons->begin();  muon != muons->end(); ++muon){
    fill("muon_pt",muon->pt());  
    fill("muon_eta",muon->eta());
    fill("muon_phi",muon->phi());
    fill("muon_dB",muon->dB());
    double muIP      = fabs(muon->dB(pat::Muon::PV3D));
    double muIPError = muon->edB(pat::Muon::PV3D); 
    double SIP =  muIP/muIPError;
    fill("muon_SIP",SIP);
    if(muon->isGlobalMuon()==true || muon->isTrackerMuon()==true){
      fill("muon_iso",(muon->trackIso() + muon->caloIso())/muon->pt());
      fill("isGlobalMuon",muon->isGlobalMuon()); 
      fill("isTrackerMuon",muon->isTrackerMuon()); 
      fill("muon_trackerhits",muon->innerTrack()->hitPattern().numberOfValidTrackerHits());
      fill("muon_pixelhits",muon->innerTrack()->hitPattern().numberOfValidPixelHits());
      if(muon->isGlobalMuon()==true){
	fill("muon_chi2",muon->globalTrack()->normalizedChi2());
	fill("muon_muonhits",muon->globalTrack()->hitPattern().numberOfValidMuonHits());
	fill("muon_numberOfMatches",muon->numberOfMatches());
      }
    }
  }
  
  //std::cout<<"after muon plots"<<std::endl;

  edm::View<pat::Muon>::const_iterator muon1,muon2;
  reco::Candidate::LorentzVector p1,p2,sum;

  //Fill histogram AllDiMuMass
  if(muons->size()>1){
    
    for (muon1 = muons->begin(); muon1 != muons->end(); muon1++){
      for (muon2 = muon1+1; muon2 != muons->end(); muon2++) {
	p1 = muon1->p4();
	p2 = muon2->p4();  
	sum = p1 + p2;
	if((muon1->charge()*muon2->charge())<0){ 
	  fill("AllDiMuMass",sum.mass());
	  fill("AllDiMuPt",sum.pt());
	  fill("AllDiMuEta",sum.eta());
	  fill("AllDiMuPhi",sum.phi());

	  double deltaR = ROOT::Math::VectorUtil::DeltaR(muon1->momentum(),muon2->momentum());
	  double deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(muon1->momentum(),muon2->momentum());
	  fill("AllDiMuDeltaPhi",deltaPhi);
	  fill("AllDiMuDeltaR",deltaR);
	}
      }						 
    }
  }
  
  //Fill histogram PreselMuMass
  if(muons->size()>1){
    
    for (muon1 = muons->begin(); muon1 != muons->end(); muon1++){
      
      if ( muon1->pt() > 10. && 
	   muon1->isGlobalMuon()==true && muon1->isTrackerMuon()==true && 
	   muon1->normChi2() < 15 && muon1->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 && muon1->innerTrack()->hitPattern().numberOfValidPixelHits() > 0  &&
	   muon1->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && 
	   muon1->dB() < 0.2 && (muon1->trackIso()+muon1->caloIso()) < 0.15*muon1->pt() && muon1->numberOfMatches() > 1 && abs(muon1->eta()) < 2.1){
	
	for (muon2 = muon1+1; muon2 != muons->end(); muon2++) {
	
	if ( muon2->pt() > 10. && 
	     muon2->isGlobalMuon()==true && muon2->isTrackerMuon()==true && 
	     muon2->normChi2() < 15 && muon2->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 && muon2->innerTrack()->hitPattern().numberOfValidPixelHits() > 0  &&
	     muon2->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && 
	     muon2->dB() < 0.2 && (muon2->trackIso()+muon2->caloIso()) < 0.15*muon2->pt() && muon2->numberOfMatches() > 1 && abs(muon2->eta()) < 2.1){
	  p1 = muon1->p4();
	  p2 = muon2->p4();  
	  sum = p1 + p2;
	  if((muon1->charge()*muon2->charge())<0){ 
	    fill("PreselMuMass",sum.mass());
	    fill("PreselDiMuPt",sum.pt());
	    fill("PreselDiMuEta",sum.eta());
	    fill("PreselDiMuPhi",sum.phi());
	    
	    // check the Delta R between the two muons
	    // (Delta_R^2 = Delta_Eta^2 + Delta_Phi^2)
	    double deltaR = ROOT::Math::VectorUtil::DeltaR(muon1->momentum(),muon2->momentum());
	    double deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(muon1->momentum(),muon2->momentum());
	    fill("PreselDiMuDeltaPhi",deltaPhi);
	    fill("PreselDiMuDeltaR",deltaR);
	  }
	}
      }
      }
    }
  }
 
  //std::cout<<"after muon mass plots"<<std::endl;

//      _     _          _     _      
//   _ | |___| |_   _ __| |___| |_ ___
//  | || / -_)  _| | '_ \ / _ \  _(_-<
//   \__/\___|\__| | .__/_\___/\__/__/
//                 |_|                

  // loop over jets
  size_t nJets=0;
  for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
    // fill basic kinematics for all jets
    fill( "jet_pt" , jet->pt());
    for(Int_t j = 0; j <  corrsize; j++) {
      h_ptJetsCorrected[j]->Fill(jet->correctedJet(corrections_[j]).pt());
    }
    fill( "jet_eta", jet->eta());
    fill( "jet_phi", jet->phi());
    fill( "jet_emf", jet->emEnergyFraction() );
    double discrTC = jet->bDiscriminator("trackCountingHighEffBJetTags");
    double discrSSVHP = jet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
    double discrSSVHE = jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    double discrCSV = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
    fill("jet_discrTC", discrTC);
    fill("jet_discrSSVHP", discrSSVHP);
    fill("jet_discrSSVHE", discrSSVHE);
    fill("jet_discrCSV", discrCSV);
  
    // multiplicity of pt>15 jets
    if(jet->pt()>15){ 
      ++nJets;
    } 
  }

  //std::cout<<"after jet loop"<<std::endl;

  if(jets->size()>0){
    // fill basic kinematics for leading jets
    fill( "leadingjet_pt" , (*jets)[0].pt());
    fill( "leadingjet_eta", (*jets)[0].eta());
    fill( "leadingjet_phi", (*jets)[0].phi());
    
    //std::cout<<"after leading jet plots"<<std::endl;
    
    // fill basic kinematics for sub-leading jets (if any)
    if(jets->size()>1){
      fill( "di_jet_mass", ((*jets)[0].p4()+(*jets)[1].p4()).mass());  
      fill( "subleadingjet_pt" , (*jets)[1].pt());
      fill( "subleadingjet_eta", (*jets)[1].eta());
      fill( "subleadingjet_phi", (*jets)[1].phi());
    }
  }
  
  //std::cout<<"after jet plots"<<std::endl;

  // do something similar for the other candidates
  //histContainer_["photons"]->Fill(photons->size() );
  //histContainer_["taus" ]->Fill(taus->size()  );
  //histContainer_["met"  ]->Fill(mets->empty() ? 0 : (*mets)[0].et());             // TO BE FIXED!!!
  histContainer_["jets"]->Fill(nJets);
  histContainer_["elecs" ]->Fill(electrons->size());
  histContainer_["muons"]->Fill(muons->size() );

  //std::cout<<"after collection plots"<<std::endl;

//   ___ _        _                       _     _      
//  | __| |___ __| |_ _ _ ___ _ _    _ __| |___| |_ ___
//  | _|| / -_) _|  _| '_/ _ \ ' \  | '_ \ / _ \  _(_-<
//  |___|_\___\__|\__|_| \___/_||_| | .__/_\___/\__/__/
//                                  |_|                

  // loop over electrons
  for(edm::View<pat::Electron>::const_iterator elec=electrons->begin(); elec !=electrons->end(); ++elec){
    // fill simple histograms
    histContainer_["electron_pt"]->Fill(elec->pt());
    histContainer_["electron_eta"]->Fill(elec->eta());
    histContainer_["electron_phi"]->Fill(elec->phi());
    histContainer_["electron_iso"]->Fill((elec->trackIso()+elec->caloIso())/elec->pt() );
    histContainer_["electron_eop"]->Fill(elec->eSeedClusterOverP() );
    histContainer_["electron_clus"]->Fill(elec->e1x5()/elec->e5x5() );
    //fill electron id histograms
    if( elec->electronID("eidVBTFCom80") > 0.5 )
      histContainer_["electron_eIDs" ]->Fill(0);
    if( elec->electronID("eidVBTFCom95") > 0.5 )
      histContainer_["electron_eIDs" ]->Fill(1);
    if( elec->electronID("eidVBTFRel80") > 0.5 )
      histContainer_["electron_eIDs" ]->Fill(2);
    if( elec->electronID("eidVBTFRel95") > 0.5 )
      histContainer_["electron_eIDs" ]->Fill(3);
  }
  
  //std::cout<<"after electron plots"<<std::endl;

  edm::View<pat::Electron>::const_iterator electron1,electron2;
  reco::Candidate::LorentzVector ep1,ep2,esum;

  //Fill histogram AllDiEleMass
  if(electrons->size()>1){
    
    for (electron1 = electrons->begin(); electron1 != electrons->end(); electron1++){
      for (electron2 = electron1+1; electron2 != electrons->end(); electron2++) {
	ep1 = electron1->p4();
	ep2 = electron2->p4();  
	esum = ep1 + ep2;
	if((electron1->charge()*electron2->charge())<0){
	  fill("AllDiEleMass",esum.mass());
	  fill("AllDiElePt",esum.pt());
	  fill("AllDiEleEta",esum.eta());
	  fill("AllDiElePhi",esum.phi());

	  double deltaR = ROOT::Math::VectorUtil::DeltaR(electron1->momentum(),electron2->momentum());
	  double deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(electron1->momentum(),electron2->momentum());
	  fill("AllDiEleDeltaPhi",deltaPhi);
	  fill("AllDiEleDeltaR",deltaR);
	}
      }						 
    }
  }
  
  //Fill histogram PreselEleMass
  if(electrons->size()>1){
    
    for (electron1 = electrons->begin(); electron1 != electrons->end(); electron1++){
      
      if ( electron1->pt() > 10.0 && abs(electron1->eta()) < 2.5 && (electron1->isEE() || electron1->isEB()) && !electron1->isEBEEGap() && electron1->electronID("eidVBTFRel80") == 7 ){
	
	for (electron2 = electron1+1; electron2 != electrons->end(); electron2++) {
	
	  if (  electron2->pt() > 10.0 && abs(electron2->eta()) < 2.5 && (electron2->isEE() || electron2->isEB()) && !electron2->isEBEEGap() && electron2->electronID("eidVBTFRel80") == 7 ){
	    ep1 = electron1->p4();
	    ep2 = electron2->p4();  
	    esum = ep1 + ep2;
	    if((electron1->charge()*electron2->charge())<0){ 
	      fill("PreselEleMass",esum.mass());
	      fill("PreselDiElePt",esum.pt());
	      fill("PreselDiEleEta",esum.eta());
	      fill("PreselDiElePhi",esum.phi());

	      double deltaR = ROOT::Math::VectorUtil::DeltaR(electron1->momentum(),electron2->momentum());
	      double deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(electron1->momentum(),electron2->momentum());
	      fill("PreselDiEleDeltaPhi",deltaPhi);
	      fill("PreselDiEleDeltaR",deltaR);
	    }
	  }
	}
      }
    }
  }

  //std::cout<<"after electron mass plots"<<std::endl;

}

void 
PatBasicAnalyzer::beginJob()
{
  // register to the TFileService
  edm::Service<TFileService> fs;
  
  // book histograms:

  // General Object Multiplicity
  // histContainer_["photons"]=fs->make<TH1F>("photons","photon multiplicity",10, 0,10);
  // histContainer_["taus"]=fs->make<TH1F>("taus","tau multiplicity",10,0, 10);
  // histContainer_["met"]=fs->make<TH1F>("met","missing E_{T}",20, 0, 100);                             // TO BE FIXED!!!!
  histContainer_["elecs"  ]=fs->make<TH1F>("n_elecs","electron multiplicity; n_{electrons};events",15,-0.5,14.5);
  histContainer_["muons"  ]=fs->make<TH1F>("n_muons","muon multiplicity; n_{muons};events",15,-0.5,14.5);
  histContainer_["jets"]=fs->make<TH1F>("n_jets","jet multiplicity; n_{jets};events",15,-0.5,14.5);
  
  // electron variables (for all electrons)
  histContainer_["electron_pt"  ]=fs->make<TH1F>("electron_pt","electron pt;p_{T} [GeV]",150,0.,150.);
  histContainer_["electron_eta" ]=fs->make<TH1F>("electron_eta","electron eta; electron #eta",30,-3.,3.);
  histContainer_["electron_phi" ]=fs->make<TH1F>("electron_phi","electron phi; electron #phi [rad];",60,-TMath::Pi(),TMath::Pi());
  histContainer_["electron_iso" ]=fs->make<TH1F>("electron_iso","electron iso;(#sum_{trk}p_{T} + #sum_{cal}E_{T})/p_{T} (combRelIso);",60,0.,20.);
  histContainer_["electron_eop" ]=fs->make<TH1F>("electron_eop","electron eop; electron E/p;",40,0.,1.);
  histContainer_["electron_clus"]=fs->make<TH1F>("electron_clus","electron clus; electron E_{1x5}/E_{5x5}",40,0.,1.);
  histContainer_["electron_eIDs"]=fs->make<TH1F>("electron_eIDs","electron eIDS; electron eleID;",4,0.,4.);
  histContainer_["AllDiEleMass"] = fs->make<TH1F>("AllDiEleMass","Invariant mass spectrum of all di-electrons; M(e^{+}e^{-}) [GeV]",100,0,120);
  histContainer_["AllDiElePt"] = fs->make<TH1F>("AllDiElePt","p_{T} of all di-electrons; p_{T}(e^{+}e^{-}) [GeV]",30,0,150);
  histContainer_["AllDiEleEta"] = fs->make<TH1F>("AllDiEleEta","#eta of all di-electrons; #eta(e^{+}e^{-})",50,-5,5);
  histContainer_["AllDiElePhi"] = fs->make<TH1F>("AllDiElePhi","#phi of all di-electrons; #phi(e^{+}e^{-}) [rad]",30,-TMath::Pi(),TMath::Pi());  
  histContainer_["AllDiEleDeltaPhi"] = fs->make<TH1F>("AllDiEleDeltaPhi","#Delta #phi of all di-electrons; #Delta#phi(e^{+}e^{-}) [rad]",60,0.,TMath::Pi());  
  histContainer_["AllDiEleDeltaR"] = fs->make<TH1F>("AllDiEleDeltaR","#Delta R of all di-electrons; #DeltaR(e^{+}e^{-}) [rad]",100,0,10); 
  histContainer_["PreselEleMass"] = fs->make<TH1F>("PreselEleMass","Invariant mass spectrum of Presel di-electrons; M(e^{+}e^{-}) [GeV]",100,0,120);
  histContainer_["PreselDiElePt"] = fs->make<TH1F>("PreselDiElePt","p_{T} of Presel di-electrons; p_{T}(e^{+}e^{-}) [GeV]",30,0,150);
  histContainer_["PreselDiEleEta"] = fs->make<TH1F>("PreselDiEleEta","#eta of Presel di-electrons; #eta(e^{+}e^{-})",50,-5,5);
  histContainer_["PreselDiElePhi"] = fs->make<TH1F>("PreselDiElePhi","#phi of Presel di-electrons; #phi(e^{+}e^{-}) [rad]",30,-TMath::Pi(),TMath::Pi());  
  histContainer_["PreselDiEleDeltaPhi"] = fs->make<TH1F>("PreselDiEleDeltaPhi","#Delta #phi of Presel di-electrons; #Delta#phi(e^{+}e^{-}) [rad]",60,0.,TMath::Pi());  
  histContainer_["PreselDiEleDeltaR"] = fs->make<TH1F>("PreselDiEleDeltaR","#Delta R of Presel di-electrons; #DeltaR(e^{+}e^{-}) [rad]",100,0,10);  

  // jet variables (for all jets)
  histContainer_["jet_pt"]=fs->make<TH1F>("alljet_pt","p_{T}(Jet) all jets;p_{T}(Jet) [GeV]",60,0.,300.);
  histContainer_["jet_eta"]=fs->make<TH1F>("alljet_eta","#eta (Jet) all jets;#eta (Jet)",60,-3.,3.);
  histContainer_["jet_phi"]=fs->make<TH1F>("alljet_phi","#phi (Jet) all jets;#phi (Jet)",60,-TMath::Pi(),TMath::Pi());
  histContainer_["jet_emf"]= fs->make<TH1F>("emf","emf;jet EM fraction",40,0.,1.);
 
  // jet pt corrected

  for(Int_t i = 0; i < corrsize ; i++) {
    h_ptJetsCorrected.push_back(fs->make<TH1F>(("jet_pt"+corrections_[i]).c_str(),("p_{T}(Jet);p_{T} (Jet "+corrections_[i]+ ")[GeV]").c_str(),60,0.,300.));
  }
 
  // leading jet variables 
  histContainer_["leadingjet_pt"]=fs->make<TH1F>("leadingjet_pt","p_{T}(Jet);p_{T}(Leading Jet) [GeV]",60,0.,300.);
  histContainer_["leadingjet_eta"]=fs->make<TH1F>("leadingjet_eta","#eta (Jet);#eta (Leading Jet)",60,-3.,3.);
  histContainer_["leadingjet_phi"]=fs->make<TH1F>("leadingjet_phi","#phi (Jet);#phi (Leading Jet)",60,-TMath::Pi(),TMath::Pi());
  // leading jet variables (if any)
  histContainer_["subleadingjet_pt"]=fs->make<TH1F>("subleadingjet_pt","p_{T}(Jet);p_{T}(Subleading Jet) [GeV]",60,0.,300.);
  histContainer_["subleadingjet_eta"]=fs->make<TH1F>("subleadingjet_eta","#eta (Jet);#eta (Subleading Jet)",60,-3.,3.);
  histContainer_["subleadingjet_phi"]=fs->make<TH1F>("subleadingjet_phi","#phi (Jet);#phi (Subleading Jet)",60,-TMath::Pi(),TMath::Pi());
  histContainer_["jet_discrTC"]=fs->make<TH1F>("jet_discrTC", "TC discriminant;TC discriminant",100,0.,20.); 
  histContainer_["jet_discrSSVHP"]=fs->make<TH1F>("jet_discrSSVHP","SSVHP discriminant; SSVHP discriminant",100,-2.,10.);
  histContainer_["jet_discrSSVHE"]=fs->make<TH1F>("jet_discrSSVHE","SSVHE discriminant;SSVHE discriminant",100,-2.,10.);
  histContainer_["jet_discrCSV"]=fs->make<TH1F>("jet_discrCSV","CSV discriminant;CSV discriminant",100,0.,1.);
						
  // dijet mass (if available)
  histContainer_["di_jet_mass" ]=fs->make<TH1F>("di_jet_mass","M_{jj}; M_{jj} [GeV]",50,0.,500.);

  // muon variables (for all muons)
  histContainer_["muon_pt"]=fs->make<TH1F>("muon_pt","muon pt;muon p_{T} [GeV]",150,0.,150.);
  histContainer_["muon_eta"]=fs->make<TH1F>("muon_eta","muon eta; muon #eta",30,-3.,3.);
  histContainer_["muon_phi"]=fs->make<TH1F>("muon_phi","muon phi; muon #phi [rad];",60,-TMath::Pi(),TMath::Pi());
  histContainer_["muon_iso"]=fs->make<TH1F>("muon_iso","muon iso;(#sum_{trk}p_{T} + #sum_{cal}E_{T})/p_{T} (combRelIso);",60,0.,20.);
  histContainer_["muon_SIP"]=fs->make<TH1F>("muon_SIP","muon IP/#sigma_{IP};muon IP/#sigma_{IP}",30,0.,10.);
  histContainer_["isGlobalMuon"]=fs->make<TH1F>("isGlobalMuon","isGlobalMuon",2,-0.5,1.5); 
  histContainer_["isTrackerMuon"]=fs->make<TH1F>("isTrackerMuon","isTrackerMuon",2,-0.5,1.5);  
  histContainer_["muon_chi2"]= fs->make<TH1F>("muon_chi2","muon #chi^{2}/ndof; muon #chi^{2}/ndof", 100,0.,10.);
  histContainer_["muon_trackerhits"]=fs->make<TH1F>("muon_trackerhits","muon Trk hits;tracker hits",40,-0.5,39.5);
  histContainer_["muon_pixelhits"]=fs->make<TH1F>("muon_pixelhits","muon Pxl hits;pixel hits",10,-0.5,9.5);
  histContainer_["muon_muonhits"]=   fs->make<TH1F>("muon_muonhits","muon Muon hits;muon hits",100,-0.,99.5);
  histContainer_["muon_dB"]= fs->make<TH1F>("muon_dB","muon dB; dB [cm]",150,0.,10.);
  histContainer_["muon_numberOfMatches"]= fs->make<TH1F>("muon_numberOfMatches","muon number of matches;n.of matches",10,-0.5,9.5);
  histContainer_["AllDiMuMass"] = fs->make<TH1F>("AllDiMuMass","Invariant mass spectrum of all dimuons; M(#mu^{+}#mu^{-}) [GeV]",100,0,120);
  histContainer_["AllDiMuPt"] = fs->make<TH1F>("AllDiMuPt","p_{T} of all di-muons; p_{T}(#mu^{+}#mu^{-}) [GeV]",30,0,150);
  histContainer_["AllDiMuEta"] = fs->make<TH1F>("AllDiMuEta","#eta of all di-muons; #eta(#mu^{+}#mu^{-})",50,-5,5);
  histContainer_["AllDiMuPhi"] = fs->make<TH1F>("AllDiMuPhi","#phi of all di-muons; #phi(#mu^{+}#mu^{-}) [rad]",30,-TMath::Pi(),TMath::Pi());
  histContainer_["AllDiMuDeltaPhi"] = fs->make<TH1F>("AllDiMuDeltaPhi","#Delta #phi of all di-muons; #Delta#phi(#mu^{+}#mu^{-}) [rad]",60,0.,TMath::Pi());  
  histContainer_["AllDiMuDeltaR"] = fs->make<TH1F>("AllDiMuDeltaR","#Delta R of all di-muons; #DeltaR(#mu^{+}#mu^{-}) [rad]",100,0,10); 
  histContainer_["PreselMuMass"] = fs->make<TH1F>("PreselMuMass","Invariant mass spectrum of Presel dimuons; M(#mu^{+}#mu^{-}) [GeV]",100,0,120);
  histContainer_["PreselDiMuPt"] = fs->make<TH1F>("PreselDiMuPt","p_{T} of Presel di-muons; p_{T}(#mu^{+}#mu^{-}) [GeV]",30,0,150);
  histContainer_["PreselDiMuEta"] = fs->make<TH1F>("PreselDiMuEta","#eta of Presel di-muons; #eta(#mu^{+}#mu^{-})",50,-5,5);
  histContainer_["PreselDiMuPhi"] = fs->make<TH1F>("PreselDiMuPhi","#phi of Presel di-muons; #phi(#mu^{+}#mu^{-}) [rad]",30,-TMath::Pi(),TMath::Pi()); 
  histContainer_["PreselDiMuDeltaPhi"] = fs->make<TH1F>("PreselDiMuDeltaPhi","#Delta #phi of Presel di-muons; #Delta#phi(#mu^{+}#mu^{-}) [rad]",60,0.,TMath::Pi());  
  histContainer_["PreselDiMuDeltaR"] = fs->make<TH1F>("PreselDiMuDeltaR","#Delta R of Presel di-muons; #DeltaR(#mu^{+}#mu^{-}) [rad]",100,0,10);  

}

void 
PatBasicAnalyzer::endJob() 
{
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PatBasicAnalyzer);
