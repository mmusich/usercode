#include <iostream>
#include "FWCore/Framework/interface/ESHandle.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbUtils.h"
#include "ZbbAnalysis/AnalysisStep/interface/BTagWeight.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"
#include "RecoBTag/PerformanceDB/interface/BtagPerformance.h" 
#include "RecoBTag/Records/interface/BTagPerformanceRecord.h"
#include "CondFormats/PhysicsToolsObjects/interface/PerformancePayloadFromTable.h"
#include "CondFormats/PhysicsToolsObjects/interface/PerformanceWorkingPoint.h"
#include "ZbbAnalysis/AnalysisStep/interface/VtxAssociatorsUtils.h"
#include "ZbbAnalysis/AnalysisStep/interface/AcceptanceCuts.h"
#include <Math/VectorUtil.h>
#include "TLorentzVector.h"
#include "PhysicsTools/JetMCUtils/interface/CandMCTag.h"

using namespace ZbbUtils;

//______________________________________________________________________________
Double_t ZbbUtils::CostCS(reco::CompositeCandidate ZCand){
  // Calculation the Collins-Soper angle (adapted from code by R. Arnaldi)
  
  Double_t fMProton = 0.93827231;
  Double_t ebeam=3500.;  //temporary
  if(ebeam<=0){
    printf("Can not compute costCS with EBeam=%f\n",ebeam);
    return -999999999;
  }
  
  Double_t mp=fMProton;
  Double_t pbeam=TMath::Sqrt(ebeam*ebeam-mp*mp);
  
  const reco::Candidate* lepton1 = ZCand.daughter(0);
  const reco::Candidate* lepton2 = ZCand.daughter(1);
  
  Double_t pla10 = lepton1->px();
  Double_t pla11 = lepton1->py();
  Double_t pla12 = lepton1->pz();
  Double_t e1 = lepton1->energy();
  Int_t lep1Charge = lepton1->charge();
  
  Double_t pla20 = lepton2->px();
  Double_t pla21 = lepton2->py();
  Double_t pla22 = lepton2->pz();
  Double_t e2 = lepton2->energy();
  Int_t lep2Charge = lepton2->charge();

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

//______________________________________________________________________________
Double_t ZbbUtils::CosThetaStar(reco::CompositeCandidate ZCand){

  const reco::Candidate* lepton1 = ZCand.daughter(0);
  const reco::Candidate* lepton2 = ZCand.daughter(1);

  TLorentzVector p4_lp;
  TLorentzVector p4_lm;

  Double_t pla10 = lepton1->px();
  Double_t pla11 = lepton1->py();
  Double_t pla12 = lepton1->pz();
  Double_t e1 = lepton1->energy();
  Int_t lep1Charge = lepton1->charge();
  
  Double_t pla20 = lepton2->px();
  Double_t pla21 = lepton2->py();
  Double_t pla22 = lepton2->pz();
  Double_t e2 = lepton2->energy();

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

//______________________________________________________________________________
Bool_t ZbbUtils::isJetIdOk(const pat::Jet& jet,TString jetcategory){

  //if(debug_) std::cout<<"isJetIdOk"<<std::endl;

  Bool_t isgoodjetID = false;
  
  Double_t nhf = ( jet.neutralHadronEnergy() + jet.HFHadronEnergy() ) / jet.energy();
  Double_t nEF = jet.neutralEmEnergyFraction();
  Double_t nconstituents = jet.numberOfDaughters();
  Double_t chf = jet.chargedHadronEnergyFraction();
  Double_t nch = jet.chargedMultiplicity();
  Double_t cef = jet.chargedEmEnergyFraction();
  
  if(jetcategory=="tight"){
    if(nhf<0.90 && nEF<0.90 && nconstituents>1 && chf>0 && nch>0 && cef<0.99) isgoodjetID=true;
  } else if(jetcategory=="medium"){
    if(nhf<0.95 && nEF<0.95 && nconstituents>1 && chf>0 && nch>0 && cef<0.99) isgoodjetID=true;
  } else if(jetcategory=="loose"){
    if(nhf<0.99 && nEF<0.99 && nconstituents>1 && chf>0 && nch>0 && cef<0.99) isgoodjetID=true;
  } else {
    std::cout<<"isJetIdOk: Error unforeseen jet category. Use thight,medium or loose."<<std::endl;
  }
  return isgoodjetID; 
}

//______________________________________________________________________________
Bool_t ZbbUtils::isGoodJet(const pat::Jet& jet,const reco::CompositeCandidate& ZCand,const AcceptanceCuts& lCuts){
  
  //if(debug_) std::cout<<"isGoodJet"<<std::endl;

  Bool_t isgoodjet = false;

  if ( ZbbUtils::hasOverlap(jet.eta(),jet.phi(),ZCand.daughter(0)->eta(),ZCand.daughter(0)->phi(),0.5) ) return false;
  if ( ZbbUtils::hasOverlap(jet.eta(),jet.phi(),ZCand.daughter(1)->eta(),ZCand.daughter(1)->phi(),0.5) ) return false;
  isgoodjet = jet.pt()>lCuts.bjetPtMin_ && TMath::Abs(jet.eta())<lCuts.bjetEtaMax_ ;
  return isgoodjet;
}
                                                 
//______________________________________________________________________________
Bool_t ZbbUtils::isBJet(const pat::Jet& jet,TString theAlgoWP){
  
  //if(debug_) std::cout<<"isBJet"<<std::endl;

  Bool_t isbjet = false;

  if( theAlgoWP=="SSVHEM" ) {
    isbjet = jet.bDiscriminator("simpleSecondaryVertexHighEffBJetTags")>1.74;
  } else if ( theAlgoWP=="SSVHPT" ) {
      isbjet = jet.bDiscriminator("simpleSecondaryVertexHighPurBJetTags")>2.00;
  } else if ( theAlgoWP=="TCHEL" ) {
    isbjet =jet.bDiscriminator("trackCountingHighEffBJetTags")>1.70;
  } else if ( theAlgoWP=="TCHEM" ) {
    isbjet =jet.bDiscriminator("trackCountingHighEffBJetTags")>3.30;
  } else if ( theAlgoWP=="TCHPL" ) {
    isbjet = jet.bDiscriminator("trackCountingHighPurBJetTags")>1.19;
  } else if ( theAlgoWP=="TCHPT" ) {
    isbjet = jet.bDiscriminator("trackCountingHighPurBJetTags")>3.41;
  } else if ( theAlgoWP=="CSVM" ) {
    isbjet = jet.bDiscriminator("combinedSecondaryVertexBJetTags")>0.679;
  } else if ( theAlgoWP=="CSVT" ) {
    isbjet = jet.bDiscriminator("combinedSecondaryVertexBJetTags")>0.898;
  } else if  ( theAlgoWP=="JPT" ) {
    isbjet = jet.bDiscriminator("jetProbabilityBJetTags")>0.790;
  } 
  else {   
    std::cout<<"isBJet: Error unforeseen algo or working point for b-tagging. Use SSVHEM/SSVHPT/TCHEM/TCHPT/CSVM/CSVT/JPT "<<std::endl; 	   
  } 	    
  return isbjet;
}

//______________________________________________________________________________
Bool_t ZbbUtils::isBJetMCMatched(const pat::Jet& jet,edm::Handle<reco::GenParticleCollection> genParticlesCollection,const reco::GenJetCollection & genJets ){
  
  using namespace CandMCTagUtils;

  Bool_t isbjetMcMatched_ = false;
  
  // std::cout<<"inside method"<<std::endl;

  for( reco::GenParticleCollection::const_iterator genp = genParticlesCollection->begin();genp != genParticlesCollection->end(); ++ genp ) {  // loop over GEN particles

    // std::cout<<"inside genParticles loop"<<std::endl;

    //do a copy of the genParticle, not to skip anythin in the loop
    const reco::Candidate* p = &(*genp);
    
    //if b-hadrons
    if( hasBottom(*p) ){

      //std::cout<<"has bottom"<<std::endl;

      // looks if all daughter are bottomless
      Bool_t hasBottomedDaughter = false;
      for(UInt_t i=0; i< p->numberOfDaughters(); i++){
	if(hasBottom(*(p->daughter(i)))){	
	  hasBottomedDaughter=true;
	}
      }
      
      // if b-hadron is sensible
      if(!hasBottomedDaughter && TMath::Abs(p->eta())<3.5 ){
	
	//std::cout<<"is good bottom jet"<<std::endl;

	math::XYZTLorentzVectorD p4HFGEN(p->px(),p->py(),p->pz(),p->energy());
	std::pair<int,double> GenGenAssociation = std::make_pair(-1,9999.);
	double minDeltaRGenGen(9999.);
	int i(0);
	for(std::vector<reco::GenJet>::const_iterator genjet_it=genJets.begin(); genjet_it!=genJets.end(); ++genjet_it){ 
	  if(ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4HFGEN.Vect())<0.5) { 
	    if(ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4HFGEN.Vect())< minDeltaRGenGen){
	      minDeltaRGenGen = ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4HFGEN.Vect());
	      GenGenAssociation.first  = i;
	      GenGenAssociation.second = minDeltaRGenGen;
	    } // end if minimum distance gen jet-parton 
	  } // ends if parton matched to gen jet
	  //std::cout<<"i: "<<i<<" minDeltaR: "<<minDeltaRGenGen<<std::endl;
	  i++;
	}// ends loop on gen jets

 	if(GenGenAssociation.first == -1 ) continue;
	
	if(ROOT::Math::VectorUtil::DeltaR(genJets.at(GenGenAssociation.first).momentum(),jet.momentum())<0.5) {
	  isbjetMcMatched_=true;
	}// end if distance gen jet - reco jet < 0.5
      }// ends if b-hadron is sensible
    }// ends if hasbottom
  }// ends loop on genParticles 

  return isbjetMcMatched_;
}

//______________________________________________________________________________
Bool_t ZbbUtils::isZLLCandidateMCMatched(const reco::CompositeCandidate& ZCand,edm::Handle<reco::GenParticleCollection> genParticlesCollection){
  
  Bool_t isZLLCandidateMatched_= false;

  const reco::Candidate* lepton1 = ZCand.daughter(0);
  const reco::Candidate* lepton2 = ZCand.daughter(1);

  std::vector<const reco::Candidate*> leptons;
  leptons.push_back(lepton1);
  leptons.push_back(lepton2);
  
  Double_t n_matched_leptons_(0.);

  Int_t nLeptonsInLoop(0);
  for(std::vector<const reco::Candidate*>::const_iterator theZlepton = leptons.begin(); theZlepton != leptons.end(); ++theZlepton) {
    nLeptonsInLoop++;
    for(reco::GenParticleCollection::const_iterator genp = genParticlesCollection->begin();genp != genParticlesCollection->end(); ++ genp ) { 
      // lepton selection of selected flavour, status 3 and daughter of a Z
      if( (fabs(genp->pdgId()) == 11 || fabs(genp->pdgId()) == 13) ){ //(ZbbUtils::getParentCode(&(*genp))==23) ){
	// std::cout<<"was a Z daughter "<<std::endl;
	// dump the 4-momentum into a LorentzVector
	math::XYZTLorentzVectorD p4LeptonGEN(genp->px(),genp->py(),genp->pz(),genp->energy());
	// gen lepton to reco matching
	if((*theZlepton)->isMuon()){
	  //std::cout<<"was a reco muon "<<std::endl;
	  if( fabs(genp->pdgId())==13 && (genp->status()==1)){
	    //std::cout<<"was a gen muon"<<std::endl;
	    if(ROOT::Math::VectorUtil::DeltaR((*theZlepton)->momentum(),p4LeptonGEN.Vect())<0.1){
	      n_matched_leptons_++;
	      break;
	    }
	  } 
	} else if((*theZlepton)->isElectron()){
	  if( fabs(genp->pdgId())==11 && (genp->status()==3)) {
	    if(ROOT::Math::VectorUtil::DeltaR((*theZlepton)->momentum(),p4LeptonGEN.Vect())<0.1){
	      n_matched_leptons_++;
	      break;
	    } // loop on reco electrons
	  } // flavour switch
	} // if is reco::Electron
      } // if gen lepton is ok
    } // ends loop on genParticles
    //   std::cout<<"nLeptonsInLoop: "<<nLeptonsInLoop<<std::endl;
  } // ends loop on reco::Leptons

  if (n_matched_leptons_>=2) isZLLCandidateMatched_ = true;
 
  if(isZLLCandidateMatched_){
    //std::cout<<"n_matched_leptons_: "<<n_matched_leptons_<<std::endl;
  }
  
  return isZLLCandidateMatched_;
  
}

//______________________________________________________________________________
Bool_t ZbbUtils::isGoodMet(const pat::MET& met, Double_t theMetCut){
  
  //if(debug_) std::cout<<"isGoodMet"<<std::endl;

  Bool_t isgoodmet = false;

  if(met.pt()<theMetCut){
    isgoodmet = true;
  }
  return isgoodmet;
}

//______________________________________________________________________________
Bool_t ZbbUtils::hasOverlap(Double_t etaJ, Double_t phiJ, Double_t etaLept,Double_t phiLept, Double_t theDeltaR){
  return  deltaR(etaJ,phiJ,etaLept,phiLept) < theDeltaR;
}

//______________________________________________________________________________
std::vector<reco::CompositeCandidate> ZbbUtils::sortCandidatesByDifference(std::vector<reco::CompositeCandidate> unsortedCands){

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
std::vector<pat::Jet> ZbbUtils::sortJetsBypT(std::vector<pat::Jet> unsortedJets){

  std::vector<pat::Jet> sortedJets = unsortedJets;
  //Double_t pT_max = 0;
  std::vector<Float_t> pTJets;
  UInt_t JetsSize = unsortedJets.size();

  for (pat::JetCollection::const_iterator jet = unsortedJets.begin(); jet!= unsortedJets.end(); ++jet){
    pTJets.push_back(jet->pt());
  }

  for (UInt_t i = 0; i < JetsSize; ++i){
    for (UInt_t j = i+1; j < JetsSize; ++j){
      if (pTJets[i] < pTJets[j]){
	pat::Jet auxJet = sortedJets[i];
	sortedJets[i] = sortedJets[j];
	sortedJets[j] = auxJet;
      }// if
    }// j loop
  }// i loop

  return sortedJets;

}


//______________________________________________________________________________
template<class T>
std::vector<T> ZbbUtils::sortObjectsBypT(std::vector<T> const& unsortedObj){

  std::vector<T> sortedObj= unsortedObj;
  std::vector<Float_t> pTObj;
  UInt_t ObjSize=unsortedObj.size();
  for (UInt_t n=0;n <ObjSize; ++n){
    pTObj.push_back(unsortedObj[n].pt());
  }

  for (UInt_t i = 0; i < ObjSize; ++i){
    for (UInt_t j = i+1; j < ObjSize; ++j){
      if (pTObj[i] < pTObj[j]){
	T auxObj = sortedObj[i];
	sortedObj[i] = sortedObj[j];
	sortedObj[j] = auxObj;
      }// if
    }// j loop
  }// i loop

  return sortedObj;

}

//______________________________________________________________________________
Bool_t ZbbUtils::isTightZCandidate(reco::CompositeCandidate ZCand, const reco::BeamSpot& beamSpot, Bool_t isMuChannel, const AcceptanceCuts& lCuts){

  Bool_t istightZcandidate = false;
  const reco::Candidate* lep0 = ZCand.daughter(0);
  const reco::Candidate* lep1 = ZCand.daughter(1);

  if(isMuChannel){
    const pat::Muon* muon0 = dynamic_cast<const pat::Muon*>(&(*lep0->masterClone()));
    const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(&(*lep1->masterClone()));    

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
    const pat::Electron* ele0 = dynamic_cast<const pat::Electron*>(&(*lep0->masterClone()));
    const pat::Electron* ele1 = dynamic_cast<const pat::Electron*>(&(*lep1->masterClone()));
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
Bool_t ZbbUtils::isPreselMu(const pat::Muon& muon, const reco::Vertex& vertex,const AcceptanceCuts& lCuts){

  Bool_t isPreselMu_ = false;
  
  if ( muon.pt()> 5. &&
       muon.isGlobalMuon()==true && 
       muon.isTrackerMuon()==true && 
       muon.normChi2() < 15 && 
       muon.innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 && 
       muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0  &&
       muon.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && 
       VtxAssociatorsUtils::tightzVertexMu(muon,0.1,vertex) &&
       //muon.dB() < 0.2 && 
       //(muon.trackIso()+muon.caloIso()) < 0.15*muon.pt() && 
       //muon.trackIso() < 0.7*(muon.pt()) &&
       muon.numberOfMatches() > 1 && 
       abs(muon.eta()) < lCuts.muEtaMax_ ) isPreselMu_=true;

  return isPreselMu_;
}

//______________________________________________________________________________
Bool_t ZbbUtils::isPreselEle(const pat::Electron& electron, const reco::Vertex& vertex,const AcceptanceCuts& lCuts){

  Bool_t isPreselEle_=false;

  if (electron.pt() > 5. &&
      (electron.isEE() || electron.isEB()) && 
      !electron.isEBEEGap() && 
      electron.electronID("eidVBTFRel95") == 5 &&
      VtxAssociatorsUtils::tightzVertexEle(electron,0.1,vertex) &&
      //electron.trackIso()<0.7*(electron.pt()) &&
      abs(electron.eta())< lCuts.eleEtaMax_ ) isPreselEle_=true;

  return isPreselEle_;

}

//______________________________________________________________________________
Bool_t ZbbUtils::ProduceTheCut(const EventCategory& ec,TString categoryname, TString extractThisCut){
  
  std::map<TString,Bool_t> theEvtCategoryCut;
  
  //define event categories and cuts:
  TString eventcategories[12] = {"AllEvts",
				 "TrigOk",
				 "LL",
				 "LLTightAndTrigMatch",
				 "ZLL",
				 "GoodJet",
				 "JetBTag",
				 "JetBTagExcl",   
				 "JetBTagPlusl",
				 "JetBTagPlus2l",
				 "JetBTagMCMatched",
				 "GoodJetAndVertexAssoc"};
  
  Bool_t  cuts[12]            = {true,
				 ec.isTrigger(),
				 ec.isLL(),
				 ec.isLLTightAndTriggerMatched(),
				 ec.isZLL(),
				 ec.isGoodJet(),
				 ec.isJetBTagged(),
				 ec.isExclusive(),
				 ec.isThreeLeptons(),
				 ec.isFourLeptons(),
				 ec.isJetBTaggedMCMatched(),
				 ec.isJetVertexAssoc()};
  
  Bool_t theCut=true;

  for(UInt_t cat=0; cat<10; cat++){
    
    if(eventcategories[cat]==extractThisCut) continue; 
    
    if(cat>=8){  // exclusive f.s. not in the sequence if additional leptons
      if(eventcategories[cat]=="JetBTagExcl") continue; 
    }

    theCut =  theCut&&cuts[cat];
    theEvtCategoryCut[eventcategories[cat]] = theCut;
  }
  
  //separate case for mc truth matching
  theEvtCategoryCut["JetBTagMCMatched"]      = (theEvtCategoryCut.find("JetBTag")->second && ec.isJetBTaggedMCMatched());
  theEvtCategoryCut["GoodJetAndVertexAssoc"] = (theEvtCategoryCut.find("GoodJet")->second && ec.isJetVertexAssoc());

  // returns the cut
  return theEvtCategoryCut.find(categoryname)->second;
  
}

//______________________________________________________________________________
Double_t ZbbUtils::getTheWeight(const EventCategory& ec,TString categoryname){

  Double_t w(1);
  
  // corrects the weight for exclusive final state  
  if(categoryname=="JetBTagExcl"){
    w =  ec.getWeight()*ec.getBTagEffWeightCF();
  } else {
    w = ec.getWeight();
  }
  
  return w;
  
}

//______________________________________________________________________________
Double_t ZbbUtils::myDeltaPhi(Double_t phi1, Double_t phi2){
  // in CMSSW phi = [0,2pi], in TLorentzVector phi = [-pi,pi].
  // With the conversion below deltaPhi works ok despite the
  // 2*pi difference in phi definitions.
  Double_t PI = 3.1415;
  if(phi1 < 0) phi1 += 2*PI;
  if(phi2 < 0) phi2 += 2*PI;
  
  Double_t dphi = fabs(phi1-phi2);
  
  if(dphi > PI) dphi = 2*PI - dphi;
  return dphi;
}

//______________________________________________________________________________
FLAVOUR ZbbUtils::getPartonFlavour(const pat::Jet& jet){

  FLAVOUR flavour;
  //flavour switch
  switch(std::abs(jet.partonFlavour())) {
  case 1:
  case 2:
  case 3:
  case 21:
    flavour = ZbbUtils::UDSG_JETS;
    break;
  case 4:
    flavour = ZbbUtils::C_JETS;
    break;
  case 5:
    flavour = ZbbUtils::B_JETS;
    break;
  default:
    flavour = ZbbUtils::NONID_JETS;
  }
  return flavour;
}

//______________________________________________________________________________
double ZbbUtils::getbEffScaleFactor(FLAVOUR jetFlavour, const pat::Jet& jet, const edm::EventSetup& iSetup, std::string beffcalibmethod, std::string bmistagcalibmethod){

  double w(1);
  
  // it would be better determine the jet flavour inside this method....
  if ( jetFlavour != B_JETS ) return w; 
  double ETjet=jet.pt();
  double AbsEtajet=fabs(jet.eta());

  edm::ESHandle<BtagPerformance> perfEFFH;
  iSetup.get<BTagPerformanceRecord>().get(beffcalibmethod,perfEFFH);
  const BtagPerformance & pbeff = *(perfEFFH.product());

  edm::ESHandle<BtagPerformance> perfMISH;
  iSetup.get<BTagPerformanceRecord>().get(bmistagcalibmethod,perfMISH);
  const BtagPerformance & pleff = *(perfMISH.product());

  //  std::cout << "ET/eta " << ETjet << " / " << AbsEtajet << std::endl;
  //  std::cout <<" Discriminant is "<<pbeff.workingPoint().discriminantName()<<std::endl;
  //  std::cout << "Working point: " << pbeff.workingPoint().cut() << std::endl;

  BinningPointByMap p;
  p.reset();
  p.insert(BinningVariables::JetAbsEta,AbsEtajet);
  p.insert(BinningVariables::JetEt,ETjet);

  // it would be better determine the jet flavour inside this method....
  if ( jetFlavour == B_JETS ) {
    if ( pbeff.isResultOk(PerformanceResult::BTAGBEFFCORR,p) )
      w=pbeff.getResult(PerformanceResult::BTAGBEFFCORR,p);
  } else {
    if ( pleff.isResultOk(PerformanceResult::BTAGLEFFCORR,p) )
      w=pleff.getResult(PerformanceResult::BTAGLEFFCORR,p);
  }

  // sanity check
  if ( w<0.25 || w>4 ) w=1;

  return w;
}

//______________________________________________________________________________
Double_t ZbbUtils::getMuonOfflineScaleFactor(double pt, double eta, bool uncert){

  // y for pt, x for eta

  double sf_muons_ID[2][2]={{0.995,0.994},{0.998,0.991}};
  //double sf_muons_ID_err[2][2]={{0.004, 0.004}, {0.004, 0.005}};

  double sf_muons_ISO[2][2]={{1.015, 0.999},{1.013, 1.002}};
  //double sf_muons_ISO_err[2][2]={{0.004, 0.004}, {0.004, 0.005}};
  
  //default
  double w(1);
  //double w_err(1); 
  
  int idx=-1;
  int idy=-1;

  int isox=-1;
  int isoy=-1;
  
  // -------- pt ----------------------------------
  if (pt<=50) {
    idx  =0;
    isox =0;
   } else if (pt>50) {
    idx  =1;
    isox =1;
  }
  
  // -------- eta ---------------------------------
  if (TMath::Abs(eta)<=0.9) {
    idy  =0;
    isoy =0;
  } else if (TMath::Abs(eta)>0.9 && TMath::Abs(eta)<=1.2) {   
    idy  =0;
    isoy =0;
  } else if (TMath::Abs(eta)>1.2 && TMath::Abs(eta)<=2.1) {
    idy  =1;
    isoy =1;
  } else if (TMath::Abs(eta)>2.1 && TMath::Abs(eta)<=2.4) {
    idy  =1;
    isoy =1;
  }
  
  // --------------- final SF computation (and error) ------------- 
  if (idx>=0 && idy>=0 && isox>=0 && isoy>=0) {
    w = sf_muons_ID[idx][idy]*sf_muons_ISO[isox][isoy];
    //w_err = TMath::Sqrt( TMath::Power(sf_muons_ID_err[idx][idy]*sf_muons_ISO[isox][isoy]*sf_muons_TRG[0][trgy],2)+ 
    //			 TMath::Power(sf_muons_ISO_err[isox][isoy]*sf_muons_ID[idx][idy]*sf_muons_TRG[0][trgy],2)+
    //			 TMath::Power(sf_muons_TRG_err[0][trgy]*sf_muons_ID[idx][idy]*sf_muons_ISO[isox][isoy],2));// errors to be propagated
  }
  return w;
}

//______________________________________________________________________________
Double_t ZbbUtils::getMuonTrgScaleFactor_H( double pt, double eta,const std::string& runPeriod, bool uncert){

  double sf_muons_TRG[1][3];  
  // double sf_muons_TRG_err[1][3];

  if(runPeriod=="2011A1"){
    sf_muons_TRG[0][0]=0.971;
    sf_muons_TRG[0][1]=0.957;
    sf_muons_TRG[0][2]=0.954;
   //double  sf_muons_TRG_err[1][3]={{0.001, 0.001,0.005}};
  }
  else if(runPeriod=="2011A2"){
    sf_muons_TRG[0][0]=0.973;
    sf_muons_TRG[0][1]=0.962;
    sf_muons_TRG[0][2]=0.946;
    //double  sf_muons_TRG_err[1][3]={{0.001, 0.001,0.005}};   
  }
  else{
    std::cout << "WRONG PERIOD!! " << std::endl;  
    sf_muons_TRG[0][0]=0.;
    sf_muons_TRG[0][1]=0.;
    sf_muons_TRG[0][2]=0.;     
     //double  sf_muons_TRG_err[1][3]={{0., 0.,0.}};
  }
  
  double w(1);
  //double w_err(1); 
  
  int trgy=-1;
  
  // --------eta ----------------------------------
  if (TMath::Abs(eta)<=0.9) {
    trgy =0;
  } else if (TMath::Abs(eta)>0.9 && TMath::Abs(eta)<=2.1) {   
    trgy =1;
  } else if (TMath::Abs(eta)>2.1 && TMath::Abs(eta)<=2.4) {
    trgy =2;
  }
  
  // --------------- final SF computation (and error) ------------- 
  if (trgy>=0) {
    w = sf_muons_TRG[0][trgy];
    //    w_err = TMath::Sqrt( TMath::Power(sf_muons_ID_err[idx][idy]*sf_muons_ISO[isox][isoy]*sf_muons_TRG[0][trgy],2)+ 
    // 			 TMath::Power(sf_muons_ISO_err[isox][isoy]*sf_muons_ID[idx][idy]*sf_muons_TRG[0][trgy],2)+
    // 			 TMath::Power(sf_muons_TRG_err[0][trgy]*sf_muons_ID[idx][idy]*sf_muons_ISO[isox][isoy],2));// errors to be propagated
  }
  return w;
}

//______________________________________________________________________________
Double_t ZbbUtils::getMuonTrgScaleFactor_L( double pt, double eta, const std::string& runPeriod, bool uncert){
  
  double sf_muons_TRG[1][3];
  //  double sf_muons_TRG_err[1][3];
  if(runPeriod=="2011A1"){
    sf_muons_TRG[0][0]=0.971;
    sf_muons_TRG[0][1]=0.957;
    sf_muons_TRG[0][2]=0.954;
    // double  sf_muons_TRG_err[1][3]={{0.001, 0.001,0.005}};
  }
  else if(runPeriod=="2011A2"){
    sf_muons_TRG[0][0]=0.973;
    sf_muons_TRG[0][1]=0.964;
    sf_muons_TRG[0][2]=0.952;
    //   double  sf_muons_TRG_err[1][3]={{0.001, 0.001,0.005}};
  }
  else{
    std::cout << "WRONG PERIOD!! " << std::endl;  
    sf_muons_TRG[0][0]=0.;
    sf_muons_TRG[0][1]=0.;
    sf_muons_TRG[0][2]=0.;
    // double  sf_muons_TRG_err[1][3]={{0., 0.,0.}};
  }
  
  double w(1);
  //double w_err(1); 
  
  int trgy=-1;
  
  // --------eta ----------------------------------
  if (TMath::Abs(eta)<=0.9) {
    trgy =0;
  } else if (TMath::Abs(eta)>0.9 && TMath::Abs(eta)<=2.1) {   
    trgy =1;
  }
  else if (TMath::Abs(eta)>2.1 && TMath::Abs(eta)<=2.4) {
    trgy =2;
  }
  // --------------- final SF computation (and error) ------------- 
  if (trgy>=0) {
    w = sf_muons_TRG[0][trgy];
    //    w_err = TMath::Sqrt( TMath::Power(sf_muons_ID_err[idx][idy]*sf_muons_ISO[isox][isoy]*sf_muons_TRG[0][trgy],2)+ 
    // 			 TMath::Power(sf_muons_ISO_err[isox][isoy]*sf_muons_ID[idx][idy]*sf_muons_TRG[0][trgy],2)+
    // 			 TMath::Power(sf_muons_TRG_err[0][trgy]*sf_muons_ID[idx][idy]*sf_muons_ISO[isox][isoy],2));// errors to be propagated
  }
  return w;
}


//______________________________________________________________________________
Double_t ZbbUtils::getElectronOfflineScaleFactor(double pt, double eta, bool uncert){
  double sf_electrons_RECO[2][4]={{1.002, 1.001, 1.009, 1.008}, {1.004, 0.984, 0.993, 0.944}};
  //double sf_electrons_ID_err[2][4]={{0.005, 0.005, 0.005, 0.005}, {0.005, 0.005, 0.007, 0.010}};
  
  double sf_electrons_IDISO[2][2]={{1.004, 1.043}, {0.988, 1.013}};
  //double sf_electrons_ISO_err[2][2]={{0.005, 0.005}, {0.005, 0.006}};

  //double sf_electrons_TRG[1][2]={{0.999, 0.999}};
  //double sf_electrons_TRG_err[1][2]={{0.01, 0.01}};

  // defaults  
  double w(1);
  //double w_err(1); 

  int idx=-1;
  int idy=-1;

  int isox=-1;
  int isoy=-1;
    
  // -------- pt ------------------------------
  if (pt<=50) {
    idx  =0;
    isox =0;
   } else if (pt>50) {
    idx  =1;
    isox =1;
  }
  
  // -------- eta -----------------------------
  if (eta<=0.8) {
    idy  =0;
    isoy =0;
  } else if (TMath::Abs(eta)>0.8 && TMath::Abs(eta)<=1.44) {   
    idy  =1;
    isoy =0;
  } else if (TMath::Abs(eta)>1.57 && TMath::Abs(eta)<=1.6) {
    idy  =2;
    isoy =0;
  } else if (TMath::Abs(eta)>1.6 && TMath::Abs(eta)<=2.0) {
    idy  =2;
    isoy =1;
  } else if (TMath::Abs(eta)>2.0 && TMath::Abs(eta)<=2.5) {
    idy =3;
    isoy =1;
  }

  // --------------- final SF computation (and error) ------------- 
  if (idx>=0 && idy>=0 && isox>=0 && isoy>=0){
    w = sf_electrons_RECO[idx][idy]*sf_electrons_IDISO[isox][isoy];
    //w_err = TMath::Sqrt( TMath::Power(sf_electrons_ID_err[idx][idy]*sf_electrons_ISO[isox][isoy]*sf_electrons_TRG[0][trgy],2)+ 
    //			 TMath::Power(sf_electrons_ISO_err[isox][isoy]*sf_electrons_ID[idx][idy]*sf_electrons_TRG[0][trgy],2)+
    //			 TMath::Power(sf_electrons_TRG_err[0][trgy]*sf_electrons_ID[idx][idy]*sf_electrons_ISO[isox][isoy],2));// errors to be propagated
  }
  return w;
}

//______________________________________________________________________________
Double_t ZbbUtils::getElectronTrgScaleFactor_H( double pt, double eta, const std::string& runPeriod, bool uncert){

  double sf_electrons_TRG[1][3];

  if(runPeriod=="2011A1"){
    sf_electrons_TRG[0][0]=0.990;
    sf_electrons_TRG[0][1]=0.994;
    sf_electrons_TRG[0][2]=0.992;
  } else if(runPeriod=="2011A2"){
    sf_electrons_TRG[0][0]=0.990;
    sf_electrons_TRG[0][1]=0.989;
    sf_electrons_TRG[0][2]=0.989;
  } 
  else{
    std::cout << "WRONG PERIOD!! " << std::endl;  
    sf_electrons_TRG[0][0]=0.;
    sf_electrons_TRG[0][1]=0.;
    sf_electrons_TRG[0][2]=0.;
  }
  
  double w(1);
  //double w_err(1); 
  
  int trgy=-1;
  
  // --------eta ----------------------------------
  if (TMath::Abs(eta)<=0.9) {
    trgy =0;
  } else if (TMath::Abs(eta)>0.9 && TMath::Abs(eta)<=2.1) {   
    trgy =1;
  } else if (TMath::Abs(eta)>2.1 && TMath::Abs(eta)<=2.4) {
    trgy =2;
  }
  // --------------- final SF computation (and error) ------------- 
  if (trgy>=0) {
    w = sf_electrons_TRG[0][trgy];
    //w_err = TMath::Sqrt( TMath::Power(sf_electrons_ID_err[idx][idy]*sf_electrons_ISO[isox][isoy]*sf_electrons_TRG[0][trgy],2)+ 
    // 			 TMath::Power(sf_electrons_ISO_err[isox][isoy]*sf_electrons_ID[idx][idy]*sf_electrons_TRG[0][trgy],2)+
    // 			 TMath::Power(sf_electrons_TRG_err[0][trgy]*sf_electrons_ID[idx][idy]*sf_electrons_ISO[isox][isoy],2));// errors to be propagated
  }
  return w;
}

//______________________________________________________________________________
Double_t ZbbUtils::getElectronTrgScaleFactor_L( double pt, double eta, const std::string& runPeriod, bool uncert){

  double sf_electrons_TRG[1][3];

  if(runPeriod=="2011A1"){
    sf_electrons_TRG[0][0]=0.987;
    sf_electrons_TRG[0][1]=0.989;
    sf_electrons_TRG[0][2]=0.989;
  } else if(runPeriod=="2011A2"){
    sf_electrons_TRG[0][0]=0.991;
    sf_electrons_TRG[0][1]=0.989;
    sf_electrons_TRG[0][2]=0.989;
  } else{
    std::cout << "WRONG PERIOD!! " << std::endl;  
    sf_electrons_TRG[0][0]=0.;
    sf_electrons_TRG[0][1]=0.;
    sf_electrons_TRG[0][2]=0.;
  }
 
  double w(1);
  //double w_err(1); 
  
  int trgy=-1;
  
  // --------eta ----------------------------------
  if (TMath::Abs(eta)<=0.9) {
    trgy =0;
  } else if (TMath::Abs(eta)>0.9 && TMath::Abs(eta)<=2.1) {   
    trgy =1;
  }
  else if (TMath::Abs(eta)>2.1 && TMath::Abs(eta)<=2.4) {
    trgy =2;
  }
  // --------------- final SF computation (and error) ------------- 
  if (trgy>=0) {
    w = sf_electrons_TRG[0][trgy];
    //w_err = TMath::Sqrt( TMath::Power(sf_electrons_ID_err[idx][idy]*sf_electrons_ISO[isox][isoy]*sf_electrons_TRG[0][trgy],2)+ 
    // 			 TMath::Power(sf_electrons_ISO_err[isox][isoy]*sf_electrons_ID[idx][idy]*sf_electrons_TRG[0][trgy],2)+
    //			 TMath::Power(sf_electrons_TRG_err[0][trgy]*sf_electrons_ID[idx][idy]*sf_electrons_ISO[isox][isoy],2));// errors to be propagated
  }
  return w;
}

//______________________________________________________________________________
double ZbbUtils::getbEffScaleFactorAR(std::map<EffMCpair,std::vector<double> > theAssociativeMapEffbMC_,
				      std::map<EffMCpair,std::vector<double> > theAssociativeMapEffcMC_,
				      std::vector<pat::Jet> btagjets, const edm::EventSetup& iSetup, std::string beffcalibmethod, std::string bmistagcalibmethod,bool isExclusive,int minNBtag,int theSFkFactor_,int theLEffkFactor_){
  
  double w(1);

  // if no btag jets just return 1
  if ( btagjets.size()<1 ) return w; 

  edm::ESHandle<BtagPerformance> perfEFFH;
  iSetup.get<BTagPerformanceRecord>().get(beffcalibmethod,perfEFFH);
  const BtagPerformance & pbeff = *(perfEFFH.product());

  edm::ESHandle<BtagPerformance> perfMISH;
  iSetup.get<BTagPerformanceRecord>().get(bmistagcalibmethod,perfMISH);
  const BtagPerformance & pleff = *(perfMISH.product());

  std::vector<BTagWeight::JetInfo> vJetInfo;

  // first iterate over bTag jets 
  // for each determine flavour/SF/effMC
  double ETjet, AbsEtajet;
  BinningPointByMap p;

  std::vector<pat::Jet>::const_iterator btagjet_it;
  for (btagjet_it=btagjets.begin(); btagjet_it!=btagjets.end(); ++btagjet_it) {
    ETjet=btagjet_it->pt();
    AbsEtajet=fabs(btagjet_it->eta());

    p.reset();
    p.insert(BinningVariables::JetAbsEta,AbsEtajet);
    p.insert(BinningVariables::JetEt,ETjet);

    FLAVOUR jetFlavour = getPartonFlavour(*(btagjet_it));    

    double mcEff(1), mcEfferror(0), correctedMcEff(1), dataSF(1), dataSFerror(0), correctedDataSF(1);
  
    if ( jetFlavour == B_JETS ) {    // ----------- b-jets
      if ( pbeff.isResultOk(PerformanceResult::BTAGBEFFCORR,p) ){
	dataSF=pbeff.getResult(PerformanceResult::BTAGBEFFCORR,p);
	dataSFerror=pbeff.getResult(PerformanceResult::BTAGBERRCORR,p);
      }

      // MC b eff from user table -> loop over the map
      EffMCpair mcEffPair_;
      std::map<EffMCpair,std::vector<Double_t> >::iterator it;  
      for ( it=theAssociativeMapEffbMC_.begin(); it !=theAssociativeMapEffbMC_.end(); it++ ){
	// if pT in the range
	if( ETjet >=(*it).second[0] && ETjet <=(*it).second[1]){
	  mcEffPair_ = (*it).first;
	  if ( AbsEtajet < 1.2 ) {
	    mcEff = mcEffPair_.first;
	  } else if ( AbsEtajet>1.2 && AbsEtajet<2.4 ) {
	    mcEff = mcEffPair_.second;
	  } else {
	    std::cout << "THIS JET SHOULD NOT BE THERE! " << std::endl;
	  }	    
	  break;
	}// ends if on pT range
      }// ends loop over b map
      // std::cout << "### B JET ###  ET "<< ETjet << " eta " << AbsEtajet << " effb MC " << mcEff<< std::endl;

    } else if ( jetFlavour == C_JETS ) { // ----------- c-jets  (same SF as b-jets)
      if ( pbeff.isResultOk(PerformanceResult::BTAGBEFFCORR,p) ){
	dataSF=pbeff.getResult(PerformanceResult::BTAGBEFFCORR,p);
	dataSFerror=pbeff.getResult(PerformanceResult::BTAGBERRCORR,p);
      }
      
      // MC c eff from user table -> loop over the map
      EffMCpair mcEffPair_;
      std::map<EffMCpair,std::vector<Double_t> >::iterator it;  
      for ( it=theAssociativeMapEffcMC_.begin() ; it !=theAssociativeMapEffcMC_.end(); it++ ){
	// if pT in the range
	if( ETjet >=(*it).second[0] && ETjet <=(*it).second[1]){
	  mcEffPair_ = (*it).first;
	  if ( AbsEtajet < 1.2 ) {
	    mcEff = mcEffPair_.first;
	  } else if ( AbsEtajet>1.2 && AbsEtajet<2.4 ) {
	    mcEff = mcEffPair_.second;
	  } else {
	    std::cout << "THIS JET SHOULD NOT BE THERE! " << std::endl;
	  }	    
	  break;
	}// ends if on pT range
      }// ends loop over map
      // std::cout << "### C JET ###  ET "<< ETjet << " eta " << AbsEtajet << " effb MC " << mcEff<< std::endl;

    } else {      

      // light-jets FIXME: treating not identified jets as light jets
      if ( pleff.isResultOk(PerformanceResult::BTAGLEFFCORR,p) ){
	dataSF=pleff.getResult(PerformanceResult::BTAGLEFFCORR,p);
	dataSFerror=pleff.getResult(PerformanceResult::BTAGLERRCORR,p);
      }

      // EM check if it exists...
      if ( pleff.isResultOk(PerformanceResult::BTAGLEFF,p) ){
	mcEff=pleff.getResult(PerformanceResult::BTAGLEFF,p);
	mcEfferror=pleff.getResult(PerformanceResult::BTAGLERR,p); 
      }
      // std::cout << "### UDSG JET ###  ET "<< ETjet << " eta " << AbsEtajet << " effb MC " << mcEff<< std::endl;
    }

    // std::cout << "### JET ### " << jetFlavour << " ET "<< ETjet << " eta " << AbsEtajet << " eff/SF MC " << mcEff<<" / " <<dataSF << std::endl;

    correctedDataSF=dataSF+(theSFkFactor_*dataSFerror);
    correctedMcEff=mcEff+(theLEffkFactor_*mcEfferror);

    // protect against crazy values
    if ( dataSF>0.25 && dataSF<4 ) {
      BTagWeight::JetInfo jetInfo(correctedMcEff,correctedDataSF);
      vJetInfo.push_back(jetInfo);
    }    
  } // end of loop on b-tag jets

  //finally get the weight
  if(isExclusive){
    BTagWeight btw1(minNBtag,minNBtag);
    w = btw1.weight(vJetInfo,minNBtag);
  } else {
    BTagWeight btw(minNBtag,vJetInfo.size());
    w = btw.weight(vJetInfo,vJetInfo.size()); 
    //  std::cout << "Number of b-tagged jets: " << vJetInfo.size() << " weight: " << w << std::endl;
  }
  return w;
}

//______________________________________________________________________________
std::vector<Double_t> ZbbUtils::generate_pu_weights(TH1D* data_npu_estimated){

  //********************* FOR 2011 ****************************
  // Flat in PU until npu=10 (for Spring 11 production) 
  // see SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py; copy and paste from there:
  const Double_t npu_probs[25] = {0.0698146584, 0.0698146584, 0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584 /* <-- 10*/,
				  0.0630151648,0.0526654164,0.0402754482,0.0292988928,0.0194384503,0.0122016783,0.007207042,0.004003637,0.0020278322,
				  0.0010739954,0.0004595759,0.0002229748,0.0001028162,4.58337152809607E-05 /* <-- 24 */};

  //********************* FOR 2010 ****************************
  // pseudo-poissonian PU (for Winter 10 production) 
  // see https://twiki.cern.ch/twiki/bin/view/CMS/PileupInformation#2010_Pileup_Scenarios
  // const Double_t npu_probs[10] = {0.145168,0.251419,0.251596,0.17943,0.10,0.05,0.02,0.01,0.005,0.002,0.001}

  Int_t maxNPU =sizeof(npu_probs);

  std::vector<Double_t> result(maxNPU);

  Double_t s = 0.0;
  for(int npu=0; npu<maxNPU; ++npu){
    Double_t npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
    result[npu] = npu_estimated / npu_probs[npu];
    s += npu_estimated;
  }

  // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
  for(int npu=0; npu<maxNPU; ++npu){
    result[npu] /= s;
  }
  return result;
}

//______________________________________________________________________________
int ZbbUtils::getParentCode(const reco::GenParticle* genLep) {

  const reco::Candidate* particle = &(*genLep);
  // Find the actual lepton parent  
  int flavor = particle->pdgId();
  //  cout << " flavor " << flavor;
  while (particle->mother()!=0 && particle->mother()->pdgId() == flavor) {
    //    cout  << " " << particle->mother()->pdgId();
    particle = particle->mother();
  }
  //  cout << endl;
  int parentId = 0;
  if (particle->mother()!=0) {
    parentId = particle->mother()->pdgId();
    //if (parentId == 23 && particle->mother()->mother()!=0 && particle->mother()->mother()->pdgId() == 25) return 25;
  }
  return parentId;
}


//______________________________________________________________________________
std::pair<Double_t,Double_t> ZbbUtils::effCalc(const double& num, const double& den){
  
  Double_t eff=0.; 
  Double_t erreff=0.;

  std::pair<Double_t,Double_t> effResult = std::make_pair(eff,erreff);

  if ( den>0 ) {
      eff=num/den;
      erreff=TMath::Sqrt(eff*(1-eff)/den);
  }

  effResult.first=eff;
  effResult.second=erreff;

  return effResult;

}

//______________________________________________________________________________
std::vector<float> ZbbUtils::getPU(Bool_t getData){

  std::vector<float> dataPU(36);
  std::vector<float> MCPU(36);

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

  copy(data_2011_to_run173692_v3,data_2011_to_run173692_v3+36,dataPU.begin());
  copy(MC_summer11_PUS4,MC_summer11_PUS4+36,MCPU.begin());

  if(getData) return dataPU;
  else return MCPU;
  
}
