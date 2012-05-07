#ifndef ZbbAnalysis_AnalysisStep_EventCategory_h
#define ZbbAnalysis_AnalysisStep_EventCategory_h

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

struct EventCategory
{
  Bool_t MC_, Trigger_, LL_, LLTight_, LLTightAndTriggerMatch_, ZLL_, GoodJet_, JetBTagged_,JetBTaggedMCMatched_,ExclFinalState_, MET_, ThreeLeptons_, FourLeptons_,VertexAssoc_;
  reco::CompositeCandidate bestZcandidate_;
  reco::Vertex theZvertex_;
  Int_t    nBTags_;
  Double_t weight_;
  Double_t weightCorrectionFromInclToExcl_;
  // constructor
  EventCategory() { 
    //  std::cout << "EventCategory::EventCategory()"<< std::endl; 
    this->reset();
  }

  void reset(){
    MC_=false; Trigger_=false; LL_=false; LLTight_=false; LLTightAndTriggerMatch_=false; ZLL_=false; GoodJet_=false; JetBTagged_=false; JetBTaggedMCMatched_=false; ExclFinalState_=false; MET_=false; ThreeLeptons_=false, FourLeptons_=false,VertexAssoc_=false;
    nBTags_ = 0.;
    weight_ = 1.;
    weightCorrectionFromInclToExcl_=1;
  }
  // setters and getters
  
  Bool_t isMC()                       const {return MC_;}
  Bool_t isTrigger()                  const {return Trigger_;}
  Bool_t isLL()                       const {return LL_;}
  Bool_t isLLTight()                  const {return LLTight_;}
  Bool_t isLLTightAndTriggerMatched() const {return LLTightAndTriggerMatch_;}
  Bool_t isZLL()                      const {return ZLL_;}
  Bool_t isGoodJet()                  const {return GoodJet_;}
  Bool_t isJetBTagged()               const {return JetBTagged_;}
  Bool_t isJetBTaggedMCMatched()      const {return JetBTaggedMCMatched_;}
  Bool_t isExclusive()                const {return ExclFinalState_;}
  Bool_t isMET()                      const {return MET_;}
  Bool_t isThreeLeptons()             const {return ThreeLeptons_;}
  Bool_t isFourLeptons()              const {return FourLeptons_;}
  Bool_t isJetVertexAssoc()           const {return VertexAssoc_;}

  Int_t    getBTagNumber()            const {return nBTags_;}
  Double_t getWeight()                const {return weight_;}
  Double_t getBTagEffWeightCF()       const {return weightCorrectionFromInclToExcl_;}
};

#endif
