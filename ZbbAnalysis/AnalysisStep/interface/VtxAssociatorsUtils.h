#ifndef ZbbAnalysis_AnalysisStep_VtxAssociatorsUtils_H
#define ZbbAnalysis_AnalysisStep_VtxAssociatorsUtils_h

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "ZbbAnalysis/AnalysisStep/interface/EventCategory.h"

#include <TString.h>
#include <vector>

namespace VtxAssociatorsUtils {

  Bool_t checkVertexAssociation(const reco::CompositeCandidate& ZCand, edm::Handle<edm::View<pat::Jet> > jets, edm::View<reco::Vertex> vertices);
  reco::Vertex findPrimaryVertex(const reco::CompositeCandidate& ZCand,edm::View<reco::Vertex> vertices);
  
  Bool_t loosezVertex(const reco::CompositeCandidate& ZCand, Double_t cut);
  Bool_t tightzVertex(const reco::CompositeCandidate& ZCand, Double_t cut,const reco::Vertex& vertex);
  Bool_t tightzVertexMu (const pat::Muon& muon,         Double_t cut,const reco::Vertex& vertex);
  Bool_t tightzVertexEle(const pat::Electron& electron, Double_t cut,const reco::Vertex& vertex);
  
  Bool_t jetVertex(const reco::Vertex& vertex, const pat::Jet& jet, Int_t algo, Double_t sigmaCut, Double_t fraction);

  Bool_t jetVertex_1 (const reco::Vertex& vertex, const pat::Jet& jet, Double_t sigcut, Double_t etcut);
  Bool_t jetVertex_2 (const reco::Vertex& vertex, const pat::Jet& jet, Double_t sigcut, Double_t ptcut);
  Bool_t jetVertex_3 (const reco::Vertex& vertex, const pat::Jet& jet, Double_t sigcut, Double_t etcut);
  Bool_t jetVertex_2b(const reco::Vertex& vertex, const pat::Jet& jet, Double_t ptcut);
  Bool_t muVertex_2b (const reco::Vertex& vertex, const pat::Muon& muon);
  Bool_t eleVertex_2b(const reco::Vertex& vertex, const pat::Electron& electron);
}
#endif
